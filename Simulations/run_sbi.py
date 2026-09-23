"""
File: run_sbi.py
Author: Rosa Grossmann
Description:
    Simulation-based inference (SBI) of SomatoModel parameters against measured
    electrical-stimulation EEG data, as an alternative to the genetic-algorithm
    fit in run_optimization.py.

    Instead of minimising a scalar error, this learns the posterior p(theta | x_obs)
    over the 7 free parameters using multi-round NPE (SNPE), where x is the
    summary-statistic vector defined in sbi_features.py.

    Quick local smoke test:
        python run_sbi.py --num-rounds 1 --num-sims 50 --workers 4

    Full run (HPC, see hpc_scripts/run_sbi.sh):
        python run_sbi.py --num-rounds 3 --num-sims 3000 --workers $SLURM_CPUS_PER_TASK

    Synthetic parameter-recovery sanity check (no measured data needed):
        python run_sbi.py --num-rounds 2 --num-sims 300 --workers 4 --synthetic-obs
"""

import argparse
import contextlib
import os
import pickle
import sys
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

# ── paths (mirror run_optimization.py) ───────────────────────────────────────────
WDDIR  = os.getenv("WDDIR")
RESDIR = os.getenv("RESDIR")

# Where finished SBI runs are stored (each run -> RESULTS_ROOT/sbi_<timestamp>/).
# Overridable via env var; --outdir still takes precedence over this default root.
RESULTS_ROOT = os.getenv("SBI_RESULTS_ROOT", "/data/pt_02989/optimization_results")

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(_HERE)                          # Simulations/ (for parameters.py)
sys.path.append(os.path.join(_HERE, "model"))   # Simulations/model/ (for somato_model.py)
if WDDIR:
    sys.path.append(os.path.join(WDDIR, "Simulations"))
    sys.path.append(os.path.join(WDDIR, "Simulations", "model"))

from somato_model import SomatoModel, read_simulation_params   # noqa: E402
from sbi_features import (                                 # noqa: E402
    extract_summary_features, feature_names, N_FEATURES,
)

from sbi.inference import NPE, simulate_for_sbi            # noqa: E402
from sbi.utils import BoxUniform, process_prior, process_simulator  # noqa: E402
from sbi.analysis import pairplot                          # noqa: E402

# ── target data paths ────────────────────────────────────────────────────────────
_roi_dir = os.path.join(
    RESDIR or "", "Figures", "Main", "eeg_results", "source_reconstruction",
    "group", "_preprestim_corrected", "roi_epochswise",
)
TF_DATA_PATH = os.path.join(_roi_dir, "group_roi_tf_morlet_ses-elec_preprestim_corrected.csv")
TC_DATA_PATH = os.path.join(_roi_dir, "group_roi_timecourse_pooled_ses-elec_preprestim_corrected.csv")

# ── parameters to infer + prior bounds (identical to the GA in run_optimization.py) ─
PARAM_NAMES = [
    "coupling_strength",
    "strength_I",
    "g_intercortical",
    "g_thalPOm",
    "Ib_strength",
    # background input onto each interneuron class, relative to the drive onto E.
    # Ib_ratio_E is not inferred: it is the reference the others are defined against, so
    # it is exactly collinear with Ib_strength. It stays at its BASE_PARAMS value.
    "Ib_ratio_PV",
    "Ib_ratio_SST",
    "Ib_ratio_VIP",
    "Iext_strength",
    "Iext_duration",
]
BOUNDS = np.array([
    [0,     50  ],   # coupling_strength
    [0,      0.5],   # strength_I
    [0,      2  ],   # g_intercortical
    [0,      2  ],   # g_thalPOm (scales POm output connectivity)
    [0,     10  ],   # Ib_strength
    [0.2,    1.5],   # Ib_ratio_PV
    [0.2,    1.5],   # Ib_ratio_SST
    [0.1,    1.5],   # Ib_ratio_VIP
    [0,    100  ],   # Iext_strength
    [0.001,  0.1],   # Iext_duration
], dtype=float)

# Defaults from Simulations/simulation_parameter.json (the single source of default
# parameters); only the deviations shared with the GA in run_optimization.py are set here,
# so both fitting entry points evaluate the same fixed model.
BASE_PARAMS = read_simulation_params()
BASE_PARAMS.update({
    "input_onset":      2.001,   # settled pre-stimulus period, as in run_optimization.py
    "simulation_dur":   3,
    "Ib_noise_std":     0.0,     # background noise OFF, as in run_optimization.py
    "Ib_noise_seed":    0,       # inert while Ib_noise_std == 0; makes the sim deterministic
})

# One model instance per process (workers each import this module afresh).
_MODEL = SomatoModel(BASE_PARAMS)

# Finite sentinel returned when a simulation blows up / yields non-finite features,
# so SBI training never sees NaNs. Far outside the plausible feature range.
_BAD_FEATURES = np.full(N_FEATURES, -999.0, dtype=float)


def simulator(theta):
    """theta (7,) -> summary-feature vector (N_FEATURES,). Robust to blow-ups."""
    theta = np.asarray(theta, dtype=float).reshape(-1)
    params = dict(zip(PARAM_NAMES, theta))
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            _MODEL.apply_params(params)
            _MODEL.initialize_state()
            _MODEL.simulate()
            sim_dip = _MODEL.compute_dipoles()
            tf = _MODEL.compute_timefreq(simulated_dip=sim_dip)
            tc = _MODEL.compute_timecourse(simulated_dip=sim_dip)
            x = extract_summary_features(tf, tc, _MODEL.step_size)
    except Exception:
        return _BAD_FEATURES.copy()
    if not np.all(np.isfinite(x)):
        return _BAD_FEATURES.copy()
    return x


def build_x_obs():
    """Build the observation vector from the measured electrical-stimulation CSVs."""
    tf_target = _MODEL.load_target_timefreq(TF_DATA_PATH)
    tc_target = _MODEL.load_target_timecourse(TC_DATA_PATH)
    return extract_summary_features(tf_target, tc_target, _MODEL.step_size)


def run(num_rounds, num_sims, workers, seed, outdir, synthetic_obs):
    os.makedirs(outdir, exist_ok=True)
    torch.manual_seed(seed)
    np.random.seed(seed)

    low  = torch.as_tensor(BOUNDS[:, 0], dtype=torch.float32)
    high = torch.as_tensor(BOUNDS[:, 1], dtype=torch.float32)
    prior = BoxUniform(low=low, high=high)
    prior, _, prior_returns_numpy = process_prior(prior)
    sim = process_simulator(simulator, prior, prior_returns_numpy)

    # ── observation ──────────────────────────────────────────────────────────────
    if synthetic_obs:
        theta_true = prior.sample((1,))
        x_obs = torch.as_tensor(simulator(theta_true.numpy()[0]), dtype=torch.float32)
        np.save(os.path.join(outdir, "theta_true.npy"), theta_true.numpy()[0])
        print("Synthetic ground-truth theta:",
              dict(zip(PARAM_NAMES, theta_true.numpy()[0].round(4))))
    else:
        x_obs = torch.as_tensor(build_x_obs(), dtype=torch.float32)
        theta_true = None
    np.save(os.path.join(outdir, "x_obs.npy"), x_obs.numpy())
    print(f"x_obs ({len(x_obs)} features):", np.array2string(x_obs.numpy(), precision=3))

    # ── multi-round NPE (SNPE) ─────────────────────────────────────────────────────
    inference = NPE(prior=prior)
    proposal = prior
    all_theta, all_x = [], []
    for r in range(num_rounds):
        print(f"\n── Round {r + 1}/{num_rounds}: simulating {num_sims} ──")
        theta, x = simulate_for_sbi(
            sim, proposal, num_simulations=num_sims,
            num_workers=workers, seed=seed + r,
        )
        all_theta.append(theta.numpy())
        all_x.append(x.numpy())
        n_bad = int((x == torch.as_tensor(_BAD_FEATURES, dtype=x.dtype)).all(dim=1).sum())
        print(f"  {n_bad}/{num_sims} simulations failed (sentinel).")

        density_estimator = inference.append_simulations(
            theta, x, proposal=proposal,
        ).train()
        posterior = inference.build_posterior(density_estimator)
        posterior.set_default_x(x_obs)
        proposal = posterior

    # ── persist ────────────────────────────────────────────────────────────────────
    with open(os.path.join(outdir, "posterior.pkl"), "wb") as f:
        pickle.dump(posterior, f)
    np.savez(
        os.path.join(outdir, "simulations.npz"),
        theta=np.concatenate(all_theta), x=np.concatenate(all_x),
        x_obs=x_obs.numpy(), param_names=np.array(PARAM_NAMES),
        feature_names=np.array(feature_names()),
    )

    # ── posterior summary ──────────────────────────────────────────────────────────
    samples = posterior.sample((10000,), x=x_obs)
    s = samples.numpy()
    # Approximate MAP = highest-log-prob posterior sample. Robust across sbi/torch
    # versions, unlike posterior.map()'s gradient ascent (which trips a transform
    # _inv bug with a BoxUniform prior on recent torch).
    log_probs = posterior.log_prob(samples, x=x_obs)
    map_est = samples[log_probs.argmax()].numpy()

    summary_path = os.path.join(outdir, "posterior_summary.csv")
    with open(summary_path, "w") as f:
        f.write("parameter,map,mean,std,p5,p50,p95,prior_low,prior_high\n")
        for i, name in enumerate(PARAM_NAMES):
            p5, p50, p95 = np.percentile(s[:, i], [5, 50, 95])
            f.write(f"{name},{map_est[i]:.6g},{s[:, i].mean():.6g},{s[:, i].std():.6g},"
                    f"{p5:.6g},{p50:.6g},{p95:.6g},{BOUNDS[i, 0]:.6g},{BOUNDS[i, 1]:.6g}\n")
    print(f"\nPosterior summary written to {summary_path}")

    # ── pairplot ─────────────────────────────────────────────────────────────────────
    points = theta_true.numpy() if theta_true is not None else None
    fig, _ = pairplot(
        samples, labels=PARAM_NAMES,
        limits=[[lo, hi] for lo, hi in BOUNDS],
        points=points, figsize=(11, 11),
    )
    fig.savefig(os.path.join(outdir, "posterior_pairplot.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Results saved to {outdir}")
    return posterior


def parse_args():
    p = argparse.ArgumentParser(description="SBI (multi-round NPE) fit of SomatoModel.")
    p.add_argument("--num-rounds", type=int, default=2)
    p.add_argument("--num-sims",   type=int, default=500, help="simulations per round")
    p.add_argument("--workers",    type=int, default=1)
    p.add_argument("--seed",       type=int, default=0)
    p.add_argument("--outdir",     type=str, default=None)
    p.add_argument("--synthetic-obs", action="store_true",
                   help="use a simulated observation (parameter-recovery sanity check)")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.outdir is None:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.outdir = os.path.join(RESULTS_ROOT, f"sbi_{stamp}")
    run(args.num_rounds, args.num_sims, args.workers, args.seed,
        args.outdir, args.synthetic_obs)
