"""
File: run_optimization.py
Author: Rosa Grossmann
Description:
    Fit SomatoModel parameters by minimising the combined time-frequency and
    time-course error against measured electrical-stimulation data.

    Uses the GA optimizer from the neuronaldynamics package.

    Quick sanity-check run (~25 min):
        Set N1=10, N2=20, N3=20, n_iter=5 below.

    Full run (~4 h):
        Set N1=20, N2=40, N3=40, n_iter=20 below.
"""

import numpy as np
import os
import sys
import json
import mne
import matplotlib.pyplot as plt

# ── paths ──────────────────────────────────────────────────────────────────────
WDDIR    = os.getenv("WDDIR")   # /data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel
RESDIR   = os.getenv("RESDIR")  # /data/p_02989/shared_workspace/results_grossmannr
SIMDIR   = os.getenv("SIMDIR")  # base output dir; optimization diagnostics go in SIMDIR/optimization

# electrical-modality subjects; compute_dipoles reads each one's forward model
# (same list as Analysis/SourceReconstruction/step002_inverse_solution_multisub_epochswise.py)
subID_elec = [15, 16, 17, 18, 23, 24, 25, 26, 27, 28, 29, 34, 35, 36, 37, 38, 39, 40,
              42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52]

# where to store optimization diagnostics (fit comparison plots use the model's own dirs)
diag_dir = os.path.join(SIMDIR, "optimization")
os.makedirs(diag_dir, exist_ok=True)

sys.path.append(os.path.join(WDDIR, "Simulations", "model"))
sys.path.insert(0, "/data/p_02989/Modelling/neuronaldynamics/src")
# Optimizer.py imports Utils via a bare `from Utils import ...`, so its containing
# directory must be on sys.path (Utils.py lives at src/neuronaldynamics/Utils.py).
sys.path.insert(0, "/data/p_02989/Modelling/neuronaldynamics/src/neuronaldynamics")

from somato_model import SomatoModel
from neuronaldynamics.Optimizers.Optimizer import GA

# neuronaldynamics.Utils forces matplotlib's interactive 'TkAgg' backend on import,
# whose plt.show() blocks on a window when run headless (e.g. over SSH on the
# cluster). For a non-interactive run, switch to the Agg backend *after* that
# import: figures are still saved via savefig but plt.show() becomes a no-op.
# Set OPT_HEADLESS=1 to force this even under a virtual display (xvfb, where
# DISPLAY is set so TkAgg can import). Interactive runs leave both unset and keep
# TkAgg so windows still pop up.
import matplotlib
if os.environ.get("OPT_HEADLESS") or not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")

# ── target data paths ──────────────────────────────────────────────────────────
_roi_dir = os.path.join(
    RESDIR, "Figures", "Main", "eeg_results", "source_reconstruction",
    "group", "_preprestim_corrected", "roi_epochswise",
)
tf_data_path = os.path.join(_roi_dir, "group_roi_tf_morlet_ses-elec_preprestim_corrected.csv")
tc_data_path = os.path.join(_roi_dir, "group_roi_timecourse_pooled_ses-elec_preprestim_corrected.csv")
ps_data_path = os.path.join(_roi_dir, "group_roi_prestim_spectrum_ses-elec_preprestim_corrected.csv")

# ── model (fixed non-optimised parameters) ─────────────────────────────────────
base_params = {
    "g_thal":           2,
    "sI_thal":          0.5,
    "extI_cellcounts":  1000,
    "bI_cellcounts":    100,
    "thal_cellcounts":  500,
    "area":             "all",
}
model = SomatoModel(base_params)

# ── objective function ─────────────────────────────────────────────────────────
# Which error(s) the GA optimizes against: "tf" (time-frequency), "tc" (time-course),
# "ps" (pre-stim spectrum) individually; "both" (tf+tc); or "all" (tf+tc+ps).
ERROR_MODE = "tc"
assert ERROR_MODE in ("tf", "tc", "ps", "both", "all"), f"invalid ERROR_MODE: {ERROR_MODE!r}"

# ── parameter-recovery mode ─────────────────────────────────────────────────────
# If SYNTHETIC_TARGET_PATH is set, the GA optimizes against a synthetic dipole trace
# saved by model.save_dipole_trace (instead of the measured CSVs), so we can check
# whether it recovers the known TRUE_PARAMS. None = normal measured-data fit.
# Generate the target file once with:  python run_optimization.py --generate-target
SYNTHETIC_TARGET_PATH = "/data/p_02989/Modelling/output_grossmannr/optimization/synthetic_target.hdf5"  # e.g. os.path.join(diag_dir, "synthetic_target.hdf5")
TRUE_PARAMS = {               # params used to generate the synthetic target
    "coupling_strength": 10, "strength_I": 0.68, "g_intercortical": 1.0, "g_thalPOm": 1.0,
    "Ib_strength": 6, "Iext_strength": 40, "Iext_duration": 0.016, "scaling_factor": 1.0,
}

# loaded once: the (9, n_times) synthetic target dipole trace, or None for measured-data fit
target_dip = model.load_dipole_trace(SYNTHETIC_TARGET_PATH) if SYNTHETIC_TARGET_PATH else None


def generate_synthetic_target():
    """Run the model at TRUE_PARAMS and save its dipole trace as a synthetic target."""
    model.apply_params(TRUE_PARAMS)
    model.initialize_state()
    model.simulate()
    sim_dip = model.compute_dipoles(subID_elec)
    path = model.save_dipole_trace(sim_dip, diag_dir, filename="synthetic_target")
    print(f"Synthetic target saved to {path}")
    print(f"TRUE_PARAMS: {TRUE_PARAMS}")
    return path


def objective(**params):
    """
    Run one full simulation and return the selected error (see ERROR_MODE):
    TF, timecourse, pre-stim spectrum, or their sum.
    The GA minimises (0 - objective)**2, i.e. the squared selected error.
    When target_dip is set, errors are computed against the synthetic target instead
    of the measured CSVs.
    """
    model.apply_params(params)
    model.initialize_state()
    model.simulate()
    sim_dip = model.compute_dipoles(subID_elec)
    err_tf = err_tc = err_ps = 0.0
    if ERROR_MODE in ("tf", "both", "all"):
        err_tf, _, _ = model.compute_error_timefreq(tf_data_path, sim_dip, target_dip=target_dip)
    if ERROR_MODE in ("tc", "both", "all"):
        err_tc, _, _ = model.compute_error_timecourse(tc_data_path, sim_dip, target_dip=target_dip)
    if ERROR_MODE in ("ps", "all"):
        err_ps, _, _ = model.compute_error_prestim_spectrum(ps_data_path, sim_dip, target_dip=target_dip)
    combined = err_tf + err_tc + err_ps
    print(f"  params={params}  →  err_tf={err_tf:.4f}  err_tc={err_tc:.4f}  err_ps={err_ps:.4f}  total={combined:.4f}")
    return combined

# ── GA setup ───────────────────────────────────────────────────────────────────
opt_config = {
    "model_parameters": [
        "coupling_strength",
        "strength_I",
        "g_intercortical",
        "g_thalPOm",
        "Ib_strength",
        "Iext_strength",
        "Iext_duration",
        "scaling_factor"
    ],
    "bounds": np.array([
        [0,     50  ],   # coupling_strength
        [0.5,      0.8],   # strength_I
        [0.5,      2  ],   # g_intercortical
        [0,      2  ],   # g_thalPOm (scales POm output connectivity)
        [3,     10  ],   # Ib_strength
        [0,    100  ],   # Iext_strength
        [0.001,  0.1],   # Iext_duration
        [0.5,      1.1],   # scaling factor
    ]),
    "reference":  0.0,
    "simulation": objective,
    "op":         -1,    # minimise
    "N1":         20,    # initial population size
    "N2":         40,    # crossover offspring per iteration
    "N3":         40,    # mutation offspring per iteration
    "n_iter":     20,
    "tolerance":  0.05,
    "verbose":    1,
}

# ── diagnostics from the fitting process ────────────────────────────────────────
def plot_fit_diagnostics(ga, config, outdir):
    """
    Visualise how the GA progressed: best combined error per iteration and the
    evolution of each fitted parameter (normalised to its [0, 1] bound range so
    the differently-scaled parameters are comparable on one axis).

    ga.errors               : list, best combined error per iteration
    ga.parameter_evolution  : (n_iter, n_param) best parameter set per iteration
    """
    names  = config["model_parameters"]
    bounds = np.asarray(config["bounds"], dtype=float)
    errors = np.asarray(ga.errors, dtype=float)
    evo    = np.asarray(ga.parameter_evolution, dtype=float)   # (n_iter, n_param)

    fig, (ax_err, ax_par) = plt.subplots(1, 2, figsize=(12, 4.2))

    # convergence
    ax_err.plot(np.arange(1, len(errors) + 1), errors, marker="o", color="C3")
    ax_err.set_xlabel("iteration")
    ax_err.set_ylabel("best combined error")
    ax_err.set_yscale("log")
    ax_err.set_title("Convergence")

    # parameter evolution, normalised within bounds
    span = (bounds[:, 1] - bounds[:, 0])
    span[span == 0] = 1.0
    evo_norm = (evo - bounds[:, 0]) / span
    iters = np.arange(1, evo.shape[0] + 1)
    for j, name in enumerate(names):
        ax_par.plot(iters, evo_norm[:, j], marker="o", label=name)
    ax_par.set_ylim(-0.05, 1.05)
    ax_par.set_xlabel("iteration")
    ax_par.set_ylabel("value (normalised to bounds)")
    ax_par.set_title("Parameter evolution")
    ax_par.legend(fontsize=8, loc="best")

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fit_diagnostics.png"), dpi=300, bbox_inches="tight")
    plt.show()


def plot_best_fit(model, best_params, outdir):
    """Re-run the model at the optimised parameters and plot/save the fit vs. data."""
    model.apply_params(best_params)
    model.initialize_state()
    model.simulate()

    sim_dip = model.compute_dipoles(subID_elec)

    err_tf, tf_sim, tf_target = model.compute_error_timefreq(tf_data_path, sim_dip, target_dip=target_dip)
    model.plot_timefreq_comparison(tf_sim, tf_target)        # saves to the model's TIMEFREQ_DIR

    err_tc, tc_sim, tc_target = model.compute_error_timecourse(tc_data_path, sim_dip, target_dip=target_dip)
    model.plot_timecourse_comparison(tc_sim, tc_target)      # saves to the model's TIMECOURSE_DIR

    err_ps, ps_sim, ps_target = model.compute_error_prestim_spectrum(ps_data_path, sim_dip, target_dip=target_dip)
    model.plot_prestim_spectrum_comparison(ps_data_path, sim_dip, target_dip=target_dip)  # saves to PRESTIM_SPECTRUM_DIR

    # also persist the comparison maps/traces for this best run
    model.save_timefreq_comparison(outdir, tf_sim, tf_target, err_tf, filename="best_tf_comparison")
    model.save_timecourse_comparison(outdir, tc_sim, tc_target, err_tc, filename="best_tc_comparison")
    return err_tf, err_tc, err_ps


# ── run ────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Parameter-recovery helper: generate the synthetic target once, then exit.
    if "--generate-target" in sys.argv:
        generate_synthetic_target()
        sys.exit(0)

    ga = GA(opt_config)
    ga.run()

    best_params = dict(zip(opt_config["model_parameters"], ga.optimum))
    print("\n── Optimised parameters ──")
    for name, val in best_params.items():
        print(f"  {name}: {val:.4f}")
    print(f"  Best combined error: {ga.errors[-1]:.4f}")

    # in parameter-recovery mode, show the known TRUE_PARAMS alongside the recovered ones
    if target_dip is not None:
        print("\n── True parameters (synthetic target) ──")
        for name in opt_config["model_parameters"]:
            true_val = TRUE_PARAMS.get(name)
            true_str = f"{true_val:.4f}" if true_val is not None else "—"
            print(f"  {name}: true={true_str}  recovered={best_params[name]:.4f}")

    # diagnostics of the optimisation itself
    ga.plot_fit()
    plot_fit_diagnostics(ga, opt_config, diag_dir)

    # visualise (and persist) the best-fit simulation against the measured data
    err_tf, err_tc, err_ps = plot_best_fit(model, best_params, diag_dir)
    print(f"  Best-fit re-evaluation: err_tf={err_tf:.4f}  err_tc={err_tc:.4f}  err_ps={err_ps:.4f}")

    # write a small machine-readable summary of the run
    summary = {
        "error_mode": ERROR_MODE,
        "best_params": {k: float(v) for k, v in best_params.items()},
        "best_combined_error": float(ga.errors[-1]),
        "error_per_iteration": [float(e) for e in ga.errors],
        "best_fit_err_tf": float(err_tf),
        "best_fit_err_tc": float(err_tc),
        "best_fit_err_ps": float(err_ps),
    }
    if target_dip is not None:
        summary["synthetic_target_path"] = SYNTHETIC_TARGET_PATH
        summary["true_params"] = {k: float(v) for k, v in TRUE_PARAMS.items()}
    with open(os.path.join(diag_dir, "optimization_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  Diagnostics written to {diag_dir}")
