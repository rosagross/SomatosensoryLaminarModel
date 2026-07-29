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

    Fitting a subset of ROIs (env var ROI, default "all"):
        ROI=A3b python run_optimization.py       # score the error on BA3b only
        ROI=A3b,S2 python run_optimization.py    # average the error over BA3b and S2
    The full network is always simulated; ROI only selects which ROIs the objective
    is scored on (use base_params["area"] to isolate part of the network instead).
"""

import numpy as np
import pandas as pd
import os
import sys
import json
import mne
import matplotlib.pyplot as plt
from datetime import datetime

# ── paths ──────────────────────────────────────────────────────────────────────
WDDIR    = os.getenv("WDDIR")   # /data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel
RESDIR   = os.getenv("RESDIR")  # /data/p_02989/shared_workspace/results_grossmannr
SIMDIR   = os.getenv("SIMDIR")  # base output dir; optimization diagnostics go in SIMDIR/optimization

# ── ROI selection (env var ROI) ─────────────────────────────────────────────────
# Which ROI(s) the objective is scored on: "all" (default) averages the error over all
# three, otherwise a comma-separated subset, e.g. ROI=A3b or ROI=A3b,S2. The measured-CSV
# names BA3b/BA1 are accepted as aliases and matching is case-insensitive. Parsed up front
# so a typo fails before the model is built.
ALL_ROIS   = ("A3b", "A1", "S2")
_ROI_ALIAS = {"A3B": "A3b", "BA3B": "A3b", "A1": "A1", "BA1": "A1", "S2": "S2"}

_roi_env = (os.getenv("ROI") or "all").strip()
if _roi_env.lower() == "all":
    FIT_ROIS = ALL_ROIS
else:
    _selected = []
    for _tok in _roi_env.split(","):
        _tok = _tok.strip()
        if not _tok:
            continue
        if _tok.upper() not in _ROI_ALIAS:
            raise SystemExit(f"invalid ROI {_tok!r} (from ROI={_roi_env!r}); "
                             f"valid: {', '.join(ALL_ROIS)} (aliases BA3b, BA1), 'all', "
                             f"or a comma-separated subset")
        _selected.append(_ROI_ALIAS[_tok.upper()])
    if not _selected:
        raise SystemExit(f"empty ROI selection (ROI={_roi_env!r})")
    # keep ALL_ROIS order and drop duplicates, so ROI=S2,A3b == ROI=A3b,S2
    FIT_ROIS = tuple(r for r in ALL_ROIS if r in _selected)

# tag for output filenames/dirs
roi_tag = "all" if FIT_ROIS == ALL_ROIS else "-".join(FIT_ROIS)

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

import somato_model as sm
from somato_model import SomatoModel
from signal_preprocessing import smooth_timecourse, remove_aperiodic
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

# preprocess_targets repoints tc_data_path / ps_data_path at the processed CSVs below; keep
# the raw measured prestim spectrum around for plot_prestim_spectrum_comparison, which does
# its own 1/f removal and therefore needs the unflattened target.
ps_data_path_raw = ps_data_path

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

# loaded once: the (9, n_times) synthetic target dipole trace, or None for measured-data fit.
# Tolerate a missing file (e.g. before --generate-target has been run, or when doing a
# measured-data / --preprocess-only run) instead of crashing at import: warn and fall back
# to the measured-data path.
if SYNTHETIC_TARGET_PATH and os.path.exists(SYNTHETIC_TARGET_PATH):
    target_dip = model.load_dipole_trace(SYNTHETIC_TARGET_PATH)
else:
    if SYNTHETIC_TARGET_PATH:
        print(f"[run_optimization] SYNTHETIC_TARGET_PATH not found ({SYNTHETIC_TARGET_PATH}); "
              f"falling back to measured-data fit. Run --generate-target to create it.")
    target_dip = None


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


def preprocess_targets(tc_path, ps_path, outdir, fmin=1.0, fmax=40.0):
    """Preprocess the measured fitting targets before optimisation.

    - Timecourse: Savitzky-Golay smoothing per ROI (avoid overfitting sample noise).
    - Prestim spectrum: FOOOF removal of the aperiodic 1/f component per ROI, so the
      oscillatory peaks (not the slope) drive the error.

    Processed CSVs are written to `outdir` with the same schema as the originals
    (so the model's load_target_* readers work unchanged); the raw measured CSVs
    in RESDIR are never mutated. A before/after diagnostic PNG is saved per target.

    Note: the prestim `fmin/fmax` must match those used in
    compute_error_prestim_spectrum so target and simulated spectra are flattened
    over the same band.

    Returns:
        (tc_out_path, ps_out_path) processed CSV paths.
    """
    rois = ["BA3b", "BA1", "S2"]

    # --- timecourse: smooth the amplitude per ROI ---
    tc_df = pd.read_csv(tc_path)
    fig_tc, axes_tc = plt.subplots(1, 3, figsize=(15, 4), sharex=True)
    for ax, roi in zip(axes_tc, rois):
        mask = (tc_df["modality"] == "elec") & (tc_df["roi"] == roi)
        sub = tc_df[mask].sort_values("time_s")
        raw = sub["amplitude"].to_numpy()
        smoothed = smooth_timecourse(raw, window_ms=60.0)
        tc_df.loc[sub.index, "amplitude"] = smoothed
        ax.plot(sub["time_s"], raw, color="0.7", label="raw")
        ax.plot(sub["time_s"], smoothed, color="C3", label="smoothed")
        ax.set_title(roi); ax.set_xlabel("time (s)")
    axes_tc[0].legend(); axes_tc[0].set_ylabel("amplitude")
    fig_tc.suptitle("Timecourse target: raw vs smoothed"); fig_tc.tight_layout()
    fig_tc.savefig(os.path.join(outdir, "target_timecourse_smoothing.png"), dpi=200, bbox_inches="tight")
    plt.close(fig_tc)
    tc_out_path = os.path.join(outdir, "group_roi_timecourse_pooled_ses-elec_smoothed.csv")
    tc_df.to_csv(tc_out_path, index=False)

    # --- prestim spectrum: FOOOF 1/f removal per ROI ---
    ps_elec = pd.read_csv(ps_path)
    ps_elec = ps_elec[ps_elec["modality"] == "elec"]
    rows = []
    fig_ps, axes_ps = plt.subplots(1, 3, figsize=(15, 4))
    for ax, roi in zip(axes_ps, rois):
        sub = ps_elec[ps_elec["roi"] == roi].sort_values("freq_hz")
        freqs = sub["freq_hz"].to_numpy()
        power = sub["power"].to_numpy()
        f_flat, flat = remove_aperiodic(freqs, power, fmin, fmax)
        n_sub = int(sub["n_subjects"].iloc[0]) if "n_subjects" in sub.columns else 0
        for f, p in zip(f_flat, flat):
            rows.append({"modality": "elec", "roi": roi, "freq_hz": f, "power": p, "n_subjects": n_sub})
        ax.plot(freqs, power, color="0.7", label="raw")
        ax.plot(f_flat, flat, color="C0", label="1/f removed")
        ax.set_title(roi); ax.set_xlabel("freq (Hz)"); ax.set_xlim(0, fmax)
    axes_ps[0].legend(); axes_ps[0].set_ylabel("power")
    fig_ps.suptitle("Prestim spectrum target: raw vs 1/f removed"); fig_ps.tight_layout()
    fig_ps.savefig(os.path.join(outdir, "target_prestim_spectrum_1fremoved.png"), dpi=200, bbox_inches="tight")
    plt.close(fig_ps)
    ps_out_path = os.path.join(outdir, "group_roi_prestim_spectrum_ses-elec_1fremoved.csv")
    pd.DataFrame(rows).to_csv(ps_out_path, index=False)

    return tc_out_path, ps_out_path


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
        err_tf, _, _ = model.compute_error_timefreq(tf_data_path, sim_dip, target_dip=target_dip, rois=FIT_ROIS)
    if ERROR_MODE in ("tc", "both", "all"):
        err_tc, _, _ = model.compute_error_timecourse(tc_data_path, sim_dip, target_dip=target_dip, rois=FIT_ROIS)
    if ERROR_MODE in ("ps", "all"):
        # flatten the simulated spectrum only for measured-data fits, where the target
        # CSV has been 1/f-removed; the synthetic target (target_dip) is left raw.
        err_ps, _, _ = model.compute_error_prestim_spectrum(
            ps_data_path, sim_dip, target_dip=target_dip, flatten_sim=(target_dip is None), rois=FIT_ROIS)
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
        "scaling_factor",
        "e3b_tau",
        "e1_tau",
        "e2_tau",
        "thal_delay_factor",
        "delay_factor",
        "receptor_thalamus_delay"
    ],
    "bounds": np.array([
        [0,     50  ],   # coupling_strength
        [0.4,      0.8],   # strength_I
        [0.5,      2  ],   # g_intercortical
        [0,      2  ],   # g_thalPOm (scales POm output connectivity)
        [3,     10  ],   # Ib_strength
        [0,    100  ],   # Iext_strength
        [0.001,  0.05],   # Iext_duration
        [0.5,      1.1],   # scaling factor
        [2,     10  ],   # e3b_tau (ms, default 6)
        [2,     10  ],   # e1_tau  (ms, default 6)
        [2,     10  ],   # e2_tau  (ms, default 6)
        [0.001,  0.005],  # thal_delay_factor (s, default 3e-3)
        [0.001,  0.008],  # delay_factor      (s, default 5e-3)
        [0.010,  0.03],   # receptor_thalamus_delay (s, default 0.050; ~N20 latency minus thalamus to cortex and synaptic delay)
    ]),
    "reference":  0.0,
    "simulation": objective,
    "op":         -1,    # minimise
    "N1":         30, #30,    # initial population size
    "N2":         30, #40,    # crossover offspring per iteration
    "N3":         30, #40,    # mutation offspring per iteration
    "n_iter":     10,
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

    # errors are scored on FIT_ROIS (as the GA did), while the comparison figures below
    # still show all three ROIs — useful to see what the fit did to the unfitted ones.
    err_tf, tf_sim, tf_target = model.compute_error_timefreq(tf_data_path, sim_dip, target_dip=target_dip, rois=FIT_ROIS)
    model.plot_timefreq_comparison(tf_sim, tf_target)        # saves to the model's TIMEFREQ_DIR

    err_tc, tc_sim, tc_target = model.compute_error_timecourse(tc_data_path, sim_dip, target_dip=target_dip, rois=FIT_ROIS)
    model.plot_timecourse_comparison(tc_sim, tc_target)      # saves to the model's TIMECOURSE_DIR

    err_ps, ps_sim, ps_target = model.compute_error_prestim_spectrum(
        ps_data_path, sim_dip, target_dip=target_dip, flatten_sim=(target_dip is None), rois=FIT_ROIS)
    # the figure flattens both sides itself, so it gets the *raw* measured CSV (ps_data_path
    # points at the already-flattened one after preprocess_targets); its bottom row then
    # reproduces the flattened comparison the error above scores.
    model.plot_prestim_spectrum_comparison(ps_data_path_raw, sim_dip, target_dip=target_dip)  # saves to PRESTIM_SPECTRUM_DIR

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

    if "--preprocess-only" in sys.argv:
        # Run only the target preprocessing (smooth timecourse + 1/f-remove prestim
        # spectrum), write the processed CSVs + diagnostic PNGs, then exit — no GA.
        preprocess_targets(tc_data_path, ps_data_path, diag_dir)
        print(f"Preprocessed targets written to {diag_dir}")
        sys.exit(0)

    # Each optimization run gets its own timestamped session folder so runs no
    # longer overwrite each other. All figures/CSVs for this run go here.
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    session_dir = os.path.join(diag_dir, f"opt_{run_id}_{ERROR_MODE}_roi-{roi_tag}")
    os.makedirs(session_dir, exist_ok=True)
    print(f"Optimization session: {session_dir}")
    print(f"Fitting ROIs: {', '.join(FIT_ROIS)}  (ROI={_roi_env!r})")

    # Log the full run configuration up front (persisted even if the run crashes).
    # target_csvs are recorded before preprocess_targets repoints them below, so
    # the log keeps the raw measured sources.
    run_config = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "error_mode": ERROR_MODE,
        "fit_rois": list(FIT_ROIS),
        "mode": "synthetic" if target_dip is not None else "measured",
        "base_params": base_params,
        "optimized_parameters": opt_config["model_parameters"],
        "bounds": np.asarray(opt_config["bounds"]).tolist(),
        "ga": {k: opt_config[k] for k in ("N1", "N2", "N3", "n_iter", "tolerance", "op", "reference")},
        "subID_elec": subID_elec,
        "target_csvs": {"tf": tf_data_path, "tc": tc_data_path, "ps": ps_data_path},
        "synthetic_target_path": SYNTHETIC_TARGET_PATH if target_dip is not None else None,
        "true_params": TRUE_PARAMS if target_dip is not None else None,
    }
    with open(os.path.join(session_dir, "run_config.json"), "w") as f:
        json.dump(run_config, f, indent=2)

    # Redirect the model's comparison-figure dirs (read as module globals at call
    # time) into this session so the best-fit comparison figures land here too.
    for _attr, _name in [("TIMEFREQ_DIR", "timefreq_comparison"),
                         ("TIMECOURSE_DIR", "timecourse_comparison"),
                         ("PRESTIM_SPECTRUM_DIR", "prestim_spectrum_comparison")]:
        _d = os.path.join(session_dir, _name)
        os.makedirs(_d, exist_ok=True)
        setattr(sm, _attr, _d)

    if target_dip is None:
        # Preprocess the measured targets before fitting: smooth the timecourse
        # and remove the aperiodic 1/f component from the prestim spectrum. Repoint
        # the (module-global) target paths at the processed CSVs so objective() and
        # plot_best_fit() fit against them. Skipped in synthetic-target mode, which
        # fits against target_dip rather than these CSVs.
        tc_data_path, ps_data_path = preprocess_targets(tc_data_path, ps_data_path, session_dir)
        print(f"Preprocessed targets written to {session_dir}")

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
    plot_fit_diagnostics(ga, opt_config, session_dir)

    # visualise (and persist) the best-fit simulation against the measured data
    err_tf, err_tc, err_ps = plot_best_fit(model, best_params, session_dir)
    print(f"  Best-fit re-evaluation: err_tf={err_tf:.4f}  err_tc={err_tc:.4f}  err_ps={err_ps:.4f}")

    # write a small machine-readable summary of the run
    summary = {
        "error_mode": ERROR_MODE,
        "fit_rois": list(FIT_ROIS),
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
    with open(os.path.join(session_dir, "optimization_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  Diagnostics written to {session_dir}")
