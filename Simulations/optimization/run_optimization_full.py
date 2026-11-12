"""
File: run_optimization_full.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-11-05
Description: Run an optimization run using the hybrid genetic algorithm. 
In this script we will be using the full 
"""

import argparse
import ast
import json
import os
import sys
import time
import types
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import nnls, lsq_linear, minimize

# ── paths ──────────────────────────────────────────────────────────────────────
WDDIR  = os.getenv("WDDIR")   # /data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel
SIMDIR = os.getenv("SIMDIR")  # /data/p_02989/Modelling/output_grossmannr
OPTDIR   = os.path.dirname(os.path.abspath(__file__))
RESDIR   = os.getenv("RESDIR")  # /data/p_02989/shared_workspace/results_grossmannr

sys.path.insert(0, OPTDIR)
sys.path.append(os.path.join(WDDIR, "Simulations", "model"))
sys.path.insert(0, "/data/p_02989/Modelling/neuronaldynamics/src")
# Optimizer.py imports Utils via a bare `from Utils import ...`, so its containing
# directory must be on sys.path (same setup as run_optimization.py).
sys.path.insert(0, "/data/p_02989/Modelling/neuronaldynamics/src/neuronaldynamics")

# import scripts from other folders
import dipole_parameterization as dp
import somato_model as sm
from neuronaldynamics.Optimizers.Optimizer import GA

ROIS = ("A3b", "A1", "S2")


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

# where to store optimization diagnostics (fit comparison plots use the model's own dirs)
diag_dir = os.path.join(SIMDIR, "optimization_timecourse")
os.makedirs(diag_dir, exist_ok=True)

# ── target data paths ──────────────────────────────────────────────────────────
_roi_dir = os.path.join(
    RESDIR, "Figures", "Main", "eeg_results", "source_reconstruction",
    "group", "_preprestim_corrected", "roi_epochswise",
)
tc_data_path = os.path.join(_roi_dir, "group_roi_timecourse_pooled_ses-elec_preprestim_corrected.csv")


# ── parameters that can be fitted ────────────────────────────────────────────────
PARAMETER_DICT = {
    # DELAYS
    "thal_delay_factor":       dict(enabled=False, abs=(0.005, 0.01), default=0.007),
    "delay_factor":            dict(enabled=True, abs=(0.005, 0.015), default=0.011),
    "delay_factor_short":      dict(enabled=False, abs=(0.001, 0.006), default=0.0015),
    "receptor_thalamus_delay": dict(enabled=False, abs=(0.005, 0.025), default=0.01),
    
    # CONNECTIVITY & COUPLING 
    "coupling_strength":       dict(enabled=True, abs=(0.5, 8), default=1), 
    "strength_I":              dict(enabled=True, abs=(0.5, 1.4), default=1), 
    "g_intercortical":         dict(enabled=True, abs=(0.2, 3), default=1),
    # intra-cortical
    "f_intra_S2":              dict(enabled=True, abs=(1, 5.0), default=1),
    "f_intra":                 dict(enabled=True, abs=(0.5, 3.0), default=1),
    "p_2PVE":                  dict(enabled=False, abs=(37.4, 37.4), default=37.6),
    "p_4PVE":                  dict(enabled=False, abs=(39.6, 39.6), default=39.4),
    # feedforward
    "ff_A1_S2":                dict(enabled=True, abs=(0.5, 5.0), default=1),
    "ff_A3b_S2":               dict(enabled=True, abs=(0.5, 2.0), default=1),
    "ff_A3b_A1":               dict(enabled=True, abs=(0.5, 2.0), default=1),
    # feedback
    "fb_A1_S2":                dict(enabled=True, abs=(0.2, 5.0), default=1),
    "fb_A3b_A1":               dict(enabled=True, abs=(0.5, 5.0), default=1),
    "fb_A3b_S2":               dict(enabled=True, abs=(0.5, 5.0), default=1),
    # thalamus
    "g_thalPOm":               dict(enabled=False, abs=(0.5, 3.0), default=1),
    "g_thal":                  dict(enabled=False, abs=(0.5, 3.0), default=2),
    "POm_drive_A3b":           dict(enabled=False, abs=(0.5, 2.0), default=1),
    "POm_drive_A1":            dict(enabled=False, abs=(0.5, 5.0), default=1),
    "POm_drive_S2":            dict(enabled=True, abs=(0.5, 2.0), default=1),
    "POm_tgt_A3b":             dict(enabled=False, abs=(0.5, 2.0), default=1),
    "POm_tgt_A1":              dict(enabled=False, abs=(0.5, 2.0), default=1),
    "POm_tgt_S2":              dict(enabled=False, abs=(0.5, 5.0), default=1),
    
    # TIME CONSTANTS
    "e3b_tau":                 dict(enabled=True, abs=(4.0, 10), default=6),
    "e1_tau":                  dict(enabled=True, abs=(4.0, 10), default=6),
    "e2_tau":                  dict(enabled=True, abs=(4.0, 10), default=6),
    # NMDA decay (s) of E -> cell class connections, by postsynaptic class
    "nmda_tau_E":              dict(enabled=True, abs=(5e-3, 20e-3), default=10e-3),
    "nmda_tau_PV":             dict(enabled=True, abs=(2e-3, 11e-3), default=5e-3),
    "nmda_tau_SST":            dict(enabled=False, abs=(5e-3, 40e-3), default=15e-3),
    "nmda_tau_VIP":            dict(enabled=False, abs=(5e-3, 40e-3), default=10e-3),

    # INPUTS 
    "Ib_strength":             dict(enabled=True, abs=(5, 15), default=10),
    "Iext_strength":           dict(enabled=True, abs=(1, 40), default=35),
    "Iext_duration":           dict(enabled=False, abs=(0.005, 0.02), default=0.013),
    "Im_strength":             dict(enabled=False, abs=(0, 10), default=0),
    "Ib_noise_std":            dict(enabled=False, abs=(0.5, 2), default=0.5),
    "Ib_ratio_E":              dict(enabled=False, abs=(1, 1), default=1),
    "Ib_ratio_PV":             dict(enabled=False, abs=(0.6, 0.9), default=0.7),
    "Ib_ratio_SST":            dict(enabled=False, abs=(0.7, 1.2), default=1),
    "Ib_ratio_VIP":            dict(enabled=False, abs=(0.4, 0.7), default=0.4),
    "scaling_factor":          dict(enabled=False, abs=(1, 1), default=1),

    # DIPOLE LENGTHS
    "F_L5": dict(enabled=False, abs=(dp.F_L5_MIN, 20.0), default=13),
    "F_L6": dict(enabled=False, abs=(1.0, 10.0), default=9),
    "L_BASE":        dict(enabled=False, abs=(0.5, 15.0), default=2),
    "F_L4_A1":       dict(enabled=False, abs=(0.01, 0.5),  default=dp.F_L4_A1),
    "F_L4_S2":       dict(enabled=False, abs=(0.01, 0.5),  default=dp.F_L4_S2),
    "L_THALE_CTX":   dict(enabled=False, abs=(0.5, 10.0), default=dp.L_THALE_CTX),
    "L_THALPOM_CTX": dict(enabled=False, abs=(0.5, 20.0), default=dp.L_THALPOM_CTX),
    "L_THALE_A3B":   dict(enabled=False, abs=(0.02, 1.0),  default=dp.L_THALE_A3B),
    "L_THALPOM_A3B": dict(enabled=False, abs=(0.5, 20.0), default=4),
    "L_A3B_CORT":    dict(enabled=False, abs=(0.5, 20.0), default=dp.L_A3B_CORT),
}

# the 10 parameters of the dipole parameterization are not model parameters: they go into
# the dipole lengths (model.dipole_length_params), see objective()
DIPOLE_PARAMS = tuple(dp.PARAM_NAMES)
_missing = [n for n in DIPOLE_PARAMS if n not in PARAMETER_DICT]
if _missing:
    raise ValueError(f"dipole parameters missing in PARAMETER_DICT: {_missing}")

# every entry needs bounds and a default; fail up front instead of inside the GA
for _n, _s in PARAMETER_DICT.items():
    if "abs" not in _s or "default" not in _s:
        raise ValueError(f"PARAMETER_DICT[{_n!r}] needs 'abs' and 'default'")

# genes = enabled entries only; these are the parameters the GA searches
fitted_params = [n for n, s in PARAMETER_DICT.items() if s["enabled"]]
bounds = [PARAMETER_DICT[n]["abs"] for n in fitted_params]
for _n in fitted_params:
    _lo, _hi = PARAMETER_DICT[_n]["abs"]
    if not _lo < _hi:
        raise ValueError(f"{_n}: enabled but its bounds are degenerate {(_lo, _hi)}")

# everything not searched is held at its default in every objective call
fixed_params = {n: s["default"] for n, s in PARAMETER_DICT.items() if not s["enabled"]}


# ── evaluation log ──────────────────
# main() sets the path of this run's evaluations.csv; every objective call appends one row.
# A crash inside the simulation is not caught here: it stops the run instead of being logged as a result.
EVAL_LOG_PATH = None
_n_eval = 0


def log_evaluation(params, dipole_params, err_tc, seconds):
    """Append one row (eval number, searched values, the 10 effective dipole parameters, error) to evaluations.csv."""
    global _n_eval
    _n_eval += 1
    if EVAL_LOG_PATH is None:
        return
    row = {"eval": _n_eval, **params, **dipole_params, "err_tc": err_tc, "seconds": round(seconds, 2)}
    pd.DataFrame([row]).to_csv(EVAL_LOG_PATH, mode="a", header=not os.path.exists(EVAL_LOG_PATH), index=False)


# ── Objective Function ──────────────
def objective(**params):
    """
    Run one full simulation and return the error. 
    """

    # de-normalize parameters
    params = denormalize_params(params, bounds)


    t0 = time.time()
    # searched values override the fixed ones; the dipole genes are not model parameters
    all_params = {**fixed_params, **params}
    dipole_params = {k: all_params.pop(k) for k in DIPOLE_PARAMS}
    # scaling_factor belongs to the error computation, not to the model
    scaling_factor = all_params.pop("scaling_factor")
    model.apply_params(all_params)
    model.dipole_length_params = dipole_params

    model.initialize_state()
    model.simulate()

    mid_rate_error = 0
    if MID_FIRING_RATE:
        # evaluate how far the firing rates are from their midpoint during baseline
        stimulation_start = model.input_onset / model.step_size 
        prestim_rates = model.rates[:,:stimulation_start]
        firing_rate_bounds = model.p.get_sigmoids
        # TODO: implement error increaseth
        




    sim_dip = model.compute_dipoles()
    err_tc = 0.0

    err_tc, _, _ = model.compute_error_timecourse(tc_data_path, sim_dip, scaling_factor=scaling_factor, target_dip=None, rois=FIT_ROIS)
    error = err_tc + mid_rate_error

    log_evaluation(params, dipole_params, error, time.time() - t0)
    return error


# ── optimization settings ───────────────────────────────────────────────────────────────────
opt_config = {
    "model_parameters": fitted_params,
    "bounds": np.array(bounds, dtype=float),
    "reference":  0.0,
    "simulation": objective,
    "op":         -1,    # minimise
    "N1":         60,    # initial population size
    "N2":         80,        # crossover offspring per iteration
    "N3":         80,       # mutation offspring per iteration
    "n_iter":     20,
    # Gradient-search tolerance (conf['gTol']): gauss_newton_slow quits once the cost improves
    # by less than this between inner loops. The cost is err**2, and a scored "ps" err is
    # ~0.013..0.049, i.e. a cost of ~2e-4..2e-3 — far below the old 0.05, so every gradient
    # search quit after two loops having spent ~80 simulations it could never use.
    "tolerance":  1e-8,   # gradient-search tolerance (conf['gTol'])
    # Early-stop threshold on the GA cost, which is err**2. Kept far below any error the
    # objectives actually reach so the GA always runs its full n_iter: the default 1e-5 sits
    # above a converged fit's cost and would break after the first iteration, returning a
    # random draw from the initial population.
    "single_run_tol": 1e-12,
    "verbose":    1,
}

# ── base simulation settings ───────────────────────────────────────────────────────────────────
base_params = {
    "step_size": 0.001,
    "resolution_tstep": 0.001,
    "input_type": "step",
    "input_onset": 4.001,
    "simulation_dur": 6,
    "thal_connect": [
        0,
        1000,
        300,
        0,
        300,
        1000
    ],
    "extI_cellcounts": 1000,
    "mI_cellcounts": 2000,
    "bI_cellcounts": 100,
    "thalE_cellcounts": 500,
    "thalI_cellcounts": 500,
    "pom_cellcounts": 500,
    "pomI_cellcounts": 500,
    "resistance_factor": 1,
}

# initialize the model
model = sm.SomatoModel(base_params)
# dipole lengths from the 10-parameter parameterization (the 10 dipole parameters are in PARAMETER_DICT) instead of the JSON file
model.load_dipole_params = model.load_dipole_params_parameterized
# every searched or fixed model parameter must exist on the model, otherwise apply_params
# would only create an unused attribute (typo in PARAMETER_DICT)
_unknown = [n for n in PARAMETER_DICT if n not in DIPOLE_PARAMS + ("scaling_factor",) and not hasattr(model, n)]
if _unknown:
    raise ValueError(f"PARAMETER_DICT entries that are not model parameters: {_unknown}")

# ── start from an earlier run ──────────────────────────────────────────────────
# Session folder name (or path) under diag_dir, e.g. "optfull_20261110_013709_roi-all"; empty = random start.
# The environment variable START_FROM overrides it. The old run must have the same full parameter set
# (the keys of PARAMETER_DICT); which of them it optimized and which it held fixed may differ.
START_FROM = "optfull_20261110_174017_roi-all_UNSATURATED" # "optfull_20261110_123724_roi-all"
START_FROM_OPT_RUN = START_FROM #os.getenv("START_FROM", "")


def load_start_params(run, names, bounds):
    """Seed {parameter: value} for the searched `names` from the finished or interrupted run `run`.

    Hard errors (SystemExit): missing folder or files, an old parameter set that differs from
    PARAMETER_DICT, a seed value outside `bounds`. A parameter the old run searched takes its best
    value (optimization_summary.json, else the best row of evaluations.csv); one it held fixed takes
    its fixed value."""
    run_dir = run if os.path.isdir(run) else os.path.join(diag_dir, run)
    cfg_path = os.path.join(run_dir, "run_config.json")
    if not os.path.isfile(cfg_path):
        raise SystemExit(f"START_FROM: {cfg_path} not found")
    with open(cfg_path) as f:
        old_cfg = json.load(f)

    old_all = set(old_cfg["optimized_parameters"]) | set(old_cfg["fixed_params"])
    missing, extra = sorted(set(PARAMETER_DICT) - old_all), sorted(old_all - set(PARAMETER_DICT))
    if missing or extra:
        raise SystemExit(f"START_FROM={run}: its parameter set differs from PARAMETER_DICT "
                         f"(not in the old run: {missing}; not in PARAMETER_DICT: {extra})")

    summary_path = os.path.join(run_dir, "optimization_summary.json")
    if os.path.isfile(summary_path):
        with open(summary_path) as f:
            old_best = json.load(f)["best_params"]
        source = "optimization_summary.json"
    else:
        eval_path = os.path.join(run_dir, "evaluations.csv")
        if not os.path.isfile(eval_path):
            raise SystemExit(f"START_FROM={run}: neither optimization_summary.json nor evaluations.csv found")
        ev = pd.read_csv(eval_path)
        best_row = ev.loc[ev["err_tc"].idxmin()]
        old_best = {n: float(best_row[n]) for n in old_cfg["optimized_parameters"]}
        source = f"best of {len(ev)} evaluations in evaluations.csv (no summary: run unfinished), error {best_row['err_tc']:.6f}"

    seed = {n: float(old_best[n] if n in old_best else old_cfg["fixed_params"][n]) for n in names}

    outside = [f"{n}={v:.6g} not in [{lo:g}, {hi:g}]" for (n, v), (lo, hi) in zip(seed.items(), bounds)
               if not lo <= v <= hi]
    if outside:
        raise SystemExit(f"START_FROM={run}: seed outside the current bounds: " + "; ".join(outside))

    print(f"Starting from {run}: {source}")
    for n, v in seed.items():
        print(f"  {n} = {v:.6g}" + ("" if n in old_best else "  (fixed in the old run)"))
    return seed


initial_params = None
if START_FROM_OPT_RUN:
    _seed = load_start_params(START_FROM_OPT_RUN, fitted_params, bounds)
    opt_config["initial_population"] = np.array([[_seed[n] for n in fitted_params]], dtype=float)
    initial_params = {"source": START_FROM_OPT_RUN, "values": _seed}

# FUNCTIONS


def plot_best_fit(model, best_params, outdir):
    """Re-simulate the best parameters like objective() and save the fit against the measured time courses.

    Writes best_tc_comparison.png (simulated vs measured per ROI, normalised as in the error) and
    best_tc_comparison.hdf5 (the traces and the error) to outdir. Returns the time-course error."""
    all_params = {**fixed_params, **best_params}
    dipole_params = {k: all_params.pop(k) for k in DIPOLE_PARAMS}
    scaling_factor = all_params.pop("scaling_factor")
    model.apply_params(all_params)
    model.dipole_length_params = dipole_params
    model.initialize_state()
    model.simulate()
    sim_dip = model.compute_dipoles()

    # scored on FIT_ROIS as in the GA; the figure shows all three ROIs
    err_tc, tc_sim, tc_target = model.compute_error_timecourse(
        tc_data_path, sim_dip, scaling_factor=scaling_factor, target_dip=None, rois=FIT_ROIS)
    model.plot_timecourse_comparison(tc_sim, tc_target, global_peak_rois=FIT_ROIS,
                                     save_path=os.path.join(outdir, "best_tc_comparison.png"))
    model.save_timecourse_comparison(outdir, tc_sim, tc_target, err_tc, filename="best_tc_comparison")
    return err_tc


def normalize_params(params, bounds):
    """Physical parameter values -> [0, 1] per parameter (column).

    `params` is one vector (n_params,) or a population (n_pop, n_params); `bounds` is (n_params, 2).
    Values outside the bounds are clipped."""
    bounds = np.asarray(bounds, dtype=float)
    lo, hi = bounds[:, 0], bounds[:, 1]
    return np.clip((np.asarray(params, dtype=float) - lo) / (hi - lo), 0.0, 1.0)


def denormalize_params(params, bounds):
    """Normalised values of the GA -> physical values.

    `params` is a dict {name: value in [0, 1]} in the order of fitted_params (the order the GA passes
    the keyword arguments in), which is the order of `bounds`; returns a new dict of physical values."""
    if len(params) != len(bounds):
        raise ValueError(f"{len(params)} parameters but {len(bounds)} bounds")
    return {name: v * (hi - lo) + lo for (name, v), (lo, hi) in zip(params.items(), bounds)}


def denormalize_vector(vec, bounds):
    """Same as denormalize_params for a vector (n_params,) in the order of `bounds`."""
    bounds = np.asarray(bounds, dtype=float)
    return np.asarray(vec, dtype=float) * (bounds[:, 1] - bounds[:, 0]) + bounds[:, 0]


# ── run ────────────────────────────────────────────────────────────────────────
def main():
    # Each optimization run gets its own timestamped session folder so runs no
    # longer overwrite each other. All figures/CSVs for this run go here.
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    session_dir = os.path.join(diag_dir, f"optfull_{run_id}_roi-{roi_tag}")
    os.makedirs(session_dir, exist_ok=True)
    global EVAL_LOG_PATH
    EVAL_LOG_PATH = os.path.join(session_dir, "evaluations.csv")
    print(f"Optimization session: {session_dir}")
    print(f"Fitting ROIs: {', '.join(FIT_ROIS)}  (ROI={_roi_env!r})")

    # Log the full run configuration up front (persisted even if the run crashes).
    # target_csvs are recorded before preprocess_targets repoints them below, so
    # the log keeps the raw measured sources.
    run_config = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "fit_rois": list(FIT_ROIS),
        "optimized_parameters": opt_config["model_parameters"],
        "bounds": np.asarray(opt_config["bounds"]).tolist(),
        "fixed_params": fixed_params,
        "parameter_defaults": {n: sp["default"] for n, sp in PARAMETER_DICT.items()},
        "base_params": base_params,   # simulation settings the model was built with (needed to re-simulate the run)
        "ga": {k: opt_config[k] for k in ("N1", "N2", "N3", "n_iter", "tolerance", "op", "reference")},
        "target_csv": tc_data_path,
        "initial_params": initial_params,
        "search_space": "normalized",   # the GA works on [0, 1]; evaluations.csv and the summary are in physical units
    }
    with open(os.path.join(session_dir, "run_config.json"), "w") as f:
        json.dump(run_config, f, indent=2)


    # the GA searches normalised values (0 to 1 for every parameter); objective() converts them back
    # to physical values with the bounds in `bounds`
    opt_config["bounds"] = np.repeat([[0.0, 1.0]], len(fitted_params), axis=0)
    if "initial_population" in opt_config:      # seed from START_FROM (physical values) -> normalised
        opt_config["initial_population"] = normalize_params(opt_config["initial_population"], bounds)

    # run GA
    ga = GA(opt_config)
    ga.run()

    # ga.optimum is normalised: convert to physical values for the model, the dipole file and the summary
    best_params = dict(zip(opt_config["model_parameters"], denormalize_vector(ga.optimum, bounds)))

    # dipole file of the best parameters: the 10 dipole parameters (searched ones from the optimum,
    # the others at their fixed value), in the schema the model reads (dipole_params_file)
    best_dipole = {k: float(best_params.get(k, fixed_params.get(k))) for k in DIPOLE_PARAMS}
    dp.build_params_json(os.path.join(session_dir, "dipole_parameters_finetuned.json"), best_dipole)

    print("\n── Optimised parameters ──")
    err_tc = plot_best_fit(model, best_params, session_dir)
    print(f"  Best error: {ga.errors[-1]}")

    # diagnostics of the optimisation itself
    ga.plot_fit()


    # write a summary of the run
    summary = {
        "fit_rois": list(FIT_ROIS),
        "best_params": {k: float(v) for k, v in best_params.items()},
        "dipole_parameters": best_dipole,   # the 10 dipole parameters of the best fit (as in the ft_* summaries)
        "best_combined_error": float(np.sqrt(ga.errors[-1])),
        "cost_per_iteration": [float(e) for e in ga.errors],   # GA cost = err**2
        "error_per_iteration": [float(np.sqrt(e)) for e in ga.errors],
        "best_fit_err_tc": float(err_tc),
    }

    with open(os.path.join(session_dir, "optimization_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  Diagnostics written to {session_dir}")

if __name__ == "__main__":
    main()
