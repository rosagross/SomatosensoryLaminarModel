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
    "thal_delay_factor":       dict(enabled=True, abs=(0.003, 0.01), default=0.005),
    "delay_factor":            dict(enabled=True, abs=(0.005, 0.015), default=0.01),
    "delay_factor_short":      dict(enabled=True, abs=(0.001, 0.006), default=0.002),
    "receptor_thalamus_delay": dict(enabled=False, abs=(0.005, 0.025), default=0.01),
    
    # CONNECTIVITY & COUPLING 
    "coupling_strength":       dict(enabled=False, abs=(0.5, 10), default=1), 
    "strength_I":              dict(enabled=False, abs=(0.5, 1.5), default=1), 
    "g_intercortical":         dict(enabled=False, abs=(0.5, 5), default=1),
    # intra-cortical
    "f_intra_S2":              dict(enabled=False, abs=(0.5, 3.0), default=1),
    "f_intra":                 dict(enabled=True, abs=(0.5, 3.0), default=1),
    "p_2PVE":                  dict(enabled=False, abs=(37.4, 37.4), default=37.4),
    "p_4PVE":                  dict(enabled=False, abs=(39.6, 39.6), default=39.6),
    # feedforward
    "ff_A1_S2":                dict(enabled=True, abs=(0.5, 2.0), default=1),
    "ff_A3b_S2":               dict(enabled=False, abs=(0.5, 2.0), default=1),
    "ff_A3b_A1":               dict(enabled=True, abs=(0.5, 2.0), default=1),
    # feedback
    "fb_A1_S2":                dict(enabled=True, abs=(0.5, 2.0), default=1),
    "fb_A3b_A1":               dict(enabled=True, abs=(0.5, 2.0), default=1),
    "fb_A3b_S2":               dict(enabled=True, abs=(0.5, 2.0), default=1),
    # thalamus
    "g_thalPOm":               dict(enabled=False, abs=(0.5, 3.0), default=1),
    "g_thal":                  dict(enabled=False, abs=(0.5, 3.0), default=1),
    "POm_drive_A3b":           dict(enabled=False, abs=(0.5, 2.0), default=1),
    "POm_drive_A1":            dict(enabled=True, abs=(0.5, 2.0), default=1),
    "POm_drive_S2":            dict(enabled=True, abs=(0.5, 2.0), default=1),
    "POm_tgt_A3b":             dict(enabled=False, abs=(0.5, 2.0), default=1),
    "POm_tgt_A1":              dict(enabled=False, abs=(0.5, 2.0), default=1),
    "POm_tgt_S2":              dict(enabled=False, abs=(0.5, 2.0), default=1),
    
    # TIME CONSTANTS
    "e3b_tau":                 dict(enabled=True, abs=(4.0, 8.0), default=6),
    "e1_tau":                  dict(enabled=True, abs=(4.0, 8.0), default=6),
    "e2_tau":                  dict(enabled=True, abs=(4.0, 8.0), default=6),

    # INPUTS 
    "Ib_strength":             dict(enabled=False, abs=(1, 60), default=5),
    "Iext_strength":           dict(enabled=False, abs=(1, 100), default=10),
    "Iext_duration":           dict(enabled=False, abs=(0.005, 0.02), default=0.01),
    "Im_strength":             dict(enabled=False, abs=(0, 10), default=0),
    "Ib_noise_std":            dict(enabled=False, abs=(0.5, 2), default=0.5),
    "Ib_ratio_E":              dict(enabled=False, abs=(1, 1), default=1),
    "Ib_ratio_PV":             dict(enabled=False, abs=(0.6, 0.9), default=0.7),
    "Ib_ratio_SST":            dict(enabled=False, abs=(0.7, 1.2), default=1),
    "Ib_ratio_VIP":            dict(enabled=False, abs=(0.4, 0.7), default=0.6),
    "scaling_factor":          dict(enabled=False, abs=(1, 1), default=1),

    # DIPOLE LENGTHS
    "F_L5": dict(enabled=True, abs=(dp.F_L5_MIN, 10.0), default=5),
    "F_L6": dict(enabled=True, abs=(1.0, 10.0), default=5),
}

fitted_params = [p for p in PARAMETER_DICT]
parameter_vals = [PARAMETER_DICT[p] for p in fitted_params]
bounds = [p['abs'] for p in parameter_vals]


# ── Objective Function ──────────────
def objective(**params):
    """
    Run one full simulation and return the error. 
    """

    model.apply_params(params)

    model.initialize_state()
    model.simulate()
    sim_dip = model.compute_dipoles()
    err_tf = err_tc = err_ps = 0.0

    err_tc, _, _ = model.compute_error_timecourse(tc_data_path, sim_dip, target_dip=None, rois=FIT_ROIS)

    return err_tc


# ── optimization settings ───────────────────────────────────────────────────────────────────
N1, N2, N3, N_ITER = 10, 20, 20, 1     # GA: initial population, crossover / mutation offspring, iterations
opt_config = {
    "model_parameters": fitted_params,
    "bounds": np.array(bounds, dtype=float),
    "reference":  0.0,
    "simulation": objective,
    "op":         -1,    # minimise
    "N1":         30,    # initial population size
    "N2":         80,        # crossover offspring per iteration
    "N3":         80,       # mutation offspring per iteration
    "n_iter":     1,
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
    # the following parameters are only place holders for the rest of the required ones, 
    # we'll use the defaults defined above 
    
     
}

# initialize the model
model = sm.SomatoModel(base_params)


# FUNCTIONS


#def plot_best_fit(model, best_params, outdir):
#    #TODO: to be implemented


# ── run ────────────────────────────────────────────────────────────────────────
def main():
    # Each optimization run gets its own timestamped session folder so runs no
    # longer overwrite each other. All figures/CSVs for this run go here.
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    session_dir = os.path.join(diag_dir, f"optfull_{run_id}_roi-{roi_tag}")
    os.makedirs(session_dir, exist_ok=True)
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
        "ga": {k: opt_config[k] for k in ("N1", "N2", "N3", "n_iter", "tolerance", "op", "reference")},
        "target_csv": tc_data_path,
        #"initial_params": initial_params, # TODO: only those with a hardcoded "init" value
    }
    with open(os.path.join(session_dir, "run_config.json"), "w") as f:
        json.dump(run_config, f, indent=2)

    # run GA
    ga = GA(opt_config)
    ga.run()

    best_params = dict(zip(opt_config["model_parameters"], ga.optimum))
    print("\n── Optimised parameters ──")
    err_tf, err_tc, err_ps = plot_best_fit(model, best_params, session_dir)
    print(f"  Best error: {ga.errors[-1]}")

    # diagnostics of the optimisation itself
    ga.plot_fit()

    # TODO: plot fitting diagnostics

    # write a summary of the run
    summary = {
        "fit_rois": list(FIT_ROIS),
        "best_params": {k: float(v) for k, v in best_params.items()},
        "best_combined_error": float(np.sqrt(ga.errors[-1])),
        "cost_per_iteration": [float(e) for e in ga.errors],   # GA cost = err**2
        "error_per_iteration": [float(np.sqrt(e)) for e in ga.errors],
        "best_fit_err_tf": float(err_tf),
        "best_fit_err_tc": float(err_tc),
        "best_fit_err_ps": float(err_ps),
    }

    with open(os.path.join(session_dir, "optimization_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  Diagnostics written to {session_dir}")

if __name__ == "__main__":
    main()
