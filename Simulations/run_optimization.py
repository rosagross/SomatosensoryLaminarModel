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

# ── paths ──────────────────────────────────────────────────────────────────────
WDDIR  = os.getenv("WDDIR")   # /data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel
RESDIR = os.getenv("RESDIR")  # /data/p_02989/shared_workspace/results_grossmannr

sys.path.append(os.path.join(WDDIR, "Simulations", "model"))
sys.path.insert(0, "/data/p_02989/Modelling/neuronaldynamics/src")
sys.path.insert(0, "/data/p_02989/Modelling/neuronaldynamics")  # for Utils.py

from somato_model import SomatoModel
from neuronaldynamics.Optimizers.Optimizer import GA

# ── target data paths ──────────────────────────────────────────────────────────
_roi_dir = os.path.join(
    RESDIR, "Figures", "Main", "eeg_results", "source_reconstruction",
    "group", "_preprestim_corrected", "roi_epochswise",
)
tf_data_path = os.path.join(_roi_dir, "group_roi_tf_morlet_ses-elec_preprestim_corrected.csv")
tc_data_path = os.path.join(_roi_dir, "group_roi_timecourse_pooled_ses-elec_preprestim_corrected.csv")

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
def objective(**params):
    """
    Run one full simulation and return the combined TF + timecourse error.
    The GA minimises (0 - objective)**2, i.e. the squared combined error.
    """
    model.apply_params(params)
    model.initialize_state()
    model.simulate()
    err_tf, _, _ = model.compute_error_timefreq(tf_data_path)
    err_tc, _, _ = model.compute_error_timecourse(tc_data_path)
    combined = err_tf + err_tc
    print(f"  params={params}  →  err_tf={err_tf:.4f}  err_tc={err_tc:.4f}  total={combined:.4f}")
    return combined

# ── GA setup ───────────────────────────────────────────────────────────────────
opt_config = {
    "model_parameters": [
        "coupling_strength",
        "strength_I",
        "g_intercortical",
        "Ib_strength",
        "Iext_strength",
        "Iext_duration",
    ],
    "bounds": np.array([
        [0,     50  ],   # coupling_strength
        [0,      0.5],   # strength_I
        [0,      2  ],   # g_intercortical
        [0,     10  ],   # Ib_strength
        [0,    100  ],   # Iext_strength
        [0.001,  0.1],   # Iext_duration
    ]),
    "reference":  0.0,
    "simulation": objective,
    "op":         -1,    # minimise
    "N1":         10,    # initial population size
    "N2":         20,    # crossover offspring per iteration
    "N3":         20,    # mutation offspring per iteration
    "n_iter":     5,
    "tolerance":  0.05,
    "verbose":    1,
}

# ── run ────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    ga = GA(opt_config)
    ga.run()

    best_params = dict(zip(opt_config["model_parameters"], ga.optimum))
    print("\n── Optimised parameters ──")
    for name, val in best_params.items():
        print(f"  {name}: {val:.4f}")
    print(f"  Best combined error: {ga.errors[-1]:.4f}")

    ga.plot_fit()
