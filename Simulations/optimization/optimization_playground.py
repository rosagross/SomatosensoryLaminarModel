"""
File: optimization_playground.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-11-10
Description: in this script we can load an optimization run and then change parameters to see what happens.
The parameters that are changed are overwritten.
It plots directly in the interactive session, no plots are saved.
"""

# %%
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# the helper functions of the diagnostics script (also puts the model folders on sys.path)
import visualize_optmization_results as viz
import dipole_parameterization as dp
import somato_model as sm
from plot_dipole_computation import plot_layer_interneuron_contributions

# Run to check out
OPT_RUN_ID = "optfull_20261110_162957_roi-all" 

# All parameters
PARAMETER_DICT = {
    # DELAYS
    "thal_delay_factor":       dict(overwrite=False, value=0.015),
    "delay_factor":            dict(overwrite=False, value=0.011),
    "delay_factor_short":      dict(overwrite=False, value=0.004),
    "receptor_thalamus_delay": dict(overwrite=False, value=0.01),
    
    # CONNECTIVITY & COUPLING 
    "coupling_strength":       dict(overwrite=True, value =1),
    "strength_I":              dict(overwrite=True, value =1),
    "g_intercortical":         dict(overwrite=True, value=1),
    # intra-cortical
    "f_intra_S2":              dict(overwrite=True, value=1),
    "f_intra":                 dict(overwrite=True, value=1),
    "p_2PVE":                  dict(overwrite=False, value=37.6),
    "p_4PVE":                  dict(overwrite=False, value=39.4),
    # feedforward
    "ff_A1_S2":                dict(overwrite=False, value=1),
    "ff_A3b_S2":               dict(overwrite=False, value=1),
    "ff_A3b_A1":               dict(overwrite=False, value=1),
    # feedback
    "fb_A1_S2":                dict(overwrite=False, value=4),
    "fb_A3b_A1":               dict(overwrite=False, value=4),
    "fb_A3b_S2":               dict(overwrite=False, value=4),
    # thalamus
    "g_thalPOm":               dict(overwrite=False, value=1),
    "g_thal":                  dict(overwrite=False, value=2),
    "POm_drive_A3b":           dict(overwrite=False, value=1),
    "POm_drive_A1":            dict(overwrite=False, value=1),
    "POm_drive_S2":            dict(overwrite=False, value=1),
    "POm_tgt_A3b":             dict(overwrite=False, value=1),
    "POm_tgt_A1":              dict(overwrite=False, value=1),
    "POm_tgt_S2":              dict(overwrite=False, value=1),
    
    # TIME CONSTANTS
    "e3b_tau":                 dict(overwrite=False, value=6),
    "e1_tau":                  dict(overwrite=False, value=6),
    "e2_tau":                  dict(overwrite=False, value=6),
    # NMDA decay (s) of E -> cell class connections, by postsynaptic class
    "nmda_tau_E":              dict(overwrite=False, value=10e-3),
    "nmda_tau_PV":             dict(overwrite=False, value=5e-3),
    "nmda_tau_SST":            dict(overwrite=False, value=25e-3),
    "nmda_tau_VIP":            dict(overwrite=False, value=20e-3),

    # INPUTS 
    "Ib_strength":             dict(overwrite=True, value=15),
    "Iext_strength":           dict(overwrite=False, value=35),
    "Iext_duration":           dict(overwrite=False, value=0.013),
    "Im_strength":             dict(overwrite=False, value=0),
    "Ib_noise_std":            dict(overwrite=False, value=0.5),
    "Ib_ratio_E":              dict(overwrite=False, value=1),
    "Ib_ratio_PV":             dict(overwrite=False, value=0.7),
    "Ib_ratio_SST":            dict(overwrite=False, value=1),
    "Ib_ratio_VIP":            dict(overwrite=False, value=0.4),
    "scaling_factor":          dict(overwrite=False, value=1),

    # DIPOLE LENGTHS
    "F_L5": dict(overwrite=True, value=5),
    "F_L6": dict(overwrite=True, value=5),
    "L_BASE":        dict(overwrite=True, value=8),
    "F_L4_A1":       dict(overwrite=False, value=0.05),
    "F_L4_S2":       dict(overwrite=False, value=0.05),
    "L_THALE_CTX":   dict(overwrite=False, value=1),
    "L_THALPOM_CTX": dict(overwrite=False, value=2),
    "L_THALE_A3B":   dict(overwrite=False, value=0.1),
    "L_THALPOM_A3B": dict(overwrite=False, value=2),
    "L_A3B_CORT":    dict(overwrite=True, value=5),
}

# SIMULATE WITH OPTIMIZED (or overwritten parameters)
# Same simulation as objective() in run_optimization_full.py: the best parameters of the run
# (optimization_summary.json, else the best row of evaluations.csv), the run's fixed parameters,
# and on top of that every PARAMETER_DICT entry with overwrite=True.
run = viz.load_run(viz.resolve_run_dir(OPT_RUN_ID))
cfg = run["config"]
best, best_logged_err, source = viz.best_parameters(run)

unknown = [n for n, s in PARAMETER_DICT.items() if s["overwrite"] and n not in cfg["parameter_defaults"]]
if unknown:
    raise ValueError(f"overwrite requested for parameters the run does not know: {unknown}")

all_params = {**cfg["fixed_params"], **best}
overwritten = {n: s["value"] for n, s in PARAMETER_DICT.items() if s["overwrite"]}
all_params.update(overwritten)
print(f"run {OPT_RUN_ID}: parameters from {source}, best logged error {best_logged_err:.6f}")
print("overwritten:", overwritten if overwritten else "none")

dipole_params = {k: all_params.pop(k) for k in dp.PARAM_NAMES}
scaling_factor = all_params.pop("scaling_factor")      # belongs to the error, not to the model

model = sm.SomatoModel(cfg["base_params"])
model.load_dipole_params = model.load_dipole_params_parameterized
model.apply_params(all_params)
model.dipole_length_params = dipole_params
model.initialize_state()
model.simulate()
sim_dip = model.compute_dipoles()
err, tc_sim, tc_target = model.compute_error_timecourse(
    cfg["target_csv"], sim_dip, scaling_factor=scaling_factor, target_dip=None, rois=tuple(cfg["fit_rois"]))
print(f"time-course error: {err:.6f}" + ("" if overwritten else " (no overwrite: should equal the best logged error)"))



plot1 = True
plot2 = False
plot3 = False
plot4 = True


if plot1:
    # PLOT 1 - simulated time-course error
    # simulated vs measured ROI time courses, normalised as in the error; per-ROI error terms below
    model.plot_timecourse_comparison(tc_sim, tc_target, show=True, global_peak_rois=tuple(cfg["fit_rois"]))
    fig = viz.plot_roi_errors(run, err, tc_sim, tc_target)
    plt.show()

if plot2:
    # Plot 2 - dipole computation contribution (full interneuron contributions)
    # which source populations (E/PV/SST/VIP) build the dipole of each layer
    for area in ("A3b", "A1", "S2"):
        plot_layer_interneuron_contributions(model, sim_dip, area)
        plt.show()

if plot3:
    # Plot 3 - dipole length plot (for each layer and area)
    # dipole length of each input cell type onto the E population of each layer
    viz.plot_dipole_length_by_celltype(model)
    plt.show()

if plot4:
    # Plot 4 - saturation of the population firing rates
    # baseline rate and post-stimulus range of every population relative to its maximum rate,
    # then the sigmoid of every population with the drive it visits from -500 to +250 ms
    viz.plot_saturation(model)
    plt.show()
    viz.plot_sigmoid_operating_points(model)
    plt.show()


