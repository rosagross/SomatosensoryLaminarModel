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
import pandas as pd
import matplotlib.pyplot as plt

# the helper functions of the diagnostics script (also puts the model folders on sys.path)
import visualize_optmization_results as viz
import dipole_parameterization as dp
import somato_model as sm
from plot_dipole_computation import plot_layer_interneuron_contributions

# Run to check out
OPT_RUN_ID = "optfull_20261110_174017_roi-all_UNSATURATED" 

# All parameters
PARAMETER_DICT = {
    # DELAYS
    "thal_delay_factor":       dict(overwrite=False, value=0.008),
    "delay_factor":            dict(overwrite=False, value=0.011),
    "delay_factor_short":      dict(overwrite=False, value=0.004),
    "receptor_thalamus_delay": dict(overwrite=False, value=0.01),
    
    # CONNECTIVITY & COUPLING 
    "coupling_strength":       dict(overwrite=False, value =1),
    "strength_I":              dict(overwrite=False, value =1),
    "g_intercortical":         dict(overwrite=False, value =1),
    # intra-cortical
    "f_intra_S2":              dict(overwrite=False, value=1),
    "f_intra":                 dict(overwrite=False, value=1),
    "p_2PVE":                  dict(overwrite=False, value=37.6),
    "p_4PVE":                  dict(overwrite=False, value=39.4),
    # feedforward
    "ff_A1_S2":                dict(overwrite=False, value=1),
    "ff_A3b_S2":               dict(overwrite=False, value=1),
    "ff_A3b_A1":               dict(overwrite=False, value=2.5),
    # feedback
    "fb_A1_S2":                dict(overwrite=False, value=4),
    "fb_A3b_A1":               dict(overwrite=False, value=7),
    "fb_A3b_S2":               dict(overwrite=False, value=4),
    # thalamus
    "g_thalPOm":               dict(overwrite=False, value=1),
    "g_thal":                  dict(overwrite=False, value=2),
    "POm_drive_A3b":           dict(overwrite=False, value=1),
    "POm_drive_A1":            dict(overwrite=False, value=1),
    "POm_drive_S2":            dict(overwrite=False, value=1),
    "POm_tgt_A3b":             dict(overwrite=False, value=1),
    "POm_tgt_A1":              dict(overwrite=False, value=5),
    "POm_tgt_S2":              dict(overwrite=False, value=1.5),
    
    # TIME CONSTANTS
    "e3b_tau":                 dict(overwrite=False, value=4),
    "e1_tau":                  dict(overwrite=False, value=6),
    "e2_tau":                  dict(overwrite=False, value=6),
    # NMDA decay (s) of E -> cell class connections, by postsynaptic class
    "nmda_tau_E":              dict(overwrite=False, value=10e-3),
    "nmda_tau_PV":             dict(overwrite=False, value=5e-3),
    "nmda_tau_SST":            dict(overwrite=False, value=25e-3),
    "nmda_tau_VIP":            dict(overwrite=False, value=20e-3),

    # INPUTS 
    "Ib_strength":             dict(overwrite=False, value=14),
    "Iext_strength":           dict(overwrite=False, value=10),
    "Iext_duration":           dict(overwrite=False, value=0.013),
    "Im_strength":             dict(overwrite=False, value=0),
    "Ib_noise_std":            dict(overwrite=False, value=0.5),
    "Ib_ratio_E":              dict(overwrite=False, value=1),
    "Ib_ratio_PV":             dict(overwrite=False, value=0.7),
    "Ib_ratio_SST":            dict(overwrite=False, value=1),
    "Ib_ratio_VIP":            dict(overwrite=False, value=0.4),
    "scaling_factor":          dict(overwrite=False, value=1),

    # DIPOLE LENGTHS
    "F_L5": dict(overwrite=False, value=5),
    "F_L6": dict(overwrite=False, value=5),
    "L_BASE":        dict(overwrite=False, value=8),
    "F_L4_A1":       dict(overwrite=False, value=0.1),
    "F_L4_S2":       dict(overwrite=False, value=0.05),
    "L_THALE_CTX":   dict(overwrite=False, value=1),
    "L_THALPOM_CTX": dict(overwrite=False, value=2),
    "L_THALE_A3B":   dict(overwrite=False, value=2),
    "L_THALPOM_A3B": dict(overwrite=False, value=2),
    "L_A3B_CORT":    dict(overwrite=False, value=5),
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



# spectrum of the pre-stimulus activity (own implementation, independent of the model and the diagnostics script)
def plot_prestim_spectrum(model, sim_dip, fmin=1.0, fmax=40.0, seg_dur=1.0, overlap=0.5, settle_s=2.0):
    """Power spectrum of the pre-stimulus activity of the three ROI dipoles (best fit).

    ROI signals: A3b = sim_dip[0], A1 = sum of sim_dip[1:5], S2 = sum of sim_dip[5:9]. Welch estimate over the
    settled pre-stimulus period (first settle_s seconds skipped): segments of seg_dur seconds with the given
    overlap, the last one ending at stimulus onset; each segment has its mean removed and a Hann window, power is
    |rfft|^2 / n^2, averaged over the segments. Row 1: linear axes, row 2: log-log axes. The strongest peak
    between 5 and 15 Hz is marked."""
    dt = model.step_size
    stim_idx = round(model.input_onset / dt) - 1
    n = round(seg_dur / dt)
    hop = max(round(n * (1 - overlap)), 1)
    first = round(settle_s / dt)
    starts = []
    s0 = stim_idx - n
    while s0 >= first:
        starts.append(s0)
        s0 -= hop
    if not starts:
        raise ValueError(f"no {seg_dur} s segment fits between {settle_s} s and stimulus onset at {stim_idx * dt:.3f} s")
    starts.sort()

    sim_dip = np.asarray(sim_dip)
    rois = {"A3b": sim_dip[0], "A1": sim_dip[1:5].sum(axis=0), "S2": sim_dip[5:9].sum(axis=0)}
    win = np.hanning(n)
    freqs = np.fft.rfftfreq(n, d=dt)
    keep = (freqs >= fmin) & (freqs <= fmax)
    band = (freqs[keep] >= 5) & (freqs[keep] <= 15)

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for col, (roi, sig) in enumerate(rois.items()):
        pre = sig[starts[0]:stim_idx]
        if np.std(pre) <= 1e-9 * max(np.abs(np.mean(pre)), 1e-30):
            print(f"prestim spectrum: WARNING the pre-stimulus signal of {roi} is flat, its spectrum is numerical residue")
        psd = np.zeros(keep.sum())
        for s0 in starts:
            seg = sig[s0:s0 + n]
            psd += (np.abs(np.fft.rfft((seg - seg.mean()) * win)) ** 2 / n ** 2)[keep]
        psd /= len(starts)
        f = freqs[keep]
        i_pk = int(np.argmax(np.where(band, psd, -np.inf)))
        for row, ax in enumerate(axes[:, col]):
            ax.plot(f, psd, color="k", linewidth=1.4)
            ax.plot(f[i_pk], psd[i_pk], "o", color="C3", markeredgecolor="k", zorder=5)
            ax.set_xlabel("frequency (Hz)")
            if row == 1:
                ax.set_xscale("log")
                ax.set_yscale("log")
        axes[0, col].set_title(f"{roi}: peak {f[i_pk]:.1f} Hz (5-15 Hz)")
        axes[0, col].set_ylabel("power")
        axes[1, col].set_ylabel("power (log)")
    fig.suptitle(f"pre-stimulus spectrum of the ROI dipoles, best fit: {len(starts)} segments of {seg_dur:g} s, "
                 f"{overlap:.0%} overlap, ending at stimulus onset")
    return fig


def plot_2D_errorlandscape(run, param1, param2):
    """Error (err_tc) over two optimized parameters: every evaluation as a point coloured by error (no interpolation),
    the best fit as a star."""
    ev = run["evaluations"]
    best = ev.loc[ev["err_tc"].idxmin()]
    data = ev[[param1, param2, "err_tc"]].astype(float).replace([np.inf, -np.inf], np.nan).dropna()

    fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)
    sc = ax.scatter(data[param1], data[param2], c=data["err_tc"], cmap="viridis_r", s=14, alpha=0.7, zorder=3)
    fig.colorbar(sc, ax=ax, label="error (err_tc)")
    ax.scatter(best[param1], best[param2], marker="*", s=220, color="red", edgecolor="black", zorder=5,
               label=f"best fit (error={best['err_tc']:.5f})")
    bounds = dict(zip(run["config"]["optimized_parameters"], run["config"]["bounds"]))
    for name, set_lim in ((param1, ax.set_xlim), (param2, ax.set_ylim)):
        if name in bounds:
            set_lim(*bounds[name])
    ax.set_xlabel(param1)
    ax.set_ylabel(param2)
    ax.set_title(f"2D error landscape: {param1} vs {param2}")
    return fig


def plot_3D_errorlandscape(run, param1, param2, param3):
    """Error (err_tc) over three optimized parameters: every evaluation as a point coloured by error, the best fit as a star."""
    ev = run["evaluations"]
    best = ev.loc[ev["err_tc"].idxmin()]
    data = ev[[param1, param2, param3, "err_tc"]].astype(float).replace([np.inf, -np.inf], np.nan).dropna()

    fig = plt.figure(figsize=(11, 9), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(data[param1], data[param2], data[param3], c=data["err_tc"], cmap="viridis_r", s=18, alpha=0.35, edgecolors="none")
    fig.colorbar(sc, ax=ax, shrink=0.65, pad=0.1, label="error (err_tc)")
    ax.scatter(best[param1], best[param2], best[param3], marker="*", s=800, color="red", edgecolor="black",
               linewidth=2.5, depthshade=False, zorder=100, label=f"best fit (error={best['err_tc']:.5f})")
    bounds = dict(zip(run["config"]["optimized_parameters"], run["config"]["bounds"]))
    for name, set_lim in ((param1, ax.set_xlim), (param2, ax.set_ylim), (param3, ax.set_zlim)):
        if name in bounds:
            set_lim(*bounds[name])
    ax.set_xlabel(param1, labelpad=10)
    ax.set_ylabel(param2, labelpad=10)
    ax.set_zlabel(param3, labelpad=10)
    ax.set_title(f"3D parameter space: {param1}, {param2}, {param3}")
    return fig


plot1 = True
plot2 = False
plot3 = False
plot4 = False
plot5 = True
plot6 = False


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

if plot5:
    # Plot 5 - spectrum of the pre-stimulus activity
    # Welch spectrum of the three ROI dipoles over the settled period before the stimulus
    plot_prestim_spectrum(model, sim_dip)
    plt.show()

if plot6:
    # Plot 6 - error landscapes of the run: err_tc over two and over three optimized parameters
    plot_2D_errorlandscape(run, "coupling_strength", "strength_I")
    plt.show()
    plot_3D_errorlandscape(run, "ff_A1_S2", "ff_A3b_A1", "ff_A3b_S2")
    plt.show()
