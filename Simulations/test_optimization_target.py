"""
File: test_optimization_target.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-09-30
Description: In this script I test how different error function behave. 
For that I load previous optimization results of simulated curves and compute the error.
We only need to load the signal from one ROI only. 
Based on that I will choose the most suitable error function. 
"""

# %%
# imports 
import numpy as np
import os
import sys
import h5py
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

DATADIR = os.getenv('DATADIR')
RECONDIR = os.getenv('SUBJECTS_DIR')
SIMDIR = os.getenv("SIMDIR")
WDDIR = os.getenv("WDDIR")
RESDIR = os.getenv("RESDIR")
# add model to datapath
sys.path.append(os.path.join(WDDIR, 'Simulations', 'model'))
from somato_model import SomatoModel, read_simulation_params, load_optimized_params

# %%

# define different error functions

def error_1(simulated_pre, simulated_post, sim_peak, target_pre, target_post, target_peak):
    """
    Pre-Error: RMS (only simulation, no error) 
    Post-Error: RMSE
    """
    error_post = float(np.sqrt(np.mean((simulated_post / sim_peak - target_post / target_peak) ** 2)))
    error_pre = float(np.sqrt(np.mean((simulated_pre / sim_peak) ** 2)))
    return error_pre, error_post

def error_2(simulated_pre, simulated_post, sim_peak, target_pre, target_post, target_peak):
    """
    Pre-Error: RMS squared (only simulation, no error) 
    Post-Error: RMSE
    """
    error_post = float(np.sqrt(np.mean((simulated_post / sim_peak - target_post / target_peak) ** 2)))
    error_pre = float(np.sqrt(np.mean((simulated_pre / sim_peak) ** 2)))**2
    return error_pre, error_post

def error_3(simulated_pre, simulated_post, sim_peak, target_pre, target_post, target_peak):
    """
    Pre-Error: RMS normalized by post (only simulation, no error) 
    Post-Error: RMSE
    """
    error_post = float(np.sqrt(np.mean((simulated_post / sim_peak - target_post / target_peak) ** 2)))
    error_pre = (float(np.sqrt(np.mean((simulated_pre / sim_peak) ** 2)))/float(np.sqrt(np.mean((simulated_post / sim_peak) ** 2))))**2
    return error_pre, error_post

def error_4(simulated_pre, simulated_post, sim_peak, target_pre, target_post, target_peak):
    """
    Pre-Error: RMS squared (only simulation, no error) 
    Post-Error: RMSE
    WEIGHTED
    """
    weight_pre = 5
    error_post = float(np.sqrt(np.mean((simulated_post / sim_peak - target_post / target_peak) ** 2)))
    error_pre = weight_pre * float(np.sqrt(np.mean((simulated_pre / sim_peak) ** 2)))**2
    return error_pre, error_post


def error_5(simulated_pre, simulated_post, sim_peak, target_pre, target_post, target_peak):
    """
    Pre-Error: pre/post RMS ratio with ratio constraint
    Post-Error: RMSE
    """
    r_max = 0.2
    error_post = float(np.sqrt(np.mean((simulated_post / sim_peak - target_post / target_peak) ** 2)))
    rms_pre = float(np.sqrt(np.mean((simulated_pre / sim_peak) ** 2)))
    rms_post = float(np.sqrt(np.mean((simulated_post / sim_peak) ** 2)))
    error_pre = (rms_pre / rms_post - r_max) ** 2
    return error_pre, error_post

def error_6(simulated_pre, simulated_post, sim_peak, target_pre, target_post, target_peak):
    """
    Pre-Error: RMS squared (only simulation, no error) 
    Post-Error: RMSE
    """
    error_post = float(np.sqrt(np.mean((simulated_post / sim_peak - target_post / target_peak) ** 2)))
    error_sim_pre = float(np.sqrt(np.mean((simulated_pre/sim_peak - np.mean(simulated_pre/sim_peak)) ** 2)))
    error_target_pre = float(np.sqrt(np.mean((target_pre/ target_peak - np.mean(target_pre/ target_peak)) ** 2)))
    error_pre = ((error_sim_pre - error_target_pre)/error_target_pre)**2
    
    return error_pre, error_post


error_functions = {
    "error_1": error_1,
    "error_2": error_2,
    "error_3": error_3,
    "error_4": error_4,
    "error_5": error_5,
    "error_6": error_6,
}


#%%

# paths of previous optimization results
OPTDIR = os.path.join(SIMDIR, 'optimization')
opt_paths = [os.path.join(OPTDIR, d) for d in [
    "opt_20260818_171608_tc_roi-S2",
    "opt_20260818_171706_tc_roi-all",
    "opt_20260819_110148_tc_roi-A1",
    "opt_20260825_173230_ps_roi-all",
    "opt_20260826_090627_ps_roi-all",
    "opt_20260826_111151_ps_roi-A3b",
    "opt_20260826_112159_ps_roi-A3b",
    "opt_20260826_120040_ps_roi-A1",
    "opt_20260827_153421_tc_roi-A1",
    "opt_20260918_183022_tc_roi-A1",
    "opt_20260930_210433_tc_roi-A1"
]]
roi = "A1"

eps = 1e-10
rows = []
traces = {}

# loop over all paths
for p in opt_paths[-6:]:

    run = os.path.basename(p).removeprefix("opt_")
    tc_file = os.path.join(p, "best_tc_comparison.hdf5")
    if not os.path.exists(tc_file):
        print(f"skipping {run}: no best_tc_comparison.hdf5")
        continue

    # load preivous optimization results
    with h5py.File(tc_file, "r") as f:
        times_ms = f["times_ms"][:]
        simulated = f["sim"][roi][:]
        target = f["target"][roi][:]

    # split at stimulus onset (0 ms)
    i0 = int(np.argmin(np.abs(times_ms)))
    sim_pre, sim_post, tgt_pre, tgt_post = simulated[:i0], simulated[i0:], target[:i0], target[i0:]
    sim_peak = np.max(np.abs(sim_post)) + eps
    tgt_peak = np.max(np.abs(tgt_post)) + eps
    traces[run] = (times_ms, simulated / sim_peak, target / tgt_peak, i0)

    # compute the errors
    for name, fn in error_functions.items():
        err_pre, err_post = fn(sim_pre, sim_post, sim_peak, tgt_pre, tgt_post, tgt_peak)
        rows.append(dict(run=run, error=name, pre=err_pre, post=err_post, total=err_pre + err_post))

# print error summary
df = pd.DataFrame(rows)
for col in ("pre", "post", "total"):
    print(f"\n{col} error, ROI {roi}")
    print(df.pivot(index="run", columns="error", values=col).round(4).to_string())


# plot the error comparison
runs = list(traces)
error_names = list(error_functions)
colors = plt.get_cmap("tab10")
fig, axes = plt.subplots(2, len(runs), figsize=(4 * len(runs), 8), squeeze=False)

for ax, run in zip(axes[0], runs):
    times_ms, sim_norm, tgt_norm, i0 = traces[run]
    ax.plot(times_ms, tgt_norm, color="k", label="target")
    ax.plot(times_ms, sim_norm, color="C3", label="simulated")
    ax.axvline(0, color="gray", ls="--", lw=0.8)
    ax.set_title(run, fontsize=8)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("peak-normalised amplitude")
axes[0, 0].legend()

# bottom row: one panel per run, per error function a pre bar (hatched, light) next to a
# post bar (solid); colour identifies the error function
width = 0.4
for ax, run in zip(axes[1], runs):
    sub = df[df.run == run].set_index("error").loc[error_names]
    x = np.arange(len(error_names))
    bar_colors = [colors(i) for i in x]
    ax.bar(x - width / 2, sub["pre"], width, color=bar_colors, alpha=0.5, hatch="//")
    ax.bar(x + width / 2, sub["post"], width, color=bar_colors)
    ax.set_xticks(x, error_names, rotation=45, ha="right")
    ax.set_ylabel("error")
axes[1, 0].legend(handles=[
    Patch(facecolor="lightgray", edgecolor="gray", hatch="//", label="pre"),
    Patch(facecolor="gray", label="post"),
])

fig.tight_layout()
plt.show()
# %%
