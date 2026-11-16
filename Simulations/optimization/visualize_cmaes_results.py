"""
File: visualize_cmaes_results.py
Description: Diagnostic plots for a session folder written by run_optimization_cmaes.py
(optimization_timecourse/optcma_<id>_roi-all/). Independent of visualize_optmization_results.py (which plots
the sessions of run_optimization_full.py); the plot functions are copies, changed where the CMA-ES run differs:
the ERP error is computed on the average of several simulations (run_config.json "sims"), the 8 dipole lengths
are solved after the simulations (F_L5 and F_L6 are fixed), and the baseline rate penalty is the linear
distance from the band. Figures are saved to <session>/diagnostics/:

  01_convergence.png                total error, its three terms and the spread per CMA-ES generation
  01b_parameter_evolution.png       the searched parameters against the evaluation number
  01c_dipole_length_evolution.png   the solved dipole lengths against the evaluation number
  01d_parameter_error_correlation.png  Spearman and Pearson correlation of every parameter with the error
  02_roi_errors.png                 error terms per ROI of the best fit (average trace) and the time courses
  02b_timecourse_comparison.png     model.plot_timecourse_comparison: simulated vs measured, global-peak normalised
  02c_simulations.png               the single simulations of the best fit and their error against the average
  03_dipole_parameters_best.png     dipole length x orientation per population (best fit)
  03b_dipole_parameters_10.png      the 10 dipole parameters inside their bounds
  04_dipole_computation_<area>.png, 04b_layer_interneuron_contributions_<area>.png, 04c_dipole_length_by_celltype.png,
  single_simulations/all_potentials_*.pdf   how the potentials become the dipole (PLOT_DETAILS)
  05_connectivity_neural.png        connectivity between the 34 populations of the best fit
  06_connectivity_inputs.png        weights of the background / external / modulatory input
  07_saturation.png                 baseline rate (with the penalty band) and post-stimulus range of every population
  07b_sigmoid_operating_points.png  sigmoid of every population with the visited drive
  07c_prestim_spectrum.png          power spectrum of the pre-stimulus activity of the three ROI dipoles
  08_error_landscape_2D_*.png, 08_error_landscape_3D_*.png   evaluations coloured by error

01 reads only the files of the run. 02-07 re-simulate the best parameter set as run_optimization_cmaes.py does
at its end (all simulations of "sims", dipole lengths as solved) and stop if the total error differs from the best
error in evaluations.csv. The single-simulation figures (04-07) use the first simulation of "sims".

Runs made before "sims" was saved (one simulation, squared distance from the band) are re-simulated that way.

Usage (WDDIR and SIMDIR must be exported):
    python visualize_cmaes_results.py                          # session OPT_RUN_ID below
    python visualize_cmaes_results.py optcma_20261114_173224_roi-all   # name or path of a session
"""
# %%
import os
os.environ.setdefault("MPLBACKEND", "Agg")      # headless by default; set MPLBACKEND to look at the figures

import json
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.colors as mcolors

# ── settings ───────────────────────────────────────────────────────────────────
OPT_RUN_ID = "optcma_20261113_173145_roi-all"   # session folder name or path
CHECK_TOL = 1e-6        # allowed difference between the re-simulated and the logged best error
PARAM_TOL = 0.01        # two evaluations count as the same point if no parameter differs by more than this fraction of its range (correlation figure)
PLOT_DETAILS = True     # figures 04 (dipole computation, per-synapse potentials)
PARAM_2D = ("strength_I", "Iext_strength")                          # error landscapes (searched parameters)
PARAM_3D = ("coupling_strength", "strength_I", "g_intercortical")

# ── paths ──────────────────────────────────────────────────────────────────────
WDDIR = os.getenv("WDDIR")      # .../SomatosensoryLaminarModel
SIMDIR = os.getenv("SIMDIR")    # output directory of the simulations
if not WDDIR or not SIMDIR:
    raise SystemExit("WDDIR and SIMDIR must be exported")
OPT_OUT_DIR = os.path.join(SIMDIR, "optimization_timecourse")

for _p in (os.path.join(WDDIR, "Simulations"), os.path.join(WDDIR, "Simulations", "model"),
           os.path.join(WDDIR, "EEGSimulation")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import dipole_parameterization as dp
import somato_model as sm
import plotting_functions as pf
from plot_dipole_computation import (plot_dipole_computation_area, plot_dipole_parameters_clean,
                                    plot_layer_interneuron_contributions, _source_population_table,
                                    _SOURCE_AREA_COLORS)

# time axis of the ROI time courses: 2 ms samples, stimulus onset at sample 250 (as in compute_error_timecourse)
ONSET_SAMPLE, DT_MS = 250, 2
# block boundaries of the 34 populations: A3b | A1 | S2 | thalamus
BLOCKS = (("A3b", 0, 4), ("A1", 4, 17), ("S2", 17, 30), ("Thalamus", 30, 34))

# The 8 dipole lengths that run_optimization_cmaes.py solves after each simulation (F_L5 and F_L6 are fixed), with
# the bounds of the solve. Copy of the entries of PARAMETER_DICT in run_optimization_cmaes.py (run_config.json
# does not store them): keep the two in step.
DIPOLE_BOUNDS = {
    "L_BASE": (0.5, 15.0), "F_L4_A1": (0.01, 0.5), "F_L4_S2": (0.01, 0.5), "L_THALE_CTX": (0.5, 10.0),
    "L_THALPOM_CTX": (0.5, 20.0), "L_THALE_A3B": (0.02, 1.0), "L_THALPOM_A3B": (0.5, 20.0), "L_A3B_CORT": (0.5, 20.0),
}
SOLVED = [n for n in dp.PARAM_NAMES if n not in ("F_L5", "F_L6")]


# ── loading ────────────────────────────────────────────────────────────────────
def resolve_run_dir(run):
    """Session folder from a folder name or a path."""
    path = run if os.path.isdir(run) else os.path.join(OPT_OUT_DIR, run)
    if not os.path.isdir(path):
        raise FileNotFoundError(f"session folder not found: {path}")
    return path


def load_run(run_dir):
    """Files of a session: dict with config, evaluations, summary (None while the run is unfinished)."""
    def _read_json(name, required):
        path = os.path.join(run_dir, name)
        if not os.path.isfile(path):
            if required:
                raise FileNotFoundError(path)
            return None
        with open(path) as f:
            return json.load(f)

    config = _read_json("run_config.json", True)
    missing = [k for k in ("base_params", "fixed_params", "parameter_defaults", "optimized_parameters", "bounds", "fit_rois",
                           "target_csv", "cma", "rate_band", "mid_rate_penalty_weight", "saturation_threshold") if k not in config]
    if missing:
        raise KeyError(f"run_config.json lacks {missing}: this is not a session of run_optimization_cmaes.py")
    eval_path = os.path.join(run_dir, "evaluations.csv")
    if not os.path.isfile(eval_path):
        raise FileNotFoundError(eval_path)
    return {"dir": run_dir, "config": config, "evaluations": pd.read_csv(eval_path),
            "summary": _read_json("optimization_summary.json", False)}


def best_parameters(run):
    """(searched parameters, the 10 dipole parameters, best logged total error, source description)."""
    ev, cfg = run["evaluations"], run["config"]
    best_row = ev.loc[ev["err_tc"].idxmin()]
    if run["summary"] is not None:
        summary = run["summary"]
        searched = {n: float(summary["best_params"][n]) for n in cfg["optimized_parameters"]}
        dipole = {k: float(v) for k, v in summary["dipole_parameters"].items()}
        source = "optimization_summary.json"
    else:
        print("optimization_summary.json not found (run unfinished): using the best row of evaluations.csv")
        searched = {n: float(best_row[n]) for n in cfg["optimized_parameters"]}
        dipole = {k: float(best_row[k]) for k in dp.PARAM_NAMES}
        source = "evaluations.csv"
    return searched, dipole, float(best_row["err_tc"]), source


# ── re-simulation of the best parameter set ────────────────────────────────────
def tc_error(tc, target):
    """ERP error of the traces tc (n_roi, 376) against target, same arithmetic as SomatoModel.compute_error_timecourse."""
    sim_post, tgt_post = tc[:, ONSET_SAMPLE:], target[:, ONSET_SAMPLE:]
    sim_peak = np.max(np.abs(sim_post)) + 1e-10
    tgt_peak = np.max(np.abs(tgt_post)) + 1e-10
    err_post = np.sqrt(np.mean((sim_post / sim_peak - tgt_post / tgt_peak) ** 2, axis=1))
    err_pre = np.mean((tc[:, :ONSET_SAMPLE] / sim_peak) ** 2, axis=1)
    return float(np.mean(err_pre + err_post))


def rate_penalties(model, cfg, linear):
    """(baseline rate penalty, thalamus saturation penalty) of the simulation the model holds, as in run_optimization_cmaes.py.

    linear=True: distance from the band relative to the band edge, 0 inside the band and 1 at a rate of 0 or of the
    maximum. linear=False: squared distance (the first runs)."""
    stim_idx = round(model.input_onset / model.step_size)
    window = round(0.5 / model.step_size)                                 # last 500 ms before the stimulus
    frac = model.rate[:30, stim_idx - window:stim_idx] / model.sigm[:30, 2, np.newaxis]
    lo, hi = cfg["rate_band"]
    if linear:
        distance = np.maximum(lo - frac, 0) / lo + np.maximum(frac - hi, 0) / (1 - hi)
    else:
        distance = np.maximum(lo - frac, 0) ** 2 + np.maximum(frac - hi, 0) ** 2
    mid_rate_error = cfg["mid_rate_penalty_weight"] * np.mean(distance)
    stim_duration = round(model.Iext_duration / model.step_size)
    stim_fraction = model.rate[30, stim_idx:stim_idx + stim_duration] / model.sigm[30, 2]
    saturation_error = np.mean(np.maximum(stim_fraction - cfg["saturation_threshold"], 0) ** 2)
    return float(mid_rate_error), float(saturation_error)


def simulate_best(run):
    """Re-simulate the best parameters the way run_optimization_cmaes.py does at its end.

    Returns a dict: model and sim_dip of the FIRST simulation of "sims", rois, sims, traces (n_sims, n_roi, 376),
    tc_avg (n_roi, 376), target (n_roi, 376), single_errs, err_erp (of the average), mid_rate_error,
    saturation_error, total, dipole. Raises if the total differs from the best logged error by more than CHECK_TOL."""
    cfg = run["config"]
    searched, dipole, logged, source = best_parameters(run)
    rois = tuple(cfg["fit_rois"])
    linear = "sims" in cfg                                    # runs without "sims": one simulation, squared band penalty
    sims = [tuple(s) for s in cfg["sims"]] if linear else [(cfg["base_params"]["input_onset"], None)]

    fixed = {k: v for k, v in cfg["fixed_params"].items() if k not in dp.PARAM_NAMES and k != "scaling_factor"}
    model = sm.SomatoModel(cfg["base_params"])
    model.load_dipole_params = model.load_dipole_params_parameterized
    model.dipole_length_params = dipole
    target_dict = model.load_target_timecourse(cfg["target_csv"])
    target = np.array([target_dict[r] for r in rois])

    traces, penalties = [], []
    for onset, seed in reversed(sims):                        # the first simulation last: the model then holds it
        extra = {} if seed is None else {"input_onset": onset, "Ib_noise_seed": seed}
        model.apply_params({**fixed, **searched, **extra})
        model.initialize_state()
        model.simulate()
        tc = model.compute_timecourse(model.compute_dipoles())
        traces.append(np.array([tc[r] for r in rois]))
        penalties.append(rate_penalties(model, cfg, linear))
    traces, penalties = np.array(traces[::-1]), np.array(penalties[::-1])
    tc_avg = traces.mean(axis=0)
    err_erp = tc_error(tc_avg, target)
    mid_rate_error, saturation_error = penalties.mean(axis=0)
    total = err_erp + mid_rate_error + saturation_error
    if abs(total - logged) > CHECK_TOL:
        raise RuntimeError(f"re-simulated error {total:.8f} differs from the best error in evaluations.csv "
                           f"{logged:.8f} (best parameters from {source}); the run cannot be reproduced with the saved settings")
    print(f"re-simulated best fit: error {total:.6f} = ERP {err_erp:.6f} + rate {mid_rate_error:.6f} + thalamus "
          f"{saturation_error:.6f} over {len(sims)} simulation(s) (matches evaluations.csv, parameters from {source})")
    return {"model": model, "sim_dip": model.compute_dipoles(), "rois": rois, "sims": sims, "traces": traces, "tc_avg": tc_avg,
            "target": target, "single_errs": [tc_error(t, target) for t in traces], "err_erp": err_erp,
            "mid_rate_error": float(mid_rate_error), "saturation_error": float(saturation_error), "total": float(total),
            "dipole": dipole}


# ── 1: convergence ─────────────────────────────────────────────────────────────
def plot_convergence(run):
    """Error against the evaluation number, the terms of the best error so far, and the spread per generation.

    Top: every evaluation, the lowest total error so far, and the ERP part of that best evaluation. Bottom left:
    the three terms of the best evaluation so far (ERP error, baseline rate penalty, thalamus saturation penalty),
    stacked to the total. Bottom right: lowest, median and highest total error of each CMA-ES generation
    (popsize evaluations); the generations narrow when the search converges."""
    ev, cfg = run["evaluations"], run["config"]
    err = ev["err_tc"].to_numpy(dtype=float)
    n, x = len(err), ev["eval"].to_numpy()
    best = np.minimum.accumulate(err)
    best_row = np.maximum.accumulate(np.where(err == best, np.arange(n), 0))       # row of the best evaluation so far
    terms = (("ERP error", "err_erp"), ("baseline rate penalty", "mid_rate_error"), ("thalamus saturation", "saturation_error"))

    fig = plt.figure(figsize=(12, 9), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.3, 1])

    ax = fig.add_subplot(grid[0, :])
    ax.plot(x, err, ".", color="0.65", markersize=3, label="evaluation")
    ax.plot(x, best, color="C3", linewidth=1.6, label="best so far (total)")
    ax.plot(x, ev["err_erp"].to_numpy(dtype=float)[best_row], color="C0", linestyle="--", linewidth=1.4,
            label="ERP part of the best so far")
    ax.set_xlabel("evaluation")
    ax.set_ylabel("error")
    ax.set_ylim(top=np.nanpercentile(err, 98))
    ax.set_title(f"{os.path.basename(run['dir'])}: {n} evaluations, best error {best[-1]:.4f} (ERP {ev['err_erp'].iloc[best_row[-1]]:.4f}, "
                 f"rate {ev['mid_rate_error'].iloc[best_row[-1]]:.4f}, thalamus {ev['saturation_error'].iloc[best_row[-1]]:.4f})")
    ax.legend(loc="upper right")

    ax = fig.add_subplot(grid[1, 0])
    ax.stackplot(x, [ev[c].to_numpy(dtype=float)[best_row] for _, c in terms], labels=[t for t, _ in terms],
                 colors=["C0", "C1", "C2"])
    ax.set_xlabel("evaluation")
    ax.set_ylabel("error of the best so far")
    ax.set_title("terms of the best error so far")
    ax.legend(loc="upper right", fontsize=8)

    ax = fig.add_subplot(grid[1, 1])
    per_gen = pd.Series(err).groupby(np.arange(n) // cfg["cma"]["popsize"])
    gen = per_gen.min().index.to_numpy() + 1
    ax.fill_between(gen, per_gen.min(), per_gen.max(), color="0.85", label="lowest to highest")
    ax.plot(gen, per_gen.median(), color="0.4", label="median")
    ax.plot(gen, per_gen.min(), color="C3", label="lowest")
    ax.set_ylim(per_gen.min().min() * 0.95, np.nanpercentile(err, 98))
    ax.set_xlabel(f"generation ({cfg['cma']['popsize']} evaluations)")
    ax.set_ylabel("error")
    ax.set_title("spread within each generation")
    ax.legend(loc="upper right", fontsize=8)
    return fig


def plot_evolution(run, names, bounds, title):
    """The columns `names` of evaluations.csv against the evaluation number, coloured by the total error.

    The best evaluation is marked, dashed lines are the bounds (dict name -> (low, high))."""
    ev = run["evaluations"]
    err = ev["err_tc"].to_numpy(dtype=float)
    best_i = int(np.argmin(err))
    norm = mcolors.Normalize(vmin=np.nanpercentile(err, 2), vmax=np.nanpercentile(err, 90))
    n_cols = 4
    n_rows = int(np.ceil(len(names) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.6 * n_cols, 2.1 * n_rows + 0.8), sharex=True, squeeze=False,
                             constrained_layout=True)
    scatter = None
    for k, (name, ax) in enumerate(zip(names, axes.ravel())):
        scatter = ax.scatter(ev["eval"], ev[name], c=err, cmap="viridis_r", norm=norm, s=6)
        ax.scatter(ev["eval"].iloc[best_i], ev[name].iloc[best_i], marker="*", s=130, color="C3", edgecolor="k", zorder=5)
        lo, hi = bounds[name]
        for b in (lo, hi):
            ax.axhline(b, color="0.4", linestyle="--", linewidth=0.8)
        pad = 0.05 * (hi - lo)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_title(name, fontsize=9)
        ax.tick_params(labelsize=7)
        if k // n_cols == n_rows - 1:
            ax.set_xlabel("evaluation", fontsize=8)
    for ax in axes.ravel()[len(names):]:
        ax.axis("off")
    fig.colorbar(scatter, ax=axes, shrink=0.5, label="error (total)", location="right")
    fig.suptitle(title)
    return fig


# ── 2: per-ROI error ───────────────────────────────────────────────────────────
def roi_error_terms(tc_sim, tc_target, rois):
    """Pre- and post-stimulus error term per ROI, with the formulas of compute_error_timecourse.

    The error of an ROI is pre + post; the model's error is the mean over rois. Both traces are normalised by
    the peak of the post-stimulus part over rois."""
    eps = 1e-10
    sim_peak = max(np.max(np.abs(tc_sim[r][ONSET_SAMPLE:])) for r in rois) + eps
    tgt_peak = max(np.max(np.abs(tc_target[r][ONSET_SAMPLE:])) for r in rois) + eps
    terms = {}
    for r in rois:
        post = float(np.sqrt(np.mean((tc_sim[r][ONSET_SAMPLE:] / sim_peak - tc_target[r][ONSET_SAMPLE:] / tgt_peak) ** 2)))
        pre = float(np.mean((tc_sim[r][:ONSET_SAMPLE] / sim_peak) ** 2))      # penalises fluctuations before the stimulus
        terms[r] = (pre, post)
    return terms, sim_peak, tgt_peak


def plot_roi_errors(fit):
    """Per-ROI error terms (bars) above the simulated (average of the simulations) and measured time courses of each ROI."""
    rois = fit["rois"]
    tc_sim, tc_target = dict(zip(rois, fit["tc_avg"])), dict(zip(rois, fit["target"]))
    terms, sim_peak, tgt_peak = roi_error_terms(tc_sim, tc_target, rois)
    mean_check = float(np.mean([sum(terms[r]) for r in rois]))
    if abs(mean_check - fit["err_erp"]) > CHECK_TOL:
        raise RuntimeError(f"per-ROI terms average {mean_check:.8f} but the ERP error is {fit['err_erp']:.8f}")

    fig = plt.figure(figsize=(13, 7), constrained_layout=True)
    grid = fig.add_gridspec(2, len(rois), height_ratios=[1, 1.3])
    ax = fig.add_subplot(grid[0, :])
    x = np.arange(len(rois))
    pre = [terms[r][0] for r in rois]
    post = [terms[r][1] for r in rois]
    ax.bar(x, post, color="C0", label="post-stimulus term")
    ax.bar(x, pre, bottom=post, color="C1", label="pre-stimulus term")
    for i in range(len(rois)):
        ax.text(i, pre[i] + post[i], f"{pre[i] + post[i]:.3f}", ha="center", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(rois)
    ax.set_ylabel("error")
    ax.set_title(f"ERP error per ROI of the average of {len(fit['sims'])} simulation(s), mean over ROIs = {fit['err_erp']:.4f}")
    ax.legend()

    t = (np.arange(len(tc_sim[rois[0]])) - ONSET_SAMPLE) * DT_MS
    for col, r in enumerate(rois):
        ax = fig.add_subplot(grid[1, col])
        ax.plot(t, tc_target[r] / tgt_peak, color="k", label="measured")
        ax.plot(t, tc_sim[r] / sim_peak, color="C3", label="simulated")
        ax.axvline(0, color="0.5", linewidth=0.8)
        ax.set_title(r)
        ax.set_xlabel("time from stimulus (ms)")
        if col == 0:
            ax.set_ylabel("normalised amplitude")
            ax.legend()
    return fig


def plot_simulations(fit):
    """The single simulations of the best fit against their average (only for runs with several simulations).

    Top: the single traces (colour = simulation) with the average (red) and the measured trace (black), divided by
    the peak of the average as in 02. Bottom: ERP error of each single simulation (own peak normalisation) and of the
    average. An average error well below the single errors means the simulations differ in the response."""
    rois, traces, tc_avg, target = fit["rois"], fit["traces"], fit["tc_avg"], fit["target"]
    sim_peak = np.max(np.abs(tc_avg[:, ONSET_SAMPLE:]))
    tgt_peak = np.max(np.abs(target[:, ONSET_SAMPLE:]))
    t = (np.arange(tc_avg.shape[1]) - ONSET_SAMPLE) * DT_MS
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(traces)))
    labels = [f"stimulus {onset:.3f} s\nseed {seed}" for onset, seed in fit["sims"]]

    fig = plt.figure(figsize=(13, 8), constrained_layout=True)
    grid = fig.add_gridspec(2, len(rois), height_ratios=[1.3, 1])
    for col, r in enumerate(rois):
        ax = fig.add_subplot(grid[0, col])
        for k, tr in enumerate(traces):
            ax.plot(t, tr[col] / sim_peak, color=colors[k], linewidth=0.9, label=labels[k].replace("\n", ", "))
        ax.plot(t, tc_avg[col] / sim_peak, color="C3", linewidth=2.0, label="average")
        ax.plot(t, target[col] / tgt_peak, color="k", linewidth=1.6, label="measured")
        ax.axvline(0, color="0.5", linewidth=0.8)
        ax.set_title(r)
        ax.set_xlabel("time from stimulus (ms)")
        if col == 0:
            ax.set_ylabel("normalised amplitude")
            ax.legend(fontsize=7)
    ax = fig.add_subplot(grid[1, :])
    x = np.arange(len(traces) + 1)
    ax.bar(x[:-1], fit["single_errs"], color=colors)
    ax.bar(x[-1], fit["err_erp"], color="C3")
    for xi, v in zip(x, fit["single_errs"] + [fit["err_erp"]]):
        ax.text(xi, v, f"{v:.3f}", ha="center", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(labels + ["average"])
    ax.set_ylabel("ERP error")
    ax.set_title("ERP error of each simulation and of the average")
    return fig


# ── 3: dipole parameters ───────────────────────────────────────────────────────
def plot_dipole_parameter_values(run, dipole_params):
    """The 10 dipole parameters: position inside the bounds of the solve, default marked; F_L5 and F_L6 are fixed."""
    cfg = run["config"]
    fig, axes = plt.subplots(len(dp.PARAM_NAMES), 1, figsize=(9, 0.62 * len(dp.PARAM_NAMES) + 1.6), constrained_layout=True)
    for ax, name in zip(axes, dp.PARAM_NAMES):
        value, default = float(dipole_params[name]), float(cfg["parameter_defaults"][name])
        if name in DIPOLE_BOUNDS:
            lo, hi = DIPOLE_BOUNDS[name]
            ax.plot([lo, hi], [0, 0], color="0.6", linewidth=6, solid_capstyle="round")
            ax.plot(value, 0, "o", color="C3", markersize=9, label="best")
            ax.plot(default, 0, "|", color="k", markersize=16, markeredgewidth=2, label="default")
            ax.set_xlim(lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))
            ax.text(1.01, 0.5, f"{value:.4g}  [{lo:g}, {hi:g}]", transform=ax.transAxes, va="center", fontsize=8)
        else:
            ax.plot(value, 0, "s", color="0.3", markersize=9)
            ax.text(1.01, 0.5, f"{value:.4g}  (fixed)", transform=ax.transAxes, va="center", fontsize=8)
        ax.set_yticks([])
        ax.set_ylabel(name, rotation=0, ha="right", va="center", fontsize=9)
        ax.tick_params(labelsize=7)
        for side in ("left", "right", "top"):
            ax.spines[side].set_visible(False)
    axes[0].legend(loc="lower left", ncol=2, fontsize=8, bbox_to_anchor=(0, 1.0))
    fig.suptitle("dipole parameters of the best fit (grey bar = bounds of the solve; only the ratios of the lengths are determined)")
    return fig


# ── 7: saturation of the population firing rates ───────────────────────────────
def saturation_table(model):
    """Operating point of every population on its sigmoid (firing rate = r_max / (1 + exp(slope * (v_thr - v)))).

    v is the summed potential of all synapses of the population, as in model.compute_rates. Returns a dict of
    (34,) arrays: slope, v_thr, r_max, and for the pre-stimulus window (the last 500 ms before the stimulus, the window of
    the rate penalty) the mean v, drive x = slope * (v - v_thr) and rate fraction r / r_max, and for the post-stimulus
    window (250 ms after the stimulus, shifted by the receptor-thalamus delay as the ROI time course) the minimum /
    maximum of x and of the rate fraction; v, x, frac over -500..+250 ms. x = 0 is the half-maximum point,
    |x| > 2.94 means the rate is below 5% or above 95% of r_max."""
    stim = round(model.input_onset / model.step_size)
    delay = round(model.receptor_thalamus_delay / model.step_size)
    pre = slice(stim - 500, stim)
    post = slice(stim - 1 - delay, stim + 250 - delay)
    span = slice(stim - 1 - 500 - delay, post.stop)
    slope, v_thr, r_max = (np.asarray(model.sigm[:, k], dtype=float) for k in range(3))
    v = model.potential.sum(axis=1)                                   # (34, n_samples)
    x = slope[:, None] * (v - v_thr[:, None])
    frac = model.rate / r_max[:, None]
    return {"slope": slope, "v_thr": v_thr, "r_max": r_max,
            "v_pre": v[:, pre].mean(axis=1), "x_pre": x[:, pre].mean(axis=1), "frac_pre": frac[:, pre].mean(axis=1),
            "x_post_min": x[:, post].min(axis=1), "x_post_max": x[:, post].max(axis=1),
            "frac_post_min": frac[:, post].min(axis=1), "frac_post_max": frac[:, post].max(axis=1),
            "v": v[:, span], "x": x[:, span], "frac": frac[:, span]}


def plot_saturation(model, band):
    """Firing-rate saturation of all 34 populations, grouped by area (best fit, first simulation).

    Rows: (1) baseline rate (last 500 ms before the stimulus) as a fraction of the population's maximum rate, with the
    min-max range during the 250 ms after the stimulus as a line and the band of the rate penalty in green; (2) the
    same on the sigmoid input axis x = slope * (v - v_thr) (symmetric log axis, dotted lines at x = +-2.94 = 5% / 95% of
    the maximum rate); (3) how much of the rate range the stimulus uses (max - min of the rate fraction after the
    stimulus). Bars: red above 95% of the maximum, blue below 5%, orange outside the band (cortical populations
    only, the rate penalty does not include the thalamus), grey otherwise."""
    tab = saturation_table(model)
    labels = list(model.get_population_labels())
    n = len(labels)
    xs = np.arange(n)
    lo, hi = band
    sat_x = float(np.log(0.95 / 0.05))
    frac = tab["frac_pre"]
    outside = (xs < 30) & ((frac < lo) | (frac > hi))
    colors = np.where(frac > 0.95, "C3", np.where(frac < 0.05, "C0", np.where(outside, "C1", "0.6")))

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True, constrained_layout=True)

    ax = axes[0]
    ax.axhspan(lo, hi, color="C2", alpha=0.15, zorder=0)
    ax.bar(xs, frac, color=colors)
    ax.vlines(xs, tab["frac_post_min"], tab["frac_post_max"], color="k", linewidth=1.8)
    ax.axhline(0.95, color="0.3", linestyle=":", linewidth=0.9)
    ax.axhline(0.05, color="0.3", linestyle=":", linewidth=0.9)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("rate / maximum rate")
    ax.set_title(f"baseline firing rate (bars), its range 0-250 ms after the stimulus (black lines), penalty band {lo:g}-{hi:g} (green)")
    for i in xs:
        ax.text(i, 1.015, f"{tab['r_max'][i]:.0f}", ha="center", va="bottom", fontsize=5.5, color="0.35")
    ax.text(-0.02, 1.015, "r_max (Hz):", transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=6, color="0.35")

    ax = axes[1]
    ax.axhspan(-sat_x, sat_x, color="0.92", zorder=0)
    ax.bar(xs, tab["x_pre"], color=colors)
    ax.vlines(xs, tab["x_post_min"], tab["x_post_max"], color="k", linewidth=1.8)
    ax.axhline(0, color="0.3", linewidth=0.6)
    ax.set_yscale("symlog", linthresh=sat_x)
    y_lim = max(float(np.max(np.abs([tab["x_post_min"], tab["x_post_max"], tab["x_pre"]]))) * 1.3, 2 * sat_x)
    ax.set_ylim(-y_lim, y_lim)
    ticks = [t for t in (-1000, -100, -10, -sat_x, 0, sat_x, 10, 100, 1000) if abs(t) <= y_lim]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t:.2g}" for t in ticks])
    ax.minorticks_off()
    ax.set_ylabel("drive x = slope * (v - v_thr)")
    ax.set_title("position on the sigmoid (grey band: rate between 5% and 95% of the maximum)")

    ax = axes[2]
    swing = tab["frac_post_max"] - tab["frac_post_min"]
    ax.bar(xs, swing, color="C2")
    ax.set_ylabel("swing after stimulus\n(fraction of max rate)")
    ax.set_title("how much of its rate range each population uses in response to the stimulus")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=90, fontsize=8)

    for ax in axes:
        for _, start, _ in BLOCKS[1:]:
            ax.axvline(start - 0.5, color="k", linewidth=0.8)
        ax.set_xlim(-0.6, n - 0.4)
    for name, start, stop in BLOCKS:
        axes[2].text((start + stop - 1) / 2, 0.97, name, ha="center", va="top", fontsize=10, fontweight="bold",
                     transform=axes[2].get_xaxis_transform())
    n_in = int(np.sum((frac[:30] >= lo) & (frac[:30] <= hi)))
    n_high, n_low = int(np.sum(frac > 0.95)), int(np.sum(frac < 0.05))
    n_resp = int(np.sum(tab["frac_post_max"] > 0.95))
    fig.suptitle(f"saturation of the population firing rates, best fit: {n_in} of 30 cortical baselines inside the band, "
                 f"{n_high} of {n} above 95% and {n_low} below 5% of the maximum rate at baseline, "
                 f"{n_resp} reach above 95% after the stimulus")
    print(f"saturation: {n_in} of 30 cortical baselines inside the band {lo:g}-{hi:g}")
    print(f"saturation: {n_low} populations below 5% of r_max at baseline: {[str(labels[i]) for i in np.where(frac < 0.05)[0]]}")
    print(f"saturation: {n_resp} populations above 95% of r_max after the stimulus: "
          f"{[str(labels[i]) for i in np.where(tab['frac_post_max'] > 0.95)[0]]}")
    return fig


# ── copied from visualize_optmization_results.py: correlation of the parameters with the error (figure 01d) ──
def correlation_with_error(values, err, method):
    """Correlation of each column of values (DataFrame) with err (array); method is "spearman" or "pearson".

    Spearman uses ranks: robust to the skewed error (a few evaluations are far above the median) and to
    monotonic non-linearity. Pearson is linear and is dominated by the large errors; it also misses
    U-shaped dependencies. A constant column gives NaN."""
    if method == "spearman":
        values, err = values.rank(), pd.Series(err, index=values.index).rank()
    elif method == "pearson":
        err = pd.Series(err, index=values.index)
    else:
        raise ValueError(f"unknown correlation method: {method}")
    return values.apply(lambda col: col.corr(err))


def distinct_evaluations(ev, names, bounds, tol):
    """Row positions of ev that are distinct in parameter space, ordered by increasing error.

    Every parameter is scaled to its search bounds (0..1) and cut into bins of width tol. Evaluations in
    the same bin in ALL parameters count as the same point; the one with the lowest error is kept. Two
    evaluations are therefore distinct when at least one parameter differs by about more than tol of its
    range (bin edges make this approximate)."""
    bounds = np.asarray(bounds, dtype=float)
    scaled = (ev[names].to_numpy(dtype=float) - bounds[:, 0]) / (bounds[:, 1] - bounds[:, 0])
    cells = pd.DataFrame(np.floor(scaled / tol).astype(np.int64))
    order = np.argsort(ev["err_tc"].to_numpy(dtype=float), kind="stable")
    keep = ~cells.iloc[order].duplicated().to_numpy()
    return order[keep]


def plot_parameter_error_correlation(run, best_fraction=0.05, tol=PARAM_TOL):
    """Spearman (top row) and Pearson (bottom row) correlation of every searched parameter with the error,
    over all evaluations (left) and over the best_fraction of the distinct evaluations with the lowest error (right).

    The left column is dominated by the bad regions the optimizer left early on; it tells which parameters
    the error depends on at all. The right column only contains the region the optimizer settled in; it
    tells which parameters still matter there. Parameters are sorted by |Spearman rho| over all
    evaluations in every panel. A missing bar means the parameter is constant in that subset (correlation
    undefined). Positive: a larger value goes with a larger error.

    Distinct means: not in the same bin as a better evaluation in every parameter, with bins of width
    tol (fraction of the search range, see distinct_evaluations). Without this the best 5% are mostly
    near-copies from the final refinement, whose correlations reflect the order of the steps."""
    ev, names = run["evaluations"], run["config"]["optimized_parameters"]
    err = ev["err_tc"].to_numpy(dtype=float)
    distinct = distinct_evaluations(ev, names, run["config"]["bounds"], tol)    # ordered by error
    n_best = int(np.ceil(best_fraction * len(distinct)))
    best_idx = distinct[:n_best]
    print(f"correlation figure: {len(distinct)} of {len(ev)} evaluations are distinct (tol {tol:g} of the range), "
          f"best {100 * best_fraction:g}% = {n_best}")

    fig, axes = plt.subplots(2, 2, figsize=(11, 2 * (0.28 * len(names) + 1.2)), sharey=True,
                             constrained_layout=True)
    order = correlation_with_error(ev[names], err, "spearman").abs().sort_values(ascending=False).index.tolist()
    y = np.arange(len(order))
    subsets = ((ev[names], err, f"all {len(ev)} evaluations"),
               (ev[names].iloc[best_idx], err[best_idx], f"lowest {100 * best_fraction:g}% of {len(distinct)} distinct (n={n_best})"))
    for row, (method, label) in enumerate((("spearman", "Spearman rho"), ("pearson", "Pearson r"))):
        for col, (values, e, title) in enumerate(subsets):
            ax = axes[row, col]
            vals = correlation_with_error(values, e, method)[order].to_numpy(dtype=float)
            ax.barh(y, vals, color=np.where(vals >= 0, "C3", "C0"))
            ax.axvline(0, color="k", linewidth=0.8)
            ax.set_xlim(-1, 1)
            ax.set_title(f"{label}: {title}", fontsize=9)
            ax.set_xlabel(f"{label} (parameter, error)")
    for row in range(2):
        axes[row, 0].set_yticks(y)
        axes[row, 0].set_yticklabels(order, fontsize=8)
    axes[0, 0].invert_yaxis()
    return fig


# ── copied: dipole file of the best fit (03), dipole length per input cell type (04c), per-synapse potentials (04d) ──
def best_dipole_json(run, dipole_params):
    """Path of the best fit's dipole file: the run's own file, otherwise one built from the best parameters."""
    path = os.path.join(run["dir"], "dipole_parameters_finetuned.json")
    if os.path.isfile(path):
        return path
    out = os.path.join(run["dir"], "diagnostics", "dipole_parameters_best.json")
    return dp.build_params_json(out, {k: float(v) for k, v in dipole_params.items()})


def plot_dipole_length_by_celltype(model):
    """Dipole length (mm) of every input cell type onto the E population of each layer, all areas in one figure.

    Rows A3b / S1 / S2, columns the target layers (A3b is unlaminated: one panel). In a panel every source
    population is a dot at its length, coloured by source area, on the x position of its cell type; the tick
    marks the largest length of the cell type. Only the length parameter is shown: no orientation, cell count
    or potential. The lengths are those the model uses (model.load_dipole_params)."""
    lengths = model.load_dipole_params()[0]
    sources = _source_population_table(model)
    cell_types = ("E", "PV", "SST", "VIP")
    layer_names = ("L2/3", "L4", "L5", "L6")
    AREA_OFFSET = {"A3b": -0.21, "S1": -0.07, "S2": 0.07, "Thalamus": 0.21}    # dots of one cell type side by side
    rows = (("A3b", "A3b"), ("A1", "S1"), ("S2", "S2"))

    fig, axes = plt.subplots(len(rows), 4, figsize=(11, 8), sharey="row", squeeze=False, constrained_layout=True)
    for r, (key, display) in enumerate(rows):
        mat = np.atleast_2d(np.asarray(lengths[key], dtype=float))      # (n_layers, 34)
        for c in range(4):
            ax = axes[r, c]
            if c >= mat.shape[0]:
                ax.axis("off")
                continue
            for x, ct in enumerate(cell_types):
                members = [(j, rec) for j, rec in enumerate(sources) if rec[3] == ct]
                vals = [mat[c, j] for j, _ in members]
                for (j, rec), v in zip(members, vals):
                    ax.plot(x + AREA_OFFSET[rec[1]], v, "o", color=_SOURCE_AREA_COLORS[rec[1]],
                            markersize=6, markeredgecolor="k", markeredgewidth=0.4, alpha=0.9)
                ax.hlines(max(vals), x - 0.3, x + 0.3, color="k", linewidth=1.5)
            ax.axhline(0, color="0.6", linewidth=0.6)
            ax.set_xticks(range(len(cell_types)))
            ax.set_xticklabels(cell_types)
            ax.set_xlim(-0.6, len(cell_types) - 0.4)
            ax.set_title(f"{display} E" if mat.shape[0] == 1 else f"{display} {layer_names[c]} E", fontsize=10)
            if c == 0:
                ax.set_ylabel("dipole length (mm)")
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=col, markeredgecolor="k", label=f"source {name}")
               for name, col in _SOURCE_AREA_COLORS.items()]
    handles.append(plt.Line2D([0], [0], color="k", linewidth=1.5, label="largest length of the cell type"))
    fig.legend(handles=handles, loc="outside lower center", ncol=5, frameon=False, fontsize=8)
    fig.suptitle("dipole length of each input cell type onto the E population of each layer (best fit)")
    return fig


def plot_full_potentials(model, out_dir, before_s=0.2, after_s=0.3):
    """Per-synapse potential decomposition of every population (pf.plot_all_potentials, as in simulation_main.py).

    Window: before_s before to after_s after the input onset. pf.plot_all_potentials writes one PDF per
    area (A3b + thalamus, S1, S2) to <out_dir>/single_simulations/."""
    onset = int(model.input_onset / model.step_size)
    start = max(0, onset - int(before_s / model.step_size))
    stop = onset + int(after_s / model.step_size)
    pf.plot_all_potentials(
        model.potential, model.Iext[-2], model.Ib[0], model.step_size, model.simulation_dur, start, out_dir,
        model.strength_I, model.coupling_strength, model.Iext_duration, model.Ib_strength, model.Iext_strength,
        pop_labels=model.get_population_labels(), stop_plot=stop)
    print("saved", os.path.join(out_dir, "single_simulations", "all_potentials_*.pdf"))


# ── copied: connectivity (05, 06) ──
def _heatmap(ax, matrix, row_labels, col_labels):
    """Diverging heatmap centred at 0 with a symmetric logarithmic colour scale.

    The weights span about three orders of magnitude (thalamic weights up to ~1000, cortical ones tens),
    so a linear scale would show only the thalamus. Values below 1% of the maximum are drawn linearly."""
    vmax = np.max(np.abs(matrix))
    norm = mcolors.SymLogNorm(linthresh=0.01 * vmax, vmin=-vmax, vmax=vmax)
    image = ax.imshow(matrix, cmap="coolwarm", norm=norm, aspect="auto", interpolation="nearest")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=6)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=6, rotation=90)
    return image


def plot_connectivity_neural(model):
    """Connectivity between the 34 populations (rows = target, columns = source), block borders between areas."""
    labels = list(model.get_population_labels())
    w = model.W[:, :len(labels)]
    fig, ax = plt.subplots(figsize=(11, 9.5), constrained_layout=True)
    image = _heatmap(ax, w, labels, labels)
    for _, start, _ in BLOCKS[1:]:
        ax.axhline(start - 0.5, color="k", linewidth=0.8)
        ax.axvline(start - 0.5, color="k", linewidth=0.8)
    ax.set_xlabel("source population")
    ax.set_ylabel("target population")
    ax.set_title("connectivity W between populations (inhibitory weights negative)")
    fig.colorbar(image, ax=ax, shrink=0.6, label="weight")
    return fig


def plot_connectivity_inputs(model):
    """Weights of the three input channels onto every population (background, external, modulatory)."""
    labels = list(model.get_population_labels())
    w = model.W[:, len(labels):]
    cols = ["background", "external", "modulatory"]
    if w.shape[1] != len(cols):
        raise ValueError(f"expected {len(cols)} input columns in W, got {w.shape[1]}")
    fig, ax = plt.subplots(figsize=(4.5, 9.5), constrained_layout=True)
    image = _heatmap(ax, w, labels, cols)
    for _, start, _ in BLOCKS[1:]:
        ax.axhline(start - 0.5, color="k", linewidth=0.8)
    ax.set_ylabel("target population")
    ax.set_title("input weights")
    fig.colorbar(image, ax=ax, shrink=0.6, label="weight")
    return fig


# ── copied: sigmoid operating points (07b), pre-stimulus spectrum (07c) ──
def plot_sigmoid_operating_points(model):
    """Sigmoid of each of the 34 populations with the trajectory of its drive over -500..+250 ms (best fit).

    x axis: drive x = slope * (v - v_thr) on a symmetric log axis, so potentials far outside the steep part of
    the sigmoid stay visible; y axis: rate / maximum rate. Black line: sigmoid. Red line: the (x, rate) pairs
    visited from 500 ms before to 250 ms after the stimulus, circle: mean over the pre-stimulus part. A trajectory
    that sits on the flat part of the sigmoid cannot transmit changes of its input."""
    tab = saturation_table(model)
    labels = list(model.get_population_labels())
    n = len(labels)
    n_cols = 6
    n_rows = int(np.ceil(n / n_cols))
    lim = max(8.0, float(np.max(np.abs(tab["x"]))) * 1.3)
    grid = np.concatenate([-np.logspace(np.log10(lim), -2, 300), [0], np.logspace(-2, np.log10(lim), 300)])
    curve = 1 / (1 + np.exp(-grid))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.6 * n_cols, 2.0 * n_rows), sharey=True, constrained_layout=True)
    for i, ax in enumerate(axes.ravel()):
        if i >= n:
            ax.axis("off")
            continue
        ax.plot(grid, curve, color="k", linewidth=1.0)
        ax.plot(tab["x"][i], tab["frac"][i], color="C3", linewidth=1.6)
        ax.plot(tab["x_pre"][i], tab["frac_pre"][i], "o", color="C3", markeredgecolor="k", markersize=5, zorder=5)
        ax.set_xscale("symlog", linthresh=1.0)
        ax.set_xlim(-lim, lim)
        ax.set_xticks([t for t in (-1000, -100, -10, -1, 0, 1, 10, 100, 1000) if abs(t) <= lim])
        ax.set_xticklabels([f"{t:g}" for t in (-1000, -100, -10, -1, 0, 1, 10, 100, 1000) if abs(t) <= lim])
        ax.minorticks_off()
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(f"{labels[i]}  (r_max {tab['r_max'][i]:.0f} Hz)", fontsize=8)
        ax.tick_params(labelsize=6)
        if i // n_cols == n_rows - 1 or i + n_cols >= n:
            ax.set_xlabel("drive x", fontsize=7)
        if i % n_cols == 0:
            ax.set_ylabel("rate / max", fontsize=7)
    fig.suptitle("sigmoid operating points of all populations, best fit (red: -500..+250 ms, circle: baseline)")
    return fig


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


# ── copied: error landscapes (08) ──
def _landscape_bounds(run, names):
    """Search bounds for the selected parameters, if they were optimized."""
    cfg = run["config"]
    bounds = dict(zip(cfg["optimized_parameters"], cfg["bounds"]))
    return {name: bounds.get(name, None) for name in names}


def plot_2D_errorlandscape(run, param1, param2, levels=30):
    """2D error landscape of two optimized parameters.

    x/y: selected parameters
    colour: error (err_tc)

    Uses triangulation to interpolate between the sampled evaluations.
    Actual evaluations are shown as points and the best fit as a star.
    Duplicate parameter combinations are represented by their lowest error.
    """
    ev = run["evaluations"]
    best = ev.loc[ev["err_tc"].idxmin()]

    # Extract sampled parameter combinations
    x = ev[param1].to_numpy(dtype=float)
    y = ev[param2].to_numpy(dtype=float)
    err = ev["err_tc"].to_numpy(dtype=float)

    # Remove invalid evaluations
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(err)
    x, y, err = x[valid], y[valid], err[valid]

    # Aggregate duplicate parameter combinations using their minimum error
    data = pd.DataFrame({param1: x, param2: y, "err": err})
    data = data.groupby([param1, param2], as_index=False)["err"].min()

    x_unique = data[param1].to_numpy()
    y_unique = data[param2].to_numpy()
    err_unique = data["err"].to_numpy()

    fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)

    # Interpolated error landscape
    if len(data) >= 3:
        try:
            triang = mtri.Triangulation(x_unique, y_unique)
            contour = ax.tricontourf(
                triang, err_unique,
                levels=levels,
                cmap="viridis_r"
            )
            fig.colorbar(contour, ax=ax, label="error (err_tc)")
        except (ValueError, RuntimeError) as exc:
            print(f"Could not interpolate 2D landscape: {exc}")

    # Actual evaluations
    scatter = ax.scatter(
        x, y,
        c=err,
        cmap="viridis_r",
        s=12,
        edgecolors="none",
        alpha=0.65,
        zorder=3
    )

    # Best fit
    ax.scatter(
        best[param1], best[param2],
        marker="*",
        s=220,
        color="red",
        edgecolor="black",
        linewidth=0.8,
        label=f"best fit (error={best['err_tc']:.5f})",
        zorder=5
    )

    # Search bounds
    bounds = _landscape_bounds(run, (param1, param2))
    for name, axis in ((param1, ax.set_xlim), (param2, ax.set_ylim)):
        if bounds[name] is not None:
            lo, hi = bounds[name]
            axis(lo, hi)

    ax.set_xlabel(param1)
    ax.set_ylabel(param2)
    ax.set_title(f"2D error landscape: {param1} vs {param2}")
    #ax.legend(loc="best")

    return fig


def plot_3D_errorlandscape(run, param1, param2, param3):
    """3D parameter-space visualization of three parameters.

    x/y/z: selected parameters
    colour: error (err_tc)

    Each point represents an actual evaluation. The best fit is marked
    with a star. Unlike the 2D plot, no interpolation is performed.
    """
    ev = run["evaluations"]
    best = ev.loc[ev["err_tc"].idxmin()]

    x = ev[param1].to_numpy(dtype=float)
    y = ev[param2].to_numpy(dtype=float)
    z = ev[param3].to_numpy(dtype=float)
    err = ev["err_tc"].to_numpy(dtype=float)

    valid = (
        np.isfinite(x) & np.isfinite(y) &
        np.isfinite(z) & np.isfinite(err)
    )
    x, y, z, err = x[valid], y[valid], z[valid], err[valid]

    fig = plt.figure(figsize=(11, 9), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")

    # Actual evaluations, coloured by error
    scatter = ax.scatter(
        x, y, z,
        c=err,
        cmap="viridis_r",
        s=18,
        alpha=0.35,
        edgecolors="none"
    )

    fig.colorbar(scatter, ax=ax, shrink=0.65, pad=0.1,
                 label="error (err_tc)")

    # Best fit
    ax.scatter(
        best[param1], best[param2], best[param3],
        marker="*",
        s=800,
        color="red",
        edgecolor="black",
        linewidth=2.5,
        depthshade=False,
        label=f"best fit (error={best['err_tc']:.5f})",
        zorder=100
    )

    # Search bounds
    bounds = _landscape_bounds(run, (param1, param2, param3))
    for name, axis in zip(
        (param1, param2, param3),
        (ax.set_xlim, ax.set_ylim, ax.set_zlim)
    ):
        if bounds[name] is not None:
            lo, hi = bounds[name]
            axis(lo, hi)

    ax.set_xlabel(param1, labelpad=10)
    ax.set_ylabel(param2, labelpad=10)
    ax.set_zlabel(param3, labelpad=10)

    ax.set_title(
        f"3D parameter space: {param1}, {param2}, {param3}"
    )
    #ax.legend(loc="upper left")

    return fig


# ── main ──
def save(fig, out_dir, name, show=False):
    
    if show:
        plt.show()
    else:
        path = os.path.join(out_dir, name)
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print("saved", path)


def main():
    run_dir = os.path.join(OPT_OUT_DIR, OPT_RUN_ID)
    print("session:", run_dir)
    run = load_run(run_dir)
    cfg = run["config"]
    out_dir = os.path.join(run_dir, "diagnostics")
    os.makedirs(out_dir, exist_ok=True)

    # 1: convergence (files only)
    names = cfg["optimized_parameters"]
    save(plot_convergence(run), out_dir, "01_convergence.png")
    save(plot_evolution(run, names, dict(zip(names, cfg["bounds"])), "searched parameters (dashed: bounds, star: best)"),
         out_dir, "01b_parameter_evolution.png")
    save(plot_evolution(run, SOLVED, DIPOLE_BOUNDS, "dipole lengths solved after each simulation (dashed: bounds, star: best); "
                        "only the ratios are determined"), out_dir, "01c_dipole_length_evolution.png")
    save(plot_parameter_error_correlation(run), out_dir, "01d_parameter_error_correlation.png")

    # the following need the best fit simulated again
    fit = simulate_best(run)
    model, sim_dip, rois = fit["model"], fit["sim_dip"], fit["rois"]

    # 2: fit quality
    save(plot_roi_errors(fit), out_dir, "02_roi_errors.png")
    model.plot_timecourse_comparison(dict(zip(rois, fit["tc_avg"])), dict(zip(rois, fit["target"])), show=False,
                                     global_peak_rois=rois, save_path=os.path.join(out_dir, "02b_timecourse_comparison.png"))
    print("saved", os.path.join(out_dir, "02b_timecourse_comparison.png"))
    if len(fit["sims"]) > 1:
        save(plot_simulations(fit), out_dir, "02c_simulations.png")

    # 3: dipole parameters found
    save(plot_dipole_parameters_clean(best_dipole_json(run, fit["dipole"])), out_dir, "03_dipole_parameters_best.png")
    save(plot_dipole_parameter_values(run, fit["dipole"]), out_dir, "03b_dipole_parameters_10.png")

    if PLOT_DETAILS:
        # 4: how the dipole is computed from the potentials (first simulation)
        for area in ("S2", "A1", "A3b"):
            save(plot_dipole_computation_area(model, sim_dip, area), out_dir, f"04_dipole_computation_{area}.png")
            save(plot_layer_interneuron_contributions(model, sim_dip, area), out_dir,
                 f"04b_layer_interneuron_contributions_{area}.png")
        save(plot_dipole_length_by_celltype(model), out_dir, "04c_dipole_length_by_celltype.png")
        plot_full_potentials(model, out_dir)

    # 5/6: network structure
    save(plot_connectivity_neural(model), out_dir, "05_connectivity_neural.png")
    save(plot_connectivity_inputs(model), out_dir, "06_connectivity_inputs.png")

    # 7: are the populations saturated, are the baselines inside the band
    save(plot_saturation(model, cfg["rate_band"]), out_dir, "07_saturation.png")
    save(plot_sigmoid_operating_points(model), out_dir, "07b_sigmoid_operating_points.png")
    save(plot_prestim_spectrum(model, sim_dip, settle_s=min(2.0, 0.4 * model.input_onset)), out_dir, "07c_prestim_spectrum.png")

    # 8: error landscapes
    save(plot_2D_errorlandscape(run, *PARAM_2D), out_dir, f"08_error_landscape_2D_{PARAM_2D[0]}_{PARAM_2D[1]}.png")
    save(plot_3D_errorlandscape(run, *PARAM_3D), out_dir, f"08_error_landscape_3D_{PARAM_3D[0]}_{PARAM_3D[1]}_{PARAM_3D[2]}.png")


if __name__ == "__main__":
    main()
