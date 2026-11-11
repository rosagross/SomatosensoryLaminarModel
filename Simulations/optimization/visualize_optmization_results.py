"""
File: vizualize_optmization_results.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-11-05
Description: Diagnostic plots for a session folder written by run_optimization_full.py
(optimization_timecourse/optfull_<id>_roi-<tag>/). The figures are numbered in the order
in which they answer the questions about a run, and saved to <session>/diagnostics/:

  01_parameter_evolution.png        did the optimization converge, which parameters sit at a bound
  01b_parameter_error_correlation.png  Spearman and Pearson correlation of every parameter with the error (all evaluations, best 5%)
  02_roi_errors.png                 how well does the best fit match the measured time courses, per ROI
  02b_timecourse_comparison.png     model.plot_timecourse_comparison: simulated vs measured, global-peak normalised as in the error
  03_dipole_parameters_best.png     dipole length x orientation per population (best fit)
  03b_dipole_parameters_10.png      the 10 dipole parameters inside their search bounds
  04_dipole_computation_<area>.png  how the best fit's membrane potentials become the dipole (S2, A1, A3b)
  04b_layer_interneuron_contributions_<area>.png  which source populations (E/PV/SST/VIP) build each layer's dipole
  04c_dipole_length_by_celltype.png dipole length of each input cell type onto the E of every layer (all areas)
  single_simulations/all_potentials_*.pdf  04d: per-synapse potentials of every population (A3b+thalamus, S1, S2)
  05_connectivity_neural.png        connectivity between the 34 populations of the best fit
  06_connectivity_inputs.png        weights of the background / external / modulatory input
  07_saturation.png                 baseline firing rate and post-stimulus range of every population relative to its maximum rate
  07b_sigmoid_operating_points.png  sigmoid of every population with the visited drive, to see whether it sits on a flat part

Figure 01 only reads the files of the run. 02-07 re-simulate the best parameter set; the
re-simulated error must match the best error in evaluations.csv, otherwise the script stops.

Files read from the session folder: run_config.json, evaluations.csv and, once the run has
finished, optimization_summary.json and dipole_parameters_finetuned.json. run_config.json must
contain base_params, fixed_params and parameter_defaults, which run_optimization_full.py saves.

Usage (WDDIR and SIMDIR must be exported):
    python visualize_optmization_results.py                      # newest optfull_* session
    python visualize_optmization_results.py --run optfull_20261109_163717_roi-all
"""
# %%
import os
os.environ.setdefault("MPLBACKEND", "Agg")      # headless by default; set MPLBACKEND to look at the figures

import argparse
import glob
import json
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ── settings ───────────────────────────────────────────────────────────────────
OPT_RUN_ID = "optfull_20261110_174017_roi-all" #optfull_20261110_233754_roi-all"    # session folder name (or path); empty = newest optfull_* session
CHECK_TOL = 1e-6        # allowed difference between the re-simulated and the logged best error
PARAM_TOL = 0.01       # two evaluations count as the same point if no parameter differs by more than this fraction of its search range (correlation figure)

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


# ── loading ────────────────────────────────────────────────────────────────────
def resolve_run_dir(run):
    """Session folder from a folder name, a path, or (empty) the newest optfull_* folder."""
    if run:
        path = run if os.path.isdir(run) else os.path.join(OPT_OUT_DIR, run)
        if not os.path.isdir(path):
            raise FileNotFoundError(f"session folder not found: {path}")
        return path
    sessions = sorted(glob.glob(os.path.join(OPT_OUT_DIR, "optfull_*")))
    if not sessions:
        raise FileNotFoundError(f"no optfull_* session in {OPT_OUT_DIR}")
    return sessions[-1]


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
    missing = [k for k in ("base_params", "fixed_params", "parameter_defaults", "optimized_parameters", "bounds",
                           "fit_rois", "target_csv") if k not in config]
    if missing:
        raise KeyError(f"run_config.json lacks {missing}: the run was made before run_optimization_full.py "
                       f"saved them, so it cannot be re-simulated")
    eval_path = os.path.join(run_dir, "evaluations.csv")
    if not os.path.isfile(eval_path):
        raise FileNotFoundError(eval_path)
    return {"dir": run_dir, "config": config, "evaluations": pd.read_csv(eval_path),
            "summary": _read_json("optimization_summary.json", False)}


def best_parameters(run):
    """(best_params dict of the searched parameters, best logged error, source description)."""
    ev, names = run["evaluations"], run["config"]["optimized_parameters"]
    best_row = ev.loc[ev["err_tc"].idxmin()]
    if run["summary"] is not None:
        params = run["summary"]["best_params"]
        return {k: float(v) for k, v in params.items()}, float(best_row["err_tc"]), "optimization_summary.json"
    print("optimization_summary.json not found (run unfinished): using the best row of evaluations.csv")
    return {n: float(best_row[n]) for n in names}, float(best_row["err_tc"]), "evaluations.csv"


# ── re-simulation of the best parameter set ────────────────────────────────────
def resimulate_best(run):
    """Re-simulate the best parameters the way objective() of run_optimization_full.py does.

    Returns (model, sim_dip, err, tc_sim, tc_target, dipole_params). Raises if the error differs from the
    best logged error by more than CHECK_TOL."""
    cfg = run["config"]
    best, best_logged_err, source = best_parameters(run)

    # searched values override the fixed ones; dipole parameters and scaling_factor are not model parameters
    all_params = {**cfg["fixed_params"], **best}
    dipole_params = {k: all_params.pop(k) for k in dp.PARAM_NAMES}
    scaling_factor = all_params.pop("scaling_factor")

    model = sm.SomatoModel(cfg["base_params"])
    model.apply_params(all_params)
    model.dipole_length_params = dipole_params
    model.load_dipole_params = model.load_dipole_params_parameterized
    model.initialize_state()
    model.simulate()
    sim_dip = model.compute_dipoles()
    err, tc_sim, tc_target = model.compute_error_timecourse(
        cfg["target_csv"], sim_dip, scaling_factor=scaling_factor, target_dip=None, rois=tuple(cfg["fit_rois"]))

    if abs(err - best_logged_err) > CHECK_TOL:
        raise RuntimeError(f"re-simulated error {err:.8f} differs from the best error in evaluations.csv "
                           f"{best_logged_err:.8f} (best parameters from {source}); the run cannot be reproduced "
                           f"with the saved settings")
    print(f"re-simulated best fit: error {err:.6f} (matches evaluations.csv, parameters from {source})")
    return model, sim_dip, err, tc_sim, tc_target, dipole_params


# ── 1: parameter evolution ─────────────────────────────────────────────────────
def plot_parameter_evolution(run):
    """Error and every searched parameter against the evaluation number.

    The x axis is the order of the objective calls (initial population, offspring, gradient steps), because
    the optimizer's own iteration count is 1 for these runs. Points are coloured by their error, the best
    evaluation is marked, dashed lines are the search bounds."""
    ev, cfg = run["evaluations"], run["config"]
    names, bounds = cfg["optimized_parameters"], np.asarray(cfg["bounds"], dtype=float)
    err = ev["err_tc"].to_numpy()
    best_i = int(np.argmin(err))
    norm = mcolors.Normalize(vmin=np.nanpercentile(err, 2), vmax=np.nanpercentile(err, 90))

    n_cols = 4
    n_rows = int(np.ceil(len(names) / n_cols))
    fig = plt.figure(figsize=(3.6 * n_cols, 2.1 * n_rows + 3.2), constrained_layout=True)
    grid = fig.add_gridspec(n_rows + 1, n_cols, height_ratios=[1.6] + [1] * n_rows)

    ax = fig.add_subplot(grid[0, :])
    ax.plot(ev["eval"], err, ".", color="0.65", markersize=3, label="evaluation")
    ax.plot(ev["eval"], np.minimum.accumulate(err), color="C3", linewidth=1.6, label="best so far")
    ax.set_xlabel("evaluation")
    ax.set_ylabel("error (time course)")
    ax.set_title(f"{os.path.basename(run['dir'])}: {len(ev)} evaluations, best error {err[best_i]:.4f} "
                 f"(evaluation {int(ev['eval'].iloc[best_i])})")
    ax.set_ylim(top=np.nanpercentile(err, 98))
    ax.legend(loc="upper right")

    scatter = None
    param_axes = []
    for k, name in enumerate(names):
        ax = fig.add_subplot(grid[1 + k // n_cols, k % n_cols])
        param_axes.append(ax)
        scatter = ax.scatter(ev["eval"], ev[name], c=err, cmap="viridis_r", norm=norm, s=6)
        ax.scatter(ev["eval"].iloc[best_i], ev[name].iloc[best_i], marker="*", s=130, color="C3",
                   edgecolor="k", zorder=5)
        lo, hi = bounds[k]
        for b in (lo, hi):
            ax.axhline(b, color="0.4", linestyle="--", linewidth=0.8)
        pad = 0.05 * (hi - lo)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_title(name, fontsize=9)
        ax.tick_params(labelsize=7)
        if k // n_cols == n_rows - 1:
            ax.set_xlabel("evaluation", fontsize=8)
    fig.colorbar(scatter, ax=param_axes, shrink=0.5, label="error (time course)", location="right")
    return fig


# ── 1b: correlation of the parameters with the error ───────────────────────────
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


# ── 2: per-ROI error ───────────────────────────────────────────────────────────
def roi_error_terms(tc_sim, tc_target, fit_rois):
    """Pre- and post-stimulus error term per ROI, with the formulas of compute_error_timecourse.

    The error of an ROI is pre + post; the model's error is the mean over fit_rois. Both traces are
    normalised by the peak of the post-stimulus part over fit_rois."""
    eps = 1e-10
    rois = ("A3b", "A1", "S2")
    tgt_pre = {r: tc_target[r][:ONSET_SAMPLE] for r in rois}
    tgt_post = {r: tc_target[r][ONSET_SAMPLE:] for r in rois}
    sim_post = {r: tc_sim[r][ONSET_SAMPLE:][:tgt_post[r].shape[0]] for r in rois}
    sim_peak = max(np.max(np.abs(sim_post[r])) for r in fit_rois) + eps
    tgt_peak = max(np.max(np.abs(tgt_post[r])) for r in fit_rois) + eps
    terms = {}
    for r in rois:
        post = float(np.sqrt(np.mean((sim_post[r] / sim_peak - tgt_post[r] / tgt_peak) ** 2)))
        #pre_tgt = float(np.sqrt(np.mean((tgt_pre[r][:ONSET_SAMPLE] / tgt_peak -  np.mean(tgt_pre[r][:ONSET_SAMPLE] / tgt_peak)) ** 2)))
        #pre_sim = float(np.sqrt(np.mean((tc_sim[r][:ONSET_SAMPLE] / sim_peak -  np.mean(tc_sim[r][:ONSET_SAMPLE] / sim_peak)) ** 2)))
        #pre = ((pre_sim - pre_tgt)/pre_tgt)**2
        
        # classical error: favors no fluctuations
        pre = float(np.mean((tc_sim[r][:ONSET_SAMPLE] / sim_peak) ** 2))
        terms[r] = (pre, post)
    return terms, sim_peak, tgt_peak




def plot_roi_errors(run, err, tc_sim, tc_target):
    """Per-ROI error terms (bars) above the simulated and measured time courses of each ROI."""
    fit_rois = tuple(run["config"]["fit_rois"])
    terms, sim_peak, tgt_peak = roi_error_terms(tc_sim, tc_target, fit_rois)
    rois = ("A3b", "A1", "S2")
    mean_check = float(np.mean([sum(terms[r]) for r in fit_rois]))
    if abs(mean_check - err) > CHECK_TOL:
        raise RuntimeError(f"per-ROI terms average {mean_check:.8f} but the model error is {err:.8f}")

    fig = plt.figure(figsize=(13, 7), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=[1, 1.3])

    ax = fig.add_subplot(grid[0, :])
    x = np.arange(len(rois))
    pre = [terms[r][0] for r in rois]
    post = [terms[r][1] for r in rois]
    ax.bar(x, post, color="C0", label="post-stimulus term")
    ax.bar(x, pre, bottom=post, color="C1", label="pre-stimulus term")
    for i in range(len(rois)):
        ax.text(i, pre[i] + post[i], f"{pre[i] + post[i]:.3f}", ha="center", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels([r if r in fit_rois else f"{r} (not fitted)" for r in rois])
    ax.set_ylabel("error")
    ax.set_title(f"error per ROI, mean over fitted ROIs = {err:.4f}")
    ax.legend()

    t = (np.arange(len(tc_sim["A3b"])) - ONSET_SAMPLE) * DT_MS
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


# ── 3: dipole parameters ───────────────────────────────────────────────────────
def best_dipole_json(run, dipole_params):
    """Path of the best fit's dipole file: the run's own file, otherwise one built from the best parameters."""
    path = os.path.join(run["dir"], "dipole_parameters_finetuned.json")
    if os.path.isfile(path):
        return path
    out = os.path.join(run["dir"], "diagnostics", "dipole_parameters_best.json")
    return dp.build_params_json(out, {k: float(v) for k, v in dipole_params.items()})


def plot_dipole_parameter_values(run, dipole_params):
    """The 10 dipole parameters: position inside the search bounds (searched) or fixed value, default marked."""
    cfg = run["config"]
    bounds = dict(zip(cfg["optimized_parameters"], cfg["bounds"]))
    fig, axes = plt.subplots(len(dp.PARAM_NAMES), 1, figsize=(9, 0.62 * len(dp.PARAM_NAMES) + 1.6),
                             constrained_layout=True)
    for ax, name in zip(axes, dp.PARAM_NAMES):
        value, default = float(dipole_params[name]), float(cfg["parameter_defaults"][name])
        if name in bounds:
            lo, hi = bounds[name]
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
    fig.suptitle("dipole parameters of the best fit (grey bar = search bounds)")
    return fig


# ── 4c: dipole length per input cell type ──────────────────────────────────────
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


# ── 4d: per-synapse potentials ─────────────────────────────────────────────────
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


# ── 5/6: connectivity ──────────────────────────────────────────────────────────
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


# ── 7: saturation of the population firing rates ───────────────────────────────
def _roi_windows(model):
    """Sample slices (pre-stimulus 500 ms, post-stimulus 250 ms) of the ROI time axis of compute_timecourse."""
    stim_idx = round(model.input_onset / model.step_size) - 1
    delay = round(model.receptor_thalamus_delay / model.step_size)
    return slice(stim_idx - 500 - delay, stim_idx - delay), slice(stim_idx - delay, stim_idx + 251 - delay)


def saturation_table(model):
    """Operating point of every population on its sigmoid (firing rate = r_max / (1 + exp(slope * (v_thr - v)))).

    v is the summed potential of all synapses of the population, as in model.compute_rates. Returns a dict of
    (34,) arrays: slope, v_thr, r_max, and for the pre-stimulus window the mean v, drive x = slope * (v - v_thr)
    and rate fraction r / r_max, and for the post-stimulus window the minimum / maximum of x and of the rate
    fraction. x = 0 is the half-maximum point, |x| > 2.94 means the rate is below 5% or above 95% of r_max."""
    pre, post = _roi_windows(model)
    slope, v_thr, r_max = (np.asarray(model.sigm[:, k], dtype=float) for k in range(3))
    v = model.potential.sum(axis=1)                                   # (34, n_samples)
    x = slope[:, None] * (v - v_thr[:, None])
    frac = model.rate / r_max[:, None]
    return {"slope": slope, "v_thr": v_thr, "r_max": r_max,
            "v_pre": v[:, pre].mean(axis=1), "x_pre": x[:, pre].mean(axis=1), "frac_pre": frac[:, pre].mean(axis=1),
            "x_post_min": x[:, post].min(axis=1), "x_post_max": x[:, post].max(axis=1),
            "frac_post_min": frac[:, post].min(axis=1), "frac_post_max": frac[:, post].max(axis=1),
            "v": v[:, pre.start:post.stop], "x": x[:, pre.start:post.stop], "frac": frac[:, pre.start:post.stop],
            "n_pre": pre.stop - pre.start}


def plot_saturation(model):
    """Firing-rate saturation of all 34 populations, grouped by area (best fit).

    Rows: (1) pre-stimulus rate as a fraction of the population's maximum rate, with the min-max range
    during the 250 ms after the stimulus as a line; (2) the same on the sigmoid input axis x = slope * (v - v_thr)
    (symmetric log axis, dotted lines at x = +-2.94 = 5% / 95% of the maximum rate); (3) how much of the rate range
    the stimulus uses (max - min of the rate fraction after the stimulus). Bars are red when the baseline rate is
    above 95% of the maximum, blue below 5%, grey otherwise. Populations with a swing close to 0 do not
    respond to anything but their own external input."""
    tab = saturation_table(model)
    labels = list(model.get_population_labels())
    n = len(labels)
    xs = np.arange(n)
    sat_x = float(np.log(0.95 / 0.05))
    colors = np.where(tab["frac_pre"] > 0.95, "C3", np.where(tab["frac_pre"] < 0.05, "C0", "0.6"))

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True, constrained_layout=True)

    ax = axes[0]
    ax.bar(xs, tab["frac_pre"], color=colors)
    ax.vlines(xs, tab["frac_post_min"], tab["frac_post_max"], color="k", linewidth=1.8)
    ax.axhline(0.95, color="0.3", linestyle=":", linewidth=0.9)
    ax.axhline(0.05, color="0.3", linestyle=":", linewidth=0.9)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("rate / maximum rate")
    ax.set_title("baseline firing rate (bars) and its range 0-250 ms after the stimulus (black lines)")
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
    n_high, n_low = int(np.sum(tab["frac_pre"] > 0.95)), int(np.sum(tab["frac_pre"] < 0.05))
    fig.suptitle(f"saturation of the population firing rates, best fit: {n_high} of {n} populations above 95% and "
                 f"{n_low} below 5% of their maximum rate at baseline")
    print(f"saturation: {n_high} populations above 95% of r_max at baseline: "
          f"{[labels[i] for i in np.where(tab['frac_pre'] > 0.95)[0]]}")
    print(f"saturation: {n_low} populations below 5% of r_max at baseline: "
          f"{[labels[i] for i in np.where(tab['frac_pre'] < 0.05)[0]]}")
    return fig


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


# ── main ───────────────────────────────────────────────────────────────────────
def save(fig, out_dir, name):
    path = os.path.join(out_dir, name)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("saved", path)


def main(run_name=None):
    run_dir = resolve_run_dir(run_name if run_name else OPT_RUN_ID)
    print("session:", run_dir)
    run = load_run(run_dir)
    out_dir = os.path.join(run_dir, "diagnostics")
    os.makedirs(out_dir, exist_ok=True)

    # 1: convergence (files only)
    save(plot_parameter_evolution(run), out_dir, "01_parameter_evolution.png")
    save(plot_parameter_error_correlation(run), out_dir, "01b_parameter_error_correlation.png")

    # the following need the best fit simulated again
    model, sim_dip, err, tc_sim, tc_target, dipole_params = resimulate_best(run)

    # 2: fit quality
    save(plot_roi_errors(run, err, tc_sim, tc_target), out_dir, "02_roi_errors.png")

    # 2b: the model's own time course comparison, normalised with the fitted scaling factor
    model.plot_timecourse_comparison(tc_sim, tc_target, show=False, global_peak_rois=tuple(run["config"]["fit_rois"]),
                                     save_path=os.path.join(out_dir, "02b_timecourse_comparison.png"))
    print("saved", os.path.join(out_dir, "02b_timecourse_comparison.png"))

    # 3: dipole parameters found
    save(plot_dipole_parameters_clean(best_dipole_json(run, dipole_params)), out_dir, "03_dipole_parameters_best.png")
    save(plot_dipole_parameter_values(run, dipole_params), out_dir, "03b_dipole_parameters_10.png")

    # 4: how the dipole is computed from the potentials
    for area in ("S2", "A1", "A3b"):
        save(plot_dipole_computation_area(model, sim_dip, area), out_dir, f"04_dipole_computation_{area}.png")
        save(plot_layer_interneuron_contributions(model, sim_dip, area), out_dir,
             f"04b_layer_interneuron_contributions_{area}.png")

    # 4c: dipole length per input cell type
    save(plot_dipole_length_by_celltype(model), out_dir, "04c_dipole_length_by_celltype.png")

    # 4d: per-synapse potentials of all populations
    plot_full_potentials(model, out_dir)

    # 5/6: network structure
    save(plot_connectivity_neural(model), out_dir, "05_connectivity_neural.png")
    save(plot_connectivity_inputs(model), out_dir, "06_connectivity_inputs.png")

    # 7: are the populations saturated (do they respond to input at all)
    save(plot_saturation(model), out_dir, "07_saturation.png")
    save(plot_sigmoid_operating_points(model), out_dir, "07b_sigmoid_operating_points.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Diagnostic plots of a run_optimization_full.py session")
    parser.add_argument("--run", default=None, help="session folder name or path (default: OPT_RUN_ID, else newest)")
    main(parser.parse_args().run)
