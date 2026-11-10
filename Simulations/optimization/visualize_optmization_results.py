"""
File: vizualize_optmization_results.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-11-05
Description: Diagnostic plots for a session folder written by run_optimization_full.py
(optimization_timecourse/optfull_<id>_roi-<tag>/). The figures are numbered in the order
in which they answer the questions about a run, and saved to <session>/diagnostics/:

  01_parameter_evolution.png        did the optimization converge, which parameters sit at a bound
  02_roi_errors.png                 how well does the best fit match the measured time courses, per ROI
  02b_timecourse_comparison.png     model.plot_timecourse_comparison: simulated vs measured, global-peak normalised as in the error
  03_dipole_parameters_best.png     dipole length x orientation per population (best fit)
  03b_dipole_parameters_10.png      the 10 dipole parameters inside their search bounds
  04_dipole_computation_<area>.png  how the best fit's membrane potentials become the dipole (S2, A1, A3b)
  04b_layer_interneuron_contributions_<area>.png  which source populations (E/PV/SST/VIP) build each layer's dipole
  single_simulations/all_potentials_*.pdf  04c: per-synapse potentials of every population (A3b+thalamus, S1, S2)
  05_connectivity_neural.png        connectivity between the 34 populations of the best fit
  06_connectivity_inputs.png        weights of the background / external / modulatory input

Figure 01 only reads the files of the run. 02-06 re-simulate the best parameter set; the
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
OPT_RUN_ID = "optfull_20261110_013733_roi-all"         # session folder name (or path); empty = newest optfull_* session
CHECK_TOL = 1e-6        # allowed difference between the re-simulated and the logged best error

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
                                    plot_layer_interneuron_contributions)

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


# ── 2: per-ROI error ───────────────────────────────────────────────────────────
def roi_error_terms(tc_sim, tc_target, fit_rois):
    """Pre- and post-stimulus error term per ROI, with the formulas of compute_error_timecourse.

    The error of an ROI is pre + post; the model's error is the mean over fit_rois. Both traces are
    normalised by the peak of the post-stimulus part over fit_rois."""
    eps = 1e-10
    rois = ("A3b", "A1", "S2")
    tgt_post = {r: tc_target[r][ONSET_SAMPLE:] for r in rois}
    sim_post = {r: tc_sim[r][ONSET_SAMPLE:][:tgt_post[r].shape[0]] for r in rois}
    sim_peak = max(np.max(np.abs(sim_post[r])) for r in fit_rois) + eps
    tgt_peak = max(np.max(np.abs(tgt_post[r])) for r in fit_rois) + eps
    terms = {}
    for r in rois:
        post = float(np.sqrt(np.mean((sim_post[r] / sim_peak - tgt_post[r] / tgt_peak) ** 2)))
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


# ── 4c: per-synapse potentials ─────────────────────────────────────────────────
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

    # 4c: per-synapse potentials of all populations
    plot_full_potentials(model, out_dir)

    # 5/6: network structure
    save(plot_connectivity_neural(model), out_dir, "05_connectivity_neural.png")
    save(plot_connectivity_inputs(model), out_dir, "06_connectivity_inputs.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Diagnostic plots of a run_optimization_full.py session")
    parser.add_argument("--run", default=None, help="session folder name or path (default: OPT_RUN_ID, else newest)")
    main(parser.parse_args().run)
