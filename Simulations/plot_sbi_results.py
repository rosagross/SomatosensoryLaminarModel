"""
File: plot_sbi_results.py
Author: Rosa Grossmann
Description:
    Detailed plotting of a finished SBI (multi-round NPE) run produced by run_sbi.py.

    Reads a single run folder (the `sbi_<timestamp>/` directory containing posterior.pkl,
    simulations.npz, posterior_summary.csv and x_obs.npy) and writes a set of diagnostic
    figures into a `plots/` subfolder of that run:

        posterior_marginals.png  1D posterior per parameter (MAP, median, 90% CI, prior)
        posterior_corner.png     corner/pairplot with the MAP (and ground truth) marked
        feature_ppc.png          simulated feature distributions vs the observation x_obs
        sim_diagnostics.png      failed-sim count, distance-to-x_obs, best-matching theta
        map_timecourse.png       MAP re-simulation: dipole time course vs measured data
        map_timefreq.png         MAP re-simulation: Morlet TF maps vs measured data

    The last two re-run the SomatoModel at the best-fit (MAP) parameters, so they need
    RESDIR set (measured target CSVs); pass --no-resim to skip them.

    Usage:
        python plot_sbi_results.py [run_dir] [--no-resim]

    run_dir defaults to the most recent sbi_* folder under SBI_RESULTS_ROOT.
"""

import argparse
import glob
import os
import pickle
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

# ── paths / imports (mirror run_sbi.py) ───────────────────────────────────────────
WDDIR  = os.getenv("WDDIR")
RESDIR = os.getenv("RESDIR")

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(_HERE)                          # Simulations/ (for sbi modules)
sys.path.append(os.path.join(_HERE, "model"))   # Simulations/model/ (for somato_model.py)
if WDDIR:
    sys.path.append(os.path.join(WDDIR, "Simulations"))
    sys.path.append(os.path.join(WDDIR, "Simulations", "model"))

# Reuse the exact configuration the run used, rather than redefining it.
from run_sbi import (                                       # noqa: E402
    PARAM_NAMES, BOUNDS, BASE_PARAMS, TF_DATA_PATH, TC_DATA_PATH, RESULTS_ROOT,
    _BAD_FEATURES,
)
from sbi_features import ROIS, feature_names                # noqa: E402
from sbi.analysis import pairplot                           # noqa: E402

N_SAMPLES = 10000
_ROI_COLORS = {"A3b": "C0", "A1": "C1", "S2": "C2"}


# ── small headless plotting style (NOT plotting_style.figure_style: it calls tk.Tk
#    which crashes on a headless node) ──────────────────────────────────────────────
def _set_style():
    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.bbox": "tight",
        "figure.dpi": 110,
    })


# ── loading ────────────────────────────────────────────────────────────────────────
def _latest_run_dir():
    cands = sorted(glob.glob(os.path.join(RESULTS_ROOT, "sbi_*")))
    if not cands:
        raise SystemExit(f"No sbi_* runs found under {RESULTS_ROOT}; pass a run_dir explicitly.")
    return cands[-1]


def load_run(run_dir):
    """Load everything a run wrote. Returns a dict of arrays / objects."""
    npz = np.load(os.path.join(run_dir, "simulations.npz"), allow_pickle=True)
    theta = npz["theta"].astype(float)            # (n_sims, 6)
    x     = npz["x"].astype(float)                # (n_sims, 30)
    x_obs = npz["x_obs"].astype(float)            # (30,)

    with open(os.path.join(run_dir, "posterior.pkl"), "rb") as f:
        posterior = pickle.load(f)

    summary = _read_summary(os.path.join(run_dir, "posterior_summary.csv"))

    theta_true_path = os.path.join(run_dir, "theta_true.npy")
    theta_true = np.load(theta_true_path) if os.path.exists(theta_true_path) else None

    # Failed simulations are the constant -999 sentinel rows; drop them for feature plots.
    bad = np.all(np.isclose(x, _BAD_FEATURES), axis=1)
    return dict(theta=theta, x=x, x_obs=x_obs, posterior=posterior, summary=summary,
                theta_true=theta_true, bad_mask=bad)


def _read_summary(path):
    """posterior_summary.csv -> {param: {col: value}}, column order from the header."""
    out = {}
    with open(path) as f:
        header = f.readline().strip().split(",")
        for line in f:
            parts = line.strip().split(",")
            name = parts[0]
            out[name] = {k: float(v) for k, v in zip(header[1:], parts[1:])}
    return out


# ── figure 1: posterior marginals ────────────────────────────────────────────────────
def plot_marginals(data, out_path):
    samples = _posterior_samples(data)
    summary = data["summary"]
    theta_true = data["theta_true"]

    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5))
    for i, (name, ax) in enumerate(zip(PARAM_NAMES, axes.ravel())):
        s = summary[name]
        lo, hi = BOUNDS[i]
        ax.hist(samples[:, i], bins=50, range=(lo, hi), density=True,
                color="0.7", edgecolor="none")
        # flat uniform-prior reference
        ax.axhline(1.0 / (hi - lo), color="k", ls=":", lw=0.8, label="prior")
        ax.axvspan(s["p5"], s["p95"], color="C0", alpha=0.15, label="90% CI")
        ax.axvline(s["map"], color="C3", lw=1.5, label="MAP")
        ax.axvline(s["p50"], color="C0", lw=1.2, ls="--", label="median")
        if theta_true is not None:
            ax.axvline(theta_true[i], color="k", lw=1.5, label="truth")
        ax.set_xlim(lo, hi)
        ax.set_title(name)
        ax.set_yticks([])
        if i == 0:
            ax.legend(loc="upper right")
    fig.suptitle("Posterior marginals", y=1.0)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


# ── figure 2: corner / pairplot ──────────────────────────────────────────────────────
def plot_corner(data, out_path):
    samples = _posterior_samples(data)
    map_est = np.array([data["summary"][n]["map"] for n in PARAM_NAMES])
    points = [map_est]
    if data["theta_true"] is not None:
        points.append(data["theta_true"])
    fig, _ = pairplot(
        torch.as_tensor(samples, dtype=torch.float32),
        labels=PARAM_NAMES,
        limits=[[lo, hi] for lo, hi in BOUNDS],
        points=np.vstack(points),
        figsize=(11, 11),
    )
    fig.suptitle("Posterior corner (red = MAP)", y=1.0)
    fig.savefig(out_path)
    plt.close(fig)


# ── figure 3: feature posterior-predictive check ─────────────────────────────────────
def plot_feature_ppc(data, out_path):
    names = feature_names()
    x = data["x"][~data["bad_mask"]]              # (n_valid, 30)
    x_obs = data["x_obs"]
    n_feat = len(names)

    # z-score each feature by its simulated spread so the 30 features fit one axis.
    mu = x.mean(axis=0)
    sd = x.std(axis=0)
    sd[sd == 0] = 1.0
    xz = (x - mu) / sd
    obsz = (x_obs - mu) / sd

    roi_of = [n.split("_")[0] for n in names]
    colors = [_ROI_COLORS.get(r, "0.5") for r in roi_of]

    fig, ax = plt.subplots(figsize=(13, 6))
    parts = ax.violinplot([xz[:, j] for j in range(n_feat)],
                          positions=np.arange(n_feat), widths=0.8,
                          showextrema=False)
    for body, c in zip(parts["bodies"], colors):
        body.set_facecolor(c)
        body.set_alpha(0.35)
    ax.scatter(np.arange(n_feat), obsz, color="k", marker="D", s=22,
               zorder=5, label="x_obs")
    ax.axhline(0, color="0.6", lw=0.6)
    # Clip the y-range to a robust window so a few heavy-tailed features don't
    # squash everything near zero (the violins are simply cropped, not removed).
    hi = max(np.nanpercentile(xz, 99), np.nanmax(np.abs(obsz)) + 1)
    ax.set_ylim(-hi, hi)
    ax.set_xticks(np.arange(n_feat))
    ax.set_xticklabels(names, rotation=90)
    ax.set_ylabel("z-scored feature (vs simulated spread)")
    ax.set_title(f"Feature posterior-predictive check  ({x.shape[0]} valid sims)")
    handles = [plt.Line2D([], [], color=c, lw=6, alpha=0.35, label=r)
               for r, c in _ROI_COLORS.items()]
    handles.append(plt.Line2D([], [], color="k", marker="D", ls="", label="x_obs"))
    ax.legend(handles=handles, ncol=4, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


# ── figure 4: simulation diagnostics ─────────────────────────────────────────────────
def plot_sim_diagnostics(data, out_path):
    x = data["x"]
    x_obs = data["x_obs"]
    bad = data["bad_mask"]
    n_total = x.shape[0]
    n_bad = int(bad.sum())

    valid = ~bad
    xv = x[valid]
    thetav = data["theta"][valid]
    # z-scored Euclidean distance of each valid sim's features to the observation.
    sd = xv.std(axis=0)
    sd[sd == 0] = 1.0
    dist = np.sqrt((((xv - x_obs) / sd) ** 2).sum(axis=1))

    n_best = max(1, int(0.05 * len(dist)))        # top 5% closest sims
    best_idx = np.argsort(dist)[:n_best]

    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 3)

    # (a) failed-sim bar
    ax0 = fig.add_subplot(gs[0, 0])
    ax0.bar(["valid", "failed"], [n_total - n_bad, n_bad],
            color=["C0", "C3"])
    ax0.set_title(f"Simulations: {n_bad}/{n_total} failed "
                  f"({100 * n_bad / n_total:.1f}%)")
    ax0.set_ylabel("count")

    # (b) distance histogram
    ax1 = fig.add_subplot(gs[0, 1])
    ax1.hist(dist, bins=40, color="0.7")
    ax1.axvline(dist[best_idx].max(), color="C3", ls="--",
                label=f"top {n_best} cutoff")
    ax1.set_xlabel("z-scored distance to x_obs")
    ax1.set_ylabel("count")
    ax1.set_title("Feature distance to observation")
    ax1.legend()

    # (c) param-value distribution of best sims vs all (normalised to bounds)
    ax2 = fig.add_subplot(gs[0, 2])
    span = (BOUNDS[:, 1] - BOUNDS[:, 0]).copy()
    span[span == 0] = 1.0
    best_norm = (thetav[best_idx] - BOUNDS[:, 0]) / span
    ax2.boxplot([best_norm[:, j] for j in range(len(PARAM_NAMES))])
    ax2.set_xticks(np.arange(1, len(PARAM_NAMES) + 1))
    ax2.set_xticklabels(PARAM_NAMES)
    ax2.set_ylim(-0.05, 1.05)
    ax2.set_ylabel("value (normalised to bounds)")
    ax2.set_title(f"Best-matching sims (top {n_best})")
    ax2.tick_params(axis="x", rotation=90)

    # (d) theta scatter of best matches for the first 3 informative param pairs
    pairs = [(0, 1), (2, 3), (4, 5)]
    for k, (a, b) in enumerate(pairs):
        ax = fig.add_subplot(gs[1, k])
        ax.scatter(thetav[:, a], thetav[:, b], s=4, color="0.8", label="all")
        sc = ax.scatter(thetav[best_idx, a], thetav[best_idx, b], s=14,
                        c=dist[best_idx], cmap="viridis_r", label="best")
        ax.set_xlabel(PARAM_NAMES[a])
        ax.set_ylabel(PARAM_NAMES[b])
        ax.set_xlim(*BOUNDS[a])
        ax.set_ylim(*BOUNDS[b])
        if k == 0:
            ax.legend(loc="upper right")
    fig.colorbar(sc, ax=fig.axes[-1], label="distance", fraction=0.046)
    fig.suptitle("Simulation diagnostics", y=1.0)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


# ── figures 5+6: MAP posterior-predictive re-simulation vs data ───────────────────────
def plot_map_resim(data, plots_dir):
    """Re-run the model at the MAP parameters and compare TC + TF to measured data."""
    if not RESDIR:
        print("  [skip] RESDIR not set -> cannot load target data for MAP re-simulation.")
        return
    if not (os.path.exists(TC_DATA_PATH) and os.path.exists(TF_DATA_PATH)):
        print("  [skip] target CSVs not found -> skipping MAP re-simulation.")
        return

    from somato_model import SomatoModel

    map_est = np.array([data["summary"][n]["map"] for n in PARAM_NAMES])
    params = dict(zip(PARAM_NAMES, map_est))
    print("  MAP params:", {k: round(v, 4) for k, v in params.items()})

    model = SomatoModel(BASE_PARAMS)
    model.apply_params(params)
    model.initialize_state()
    model.simulate()

    err_tc, tc_sim, tc_target = model.compute_error_timecourse(TC_DATA_PATH)
    err_tf, tf_sim, tf_target = model.compute_error_timefreq(TF_DATA_PATH)

    # ── time course ────────────────────────────────────────────────────────────────
    # Sim is in arbitrary model-dipole units, measured is in nAm — different scales.
    # Peak-normalize each trace by its max |amplitude| over the analysis window
    # (-200 ms onward, idx 150), exactly as compute_error_timecourse does, so the
    # overlay reflects waveform shape/timing rather than the unit mismatch.
    n = np.asarray(tc_sim["A3b"]).shape[0]
    times_ms = np.linspace(-500, 400, n)
    eps = 1e-10
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), sharex=True)
    for ax, roi in zip(axes, ROIS):
        tgt = np.asarray(tc_target[roi])[:n]
        sim = np.asarray(tc_sim[roi])[:n]
        tgt = tgt / (np.max(np.abs(tgt[150:])) + eps)
        sim = sim / (np.max(np.abs(sim[150:])) + eps)
        ax.plot(times_ms, tgt, color="k", label="measured")
        ax.plot(times_ms, sim, color="C3", label="MAP sim")
        ax.axvline(0, color="0.6", lw=0.6)
        ax.set_title(roi)
        ax.set_xlabel("time (ms)")
    axes[0].set_ylabel("dipole (peak-normalized)")
    axes[0].legend(loc="best")
    fig.suptitle(f"MAP re-simulation — time course  (err={err_tc:.3f})", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "map_timecourse.png"))
    plt.close(fig)

    # ── time-frequency ────────────────────────────────────────────────────────────────
    # Both sim power and measured (beamformer) power are arbitrary units on different
    # scales, so show the quantity compute_error_timefreq actually compares: log10 of
    # each frequency's power normalized by its own baseline (idx 150:200 = -200..0 ms),
    # over the analysis window (idx 200:400 = 0..+400 ms). Diverging scale centered at 0.
    freqs = np.arange(1, 41, 1)
    baseline_slice = slice(150, 200)
    analysis_slice = slice(200, 400)
    eps = 1e-10
    times_full = np.linspace(-500, 400, 451)
    t_ms = times_full[analysis_slice]
    extent = [t_ms[0], t_ms[-1], freqs[0], freqs[-1]]

    def _norm_log(P):
        bl = P[:, baseline_slice].mean(axis=1, keepdims=True)
        return np.log10(P[:, analysis_slice] / (bl + eps) + eps)

    fig, axes = plt.subplots(3, 2, figsize=(10, 9), sharex=True, sharey=True)
    for r, roi in enumerate(ROIS):
        P_sim = np.asarray(tf_sim[roi])
        P_tgt = np.asarray(tf_target[roi])
        if P_sim.shape != P_tgt.shape:                 # match error code: trim sim to target
            P_sim = P_sim[:, :P_tgt.shape[1]]
        L_sim, L_tgt = _norm_log(P_sim), _norm_log(P_tgt)
        vmax = float(np.nanmax(np.abs([L_sim, L_tgt])))   # symmetric, shared per ROI
        for c, (L, label) in enumerate([(L_sim, "MAP sim"), (L_tgt, "measured")]):
            ax = axes[r, c]
            im = ax.imshow(L, aspect="auto", origin="lower", extent=extent,
                           vmin=-vmax, vmax=vmax, cmap="RdBu_r")
            ax.axvline(0, color="k", lw=0.5)
            if r == 0:
                ax.set_title(label)
            if c == 0:
                ax.set_ylabel(f"{roi}\nfreq (Hz)")
            if r == len(ROIS) - 1:
                ax.set_xlabel("time (ms)")
        fig.colorbar(im, ax=axes[r, :].tolist(), fraction=0.025, pad=0.02,
                     label="log10(power / baseline)")
    fig.suptitle(f"MAP re-simulation — Morlet TF, baseline-normalized  (err={err_tf:.3f})", y=1.0)
    fig.savefig(os.path.join(plots_dir, "map_timefreq.png"))
    plt.close(fig)


# ── shared posterior sampling (cached on the data dict) ───────────────────────────────
def _posterior_samples(data):
    if "_samples" not in data:
        x_obs = torch.as_tensor(data["x_obs"], dtype=torch.float32)
        data["_samples"] = data["posterior"].sample((N_SAMPLES,), x=x_obs).numpy()
    return data["_samples"]


# ── driver ────────────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description="Detailed plots for an SBI (NPE) run.")
    p.add_argument("run_dir", nargs="?", default=None,
                   help="sbi_<timestamp> run folder (default: latest under SBI_RESULTS_ROOT)")
    p.add_argument("--no-resim", action="store_true",
                   help="skip the MAP re-simulation TC/TF comparison (no model run / RESDIR)")
    args = p.parse_args()

    run_dir = args.run_dir or _latest_run_dir()
    run_dir = os.path.abspath(run_dir)
    plots_dir = os.path.join(run_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    print(f"Plotting SBI run: {run_dir}")

    _set_style()
    data = load_run(run_dir)

    written = []
    for fname, fn in [
        ("posterior_marginals.png", plot_marginals),
        ("posterior_corner.png",    plot_corner),
        ("feature_ppc.png",         plot_feature_ppc),
        ("sim_diagnostics.png",     plot_sim_diagnostics),
    ]:
        out = os.path.join(plots_dir, fname)
        fn(data, out)
        written.append(out)
        print(f"  wrote {fname}")

    if not args.no_resim:
        print("  MAP re-simulation ...")
        plot_map_resim(data, plots_dir)
        for f in ("map_timecourse.png", "map_timefreq.png"):
            if os.path.exists(os.path.join(plots_dir, f)):
                written.append(os.path.join(plots_dir, f))

    print("\nDone. Figures written to:")
    for w in written:
        print(f"  {w}")


if __name__ == "__main__":
    main()
