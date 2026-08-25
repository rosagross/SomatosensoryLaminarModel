"""
File: plotting_functions_parameter_space.py
Description:
    Figures for the parameter-space exploration produced by
    Simulations/run_parameter_sweep.py - where in parameter space the model
    oscillates, in which band, and how strongly.

    Driven by Analysis/step003_plot_parameter_space.py. Kept separate from
    plotting_functions_analysis.py because it consumes a different table: the sweep's
    tidy `sweep_features.csv` (one row per parameter set x signal), not the per-run
    `processed.csv` that step001 writes.

    The three sweep modes each get their own family of figures:
      grid    2-D heatmaps: regime, dominant band, peak frequency, alpha prominence
      line    one panel per parameter: how frequency / alpha prominence respond
      screen  marginals of alpha prominence against every parameter, plus a ranking
              of the strongest oscillatory parameter sets
    plus a spectral waterfall that works on any 1-D sweep.
"""

import json
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm

# ── palette ────────────────────────────────────────────────────────────────────
# Okabe-Ito steps, checked with the dataviz palette validator: every pair clears the
# colour-vision-deficiency separation floor and the normal-vision floor at --pairs all.
# The two "nothing is happening here" states are deliberately neutral rather than a
# hue, so a map that is mostly inert reads as mostly blank.
REGIME_COLORS = {
    "diverged": "#000000",      # the integration blew up
    "pinned": "#D9D9D6",        # firing rate frozen at a rail of the sigmoid
    "fixed_point": "#F2F2EF",   # settled, no fluctuation
    "damped": "#E69F00",        # still ringing down
    "noise_driven": "#56B4E9",  # fluctuating, spectrally featureless
    "oscillation": "#D55E00",   # sustained, with a spectral peak
}
# Ordered from least to most alive, so the colour bar reads as a progression.
REGIME_ORDER = ["diverged", "pinned", "fixed_point", "damped", "noise_driven",
                "oscillation"]

BAND_COLORS = {
    "none": "#F2F2EF",
    "theta": "#56B4E9",
    "alpha": "#D55E00",      # the band of interest gets the most salient hue
    "beta": "#009E73",
    "gamma": "#7A5195",
}
BAND_ORDER = ["none", "theta", "alpha", "beta", "gamma"]

# Frequency extent of each band, for shading the frequency axis (matches BANDS in
# Simulations/model/oscillation_metrics.py).
_BAND_SPANS = {"theta": (4.0, 7.0), "alpha": (8.0, 13.0),
               "beta": (17.5, 30.0), "gamma": (31.0, 50.0)}

# The signals plotted by default: the three dipole ROIs, which are what the measured
# EEG is compared against, plus the excitatory population of each S1 layer.
DEFAULT_SIGNALS = ["roi_A3b", "roi_A1", "roi_S2"]
SIGNAL_LABELS = {"roi_A3b": "A3b (dipole)", "roi_A1": "S1 (dipole)",
                 "roi_S2": "S2 (dipole)"}

# Layer x cell-type layout per area, as in plotting_functions_analysis
AREA_GRIDS = {
    "A3b": ([["E3b", "PV3b", "SST3b", "VIP3b"]], ["A3b"]),
    "S1": ([["E1", "PV1", "SST1", "VIP1"],
            ["E2", "PV2", "SST2", None],
            ["E3", "PV3", "SST3", None],
            ["E4", "PV4", "SST4", None]],
           ["Layer 2/3", "Layer 4", "Layer 5", "Layer 6"]),
    "S2": ([["E1S2", "PV1S2", "SST1S2", "VIP1S2"],
            ["E2S2", "PV2S2", "SST2S2", None],
            ["E3S2", "PV3S2", "SST3S2", None],
            ["E4S2", "PV4S2", "SST4S2", None]],
           ["Layer 2/3", "Layer 4", "Layer 5", "Layer 6"]),
}
CELLTYPE_LABELS = ["Excitatory", "PV", "SST", "VIP"]

# Readable axis labels for the swept parameters
PARAM_LABELS = {
    "coupling_strength": r"global coupling $g$",
    "strength_I": r"E/I balance $s_I$",
    "Ib_strength": "background drive",
    "g_intercortical": "inter-areal gain",
    "g_thal": "VPM gain",
    "sI_thal": "reticular gain ratio",
    "g_thalPOm": "POm gain",
    "thal_EtoE": "VPM→VPM",
    "thal_EtoI": "VPM→reticular",
    "thal_ItoE": "reticular→VPM",
    "thal_ItoI": "reticular→reticular",
    "thal_POmtoI": "POm→reticular",
    "e3b_tau": r"A3b $\tau_E$ (ms)",
    "e1_tau": r"S1 $\tau_E$ (ms)",
    "e2_tau": r"S2 $\tau_E$ (ms)",
    "delay_factor": "long-range delay (s)",
    "delay_factor_short": "A3b↔S1 delay (s)",
    "thal_delay_factor": "thalamic delay (s)",
    "p_2PVE": "L4 PV←E prob. (%)",
    "p_4PVE": "L6 PV←E prob. (%)",
}


def param_label(name):
    return PARAM_LABELS.get(name, name)


def signal_label(name):
    return SIGNAL_LABELS.get(name, name)


def set_style():
    """Repo figure style, but safe on a headless node.

    plotting_style.figure_style() sizes figures from the screen width via tkinter and
    raises without a display, which is how these sweeps are usually run - so fall back
    to the same rcParams that Simulations/plot_sbi_results.py uses there.
    """
    try:
        from plotting_style import figure_style
        colors, _ = figure_style()
        return colors
    except Exception:
        sns.set_theme(style="ticks", context="paper")
        plt.rcParams.update({"font.size": 9, "axes.titlesize": 10,
                             "axes.labelsize": 9, "savefig.bbox": "tight",
                             "figure.dpi": 110, "pdf.fonttype": 42,
                             "ps.fonttype": 42})
        return {}


def save_figure(fig, figure_dir, name):
    """Save as PDF + 300 dpi PNG under `name`, the repo's convention."""
    os.makedirs(figure_dir, exist_ok=True)
    fig.savefig(os.path.join(figure_dir, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(figure_dir, name + ".png"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    return os.path.join(figure_dir, name + ".png")


# ── loading ────────────────────────────────────────────────────────────────────
def load_sweep(sweep_dir):
    """(features DataFrame, config dict) of one sweep directory."""
    df = pd.read_csv(os.path.join(sweep_dir, "sweep_features.csv"))
    with open(os.path.join(sweep_dir, "sweep_config.json")) as f:
        cfg = json.load(f)
    return df, cfg


def load_spectra(sweep_dir):
    """The stored spectra of one sweep: freqs, psd, psd_null, signal names, theta."""
    path = os.path.join(sweep_dir, "sweep_spectra.hdf5")
    with h5py.File(path, "r") as f:
        return {
            "freqs": f["freqs"][:],
            "psd": f["psd"][:],
            "psd_null": f["psd_null"][:],
            "signals": [s.decode() for s in f["signals"][:]],
            "param_names": [s.decode() for s in f["param_names"][:]],
            "theta": f["theta"][:],
        }


def gain_spectra(spec):
    """psd / psd_null, the network gain spectrum stored per simulation and signal."""
    null = spec["psd_null"]
    with np.errstate(divide="ignore", invalid="ignore"):
        gain = spec["psd"] / np.where(null > 0, null, np.nan)
    return gain


# ── shared heatmap helpers ─────────────────────────────────────────────────────
def _pivot(df, x, y, value, xs=None, ys=None):
    """Wide matrix of `value` over the (x, y) parameter grid, rows = y.

    `xs` / `ys` force the full grid axes. Without them a row or column that is entirely
    masked out - a band of the grid where nothing oscillates - would be dropped, and
    panels of the same figure would end up spanning different parameter ranges, which
    makes them impossible to compare by eye.
    """
    sub = df[[x, y, value]].dropna(subset=[x, y])
    if sub.empty:
        return pd.DataFrame()
    matrix = sub.pivot_table(index=y, columns=x, values=value, aggfunc="first",
                             dropna=False)
    if xs is not None:
        matrix = matrix.reindex(columns=xs)
    if ys is not None:
        matrix = matrix.reindex(index=ys)
    return matrix


def _grid_axes(df, x, y):
    """The full sorted set of values each grid axis takes across the whole sweep."""
    return sorted(df[x].dropna().unique()), sorted(df[y].dropna().unique())


def _tick_labels(values, max_ticks=8):
    """Show at most `max_ticks` of the axis values, blanking the rest."""
    values = np.asarray(values, dtype=float)
    step = max(int(np.ceil(len(values) / max_ticks)), 1)
    return [f"{v:g}" if i % step == 0 else "" for i, v in enumerate(values)]


def _categorical_heatmap(ax, matrix, order, colors, show_cbar=False, cbar_ax=None):
    """Draw a heatmap of category labels using a fixed category order and palette."""
    codes = matrix.replace({name: i for i, name in enumerate(order)})
    codes = codes.apply(pd.to_numeric, errors="coerce")
    cmap = ListedColormap([colors[name] for name in order])
    norm = BoundaryNorm(np.arange(-0.5, len(order) + 0.5), ncolors=cmap.N)
    hm = sns.heatmap(codes, cmap=cmap, norm=norm, ax=ax, cbar=show_cbar,
                     cbar_ax=cbar_ax,
                     cbar_kws={"ticks": np.arange(len(order))} if show_cbar else None)
    if show_cbar:
        cbar = hm.collections[0].colorbar
        cbar.set_ticklabels(order)
        cbar.ax.tick_params(length=0)
    return hm


def _format_grid_axes(ax, matrix, x, y, show_x, show_y):
    ax.invert_yaxis()
    ax.tick_params(axis="both", length=0)
    ax.set_xticks(np.arange(len(matrix.columns)) + 0.5)
    ax.set_yticks(np.arange(len(matrix.index)) + 0.5)
    ax.set_xticklabels(_tick_labels(matrix.columns), rotation=0)
    ax.set_yticklabels(_tick_labels(matrix.index), rotation=0)
    ax.set_xlabel(param_label(x) if show_x else "")
    ax.set_ylabel(param_label(y) if show_y else "")


# ── grid-mode figures ──────────────────────────────────────────────────────────
def plot_regime_map(df, x, y, figure_dir, name, signals=None, kind="regime"):
    """One categorical map per signal: which dynamical regime each parameter set is in.

    `kind` is 'regime' (dead / pinned / ringing / noisy / oscillating) or 'band' (which
    band the dominant peak falls in). This is the figure to read first: a continuous
    frequency map is meaningless wherever the regime is not 'oscillation', so the
    frequency figures below mask exactly the cells this one marks as inert.
    """
    signals = signals or DEFAULT_SIGNALS
    order, colors = ((REGIME_ORDER, REGIME_COLORS) if kind == "regime"
                     else (BAND_ORDER, BAND_COLORS))
    xs, ys = _grid_axes(df, x, y)
    fig, axes = plt.subplots(1, len(signals), figsize=(3.1 * len(signals) + 1.4, 3.0),
                             squeeze=False)
    axes = axes[0]
    cbar_ax = fig.add_axes([0.92, 0.18, 0.015, 0.64])
    for i, (ax, sig) in enumerate(zip(axes, signals)):
        matrix = _pivot(df[df["signal"] == sig], x, y, kind, xs, ys)
        if matrix.empty:
            ax.axis("off")
            continue
        _categorical_heatmap(ax, matrix, order, colors,
                             show_cbar=(i == len(signals) - 1), cbar_ax=cbar_ax)
        _format_grid_axes(ax, matrix, x, y, True, i == 0)
        ax.set_title(signal_label(sig))
    fig.suptitle(f"{'Dynamical regime' if kind == 'regime' else 'Dominant band'}: "
                 f"{param_label(x)} vs {param_label(y)}", fontweight="bold")
    fig.tight_layout(rect=[0, 0, 0.9, 0.97])
    return save_figure(fig, figure_dir, f"{name}_{kind}map")


def plot_value_map(df, x, y, value, figure_dir, name, signals=None,
                   mask_non_oscillating=True, cmap=None, vmin=None, vmax=None,
                   label=None, diverging=False):
    """One continuous map per signal (peak frequency, alpha prominence, amplitude...).

    Cells whose regime is not 'oscillation' are blanked by default: a peak frequency
    read off a pinned or settled network is the argmax of numerical residue, and
    plotting it would fill the map with structure that is not there.
    """
    signals = signals or DEFAULT_SIGNALS
    label = label or value
    plot_df = df.copy()
    if mask_non_oscillating:
        plot_df.loc[plot_df["regime"] != "oscillation", value] = np.nan

    xs, ys = _grid_axes(df, x, y)
    mats = {sig: _pivot(plot_df[plot_df["signal"] == sig], x, y, value, xs, ys)
            for sig in signals}
    finite = np.concatenate([m.to_numpy(dtype=float).ravel() for m in mats.values()
                             if not m.empty] or [np.array([np.nan])])
    finite = finite[np.isfinite(finite)]
    if vmin is None or vmax is None:
        if finite.size == 0:
            vmin, vmax = 0.0, 1.0
        elif diverging:
            lim = float(np.nanmax(np.abs(finite))) or 1.0
            vmin, vmax = -lim, lim
        else:
            vmin, vmax = float(np.nanmin(finite)), float(np.nanmax(finite))
    if cmap is None:
        cmap = "RdBu_r" if diverging else "magma"
    norm = (TwoSlopeNorm(vcenter=0.0, vmin=min(vmin, -1e-9), vmax=max(vmax, 1e-9))
            if diverging else plt.Normalize(vmin=vmin, vmax=vmax))

    fig, axes = plt.subplots(1, len(signals), figsize=(3.1 * len(signals) + 1.4, 3.0),
                             squeeze=False)
    axes = axes[0]
    for i, (ax, sig) in enumerate(zip(axes, signals)):
        matrix = mats[sig]
        if matrix.empty or not np.isfinite(matrix.to_numpy(dtype=float)).any():
            ax.axis("off")
            ax.text(0.5, 0.5, f"{signal_label(sig)}\n(nothing oscillating)",
                    ha="center", va="center", transform=ax.transAxes, fontsize=8)
            continue
        sns.heatmap(matrix, cmap=cmap, norm=norm, ax=ax, cbar=False)
        _format_grid_axes(ax, matrix, x, y, True, i == 0)
        ax.set_title(signal_label(sig))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.02)
    cbar.set_label(label)
    fig.suptitle(f"{label}: {param_label(x)} vs {param_label(y)}", fontweight="bold",
                 y=1.02)
    return save_figure(fig, figure_dir, f"{name}_{value}")


def plot_population_grid_map(df, x, y, value, area, figure_dir, name,
                             mask_non_oscillating=True, cmap="magma",
                             vmin=None, vmax=None, label=None, categorical=None):
    """Layer x cell-type panel grid of one area, same layout as the step002 heatmaps.

    Shows whether an oscillation is a property of the whole area or of particular
    layers / interneuron classes - the laminar detail the three ROI dipoles average
    away.
    """
    layers, row_labels = AREA_GRIDS[area]
    label = label or value
    xs, ys = _grid_axes(df, x, y)
    plot_df = df.copy()
    if mask_non_oscillating and categorical is None:
        plot_df.loc[plot_df["regime"] != "oscillation", value] = np.nan

    nrows, ncols = len(layers), 4
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.4 * ncols, 2.2 * nrows),
                             squeeze=False)
    if categorical is None:
        vals = plot_df[plot_df["signal"].isin(
            [p for row in layers for p in row if p])][value].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        vmin = vmin if vmin is not None else (float(vals.min()) if vals.size else 0.0)
        vmax = vmax if vmax is not None else (float(vals.max()) if vals.size else 1.0)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        cbar_ax = None
    else:
        order, colors = ((REGIME_ORDER, REGIME_COLORS) if categorical == "regime"
                         else (BAND_ORDER, BAND_COLORS))
        cbar_ax = fig.add_axes([0.94, 0.25, 0.012, 0.5])

    for r, row in enumerate(layers):
        for c in range(ncols):
            ax = axes[r][c]
            pop = row[c] if c < len(row) else None
            if pop is None:
                ax.axis("off")
                continue
            matrix = _pivot(plot_df[plot_df["signal"] == pop], x, y, value, xs, ys)
            if matrix.empty or (categorical is None and
                                not np.isfinite(matrix.to_numpy(dtype=float)).any()):
                ax.axis("off")
                ax.text(0.5, 0.5, f"{pop}\n(none)", ha="center", va="center",
                        transform=ax.transAxes, fontsize=7)
                continue
            if categorical is None:
                sns.heatmap(matrix, cmap=cmap, norm=norm, ax=ax, cbar=False)
            else:
                _categorical_heatmap(ax, matrix, order, colors,
                                     show_cbar=(r == 0 and c == 0), cbar_ax=cbar_ax)
            _format_grid_axes(ax, matrix, x, y, r == nrows - 1, c == 0)
            if r == 0:
                ax.set_title(CELLTYPE_LABELS[c], fontweight="bold")
            if c == 0:
                ax.text(-0.42, 0.5, row_labels[r], transform=ax.transAxes,
                        rotation=90, va="center", ha="center", fontweight="bold")

    if categorical is None:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.02)
        cbar.set_label(label)
    fig.suptitle(f"{label} - {area}: {param_label(x)} vs {param_label(y)}",
                 fontweight="bold", y=1.0 + 0.06 / nrows)
    return save_figure(fig, figure_dir, f"{name}_{area}_{value}")


# ── line-mode figures ──────────────────────────────────────────────────────────
def plot_line_sensitivity(df, value, figure_dir, name, signals=None, label=None,
                          mask_non_oscillating=False, reference=None):
    """One panel per parameter: `value` against that parameter, one line per signal.

    All other parameters stay at the reference point, so each panel is a
    one-at-a-time sensitivity curve - what changing this parameter alone does to the
    oscillation. Points whose regime is not 'oscillation' are drawn hollow, so a curve
    running through an inert region cannot be mistaken for a moving rhythm.
    """
    signals = signals or DEFAULT_SIGNALS
    label = label or value
    params = sorted(df["varied_param"].dropna().unique())
    ncols = 4
    nrows = int(np.ceil(len(params) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.1 * ncols, 2.5 * nrows),
                             squeeze=False, sharey=True)
    colors = sns.color_palette("colorblind", len(signals))

    for i, param in enumerate(params):
        ax = axes[i // ncols][i % ncols]
        for sig, color in zip(signals, colors):
            sub = df[(df["varied_param"] == param) & (df["signal"] == sig)]
            sub = sub.sort_values("varied_value")
            if sub.empty:
                continue
            vals = sub[value].to_numpy(dtype=float).copy()
            osc = (sub["regime"] == "oscillation").to_numpy()
            if mask_non_oscillating:
                vals[~osc] = np.nan
            ax.plot(sub["varied_value"], vals, color=color, lw=1.2,
                    label=signal_label(sig) if i == 0 else None)
            ax.scatter(sub["varied_value"][osc], vals[osc], s=9, color=color, zorder=3)
            ax.scatter(sub["varied_value"][~osc], vals[~osc], s=9,
                       facecolors="none", edgecolors=color, linewidths=0.6, zorder=3)
        if reference is not None and param in reference:
            ax.axvline(reference[param], color="0.6", lw=0.8, ls=":", zorder=0)
        ax.set_xlabel(param_label(param))
        if i % ncols == 0:
            ax.set_ylabel(label)
    for j in range(len(params), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(signals),
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(f"{label} against each parameter (all others at the reference point; "
                 f"hollow = not oscillating)", fontweight="bold", y=1.0)
    sns.despine(fig=fig)
    fig.tight_layout()
    return save_figure(fig, figure_dir, f"{name}_line_{value}")


def plot_line_regimes(df, figure_dir, name, signal="roi_A1"):
    """A regime strip per parameter: which regime the model is in along each line."""
    params = sorted(df["varied_param"].dropna().unique())
    sub_all = df[df["signal"] == signal]
    # every parameter is swept over the same number of points, but over its own range,
    # so the strip is indexed by position along the sweep rather than by value
    n_points = int(sub_all.groupby("varied_param")["varied_value"].size().max())
    fig, ax = plt.subplots(figsize=(8, 0.32 * len(params) + 1.6))
    codes = np.full((len(params), n_points), np.nan)
    for i, param in enumerate(params):
        sub = sub_all[sub_all["varied_param"] == param].sort_values("varied_value")
        codes[i, :len(sub)] = [REGIME_ORDER.index(r) if r in REGIME_ORDER else np.nan
                               for r in sub["regime"]]
    cmap = ListedColormap([REGIME_COLORS[r] for r in REGIME_ORDER])
    norm = BoundaryNorm(np.arange(-0.5, len(REGIME_ORDER) + 0.5), ncolors=cmap.N)
    im = ax.imshow(codes, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
    ax.set_yticks(range(len(params)))
    # each row spans that parameter's own range, so the range belongs in its label
    ranges = sub_all.groupby("varied_param")["varied_value"].agg(["min", "max"])
    ax.set_yticklabels([f"{param_label(p)}  [{ranges.loc[p, 'min']:g}–"
                        f"{ranges.loc[p, 'max']:g}]" for p in params])
    ax.set_xlabel("parameter value (low → high, across its sweep range)")
    ax.set_xticks([])
    cbar = fig.colorbar(im, ax=ax, ticks=range(len(REGIME_ORDER)), fraction=0.04)
    cbar.set_ticklabels(REGIME_ORDER)
    cbar.ax.tick_params(length=0)
    ax.set_title(f"Regime along each parameter - {signal_label(signal)}",
                 fontweight="bold")
    fig.tight_layout()
    return save_figure(fig, figure_dir, f"{name}_line_regimes_{signal}")


# ── screen-mode figures ────────────────────────────────────────────────────────
def plot_screen_marginals(df, value, figure_dir, name, signal="roi_A1", label=None,
                          n_bins=12, min_per_bin=3):
    """`value` against every parameter of a Latin-hypercube screen, one panel each.

    Every other parameter varies freely between the points of a panel, so the scatter
    shows the *marginal* effect and the binned median the trend through it. A
    parameter that matters on its own tilts the median line; one that only matters in
    combination shows a wide scatter with a flat median - which is itself the answer.
    """
    label = label or value
    sub = df[df["signal"] == signal].copy()
    params = [c for c in sub.columns if c in PARAM_LABELS]
    ncols = 4
    nrows = int(np.ceil(len(params) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.1 * ncols, 2.4 * nrows),
                             squeeze=False, sharey=True)
    osc = sub["regime"] == "oscillation"
    for i, param in enumerate(params):
        ax = axes[i // ncols][i % ncols]
        ax.scatter(sub.loc[~osc, param], sub.loc[~osc, value], s=5, color="0.75",
                   linewidths=0, label="not oscillating" if i == 0 else None)
        ax.scatter(sub.loc[osc, param], sub.loc[osc, value], s=6,
                   color=REGIME_COLORS["oscillation"], linewidths=0,
                   label="oscillating" if i == 0 else None)
        edges = np.linspace(sub[param].min(), sub[param].max(), n_bins + 1)
        centres = 0.5 * (edges[:-1] + edges[1:])
        # The trend is taken over the oscillating sets only: "how prominent is alpha"
        # is a question about parameter sets that have a rhythm at all, and including
        # the pinned majority would drag the line towards a value that describes
        # nothing. A bin holding one or two points has a median that is just that
        # point, so those are left as gaps rather than drawn as a spike.
        osc_sub = sub[osc]
        med = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            vals = osc_sub.loc[(osc_sub[param] >= lo) & (osc_sub[param] < hi),
                               value].dropna()
            med.append(vals.median() if len(vals) >= min_per_bin else np.nan)
        ax.plot(centres, med, color="#0072B2", lw=1.6, zorder=4,
                label="binned median (oscillating)" if i == 0 else None)
        ax.set_xlabel(param_label(param))
        if i % ncols == 0:
            ax.set_ylabel(label)
    # The tail of near-unstable parameter sets runs to gains of 10^4, which would
    # compress the range everything else lives in to a couple of pixels. Clipping the
    # view keeps the bulk readable; the extreme sets are still listed in full by
    # `top_oscillatory_points`.
    finite = sub[value].to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    clipped = 0
    if finite.size:
        lo, hi = np.percentile(finite, [1, 95])
        pad = 0.08 * (hi - lo) if hi > lo else 0.1
        lo, hi = lo - pad, hi + pad
        clipped = int(((finite < lo) | (finite > hi)).sum())
        axes[0][0].set_ylim(lo, hi)
    for j in range(len(params), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    note = f", {clipped} outside the y-range" if clipped else ""
    fig.suptitle(f"{label} against each parameter - {signal_label(signal)} "
                 f"({len(sub)} random parameter sets{note})", fontweight="bold", y=1.0)
    sns.despine(fig=fig)
    fig.tight_layout()
    return save_figure(fig, figure_dir, f"{name}_screen_{value}_{signal}")


def plot_screen_summary(df, figure_dir, name, signals=None):
    """How much of the space is in each regime, and which bands the peaks land in."""
    signals = signals or DEFAULT_SIGNALS
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))

    frac = (df[df["signal"].isin(signals)]
            .groupby("signal")["regime"].value_counts(normalize=True)
            .unstack(fill_value=0).reindex(columns=REGIME_ORDER, fill_value=0)
            .reindex(signals))
    bottom = np.zeros(len(frac))
    for regime in REGIME_ORDER:
        ax = axes[0]
        ax.bar([signal_label(s) for s in frac.index], frac[regime], bottom=bottom,
               color=REGIME_COLORS[regime], label=regime, width=0.6,
               edgecolor="white", linewidth=0.6)
        bottom += frac[regime].to_numpy()
    axes[0].set_ylabel("fraction of parameter sets")
    axes[0].set_title("Regimes across the sampled space")
    # the bars reach 1.0, so the legend needs its own headroom rather than overlapping
    axes[0].set_ylim(0, 1.28)
    axes[0].legend(frameon=False, fontsize=7, ncol=3, loc="upper center")

    bands = (df[(df["signal"].isin(signals)) & (df["regime"] == "oscillation")]
             .groupby("signal")["band"].value_counts(normalize=True)
             .unstack(fill_value=0).reindex(columns=BAND_ORDER, fill_value=0)
             .reindex(signals).fillna(0))
    bottom = np.zeros(len(bands))
    for band in BAND_ORDER:
        axes[1].bar([signal_label(s) for s in bands.index], bands[band], bottom=bottom,
                    color=BAND_COLORS[band], label=band, width=0.6,
                    edgecolor="white", linewidth=0.6)
        bottom += bands[band].to_numpy()
    axes[1].set_ylabel("fraction of oscillating sets")
    axes[1].set_title("Band of the dominant peak")
    axes[1].set_ylim(0, 1.28)
    axes[1].legend(frameon=False, fontsize=7, ncol=3, loc="upper center")

    sns.despine(fig=fig)
    fig.tight_layout()
    return save_figure(fig, figure_dir, f"{name}_screen_summary")


def top_oscillatory_points(df, signal="roi_A1", value="alpha_prom", n=20):
    """The `n` parameter sets with the strongest peak, as a table to inspect or re-run."""
    sub = df[(df["signal"] == signal) & (df["regime"] == "oscillation")]
    cols = [c for c in sub.columns if c in PARAM_LABELS]
    return (sub.sort_values(value, ascending=False)
            .head(n)[cols + [value, "peak_freq", "band", "amplitude"]]
            .reset_index(drop=True))


# ── spectra ────────────────────────────────────────────────────────────────────
def plot_gain_waterfall(spec, param, figure_dir, name, signal="roi_A1", fmax=45.0):
    """Network gain spectra along a 1-D sweep, stacked into one image.

    x = frequency, y = the swept parameter, colour = log10 gain over the loop-free
    control. A rhythm appears as a bright vertical ridge, and whether it *moves* with
    the parameter - a ridge that bends - is exactly the question a peak-frequency
    heatmap answers one number at a time.
    """
    sig_idx = spec["signals"].index(signal)
    par_idx = spec["param_names"].index(param)
    order = np.argsort(spec["theta"][:, par_idx])
    values = spec["theta"][order, par_idx]
    gain = gain_spectra(spec)[order, sig_idx, :]
    fmask = spec["freqs"] <= fmax
    with np.errstate(divide="ignore", invalid="ignore"):
        img = np.log10(gain[:, fmask])

    fig, ax = plt.subplots(figsize=(6, 4))
    lim = np.nanpercentile(np.abs(img), 98) if np.isfinite(img).any() else 1.0
    mesh = ax.pcolormesh(spec["freqs"][fmask], values, img, cmap="RdBu_r",
                         vmin=-lim, vmax=lim, shading="nearest")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel(param_label(param))
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("log10 gain over loop-free control")
    ax.set_title(f"Network gain spectrum - {signal_label(signal)}", fontweight="bold")
    fig.tight_layout()
    return save_figure(fig, figure_dir, f"{name}_gainspectra_{param}_{signal}")


# ── parameter-influence figures (seed-averaged screens) ────────────────────────
def plot_param_vs_frequency(df, param, figure_dir, name, signals=None, suffix='',
                            agreement_min=0.5, n_bins=12, min_per_bin=4):
    """Peak frequency against one parameter, one panel per ROI.

    The figure the `e1_tau` question needs: it puts the swept parameter on x, the
    dominant peak frequency on y, and shades the four frequency bands across the panel
    so a point's band can be read off its height. Colour repeats the band, so the
    scatter and the shading agree by construction.

    Points whose `band_agreement` is below `agreement_min` - i.e. the individual noise
    realisations disagreed with the seed-averaged spectrum about where the peak sits -
    are drawn hollow. That distinction is the whole reason this sweep exists: an earlier
    single-realisation screen reported alpha at points where the peak frequency in fact
    scattered from 4 to 41 Hz across seeds, and a solid-looking trend line through such
    points is worse than no plot at all.

    Parameters
    ----------
    df : DataFrame
        A screen's sweep_features.csv, with the consistency columns.
    param : str
        Parameter for the x axis, e.g. 'e1_tau'.
    signals : sequence of str
        ROIs to show, one panel each (default the three dipole ROIs).
    agreement_min : float
        band_agreement below which a point is drawn hollow.
    """
    signals = signals or DEFAULT_SIGNALS
    fig, axes = plt.subplots(1, len(signals), figsize=(4.0 * len(signals) + 1.6, 3.6),
                             sharey=True, squeeze=False)
    axes = axes[0]
    has_agreement = 'band_agreement' in df.columns

    for i, (ax, sig) in enumerate(zip(axes, signals)):
        sub = df[(df['signal'] == sig) & (df['regime'] == 'oscillation')].copy()
        # band shading first, so the points sit on top of it
        for bname, (lo, hi) in _BAND_SPANS.items():
            ax.axhspan(lo, hi, color=BAND_COLORS[bname], alpha=0.16, linewidth=0,
                       zorder=0)
            if i == 0:
                ax.text(0.012, (lo + hi) / 2, bname, transform=ax.get_yaxis_transform(),
                        va='center', ha='left', fontsize=7, color='0.35', zorder=1)
        if sub.empty:
            ax.text(0.5, 0.5, f'{signal_label(sig)}\n(nothing oscillating)',
                    ha='center', va='center', transform=ax.transAxes, fontsize=8)
            continue

        firm = (sub['band_agreement'] >= agreement_min) if has_agreement \
            else np.ones(len(sub), dtype=bool)
        for bname in BAND_ORDER:
            m = sub['band'] == bname
            if not m.any():
                continue
            ax.scatter(sub.loc[m & firm, param], sub.loc[m & firm, 'peak_freq'],
                       s=9, color=BAND_COLORS[bname], linewidths=0, zorder=3,
                       label=bname if i == len(signals) - 1 else None)
            ax.scatter(sub.loc[m & ~firm, param], sub.loc[m & ~firm, 'peak_freq'],
                       s=9, facecolors='none', edgecolors=BAND_COLORS[bname],
                       linewidths=0.5, zorder=2)

        # binned median and interquartile ribbon over the reliable points only
        rel = sub[firm]
        if len(rel) >= min_per_bin:
            edges = np.linspace(sub[param].min(), sub[param].max(), n_bins + 1)
            centres, med, q1, q3 = [], [], [], []
            for lo, hi in zip(edges[:-1], edges[1:]):
                v = rel.loc[rel[param].between(lo, hi), 'peak_freq'].dropna()
                if len(v) >= min_per_bin:
                    centres.append(0.5 * (lo + hi))
                    med.append(v.median()); q1.append(v.quantile(.25))
                    q3.append(v.quantile(.75))
            if centres:
                ax.fill_between(centres, q1, q3, color='0.35', alpha=0.18, linewidth=0,
                                zorder=4)
                ax.plot(centres, med, color='0.15', linewidth=1.8, zorder=5)
        ax.set_xlabel(param_label(param))
        ax.set_title(signal_label(sig))
        if i == 0:
            ax.set_ylabel('peak frequency (Hz)')

    handles, labs = axes[-1].get_legend_handles_labels()
    if handles:
        axes[-1].legend(handles, labs, title='dominant band', frameon=False, fontsize=7,
                        title_fontsize=8, loc='center left', bbox_to_anchor=(1.02, 0.5))
    note = (f'hollow = band_agreement < {agreement_min:g} (the noise realisations '
            f'disagreed); line = median, band = IQR')
    fig.suptitle(f'Influence of {param_label(param)} on the oscillation frequency',
                 fontweight='bold')
    fig.text(0.5, -0.04, note, ha='center', fontsize=7.5, color='0.35')
    sns.despine(fig=fig, trim=True)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return save_figure(fig, figure_dir, f'{name}_freqVs_{param}{suffix}')


def plot_param_influence(df, figure_dir, name, signals=None, suffix='',
                         values=(('peak_freq', 'peak frequency'),
                                 ('alpha_prom', 'alpha prominence'))):
    """Rank every swept parameter by how strongly it moves each measure.

    Spearman rho (monotone, so it does not assume the relationship is linear) between
    each parameter and each measure, over the oscillating points, one bar per ROI. The
    bars are sorted by the largest |rho| across ROIs, so the parameters that matter come
    out on top; the sign is kept on the bar so the direction is not lost.
    """
    from scipy.stats import spearmanr
    signals = signals or DEFAULT_SIGNALS
    params = [c for c in df.columns if c in PARAM_LABELS and df[c].nunique() > 1]

    fig, axes = plt.subplots(1, len(values), figsize=(5.2 * len(values), 
                                                      0.34 * len(params) + 2.0),
                             sharey=True, squeeze=False)
    axes = axes[0]
    colors = sns.color_palette('colorblind', len(signals))

    rho = {}
    for value, _ in values:
        for sig in signals:
            sub = df[(df['signal'] == sig) & (df['regime'] == 'oscillation')]
            for p in params:
                v = sub[[p, value]].dropna()
                rho[(value, sig, p)] = (spearmanr(v[p], v[value])[0]
                                        if len(v) > 10 else np.nan)
    order = sorted(params, key=lambda p: max(
        abs(rho.get((v, s, p), 0)) if np.isfinite(rho.get((v, s, p), np.nan)) else 0
        for v, _ in values for s in signals))
    y = np.arange(len(order))
    h = 0.8 / len(signals)

    for ax, (value, label) in zip(axes, values):
        for k, (sig, col) in enumerate(zip(signals, colors)):
            vals = [rho.get((value, sig, p), np.nan) for p in order]
            ax.barh(y + (k - (len(signals) - 1) / 2) * h, vals, height=h, color=col,
                    label=signal_label(sig) if ax is axes[0] else None)
        ax.axvline(0, color='0.5', linewidth=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([param_label(p) for p in order])
        ax.set_xlabel(r'Spearman $\rho$')
        ax.set_title(label)
    handles, labs = axes[0].get_legend_handles_labels()
    fig.legend(handles, labs, frameon=False, fontsize=8, loc='lower center',
               ncol=len(signals), bbox_to_anchor=(0.5, -0.02))
    fig.suptitle('Which parameters move the oscillation, ranked', fontweight='bold')
    sns.despine(fig=fig, trim=True)
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    return save_figure(fig, figure_dir, f'{name}_paramInfluence{suffix}')


# ── synaptic time-constant figures ─────────────────────────────────────────────
# Distinct qualitative colours (Paul Tol's bright set): the three curves have to be
# told apart at a glance, which a sequential colormap does not do.
TAU_COLORS = {"e1_tau": "#4477AA", "e3b_tau": "#EE6677", "e2_tau": "#228833"}


def _tau_line(df, param, signal):
    """The rows of a `line` sweep for one swept parameter and one signal, sorted by x.

    Sorting matters twice over: the x axis has to be monotone, and `sim_id` is what
    indexes the spectra arrays, so the two must be put in the same order or the ridge
    ends up scrambled against its own peak markers.
    """
    sub = df[(df["varied_param"] == param) & (df["signal"] == signal)]
    return sub.sort_values("varied_value")


def plot_tau_ridge_grid(df, spec, params, figure_dir, name, signals=None, fmax=45.0,
                        agreement_min=0.5, suffix=""):
    """Where the oscillation sits in frequency as each time constant is swept.

    One panel per (time constant x ROI). The image is the seed-averaged network gain
    spectrum - the run divided by the same parameter set with its loop gains scaled to
    1e-3 - so a resonance is a bright horizontal ridge and the question "does the
    frequency move with tau?" is answered by whether that ridge *bends*.

    Each column of the image (one tau value) has its own median over frequency
    subtracted, so what is coloured is the spectrum's *shape* and not its level. Without
    that the figure is unreadable: the ROI dipoles carry a roughly constant 1e6 gain
    over the loop-free control - the control's dipole is tiny - so the panel spans 5.9
    to 7.0 in log10 and a scale wide enough to hold it flattens the ~0.5 decades of
    actual spectral structure to one colour. The overall level is what the companion
    amplitude figure reports; this one is about where in frequency the power sits.

    Drawn on top is the scored `peak_freq`, and the two are deliberately different
    measurements. The image's brightest row is the raw argmax, which at low loop gain
    sits on the 1/f corner at 1-4 Hz and says nothing; `peak_freq` is the flank-
    referenced band prominence, which is the number worth quoting (it reports 27 Hz at
    points where the raw argmax reports 1 Hz). The image shows the shape, the markers
    make the claim.

    Seed disagreement is drawn rather than hidden: the error bar is `peak_freq_sd`, the
    spread of the per-seed peak, and a marker is hollow where `band_agreement` falls
    below `agreement_min`, i.e. where the individual noise realisations did not agree
    with the averaged spectrum about which band dominates. Points that are not in the
    `oscillation` regime are left out entirely.
    """
    signals = signals or DEFAULT_SIGNALS
    gain = gain_spectra(spec)
    freqs = spec["freqs"]
    fmask = freqs <= fmax
    has_agreement = "band_agreement" in df.columns
    has_sd = "peak_freq_sd" in df.columns

    # One colour scale for the whole figure, so panels are comparable by eye. A
    # percentile rather than the extremes: a single hot corner would otherwise flatten
    # every other panel to the midpoint.
    imgs = {}
    for p in params:
        for sg in signals:
            sub = _tau_line(df, p, sg)
            if sub.empty or sg not in spec["signals"]:
                continue
            ids = sub["sim_id"].to_numpy()
            with np.errstate(divide="ignore", invalid="ignore"):
                img = np.log10(gain[ids, spec["signals"].index(sg), :][:, fmask])
            img = img - np.nanmedian(img, axis=1, keepdims=True)
            imgs[(p, sg)] = (sub["varied_value"].to_numpy(), img, sub)
    finite = np.concatenate([np.abs(v[1][np.isfinite(v[1])].ravel())
                             for v in imgs.values()]) if imgs else np.array([1.0])
    lim = float(np.nanpercentile(finite, 98)) if finite.size else 1.0
    lim = lim if np.isfinite(lim) and lim > 0 else 1.0

    # constrained layout: every row carries its own x label (each sweeps a different
    # parameter), and without it those labels collide with the titles of the row below
    fig, axes = plt.subplots(len(params), len(signals),
                             figsize=(3.5 * len(signals) + 1.4, 2.9 * len(params) + 1.1),
                             sharex="row", sharey=True, squeeze=False,
                             layout="constrained")
    mesh = None
    for r, p in enumerate(params):
        for c, sg in enumerate(signals):
            ax = axes[r][c]
            entry = imgs.get((p, sg))
            if entry is None:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8)
                continue
            values, img, sub = entry
            mesh = ax.pcolormesh(values, freqs[fmask], img.T, cmap="RdBu_r",
                                 vmin=-lim, vmax=lim, shading="nearest")
            for lo, hi in _BAND_SPANS.values():
                for edge in (lo, hi):
                    if edge <= fmax:
                        ax.axhline(edge, color="0.25", linewidth=0.4, alpha=0.35,
                                   zorder=2)

            osc = sub[sub["regime"] == "oscillation"]
            if not osc.empty:
                firm = (osc["band_agreement"] >= agreement_min) if has_agreement \
                    else np.ones(len(osc), dtype=bool)
                err = osc["peak_freq_sd"].to_numpy() if has_sd else None
                if err is not None:
                    ax.errorbar(osc["varied_value"], osc["peak_freq"], yerr=err,
                                fmt="none", ecolor="0.15", elinewidth=0.5, alpha=0.32,
                                capsize=1.2, zorder=3)
                ax.scatter(osc.loc[firm, "varied_value"], osc.loc[firm, "peak_freq"],
                           s=16, color="k", linewidths=0, zorder=5)
                ax.scatter(osc.loc[~firm, "varied_value"], osc.loc[~firm, "peak_freq"],
                           s=16, facecolors="none", edgecolors="k", linewidths=0.7,
                           zorder=4)
            # the error bars run past the plotted band, and sharey would otherwise let
            # one of them stretch every panel down to negative frequencies
            ax.set_ylim(float(freqs[fmask].min()), float(freqs[fmask].max()))
            # every row sweeps a different parameter, so the x label belongs on all of
            # them, not only the bottom one
            ax.set_xlabel(param_label(p))
            if r == 0:
                ax.set_title(signal_label(sg))
            if c == 0:
                ax.set_ylabel("frequency (Hz)")
    if mesh is not None:
        cbar = fig.colorbar(mesh, ax=axes.ravel().tolist(), fraction=0.02, pad=0.015)
        cbar.set_label("log10 gain, re. the median over frequency")
    fig.suptitle("Influence of the excitatory time constants on the oscillation "
                 "frequency", fontweight="bold")
    note = ("image = seed-averaged network gain spectrum, level removed per tau value; "
            "markers = flank-referenced peak +/- across-seed SD; "
            f"hollow = band_agreement < {agreement_min:g}")
    fig.supxlabel(note, fontsize=7.5, color="0.35")
    return save_figure(fig, figure_dir, f"{name}_tauRidgePeak{suffix}")


def plot_tau_amplitude(df, params, figure_dir, name, signals=None, value="amplitude",
                       reference=None, suffix=""):
    """Oscillation amplitude against each excitatory time constant, one panel per ROI.

    Amplitude is the settled signal's standard deviation, averaged over noise
    realisations; the ribbon is its across-seed SD (`<value>_sd`, written by
    run_parameter_sweep when --seeds > 1). Unlike the peak frequency this measure does
    not depend on picking a peak out of a noisy spectrum, so it is usually the cleaner
    read on what a time constant does - but a rise here is a rise in *fluctuation*,
    which is only an oscillation amplitude where the regime is oscillatory, hence the
    open markers on non-oscillating points.

    `reference` (a {param: value} dict) marks where the model itself sits on each curve.

    The y axis is logarithmic. A state transition - e3b_tau crossing ~9 ms flips A3b from
    a low-rate beta regime to a near-saturated alpha one - throws a single point an order
    of magnitude above the rest, and on a linear axis that one point flattens every other
    curve to the baseline. Such a point is also bistable across seeds (its SD exceeds its
    mean), so the ribbon is clipped to stay positive rather than being dropped.
    """
    signals = signals or DEFAULT_SIGNALS
    sd_col = f"{value}_sd"
    has_sd = sd_col in df.columns

    fig, axes = plt.subplots(1, len(signals),
                             figsize=(3.6 * len(signals) + 0.8, 3.2),
                             sharey=True, squeeze=False)
    axes = axes[0]
    for i, (ax, sg) in enumerate(zip(axes, signals)):
        for p in params:
            sub = _tau_line(df, p, sg)
            if sub.empty:
                continue
            col = TAU_COLORS.get(p, "0.4")
            x = sub["varied_value"].to_numpy()
            y = sub[value].to_numpy()
            if has_sd:
                sd = np.nan_to_num(sub[sd_col].to_numpy())
                floor = np.nanmin(y[y > 0]) * 0.5 if np.any(y > 0) else 1e-12
                ax.fill_between(x, np.maximum(y - sd, floor), y + sd, color=col,
                                alpha=0.18, linewidth=0, zorder=1)
            ax.plot(x, y, color=col, linewidth=1.4, zorder=3,
                    label=param_label(p) if i == 0 else None)
            osc = sub["regime"] == "oscillation"
            ax.scatter(x[osc.to_numpy()], y[osc.to_numpy()], s=10, color=col,
                       linewidths=0, zorder=4)
            ax.scatter(x[~osc.to_numpy()], y[~osc.to_numpy()], s=10, facecolors="none",
                       edgecolors=col, linewidths=0.6, zorder=4)
            if reference and p in reference:
                ax.axvline(float(reference[p]), color=col, linewidth=0.7, alpha=0.45,
                           linestyle=(0, (3, 2)), zorder=0)
        ax.set_yscale("log")
        ax.set_xlabel("synaptic time constant (ms)")
        ax.set_title(signal_label(sg))
        if i == 0:
            ax.set_ylabel("oscillation amplitude")
    handles, labs = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(handles, labs, frameon=False, fontsize=8, loc="best")
    fig.suptitle("Influence of the excitatory time constants on the oscillation "
                 "amplitude", fontweight="bold")
    note = ("log y; ribbon = across-seed SD, clipped positive; open markers = not "
            "classified as oscillating; dashed line = the reference value of that "
            "time constant")
    fig.text(0.5, -0.05, note, ha="center", fontsize=7.5, color="0.35")
    sns.despine(fig=fig, trim=True)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    return save_figure(fig, figure_dir, f"{name}_tauAmplitude{suffix}")
