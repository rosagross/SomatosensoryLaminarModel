"""
File: plotting_functions_step004.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Description:
    Plotting and IO for the step004 routine analysis of the model's parameter space.

    Consumes the tidy tables written by Simulations/run_parameter_sweep.py - one row per
    (parameter set x signal) in sweep_features.csv, the spectra themselves in
    sweep_spectra.hdf5 - and turns a 2-D grid sweep into figures. It knows nothing about
    the model or about which pairs were swept; step004_run_model_analysis.py decides that.

    Four kinds of figure:
      plot_roi_map            1 x 3 heatmaps, one per ROI dipole
      plot_population_grid_map  layer x cell-type panel grid of heatmaps, one area
      plot_grid_slices        the line companion: a few representative traces
      plot_gain_spectra       the network-gain spectra themselves at a few grid cells
"""

import json
import os

import h5py
import matplotlib
import numpy as np
import pandas as pd

if not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")     # these sweeps are usually plotted on a headless node

import matplotlib.pyplot as plt                            # noqa: E402
import seaborn as sns                                      # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm, TwoSlopeNorm  # noqa: E402
from matplotlib.ticker import ScalarFormatter                    # noqa: E402


# ── labels ─────────────────────────────────────────────────────────────────────
# The model calls the area A1; the figures call it S1. Both names are in the codebase
# (SomatoModel.DIPOLE_SOURCE_ALIASES reconciles them), so the mapping is spelled out here
# rather than left to whoever reads the axis.
DEFAULT_SIGNALS = ["roi_A3b", "roi_A1", "roi_S2"]
SIGNAL_LABELS = {"roi_A3b": "A3b (dipole)",
                 "roi_A1": "S1 (dipole)",
                 "roi_S2": "S2 (dipole)"}

PARAM_LABELS = {
    "coupling_strength": r"global coupling $g$",
    "strength_I": r"E/I balance $s_I$",
    "Ib_strength": "background input strength",
    "Ib_ratio_E": "background ratio E",
    "Ib_ratio_PV": "background ratio PV",
    "Ib_ratio_SST": "background ratio SST",
    "Ib_ratio_VIP": "background ratio VIP",
    "Ib_noise_std": "background noise SD",
    "Ib_noise_tau": "background noise " + r"$\tau$ (s)",
    "g_intercortical": "inter-cortical coupling",
    "g_thal": "thalamic gain " + r"$g_{thal}$",
    "sI_thal": "reticular ratio " + r"$s_{I,thal}$",
    "g_thalPOm": "POm coupling",
    "thal_EtoE": "VPM " + r"$\rightarrow$" + " VPM",
    "thal_EtoI": "VPM " + r"$\rightarrow$" + " reticular",
    "thal_ItoE": "reticular " + r"$\rightarrow$" + " VPM",
    "thal_ItoI": "reticular " + r"$\rightarrow$" + " reticular",
    "thal_POmtoI": "POm " + r"$\rightarrow$" + " reticular",
    "e3b_tau": r"$\tau$ A3b (ms)",
    "e1_tau": r"$\tau$ S1 (ms)",
    "e2_tau": r"$\tau$ S2 (ms)",
    "delay_factor": "long-range delay (s)",
    "delay_factor_short": "A3b" + r"$\leftrightarrow$" + "S1 delay (s)",
    "thal_delay_factor": "thalamic delay (s)",
    "p_2PVE": r"L4 PV$\leftarrow$E prob. (%)",
    "p_4PVE": r"L6 PV$\leftarrow$E prob. (%)",
    "Im_strength": "modulatory input strength",
}

# Layer x cell-type layout per area. Built from the canonical population order in
# Simulations/parameters.POPULATION_LABELS; VIP exists only in L2/3 and A3b, so the
# deeper layers leave that column empty.
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
    # the thalamus has no laminar structure; one row, the third column holding POm
    "Thal": ([["ThalE", "ThalI", "ThalPOm", None]], ["Thalamus"]),
    # Excitatory populations only, laid out layer x area instead of layer x cell type.
    # This is the layout the disinhibition analysis wants: "the E potential in all layers
    # and areas" is one figure, not three. ThalE is in it as a negative control - the
    # thalamus has no VIP and receives no modulatory input, so a delta there is
    # cortico-thalamic feedback and nothing else.
    "E": ([["E3b", "E1", "E1S2", "ThalE"],
           [None, "E2", "E2S2", None],
           [None, "E3", "E3S2", None],
           [None, "E4", "E4S2", None]],
          ["Layer 2/3", "Layer 4", "Layer 5", "Layer 6"]),
}
AREA_COL_LABELS = {"Thal": ["VPM (E)", "reticular (I)", "POm", ""],
                   "E": ["A3b", "S1", "S2", "VPM (E)"]}
# What the suptitle calls the layout; the key itself is not always a readable name.
AREA_TITLES = {"E": "excitatory populations"}
CELLTYPE_LABELS = ["Excitatory", "PV", "SST", "VIP"]
# The areas analyses 2 and 3 loop over. Deliberately not every key of AREA_GRIDS: "E"
# is a cross-area layout for one analysis, not a fifth area to plot everything for.
AREAS = ("A3b", "S1", "S2", "Thal")

# Okabe-Ito, colourblind-safe. Order matters: it is the order of the colour bar.
REGIME_ORDER = ["diverged", "pinned", "fixed_point", "damped", "noise_driven",
                "oscillation"]
REGIME_COLORS = {"diverged": "#000000", "pinned": "#D55E00",
                 "fixed_point": "#999999", "damped": "#F0E442",
                 "noise_driven": "#56B4E9", "oscillation": "#009E73"}
BAND_ORDER = ["none", "theta", "alpha", "beta", "gamma"]
BAND_COLORS = {"none": "#DDDDDD", "theta": "#332288", "alpha": "#EE6677",
               "beta": "#228833", "gamma": "#CCBB44"}

# Maximally distinct qualitative colours for the line companions - the traces have to be
# told apart, which a gradient palette does not do.
LINE_COLORS = ["#4477AA", "#EE6677", "#228833", "#AA3377", "#CCBB44",
               "#66CCEE", "#BBBBBB"]


def param_label(name):
    return PARAM_LABELS.get(name, name.replace("_", " "))


def signal_label(name):
    return SIGNAL_LABELS.get(name, name)


# ── style / IO ─────────────────────────────────────────────────────────────────
def set_style():
    """The repo figure style, without the tkinter screen-width query.

    plotting_style.figure_style() calls tk.Tk().winfo_screenwidth() to pick a dpi, which
    raises on a headless node. Everything else about the style is reproduced here.
    """
    sns.set_theme(style="ticks", context="paper", font="Arial",
                  rc={"font.size": 11, "figure.titlesize": 11, "axes.titlesize": 11,
                      "axes.labelsize": 11, "axes.linewidth": 0.5, "lines.linewidth": 1,
                      "lines.markersize": 3, "xtick.labelsize": 8, "ytick.labelsize": 8,
                      "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                      "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                      "legend.fontsize": 8, "legend.title_fontsize": 9,
                      "legend.frameon": False})
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42


def save_figure(fig, figure_dir, name):
    """Write one figure as both .pdf and .png (dpi 300); returns the pdf path."""
    os.makedirs(figure_dir, exist_ok=True)
    path = os.path.join(figure_dir, f"{name}.pdf")
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(os.path.join(figure_dir, f"{name}.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def load_sweep(sweep_dir):
    """(features DataFrame, sweep config dict) of one sweep directory."""
    df = pd.read_csv(os.path.join(sweep_dir, "sweep_features.csv"))
    with open(os.path.join(sweep_dir, "sweep_config.json")) as f:
        cfg = json.load(f)
    return df, cfg


def load_spectra(sweep_dir):
    """The stored spectra of one sweep, or None if the file is missing."""
    path = os.path.join(sweep_dir, "sweep_spectra.hdf5")
    if not os.path.exists(path):
        return None
    with h5py.File(path, "r") as f:
        return {"freqs": f["freqs"][:], "psd": f["psd"][:], "psd_null": f["psd_null"][:],
                "signals": [s.decode() for s in f["signals"][:]],
                "param_names": [s.decode() for s in f["param_names"][:]],
                "theta": f["theta"][:]}


def gain_spectra(spec):
    """The network gain: the run's spectrum divided by its loop-free control.

    Never read a peak off the raw spectrum. A network with its loops cut already scores
    an 'alpha peak' of 0.22 and a 'gamma peak' of 0.30 from the Ornstein-Uhlenbeck corner
    at 1/(2*pi*Ib_noise_tau) ~ 9.9 Hz plus the synaptic kernel roll-off - both inside the
    bands being measured. Dividing by the control removes exactly that.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        gain = spec["psd"] / spec["psd_null"]
    gain[~np.isfinite(gain)] = np.nan
    return gain


# ── grid helpers ───────────────────────────────────────────────────────────────
def _grid_axes(df, x, y):
    """The sorted unique values of both sweep axes."""
    return np.sort(df[x].unique()), np.sort(df[y].unique())


def _pivot(df, x, y, value, xs=None, ys=None):
    """One signal's values as a (len(ys), len(xs)) matrix, NaN where absent.

    The axes are forced rather than taken from the data so a masked or failed row does
    not silently shrink the grid and shift every cell.
    """
    xs = _grid_axes(df, x, y)[0] if xs is None else xs
    ys = _grid_axes(df, x, y)[1] if ys is None else ys
    wide = df.pivot_table(index=y, columns=x, values=value, aggfunc="first",
                          dropna=False)
    return wide.reindex(index=ys, columns=xs).to_numpy(dtype=float)


def _pivot_categorical(df, x, y, value, order, xs, ys):
    """The same, for a string column: returns indices into `order` (NaN if unknown)."""
    codes = df[[x, y, value]].copy()
    codes[value] = codes[value].map({v: i for i, v in enumerate(order)})
    wide = codes.pivot_table(index=y, columns=x, values=value, aggfunc="first",
                             dropna=False)
    return wide.reindex(index=ys, columns=xs).to_numpy(dtype=float)


def _tick_labels(values, max_ticks=6):
    """Indices and labels for at most `max_ticks` ticks along a swept axis."""
    n = len(values)
    step = max(1, int(np.ceil(n / max_ticks)))
    idx = list(range(0, n, step))
    if idx[-1] != n - 1:
        idx.append(n - 1)
    return idx, [f"{values[i]:.3g}" for i in idx]


def _format_axes(ax, xs, ys, x, y, show_x, show_y):
    xi, xl = _tick_labels(xs)
    yi, yl = _tick_labels(ys)
    ax.set_xticks(np.array(xi) + 0.5)
    ax.set_yticks(np.array(yi) + 0.5)
    ax.set_xticklabels(xl if show_x else [], rotation=45, ha="right")
    ax.set_yticklabels(yl if show_y else [])
    ax.set_xlabel(param_label(x) if show_x else "")
    ax.set_ylabel(param_label(y) if show_y else "")
    ax.set_xlim(0, len(xs))
    ax.set_ylim(0, len(ys))


def _mask(df, value, mask_non_oscillating):
    """The rows to plot: everything, or only the cells classified as oscillating.

    peak_freq and the band prominences are only meaningful where there is a peak; the
    rate columns are meaningful everywhere.
    """
    if not mask_non_oscillating:
        return df
    return df.assign(**{value: df[value].where(df["regime"] == "oscillation")})


# Cells whose value is a property of the sigmoid rather than of the parameters: a pinned
# population sits at m_max by definition (31.4 Hz for E, 166.8 Hz for PV), and a diverged
# one is a blow-up. They are still drawn - they clip to the top of the scale, and the
# colour bar is marked as extending - but they are kept out of the *limits*, because one
# pinned corner is 50-100x the rest of the map and flattens every other cell to one
# colour.
RAIL_REGIMES = ("pinned", "diverged")


def _limit_matrix(sub, x, y, value, xs, ys, mask_non_oscillating):
    """The same matrix as _pivot, with the rail cells blanked out."""
    keep = sub[~sub["regime"].isin(RAIL_REGIMES)] if "regime" in sub.columns else sub
    if keep.empty:
        keep = sub
    return _pivot(_mask(keep, value, mask_non_oscillating), x, y, value, xs, ys)


def _robust_extent(values, log_scale, lo_pct=2.0, hi_pct=98.0):
    """(vmin, vmax) over a list of matrices, ignoring the extreme tails."""
    finite = [v[np.isfinite(v)].ravel() for v in values if v is not None]
    finite = np.concatenate(finite) if finite else np.array([])
    if log_scale:
        finite = finite[finite > 0]
    if not finite.size:
        return None, None
    lo, hi = np.nanpercentile(finite, [lo_pct, hi_pct])
    if not hi > lo:                       # a flat or near-flat panel
        lo, hi = float(finite.min()), float(finite.max())
    return float(lo), float(hi)


def _draw_continuous(ax, matrix, cmap, vmin, vmax, log_scale, diverging):
    finite = matrix[np.isfinite(matrix)]
    if not finite.size:
        ax.set_facecolor("#F5F5F5")
        return None
    if log_scale:
        pos = finite[finite > 0]
        if not pos.size:
            ax.set_facecolor("#F5F5F5")
            return None
        lo = vmin if vmin is not None else float(pos.min())
        hi = vmax if vmax is not None else float(pos.max())
        if not hi > lo:
            hi = lo * 10
        norm = LogNorm(vmin=lo, vmax=hi)
        return ax.pcolormesh(np.ma.masked_invalid(matrix), cmap=cmap, norm=norm)
    lo = float(finite.min()) if vmin is None else vmin
    hi = float(finite.max()) if vmax is None else vmax
    if diverging:
        # a signed measure - a negative prominence is a dip where a peak was looked for -
        # so the colour scale has to be centred on zero, not on the data
        span = max(abs(lo), abs(hi)) or 1.0
        norm = TwoSlopeNorm(vmin=-span, vcenter=0.0, vmax=span)
        return ax.pcolormesh(np.ma.masked_invalid(matrix), cmap=cmap, norm=norm)
    if not hi > lo:
        hi = lo + 1e-12
    return ax.pcolormesh(np.ma.masked_invalid(matrix), cmap=cmap, vmin=lo, vmax=hi)


def _cbar_extend(meshes):
    """Whether the drawn data runs past the colour limits, so the bar can say so.

    The limits deliberately exclude the pinned/diverged cells (see RAIL_REGIMES); those
    cells are still drawn and clip to the end of the scale, and an arrow on the colour bar
    is what keeps that honest rather than silently reading as "equal to the maximum".
    """
    below = above = False
    for mesh in meshes:
        arr = mesh.get_array()
        if arr is None:
            continue
        data = np.asarray(arr.compressed() if np.ma.isMaskedArray(arr) else arr,
                          dtype=float)
        data = data[np.isfinite(data)]
        if not data.size:
            continue
        if mesh.norm.vmin is not None and (data < mesh.norm.vmin).any():
            below = True
        if mesh.norm.vmax is not None and (data > mesh.norm.vmax).any():
            above = True
    return {(False, False): "neither", (True, False): "min",
            (False, True): "max", (True, True): "both"}[(below, above)]


def _log_colorbar(fig, mesh, ax, label, extend="neither"):
    """A log colour bar with explicit ticks at both ends and the geometric middle.

    Several panels span less than a decade, and LogNorm's default locator then puts no
    tick inside the range at all - the bar comes out completely unlabelled.
    """
    cb = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.03, extend=extend)
    lo, hi = mesh.norm.vmin, mesh.norm.vmax
    ticks = [lo, float(np.sqrt(lo * hi)), hi]
    # LogNorm's default minor locator labels its own decades on top of these, which on a
    # narrow bar prints "6 x 10^-1" straight through "0.67". Turn the minor layer off.
    cb.ax.minorticks_off()
    cb.set_ticks(ticks)
    cb.set_ticklabels([f"{t:.2g}" for t in ticks])
    cb.ax.tick_params(labelsize=7)
    if label:
        cb.set_label(label, fontsize=8)
    return cb


def _categorical_legend(fig, order, colors, title):
    handles = [plt.Line2D([], [], marker="s", linestyle="none", markersize=7,
                          markerfacecolor=colors[v], markeredgecolor="none",
                          label=v.replace("_", " "))
               for v in order]
    fig.legend(handles=handles, title=title, loc="center left",
               bbox_to_anchor=(1.0, 0.5))


def _despine(ax):
    """Drop the top/right spines and trim the rest to the ticks that are actually in view.

    seaborn's own trim=True is not usable here, for two reasons:

    - It raises on a panel with no data to put ticks on - a masked measure (peak_freq
      where nothing oscillates) leaves an axis with no finite points and therefore no
      ticks, which trim=True indexes into blindly.
    - It clips each spine to the outermost *major* tick without checking that any major
      tick is in view. On a log axis spanning less than a decade there is none: the
      sig_alpha_gain slices sit at 1.30e3-1.67e3, whose major ticks are 1e3 and 1e4, so
      seaborn picks firsttick=1e4 and lasttick=1e3 and draws the left spine from 1e3 to
      1e4 - far above the axes box. bbox_inches="tight" then grows the saved figure to
      contain that spine, which is how a 2.8 in figure came out 16 in tall.

    Trimming to the in-view ticks does what trim was for and cannot leave the axes.
    """
    sns.despine(ax=ax)
    for side, axis, lim in (("bottom", ax.xaxis, ax.get_xlim()),
                            ("left", ax.yaxis, ax.get_ylim())):
        lo, hi = min(lim), max(lim)
        ticks = [t for t in axis.get_majorticklocs() if lo <= t <= hi]
        # fewer than two ticks in view leaves nothing to trim to; the view limits are
        # then the honest bounds
        ax.spines[side].set_bounds(*((ticks[0], ticks[-1]) if len(ticks) > 1 else (lo, hi)))


# ── the ROI dipole maps ────────────────────────────────────────────────────────
def plot_roi_map(df, x, y, value, figure_dir, name, signals=None, categorical=None,
                 cmap="magma", vmin=None, vmax=None, label=None, log_scale=False,
                 diverging=False, mask_non_oscillating=True):
    """A 1 x 3 row of heatmaps over the (x, y) grid, one per ROI dipole.

    `categorical` is 'regime' or 'band' (the column is then read as a string and drawn
    with the discrete colour map); otherwise `value` is a continuous column.
    """
    signals = signals or DEFAULT_SIGNALS
    xs, ys = _grid_axes(df, x, y)
    present = [s for s in signals if (df["signal"] == s).any()]
    if not present:
        return None

    fig, axes = plt.subplots(1, len(present), figsize=(3.0 * len(present) + 1.2, 3.0),
                             squeeze=False)
    axes = axes[0]
    # startswith, not ==, so a derived column (regime_lo / regime_hi in the disinhibition
    # analysis) still gets the regime palette rather than silently being drawn with the
    # band one
    regime_like = bool(categorical) and categorical.startswith("regime")
    order = REGIME_ORDER if regime_like else BAND_ORDER
    colors = REGIME_COLORS if regime_like else BAND_COLORS

    # One scale across the three panels, so the single colour bar below describes all of
    # them - and taken over the rail-free cells, so a pinned corner does not flatten the
    # rest of the map.
    if not categorical:
        subs = [df[df["signal"] == s] for s in present]
        lo, hi = _robust_extent([_limit_matrix(sb, x, y, value, xs, ys,
                                               mask_non_oscillating) for sb in subs],
                                log_scale)
        vmin = lo if vmin is None else vmin
        vmax = hi if vmax is None else vmax

    mesh = None
    for ax, sig in zip(axes, present):
        sub = df[df["signal"] == sig]
        if categorical:
            matrix = _pivot_categorical(sub, x, y, categorical, order, xs, ys)
            cmap_c = ListedColormap([colors[v] for v in order])
            norm = BoundaryNorm(np.arange(-0.5, len(order)), cmap_c.N)
            ax.pcolormesh(np.ma.masked_invalid(matrix), cmap=cmap_c, norm=norm)
        else:
            matrix = _pivot(_mask(sub, value, mask_non_oscillating), x, y, value, xs, ys)
            mesh = _draw_continuous(ax, matrix, cmap, vmin, vmax, log_scale, diverging)
        ax.set_title(signal_label(sig))
        _format_axes(ax, xs, ys, x, y, True, ax is axes[0])

    if categorical:
        _categorical_legend(fig, order, colors, categorical)
    elif mesh is not None:
        extend = _cbar_extend([mesh])
        if log_scale:
            _log_colorbar(fig, mesh, list(axes), label or value, extend)
        else:
            cb = fig.colorbar(mesh, ax=list(axes), fraction=0.046, pad=0.03,
                              extend=extend)
            cb.set_label(label or value, fontsize=8)
            cb.ax.tick_params(labelsize=7)
    return save_figure(fig, figure_dir, name)


# ── the layer x cell-type panel grids ──────────────────────────────────────────
def _panel_limits(matrices, rows, share_scale, log_scale):
    """(vmin, vmax) per panel, per column or per figure.

    A firing-rate panel grid needs its colour scale per panel: at one operating point E1
    sits near 27 Hz and E2 near 2 Hz, and the sigmoid ceiling alone is 31.4 Hz for E
    against 166.8 Hz for PV, so any shared scale flattens most panels to one colour.
    'figure' is right for a normalised quantity like rate_level, which is 0-1 everywhere
    by construction and where the comparison between panels is the whole point.
    """
    def extent(values):
        return _robust_extent(values, log_scale)

    if share_scale == "figure":
        lo, hi = extent([m for r in matrices for m in r if m is not None])
        return {(i, j): (lo, hi) for i, r in enumerate(matrices)
                for j, _ in enumerate(r)}
    if share_scale == "column":
        limits = {}
        n_col = max(len(r) for r in matrices)
        for j in range(n_col):
            lo, hi = extent([r[j] for r in matrices
                             if j < len(r) and r[j] is not None])
            for i in range(len(matrices)):
                limits[(i, j)] = (lo, hi)
        return limits
    # 'panel': each panel gets its own limits, still taken over its rail-free cells -
    # returning (None, None) here would let _draw_continuous auto-scale from the full
    # matrix, and a single pinned cell at m_max then flattens the whole panel to black.
    return {(i, j): extent([m]) for i, r in enumerate(matrices)
            for j, m in enumerate(r)}


def plot_population_grid_map(df, x, y, value, area, figure_dir, name, categorical=None,
                             cmap="magma", label=None, share_scale="panel",
                             log_scale=False, vmin=None, vmax=None,
                             mask_non_oscillating=True, diverging=False):
    """Layer x cell-type panel grid of 2-D maps for one area."""
    if area not in AREA_GRIDS:
        raise ValueError(f"unknown area {area!r}; choose from {sorted(AREA_GRIDS)}")
    pop_grid, row_labels = AREA_GRIDS[area]
    col_labels = AREA_COL_LABELS.get(area, CELLTYPE_LABELS)
    xs, ys = _grid_axes(df, x, y)
    # startswith, not ==, so a derived column (regime_lo / regime_hi in the disinhibition
    # analysis) still gets the regime palette rather than silently being drawn with the
    # band one
    regime_like = bool(categorical) and categorical.startswith("regime")
    order = REGIME_ORDER if regime_like else BAND_ORDER
    colors = REGIME_COLORS if regime_like else BAND_COLORS

    matrices, for_limits = [], []
    for row in pop_grid:
        line, line_lim = [], []
        for pop in row:
            sub = df[df["signal"] == pop] if pop else df.iloc[0:0]
            if pop is None or sub.empty:
                line.append(None)
                line_lim.append(None)
            elif categorical:
                line.append(_pivot_categorical(sub, x, y, categorical, order, xs, ys))
                line_lim.append(None)
            else:
                line.append(_pivot(_mask(sub, value, mask_non_oscillating), x, y, value,
                                   xs, ys))
                line_lim.append(_limit_matrix(sub, x, y, value, xs, ys,
                                              mask_non_oscillating))
        matrices.append(line)
        for_limits.append(line_lim)
    if all(m is None for row in matrices for m in row):
        return None

    limits = ({} if categorical
              else _panel_limits(for_limits, row_labels, share_scale, log_scale))
    anchors = _panel_anchors([[m is not None for m in row] for row in matrices])
    n_row, n_col = len(pop_grid), len(pop_grid[0])
    fig, axes = plt.subplots(n_row, n_col, figsize=(2.5 * n_col + 1.5, 2.4 * n_row + 0.6),
                             squeeze=False)
    meshes = []
    for i in range(n_row):
        for j in range(n_col):
            ax = axes[i][j]
            matrix = matrices[i][j]
            if matrix is None:
                ax.axis("off")
                continue
            if categorical:
                cmap_c = ListedColormap([colors[v] for v in order])
                norm = BoundaryNorm(np.arange(-0.5, len(order)), cmap_c.N)
                ax.pcolormesh(np.ma.masked_invalid(matrix), cmap=cmap_c, norm=norm)
            else:
                lo, hi = limits[(i, j)]
                lo = vmin if vmin is not None else lo
                hi = vmax if vmax is not None else hi
                mesh = _draw_continuous(ax, matrix, cmap, lo, hi, log_scale, diverging)
                if mesh is not None:
                    meshes.append((mesh, ax))
            row_anchor, col_anchor = anchors
            if col_anchor.get(j) == i:
                ax.set_title(col_labels[j] if j < len(col_labels) else "")
            if row_anchor.get(i) == j:
                ax.text(-0.35, 0.5, row_labels[i], transform=ax.transAxes,
                        rotation=90, va="center", ha="center", fontsize=10)
            # the x labels go on the lowest drawn panel of the column, not on row n-1,
            # which in a ragged layout is often turned off
            last_row = max((r for r in range(n_row) if matrices[r][j] is not None),
                           default=n_row - 1)
            _format_axes(ax, xs, ys, x, y, i == last_row, row_anchor.get(i) == j)

    if categorical:
        _categorical_legend(fig, order, colors, categorical)
    elif share_scale == "figure" and meshes:
        cb = fig.colorbar(meshes[0][0], ax=axes.ravel().tolist(), fraction=0.03,
                          pad=0.02, extend=_cbar_extend([m for m, _ in meshes]))
        cb.set_label(label or value, fontsize=8)
        cb.ax.tick_params(labelsize=7)
    else:
        for mesh, ax in meshes:
            extend = _cbar_extend([mesh])
            if log_scale:
                _log_colorbar(fig, mesh, ax, None, extend)
            else:
                cb = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.03, extend=extend)
                cb.ax.tick_params(labelsize=6)
        if label:
            fig.suptitle(f"{AREA_TITLES.get(area, area)} - {label}", y=1.01)
    if not categorical and share_scale == "figure" and label:
        fig.suptitle(f"{AREA_TITLES.get(area, area)} - {label}", y=1.01)
    return save_figure(fig, figure_dir, name)


def _plain_log_ticks(ax):
    """Plain numbers on a log y axis whose view spans less than a decade.

    Under a decade there is no major tick in view, so matplotlib labels the minor ones
    with LogFormatterSciNotation: an alpha gain of 1.3e3-1.7e3 reads "1.3 x 10^3,
    1.4 x 10^3, ..." - five superscripts carrying nothing the mantissa does not already
    say. Only the sub-decade case is touched; over a wider range the exponents are the
    point of the log scale.
    """
    lo, hi = ax.get_ylim()
    if not (lo > 0 and hi / lo < 10):
        return
    for setter in (ax.yaxis.set_major_formatter, ax.yaxis.set_minor_formatter):
        fmt = ScalarFormatter()
        fmt.set_scientific(False)
        setter(fmt)


# ── the line companions ────────────────────────────────────────────────────────
def _representative(values, n_lines):
    """A few values of a swept axis, spread evenly across it including both ends.

    A small hand-picked subset with a legend reads better than the whole dense sweep
    under a colour bar - the point is to see where the behaviour diverges.
    """
    if len(values) <= n_lines:
        return list(values)
    idx = np.unique(np.linspace(0, len(values) - 1, n_lines).round().astype(int))
    return [values[i] for i in idx]


def plot_grid_slices(df, x, y, value, figure_dir, name, signals=None, n_lines=5,
                     label=None, log_scale=False, mask_non_oscillating=False):
    """`value` against x, one line per representative value of y, one panel per signal."""
    signals = signals or DEFAULT_SIGNALS
    present = [s for s in signals if (df["signal"] == s).any()]
    if not present:
        return None
    xs, ys = _grid_axes(df, x, y)
    picks = _representative(list(ys), n_lines)

    fig, axes = plt.subplots(1, len(present), figsize=(3.2 * len(present) + 1.0, 2.8),
                             squeeze=False, sharey=True)
    axes = axes[0]
    for ax, sig in zip(axes, present):
        sub = _mask(df[df["signal"] == sig], value, mask_non_oscillating)
        for k, yv in enumerate(picks):
            line = sub[np.isclose(sub[y], yv)].sort_values(x)
            ax.plot(line[x], line[value], marker="o", markersize=3,
                    color=LINE_COLORS[k % len(LINE_COLORS)],
                    label=f"{param_label(y)} = {yv:.3g}")
        ax.set_title(signal_label(sig))
        ax.set_xlabel(param_label(x))
        if log_scale:
            ax.set_yscale("log")
    for ax in axes:
        # after every line is drawn, so the shared view limits are final
        if log_scale:
            _plain_log_ticks(ax)
        _despine(ax)
    axes[0].set_ylabel(label or value)
    axes[-1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    return save_figure(fig, figure_dir, name)


# ── the population panel grid, as lines ────────────────────────────────────────
def _panel_grid(area):
    """(pop_grid, row_labels, col_labels) of one AREA_GRIDS entry."""
    if area not in AREA_GRIDS:
        raise ValueError(f"unknown area {area!r}; choose from {sorted(AREA_GRIDS)}")
    pop_grid, row_labels = AREA_GRIDS[area]
    return pop_grid, row_labels, AREA_COL_LABELS.get(area, CELLTYPE_LABELS)


def _bottom_panels(present):
    """The lowest drawn panel of each column: {j: i}.

    sharex=True hides the x tick labels of every row but the last, which in a ragged
    layout is the wrong row - the "E" grid has A3b in row 0 only, so its single panel
    would come out with no tick labels and no axis label at all.
    """
    return {j: max(i for i in range(len(present)) if present[i][j])
            for j in range(len(present[0]))
            if any(present[i][j] for i in range(len(present)))}


def _panel_anchors(present):
    """Which panel of each row carries the row label, and of each column the title.

    Not simply (i == 0) and (j == 0). The "E" layout has A3b in column 0 for the L2/3 row
    only, so anchoring the row label to column 0 silently drops the labels of every other
    row - the figure then has four unlabelled rows of layers. Anchor to the first panel
    that is actually drawn instead.
    """
    row_anchor, col_anchor = {}, {}
    for i, row in enumerate(present):
        for j, ok in enumerate(row):
            if not ok:
                continue
            row_anchor.setdefault(i, j)
            col_anchor.setdefault(j, i)
    return row_anchor, col_anchor


def _panel_frame(ax, i, j, row_labels, col_labels, anchors):
    """The title / row label / spine treatment shared by the panel-grid plotters."""
    row_anchor, col_anchor = anchors
    if col_anchor.get(j) == i:
        ax.set_title(col_labels[j] if j < len(col_labels) else "")
    if row_anchor.get(i) == j:
        ax.text(-0.35, 0.5, row_labels[i], transform=ax.transAxes,
                rotation=90, va="center", ha="center", fontsize=10)
    _despine(ax)


def plot_population_grid_lines(df, x, value, area, figure_dir, name, hue=None,
                               n_lines=5, label=None, log_scale=False, share_y=False,
                               mark=None):
    """`value` against `x` for every population of one area, in the same panel grid.

    The line companion to plot_population_grid_map, for a 1-D sweep: one panel per
    population, one line (or one per representative level of `hue`) in each. `mark` draws
    a vertical reference line, e.g. the two levels a delta was taken between.
    """
    pop_grid, row_labels, col_labels = _panel_grid(area)
    n_row, n_col = len(pop_grid), len(pop_grid[0])
    if value not in df.columns:
        return None
    levels = _representative(list(np.sort(df[hue].unique())), n_lines) if hue else [None]

    present = [[bool(pop) and (df["signal"] == pop).any() for pop in row]
               for row in pop_grid]
    anchors, bottom = _panel_anchors(present), _bottom_panels(present)
    fig, axes = plt.subplots(n_row, n_col, figsize=(2.6 * n_col + 1.0, 2.1 * n_row + 0.6),
                             squeeze=False, sharex=True, sharey=share_y)
    drawn = False
    for i in range(n_row):
        for j in range(n_col):
            ax = axes[i][j]
            pop = pop_grid[i][j]
            sub = df[df["signal"] == pop] if pop else df.iloc[0:0]
            if pop is None or sub.empty:
                ax.axis("off")
                continue
            for k, lv in enumerate(levels):
                line = (sub if lv is None else sub[np.isclose(sub[hue], lv)]).sort_values(x)
                if line.empty:
                    continue
                ax.plot(line[x], line[value], marker="o", markersize=2.5,
                        color=LINE_COLORS[k % len(LINE_COLORS)],
                        label=None if lv is None else f"{param_label(hue)} = {lv:.3g}")
                drawn = True
            for m in (mark or []):
                ax.axvline(float(m), color="0.6", linewidth=0.6, linestyle="--",
                           zorder=0)
            if log_scale:
                ax.set_yscale("log")
                _plain_log_ticks(ax)
            if bottom.get(j) == i:
                ax.set_xlabel(param_label(x))
                ax.tick_params(labelbottom=True)
            # the population's own name, inside the panel: the column header names the
            # area and the row label the layer, neither of which pins down the population
            ax.text(0.03, 0.95, pop, transform=ax.transAxes, va="top", ha="left",
                    fontsize=8, color="0.35")
            _panel_frame(ax, i, j, row_labels, col_labels, anchors)
    if not drawn:
        plt.close(fig)
        return None
    axes[0][0].set_ylabel(label or value)
    if hue:
        handles, labels_ = axes[0][0].get_legend_handles_labels()
        if handles:
            fig.legend(handles=handles, labels=labels_, loc="center left",
                       bbox_to_anchor=(1.0, 0.5))
    return save_figure(fig, figure_dir, name)


def plot_traces(traces, figure_dir, name, area="E", kind="potentials", label=None,
                window=None):
    """Time courses of every population of one area, one coloured line per input level.

    `traces` is what step004_run_model_analysis.run_disinh_traces wrote: `levels`, `t`,
    `potentials`/`rates` as (n_levels, n_pop, n_steps), `labels`, `settle_s`. `window`
    limits the plotted span in seconds after the settle point - the whole 6 s of settled
    run is unreadable at this panel size, and a rhythm is judged on a few cycles.
    """
    pop_grid, row_labels, col_labels = _panel_grid(area)
    n_row, n_col = len(pop_grid), len(pop_grid[0])
    labels = list(traces["labels"])
    data, t, levels = traces[kind], np.asarray(traces["t"]), np.asarray(traces["levels"])

    keep = t >= float(traces["settle_s"])
    if window:
        keep &= t <= float(traces["settle_s"]) + float(window)

    present = [[bool(pop) and pop in labels for pop in row] for row in pop_grid]
    anchors, bottom = _panel_anchors(present), _bottom_panels(present)
    fig, axes = plt.subplots(n_row, n_col, figsize=(2.6 * n_col + 1.0, 2.1 * n_row + 0.6),
                             squeeze=False, sharex=True)
    for i in range(n_row):
        for j in range(n_col):
            ax = axes[i][j]
            pop = pop_grid[i][j]
            if pop is None or pop not in labels:
                ax.axis("off")
                continue
            idx = labels.index(pop)
            for k, lv in enumerate(levels):
                ax.plot(t[keep], data[k, idx][keep],
                        color=LINE_COLORS[k % len(LINE_COLORS)], linewidth=0.8,
                        label=f"{param_label('Im_strength')} = {lv:.3g}")
            if bottom.get(j) == i:
                ax.set_xlabel("time (s)")
                ax.tick_params(labelbottom=True)
            ax.text(0.03, 0.95, pop, transform=ax.transAxes, va="top", ha="left",
                    fontsize=8, color="0.35")
            _panel_frame(ax, i, j, row_labels, col_labels, anchors)
    axes[0][0].set_ylabel(label or kind)
    handles, labels_ = axes[0][0].get_legend_handles_labels()
    fig.legend(handles=handles, labels=labels_, loc="center left",
               bbox_to_anchor=(1.0, 0.5))
    return save_figure(fig, figure_dir, name)


# ── the spectra themselves ─────────────────────────────────────────────────────
def plot_gain_spectra(spec, df, x, y, figure_dir, name, signals=None, n_lines=5,
                      fmax=45.0):
    """Network-gain spectra at a few representative cells along x, one panel per signal.

    Each trace has its own median over frequency subtracted. The ROI dipoles sit at a
    near-constant ~1e6 gain over the loop-free control (the control's dipole is tiny), so
    a raw log10 scale spans 5.9-7.0 and flattens the ~0.5 decades of real structure to a
    flat line. Removing the level is what makes the shape visible.
    """
    if spec is None or not len(spec["freqs"]):
        return None
    signals = signals or DEFAULT_SIGNALS
    present = [s for s in signals if s in spec["signals"]]
    if not present:
        return None

    gain = gain_spectra(spec)
    freqs = spec["freqs"]
    keep = freqs <= fmax
    # sim_id indexes the first axis of psd; the y value is held at the grid centre so the
    # traces differ in x alone
    ids = df[["sim_id", x, y]].drop_duplicates()
    ys = np.sort(ids[y].unique())
    y_hold = ys[len(ys) // 2]
    at_y = ids[np.isclose(ids[y], y_hold)].sort_values(x)
    picks = _representative(list(at_y[x].to_numpy()), n_lines)

    fig, axes = plt.subplots(1, len(present), figsize=(3.2 * len(present) + 1.0, 2.8),
                             squeeze=False)
    axes = axes[0]
    for ax, sig in zip(axes, present):
        j = spec["signals"].index(sig)
        for k, xv in enumerate(picks):
            row = at_y[np.isclose(at_y[x], xv)]
            if row.empty:
                continue
            sim_id = int(row["sim_id"].iloc[0])
            trace = gain[sim_id, j][keep]
            if not np.isfinite(trace).any():
                continue
            with np.errstate(divide="ignore", invalid="ignore"):
                log_gain = np.log10(np.where(trace > 0, trace, np.nan))
            ax.plot(freqs[keep], log_gain - np.nanmedian(log_gain),
                    color=LINE_COLORS[k % len(LINE_COLORS)],
                    label=f"{param_label(x)} = {xv:.3g}")
        ax.axvspan(8, 13, color="0.9", zorder=0)     # the alpha band, for reference
        ax.set_title(signal_label(sig))
        ax.set_xlabel("frequency (Hz)")
        _despine(ax)
    axes[0].set_ylabel("network gain (log10, level removed)")
    axes[-1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
                    title=f"{param_label(y)} = {y_hold:.3g}")
    return save_figure(fig, figure_dir, name)
