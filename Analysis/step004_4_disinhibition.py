"""
File: step004_4_disinhibition.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-09-01
Description: Script for analyzing what the modulatory (VIP) input does to the excitatory
             populations, across parameter space.

Simulations/parameters.py builds the modulatory weight column as

    Wm[classes == "VIP", 0] = 1 * mI_cellcount

so Im_strength lands on the three VIP populations (VIP3b, VIP1, VIP1S2) and on nothing
else. Raising it drives VIP, which suppresses SST, which releases E: disinhibition.

WHETHER THAT ACTUALLY HAPPENS DEPENDS ON THE BACKGROUND INPUT, and the analysis is
centred at Ib_strength = 80 for that reason (--ib-centre; the model's own baseline is 7).
The sign of the E delta flips between Ib 40 and 50:

    Ib_strength      7     20     40  |     50     60     80    100
    dE1 (Hz)     -0.24  -0.40  -1.60  |  +3.99  +4.27  +4.44  +4.16
    SST1 at Im=0  2.04   2.86   5.79  |  22.01  25.81  33.05  39.45

The reason is SST's headroom, which is what the disinhibitory route has to spend. At the
baseline background SST1 sits at 2 Hz, so releasing it entirely is worth about
2.04 * W[E1,SST1] = 2.04 * 29.01 ~= +59 on E1's potential, while VIP's own 36 Hz rise
delivers 36 * W[E1,VIP1] = 36 * 5.12 ~= -185 straight onto E. VIP inhibits E directly as
well as via SST, and at low background the direct term wins about 3:1 - the modulatory
input then SUPPRESSES the excitatory populations. Raise the background and SST1 sits at
25-40 Hz instead: the same release is now worth ten times more, and it wins.

Two questions:

1. What does the E potential do, in every layer and area, as the VIP drive increases?
   -> the dose-response figures, from a 1-D sweep of Im_strength (0-20), plus the
      potential time courses themselves at a few levels.
2. How large is the E delta - E at high Im_strength minus E at low - and how does it move
   across the same parameter planes as the other step004 analyses?
   -> the delta maps, from grids run twice, once per Im level.

Driven by step004_run_model_analysis.py, which runs those sweeps; like the other three
analyses this file only reads what the sweep wrote and plots it. It never simulates.

Why the delta is a paired difference and not a difference of noise draws:
SomatoModel.add_background_noise does `np.random.default_rng(self.Ib_noise_seed)` on every
apply_params, and run_parameter_sweep.make_base_params pins that seed for a whole sweep.
The low-Im and high-Im runs of the same grid cell therefore see the *identical*
Ornstein-Uhlenbeck realisation, and subtracting them cancels it exactly. Break that (a
per-point seed, say) and the delta measures noise.

Four ways to misread these figures, all of them worth checking before believing one:

- VIP exists only in L2/3 and in A3b. The modulatory input reaches exactly three
  populations, so a delta in L4/L5/L6 is a network effect propagated through the column,
  never a direct one. That propagation is the interesting part, but it is not evidence of
  modulatory input to those layers.
- Saturation fakes a plateau. If the high Im level pins VIP against its sigmoid, the delta
  stops growing and the map flattens - which looks exactly like a parameter region that
  does not care. vip_dose_rate_level and the *_E_regime_hi maps are how the two are told
  apart.
- At Ib = 80 the default operating point is itself close to the ceiling: the sign flip
  between Ib 40 and 50 coincides with PV1 pinning (rate_level 0.05 -> 0.997), and E1 runs
  at 0.82-0.96 of m_max. Part of why E rises there is that PV is railed and can no longer
  follow it. Disinhibition off the rail does exist - a coarse scan at Ib = 80 found it at
  coupling_strength ~ 0.5 with strength_I 0.7-1.3, and at coupling_strength 1 with
  strength_I 1.3 (dE1 = +5.5 with E1 at 0.58 of m_max and PV1 at 0.83) - and finding it is
  exactly what the delta maps plus the `saturated` column are for. Read a positive delta
  together with rate_level before calling it a mechanism.
- A negative delta is a result, not a bug. VIP contacts PV as well as SST, so net
  inhibition in some layers is a real prediction of this wiring. Every delta map is drawn
  on a diverging scale centred at zero so the sign reads straight off the figure.
- ThalE is the control. The thalamus has no VIP and no modulatory input, so its delta
  should be near zero; a large one means the interpretation, not the model, is wrong.

Parameter combinations are defined in step004_run_model_analysis.py.
"""

import os
import sys

import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import plotting_functions_step004 as pf                            # noqa: E402


# The measures differenced between the two Im levels. (column, delta name, cmap, label,
# kwargs for plot_population_grid_map).
#
# share_scale is 'panel' for the two absolute deltas on purpose: the question is how the
# delta moves *across a plane*, and L2/3 - which is where VIP actually is - runs an order
# of magnitude above L5/L6, so a figure-wide scale would leave the deeper layers blank.
# rel_delta is the normalised counterpart and shares one scale, because comparing layers
# is the whole point of normalising.
DELTA_VALUES = [
    ("delta_mean_level", "RdBu_r", "E potential, high Im - low Im",
     dict(share_scale="panel")),
    ("delta_mean_rate", "RdBu_r", "firing rate, high Im - low Im (Hz)",
     dict(share_scale="panel")),
    ("rel_delta", "PuOr_r", "potential change / |potential at low Im|",
     dict(share_scale="figure")),
]

# Below this the low-Im potential is too small for a ratio to mean anything: a population
# sitting at ~0 turns any absolute change into an enormous relative one.
REL_DELTA_FLOOR = 1e-3

# mean_rate / m_max at or above which a population counts as against its sigmoid ceiling.
# This, and not the regime label, is the saturation test for a delta. classify_regime asks
# what the *spectrum* looks like, so a cell sitting at 31.39 Hz against m_max = 31.4 is
# routinely labelled 'oscillation' or 'noise_driven' - it still fluctuates, it just does
# so on the rail. In the coupling x inhibition plane nine of 121 cells are there at low
# Im and near zero at high Im, giving a delta of the full 31 Hz height of the sigmoid;
# with those inside the colour limits every other cell in the panel is one flat colour.
CEILING_LEVEL = 0.95

# What goes in the per-pair summary table.
SUMMARY_COLUMNS = ["signal", "signal_kind", "mean_level_lo", "mean_level_hi",
                   "delta_mean_level", "rel_delta", "mean_rate_lo", "mean_rate_hi",
                   "delta_mean_rate", "delta_rate_level", "rate_level_hi",
                   "regime_lo", "regime_hi", "saturated"]

# The E populations the line companion is drawn for - one per area, all L2/3, which is
# where VIP is.
E_SLICE_POPS = ["E3b", "E1", "E1S2"]


def im_levels(df):
    """The Im_strength levels present in a paired sweep, low first."""
    return sorted(float(v) for v in df["Im_strength"].unique())


def delta_table(df, x, y):
    """One row per (x, y, signal): the high-Im measures minus the low-Im ones.

    The two halves come from the same np.linspace, so the axis values are bit-identical
    and the merge is exact - no tolerance join is needed or wanted here.
    """
    if "Im_strength" not in df.columns:
        raise ValueError("no Im_strength column - this is not a paired disinhibition "
                         "sweep")
    if x == "Im_strength" or y == "Im_strength":
        raise ValueError(f"the grid already varies Im_strength ({x} x {y}); a delta "
                         "between two Im levels is undefined there")
    levels = im_levels(df)
    if len(levels) != 2:
        raise ValueError(f"expected exactly two Im_strength levels, found {levels}")
    lo, hi = levels

    keys = [x, y, "signal", "signal_kind"]
    carry = [c for c in ("mean_level", "mean_rate", "rate_level", "regime", "amplitude")
             if c in df.columns]
    a = df[np.isclose(df["Im_strength"], lo)][keys + carry]
    b = df[np.isclose(df["Im_strength"], hi)][keys + carry]
    out = a.merge(b, on=keys, suffixes=("_lo", "_hi"), how="inner")

    out["delta_mean_level"] = out["mean_level_hi"] - out["mean_level_lo"]
    if "mean_rate_lo" in out.columns:
        out["delta_mean_rate"] = out["mean_rate_hi"] - out["mean_rate_lo"]
    if "rate_level_lo" in out.columns:
        out["delta_rate_level"] = out["rate_level_hi"] - out["rate_level_lo"]
    denom = out["mean_level_lo"].abs()
    out["rel_delta"] = np.where(denom > REL_DELTA_FLOOR,
                                out["delta_mean_level"] / denom.where(denom > 0, 1.0),
                                np.nan)

    # _mask and _limit_matrix both read a column literally called 'regime', and use it to
    # keep sigmoid-property cells out of the colour limits. For a delta a cell counts as
    # one if *either* endpoint is a rail regime or at the ceiling (CEILING_LEVEL): the
    # corner at high coupling / low inhibition sits at m_max with the modulation off and
    # is driven to ~0 by it, so its delta is the height of the sigmoid and says nothing
    # about the parameters. Those cells are still drawn - they clip, and the colour bar
    # grows an arrow - but they do not set the scale. `saturated` keeps the flag readable
    # in the summary table.
    ceiling = np.zeros(len(out), dtype=bool)
    for col in ("rate_level_lo", "rate_level_hi"):
        if col in out.columns:
            ceiling |= out[col].fillna(0.0).to_numpy() >= CEILING_LEVEL
    rail = (out["regime_lo"].isin(pf.RAIL_REGIMES)
            | out["regime_hi"].isin(pf.RAIL_REGIMES)
            | ceiling)
    # An ROI dipole has no firing rate of its own, so it inherits the flag from the cell:
    # if any population in that grid cell was on the rail, the dipole's delta is the
    # sigmoid's height too.
    on_rail = out[rail][[x, y]].drop_duplicates().assign(_rail=True)
    rail |= out.merge(on_rail, on=[x, y], how="left")["_rail"].fillna(False).to_numpy()

    out["saturated"] = rail
    out["regime"] = np.where(rail, "pinned", out["regime_hi"])
    out["im_low"], out["im_high"] = lo, hi
    out["x_param"], out["y_param"] = x, y
    return out


def run_disinhibition_analysis(sweep_dir, figure_dir, name):
    """Every delta figure of one paired (two Im level) grid sweep."""
    df, _ = pf.load_sweep(sweep_dir)
    if df.empty:
        print(f"    {name}: empty sweep - skipped")
        return []
    x, y = df["x_param"].iloc[0], df["y_param"].iloc[0]
    delta = delta_table(df, x, y)
    pops = delta[delta["signal_kind"] == "potential"]
    rois = delta[delta["signal_kind"] == "dipole"]

    os.makedirs(figure_dir, exist_ok=True)
    made = []

    # the headline figures: the E delta over the plane, all layers and areas
    for value, cmap, label, kwargs in DELTA_VALUES:
        if value not in pops.columns:
            continue
        made.append(pf.plot_population_grid_map(
            pops, x, y, value, "E", figure_dir, f"{name}_E_{value}",
            cmap=cmap, label=label, diverging=True, mask_non_oscillating=False,
            **kwargs))

    # both endpoints' regimes: a big delta between a pinned cell and a silent one is the
    # height of the sigmoid, not a mechanism, and these two maps are where that is seen.
    # Those cells are kept out of the colour limits above and clip with an arrow instead.
    for endpoint in ("regime_lo", "regime_hi"):
        made.append(pf.plot_population_grid_map(
            pops, x, y, endpoint, "E", figure_dir, f"{name}_E_{endpoint}",
            categorical=endpoint))

    # the same delta on the interneurons, which is where the chain can be checked:
    # Im up -> VIP up -> SST down -> E up
    for area in ("A3b", "S1", "S2"):
        made.append(pf.plot_population_grid_map(
            pops, x, y, "delta_mean_rate", area, figure_dir,
            f"{name}_{area}_delta_rate", cmap="RdBu_r",
            label="firing rate, high Im - low Im (Hz)", diverging=True,
            share_scale="panel", mask_non_oscillating=False))

    # and at the ROI level
    if not rois.empty:
        made.append(pf.plot_roi_map(
            rois, x, y, "delta_mean_level", figure_dir, f"{name}_roi_delta",
            cmap="RdBu_r", label="dipole, high Im - low Im", diverging=True,
            mask_non_oscillating=False))

    made.append(pf.plot_grid_slices(
        pops, x, y, "delta_mean_level", figure_dir, f"{name}_E_slices_delta",
        signals=E_SLICE_POPS, label="E potential, high Im - low Im"))

    cols = [c for c in SUMMARY_COLUMNS if c in delta.columns]
    path = os.path.join(figure_dir, f"{name}_delta_summary.csv")
    delta[[x, y] + cols].sort_values([y, x, "signal"]).to_csv(path, index=False)
    made.append(path)
    return [m for m in made if m]


def run_vip_dose_response(line_dir, trace_path, figure_dir, mark=None):
    """The Im_strength dose-response and the potential time courses.

    Run once per invocation rather than once per pair: it is a 1-D sweep through the
    centre point, not a property of any parameter plane.
    """
    os.makedirs(figure_dir, exist_ok=True)
    made = []

    if line_dir and os.path.exists(os.path.join(line_dir, "sweep_features.csv")):
        df, _ = pf.load_sweep(line_dir)
        df = df[df["signal_kind"] == "potential"]
        # One line per background level wherever the sweep carries several. This is the
        # figure that answers "does it disinhibit?", because the answer is not a property
        # of the modulatory input at all - it is set by Ib_strength, and the sign of the
        # slope flips between Ib 40 and 50.
        hue = "Ib_strength" if df["Ib_strength"].nunique() > 1 else None
        n_lines = max(df["Ib_strength"].nunique(), 5) if hue else 5

        # question 1: the E potential of every layer and area against increasing VIP drive
        for value, label, log in (
                ("mean_level", "mean E potential", False),
                ("mean_rate", "mean firing rate (Hz)", False),
                ("rate_level", "mean rate / sigmoid ceiling", False),
                ("amplitude", "ongoing fluctuation (SD)", False)):
            if value not in df.columns:
                continue
            made.append(pf.plot_population_grid_lines(
                df, "Im_strength", value, "E", figure_dir, f"vip_dose_E_{value}",
                hue=hue, n_lines=n_lines, label=label, log_scale=log, mark=mark))

        # the whole chain per area, so the mechanism is shown rather than asserted:
        # VIP rises, SST falls, E is released - but only where SST had headroom to lose
        for area in ("A3b", "S1", "S2", "Thal"):
            made.append(pf.plot_population_grid_lines(
                df, "Im_strength", "mean_rate", area, figure_dir,
                f"vip_dose_{area}_rate", hue=hue, n_lines=n_lines,
                label="mean firing rate (Hz)", mark=mark))
            made.append(pf.plot_population_grid_lines(
                df, "Im_strength", "rate_level", area, figure_dir,
                f"vip_dose_{area}_rate_level", hue=hue, n_lines=n_lines,
                label="mean rate / sigmoid ceiling", mark=mark))

        path = os.path.join(figure_dir, "vip_dose_response.csv")
        keep = [c for c in ("Im_strength", "signal", "signal_kind", "mean_level",
                            "mean_rate", "rate_level", "peak_saturation", "amplitude",
                            "regime", "band", "peak_freq") if c in df.columns]
        df[keep].sort_values(["signal", "Im_strength"]).to_csv(path, index=False)
        made.append(path)
    else:
        print("    no Im_strength line sweep on disk - dose-response skipped")

    if trace_path and os.path.exists(trace_path):
        with np.load(trace_path, allow_pickle=False) as tr:
            traces = {k: tr[k] for k in tr.files}
        traces["labels"] = [s.decode() if isinstance(s, bytes) else str(s)
                            for s in traces["labels"]]
        traces["settle_s"] = float(traces["settle_s"])
        # a 1 s window: the whole settled run at this panel size is a solid block, and
        # what these traces are for is seeing the level shift and the rhythm on it
        for area in ("E", "S1"):
            made.append(pf.plot_traces(traces, figure_dir, f"vip_traces_{area}",
                                       area=area, kind="potentials",
                                       label="population potential", window=1.0))
        made.append(pf.plot_traces(traces, figure_dir, "vip_traces_S1_rate",
                                   area="S1", kind="rates",
                                   label="firing rate (Hz)", window=1.0))
    else:
        print("    no trace file on disk - time courses skipped")

    return [m for m in made if m]


def report_delta(delta, pops=("VIP1", "SST1", "PV1", "E1", "E1S2", "E3b", "ThalE")):
    """Print the delta at the centre of the plane, as a sign check.

    Disinhibition has a signature: VIP up, SST down, E up. If E does not move at all then
    either the high Im level is too low or VIP is already saturated at the low one - read
    vip_dose_rate_level before touching any parameter. ThalE should be ~0.
    """
    if delta.empty:
        return
    x, y = delta["x_param"].iloc[0], delta["y_param"].iloc[0]
    xs, ys = np.sort(delta[x].unique()), np.sort(delta[y].unique())
    cell = delta[np.isclose(delta[x], xs[len(xs) // 2])
                 & np.isclose(delta[y], ys[len(ys) // 2])]
    lo, hi = delta["im_low"].iloc[0], delta["im_high"].iloc[0]
    print(f"    delta at the centre cell (Im {lo:g} -> {hi:g}):")
    for pop in pops:
        row = cell[cell["signal"] == pop]
        if row.empty:
            continue
        r = row.iloc[0]
        rate = f"{r['delta_mean_rate']:+8.3f} Hz" if "delta_mean_rate" in r else ""
        print(f"      {pop:<6} potential {r['delta_mean_level']:+10.4f}   "
              f"rate {rate}   {r['regime_lo']} -> {r['regime_hi']}")
