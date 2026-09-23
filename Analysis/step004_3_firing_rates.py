"""
File: step004_3_firing_rates.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-08-31
Description: Script for analyzing the behaviour of firing rates upon parameter changes.

Looking at each population:
1. How high is the firing rate?           -> mean_rate (Hz)
2. How far is it from the sigmoid ceiling? -> rate_level (mean/m_max) and
                                              peak_saturation (max/m_max)
3. Did it change after the input onset?    -> NOT ANSWERED HERE, see below.

Why question 3 is open: every run behind these figures is stimulus-free by design.
run_parameter_sweep.make_base_params sets Iext_strength = 0 and pushes input_onset past
the end of the run, so what is measured is the network's ongoing dynamics and there is no
onset to compare against. This is deliberate - it is what makes the spectra in step004_1
interpretable - not an oversight. Adding it later is contained: a --stim option on
make_base_params that restores input_onset = 2.001 / simulation_dur = 3, scores the
spectrum on the pre-stimulus window only, and adds a pre/post rate-difference column.

Why the ceiling matters more than it looks: above roughly coupling_strength 4 the S1 and
A3b populations pin against the sigmoid - E at m_max = 31.4 Hz, PV at 166.8 Hz, the rest
silenced at 0 - the loop through them is open, and nothing responds to any parameter at
all. rate_level and peak_saturation are how that is spotted; a map that looks flat is
usually a map of saturated cells, not of a parameter that does nothing.

Parameter combinations are defined in step004_run_model_analysis.py.
"""

import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import plotting_functions_step004 as pf                            # noqa: E402


# (column, cmap, label, kwargs) for the per-population rate panel grids.
#
# The share_scale choices are not cosmetic. mean_rate needs 'panel': at one operating
# point E1 sits near 27 Hz and E2 near 2 Hz, and the ceilings differ by 5x between E and
# PV, so a shared scale flattens most panels to one colour. rate_level and
# peak_saturation are 0-1 by construction and want 'figure', because comparing panels is
# the entire point of them.
RATE_VALUES = [
    # log, because the rate spans three orders of magnitude within a single panel: at the
    # low-inhibition / high-coupling corner PV2 reaches ~78 Hz against ~0.5 Hz over the
    # rest of the plane, and that corner is real structure rather than a rail to clip
    # away (its peak_saturation is only ~0.5). On a linear scale it turns every other
    # cell black.
    ("mean_rate", "magma", "baseline firing rate (Hz)",
     dict(share_scale="panel", log_scale=True)),
    ("rate_level", "viridis", "mean rate / sigmoid ceiling",
     dict(share_scale="figure", vmin=0.0, vmax=1.0)),
    ("peak_saturation", "viridis", "peak rate / sigmoid ceiling",
     dict(share_scale="figure", vmin=0.0, vmax=1.0)),
    ("rate_drive", "mako", "rate fluctuation / sigmoid ceiling",
     dict(share_scale="panel", log_scale=True)),
    ("alpha_gain", "rocket", "rate alpha / loop-free control",
     dict(share_scale="panel", log_scale=True)),
]

SUMMARY_COLUMNS = ["signal", "mean_rate", "max_rate", "min_rate", "rate_level",
                   "peak_saturation", "rate_drive", "alpha_amp", "alpha_gain",
                   "regime", "band", "peak_freq", "amplitude"]


def run_firing_rate_analysis(sweep_dir, figure_dir, name, areas=pf.AREAS):
    """Every firing-rate figure of one 2-D grid sweep, plus a summary table."""
    df, _ = pf.load_sweep(sweep_dir)
    df = df[df["signal_kind"] == "potential"]
    if df.empty:
        print(f"    {name}: no population rows - skipped")
        return []

    x, y = df["x_param"].iloc[0], df["y_param"].iloc[0]
    os.makedirs(figure_dir, exist_ok=True)
    made = []

    for area in areas:
        # which cells are pinned rather than merely quiet - read this alongside the rates
        made.append(pf.plot_population_grid_map(
            df, x, y, "regime", area, figure_dir, f"{name}_{area}_regime",
            categorical="regime"))
        for value, cmap, label, kwargs in RATE_VALUES:
            if value not in df.columns:
                continue
            made.append(pf.plot_population_grid_map(
                df, x, y, value, area, figure_dir, f"{name}_{area}_{value}",
                cmap=cmap, label=label, mask_non_oscillating=False, **kwargs))

    # line companions, over the excitatory populations of the three ROIs
    e_pops = ["E3b", "E1", "E1S2"]
    for value, label in (("mean_rate", "baseline firing rate (Hz)"),
                         ("rate_level", "mean rate / sigmoid ceiling")):
        made.append(pf.plot_grid_slices(df, x, y, value, figure_dir,
                                        f"{name}_slices_{value}", signals=e_pops,
                                        label=label))

    cols = [c for c in SUMMARY_COLUMNS if c in df.columns]
    summary = df[[x, y] + cols].sort_values([y, x, "signal"])
    path = os.path.join(figure_dir, f"{name}_rate_summary.csv")
    summary.to_csv(path, index=False)
    made.append(path)
    return [m for m in made if m]
