"""
File: step004_2_potentials_oscillation_peakfreq.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-08-31
Description: Script for implementation of analyzing oscillatory behaviour.

The same questions as step004_1, asked of the 33 population potentials instead of the
three ROI dipoles: where does each population oscillate, in which band, at what
frequency, and how large is its alpha envelope relative to the loop-free control.

This costs no extra simulations - run_parameter_sweep scores every population alongside
the dipoles, so this file only reads columns that are already on disk. It is the laminar
detail behind step004_1: an ROI dipole averages its layers together, and a layer that
oscillates out of phase with the others disappears from the dipole entirely.

The measurement rules of step004_1 apply unchanged (network gain spectrum, local band
flanks, never a global 1/f fit) - see that file's docstring.

Parameter combinations are defined in step004_run_model_analysis.py.
"""

import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import plotting_functions_step004 as pf                            # noqa: E402


# (column, cmap, label, kwargs) for the per-population panel grids.
POP_VALUES = [
    ("peak_freq", "magma", "peak frequency (Hz)",
     dict(mask_non_oscillating=True, share_scale="figure")),
    ("alpha_prom", "RdBu_r", "alpha prominence (log10 gain)",
     dict(mask_non_oscillating=True, diverging=True, share_scale="figure")),
    ("sig_alpha_gain", "rocket", "alpha envelope / loop-free control",
     dict(mask_non_oscillating=False, log_scale=True, share_scale="panel")),
    ("amplitude", "viridis", "ongoing fluctuation (SD)",
     dict(mask_non_oscillating=False, share_scale="panel")),
]


def run_potential_oscillation_analysis(sweep_dir, figure_dir, name, areas=pf.AREAS):
    """Every per-population oscillation figure of one 2-D grid sweep."""
    df, _ = pf.load_sweep(sweep_dir)
    df = df[df["signal_kind"] == "potential"]
    if df.empty:
        print(f"    {name}: no population rows - skipped")
        return []

    x, y = df["x_param"].iloc[0], df["y_param"].iloc[0]
    os.makedirs(figure_dir, exist_ok=True)
    made = []

    for area in areas:
        for kind in ("regime", "band"):
            made.append(pf.plot_population_grid_map(
                df, x, y, kind, area, figure_dir, f"{name}_{area}_{kind}",
                categorical=kind))
        for value, cmap, label, kwargs in POP_VALUES:
            if value not in df.columns:
                continue
            made.append(pf.plot_population_grid_map(
                df, x, y, value, area, figure_dir, f"{name}_{area}_{value}",
                cmap=cmap, label=label, **kwargs))
    return [m for m in made if m]
