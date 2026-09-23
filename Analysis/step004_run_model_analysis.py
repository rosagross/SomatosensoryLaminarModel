"""
File: step004_run_model_analysis.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-08-31
Description: Run selected analysis scripts that start specific simulations and create plots.

Analyses:
1. Oscillations (dipole signal)         -> step004_1_dipole_oscillation_peakfreq.py
2. Oscillations (population potentials) -> step004_2_potentials_oscillation_peakfreq.py
3. Firing rates                         -> step004_3_firing_rates.py
4. VIP disinhibition                    -> step004_4_disinhibition.py

Each of those analyses plots its own figures, into its own folder.

This file owns the *parameter ranges*: which pairs of parameters are swept, over what
extent, around which centre. The baseline is Simulations/simulation_parameter.json, read
through somato_model.read_simulation_params and used unchanged; the centre is that file
with --ib-centre and --centre applied on top, so the json stays the single source of
defaults for simulation_main / run_optimization / run_sbi while this suite slices
somewhere else.

The default centre is Ib_strength 60 with coupling_strength 0.5, which is a pair, not two
independent choices. Measured on the current model: at Ib 60 and the file's
coupling_strength of 1, six of 33 populations sit at >= 0.95 of their sigmoid ceiling
(E3b 0.993, E1S2 0.983) and no parameter on any axis moves them; dropping the coupling to
0.5 pins nothing, leaves the excitatory populations at 0.10-0.30 of the ceiling, and
raises the A3b alpha prominence from 0.105 to 0.222. Analysis 4 is centred separately
(--disinh-ib-centre 80, --disinh-centre) because it is a different experiment: at Ib 80
the modulatory input releases the L1 excitatory populations of A1 and S2 while at 60 it
suppresses every excitatory population.

Every run is stimulus-free (run_parameter_sweep switches Iext off entirely), so what is
measured throughout is the model's ongoing dynamics, not its evoked response.

Where things go:

    $SIMDIR/parameter_space/<tag>/grid_<x>_vs_<y>_<centre>/   simulation output, NOT
        sweep_features.csv  sweep_spectra.hdf5  sweep_config.json      timestamped

        <centre> names every parameter moved off the json, e.g. `_Ib60_g0.5`. It is in
        the directory name because a slice only says something about the point it passes
        through, so two centres are two experiments and must not collide on disk.

    <figure root>/<timestamp>/                           figures, one folder per run
        run_config.json
        Analysis1_Oscillations_dipole/
        Analysis2_Oscillations_potentials/
        Analysis3_FiringRates/
        Analysis4_Disinhibition/

The split is deliberate: sweeps are the expensive half and re-plotting is seconds, so
re-running the plotting under a new timestamp must never force a re-simulation. Bare
invocation re-plots what is already on disk; --run-sweep simulates first.

Usage:
    # validate the whole path first (~1 min, one pair, coarse grid)
    python Analysis/step004_run_model_analysis.py --run-sweep --smoke

    # the full set
    python Analysis/step004_run_model_analysis.py --run-sweep

    # somewhere else entirely: any parameter may be centred, not just the background
    python Analysis/step004_run_model_analysis.py --run-sweep \
        --ib-centre 50 --centre coupling_strength=0.7 strength_I=0.8

    # re-plot only the firing rates from sweeps already on disk
    python Analysis/step004_run_model_analysis.py --analysis rates

    # the VIP disinhibition analysis: every grid run twice, once per Im_strength level
    python Analysis/step004_run_model_analysis.py --run-sweep --analysis disinhibition

    # the additional pairs suggested in the plan
    python Analysis/step004_run_model_analysis.py --run-sweep --pairs extra

Requires WDDIR (repository root) and SIMDIR to be set.
"""

import argparse
import json
import os
import subprocess
import sys
from collections import OrderedDict
from datetime import datetime

import numpy as np

WDDIR = os.getenv("WDDIR")
if WDDIR is None:
    raise RuntimeError("WDDIR is not set - it must point at the repository root.")
SIMDIR = os.getenv("SIMDIR", "/data/pt_02989")

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(WDDIR, "Simulations"))
sys.path.append(os.path.join(WDDIR, "Simulations", "model"))

import plotting_functions_step004 as pf                            # noqa: E402
import run_parameter_sweep as rps                                  # noqa: E402
import step004_1_dipole_oscillation_peakfreq as a1                 # noqa: E402
import step004_2_potentials_oscillation_peakfreq as a2             # noqa: E402
import step004_3_firing_rates as a3                                # noqa: E402
import step004_4_disinhibition as a4                               # noqa: E402
from somato_model import read_simulation_params                    # noqa: E402

SWEEP_ROOT = os.path.join(SIMDIR, "parameter_space")
FIGURE_ROOT = "/data/pt_02989/Modelling/Figures/routine_analysis"

# The centre every grid passes through: the defaults, used as-is.
BASELINE = read_simulation_params(WDDIR)

# ── the swept extents ──────────────────────────────────────────────────────────
# Explicit, and NOT run_parameter_sweep.PARAM_RANGES. Those are the global search space
# (coupling_strength 0-20, for instance) and the centre now sits at coupling_strength = 1;
# an axis running to 20 spends almost every point above the saturation threshold, where
# every population is pinned against the sigmoid and no parameter does anything, and the
# whole map comes out one flat colour. These brackets are chosen around the centre
# instead. Widen them once the centre-point report says there is headroom.
RANGES = {
    "coupling_strength": (0.2, 4.0),
    "strength_I": (0.4, 1.0),
    "g_intercortical": (0.0, 2.0),
    "g_thalPOm": (0.0, 3.0),
    # 3-100, not 3-20: the suite is centred at Ib_strength 60 (--ib-centre) and analysis
    # 4 at 80, so a 20 ceiling would neither bracket either centre nor reach the part of
    # the axis where the sign of the VIP effect changes.
    "Ib_strength": (3.0, 100.0),
    "Ib_ratio_E": (0.1, 1.5),
    "Ib_ratio_PV": (0.1, 1.5),
    "Ib_ratio_SST": (0.1, 1.5),
    "Ib_ratio_VIP": (0.1, 1.5),
    "Ib_noise_std": (0.05, 1.5),
    "Ib_noise_tau": (0.004, 0.06),
    "e3b_tau": (2.0, 10.0),
    "e1_tau": (2.0, 10.0),
    "e2_tau": (2.0, 10.0),
    "p_2PVE": (10.0, 40.0),
    "p_4PVE": (10.0, 40.0),
    "g_thal": (0.0, 5.0),
    "sI_thal": (0.0, 1.5),
    "thal_EtoI": (0.0, 20.0),
    "delay_factor": (0.001, 0.008),
    "delay_factor_short": (0.001, 0.008),
    "thal_delay_factor": (0.001, 0.005),
    "Im_strength": (0.0, 20.0),
}

# The pairs asked for in the step004_1 / step004_3 docstrings.
PAIRS = OrderedDict([
    ("gxsI",       ("coupling_strength", "strength_I")),
    ("gxinter",    ("coupling_strength", "g_intercortical")),
    ("POmxinter",  ("g_thalPOm", "g_intercortical")),
    ("POmxIb",     ("g_thalPOm", "Ib_strength")),
    ("POmxg",      ("g_thalPOm", "coupling_strength")),
    ("POmxsI",     ("g_thalPOm", "strength_I")),
    ("IbxPV",      ("Ib_strength", "Ib_ratio_PV")),
    ("IbxSST",     ("Ib_strength", "Ib_ratio_SST")),
    ("IbxVIP",     ("Ib_strength", "Ib_ratio_VIP")),
    # the ratios against each other, at a fixed Ib_strength. The three IbxN planes above
    # move every class at once, so the ratio there is confounded with a global scale
    # change; here the classes not on an axis stay put and the comparison is clean.
    ("PVxSST",     ("Ib_ratio_PV", "Ib_ratio_SST")),
    ("ExPV",       ("Ib_ratio_E",  "Ib_ratio_PV")),
    ("ExSST",      ("Ib_ratio_E",  "Ib_ratio_SST")),
    ("tau3bxtau1", ("e3b_tau", "e1_tau")),
    ("tau3bxtau2", ("e3b_tau", "e2_tau")),
    ("tau1xtau2",  ("e1_tau", "e2_tau")),
    ("gxtau3b",    ("coupling_strength", "e3b_tau")),
    ("p4xp2",      ("p_4PVE", "p_2PVE")),
    ("p4xg",       ("p_4PVE", "coupling_strength")),
    ("p4xsI",      ("p_4PVE", "strength_I")),
    ("p2xg",       ("p_2PVE", "coupling_strength")),
    ("p2xsI",      ("p_2PVE", "strength_I")),
    ("noisexg",    ("Ib_noise_std", "coupling_strength")),
    ("noisexsI",   ("Ib_noise_std", "strength_I")),
])

# Additional planes worth having; run with --pairs extra.
EXTRA_PAIRS = OrderedDict([
    # gain vs. operating point - the best single predictor of whether the model is pinned
    ("gxIb",        ("coupling_strength", "Ib_strength")),
    ("sIxIb",       ("strength_I", "Ib_strength")),
    # the other half of the alpha loop
    ("gthalxsIthal", ("g_thal", "sI_thal")),
    ("gthalxEtoI",  ("g_thal", "thal_EtoI")),
    ("interxsI",    ("g_intercortical", "strength_I")),
    # delays set the period more directly than any gain
    ("thaldelxdel", ("thal_delay_factor", "delay_factor")),
    ("shortdelxtau3b", ("delay_factor_short", "e3b_tau")),
    # the OU corner sits at 1/(2*pi*tau) ~ 9.9 Hz, i.e. inside the alpha band; moving tau
    # is the cleanest way to tell a network alpha from noise leaking into the band
    ("noisetauxg",  ("Ib_noise_tau", "coupling_strength")),
    ("Imxg",        ("Im_strength", "coupling_strength")),
])

# key -> (figure folder, entry point, which sweep set it reads).
#
# The sweep kind matters: analysis 4 needs each grid run twice, once per Im_strength
# level, so it cannot read the same directories as 1-3. Its sweeps carry a `disinh_`
# prefix, and the prefix is load-bearing - a paired CSV holds every (x, y) cell twice,
# which would make _pivot in the other analyses silently keep one arbitrary duplicate.
ANALYSES = OrderedDict([
    ("dipole",     ("Analysis1_Oscillations_dipole", a1.run_dipole_oscillation_analysis,
                    "grid")),
    ("potentials", ("Analysis2_Oscillations_potentials",
                    a2.run_potential_oscillation_analysis, "grid")),
    ("rates",      ("Analysis3_FiringRates", a3.run_firing_rate_analysis, "grid")),
    ("disinhibition", ("Analysis4_Disinhibition", a4.run_disinhibition_analysis,
                       "disinh")),
])

# The 1-D Im_strength sweep and the time courses behind the dose-response figures. Not
# per pair: they are a cut through the centre point, not a property of any plane.
LINE_SWEEP = "disinh_line_Im_strength"
TRACE_FILE = "disinh_traces.npz"

# Range overrides for the disinhibition analysis only, on top of RANGES. Empty now that
# RANGES carries the wide Ib axis; kept as the hook for the next one.
DISINH_RANGES = {}


# Short names for the centre parameters that end up in a sweep directory name. A
# parameter not listed here keeps its full name - the suffix has to be unique and
# readable, not short.
CENTRE_ABBREV = {
    "Ib_strength": "Ib", "coupling_strength": "g", "strength_I": "sI",
    "g_intercortical": "ginter", "g_thal": "gthal", "g_thalPOm": "POm",
    "Im_strength": "Im", "Ib_noise_std": "noise", "Ib_noise_tau": "noisetau",
}


def centre_suffix(centre):
    """The part of a sweep directory name that records where the grid is centred.

    Only parameters actually moved off simulation_parameter.json appear, sorted by name.
    A centre that moves nothing gives "", and one that moves only Ib_strength gives
    "_Ib60" - exactly the spelling the sweeps written before this function existed carry,
    which is what keeps those directories loadable in plot-only mode.
    """
    parts = []
    for name in sorted(centre or {}):
        value = float(centre[name])
        if name in BASELINE and value == float(BASELINE[name]):
            continue
        parts.append(f"_{CENTRE_ABBREV.get(name, name)}{value:g}")
    return "".join(parts)


def fmt_centre(centre):
    """The centre as one readable line, for the log."""
    return ", ".join(f"{k}={float(v):g}" for k, v in sorted((centre or {}).items()))


def sweep_name(pair, kind="grid", centre=None):
    """The sweep directory of one pair.

    The centre goes in the name. The whole premise of the suite is that a slice says
    something about the point it passes through, so a grid centred at coupling_strength
    0.5 and one centred at 1 are different experiments and must not overwrite each other
    on disk. For the disinhibition sweeps this is load-bearing twice over: the sign of
    the E delta depends on the background level as well.
    """
    prefix = "disinh_" if kind == "disinh" else ""
    return f"{prefix}grid_{pair[0]}_vs_{pair[1]}{centre_suffix(centre)}"


def build_centre(args, kind="grid"):
    """The full {parameter: value} centre of one sweep family.

    --ib-centre / --disinh-ib-centre stay the shorthand they always were; --centre /
    --disinh-centre add anything else on top and win where the two overlap.

    A centre value outside RANGES for a parameter that is also swept is rejected rather
    than silently accepted: grid_points lays its axis out over RANGES, so such a centre
    is a point no grid over that axis passes through, and every slice would then be
    reported as centred somewhere it never visits.
    """
    ib = args.ib_centre if kind == "grid" else args.disinh_ib_centre
    spec = args.centre if kind == "grid" else args.disinh_centre
    centre = {"Ib_strength": float(ib)}
    centre.update(rps.parse_center(spec))
    for name, value in centre.items():
        if name not in rps.PARAM_RANGES:
            raise SystemExit(f"unknown centre parameter {name!r}; choose from "
                             f"{sorted(rps.PARAM_RANGES)}")
        if name in RANGES:
            lo, hi = RANGES[name]
            if not lo <= value <= hi:
                raise SystemExit(
                    f"centre {name}={value:g} lies outside RANGES[{name!r}] = "
                    f"({lo:g}, {hi:g}) - no grid over that axis would pass through it")
    return centre


# ── running the sweeps ─────────────────────────────────────────────────────────
def sweep_args(args, mode="grid"):
    """The run_parameter_sweep argument namespace shared by every sweep.

    Built through its own parser so every default (seg_dur, fmin/fmax, subjects, ...)
    comes from one place, and so a new flag there does not silently go missing here.
    """
    argv = [mode,
            "--n", str(args.n),
            "--seeds", str(args.seeds),
            "--n-jobs", str(args.n_jobs),
            "--sim-dur", str(args.sim_dur),
            "--settle", str(args.settle),
            # the noise level of the baseline file, not run_parameter_sweep's own default
            "--noise", str(BASELINE["Ib_noise_std"]),
            "--outdir", args.sweep_root,
            "--tag", args.tag]
    if args.no_dipoles:
        argv.append("--no-dipoles")
    return rps.parse_args(argv)


def run_grids(selected, args):
    """Simulate every selected pair. Returns {name: sweep_dir}.

    The ranges are set on the module directly rather than passed as --range strings:
    that flag is nargs="+" into a single dest, so two --range flags silently drop the
    first and the sweep quietly runs against the module default with no warning anywhere.
    """
    sw_args = sweep_args(args)
    base = rps.make_base_params(sw_args)
    # The centre is set on the base as well as passed to reference_point: on the base so
    # it reaches anything read from there directly, and to reference_point because that
    # is what resolves the thal_* names out of thal_connect and what every theta is
    # built from.
    centre = build_centre(args, "grid")
    base.update({k: v for k, v in centre.items() if k in base})
    out_root = os.path.join(args.sweep_root, args.tag)
    os.makedirs(out_root, exist_ok=True)

    dirs = {}
    for name, pair in selected.items():
        for p in pair:
            if p not in RANGES:
                raise ValueError(f"no range defined for {p!r} (pair {name})")
            rps.PARAM_RANGES[p] = RANGES[p]
        if {"Ib_strength", "Ib_ratio_E"} <= set(pair):
            # the background onto E is Ib_strength * Ib_ratio_E, so a grid over both is
            # constant along a diagonal and reads as a smooth gradient that means nothing
            raise ValueError(
                f"pair {name} sweeps Ib_strength and Ib_ratio_E together; only their "
                f"product reaches the excitatory populations, so the grid is degenerate "
                f"along a diagonal. Vary one and hold the other at the centre.")
        ref = rps.reference_point(base, centre)
        sname = sweep_name(pair, "grid", centre)
        print(f"\n### {name}: {pair[0]} {RANGES[pair[0]]} x {pair[1]} {RANGES[pair[1]]}"
              f", centre {fmt_centre(centre)}")
        rps.run_sweep(sname, rps.grid_points(pair, ref, args.n),
                      base, sw_args, out_root, reference=ref)
        dirs[name] = os.path.join(out_root, sname)
    return dirs


# ── the disinhibition sweeps ───────────────────────────────────────────────────
def disinh_points(pair, ref, n, levels):
    """Every grid point of `pair`, once per Im_strength level.

    Built on rps.grid_points so the two halves land on exactly the same axis values -
    they come from the same np.linspace, which is what lets the delta be taken by an
    exact merge rather than a tolerance join.

    The two runs of a cell also see the *same* background noise: add_background_noise
    reseeds from Ib_noise_seed on every apply_params and make_base_params pins it for the
    sweep, so subtracting them cancels the Ornstein-Uhlenbeck realisation exactly. This is
    the property that makes the delta a paired difference; a per-point seed would turn it
    into a difference of noise draws.
    """
    points = []
    for level in levels:
        tag = "low" if float(level) == float(min(levels)) else "high"
        for theta, extra in rps.grid_points(pair, ref, n):
            theta = dict(theta)
            theta["Im_strength"] = float(level)
            points.append((theta, {**extra, "im_level": tag}))
    return points


def disinh_base(sw_args, args, centre=None):
    """Base parameters for the disinhibition sweeps, at the chosen centre.

    The centre is set on the *base* rather than on each theta so it also reaches the
    pairs that do not sweep those parameters; it is passed to reference_point separately
    (see run_grids) for the names that do not live in the base dict.
    """
    base = rps.make_base_params(sw_args)
    centre = build_centre(args, "disinh") if centre is None else centre
    base.update({k: v for k, v in centre.items() if k in base})
    rps.PARAM_RANGES["Im_strength"] = RANGES["Im_strength"]
    rps.PARAM_RANGES.update(DISINH_RANGES)
    return base


def run_disinh_grids(selected, args):
    """Simulate every selected pair twice, at the low and the high Im level."""
    sw_args = sweep_args(args)
    centre = build_centre(args, "disinh")
    base = disinh_base(sw_args, args, centre)
    out_root = os.path.join(args.sweep_root, args.tag)
    os.makedirs(out_root, exist_ok=True)

    dirs = {}
    for name, pair in selected.items():
        if "Im_strength" in pair:
            print(f"### {name}: skipped - the grid already varies Im_strength")
            continue
        for prm in pair:
            if prm not in RANGES:
                raise ValueError(f"no range defined for {prm!r} (pair {name})")
            rps.PARAM_RANGES[prm] = DISINH_RANGES.get(prm, RANGES[prm])
        ref = rps.reference_point(base, centre)
        sname = sweep_name(pair, "disinh", centre)
        print(f"\n### {name} (disinhibition): {pair[0]} x {pair[1]}, "
              f"Im_strength {args.im_levels}, centre {fmt_centre(centre)}")
        rps.run_sweep(sname, disinh_points(pair, ref, args.n, args.im_levels),
                      base, sw_args, out_root, reference=ref)
        dirs[name] = os.path.join(out_root, sname)
    return dirs


def run_disinh_line(args):
    """The Im_strength dose-response, repeated at each background level.

    One sweep holding every (Im_strength, Ib_strength) combination, so the figures can put
    Ib on the hue and show the dependence directly - which is the whole point: below
    Ib ~ 45 the modulatory input suppresses E and above it releases E.
    """
    sw_args = sweep_args(args, mode="line")
    sw_args.n = args.im_n
    centre = build_centre(args, "disinh")
    base = disinh_base(sw_args, args, centre)
    out_root = os.path.join(args.sweep_root, args.tag)
    ref = rps.reference_point(base, centre)

    points = []
    for ib in args.ib_levels:
        for theta, extra in rps.line_points(["Im_strength"], ref, args.im_n):
            theta = dict(theta)
            theta["Ib_strength"] = float(ib)
            points.append((theta, {**extra, "ib_level": float(ib)}))
    print(f"\n### dose-response: Im_strength {RANGES['Im_strength']} in {args.im_n} "
          f"steps, at Ib_strength {args.ib_levels}")
    rps.run_sweep(LINE_SWEEP, points, base, sw_args, out_root, reference=ref)
    return os.path.join(out_root, LINE_SWEEP)


def run_disinh_traces(args):
    """The potential / rate time courses at a few Im levels.

    The only thing outside run_sweep that simulates. It lives here rather than in
    step004_4 so that file keeps the property the other three analyses have - it reads
    what the sweep wrote and never runs the model - and so re-plotting stays free.
    """
    sw_args = sweep_args(args)
    centre = build_centre(args, "disinh")
    base = disinh_base(sw_args, args, centre)
    model = rps._get_model(base, None)
    labels = list(model.get_population_labels())
    ref = rps.reference_point(base, centre)

    pots, rates = [], []
    for level in args.im_trace_levels:
        signals, rate = rps._simulate_signals(model, {**ref, "Im_strength": float(level)},
                                              base, None, labels)
        pots.append(np.asarray(signals, dtype=np.float32))
        rates.append(np.asarray(rate, dtype=np.float32))

    n_steps = pots[0].shape[1]
    path = os.path.join(args.sweep_root, args.tag, TRACE_FILE)
    np.savez_compressed(
        path,
        levels=np.asarray(args.im_trace_levels, dtype=float),
        t=np.arange(n_steps) * model.step_size,
        potentials=np.stack(pots), rates=np.stack(rates),
        labels=np.array(labels, dtype="S32"),
        settle_s=float(args.settle))
    print(f"    -> {path}  ({len(args.im_trace_levels)} levels, {n_steps} steps)")
    return path


# ── the centre-point report ────────────────────────────────────────────────────
def report_centre(df, x, y, centre=None, pops=("E3b", "E1", "E2", "E1S2", "ThalE")):
    """Print the grid cell nearest the reference, as a saturation check.

    Two uses. It says immediately whether the centre is pinned against the sigmoid - if
    rate_level sits at 0.9-1.0 for the excitatory populations then the loop through them
    is open and no map from this centre will show a parameter doing anything. And every
    pair that shares an axis passes through the same cell, so their numbers here must
    agree; a disagreement beyond ~1% means a grid did not centre where it claims.
    """
    centre = {**BASELINE, **(centre or {})}
    xs, ys = pf._grid_axes(df, x, y)
    xc = min(xs, key=lambda v: abs(float(v) - float(centre.get(x, xs[0]))))
    yc = min(ys, key=lambda v: abs(float(v) - float(centre.get(y, ys[0]))))
    cell = df[(df["signal_kind"] == "potential")
              & np.isclose(df[x], xc) & np.isclose(df[y], yc)]
    if cell.empty:
        return
    print(f"    centre cell {x}={xc:.4g}, {y}={yc:.4g}:")
    for pop in pops:
        row = cell[cell["signal"] == pop]
        if row.empty:
            continue
        r = row.iloc[0]
        print(f"      {pop:<6} rate {r['mean_rate']:7.3f} Hz   "
              f"level {r['rate_level']:.3f}   {r['regime']:<12} {r['band']}")


# ── driver ─────────────────────────────────────────────────────────────────────
def git_revision():
    try:
        return subprocess.check_output(["git", "-C", WDDIR, "rev-parse", "HEAD"],
                                       text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:                                              # noqa: BLE001
        return None


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--analysis", nargs="+", default=["all"],
                   choices=["all"] + list(ANALYSES),
                   help="which analyses to plot (default: all). Several may be named - "
                        "'--analysis dipole potentials rates' is the oscillation and rate "
                        "suite without the disinhibition sweeps, which are ten minutes of "
                        "simulation this does not need")
    p.add_argument("--pairs", nargs="+", default=None,
                   help="pair names, or 'all' / 'extra' (default: all of PAIRS)")
    p.add_argument("--run-sweep", action="store_true",
                   help="simulate first; without it, re-plot what is already on disk")
    p.add_argument("--n", type=int, default=11, help="points per axis")
    p.add_argument("--im-levels", type=float, nargs=2, default=[0.0, 10.0],
                   metavar=("LOW", "HIGH"),
                   help="the two Im_strength levels the E delta is taken between "
                        "(default: 0 10, i.e. modulation off vs twice the baseline of 5)")
    p.add_argument("--im-trace-levels", type=float, nargs="+",
                   default=[0.0, 2.5, 5.0, 10.0, 20.0],
                   help="Im_strength levels for the time-course figures")
    p.add_argument("--im-n", type=int, default=21,
                   help="points in the 1-D Im_strength dose-response sweep")
    p.add_argument("--ib-centre", type=float, default=60.0,
                   help="Ib_strength every grid is centred on (analyses 1-3). The default "
                        "is 60, not simulation_parameter.json's 7. At 60 the model is "
                        "only off its sigmoid ceiling because --centre lowers "
                        "coupling_strength with it; the two defaults belong together")
    p.add_argument("--disinh-ib-centre", type=float, default=80.0,
                   help="Ib_strength the paired disinhibition grids are centred on. "
                        "Separate from --ib-centre because it is a different experiment: "
                        "at 80 the modulatory input releases the L1 excitatory "
                        "populations of A1 and S2, at 60 it suppresses every one of them")
    p.add_argument("--centre", nargs="+", default=["coupling_strength=0.5"],
                   metavar="PARAM=VALUE",
                   help="anything else the grids of analyses 1-3 are centred on, on top "
                        "of --ib-centre. The default lowers the coupling because "
                        "Ib_strength 60 at the committed coupling_strength of 1 pins 6 of "
                        "33 populations against the sigmoid (E3b 0.993, E1S2 0.983), "
                        "where no parameter on any axis does anything; at 0.5 nothing is "
                        "pinned, E sits at 0.10-0.30 of its ceiling and the A3b alpha "
                        "peak is three times as prominent. Pass an empty value list to "
                        "centre on the file defaults instead")
    p.add_argument("--disinh-centre", nargs="+", default=["coupling_strength=0.5"],
                   metavar="PARAM=VALUE",
                   help="the same for the disinhibition sweeps, on top of "
                        "--disinh-ib-centre")
    p.add_argument("--ib-levels", type=float, nargs="+", default=[7.0, 40.0, 60.0, 80.0,
                                                                  100.0],
                   help="background levels the dose-response is repeated at, so the sign "
                        "flip is visible in one figure")
    p.add_argument("--seeds", type=int, default=1,
                   help="noise realisations; the per-population TIME COURSES are "
                        "averaged and every measure is then taken once on the average")
    p.add_argument("--n-jobs", type=int, default=max(os.cpu_count() - 1, 1))
    p.add_argument("--sim-dur", type=float, default=8.0)
    p.add_argument("--settle", type=float, default=2.0)
    p.add_argument("--no-dipoles", action="store_true",
                   help="score population potentials only (skips analysis 1)")
    p.add_argument("--tag", default="step004",
                   help="sweep sub-directory under the sweep root")
    p.add_argument("--sweep-root", default=SWEEP_ROOT)
    p.add_argument("--figure-root", default=FIGURE_ROOT)
    p.add_argument("--timestamp", default=None,
                   help="write into an existing output folder instead of a new one")
    p.add_argument("--smoke", action="store_true",
                   help="one pair, n=5, seeds=1 - validates the path end to end")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if args.smoke:
        args.n, args.seeds = 5, 1
        args.sim_dur, args.settle = 4.0, 1.0
        args.im_n, args.im_trace_levels = 5, [0.0, 5.0, 10.0]
        args.ib_levels = [7.0, 80.0]
        args.tag = f"{args.tag}_smoke"
        selected = OrderedDict([("gxsI", PAIRS["gxsI"])])
    elif args.pairs == ["extra"]:
        selected = EXTRA_PAIRS
    elif args.pairs and args.pairs != ["all"]:
        unknown = [n for n in args.pairs if n not in PAIRS and n not in EXTRA_PAIRS]
        if unknown:
            raise SystemExit(f"unknown pair(s) {unknown}; choose from "
                             f"{list(PAIRS)} or {list(EXTRA_PAIRS)}")
        selected = OrderedDict((n, {**PAIRS, **EXTRA_PAIRS}[n]) for n in args.pairs)
    else:
        selected = PAIRS

    wanted = (list(ANALYSES) if "all" in args.analysis
              else [k for k in ANALYSES if k in args.analysis])
    kinds = {ANALYSES[k][2] for k in wanted}

    grid_centre = build_centre(args, "grid")
    disinh_centre = build_centre(args, "disinh")

    out_root = os.path.join(args.sweep_root, args.tag)
    line_dir = os.path.join(out_root, LINE_SWEEP)
    trace_path = os.path.join(out_root, TRACE_FILE)

    def on_disk(paths, label):
        missing = [n for n, d in paths.items()
                   if not os.path.exists(os.path.join(d, "sweep_features.csv"))]
        if missing:
            raise SystemExit(f"no {label} sweep on disk for {missing} under {out_root} - "
                             f"run with --run-sweep first")
        return paths

    dirs, disinh_dirs = {}, {}
    if "grid" in kinds:
        dirs = (run_grids(selected, args) if args.run_sweep else
                on_disk({n: os.path.join(out_root,
                                              sweep_name(p, "grid", grid_centre))
                         for n, p in selected.items()}, "grid"))
    if "disinh" in kinds:
        if args.run_sweep:
            disinh_dirs = run_disinh_grids(selected, args)
            line_dir = run_disinh_line(args)
            trace_path = run_disinh_traces(args)
        else:
            # disinh_centre, not grid_centre: run_disinh_grids writes these directories
            # with the disinhibition centre in the name, so looking them up with the
            # analyses-1-3 centre never matched and plot-only mode could not find them
            disinh_dirs = on_disk(
                {n: os.path.join(out_root,
                                 sweep_name(p, "disinh", disinh_centre))
                 for n, p in selected.items() if "Im_strength" not in p},
                "disinhibition")

    timestamp = args.timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
    figure_base = os.path.join(args.figure_root, timestamp)
    os.makedirs(figure_base, exist_ok=True)
    pf.set_style()

    print(f"\nfigures -> {figure_base}")
    print(f"centre (analyses 1-3): {fmt_centre(grid_centre)}")
    print(f"centre (analysis 4):   {fmt_centre(disinh_centre)}")
    made = {}
    for name, pair in selected.items():
        print(f"\n--- {name} ({pair[0]} x {pair[1]}) ---")
        if name in dirs:
            report_centre(pf.load_sweep(dirs[name])[0], pair[0], pair[1],
                          centre=grid_centre)
        for key in wanted:
            folder, fn, kind = ANALYSES[key]
            source = dirs if kind == "grid" else disinh_dirs
            if name not in source:
                continue
            if kind == "disinh":
                a4.report_delta(a4.delta_table(pf.load_sweep(source[name])[0],
                                               pair[0], pair[1]))
            centre = disinh_centre if kind == "disinh" else grid_centre
            paths = fn(source[name], os.path.join(figure_base, folder),
                       sweep_name(pair, kind, centre))
            made.setdefault(key, []).extend(paths)
            print(f"    {folder}: {len(paths)} file(s)")

    # the dose-response, once rather than per pair
    if "disinh" in kinds:
        folder = ANALYSES["disinhibition"][0]
        paths = a4.run_vip_dose_response(line_dir, trace_path,
                                         os.path.join(figure_base, folder),
                                         mark=args.im_levels)
        made.setdefault("disinhibition", []).extend(paths)
        print(f"\n    {folder} (dose-response): {len(paths)} file(s)")

    with open(os.path.join(figure_base, "run_config.json"), "w") as f:
        json.dump({"timestamp": timestamp,
                   "pairs": {n: list(p) for n, p in selected.items()},
                   "ranges": {p: list(RANGES[p]) for pair in selected.values()
                              for p in pair},
                   "centre": BASELINE,
                   "n": args.n, "seeds": args.seeds,
                   "seed_averaging": "per-population time courses, measures taken once "
                                     "on the average",
                   "sim_dur": args.sim_dur, "settle": args.settle,
                   "stimulus": "off (ongoing dynamics only)",
                   "sweep_dirs": dirs,
                   "disinh_sweep_dirs": disinh_dirs,
                   "im_levels": args.im_levels,
                   "im_trace_levels": args.im_trace_levels,
                   "im_n": args.im_n,
                   "ib_centre": args.ib_centre,
                   "disinh_ib_centre": args.disinh_ib_centre,
                   "grid_centre": grid_centre,
                   "disinh_centre": disinh_centre,
                   "ib_levels": args.ib_levels,
                   "analyses": wanted,
                   "n_figures": {k: len(v) for k, v in made.items()},
                   "git_revision": git_revision()}, f, indent=2, default=str)
    print(f"\nwrote {sum(len(v) for v in made.values())} file(s) under {figure_base}")
    return figure_base


if __name__ == "__main__":
    main()
