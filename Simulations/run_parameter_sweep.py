"""
File: run_parameter_sweep.py
Author: Rosa Grossmann
Description:
    Explore the model's intrinsic parameter space and record where it oscillates and
    at which frequency.

    Every run here is *stimulus-free*: the external input is switched off entirely
    (Iext_strength = 0 and the onset pushed past the end of the run), so what is
    measured is the network's ongoing dynamics - a fixed point, a noise-driven
    fluctuation, or a rhythm of its own - and not its evoked response. The background
    input (Ib_strength) and its Ornstein-Uhlenbeck noise stay on, because without noise
    almost the whole parameter space settles to a fixed point and only self-sustained
    limit cycles remain visible (see run_optimization.py, "pre-stimulus rejection
    penalties"). Run once with `--noise 1` for the noise-driven picture and once with
    `--noise 0` to see which of those regions are genuine limit cycles. Note that with
    the noise off the loop-free control has no ongoing drive either - what is left of it
    is the ring-down of its own initialisation transient, which still carries a usable
    spectrum for about three quarters of the signals and none at all for the rest
    (those fall back to scoring the raw spectrum, scored_vs = 'flat'). Either way a
    limit cycle's peak dwarfs the background, so the classification holds; but the
    prominences of a noise-free and a noisy sweep rest on different references and
    should not be compared against each other numerically.

    Three sweep modes:
      grid    dense 2-D grids over parameter pairs, all other parameters held at the
              reference point                      -> readable heatmaps
      line    1-D one-at-a-time sweeps of every parameter around the same reference
              point                                -> frequency/prominence sensitivity
      screen  a Latin-hypercube sample over all parameters jointly -> finds oscillatory
              regions that no 2-D slice through the reference point passes through

    `grid` and `line` are only as informative as the point they slice through. The
    reference point defaults to the simulation_parameter.json values, but the model is
    pinned against the sigmoid there over most of the space; `--center PARAM=VALUE ...`
    moves it somewhere the model oscillates (a `screen` run is how you find such a
    point).

    Each sweep writes into <outdir>/<tag>/<sweep_name>/:
      sweep_features.csv   one row per (parameter set x signal): regime, peak
                           frequency, alpha/beta prominence, amplitude
      sweep_spectra.hdf5   the Welch spectra themselves (n_sims x n_signals x n_freqs)
      sweep_config.json    everything needed to reproduce the sweep

    Plotted by Analysis/step003_plot_parameter_space.py.

Usage:
    python Simulations/run_parameter_sweep.py grid   --n 21   --n-jobs 8
    python Simulations/run_parameter_sweep.py line   --n 25   --n-jobs 8
    python Simulations/run_parameter_sweep.py screen --n-samples 2000 --n-jobs 8
    python Simulations/run_parameter_sweep.py grid   --noise 0 --tag noise0
    python Simulations/run_parameter_sweep.py grid   --noise 0 --tag noise0_s2alpha \
        --center-json <a {parameter: value} file>
"""

import argparse
import json
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

WDDIR = os.getenv("WDDIR")
SIMDIR = os.getenv("SIMDIR")
if WDDIR is None:
    raise RuntimeError("WDDIR is not set - it must point at the repository root.")
sys.path.append(os.path.join(WDDIR, "Simulations"))
sys.path.append(os.path.join(WDDIR, "Simulations", "model"))

from somato_model import SomatoModel, read_simulation_params  # noqa: E402
from oscillation_metrics import oscillation_features, welch_spectrum  # noqa: E402


# ── the searchable parameter space ─────────────────────────────────────────────
# name -> (low, high). The gain/tau/delay/probability bounds follow
# run_optimization.SEARCH_SPACE where they overlap, so this sweep and the GA explore
# the same space. `thal_*` are the entries of the thal_connect list (see below) and
# are new here; their upper bounds bracket the defaults [tEE, tEI, tIE, tII, tPOmI] =
# [0, 10, 0, 0, 5].
PARAM_RANGES = {
    # global gains
    "coupling_strength": (0.0, 20.0),      # gE
    "strength_I": (0.4, 0.9),              # gI = coupling_strength * strength_I
    "Ib_strength": (3.0, 20.0),            # tonic background drive
    "g_intercortical": (0.0, 2.0),         # A3b <-> S1 <-> S2 long-range gain
    # thalamo-cortical loop
    "g_thal": (0.0, 5.0),                  # VPM output gain
    "sI_thal": (0.0, 1.5),                 # gIthal = g_thal * sI_thal (reticular)
    "g_thalPOm": (0.0, 5.0),               # POm output gain
    "thal_EtoE": (0.0, 20.0),              # thal_connect[0], VPM -> VPM
    "thal_EtoI": (0.0, 20.0),              # thal_connect[1], VPM -> reticular
    "thal_ItoE": (0.0, 20.0),              # thal_connect[2], reticular -> VPM
    "thal_ItoI": (0.0, 20.0),              # thal_connect[3], reticular -> reticular
    "thal_POmtoI": (0.0, 20.0),            # thal_connect[4], POm -> reticular
    # synaptic time constants (ms)
    "e3b_tau": (2.0, 10.0),
    "e1_tau": (2.0, 10.0),
    "e2_tau": (2.0, 10.0),
    # conduction delays (s)
    "delay_factor": (0.001, 0.008),        # long-range S1<->S2, A3b<->S2
    "delay_factor_short": (0.001, 0.008),  # A3b <-> A1
    "thal_delay_factor": (0.001, 0.005),   # thalamus -> cortex
    # PV connection probabilities (%)
    "p_2PVE": (10.0, 40.0),                # L4 PV <- E
    "p_4PVE": (10.0, 40.0),                # L6 PV <- E
}

# thal_connect is a list in the parameter file, so its entries get their own scalar
# sweep names and are folded back into the list before the model sees them.
THAL_CONNECT_NAMES = ("thal_EtoE", "thal_EtoI", "thal_ItoE", "thal_ItoI", "thal_POmtoI")

# The 2-D slices taken by `grid` mode. Each pair puts two parameters that interact
# mechanistically on the two axes; everything else stays at the reference point.
DEFAULT_PAIRS = [
    ("coupling_strength", "strength_I"),   # the classic E/I plane
    ("coupling_strength", "Ib_strength"),  # gain vs. operating point
    ("strength_I", "Ib_strength"),
    ("e1_tau", "e2_tau"),                  # the two S1 synaptic time constants
    ("e1_tau", "strength_I"),              # time constant vs. inhibition: sets the period
    ("g_thal", "sI_thal"),                 # thalamo-cortical vs. reticular gain
    ("g_thal", "thal_EtoI"),               # the VPM <-> reticular loop
    ("g_thalPOm", "g_intercortical"),      # the two long-range routes
    ("p_2PVE", "p_4PVE"),                  # L4 / L6 PV feedback
    ("thal_delay_factor", "delay_factor"), # loop delays
]

# Populations whose potentials are recorded (all of them) plus the three dipole ROIs.
DIPOLE_ROIS = ("A3b", "A1", "S2")


# ── parameter handling ─────────────────────────────────────────────────────────
def reference_point(base_params, center=None):
    """The value each swept parameter takes when it is not being varied.

    Taken from Simulations/simulation_parameter.json (via `base_params`), so by default
    the 2-D slices and the 1-D lines are all centred on the committed default model.

    `center` moves that point. This matters because the default model is not a
    representative place to slice through: at those values the network is pinned
    against the sigmoid over most of the space, so most slices through it are inert and
    say more about the saturation than about the dynamics. Centring on a parameter set
    that is known to oscillate - e.g. the strongest S2 alpha limit cycle found by a
    `screen` run - puts the slices somewhere the model is actually alive.
    """
    ref = {}
    thal = list(base_params["thal_connect"])
    for name in PARAM_RANGES:
        if name in THAL_CONNECT_NAMES:
            ref[name] = float(thal[THAL_CONNECT_NAMES.index(name)])
        else:
            ref[name] = float(base_params[name])
    for name, value in (center or {}).items():
        if name not in PARAM_RANGES:
            raise ValueError(f"unknown parameter {name!r} in the reference point; "
                             f"choose from {sorted(PARAM_RANGES)}")
        ref[name] = float(value)
    return ref


def parse_center(specs, path=None):
    """{param: value} from `--center PARAM=VALUE ...` and/or a `--center-json` file."""
    center = {}
    if path:
        with open(path) as f:
            center.update(json.load(f))
    for spec in specs or []:
        name, _, value = spec.partition("=")
        if not value:
            raise ValueError(f"--center expects PARAM=VALUE, got {spec!r}")
        center[name] = float(value)
    return center


def expand_theta(theta, base_params):
    """Turn a {sweep name: value} dict into the kwargs SomatoModel.apply_params takes.

    The only translation needed is thal_connect: its five entries are swept as
    separate scalars but the model reads them as one list.
    """
    params = {k: v for k, v in theta.items() if k not in THAL_CONNECT_NAMES}
    if any(k in theta for k in THAL_CONNECT_NAMES):
        thal = list(base_params["thal_connect"])
        for i, name in enumerate(THAL_CONNECT_NAMES):
            if name in theta:
                thal[i] = float(theta[name])
        params["thal_connect"] = thal
    return params


def make_base_params(args):
    """simulation_parameter.json with the stimulus switched off and the run lengthened."""
    params = read_simulation_params(WDDIR)
    params.update({
        "simulation_dur": args.sim_dur,
        # no stimulus at all: zero amplitude, and an onset past the end of the run so
        # create_Iext writes into an empty slice whatever the duration is
        "Iext_strength": 0.0,
        "Iext_duration": 0.0,
        "input_type": "step",
        "input_onset": args.sim_dur + 1.0,
        "Ib_noise_std": args.noise,
        # one fixed noise realisation across the whole sweep, so differences between
        # parameter sets are differences in dynamics and not in the noise draw
        "Ib_noise_seed": args.seed,
        "save_params": False,
        "save_results": False,
        "save_connectivity": False,
        "save_full_potentials": False,
        "plot_rates": False,
        "plot_potentials": False,
        "plot_all_potentials": False,
    })
    return params


# ── the simulator ──────────────────────────────────────────────────────────────
# One model instance per worker process: building it costs ~0.1 s and the dipole
# forward projections ~2 s, both of which are invariant across parameter sets, while
# apply_params + simulate is ~0.9 s for an 8 s run.
_MODEL = None
_MODEL_KEY = None
_NULL_CACHE = {}

# The loop-free control: the same network with every coupling gain scaled down by
# NULL_GAIN_SCALE, so a path of n synapses is attenuated by scale**n. The one-relay
# feedforward path (background -> population -> its target's synapse) survives and sets
# the reference shape; every recurrent loop, needing at least two relays, is gone. Its
# spectrum is what the real run is divided by (oscillation_metrics.network_gain_spectrum).
#
# Scaling rather than zeroing matters for two signals that would otherwise have no
# reference at all:
#   - the ROI dipoles, which are projections of self.potential[:, :-2] and therefore
#     exclude the background column entirely: with the gains at zero the dipole is
#     identically zero, not merely small.
#   - the thalamic populations, which receive no background input at all and in a
#     stimulus-free run are driven only by the cortex.
# Scale-free prominence means the absolute size of the control signal is irrelevant, so
# the attenuation costs nothing.
NULL_GAIN_SCALE = 1e-3
NULL_GAIN_PARAMS = ("coupling_strength", "g_thal", "g_thalPOm", "g_intercortical") \
    + THAL_CONNECT_NAMES
# strength_I and sI_thal are ratios (gI = coupling_strength * strength_I), so they are
# left alone - scaling coupling_strength already scales gI with it.


def _get_model(base_params, subjects):
    global _MODEL, _MODEL_KEY
    key = (json.dumps(base_params, sort_keys=True, default=str), tuple(subjects or ()))
    if _MODEL is None or _MODEL_KEY != key:
        _MODEL = SomatoModel(base_params, WDDIR=WDDIR)
        if subjects:
            # prime the forward-model cache once (~2 s), before the first timed run;
            # self.potential is already allocated by __init__, so no simulation is
            # needed to build the projections
            _MODEL.compute_dipoles(list(subjects))
        _MODEL_KEY = key
    return _MODEL


def _simulate_signals(model, params, base_params, subjects, labels):
    """Run the model once with `params` and return its recorded signals and rates."""
    model.apply_params(expand_theta(params, base_params))
    model.initialize_state()
    model.simulate()
    potentials = np.sum(model.potential, axis=1)               # (33, n_steps)
    signals = [potentials[i] for i in range(len(labels))]
    if subjects:
        # same grouping as SomatoModel._roi_dipoles: A3b, the four A1 layers, the four
        # S2 layers
        sim_dip = model.compute_dipoles(list(subjects))
        signals += [sim_dip[0], sim_dip[1:5].sum(axis=0), sim_dip[5:9].sum(axis=0)]
    return signals, model.rate.copy()


def _null_spectra(model, theta, base_params, subjects, labels, feat_kwargs):
    """Welch spectra of every signal in the loop-free control run of this parameter set.

    The control depends on the gains too (they set which one-relay paths dominate the
    reference), so it is keyed on the whole parameter set and only repeats where a
    sweep genuinely revisits a point - as `line` mode does at the reference point.
    """
    key = tuple(round(float(theta[p]), 12) for p in sorted(theta))
    if key not in _NULL_CACHE:
        null_theta = dict(theta)
        for name in NULL_GAIN_PARAMS:
            null_theta[name] = float(theta[name]) * NULL_GAIN_SCALE
        signals, _ = _simulate_signals(model, null_theta, base_params, subjects, labels)
        _NULL_CACHE[key] = [welch_spectrum(s, model.step_size,
                                           seg_dur=feat_kwargs["seg_dur"],
                                           overlap=feat_kwargs["overlap"],
                                           settle_s=feat_kwargs["settle_s"],
                                           fmin=feat_kwargs["fmin"],
                                           fmax=feat_kwargs["fmax"])[1]
                            for s in signals]
        # a sweep can visit thousands of parameter sets; keep the cache bounded
        if len(_NULL_CACHE) > 64:
            _NULL_CACHE.pop(next(iter(_NULL_CACHE)))
    return _NULL_CACHE[key]


def run_point(sim_id, theta, base_params, subjects, feat_kwargs, extra_cols=None):
    """Simulate one parameter set and score every recorded signal.

    Returns:
        (rows, psd, freqs, signal_names, psd_null) where `rows` is a list of one dict
        per signal and `psd`/`psd_null` are (n_signals, n_freqs). On a failure (blow-up,
        singular matrix, ...) the rows are still returned, marked regime='diverged', so
        a broken corner of the space appears in the map rather than aborting the sweep.
    """
    model = _get_model(base_params, subjects)
    labels = list(model.get_population_labels())
    kinds = ["potential"] * len(labels)
    signal_names = list(labels)
    if subjects:
        signal_names += [f"roi_{r}" for r in DIPOLE_ROIS]
        kinds += ["dipole"] * len(DIPOLE_ROIS)

    common = dict(theta)
    common.update(extra_cols or {})
    common["sim_id"] = sim_id

    try:
        # the loop-free reference for this parameter set first, then the real run
        nulls = _null_spectra(model, theta, base_params, subjects, labels, feat_kwargs)
        signals, rates = _simulate_signals(model, theta, base_params, subjects, labels)
    except Exception as err:                                   # noqa: BLE001
        rows = [{**common, "signal": s, "signal_kind": k, "regime": "diverged",
                 "band": "none", "error": repr(err)}
                for s, k in zip(signal_names, kinds)]
        empty = np.zeros((len(signal_names), 0))
        return rows, empty, np.array([]), signal_names, empty

    # Whether each population's firing rate still moves (rate_drive), and where on its
    # sigmoid it sits (rate_level). A population whose rate is frozen - at the ceiling
    # or silenced - has opened every loop through it, so it is reported as its own
    # regime rather than scored as an oscillation.
    settle_i = int(round(feat_kwargs["settle_s"] / model.step_size))
    m_max = np.maximum(model.sigm[:, 2], np.finfo(float).tiny)
    mean_rates = rates[:, settle_i:].mean(axis=1)
    rate_level = mean_rates / m_max
    rate_drive = rates[:, settle_i:].std(axis=1) / m_max
    roi_pops = {"roi_A3b": ("E3b",),
                "roi_A1": ("E1", "E2", "E3", "E4"),
                "roi_S2": ("E1S2", "E2S2", "E3S2", "E4S2")}
    # an ROI is only pinned once every excitatory population feeding it is
    roi_drive = {roi: rate_drive[[labels.index(p) for p in pops]].max()
                 for roi, pops in roi_pops.items()}

    rows, psds, freqs = [], [], np.array([])
    for name, kind, sig, psd_null in zip(signal_names, kinds, signals, nulls):
        drive = (rate_drive[labels.index(name)] if kind == "potential"
                 else roi_drive[name])
        feats, freqs, psd = oscillation_features(sig, model.step_size,
                                                 rate_drive=float(drive),
                                                 psd_null=psd_null,
                                                 return_psd=True, **feat_kwargs)
        row = {**common, "signal": name, "signal_kind": kind,
               "rate_drive": float(drive), **feats}
        if kind == "potential":
            i = labels.index(name)
            row["mean_rate"] = float(mean_rates[i])
            # peak of the *settled* oscillation, not of the initialisation transient,
            # which overshoots and would make every population look saturated
            row["max_rate"] = float(rates[i, settle_i:].max())
            row["min_rate"] = float(rates[i, settle_i:].min())
            row["rate_level"] = float(rate_level[i])
            # how much of the sigmoid's output range is left above the oscillation
            # peak: 1.0 means the population tops out at its ceiling on every cycle and
            # an external input can no longer raise its rate
            row["peak_saturation"] = float(rates[i, settle_i:].max() / m_max[i])
        rows.append(row)
        psds.append(psd)

    return rows, np.asarray(psds), freqs, signal_names, np.asarray(nulls)


# ── sweep construction ─────────────────────────────────────────────────────────
def grid_points(pair, ref, n):
    """Every (x, y) of a 2-D grid over `pair`, with all other parameters at `ref`."""
    (px, py) = pair
    xs = np.linspace(*PARAM_RANGES[px], n)
    ys = np.linspace(*PARAM_RANGES[py], n)
    points = []
    for y in ys:
        for x in xs:
            theta = dict(ref)
            theta[px], theta[py] = float(x), float(y)
            points.append((theta, {"x_param": px, "y_param": py}))
    return points


def line_points(params, ref, n):
    """One-at-a-time sweeps: each parameter varied across its range in turn."""
    points = []
    for name in params:
        for v in np.linspace(*PARAM_RANGES[name], n):
            theta = dict(ref)
            theta[name] = float(v)
            points.append((theta, {"varied_param": name, "varied_value": float(v)}))
    return points


def screen_points(params, n_samples, seed):
    """Latin-hypercube sample over all `params` jointly."""
    from scipy.stats import qmc

    sampler = qmc.LatinHypercube(d=len(params), seed=seed)
    unit = sampler.random(n_samples)
    lo = np.array([PARAM_RANGES[p][0] for p in params])
    hi = np.array([PARAM_RANGES[p][1] for p in params])
    scaled = qmc.scale(unit, lo, hi)
    return [(dict(zip(params, map(float, row))), {}) for row in scaled]


def run_sweep(name, points, base_params, args, out_root, reference=None):
    """Run every point of one sweep in parallel and write its outputs."""
    out_dir = os.path.join(out_root, name)
    os.makedirs(out_dir, exist_ok=True)
    subjects = None if args.no_dipoles else args.subjects
    feat_kwargs = dict(seg_dur=args.seg_dur, overlap=args.overlap,
                       settle_s=args.settle, fmin=args.fmin, fmax=args.fmax)

    print(f"\n=== {name}: {len(points)} simulations on {args.n_jobs} core(s) ===",
          flush=True)
    t0 = time.time()
    results = Parallel(n_jobs=args.n_jobs, verbose=5, batch_size=8)(
        delayed(run_point)(i, theta, base_params, subjects, feat_kwargs, extra)
        for i, (theta, extra) in enumerate(points)
    )
    elapsed = time.time() - t0
    print(f"    {elapsed/60:.1f} min ({elapsed/max(len(points),1):.2f} s per simulation)",
          flush=True)

    rows = [r for res in results for r in res[0]]
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out_dir, "sweep_features.csv"), index=False)

    # the spectra, for the waterfall plots and for re-scoring without re-simulating.
    # Both the run's own spectrum and the loop-free control it was scored against are
    # kept, so the network gain spectrum can be rebuilt downstream.
    freqs = next((res[2] for res in results if len(res[2])), np.array([]))
    signal_names = results[0][3]
    psd = np.full((len(results), len(signal_names), len(freqs)), np.nan)
    psd_null = np.full_like(psd, np.nan)
    for i, res in enumerate(results):
        if res[1].size and res[1].shape[1] == len(freqs):
            psd[i] = res[1]
        if res[4].size and res[4].shape[1] == len(freqs):
            psd_null[i] = res[4]
    with h5py.File(os.path.join(out_dir, "sweep_spectra.hdf5"), "w") as f:
        f.create_dataset("freqs", data=freqs)
        f.create_dataset("psd", data=psd, compression="gzip")
        f.create_dataset("psd_null", data=psd_null, compression="gzip")
        f.create_dataset("signals", data=np.array(signal_names, dtype="S32"))
        param_names = sorted(points[0][0])
        f.create_dataset("param_names", data=np.array(param_names, dtype="S32"))
        f.create_dataset("theta", data=np.array([[p[0][k] for k in param_names]
                                                 for p in points], dtype=float))

    with open(os.path.join(out_dir, "sweep_config.json"), "w") as f:
        json.dump({"sweep": name, "n_points": len(points), "mode": args.mode,
                   "base_params": {k: v for k, v in base_params.items()
                                   if not isinstance(v, np.ndarray)},
                   # the point every slice passes through; what the plots mark
                   "reference": reference,
                   "feature_kwargs": feat_kwargs,
                   "param_ranges": PARAM_RANGES,
                   "subjects": subjects, "elapsed_s": elapsed}, f, indent=2,
                  default=str)

    n_osc = int((df["regime"] == "oscillation").sum())
    print(f"    -> {out_dir}\n       {n_osc}/{len(df)} signal rows oscillating; "
          f"regimes: {df['regime'].value_counts().to_dict()}", flush=True)
    return df


# ── CLI ────────────────────────────────────────────────────────────────────────
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=("grid", "line", "screen"))
    p.add_argument("--n", type=int, default=21,
                   help="points per axis (grid) or per parameter (line)")
    p.add_argument("--n-samples", type=int, default=2000,
                   help="Latin-hypercube samples (screen mode)")
    p.add_argument("--pairs", nargs="+", default=None,
                   help="grid pairs as 'x:y' (default: the DEFAULT_PAIRS list)")
    p.add_argument("--params", nargs="+", default=None,
                   help="restrict line/screen mode to these parameters")
    p.add_argument("--n-jobs", type=int, default=max(os.cpu_count() - 1, 1))
    p.add_argument("--sim-dur", type=float, default=8.0, help="run length (s)")
    p.add_argument("--settle", type=float, default=2.0,
                   help="seconds skipped before the analysis window (s)")
    p.add_argument("--seg-dur", type=float, default=1.0,
                   help="Welch segment length (s); 1/seg_dur is the frequency resolution")
    p.add_argument("--overlap", type=float, default=0.5)
    p.add_argument("--fmin", type=float, default=1.0)
    p.add_argument("--fmax", type=float, default=45.0)
    p.add_argument("--noise", type=float, default=1.0,
                   help="Ib_noise_std; 0 leaves only self-sustained limit cycles visible")
    p.add_argument("--seed", type=int, default=0, help="noise / sampling seed")
    p.add_argument("--subjects", type=int, nargs="+", default=[15],
                   help="subject IDs whose forward models the dipole ROIs use")
    p.add_argument("--no-dipoles", action="store_true",
                   help="score population potentials only (no DATADIR needed)")
    p.add_argument("--range", nargs="+", default=None, dest="ranges",
                   metavar="PARAM:LO:HI",
                   help="override a parameter's sweep range, e.g. coupling_strength:0:8")
    p.add_argument("--center", nargs="+", default=None, metavar="PARAM=VALUE",
                   help="move the reference point the slices pass through (default: the "
                        "simulation_parameter.json values)")
    p.add_argument("--center-json", default=None,
                   help="JSON file of {parameter: value} for the reference point")
    p.add_argument("--outdir", default=None,
                   help="default: $SIMDIR/parameter_space")
    p.add_argument("--tag", default=None,
                   help="sub-directory name (default: noise<Ib_noise_std>)")
    return p.parse_args(argv)


def apply_range_overrides(specs):
    """Apply `--range PARAM:LO:HI` overrides to PARAM_RANGES (in place)."""
    for spec in specs or []:
        name, lo, hi = spec.split(":")
        if name not in PARAM_RANGES:
            raise ValueError(f"unknown parameter {name!r}; choose from "
                             f"{sorted(PARAM_RANGES)}")
        PARAM_RANGES[name] = (float(lo), float(hi))
        print(f"range override: {name} -> ({lo}, {hi})")


def main(argv=None):
    args = parse_args(argv)
    apply_range_overrides(args.ranges)
    base_params = make_base_params(args)
    center = parse_center(args.center, args.center_json)
    ref = reference_point(base_params, center)

    out_root = args.outdir or os.path.join(SIMDIR or ".", "parameter_space")
    out_root = os.path.join(out_root, args.tag or f"noise{args.noise:g}")
    os.makedirs(out_root, exist_ok=True)
    origin = ("simulation_parameter.json defaults" if not center
              else f"{len(center)} parameter(s) moved off the defaults")
    print(f"reference point ({origin}):\n  "
          + ", ".join(f"{k}={v:g}" for k, v in ref.items()))
    print(f"output root: {out_root}")

    if args.mode == "grid":
        pairs = ([tuple(s.split(":")) for s in args.pairs] if args.pairs
                 else DEFAULT_PAIRS)
        for pair in pairs:
            unknown = [p for p in pair if p not in PARAM_RANGES]
            if unknown:
                raise ValueError(f"unknown parameter(s) {unknown}; "
                                 f"choose from {sorted(PARAM_RANGES)}")
            run_sweep(f"grid_{pair[0]}_vs_{pair[1]}", grid_points(pair, ref, args.n),
                      base_params, args, out_root, reference=ref)
    else:
        params = args.params or list(PARAM_RANGES)
        unknown = [p for p in params if p not in PARAM_RANGES]
        if unknown:
            raise ValueError(f"unknown parameter(s) {unknown}; "
                             f"choose from {sorted(PARAM_RANGES)}")
        if args.mode == "line":
            run_sweep("line_all", line_points(params, ref, args.n),
                      base_params, args, out_root, reference=ref)
        else:
            run_sweep("screen", screen_points(params, args.n_samples, args.seed),
                      base_params, args, out_root, reference=ref)


if __name__ == "__main__":
    main()
