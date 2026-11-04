"""
File: run_optimization_dipole.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-10-13
Description: Fine tuning of the evoked responses (time course error, all three ROIs) starting
from the current settings of simulation_main.py.

Two parts:

1. Dipole lengths (10 parameters of dipole_parameterization.py: L_BASE, F_L5, F_L6, F_L4_A1, F_L4_S2, L_THALE_CTX,
   L_THALPOM_CTX, L_THALE_A3B, L_THALPOM_A3B, L_A3B_CORT). The dipole is linear in the simulated PSPs:
   dip_r(t) = sum_s w[r, s] * P[r, s](t)  with w = length * orientation * relative E count, and the length matrix is built
   from the 10 parameters (layer of the target, source population, source cell type). For fixed F_L5 and F_L6 it is linear in
   8 lengths (L_BASE, L_L4_A1 = L_BASE*F_L4_A1, L_L4_S2, L_THALE_CTX, L_THALPOM_CTX, L_THALE_A3B, L_THALPOM_A3B, L_A3B_CORT).
   Windowing, decimation and baseline subtraction are linear as well, so after one simulation the ROI traces are
   sum_g c_g * X_g(t) with X_g computed once, and the 8 lengths follow from a bounded linear least squares problem (post-stimulus
   window of all three ROIs stacked, target normalised by its global peak), about 20 ms instead of one more simulation.
   F_L5 and F_L6 act on several lengths at once (bilinear) and are GA genes. The error is invariant to a common factor, so only
   ratios are identified: L_BASE is kept at its warm value (2 mm) and the other lengths are relative to it.

2. Optional dynamic parameters (searched with the GA from neuronaldynamics, as in
   run_optimization.py). Enabled in OPTIONAL below or with --enable / --disable.
   Model parameters (delays, f_intra, e3b_tau / e1_tau / e2_tau, ...) and pathway scale factors (multiply a block of W after
   apply_params: ff/fb between A3b, A1, S2 and the POm drive/targets). With at least one of them enabled the GA
   also searches F_L5 and F_L6 (DIPOLE_GENES); without any, only one simulation runs (dipole-only mode).
   The warm start is read from simulation_main.py (module-level assignments, no import).

Results: $SIMDIR/optimization_finetune/ft_<timestamp>_<tag>/

Usage:
    python run_optimization_dipole.py --selftest
    python run_optimization_dipole.py --dipole-only            # one simulation + analytic lengths, F_L5 / F_L6 by Nelder-Mead
    python run_optimization_dipole.py                          # GA over the enabled parameters
    python run_optimization_dipole.py --enable POm_tgt_A1 --disable f_intra --n-iter 3
    python run_optimization_dipole.py --resume ft_20261014_182908_ga-16p              # start from a previous run's optimum
    python run_optimization_dipole.py --resume ft_20261014_182908_ga-16p --disable e1_tau,e2_tau --resume-top 5

--resume FOLDER (a path or a folder name under optimization_finetune/) seeds the GA by parameter name from the best parameters
of that run, so the set of searched parameters may differ from the previous run (more or fewer): parameters searched now and
before start at the old optimum, new ones at the simulation_main.py value, and parameters that are no longer searched stay
fixed at the old optimum (--revert-dropped: back to simulation_main.py). --resume-top K adds the next K-1 best evaluations
of the old run to the initial population.
"""

import argparse
import ast
import json
import os
import sys
import time
import types
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.optimize import nnls, lsq_linear, minimize

# ── paths ──────────────────────────────────────────────────────────────────────
WDDIR  = os.getenv("WDDIR")   # /data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel
SIMDIR = os.getenv("SIMDIR")  # /data/p_02989/Modelling/output_grossmannr
HERE   = os.path.dirname(os.path.abspath(__file__))
SIM_DIR = os.path.dirname(HERE)                               # .../Simulations
OUT_ROOT = os.path.join(SIMDIR, "optimization_finetune")
WARM_START_SOURCE = os.path.join(SIM_DIR, "simulation_main.py")

sys.path.insert(0, SIM_DIR)
sys.path.append(os.path.join(WDDIR, "Simulations", "model"))
sys.path.insert(0, "/data/p_02989/Modelling/neuronaldynamics/src")
# Optimizer.py imports Utils via a bare `from Utils import ...`, so its containing
# directory must be on sys.path (same setup as run_optimization.py).
sys.path.insert(0, "/data/p_02989/Modelling/neuronaldynamics/src/neuronaldynamics")

import matplotlib
sys.path.insert(0, HERE)
import dipole_parameterization as dp                # the 10 dipole parameters -> length matrix
import run_optimization as ro                  # builds the model, target paths, preprocessing
import somato_model as sm
from neuronaldynamics.Optimizers.Optimizer import GA

# neuronaldynamics.Utils forces TkAgg; switch to Agg when headless (as in run_optimization.py)
if os.environ.get("OPT_HEADLESS") or not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

model = ro.model
ROIS = ("A3b", "A1", "S2")

# ── settings ───────────────────────────────────────────────────────────────────
N1, N2, N3, N_ITER = 50, 80, 80, 5     # GA: initial population, crossover / mutation offspring, iterations
PENALTY = 5.0                          # objective value of a failed / non-finite evaluation
R_MAX = 4.0                            # upper bound of a dipole length, in units of its warm value (after the scale of pass 1)
LO_CORTEX = 0.25                       # lower bound of L_BASE and L_A3B_CORT (same units); the other lengths may go to 0.
                                       # Without it the fit can switch the main cortical dipoles off entirely.
POLISH = True                          # polish the weights on the exact error (pre + post, global peak) in the final outputs
POLISH_IN_GA = False                   # ... and inside the GA objective (False: smooth objective for the GA's finite-difference gradients)
POLISH_MAXFEV = 300
FIXED = {}                             # values carried over from a resumed run for parameters that are not searched now

# ── optional dynamic parameters ────────────────────────────────────────────────
# kind "param": model parameter, "path": factor on a block of W [rows, cols] (1 = unchanged).
# bounds: "rel" = factor on the warm start value, "abs" = absolute; "clip" limits the result.
PATHWAYS = sm.PATHWAYS                  # shared with simulation_main.py (somato_model.py)
_PATH_B = dict(kind="path", abs=(0.5, 2.0))
OPTIONAL = {
    # ranked by the scouting scan (raw and profiled error range), see the plan
    "thal_delay_factor":       dict(enabled=True,  kind="param", rel=(0.6, 1.6), clip=(0.001, 0.012)),
    "delay_factor":            dict(enabled=True,  kind="param", rel=(0.5, 1.5), clip=(0.002, 0.012)),
    "f_intra":                 dict(enabled=True,  kind="param", abs=(0.5, 3.0)),
    "ff_A1_S2":                dict(enabled=True,  **_PATH_B),
    "fb_A1_S2":                dict(enabled=True,  **_PATH_B),
    "ff_A3b_A1":               dict(enabled=True, **_PATH_B),
    "fb_A3b_A1":               dict(enabled=True, **_PATH_B),
    "fb_A3b_S2":               dict(enabled=True, **_PATH_B),
    "delay_factor_short":      dict(enabled=True, kind="param", rel=(0.5, 2.0), clip=(0.0005, 0.008)),
    # weak in the scouting scan
    "POm_drive_A3b":           dict(enabled=False, **_PATH_B),
    "POm_drive_A1":            dict(enabled=True, **_PATH_B),
    "POm_drive_S2":            dict(enabled=True, **_PATH_B),
    "POm_tgt_A3b":             dict(enabled=False, **_PATH_B),
    "POm_tgt_A1":              dict(enabled=False, **_PATH_B),
    "POm_tgt_S2":              dict(enabled=False, **_PATH_B),
    "ff_A3b_S2":               dict(enabled=False, **_PATH_B),
    # excitatory synaptic time constants (ms)
    "e3b_tau":                 dict(enabled=True, kind="param", abs=(4.0, 6.0)),
    "e1_tau":                  dict(enabled=True, kind="param", abs=(4.0, 6.0)),
    "e2_tau":                  dict(enabled=True, kind="param", abs=(4.0, 6.0)),
    "g_thalPOm":               dict(enabled=False, kind="param", abs=(0.0, 3.0)),
    "f_intra_S2":              dict(enabled=False, kind="param", abs=(0.5, 3.0)),
    "receptor_thalamus_delay": dict(enabled=False, kind="param", abs=(0.005, 0.03)),
}


# the two layer factors of the dipole parameterization, searched together with the enabled parameters above
# (not passed to the model): L5 and L6 target lengths as multiples of the L2/3 length. F_L5 >= 2 by construction.
DIPOLE_GENES = {
    "F_L5": dict(enabled=True, abs=(dp.F_L5_MIN, 10.0)),
    "F_L6": dict(enabled=True, abs=(1.0, 10.0)),
}


# ── warm start from simulation_main.py ─────────────────────────────────────────
# names of the module-level variables of simulation_main.py -> model parameter
MAIN_NAMES = {
    "coupling_strengths": "coupling_strength", "strength_Is": "strength_I",
    "backgrndI_strengths": "Ib_strength", "input_strengths": "Iext_strength",
    "input_durations": "Iext_duration", "ginters": "g_intercortical", "g_thalPOms": "g_thalPOm",
    "modulI_strengths": "Im_strength",
    "Ib_noise_std": "Ib_noise_std", "Ib_ratio_E": "Ib_ratio_E", "Ib_ratio_PV": "Ib_ratio_PV",
    "Ib_ratio_SST": "Ib_ratio_SST", "Ib_ratio_VIP": "Ib_ratio_VIP",
    "delay_factor": "delay_factor", "delay_factor_short": "delay_factor_short",
    "receptor_thalamus_delay": "receptor_thalamus_delay", "thal_delay_factor": "thal_delay_factor",
    "e3b_tau": "e3b_tau", "e1_tau": "e1_tau", "e2_tau": "e2_tau",
    "p_2PVE": "p_2PVE", "p_4PVE": "p_4PVE", "f_intra": "f_intra", "f_intra_S2": "f_intra_S2",
    "resistance_factor": "resistance_factor",
}
MAIN_PARAMS_KEYS = {"g_thal": "g_thal"}   # module-level `params['g_thal'] = ...`; onset / duration are not taken


def read_warm_start(path=WARM_START_SOURCE):
    """Model parameters set by the module-level assignments of simulation_main.py.

    The file runs simulations on import, so it is parsed instead. Only top-level statements are
    read (assignments inside the loops / `if use_opt` are not); a one-element list is unwrapped,
    a longer one (a sweep) is an error. input_onset / simulation_dur are not taken: this script
    keeps the timing of run_optimization.base_params.
    """
    with open(path) as f:
        tree = ast.parse(f.read())
    warm = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        tgt = node.targets[0]
        if isinstance(tgt, ast.Name) and tgt.id in MAIN_NAMES:
            name = MAIN_NAMES[tgt.id]
        elif (isinstance(tgt, ast.Subscript) and isinstance(tgt.value, ast.Name)
              and tgt.value.id == "params" and isinstance(tgt.slice, ast.Constant)
              and tgt.slice.value in MAIN_PARAMS_KEYS):
            name = MAIN_PARAMS_KEYS[tgt.slice.value]
        else:
            continue
        try:
            value = ast.literal_eval(node.value)
        except ValueError:
            raise SystemExit(f"{path}:{node.lineno}: {name} is not a literal, cannot use it as warm start")
        if isinstance(value, (list, tuple)):
            if len(value) != 1:
                raise SystemExit(f"{path}:{node.lineno}: {name} is a sweep {value}, set one value")
            value = value[0]
        warm[name] = float(value)
    missing = [n for n in set(MAIN_NAMES.values()) | set(MAIN_PARAMS_KEYS.values()) if n not in warm]
    if missing:
        raise SystemExit(f"{path}: no module-level value found for {sorted(missing)}")
    return warm


def search_space(warm):
    """Enabled optional parameters -> (names, bounds array, warm start values dict)."""
    names, bounds, seed = [], [], {}
    for name, spec in OPTIONAL.items():
        if not spec["enabled"]:
            continue
        v0 = 1.0 if spec["kind"] == "path" else warm[name]
        lo, hi = spec["abs"] if "abs" in spec else (v0 * spec["rel"][0], v0 * spec["rel"][1])
        if "clip" in spec:
            lo, hi = max(lo, spec["clip"][0]), min(hi, spec["clip"][1])
        names.append(name)
        bounds.append((lo, hi))
        seed[name] = v0
    if names:                               # dynamic parameters are searched: the layer factors join the GA
        for name, spec in DIPOLE_GENES.items():
            if spec["enabled"]:
                names.append(name)
                bounds.append(spec["abs"])
                seed[name] = float(getattr(dp, name))
    return names, np.array(bounds, dtype=float), seed


# ── simulation ─────────────────────────────────────────────────────────────────
def simulate(warm, searched=None, use_fixed=True):
    """Apply warm start + FIXED (resumed run) + searched values (model parameters and W block scales) and simulate."""
    searched = {**(FIXED if use_fixed else {}), **(searched or {})}
    p = dict(warm)
    scales = {}
    for k, v in searched.items():
        if k in DIPOLE_GENES:
            continue                        # used by the dipole solve, not by the model
        if k in PATHWAYS:
            scales[k] = float(v)
        else:
            p[k] = float(v)
    model.apply_params(p)                   # rebuilds W, so the block scales never accumulate
    for k, f in scales.items():
        rows, cols = PATHWAYS[k]
        model.W[rows, cols] = model.W[rows, cols] * f
    model.initialize_state()
    model.simulate()


# ── dipole basis: 8 lengths (linear) + 2 layer factors ──────────────────────────
# dipole_parameterization.py maps the 10 dipole parameters to the (9, 34) length matrix. For fixed F_L5 and F_L6 it is linear in
# 8 coefficients (dp.COEFF_NAMES, mm), so after one simulation the ROI traces are  sum_g c_g * X_g(t)  and the 8 lengths follow from
# a bounded linear least squares fit. F_L5 and F_L6 multiply L_BASE, L_THALE_CTX and L_THALPOM_CTX jointly (bilinear): they are GA genes.
COEF = dp.COEFF_NAMES
WARM_C = dp.coefficients_from_params()          # warm lengths (mm), = dipole_parameters_flippedPVSST.json
WARM_F = {"F_L5": dp.F_L5, "F_L6": dp.F_L6}
POST = 250                    # index of stimulus onset on the 2 ms / -500..+250 ms axis
T_MIN_MS, T_MAX_MS = -100, 200                 # fitted window (ms)
PRE_START = POST + T_MIN_MS // 2               # first sample of the pre-stimulus term (200)
POST_END = POST + T_MAX_MS // 2                # end of the post-stimulus window (350, exclusive)
_target_cache = {}


def get_target():
    """Smoothed measured time course per ROI (loaded after ro.tc_data_path points at it)."""
    key = ro.tc_data_path
    if key not in _target_cache:
        _target_cache[key] = model.load_target_timecourse(key)
    return _target_cache[key]


def dipole_weights():
    """(9, 34) weights w[r, s] of the dipole rows and the potential indices of the 9 target E populations,
    as model.compute_dipoles uses them (rows: A3b, A1 L2/3-L6, S2 L2/3-L6)."""
    proj = model.build_dipole_projection()
    p = proj["dipoles"][0]
    w = np.stack([np.asarray(p["A3b"])] + [np.asarray(v) for v in p["A1"]] + [np.asarray(v) for v in p["S2"]])
    return w, proj["exc_pops"]


def unit_weights():
    """(9, 34) weight per mm of length: orientation * resistance * relative E count, = w / length of the parameter file.

    The parameterized lengths are non-zero exactly where the file's are (checked), so w = u * length for any 10 parameters."""
    w, exc = dipole_weights()
    L = dp.reference_matrix()
    u = np.zeros_like(w)
    nz = L != 0
    u[nz] = w[nz] / L[nz]
    stray = np.abs(dp.coefficient_basis(WARM_F["F_L5"], WARM_F["F_L6"])).sum(axis=0)[~nz]
    if stray.size and stray.max() > 0:
        raise SystemExit("the parameterization gives a non-zero length where the dipole file has length 0")
    return u, exc


def build_groups(u, exc, F_L5, F_L6):
    """After a simulation: per length coefficient g the 9 x T dipole array z_g (per mm) and its windowed ROI traces
    X_g (dict roi -> 376 samples). sum_g c_g z_g is the model dipole for lengths c (mm), sum_g c_g X_g the ROI traces."""
    B = dp.coefficient_basis(F_L5, F_L6)
    n_t = model.potential.shape[2]
    z_all, x_all = {}, {}
    for g, name in enumerate(COEF):
        wg = u * B[g]
        z = np.zeros((9, n_t))
        for r in range(9):
            if np.any(wg[r]):
                z[r] = wg[r] @ model.potential[exc[r], :-3, :]
        if np.abs(z).max() > 0:
            z_all[name] = z
            x_all[name] = model.compute_timecourse(z)
    return z_all, x_all


def traces_from_lengths(x_all, c):
    """ROI traces (dict, 376 samples) for lengths c (dict name -> mm; missing = 0)."""
    tc = {roi: np.zeros(376) for roi in ROIS}
    for name, x in x_all.items():
        for roi in ROIS:
            tc[roi] = tc[roi] + c.get(name, 0.0) * x[roi]
    return tc


def tc_error(tc, target):
    """err_tc exactly as model.compute_error_timecourse (scaling_factor=None, all ROIs), without its figure.
    Returns (mean error, {roi: (pre, post)})."""
    eps = 1e-10
    sp = {r: tc[r][POST:POST_END] for r in ROIS}
    tp = {r: target[r][POST:POST_END] for r in ROIS}
    sp = {r: sp[r][:tp[r].shape[0]] for r in ROIS}
    gs = max(np.max(np.abs(sp[r])) for r in ROIS) + eps
    gt = max(np.max(np.abs(tp[r])) for r in ROIS) + eps
    parts = {}
    for r in ROIS:
        post = float(np.sqrt(np.mean((sp[r] / gs - tp[r] / gt) ** 2)))
        pre = float(np.sqrt(np.mean((tc[r][PRE_START:POST] / gs) ** 2))) ** 2
        parts[r] = (pre, post)
    return float(np.mean([a + b for a, b in parts.values()])), parts


def solve_lengths(x_all, target, polish=POLISH):
    """Analytic dipole lengths for the coefficients in x_all.

    1. Bounded linear least squares on the stacked post-stimulus windows,  min ||A c - y||,
       A = windowed group traces of the three ROIs, y = target / global target peak.
       Pass 1 (nnls) fixes the scale: sigma = c_L_BASE / warm L_BASE. Pass 2 bounds c in sigma * warm * [lo, R_MAX], with lo =
       LO_CORTEX for L_BASE and L_A3B_CORT (keeps the main cortical dipoles on) and 0 for the others.
    2. Optional: polish on the exact error (includes the pre-stimulus term), never worse than step 1.
    3. The error ignores a common factor: the result is rescaled so that L_BASE = its warm value (mm).
    Returns (c dict in mm, err, parts).
    """
    names = list(x_all)
    L = target[ROIS[0]][POST:POST_END].shape[0]
    gt = max(np.max(np.abs(target[r][POST:POST_END])) for r in ROIS) + 1e-10
    y = np.concatenate([target[r][POST:POST_END] / gt for r in ROIS])
    A = np.zeros((len(ROIS) * L, len(names)))
    for j, g in enumerate(names):
        for i, roi in enumerate(ROIS):
            A[i * L:(i + 1) * L, j] = x_all[g][roi][POST:POST + L]
    warm = np.array([WARM_C[COEF.index(g)] for g in names])
    c = nnls(A, y, maxiter=5000)[0]
    ib = names.index("L_BASE")
    if c[ib] > 0:
        sigma = c[ib] / warm[ib]
    elif np.any(c > 0):
        sigma = np.median(c[c > 0] / warm[c > 0])
    else:
        sigma = 1.0
    lo_frac = np.array([LO_CORTEX if g in ("L_BASE", "L_A3B_CORT") else 0.0 for g in names])
    lo, hi = lo_frac * warm * sigma, R_MAX * warm * sigma
    c = lsq_linear(A, y, bounds=(lo, hi)).x

    def err_of(vec):
        return tc_error(traces_from_lengths(x_all, dict(zip(names, vec))), target)[0]

    e0 = err_of(c)
    if polish:
        r = minimize(err_of, c, method="Powell", bounds=list(zip(lo, hi)),
                     options=dict(maxfev=POLISH_MAXFEV, xtol=1e-3, ftol=1e-6))
        if r.fun < e0:
            c = np.clip(r.x, lo, hi)
    c = c * (warm[ib] / c[ib])
    cdict = dict(zip(names, c))
    err, parts = tc_error(traces_from_lengths(x_all, cdict), target)
    return cdict, err, parts


def p10_from(c, F_L5, F_L6):
    """The 10 dipole parameters from solved lengths (dict) and the two layer factors."""
    return dp.params_from_coefficients([c.get(n, 0.0) for n in COEF], F_L5, F_L6)


# ── objective ──────────────────────────────────────────────────────────────────
class Run:
    """State shared by the GA objective: warm start, unit weights, evaluation log."""

    def __init__(self, warm, outdir):
        self.warm = warm
        self.outdir = outdir
        self.u, self.exc = unit_weights()
        self.n_eval = 0
        self.log_path = os.path.join(outdir, "evaluations.csv")
        self.header_done = False

    def factors(self, searched):
        """(F_L5, F_L6) of a searched dict; the warm value for a gene that is not searched."""
        return tuple(float(searched.get(k, FIXED.get(k, WARM_F[k]))) for k in ("F_L5", "F_L6"))

    def evaluate(self, searched, polish=POLISH_IN_GA):
        """Simulate + analytic lengths. Returns (err, parts, c, x_all, z_all, (F_L5, F_L6))."""
        simulate(self.warm, searched)
        F = self.factors(searched)
        z_all, x_all = build_groups(self.u, self.exc, *F)
        c, err, parts = solve_lengths(x_all, get_target(), polish=polish)
        return err, parts, c, x_all, z_all, F

    def fit_factors(self, searched):
        """Dipole-only mode: F_L5, F_L6 for the current simulation by Nelder-Mead on the cheap inner solve
        (no new simulation). Returns (F_L5, F_L6)."""
        simulate(self.warm, searched)
        target = get_target()
        lo = np.array([DIPOLE_GENES[k]["abs"][0] for k in ("F_L5", "F_L6")])
        hi = np.array([DIPOLE_GENES[k]["abs"][1] for k in ("F_L5", "F_L6")])

        def f(F):
            F = np.clip(F, lo, hi)
            x_all = build_groups(self.u, self.exc, *F)[1]
            return solve_lengths(x_all, target, polish=False)[1]

        best = None
        for start in ([WARM_F["F_L5"], WARM_F["F_L6"]], [2.5, 2.0], [8.0, 6.0], [3.0, 1.5]):
            r = minimize(f, np.clip(start, lo, hi), method="Nelder-Mead", bounds=list(zip(lo, hi)),
                         options=dict(xatol=1e-3, fatol=1e-7, maxiter=80))
            if best is None or r.fun < best.fun:
                best = r
        return tuple(float(v) for v in np.clip(best.x, lo, hi))

    def objective(self, **params):
        t0 = time.time()
        try:
            err, parts, c, _, _, F = self.evaluate(params)
            if not np.isfinite(err):
                raise FloatingPointError("non-finite error")
            p10 = p10_from(c, *F)
            status = "ok"
        except Exception as exc:               # a cliff in the parameter space is a result, not a crash
            err, parts, p10, status = PENALTY, {}, {}, f"{type(exc).__name__}: {exc}"[:120]
        self.n_eval += 1
        row = {"eval": self.n_eval, **params, **p10, "err_tc": err, "status": status,
               **{f"{k}_{r}": v[i] for r, v in parts.items() for i, k in enumerate(("pre", "post"))},
               "seconds": round(time.time() - t0, 2)}
        pd.DataFrame([row]).to_csv(self.log_path, mode="a", header=not self.header_done, index=False)
        self.header_done = True
        print(f"  [{self.n_eval}] {params}  ->  err_tc={err:.4f}  ({status})")
        return err


# ── outputs ────────────────────────────────────────────────────────────────────
def plot_before_after(tc_before, tc_after, target, err_before, err_after, outpath):
    """Target vs warm start (file lengths) vs fine-tuned traces, each normalised by its global post peak."""
    t = np.linspace(-500, 250, 376)
    gt = max(np.max(np.abs(target[r][POST:POST_END])) for r in ROIS)
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.8), sharex=True)
    for ax, roi in zip(axes, ROIS):
        ax.plot(t, target[roi] / gt, color="k", label="measured")
        for tc, label, c in ((tc_before, f"warm start, file lengths ({err_before:.3f})", "C0"),
                             (tc_after, f"fine-tuned ({err_after:.3f})", "C3")):
            gs = max(np.max(np.abs(tc[r][POST:POST_END])) for r in ROIS)
            ax.plot(t, tc[roi] / gs, color=c, label=label)
        ax.axvline(0, color="gray", lw=0.8)
        ax.set_xlim(T_MIN_MS, T_MAX_MS)
        ax.set_title(roi)
        ax.set_xlabel("time (ms)")
    axes[0].set_ylabel("amplitude / global peak")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_outputs(run, best_searched, session_dir, extra_summary=None, fit_factors=False):
    """Re-simulate at the best parameters, write traces, figures, dipole parameters and the summary.
    fit_factors: dipole-only mode, F_L5 / F_L6 are found by Nelder-Mead instead of coming from the GA."""
    target = get_target()
    warm_c = dict(zip(COEF, WARM_C))
    # warm start: dynamics of simulation_main.py; file lengths, then lengths fitted at the warm factors
    simulate(run.warm, use_fixed=False)
    x_warm = build_groups(run.u, run.exc, WARM_F["F_L5"], WARM_F["F_L6"])[1]
    tc_warm = traces_from_lengths(x_warm, warm_c)
    err_before = tc_error(tc_warm, target)[0]
    c_warm, err_warm_w, _ = solve_lengths(x_warm, target, polish=POLISH)
    # best dynamics (and factors)
    searched = dict(best_searched)
    if fit_factors:
        searched["F_L5"], searched["F_L6"] = run.fit_factors(best_searched)
    err, parts, c, x_all, _, F = run.evaluate(searched, polish=POLISH)
    tc_after = traces_from_lengths(x_all, c)
    p10 = p10_from(c, *F)
    err_best_file = tc_error(traces_from_lengths(x_all, warm_c), target)[0]

    sm.TIMECOURSE_DIR = os.path.join(session_dir, "timecourse_comparison")
    os.makedirs(sm.TIMECOURSE_DIR, exist_ok=True)
    model.save_timecourse_comparison(session_dir, tc_after, target, err, filename="best_tc_comparison")
    model.plot_timecourse_comparison(tc_after, target)
    plot_before_after(tc_warm, tc_after, target, err_before, err, os.path.join(session_dir, "warm_vs_finetuned.png"))
    dp.build_params_json(os.path.join(session_dir, "dipole_parameters_finetuned.json"), p10)

    summary = {
        "best_params": {k: float(v) for k, v in best_searched.items()},
        "best_params_full": {**{k: float(v) for k, v in FIXED.items()}, **{k: float(v) for k, v in best_searched.items()},
                             **({"F_L5": float(F[0]), "F_L6": float(F[1])} if fit_factors else {})},
        "dipole_parameters": p10,
        "dipole_parameters_warm_fit": p10_from(c_warm, WARM_F["F_L5"], WARM_F["F_L6"]),
        "dipole_parameters_note": "10 parameters of dipole_parameterization.py; the error is invariant to a common factor, so L_BASE is "
                                  "kept at its warm value and the rest is relative to it",
        "err_tc_warm_start_file_lengths": err_before,
        "err_tc_warm_start_lengths_fitted": err_warm_w,
        "err_tc_best_file_lengths": err_best_file,
        "err_tc_best": err,
        "err_tc_best_per_roi_pre_post": {r: list(v) for r, v in parts.items()},
        **(extra_summary or {}),
    }
    with open(os.path.join(session_dir, "optimization_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  err_tc: warm start {err_before:.4f} -> lengths only {err_warm_w:.4f} -> best {err:.4f}")
    print("  dipole parameters: " + ", ".join(f"{k}={v:.3g}" for k, v in p10.items()))
    return summary


# ── resume from a previous run ─────────────────────────────────────────────────
def resolve_run_dir(name):
    d = name if os.path.isdir(name) else os.path.join(OUT_ROOT, name)
    if not os.path.isdir(d):
        raise SystemExit(f"--resume {name}: no such folder (also tried {d})")
    return d


def load_previous(run_dir, top_k):
    """Previous fine-tune run -> (config, summary, rows). rows: list of {parameter: value}, the best parameters first, then
    up to top_k - 1 further distinct evaluations of evaluations.csv (status ok, lowest err_tc)."""
    paths = {n: os.path.join(run_dir, n) for n in ("run_config.json", "optimization_summary.json", "evaluations.csv")}
    for n, pth in paths.items():
        if not os.path.isfile(pth) and (n != "evaluations.csv" or top_k > 1):
            raise SystemExit(f"--resume: {pth} is missing")
    with open(paths["run_config.json"]) as f:
        cfg = json.load(f)
    with open(paths["optimization_summary.json"]) as f:
        summ = json.load(f)
    best = {k: float(v) for k, v in summ.get("best_params_full", summ["best_params"]).items()}
    rows = [best]
    if top_k > 1:
        ev = pd.read_csv(paths["evaluations.csv"])
        ev = ev[ev["status"] == "ok"].sort_values("err_tc")
        old = [n for n in cfg["searched"] if n in ev.columns]
        for _, r in ev.iterrows():
            row = {**best, **{n: float(r[n]) for n in old}}
            if all(not np.allclose([row[k] for k in old], [x[k] for k in old]) for x in rows):
                rows.append(row)
            if len(rows) >= top_k:
                break
    return cfg, summ, rows


def build_resume(rows, old_cfg, names, seed, warm, revert_dropped):
    """Seed rows (len(rows) x len(names)) for the new search space and the FIXED dict, both by parameter name.
    names/seed: new search space with its default seed (warm start). Returns (population, fixed, new, dropped)."""
    best = rows[0]
    old_searched = set(old_cfg["searched"])
    pop = np.array([[row.get(n, seed[n]) if n in old_searched else seed[n] for n in names] for row in rows], dtype=float)
    new = [n for n in names if n not in old_searched]
    dropped = [n for n in old_searched if n not in names]
    fixed = {}
    for n in best:
        if n in names:
            continue
        if n in old_searched and revert_dropped:
            continue
        fixed[n] = best[n]
    return pop, fixed, new, dropped


def warn_changed_context(old_cfg, warm):
    diff = sorted(k for k in set(warm) | set(old_cfg["warm_start"])
                  if k not in OPTIONAL and not np.isclose(warm.get(k, np.nan), old_cfg["warm_start"].get(k, np.nan), equal_nan=True))
    base_new = json.loads(json.dumps(ro.base_params))
    diff += sorted(f"base_params.{k}" for k in set(base_new) | set(old_cfg["base_params"]) if base_new.get(k) != old_cfg["base_params"].get(k))
    if diff:
        print(f"  WARNING: settings differ from the resumed run (its optimum was found under the old ones): {diff}")


# ── self test ──────────────────────────────────────────────────────────────────
def selftest(warm, tmpdir):
    """The identities the analytic solve relies on."""
    ro.tc_data_path, _ = ro.preprocess_targets(ro.tc_data_path, ro.ps_data_path, tmpdir)
    u, exc = unit_weights()
    simulate(warm)
    d0 = model.compute_dipoles()
    z_all, x_all = build_groups(u, exc, WARM_F["F_L5"], WARM_F["F_L6"])
    warm_c = dict(zip(COEF, WARM_C))
    dsum = sum(warm_c[g] * z for g, z in z_all.items())
    ok_b = np.abs(dsum - d0).max() < 1e-10 * max(1.0, np.abs(d0).max())
    print(f"(b) sum of group dipoles at the warm lengths == compute_dipoles: max diff {np.abs(dsum - d0).max():.2e}  {'PASS' if ok_b else 'FAIL'}")

    target = get_target()
    e_mine, _ = tc_error(traces_from_lengths(x_all, warm_c), target)
    e_ref, _, _ = model.compute_error_timecourse(ro.tc_data_path, d0, rois=ROIS, t_window=(PRE_START, POST_END))
    plt.close("all")
    ok_a = abs(e_mine - e_ref) < 1e-8
    print(f"(a) tc_error at the warm lengths == compute_error_timecourse: {e_mine:.8f} vs {e_ref:.8f}  {'PASS' if ok_a else 'FAIL'}")

    c, e_fit, _ = solve_lengths(x_all, target)
    p10 = p10_from(c, WARM_F["F_L5"], WARM_F["F_L6"])
    print(f"    solved: {', '.join(f'{k}={v:.3g}' for k, v in p10.items())}  err {e_fit:.4f}")
    # (c) the JSON written from the parameters reproduces the error through the model's own dipole code
    jpath = dp.build_params_json(os.path.join(tmpdir, "dipole_parameters_finetuned.json"), p10)
    old_load = model.load_dipole_params

    def _load(self_):
        dpar = sm.read_dipole_params(jpath)
        return dpar["dipole_lengths"], dpar["dipole_orientation"], model.resistance_factor, model.p.get_cellcounts(return_A3b=True)
    try:
        model.load_dipole_params = types.MethodType(_load, model)
        d_json = model.compute_dipoles()
    finally:
        model.load_dipole_params = old_load
    e_json, _ = tc_error(model.compute_timecourse(d_json), target)
    ok_c = abs(e_json - e_fit) < 1e-6
    print(f"(c) error from the fine-tuned JSON {e_json:.8f} vs in-memory {e_fit:.8f}  {'PASS' if ok_c else 'FAIL'}")

    # (d) coefficient basis == build_lengths for random parameters
    try:
        dp.check_basis()
        ok_d = True
    except AssertionError as exc:
        print(f"    {exc}")
        ok_d = False
    print(f"(d) coefficient basis reproduces build_lengths  {'PASS' if ok_d else 'FAIL'}")

    # (e) solved lengths: non-negative, L_BASE at its warm value, factors valid
    ok_e = all(v >= 0 for v in c.values()) and np.isclose(c["L_BASE"], WARM_C[0]) and p10["F_L5"] >= dp.F_L5_MIN
    print(f"(e) lengths >= 0, L_BASE fixed, F_L5 >= {dp.F_L5_MIN}  {'PASS' if ok_e else 'FAIL'}")
    return ok_a and ok_b and ok_c and ok_d and ok_e


# ── run ────────────────────────────────────────────────────────────────────────
def main():
    global N1, N2, N3, N_ITER
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--dipole-only", action="store_true", help="one simulation at the warm start + analytic weights")
    ap.add_argument("--enable", default="", help="comma-separated optional parameters to switch on")
    ap.add_argument("--disable", default="", help="comma-separated optional parameters to switch off")
    ap.add_argument("--resume", default="", help="folder (path or name under optimization_finetune/) of a previous run to start from")
    ap.add_argument("--resume-top", type=int, default=1, help="seed the population with this many best evaluations of the resumed run")
    ap.add_argument("--revert-dropped", action="store_true",
                    help="parameters searched in the resumed run but not now go back to the simulation_main.py value (default: stay at the old optimum)")
    ap.add_argument("--n-iter", type=int, default=N_ITER)
    ap.add_argument("--n1", type=int, default=N1)
    ap.add_argument("--n2", type=int, default=N2)
    ap.add_argument("--n3", type=int, default=N3)
    args = ap.parse_args()
    N1, N2, N3, N_ITER = args.n1, args.n2, args.n3, args.n_iter

    for flag, value in (("enable", True), ("disable", False)):
        for name in [s.strip() for s in getattr(args, flag).split(",") if s.strip()]:
            table = OPTIONAL if name in OPTIONAL else DIPOLE_GENES
            if name not in table:
                raise SystemExit(f"--{flag} {name}: unknown, choose from {sorted(OPTIONAL) + sorted(DIPOLE_GENES)}")
            table[name]["enabled"] = value

    warm = read_warm_start()
    names, bounds, seed = search_space(warm)
    dipole_only = args.dipole_only or not names

    resume = None
    if args.resume:
        run_dir = resolve_run_dir(args.resume)
        old_cfg, old_summ, rows = load_previous(run_dir, max(1, args.resume_top))
        # dipole-only: nothing is searched, so the whole old optimum is carried over as fixed
        pop, fixed, new, dropped = build_resume(rows, old_cfg, [] if dipole_only else names, seed, warm, args.revert_dropped)
        FIXED.update(fixed)
        warn_changed_context(old_cfg, warm)
        resume = {"resumed_from": run_dir, "resume_top": len(rows), "fixed_from_previous": fixed, "new_params": new,
                  "dropped_params": dropped, "revert_dropped": args.revert_dropped,
                  "err_tc_best_previous": old_summ.get("err_tc_best")}
        print(f"Resuming from {run_dir}: {len(rows)} seed row(s); new {new}; dropped {dropped}; fixed {fixed}")

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = "dipole-only" if dipole_only else f"ga-{len(names)}p"
    if resume:
        tag += "_from-" + os.path.basename(os.path.normpath(resume["resumed_from"])).split("_")[2]
    session_dir = os.path.join(OUT_ROOT, f"ft_{run_id}_{tag}")
    os.makedirs(session_dir, exist_ok=True)
    print(f"Fine-tune session: {session_dir}")

    if args.selftest:
        sys.exit(0 if selftest(warm, session_dir) else 1)

    # smoothed measured time course, as in run_optimization.py
    ro.target_dip = None
    ro.tc_data_path, _ = ro.preprocess_targets(ro.tc_data_path, ro.ps_data_path, session_dir)

    run = Run(warm, session_dir)
    config = {
        "run_id": run_id, "timestamp": datetime.now().isoformat(timespec="seconds"),
        "mode": tag, "warm_start_source": WARM_START_SOURCE, "warm_start": warm,
        "base_params": json.loads(json.dumps(ro.base_params)),
        "optional": {k: v for k, v in OPTIONAL.items()},
        "searched": names, "bounds": bounds.tolist(), "seed": seed,
        "ga": {"N1": N1, "N2": N2, "N3": N3, "n_iter": N_ITER},
        "dipole_genes": DIPOLE_GENES, "dipole_parameters_warm": dp.current_params(),
        "dipole_lengths_warm": dict(zip(COEF, WARM_C.tolist())), "fit_window_ms": [T_MIN_MS, T_MAX_MS], "r_max": R_MAX, "lo_cortex": LO_CORTEX,
        "polish": POLISH, "polish_in_ga": POLISH_IN_GA,
        "tc_target": ro.tc_data_path,
        "resume": resume,
    }
    with open(os.path.join(session_dir, "run_config.json"), "w") as f:
        json.dump(config, f, indent=2)

    if dipole_only:
        write_outputs(run, {}, session_dir, fit_factors=True, extra_summary=resume)
        print(f"  written to {session_dir}")
        return

    opt_config = {
        "model_parameters": names,
        "bounds": bounds,
        "reference": 0.0,
        "simulation": run.objective,
        "op": -1,                       # minimise
        "N1": N1, "N2": N2, "N3": N3, "n_iter": N_ITER,
        "tolerance": 1e-8,              # gradient-search tolerance (conf['gTol'])
        "single_run_tol": 1e-12,        # keep far below any reachable cost so all n_iter run (see run_optimization.py)
        "verbose": 1,
    }
    if resume:
        bounds_w = opt_config["bounds"].copy()
        bounds_w[:, 0] = np.minimum(bounds_w[:, 0], pop.min(axis=0))
        bounds_w[:, 1] = np.maximum(bounds_w[:, 1], pop.max(axis=0))
        opt_config["bounds"] = bounds_w
        opt_config["initial_population"] = pop
        print("Seed from resumed run: " + ", ".join(f"{n}={v:.6g}" for n, v in zip(names, pop[0])))
    else:
        ro._seed_population(opt_config, seed, "Warm start from simulation_main.py")

    ga = GA(opt_config)
    ga.run()
    best = dict(zip(names, ga.optimum))
    print("\n── Optimised parameters ──")
    for k, v in best.items():
        print(f"  {k}: {v:.5g}  (warm start {seed[k]:.5g})")
    print(f"  Best err_tc: {np.sqrt(ga.errors[-1]):.6f}  (cost {ga.errors[-1]:.3e})")

    ga.plot_fit()
    ro.plot_fit_diagnostics(ga, opt_config, session_dir)
    write_outputs(run, best, session_dir, extra_summary={
        "best_combined_error_ga": float(np.sqrt(ga.errors[-1])),
        "error_per_iteration": [float(np.sqrt(e)) for e in ga.errors],
        "n_evaluations": run.n_eval,
        **(resume or {}),
    })
    print(f"  written to {session_dir}")


if __name__ == "__main__":
    main()
