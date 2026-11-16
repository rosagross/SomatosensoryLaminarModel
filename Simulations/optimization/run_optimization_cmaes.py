"""
File: run_optimization_cmaes.py
Description: Fit the model parameters to the target ERP with CMA-ES.

    CMA-ES searches the enabled model parameters of PARAMETER_DICT. The dipole lengths are not searched:
    after the simulations the 8 lengths are fitted to the target without simulating again (F_L5 and F_L6
    are fixed). The populations run in parallel, one model per process.

    Every parameter set is simulated len(SIMS) times, with different stimulus times and noise seeds.
    The ERP error is computed on the average trace, the penalties are averaged over the simulations.

    error = ERP error + baseline rate penalty + thalamus saturation penalty

    Usage:
        python run_optimization_cmaes.py
        SEED=2 MAX_EVALS=10000 N_PROC=24 python run_optimization_cmaes.py
"""

import json
import multiprocessing as mp
import os
import sys
import time
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import cma
from scipy.optimize import minimize

# ── paths ──────────────────────────────────────────────────────────────────────
WDDIR  = os.getenv("WDDIR")   # /data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel
SIMDIR = os.getenv("SIMDIR")  # /data/p_02989/Modelling/output_grossmannr
RESDIR = os.getenv("RESDIR")  # /data/p_02989/shared_workspace/results_grossmannr

sys.path.append(os.path.join(WDDIR, "Simulations", "model"))
import dipole_parameterization as dp
import somato_model as sm

diag_dir = os.path.join(SIMDIR, "optimization_timecourse")
tc_data_path = os.path.join(
    RESDIR, "Figures", "Main", "eeg_results", "source_reconstruction", "group", "_preprestim_corrected",
    "roi_epochswise", "group_roi_timecourse_pooled_ses-elec_preprestim_corrected.csv")

FIT_ROIS = ("A3b", "A1", "S2")

# ── settings ───────────────────────────────────────────────────────────────────
POPSIZE   = 40                                    # simulations per CMA-ES generation
SIGMA0    = 0.3                                   # initial step size, in units of the search range
MAX_EVALS = int(os.getenv("MAX_EVALS", 20000))    # simulations in total
N_PROC    = int(os.getenv("N_PROC", 20))          # parallel processes
SEED      = int(os.getenv("SEED", 1))             # CMA-ES seed (also the random start if START_FROM is empty)

# (stimulus time in s, noise seed) of the simulations of one parameter set
SIMS = [(6.001 + 0.04 * k, k + 1) for k in range(5)]

MID_RATE_WEIGHT = 1           # the penalty of one population is 0 inside RATE_BAND and 1 at a rate of 0 or of the maximum
RATE_BAND = (0.2, 0.6)        # tolerated baseline rate, as a fraction of the maximum rate
SATURATION_THRESHOLD = 0.95   # ThalE rate (fraction of the maximum) above which the stimulus response is penalised

# Session folder (or path) of an earlier run whose best parameters are the start; "" = random start.
START_FROM = ""

# ── parameters ─────────────────────────────────────────────────────────────────
# enabled=True: searched by CMA-ES within `abs`. enabled=False: fixed at `default`.
# Dipole parameters: F_L5 and F_L6 are fixed at `default`. The other 8 are fitted after every
# simulation within `abs`, whatever their `enabled` flag says.
PARAMETER_DICT = {
    # DELAYS
    "thal_delay_factor":       dict(enabled=False, abs=(0.005, 0.01), default=0.007),
    "delay_factor":            dict(enabled=False, abs=(0.005, 0.015), default=0.011),
    "delay_factor_short":      dict(enabled=False, abs=(0.001, 0.006), default=0.0015),
    "receptor_thalamus_delay": dict(enabled=False, abs=(0.005, 0.025), default=0.01),

    # CONNECTIVITY & COUPLING
    "coupling_strength":       dict(enabled=True, abs=(0.5, 4), default=6.85925),
    "strength_I":              dict(enabled=True, abs=(0.5, 2.5), default=0.6),
    "g_intercortical":         dict(enabled=True, abs=(0.5, 3), default=1.75177),
    # intra-cortical
    "f_intra_S2":              dict(enabled=True, abs=(0.8, 2.0), default=1),
    "f_intra":                 dict(enabled=True, abs=(0.8, 2.0), default=1),
    "p_2PVE":                  dict(enabled=False, abs=(37.4, 37.4), default=37.6),
    "p_4PVE":                  dict(enabled=False, abs=(39.6, 39.6), default=39.4),
    # feedforward
    "ff_A1_S2":                dict(enabled=True, abs=(0.5, 5.0), default=1),
    "ff_A3b_S2":               dict(enabled=True, abs=(0.5, 2.0), default=1),
    "ff_A3b_A1":               dict(enabled=True, abs=(0.5, 2.0), default=1),
    # feedback
    "fb_A1_S2":                dict(enabled=True, abs=(0.5, 5.0), default=1),
    "fb_A3b_A1":               dict(enabled=True, abs=(0.5, 5.0), default=1),
    "fb_A3b_S2":               dict(enabled=True, abs=(0.5, 5.0), default=1),
    # thalamus
    "g_thalPOm":               dict(enabled=True, abs=(0.5, 3.0), default=1),
    "g_thal":                  dict(enabled=False, abs=(0.5, 3.0), default=2),
    "POm_drive_A3b":           dict(enabled=False, abs=(0.5, 2.0), default=1),
    "POm_drive_A1":            dict(enabled=False, abs=(0.5, 5.0), default=1),
    "POm_drive_S2":            dict(enabled=False, abs=(0.5, 2.0), default=1),
    "POm_tgt_A3b":             dict(enabled=True, abs=(0.1, 5.0), default=1),
    "POm_tgt_A1":              dict(enabled=True, abs=(0.1, 5.0), default=1),
    "POm_tgt_S2":              dict(enabled=True, abs=(0.1, 5.0), default=1),

    # TIME CONSTANTS
    "e3b_tau":                 dict(enabled=True, abs=(2.5, 10), default=6),
    "e1_tau":                  dict(enabled=True, abs=(4.0, 10), default=6),
    "e2_tau":                  dict(enabled=False, abs=(4.0, 10), default=5.54064),
    # NMDA decay (s) of E -> cell class connections, by postsynaptic class
    "nmda_tau_E":              dict(enabled=True, abs=(5e-3, 20e-3), default=10e-3),
    "nmda_tau_PV":             dict(enabled=True, abs=(2e-3, 20e-3), default=5e-3),
    "nmda_tau_SST":            dict(enabled=False, abs=(5e-3, 40e-3), default=15e-3),
    "nmda_tau_VIP":            dict(enabled=False, abs=(5e-3, 40e-3), default=10e-3),

    # INPUTS
    "Ib_strength":             dict(enabled=True, abs=(5, 100), default=9.61432),
    "Iext_strength":           dict(enabled=True, abs=(1, 40), default=35),
    "Iext_duration":           dict(enabled=False, abs=(0.005, 0.02), default=0.013),
    "Im_strength":             dict(enabled=False, abs=(0, 10), default=7.35317),
    "Ib_noise_std":            dict(enabled=False, abs=(0.5, 2), default=0.5),
    "Ib_ratio_E":              dict(enabled=False, abs=(1, 1), default=1),
    "Ib_ratio_PV":             dict(enabled=True, abs=(0.5, 1.2), default=0.845073),
    "Ib_ratio_SST":            dict(enabled=True, abs=(0.5, 2.0), default=1.19679),
    "Ib_ratio_VIP":            dict(enabled=False, abs=(0.4, 0.7), default=0.4),

    # DIPOLE LENGTHS
    "F_L5":          dict(enabled=False, abs=(dp.F_L5_MIN, 20.0), default=13),   # fixed
    "F_L6":          dict(enabled=False, abs=(1.0, 10.0), default=9),            # fixed
    "L_BASE":        dict(enabled=True, abs=(0.5, 15.0), default=2),
    "F_L4_A1":       dict(enabled=True, abs=(0.01, 0.5),  default=dp.F_L4_A1),
    "F_L4_S2":       dict(enabled=True, abs=(0.01, 0.5),  default=dp.F_L4_S2),
    "L_THALE_CTX":   dict(enabled=False, abs=(0.5, 5.0), default=dp.L_THALE_CTX),
    "L_THALPOM_CTX": dict(enabled=False, abs=(0.5, 5.0), default=dp.L_THALPOM_CTX),
    "L_THALE_A3B":   dict(enabled=True, abs=(0.02, 1.0),  default=dp.L_THALE_A3B),
    "L_THALPOM_A3B": dict(enabled=False, abs=(0.5, 5.0), default=4),
    "L_A3B_CORT":    dict(enabled=False, abs=(0.5, 5.0), default=dp.L_A3B_CORT),
}

F_L5, F_L6 = PARAMETER_DICT["F_L5"]["default"], PARAMETER_DICT["F_L6"]["default"]
SOLVED   = [n for n in dp.PARAM_NAMES if n not in ("F_L5", "F_L6")]                           # the 8 fitted dipole parameters
SEARCHED = [n for n, s in PARAMETER_DICT.items() if s["enabled"] and n not in dp.PARAM_NAMES]  # CMA-ES parameters
FIXED    = {n: s["default"] for n, s in PARAMETER_DICT.items()                                 # fixed model parameters
            if not s["enabled"] and n not in dp.PARAM_NAMES}
LO, HI = np.array([PARAMETER_DICT[n]["abs"] for n in SEARCHED], dtype=float).T
DIP_LO, DIP_HI = np.log(np.array([PARAMETER_DICT[n]["abs"] for n in SOLVED], dtype=float)).T   # dipole fit runs on log lengths
DIP_START = np.log([PARAMETER_DICT[n]["default"] for n in SOLVED])

# ── base simulation settings ───────────────────────────────────────────────────
base_params = {
    "step_size": 0.001,
    "resolution_tstep": 0.001,
    "input_type": "step",
    "input_onset": SIMS[0][0],                    # set per simulation, see SIMS
    "simulation_dur": 6.5,                        # at least the last stimulus time + 0.3 s
    "thal_connect": [0, 1000, 300, 0, 300, 1000],
    "extI_cellcounts": 1000,
    "mI_cellcounts": 2000,
    "bI_cellcounts": 100,
    "thalE_cellcounts": 500,
    "thalI_cellcounts": 500,
    "pom_cellcounts": 500,
    "pomI_cellcounts": 500,
    "resistance_factor": 1,
}



# ── one model per process ──────────────────────────────────────────────────────
def init_worker():
    """Build the model, load the target and the dipole basis. Runs once in every process."""
    global model, TARGET, G, EXC
    model = sm.SomatoModel(base_params)
    target = model.load_target_timecourse(tc_data_path)
    TARGET = np.array([target[r] for r in FIT_ROIS])                     # (3, 376)

    # dipole weight per mm of length (9, 34): the model's projection with every length set to 1
    _, orientation, resistance, cellcounts = model.load_dipole_params_parameterized()
    ones = {"A3b": np.ones(34), "A1": np.ones((4, 34)), "S2": np.ones((4, 34))}
    model.load_dipole_params = lambda: (ones, orientation, resistance, cellcounts)
    proj = model.build_dipole_projection()
    model.load_dipole_params = model.load_dipole_params_parameterized
    d = proj["dipoles"][0]
    unit = np.vstack([d["A3b"]] + list(d["A1"]) + list(d["S2"]))
    EXC = proj["exc_pops"]                                               # the 9 E populations that carry a dipole
    # for fixed F_L5, F_L6 the dipole weights are  sum_g coefficient_g * G[g]  (8 coefficients, dp.COEFF_NAMES)
    G = unit * dp.coefficient_basis(F_L5, F_L6)                          # (8, 9, 34)


# ── objective ──────────────────────────────────────────────────────────────────
def tc_error(tc):
    """ERP error of the traces tc (3, 376), same arithmetic as SomatoModel.compute_error_timecourse."""
    sim_post, tgt_post = tc[:, 250:], TARGET[:, 250:]
    sim_peak = np.max(np.abs(sim_post)) + 1e-10
    tgt_peak = np.max(np.abs(tgt_post)) + 1e-10
    err_post = np.sqrt(np.mean((sim_post / sim_peak - tgt_post / tgt_peak) ** 2, axis=1))
    err_pre = np.mean((tc[:, :250] / sim_peak) ** 2, axis=1)
    return float(np.mean(err_pre + err_post))


def basis_traces():
    """ROI traces per unit coefficient (8, 3, 376) of the simulation the model holds: the traces for any lengths are  coefficients @ X."""
    z = np.einsum("grs,rst->grt", G, model.potential[EXC, :-3, :])                       # (8, 9, n_steps)
    return np.array([[model.compute_timecourse(zg)[r] for r in FIT_ROIS] for zg in z])


def solve_dipoles(X):
    """Fit the 8 dipole lengths to the target for the basis traces X. Returns (10 dipole parameters, ERP error)."""

    def dipole_params(v):
        return {**dict(zip(SOLVED, np.exp(v))), "F_L5": F_L5, "F_L6": F_L6}

    def err(v):
        return tc_error(np.tensordot(dp.coefficients_from_params(dipole_params(v)), X, axes=1))

    res = minimize(err, DIP_START, method="Powell", bounds=list(zip(DIP_LO, DIP_HI)), options=dict(xtol=1e-3, ftol=1e-6))
    return dipole_params(np.clip(res.x, DIP_LO, DIP_HI)), float(res.fun)


def rate_penalties():
    """(baseline rate penalty, thalamus saturation penalty) of the simulation the model holds."""
    stim_idx = round(model.input_onset / model.step_size)

    # baseline: rate fraction of every cortical population at every step of the last 500 ms before the stimulus
    window = round(0.5 / model.step_size)
    frac = model.rate[:30, stim_idx - window:stim_idx] / model.sigm[:30, 2, np.newaxis]
    lo, hi = RATE_BAND
    mid_rate_error = MID_RATE_WEIGHT * np.mean(np.maximum(lo - frac, 0) / lo + np.maximum(frac - hi, 0) / (1 - hi))

    # ThalE must not saturate during the stimulus
    stim_duration = round(model.Iext_duration / model.step_size)
    stim_fraction = model.rate[30, stim_idx:stim_idx + stim_duration] / model.sigm[30, 2]
    saturation_error = np.mean(np.maximum(stim_fraction - SATURATION_THRESHOLD, 0) ** 2)
    return float(mid_rate_error), float(saturation_error)


def evaluate(x):
    """The simulations for the normalised parameter vector x (0..1 per searched parameter). Returns the log row."""
    params = dict(zip(SEARCHED, LO + np.asarray(x) * (HI - LO)))
    X, penalties = [], []
    for onset, seed in SIMS:
        model.apply_params({**FIXED, **params, "input_onset": onset, "Ib_noise_seed": seed})
        model.initialize_state()
        model.simulate()
        X.append(basis_traces())
        penalties.append(rate_penalties())
    mid_rate_error, saturation_error = np.mean(penalties, axis=0)
    dipole, err_erp = solve_dipoles(np.mean(X, axis=0))            # dipole fit and ERP error on the average trace
    error = err_erp + mid_rate_error + saturation_error
    if not np.isfinite(error):
        error = 10.0
    return {**params, **dipole, "err_tc": error, "err_erp": err_erp,
            "mid_rate_error": float(mid_rate_error), "saturation_error": float(saturation_error)}


# ── start values ───────────────────────────────────────────────────────────────
def load_start(run):
    """Values of the searched parameters in the earlier run `run` (its best evaluation)."""
    run_dir = run if os.path.isdir(run) else os.path.join(diag_dir, run)
    with open(os.path.join(run_dir, "run_config.json")) as f:
        cfg = json.load(f)
    summary = os.path.join(run_dir, "optimization_summary.json")
    if os.path.isfile(summary):
        with open(summary) as f:
            best = json.load(f)["best_params"]
    else:   # unfinished run: best row of evaluations.csv
        ev = pd.read_csv(os.path.join(run_dir, "evaluations.csv"))
        best = ev.loc[ev["err_tc"].idxmin(), cfg["optimized_parameters"]].to_dict()
    old = {**cfg["fixed_params"], **best}
    return np.array([float(old[n]) for n in SEARCHED])


# ── run ────────────────────────────────────────────────────────────────────────
def main():
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    session_dir = os.path.join(diag_dir, f"optcma_{run_id}_roi-all")
    os.makedirs(session_dir)
    eval_path = os.path.join(session_dir, "evaluations.csv")
    print(f"Optimization session: {session_dir}")

    # start point in normalised units
    if START_FROM:
        start = load_start(START_FROM)
        for n, v, lo, hi in zip(SEARCHED, start, LO, HI):
            print(f"  start {n} = {v:.6g}" + ("" if lo <= v <= hi else f"  (outside [{lo:g}, {hi:g}]: clipped)"))
        x0 = np.clip((start - LO) / (HI - LO), 0, 1)
    else:
        x0 = np.random.default_rng(SEED).uniform(0.1, 0.9, len(SEARCHED))

    cma_options = {"bounds": [0, 1], "popsize": POPSIZE, "seed": SEED, "maxfevals": MAX_EVALS,
                   "tolfun": 1e-7, "tolx": 1e-6, "verbose": -9}
    with open(os.path.join(session_dir, "run_config.json"), "w") as f:
        json.dump({
            "run_id": run_id,
            "fit_rois": list(FIT_ROIS),
            "optimized_parameters": SEARCHED,
            "bounds": np.column_stack([LO, HI]).tolist(),
            "fixed_params": {**FIXED, "F_L5": F_L5, "F_L6": F_L6, "scaling_factor": 1},
            "solved_dipole_parameters": SOLVED,
            "parameter_defaults": {n: s["default"] for n, s in PARAMETER_DICT.items()},
            "base_params": base_params,
            "sims": SIMS,
            "cma": {**cma_options, "sigma0": SIGMA0},
            "target_csv": tc_data_path,
            "start_from": START_FROM,
            "mid_rate_penalty": True,
            "mid_rate_penalty_weight": MID_RATE_WEIGHT,
            "rate_band": RATE_BAND,
            "saturation_threshold": SATURATION_THRESHOLD,
        }, f, indent=2)

    es = cma.CMAEvolutionStrategy(x0, SIGMA0, cma_options)
    best, best_per_iteration, n_eval, t0 = None, [], 0, time.time()
    with mp.Pool(N_PROC, initializer=init_worker) as pool:
        while not es.stop():
            X = es.ask()
            rows = pool.map(evaluate, X)
            es.tell(X, [r["err_tc"] for r in rows])

            for r in rows:
                n_eval += 1
                r["eval"] = n_eval
                if best is None or r["err_tc"] < best["err_tc"]:
                    best = r
            pd.DataFrame(rows).to_csv(eval_path, mode="a", header=not os.path.exists(eval_path), index=False)
            best_per_iteration.append(best["err_tc"])
            print(f"iteration {es.countiter:4d}  evaluations {n_eval:6d}  best {best['err_tc']:.4f}  "
                  f"(ERP {best['err_erp']:.4f}, rate {best['mid_rate_error']:.4f}, thalamus {best['saturation_error']:.4f})  "
                  f"sigma {es.sigma:.3f}  {time.time() - t0:.0f} s", flush=True)
    print(f"stopped: {es.stop()}")

    # re-simulate the best parameters the normal way (dipole lengths set on the model) and check the logged error
    init_worker()
    best_dipole = {k: float(best[k]) for k in dp.PARAM_NAMES}
    model.dipole_length_params = best_dipole
    traces, penalties = [], []
    for onset, seed in SIMS:
        model.apply_params({**FIXED, **{n: best[n] for n in SEARCHED}, "input_onset": onset, "Ib_noise_seed": seed})
        model.initialize_state()
        model.simulate()
        tc = model.compute_timecourse(model.compute_dipoles())
        traces.append([tc[r] for r in FIT_ROIS])
        penalties.append(rate_penalties())
    tc_avg = np.mean(traces, axis=0)
    err_erp = tc_error(tc_avg)
    error = err_erp + np.sum(np.mean(penalties, axis=0))
    tc_sim, tc_target = dict(zip(FIT_ROIS, tc_avg)), dict(zip(FIT_ROIS, TARGET))
    if abs(error - best["err_tc"]) > 1e-6:
        raise RuntimeError(f"re-simulated error {error:.8f} differs from the logged best error {best['err_tc']:.8f}")

    model.plot_timecourse_comparison(tc_sim, tc_target, global_peak_rois=FIT_ROIS,
                                     save_path=os.path.join(session_dir, "best_tc_comparison.png"))
    model.save_timecourse_comparison(session_dir, tc_sim, tc_target, err_erp, filename="best_tc_comparison")
    dp.build_params_json(os.path.join(session_dir, "dipole_parameters_finetuned.json"), best_dipole)

    with open(os.path.join(session_dir, "optimization_summary.json"), "w") as f:
        json.dump({
            "fit_rois": list(FIT_ROIS),
            "best_params": {n: float(best[n]) for n in SEARCHED + SOLVED},
            "dipole_parameters": best_dipole,
            "best_combined_error": float(best["err_tc"]),
            "best_fit_err_tc": float(err_erp),
            "best_mid_rate_error": float(best["mid_rate_error"]),
            "best_saturation_error": float(best["saturation_error"]),
            "error_per_iteration": [float(e) for e in best_per_iteration],
            "n_evaluations": n_eval,
            "stop": {k: str(v) for k, v in es.stop().items()},
        }, f, indent=2)
    print(f"best error {best['err_tc']:.6f} (ERP {err_erp:.6f}); results in {session_dir}")


if __name__ == "__main__":
    main()
