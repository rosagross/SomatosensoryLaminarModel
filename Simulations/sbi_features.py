"""
File: sbi_features.py
Author: Rosa Grossmann
Description:
    Summary-statistic extraction for simulation-based inference (SBI) of the
    SomatoModel. Reduces the per-ROI dipole observables (time course + Morlet
    time-frequency power, as produced by SomatoModel.compute_timecourse /
    compute_timefreq and the matching load_target_* loaders) to a fixed-length
    feature vector `x`.

    The exact same function is applied to simulated and to measured data, so the
    observation `x_obs` and every simulated `x` are directly comparable.

    Feature layout per ROI (A3b, A1, S2), 10 features each -> 30 total:
        time course : top-3 peak amplitudes  + their 3 latencies (ms)   -> 6
        time-freq   : top-2 dominant (most baseline-normalised) freqs   -> 2
        power       : mean baseline power (-200..0 ms) + late (0..400)   -> 2
"""

import numpy as np
from scipy.signal import find_peaks

# ── shared axes (compute_timecourse / compute_timefreq use this exact grid) ──────
ROIS      = ("A3b", "A1", "S2")
TIME_MS   = np.linspace(-500, 400, 451)          # 2 ms steps, stimulus at t=0 (idx 250)
TF_FREQS  = np.arange(1, 41, 1).astype(float)    # 1..40 Hz, 40 bins

N_TC_PEAKS = 3
N_TF_FREQS = 2

_BASELINE_MASK = (TIME_MS >= -200) & (TIME_MS < 0)     # -200..0 ms
_LATE_MASK     = (TIME_MS >= 0) & (TIME_MS <= 400)     # 0..+400 ms (post-stimulus)
_PEAK_MASK     = (TIME_MS >= 0) & (TIME_MS <= 400)     # detect evoked peaks post-stimulus
_EPS           = 1e-10

N_FEATURES_PER_ROI = 2 * N_TC_PEAKS + N_TF_FREQS + 2
N_FEATURES         = len(ROIS) * N_FEATURES_PER_ROI    # 30


def _tc_peak_features(tc, n_peaks=N_TC_PEAKS):
    """Top-`n_peaks` |peaks| of a baseline-corrected trace and their latencies.

    Returns (amplitudes, latencies_ms), each length `n_peaks`, ranked by |amplitude|
    descending. Padded with (0.0 amplitude, 0.0 latency) when fewer peaks are found.
    """
    win_idx = np.where(_PEAK_MASK)[0]
    seg     = tc[win_idx]
    peak_rel, _ = find_peaks(np.abs(seg))
    if peak_rel.size == 0:
        # fall back to the single largest deflection in the window
        peak_rel = np.array([int(np.argmax(np.abs(seg)))])

    order   = np.argsort(np.abs(seg[peak_rel]))[::-1][:n_peaks]
    sel     = peak_rel[order]
    amps    = seg[sel]
    lats    = TIME_MS[win_idx[sel]]

    pad = n_peaks - sel.size
    if pad > 0:
        amps = np.concatenate([amps, np.zeros(pad)])
        lats = np.concatenate([lats, np.zeros(pad)])
    return amps, lats


def _tf_dominant_freqs(tf, n_freqs=N_TF_FREQS):
    """The `n_freqs` frequencies with the largest post-stimulus / baseline power ratio."""
    bl    = tf[:, _BASELINE_MASK].mean(axis=1)
    post  = tf[:, _LATE_MASK].mean(axis=1)
    ratio = post / (bl + _EPS)
    top   = np.argsort(ratio)[::-1][:n_freqs]
    return TF_FREQS[top]


def _power_features(tf):
    """Mean baseline (-200..0 ms) and late (0..+400 ms) power, averaged over all freqs."""
    return tf[:, _BASELINE_MASK].mean(), tf[:, _LATE_MASK].mean()


def extract_summary_features(tf_dict, tc_dict, step_size=0.001):
    """Reduce per-ROI TF + TC observables to the fixed-length SBI feature vector.

    Args:
        tf_dict: roi label -> (n_freqs=40, n_times>=451) Morlet power array.
        tc_dict: roi label -> (n_times>=451,) baseline-corrected trace.
        step_size: kept for interface symmetry (axes are fixed at 2 ms).

    Returns:
        np.ndarray of shape (N_FEATURES,) == (30,). Order matches feature_names().
    """
    n = TIME_MS.size

    # Time course: one shared peak across areas (this whole dict is either all-sim or
    # all-obs). Dividing every area by the same peak preserves the amplitude ratios
    # between areas while keeping the simulated and observed vectors comparable.
    tc_peak = max(np.max(np.abs(np.asarray(tc_dict[roi], float)[:n][_PEAK_MASK])) for roi in ROIS) + _EPS

    feats = []
    for roi in ROIS:
        tc = np.asarray(tc_dict[roi], dtype=float)[:n]
        tf = np.asarray(tf_dict[roi], dtype=float)[:, :n]

        amps, lats = _tc_peak_features(tc)
        feats.extend((amps / tc_peak).tolist())     # joint peak norm (cross-area ratios preserved)
        feats.extend(lats.tolist())                 # latencies are scale-free
        feats.extend(_tf_dominant_freqs(tf).tolist())  # already a baseline ratio
        # time-freq power: each area normalized by ITS OWN baseline (baseline feature
        # -> ~1, late feature -> fold-change), mirroring compute_error_timefreq.
        bl_pow, late_pow = _power_features(tf)
        own_bl = bl_pow + _EPS
        feats.extend([bl_pow / own_bl, late_pow / own_bl])

    return np.asarray(feats, dtype=float)


def feature_names():
    """Human-readable name for each of the N_FEATURES features (same order)."""
    names = []
    for roi in ROIS:
        names += [f"{roi}_tc_amp{i+1}"  for i in range(N_TC_PEAKS)]
        names += [f"{roi}_tc_lat{i+1}"  for i in range(N_TC_PEAKS)]
        names += [f"{roi}_tf_domfreq{i+1}" for i in range(N_TF_FREQS)]
        names += [f"{roi}_pow_baseline", f"{roi}_pow_late"]
    return names
