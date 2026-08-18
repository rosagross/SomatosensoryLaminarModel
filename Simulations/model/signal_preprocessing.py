"""
File: signal_preprocessing.py
Description:
    Small, reusable signal-preprocessing helpers for the optimization pipeline.

    Used to preprocess the measured fitting targets (and, symmetrically, the
    simulated signals) so the GA fits the features that matter rather than
    measurement noise / the aperiodic 1/f slope:
      - smooth_timecourse : Savitzky-Golay smoothing of an ERP time course.
      - remove_aperiodic  : FOOOF-based removal of the 1/f background from a
                            power spectrum, returning the flattened (peaks-only)
                            spectrum.
      - log_aperiodic_residual / spectral_prominence : local, FOOOF-free
                            measures of how far a band sticks out above the
                            surrounding 1/f background (used by the pre-stim
                            peak objective).
"""

import numpy as np
from scipy.signal import savgol_filter


def smooth_timecourse(y, dt_ms=2.0, window_ms=30.0, polyorder=3):
    """Savitzky-Golay smoothing of a 1-D time course.

    Preserves ERP peak shape better than a moving average. The window length is
    derived from `window_ms / dt_ms`, forced odd and strictly greater than
    `polyorder` (as required by savgol_filter).

    Args:
        y:         1-D array, the time course to smooth.
        dt_ms:     sample spacing in ms (2 ms for the 500 Hz analysis grid).
        window_ms: smoothing window length in ms.
        polyorder: polynomial order of the Savitzky-Golay fit.

    Returns:
        np.ndarray of the same shape as `y`.
    """
    y = np.asarray(y, dtype=float)
    window = int(round(window_ms / dt_ms))
    if window % 2 == 0:
        window += 1                      # savgol requires an odd window length
    window = max(window, polyorder + 1)
    if window % 2 == 0:
        window += 1
    window = min(window, len(y) if len(y) % 2 == 1 else len(y) - 1)
    if window <= polyorder:
        return y                         # too short to smooth meaningfully
    return savgol_filter(y, window_length=window, polyorder=polyorder)


def remove_aperiodic(freqs, power, fmin=1.0, fmax=40.0, peak_width_limits=(1.0, 12.0)):
    """Remove the aperiodic (1/f) component of a power spectrum with FOOOF.

    Fits a FOOOF model over ``[fmin, fmax]`` and returns the flattened spectrum
    (the measured power with the aperiodic fit subtracted, i.e. the oscillatory
    peaks) **in log10 units** -- FOOOF's ``_spectrum_flat`` is
    ``power_spectrum - _ap_fit`` in log10, and that is what is returned.

    This used to return ``10 ** _spectrum_flat`` to get back to a linear,
    strictly-positive spectrum. That is what made the pre-stim objective useless:
    the log residual is ~0 wherever the spectrum follows its own aperiodic fit, so
    exponentiating maps *any* spectrum -- a real rhythm, a slow drift, or numerical
    residue -- onto a curve with baseline 1.0. The subsequent unit-sum normalisation
    then pushed every candidate onto ~1/n_bins and the GA's error range collapsed to
    a 0.1%-wide plateau. In log10 units a missing peak still costs what it should.
    Callers that need a positive spectrum must exponentiate themselves.

    Args:
        freqs:  1-D array of frequencies (Hz); values <= 0 are excluded.
        power:  1-D array of (linear) power at `freqs`.
        fmin, fmax: fit band in Hz.
        peak_width_limits: FOOOF peak width limits (Hz).

    Returns:
        (fit_freqs, flattened) where `fit_freqs` is the cropped `[fmin, fmax]`
        frequency grid FOOOF fit over and `flattened` is the aperiodic-removed
        power on that grid, in log10 units (0 = on the aperiodic fit).
    """
    from fooof import FOOOF

    freqs = np.asarray(freqs, dtype=float)
    power = np.asarray(power, dtype=float)

    # FOOOF requires strictly positive frequencies and power.
    valid = (freqs > 0) & np.isfinite(power) & (power > 0)
    freqs_valid, power_valid = freqs[valid], power[valid]

    if len(freqs_valid) > 0:
        fm = FOOOF(peak_width_limits=peak_width_limits, verbose=False)
        fm.fit(freqs_valid, power_valid, (fmin, fmax))
        return fm.freqs, np.asarray(fm._spectrum_flat, dtype=float)
    else:
        return freqs, np.full_like(power, np.nan)


def log_aperiodic_residual(freqs, power, ref_bands=None):
    """Log10 power minus a straight-line (log-log) aperiodic reference.

    The reference is a least-squares line through log10(power) vs log10(freq),
    fitted only on the bins inside `ref_bands`. Restricting the fit to bands that
    flank the feature of interest keeps the reference *local*: any smooth 1/f-like
    background is then removed exactly where it matters, without a global spectral
    model having to describe the whole 1-40 Hz range.

    That locality is the point. A global power law (what FOOOF fits) is a poor
    description of these measured spectra over 1-40 Hz -- they are visibly curved
    and would need a knee -- so its residual is a broad positive bow of ~0.4 log10
    units across 5-32 Hz, which swamps the alpha peak and merges it with the beta
    shoulder. A line through the immediate flanks does not have that problem.

    The residual is scale-free: multiplying `power` by a constant only shifts the
    fitted intercept, so simulated and measured spectra are directly comparable
    without any normalisation.

    Args:
        freqs:     1-D array of frequencies (Hz).
        power:     1-D array of (linear) power at `freqs`.
        ref_bands: sequence of (lo, hi) Hz windows the reference line is fitted on.
                   None (default) fits on every valid bin, giving the plain
                   broadband 1/f residual.

    Returns:
        1-D array of the same length as `freqs`; NaN at bins where `power` is not
        strictly positive, and all-NaN if the reference could not be fitted (fewer
        than two usable reference bins).
    """
    freqs = np.asarray(freqs, dtype=float)
    power = np.asarray(power, dtype=float)

    resid = np.full(freqs.shape, np.nan)
    valid = (freqs > 0) & np.isfinite(power) & (power > 0)
    if not valid.any():
        return resid

    if ref_bands is None:
        ref = valid.copy()
    else:
        ref = np.zeros(freqs.shape, dtype=bool)
        for lo, hi in ref_bands:
            ref |= (freqs >= lo) & (freqs <= hi)
        ref &= valid

    if ref.sum() < 2:
        return resid

    log_f = np.log10(freqs[valid])
    log_p = np.log10(power[valid])
    slope, intercept = np.polyfit(np.log10(freqs[ref]), np.log10(power[ref]), 1)
    resid[valid] = log_p - (intercept + slope * log_f)
    return resid


def spectral_prominence(freqs, power, band, ref_bands):
    """How far the strongest bin in `band` rises above the local 1/f reference.

    Measures the feature the pre-stim fit is actually after: "is there a bump here,
    and how big is it". The unit is log10 power ratio, so 0.2 means the peak sits
    ~1.6x above the background the flanking bands imply. Both the simulated and the
    measured spectrum are reduced to this number (plus the peak frequency), and the
    objective compares those instead of comparing 16-bin spectral shapes.

    Args:
        freqs:     1-D array of frequencies (Hz).
        power:     1-D array of (linear) power at `freqs`.
        band:      (lo, hi) Hz window the peak is searched in.
        ref_bands: sequence of (lo, hi) Hz windows flanking `band`, defining the
                   local aperiodic reference (see log_aperiodic_residual).

    Returns:
        (prominence, peak_freq, residual). `prominence` and `peak_freq` are NaN if
        the band holds no bin with a finite residual; `residual` is the full-length
        log10 residual array (useful for plotting what was measured).
    """
    freqs = np.asarray(freqs, dtype=float)
    resid = log_aperiodic_residual(freqs, power, ref_bands)

    in_band = (freqs >= band[0]) & (freqs <= band[1]) & np.isfinite(resid)
    if not in_band.any():
        return np.nan, np.nan, resid

    idx = np.flatnonzero(in_band)[np.argmax(resid[in_band])]
    return float(resid[idx]), float(freqs[idx]), resid    


    
