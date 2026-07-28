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
    peaks). FOOOF works in log10 space internally; ``_spectrum_flat`` is
    ``power_spectrum - _ap_fit`` in log10, so exponentiating brings it back to a
    linear, strictly-positive spectrum whose downstream unit-sum normalisation
    stays valid.

    Args:
        freqs:  1-D array of frequencies (Hz); values <= 0 are excluded.
        power:  1-D array of (linear) power at `freqs`.
        fmin, fmax: fit band in Hz.
        peak_width_limits: FOOOF peak width limits (Hz).

    Returns:
        (fit_freqs, flattened) where `fit_freqs` is the cropped `[fmin, fmax]`
        frequency grid FOOOF fit over and `flattened` is the aperiodic-removed
        (linear, positive) power on that grid.
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
        flattened = np.power(10.0, fm._spectrum_flat)   # aperiodic-removed, back to linear
        return fm.freqs, flattened
    else:
        return freqs, power    


    
