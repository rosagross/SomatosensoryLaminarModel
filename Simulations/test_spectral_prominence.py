"""
File: test_spectral_prominence.py
Author: Rosa Grossmann
Description: Tests for the pre-stimulus peak features the "ps" optimization is scored on
    (signal_preprocessing.spectral_prominence / log_aperiodic_residual).

    These guard the properties the objective depends on:
      - a spectrum with no oscillation scores 0, whatever its 1/f slope, so the GA cannot
        be rewarded for the slope;
      - a bump in a band scores positive and is located at the right bin;
      - a slow drift piling into the lowest bin -- the failure mode the old objective kept
        selecting for -- does *not* read as an alpha peak;
      - the measured group spectra come out with alpha at 10 Hz in every ROI.

    Run with:  SIMDIR=... RESDIR=... WDDIR=... pytest Simulations/test_spectral_prominence.py
"""

import os
import sys

import numpy as np
import pandas as pd
import pytest

WDDIR = os.getenv("WDDIR")
RESDIR = os.getenv("RESDIR")
sys.path.append(os.path.join(WDDIR, "Simulations", "model"))
from signal_preprocessing import log_aperiodic_residual, spectral_prominence

# same bands/flanks as somato_model.PRESTIM_* (imported from there would pull in the whole
# model, which needs a full environment; kept in sync deliberately)
ALPHA_BAND, ALPHA_FLANKS = (8.0, 13.0), ((5.0, 7.5), (14.0, 17.5))
BETA_BAND, BETA_FLANKS = (17.5, 30.0), ((12.5, 15.0), (32.5, 40.0))

FREQS = np.arange(2.5, 40.1, 2.5)   # the 2.5 Hz measured/simulated grid


def _spectrum(exponent=1.0, peaks=()):
    """A 1/f**exponent spectrum with optional multiplicative Gaussian bumps."""
    power = FREQS ** -exponent
    for centre, amplitude, width in peaks:
        power = power * (1 + amplitude * np.exp(-0.5 * ((FREQS - centre) / width) ** 2))
    return power


@pytest.mark.parametrize("exponent", [0.0, 0.5, 1.0, 2.0, 4.0])
def test_pure_power_law_has_no_prominence(exponent):
    """Any smooth 1/f background must score 0 -- the objective scores peaks, not slope."""
    power = _spectrum(exponent)
    for band, flanks in ((ALPHA_BAND, ALPHA_FLANKS), (BETA_BAND, BETA_FLANKS)):
        prominence, _, _ = spectral_prominence(FREQS, power, band, flanks)
        assert prominence == pytest.approx(0.0, abs=1e-9)


def test_alpha_peak_is_found_and_beta_is_not():
    prominence, freq, _ = spectral_prominence(
        FREQS, _spectrum(peaks=[(10, 2.0, 2.0)]), ALPHA_BAND, ALPHA_FLANKS)
    assert freq == pytest.approx(10.0)
    assert prominence > 0.3

    beta, _, _ = spectral_prominence(
        FREQS, _spectrum(peaks=[(10, 2.0, 2.0)]), BETA_BAND, BETA_FLANKS)
    assert beta < 0.05


def test_alpha_peak_with_smaller_beta_peak():
    """The shape the fit is aiming for: alpha dominant, beta present but weaker."""
    power = _spectrum(peaks=[(10, 2.0, 2.0), (20, 0.8, 3.0)])
    alpha, alpha_freq, _ = spectral_prominence(FREQS, power, ALPHA_BAND, ALPHA_FLANKS)
    beta, beta_freq, _ = spectral_prominence(FREQS, power, BETA_BAND, BETA_FLANKS)
    assert alpha_freq == pytest.approx(10.0)
    assert beta_freq == pytest.approx(20.0)
    assert alpha > beta > 0.05


def test_slow_drift_is_not_read_as_an_alpha_peak():
    """A spike in the lowest bin must not look like a rhythm.

    With the background noise off, a parameter set that has merely not settled puts nearly
    all its pre-stimulus power in the 2.5 Hz bin. That is what the old whole-spectrum
    objective kept selecting. Letting the 2.5 Hz bin anchor the reference line tilts it
    enough to report a ~0.5 alpha prominence for exactly this signal, which is why the alpha
    reference starts at 5 Hz.
    """
    drift = np.full_like(FREQS, 1e-30)
    drift[0] = 1.0
    prominence, _, _ = spectral_prominence(FREQS, drift, ALPHA_BAND, ALPHA_FLANKS)
    assert prominence == pytest.approx(0.0, abs=1e-9)


def test_prominence_is_scale_free():
    """Simulated and measured spectra are compared without normalisation, so scale must not matter."""
    power = _spectrum(peaks=[(10, 2.0, 2.0)])
    a1, f1, _ = spectral_prominence(FREQS, power, ALPHA_BAND, ALPHA_FLANKS)
    a2, f2, _ = spectral_prominence(FREQS, power * 1e7, ALPHA_BAND, ALPHA_FLANKS)
    assert a1 == pytest.approx(a2)
    assert f1 == pytest.approx(f2)


def test_dead_spectrum_returns_nan():
    """An all-zero ROI must be reported as unusable, not silently scored."""
    prominence, freq, resid = spectral_prominence(
        FREQS, np.zeros_like(FREQS), ALPHA_BAND, ALPHA_FLANKS)
    assert np.isnan(prominence) and np.isnan(freq)
    assert np.isnan(resid).all()


def test_broadband_residual_is_zero_for_a_power_law():
    resid = log_aperiodic_residual(FREQS, _spectrum(1.3))
    assert np.allclose(resid, 0.0, atol=1e-9)


@pytest.mark.skipif(RESDIR is None, reason="needs RESDIR for the measured target CSV")
@pytest.mark.parametrize("roi", ["BA3b", "BA1", "S2"])
def test_measured_targets_have_an_alpha_peak_at_10hz(roi):
    path = os.path.join(
        RESDIR, "Figures", "Main", "eeg_results", "source_reconstruction", "group",
        "_preprestim_corrected", "roi_epochswise",
        "group_roi_prestim_spectrum_ses-elec_preprestim_corrected.csv")
    df = pd.read_csv(path)
    df = df[(df["modality"] == "elec") & (df["roi"] == roi)].sort_values("freq_hz")
    freqs = df["freq_hz"].to_numpy()
    power = df["power"].to_numpy()
    band = (freqs >= 1) & (freqs <= 40)

    alpha, alpha_freq, _ = spectral_prominence(
        freqs[band], power[band], ALPHA_BAND, ALPHA_FLANKS)
    assert alpha_freq == pytest.approx(10.0)
    assert alpha > 0.05
