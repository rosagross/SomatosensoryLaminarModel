"""
File: test_parameter_sweep.py
Description:
    Tests for the parameter-space exploration: the oscillation measures in
    model/oscillation_metrics.py and the parameter plumbing in run_parameter_sweep.py.

    Two failure modes these guard against, both of which happened while the sweep was
    being written and neither of which is visible in the output figures:

      1. A swept parameter that never reaches the model. The sweep runs, the figures
         render, and the curve for that parameter is simply flat - indistinguishable
         from a parameter the model genuinely does not care about. (apply_params did
         silently drop delay_factor_short, which is exactly this.)
      2. A peak measure that fires on its own reference. The raw spectrum of a
         *disconnected* network scores an "alpha peak" of 0.22 and a "gamma peak" of
         0.30 purely from the shape of its input, so any sweep scored that way would
         report oscillations everywhere.

Run with:
    WDDIR=<repo root> pytest Simulations/test_parameter_sweep.py
"""

import json
import os
import sys

import numpy as np
import pytest

WDDIR = os.getenv("WDDIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(WDDIR, "Simulations"))
sys.path.insert(0, os.path.join(WDDIR, "Simulations", "model"))

from oscillation_metrics import (  # noqa: E402
    ACTIVITY_FLOOR, BANDS, PEAK_PROMINENCE_MIN, amplitude_ratio, band_prominences,
    classify_regime, dominant_peak, network_gain_spectrum, oscillation_features,
    welch_spectrum)

FS = 1000.0
DT = 1.0 / FS


def _signal(freqs_amps, dur=8.0, noise=0.0, seed=0, decay=None, offset=0.0):
    """A test signal: sum of sinusoids, optional exponential decay, noise and offset."""
    rng = np.random.default_rng(seed)
    t = np.arange(0, dur, DT)
    x = np.full_like(t, offset)
    for f, a in freqs_amps:
        x = x + a * np.sin(2 * np.pi * f * t + rng.uniform(0, 2 * np.pi))
    if decay is not None:
        x = offset + (x - offset) * np.exp(-t / decay)
    if noise:
        x = x + noise * rng.standard_normal(t.size)
    return x


# ── the spectral measures ──────────────────────────────────────────────────────
def test_welch_grid_resolves_the_alpha_band():
    """A 1 s segment gives 1 Hz bins, so a peak can be placed inside 8-13 Hz."""
    freqs, psd, seg_stds = welch_spectrum(_signal([(10.0, 1.0)]), DT, seg_dur=1.0)
    assert np.isclose(freqs[1] - freqs[0], 1.0)
    assert freqs[int(np.argmax(psd))] == pytest.approx(10.0)
    assert seg_stds.size > 1


def test_a_pure_tone_is_found_in_its_own_band_only():
    freqs, psd, _ = welch_spectrum(_signal([(10.0, 1.0)], noise=0.1), DT)
    proms = band_prominences(freqs, psd)
    band, freq, prom = dominant_peak(proms)
    assert band == "alpha"
    assert freq == pytest.approx(10.0)
    assert prom > PEAK_PROMINENCE_MIN
    assert proms["theta"][0] < prom and proms["beta"][0] < prom


def test_gain_spectrum_of_the_control_against_itself_is_flat():
    """The core guarantee: the reference cannot score as a peak against itself.

    Scored raw, the same spectrum has a large apparent peak - which is why the raw
    spectrum is not what the sweep scores.
    """
    freqs, psd, _ = welch_spectrum(_signal([], noise=1.0, seed=3), DT)
    # give it the shaped, curved background a real control run has
    shaped = psd / (1.0 + (freqs / 10.0) ** 2) ** 2
    gain = network_gain_spectrum(shaped, shaped)
    for name, (prom, _) in band_prominences(freqs, gain).items():
        assert abs(prom) < 1e-9, f"{name} scored {prom} against its own reference"


def test_a_peak_survives_division_by_the_control():
    freqs, psd, _ = welch_spectrum(_signal([], noise=1.0, seed=4), DT)
    control = psd / (1.0 + (freqs / 10.0) ** 2) ** 2
    with_peak = control * (1.0 + 3.0 * np.exp(-0.5 * ((freqs - 10.0) / 1.5) ** 2))
    gain = network_gain_spectrum(with_peak, control)
    band, freq, prom = dominant_peak(band_prominences(freqs, gain))
    assert band == "alpha"
    assert freq == pytest.approx(10.0, abs=1.0)
    assert prom > PEAK_PROMINENCE_MIN


def test_dead_control_is_rejected_rather_than_divided_by():
    freqs, psd, _ = welch_spectrum(_signal([(10.0, 1.0)]), DT)
    assert network_gain_spectrum(psd, np.zeros_like(psd)) is None
    assert network_gain_spectrum(psd, None) is None
    assert network_gain_spectrum(psd, psd[:-3]) is None


# ── stationarity and regimes ───────────────────────────────────────────────────
def test_amplitude_ratio_separates_a_ring_down_from_stationary_noise():
    _, _, steady = welch_spectrum(_signal([], noise=1.0, seed=5), DT)
    _, _, ringing = welch_spectrum(_signal([(10.0, 1.0)], decay=1.0), DT)
    assert amplitude_ratio(steady) > 0.8, "stationary noise must not read as decaying"
    assert amplitude_ratio(ringing) < 0.2


@pytest.mark.parametrize("fluctuation, amp_ratio, prominence, rate_drive, expected", [
    (np.nan, 1.0, 0.5, 1e-2, "diverged"),
    (0.3, 1.0, 0.5, 1e-9, "pinned"),          # rate frozen: nothing else matters
    (ACTIVITY_FLOOR / 10, 1.0, 0.5, 1e-2, "fixed_point"),
    (0.3, 0.2, 0.5, 1e-2, "damped"),
    (0.3, 1.0, 0.5, 1e-2, "oscillation"),
    (0.3, 1.0, 0.0, 1e-2, "noise_driven"),
])
def test_regime_classification(fluctuation, amp_ratio, prominence, rate_drive, expected):
    assert classify_regime(fluctuation, amp_ratio, prominence, rate_drive) == expected


def test_features_of_a_noisy_rhythm_over_its_own_control():
    control_sig = _signal([], noise=1.0, seed=7, offset=5.0)
    rhythm = control_sig + 0.6 * np.sin(2 * np.pi * 10.0 * np.arange(0, 8.0, DT))
    _, psd_null, _ = welch_spectrum(control_sig, DT)
    feats = oscillation_features(rhythm, DT, rate_drive=1e-2, psd_null=psd_null)
    assert feats["scored_vs"] == "loopfree"
    assert feats["regime"] == "oscillation"
    assert feats["band"] == "alpha"
    assert feats["alpha_freq"] == pytest.approx(10.0)
    assert feats["mean_level"] == pytest.approx(5.0, abs=0.1)
    # every band named in BANDS gets its own columns
    for name, _, _ in BANDS:
        assert f"{name}_prom" in feats and f"{name}_freq" in feats


def test_pinned_rate_is_not_reported_as_an_oscillation():
    """A frozen firing rate outranks a clean spectral peak: the loop through that
    population is open, so the peak cannot be something the network is generating."""
    sig = _signal([(10.0, 1.0)], noise=0.05)
    feats = oscillation_features(sig, DT, rate_drive=1e-12)
    assert feats["regime"] == "pinned"
    assert feats["band"] == "none"


# ── parameter plumbing ─────────────────────────────────────────────────────────
@pytest.fixture(scope="module")
def model_and_reference():
    from run_parameter_sweep import reference_point
    from somato_model import SomatoModel, read_simulation_params
    base = read_simulation_params(WDDIR)
    base.update({"simulation_dur": 2.0, "input_onset": 3.0, "Iext_strength": 0.0,
                 "Iext_duration": 0.0, "Ib_noise_seed": 0})
    return SomatoModel(base, WDDIR=WDDIR), reference_point(base), base


def test_expand_theta_folds_the_thalamic_weights_back_into_a_list(model_and_reference):
    from run_parameter_sweep import THAL_CONNECT_NAMES, expand_theta
    _, ref, base = model_and_reference
    theta = dict(ref)
    theta["thal_ItoE"] = 7.0
    params = expand_theta(theta, base)
    assert "thal_ItoE" not in params, "the scalar name must not reach the model"
    assert params["thal_connect"][THAL_CONNECT_NAMES.index("thal_ItoE")] == 7.0
    # the entries that were not swept keep their base values
    assert params["thal_connect"][1] == base["thal_connect"][1]


def test_center_moves_the_reference_point(model_and_reference, tmp_path):
    """--center / --center-json must move only the parameters they name."""
    from run_parameter_sweep import parse_center, reference_point
    _, ref, base = model_and_reference

    center = parse_center(["coupling_strength=7.5", "e2_tau=3"])
    assert center == {"coupling_strength": 7.5, "e2_tau": 3.0}

    path = tmp_path / "center.json"
    path.write_text(json.dumps({"strength_I": 0.81, "coupling_strength": 1.0}))
    # the explicit flag wins over the file for a parameter named in both
    center = parse_center(["coupling_strength=7.5"], str(path))
    assert center == {"strength_I": 0.81, "coupling_strength": 7.5}

    moved = reference_point(base, center)
    assert moved["coupling_strength"] == 7.5
    assert moved["strength_I"] == 0.81
    untouched = [k for k in ref if k not in center]
    assert all(moved[k] == ref[k] for k in untouched)

    with pytest.raises(ValueError):
        reference_point(base, {"not_a_parameter": 1.0})
    with pytest.raises(ValueError):
        parse_center(["coupling_strength"])


def test_every_swept_parameter_reaches_the_model(model_and_reference):
    """Each parameter in PARAM_RANGES must change the connectivity, the time constants,
    or the background input - otherwise its sweep axis is silently inert."""
    from run_parameter_sweep import PARAM_RANGES, expand_theta
    model, ref, base = model_and_reference

    model.apply_params(expand_theta(ref, base))
    W0, tau0, Ib0 = model.W.copy(), model.tau.copy(), model.Ib.copy()

    inert = []
    for name, (lo, hi) in PARAM_RANGES.items():
        theta = dict(ref)
        # move to whichever end of the range is further from the reference value
        theta[name] = hi if abs(hi - ref[name]) > abs(lo - ref[name]) else lo
        model.apply_params(expand_theta(theta, base))
        changed = (not np.array_equal(model.W, W0)
                   or not np.array_equal(model.tau, tau0)
                   or not np.array_equal(model.Ib, Ib0))
        if not changed:
            inert.append(name)
        model.apply_params(expand_theta(ref, base))
    assert not inert, f"parameters that never reach the model: {inert}"


def test_the_loop_free_control_keeps_every_signal_alive(model_and_reference):
    """The control run must leave a usable reference for every signal.

    Zeroing the gains outright does not: the ROI dipoles project
    self.potential[:, :-2] and so exclude the background column, leaving them
    identically zero, and the thalamus receives no background input at all. Scaling the
    gains down instead keeps the one-relay feedforward path, and with it a reference.
    """
    from run_parameter_sweep import NULL_GAIN_PARAMS, NULL_GAIN_SCALE, expand_theta
    model, ref, base = model_and_reference

    null_theta = dict(ref)
    for name in NULL_GAIN_PARAMS:
        null_theta[name] = float(ref[name]) * NULL_GAIN_SCALE
    model.apply_params(expand_theta(null_theta, base))
    model.initialize_state()
    model.simulate()

    potentials = np.sum(model.potential, axis=1)
    labels = list(model.get_population_labels())
    for pop in ("E1", "E1S2", "E3b", "ThalE"):
        sig = potentials[labels.index(pop)]
        assert np.all(np.isfinite(sig))
        assert sig.std() > 0, f"{pop} has no fluctuation left to reference against"
