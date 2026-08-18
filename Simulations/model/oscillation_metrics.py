"""
File: oscillation_metrics.py
Description:
    Signal-level measures used to characterise *ongoing* (stimulus-free) oscillations
    of the model, for the parameter-space exploration in
    Simulations/run_parameter_sweep.py and Analysis/step003_plot_parameter_space.py.

    The question these answer, per simulated signal, is:
      1. is the population still taking part in the dynamics, or is its firing rate
         pinned against a rail of the sigmoid?
      2. does anything fluctuate at all, or has the network settled to a fixed point?
      3. is the fluctuation sustained, or is it the ring-down of the initial transient?
      4. is there a spectral peak the *network* is producing, and at which frequency?
      5. specifically, how prominent is the alpha (8-13 Hz) peak?

    (4) and (5) reuse `signal_preprocessing.spectral_prominence`, the same scale-free
    measure the pre-stimulus objective of the GA uses (see run_optimization.py), so a
    "peak" here means the same thing it means there - but applied to the *network gain*
    spectrum rather than the raw one. With the stimulus off, everything a population
    sees is background noise through a chain of synaptic kernels, and that input
    shaping alone scores a large apparent peak; dividing by the spectrum of the same
    parameter set with its loops opened removes it. See `network_gain_spectrum`.
    (2) and (3) mirror SomatoModel._prestim_fluctuation / _prestim_amplitude_ratio, but
    act on an arbitrary 1-D signal instead of the ROI dipoles, so population potentials
    can be scored with the same yardstick.

    This module deliberately does not import somato_model (which pulls in mne): the
    band constants below are the same values as its PRESTIM_* constants, duplicated
    the way Simulations/test_spectral_prominence.py duplicates them.
"""

import numpy as np

from signal_preprocessing import spectral_prominence


# ── bands ──────────────────────────────────────────────────────────────────────
# Same values as somato_model.PRESTIM_ALPHA_BAND / _FLANKS / BETA_*.
ALPHA_BAND, ALPHA_FLANKS = (8.0, 13.0), ((5.0, 7.5), (14.0, 17.5))
BETA_BAND, BETA_FLANKS = (17.5, 30.0), ((12.5, 15.0), (32.5, 40.0))

# Every band is scored against its *own* flanking bands, never against one global
# reference line. A global log-log fit is only unbiased for a true power law, and
# these spectra are not: the background is Ornstein-Uhlenbeck noise (a Lorentzian,
# corner at 1/(2*pi*Ib_noise_tau) ~ 10 Hz) shaped further by the second-order synaptic
# kernel, so log10 power vs log10 frequency is visibly curved. Fitting one line
# through it leaves a broad positive bow that reads as a 13 Hz "peak" in every single
# run, whatever the parameters - the exact failure log_aperiodic_residual's docstring
# warns about. Local flanks remove the curvature where it matters.
BANDS = (
    ("theta", (4.0, 7.0), ((2.0, 3.0), (8.0, 12.0))),
    ("alpha", ALPHA_BAND, ALPHA_FLANKS),
    ("beta", BETA_BAND, BETA_FLANKS),
    ("gamma", (31.0, 50.0), ((22.0, 28.0), (52.0, 60.0))),
)

# The spectrum is computed a little wider than the highest band so gamma has an upper
# flank to be referenced against.
FMIN, FMAX = 1.0, 60.0

# ── regime thresholds ──────────────────────────────────────────────────────────
# Relative fluctuation (std / max|x| of the settled signal) below which the signal is
# a settled fixed point and its spectrum is numerical residue. Same floor as
# run_optimization.PS_ACTIVITY_FLOOR.
ACTIVITY_FLOOR = 1e-6
# Amplitude ratio (late/early segment amplitude, see `amplitude_ratio`) a signal must
# hold to count as sustained rather than as a decaying transient. Same value as
# run_optimization.PS_STATIONARY_MIN_RATIO, though that one measures min/max.
STATIONARY_MIN_RATIO = 0.8
# Log10 power ratio the strongest bin must rise above the reference to count as an
# oscillatory peak rather than broadband, noise-driven fluctuation. Applied to the
# network gain spectrum (see `network_gain_spectrum`), whose floor is 0 by
# construction - the control run divided by itself scores exactly 0.000 in every band -
# so a low threshold is safe here. It is NOT safe on the raw spectrum, where a network
# with its loops opened already scores 0.2-0.3 from its own input shaping; rows with
# scored_vs == 'flat' (no usable control, e.g. every run of a noise-free sweep) carry
# that offset and should be read with it in mind.
PEAK_PROMINENCE_MIN = 0.05
# How much a population's firing rate must still move, as a fraction of its own
# sigmoid ceiling (std(rate) / m_max), for it to be participating in the dynamics at
# all. Below this it is pinned against one of the two rails of the sigmoid - jammed at
# m_max, or silenced at ~0 - and its output no longer responds to its input, so the
# loop through it is open however large the coupling gains are. Both rails matter: at
# the simulation_parameter.json defaults, coupling_strength > 0 drives some populations
# to the ceiling and silences others, and a test on the mean rate alone would only see
# the first kind.
RATE_DRIVE_MIN = 1e-6

REGIMES = ("diverged", "pinned", "fixed_point", "damped", "noise_driven",
           "oscillation")


def _settled(signal, step_size, settle_s):
    """The part of `signal` after the initialisation transient."""
    return np.asarray(signal, dtype=float)[int(round(settle_s / step_size)):]


def welch_spectrum(signal, step_size, seg_dur=1.0, overlap=0.5, settle_s=2.0,
                   fmin=FMIN, fmax=FMAX):
    """Welch-averaged power spectrum of the settled part of a 1-D signal.

    Mirrors helper_functions.compute_freq_spectrum per segment (mean-detrend + Hann +
    |rfft|^2 / n^2) and averages the periodograms, exactly like
    SomatoModel.compute_prestim_spectrum, but over a stimulus-free run rather than the
    pre-stimulus period. Averaging matters: a single periodogram of one noise
    realisation has ~100% variance per bin, so a random low bin can look like a peak.

    `seg_dur` sets the frequency resolution (1 s -> 1 Hz, enough to place a peak
    inside the 8-13 Hz alpha band; the 400 ms of the measured epochs would give only
    two bins there).

    Args:
        signal:    1-D array, the simulated signal (potential, rate or dipole).
        step_size: sample spacing in s.
        seg_dur:   segment length in s (1/seg_dur = frequency resolution).
        overlap:   fractional overlap between segments.
        settle_s:  seconds skipped at the start of the run (initialisation transient).
        fmin, fmax: returned frequency range in Hz.

    Returns:
        (freqs, psd, seg_stds) - `seg_stds` is the std of each (detrended) segment,
        used by `amplitude_ratio` to tell a sustained rhythm from a ring-down.
        Returns empty arrays if the settled part is shorter than one segment.
    """
    x = _settled(signal, step_size, settle_s)
    seg_len = int(round(seg_dur / step_size))
    if len(x) < seg_len or seg_len < 4:
        return np.array([]), np.array([]), np.array([])

    hop = max(int(round(seg_len * (1.0 - overlap))), 1)
    starts = list(range(0, len(x) - seg_len + 1, hop))
    win = np.hanning(seg_len)
    freqs = np.fft.rfftfreq(seg_len, d=step_size)
    fmask = (freqs >= fmin) & (freqs <= fmax)

    psd = np.zeros(int(fmask.sum()))
    seg_stds = np.empty(len(starts))
    for i, s0 in enumerate(starts):
        seg = x[s0:s0 + seg_len]
        seg = seg - seg.mean()
        seg_stds[i] = seg.std()
        psd += ((np.abs(np.fft.rfft(seg * win)) ** 2) / seg_len ** 2)[fmask]

    return freqs[fmask], psd / len(starts), seg_stds


def relative_fluctuation(signal, step_size, settle_s=2.0):
    """std / max|x| of the settled signal, in [0, ~1].

    Referenced to the signal's own level rather than to an evoked response, so it is
    comparable across populations with very different DC offsets and reads ~1e-9 for a
    settled fixed point regardless of how large that fixed point is. Same definition as
    SomatoModel._prestim_fluctuation.
    """
    x = _settled(signal, step_size, settle_s)
    if len(x) == 0 or not np.all(np.isfinite(x)):
        return np.nan
    level = max(float(np.abs(x).max()), np.finfo(float).tiny)
    return float(x.std()) / level


def amplitude_ratio(seg_stds):
    """Late / early segment amplitude: ~1 sustained, <<1 ringing down, >1 growing.

    SomatoModel._prestim_amplitude_ratio uses min/max over the segments, which works
    for the noise-free runs it was written for but is far too strict once the
    background noise is on: the per-segment std of a stationary noise-driven signal
    scatters by tens of percent, so min/max drops below any sensible threshold and
    every noisy run would be called a decaying transient. Comparing the *mean* of the
    last third of the segments to the mean of the first third averages that scatter
    out while still falling well below 1 for a genuine ring-down, which decays
    monotonically.
    """
    seg_stds = np.asarray(seg_stds, dtype=float)
    if seg_stds.size == 0 or not np.all(np.isfinite(seg_stds)):
        return np.nan
    k = max(seg_stds.size // 3, 1)
    early = float(seg_stds[:k].mean())
    if early <= 0:
        return 0.0
    return float(seg_stds[-k:].mean()) / early


def network_gain_spectrum(psd, psd_null):
    """psd / psd_null - what the recurrent network adds on top of its own input.

    Why this and not the raw spectrum: with the stimulus off, everything that reaches a
    population is background noise poured through the synaptic kernels, and that route
    alone produces a strongly shaped spectrum - an Ornstein-Uhlenbeck corner at
    1/(2*pi*Ib_noise_tau) (~10 Hz at the default 16 ms, i.e. inside the alpha band)
    followed by the second-order synaptic kernel's own roll-off (~26 Hz at tau = 6 ms).
    Scoring the raw spectrum against any smooth reference - a local line, or FOOOF with
    a knee - reads that shaping, plus the per-bin scatter of the finite Welch average,
    as peaks: the *disconnected* network scores a "gamma peak" of 0.30 and an "alpha
    peak" of 0.22, with no recurrent connectivity at all to produce them.

    Dividing by the spectrum of the same parameter set with every coupling gain set to
    zero - same taus, same delays, same noise realisation - cancels all of it: the
    input shaping, the kernel roll-off, and most of the realisation-specific bin
    scatter. What is left is the network's own transfer function, and the control case
    (null divided by itself) scores exactly 0 in every band.

    Args:
        psd:      power spectrum of the simulation.
        psd_null: power spectrum of the disconnected run, same frequency grid.

    Returns:
        The ratio, or None if `psd_null` is unusable (all-zero or wrong length), which
        is the case for populations the disconnected network leaves with no input at
        all - the thalamus, which in a stimulus-free run is driven only by the cortex.
        Callers should fall back to scoring the raw spectrum there.
    """
    if psd_null is None or len(psd_null) != len(psd):
        return None
    if not np.all(np.isfinite(psd_null)) or np.nanmax(psd_null) <= 0:
        return None
    # only bins where the null carries real power can be corrected
    usable = psd_null > np.nanmax(psd_null) * 1e-12
    if usable.sum() < 4:
        return None
    ratio = np.full(psd.shape, np.nan)
    ratio[usable] = psd[usable] / psd_null[usable]
    return ratio


def band_prominences(freqs, psd):
    """Prominence and peak frequency of every band in BANDS.

    Returns:
        dict band name -> (prominence, peak_freq), each in (log10 power ratio, Hz) and
        nan where the band or its flanks hold no usable bin.
    """
    out = {}
    for name, band, flanks in BANDS:
        if len(freqs) == 0:
            out[name] = (np.nan, np.nan)
            continue
        prom, freq, _ = spectral_prominence(freqs, psd, band, flanks)
        out[name] = (prom, freq)
    return out


def dominant_peak(band_proms):
    """The strongest of the per-band peaks: (band name, peak_freq, prominence).

    Taking the maximum over locally-referenced band prominences rather than the argmax
    of one global residual keeps the comparison fair across bands (each is measured
    against its own flanks) and keeps the answer meaningful on a curved spectrum.
    """
    best = ("none", np.nan, -np.inf)
    for name, (prom, freq) in band_proms.items():
        if np.isfinite(prom) and prom > best[2]:
            best = (name, freq, prom)
    if not np.isfinite(best[2]):
        return "none", np.nan, np.nan
    return best


def classify_regime(fluctuation, amp_ratio, prominence, rate_drive=np.nan,
                    activity_floor=ACTIVITY_FLOOR,
                    stationary_min_ratio=STATIONARY_MIN_RATIO,
                    prominence_min=PEAK_PROMINENCE_MIN,
                    rate_drive_min=RATE_DRIVE_MIN):
    """Label the dynamical regime a signal is in.

    The order of the tests is the point - each one only makes sense once the previous
    has passed:

      diverged     non-finite values (the Euler integration blew up)
      pinned       the firing rate no longer moves - jammed at the sigmoid ceiling or
                   silenced at zero - so the population cannot take part in any loop
      fixed_point  no fluctuation at all; the spectrum below would be numerical residue
      damped       fluctuating, but the amplitude is still decaying across the run, so
                   what is measured is the ring-down of the initial transient
      oscillation  sustained, with a spectral peak above the 1/f background
      noise_driven sustained, but spectrally featureless: the background noise driving
                   a stable network, not a rhythm of its own

    With the background noise on (Ib_noise_std > 0), 'oscillation' means the network
    has a *resonance* the noise is exciting - which is what an EEG rhythm looks like -
    not necessarily a self-sustained one. Re-running the same sweep with
    `--noise 0` separates the two: only a genuine limit cycle survives there, and
    everything else falls back to 'fixed_point' or 'damped'.
    """
    if not np.isfinite(fluctuation):
        return "diverged"
    if np.isfinite(rate_drive) and rate_drive < rate_drive_min:
        return "pinned"
    if fluctuation <= activity_floor:
        return "fixed_point"
    if np.isfinite(amp_ratio) and amp_ratio < stationary_min_ratio:
        return "damped"
    if np.isfinite(prominence) and prominence >= prominence_min:
        return "oscillation"
    return "noise_driven"


def oscillation_features(signal, step_size, seg_dur=1.0, overlap=0.5, settle_s=2.0,
                         fmin=FMIN, fmax=FMAX, rate_drive=np.nan, psd_null=None,
                         return_psd=False):
    """Full set of ongoing-oscillation measures for one simulated signal.

    Args:
        signal:    1-D array (population potential, firing rate, or ROI dipole).
        step_size: sample spacing in s.
        seg_dur, overlap, settle_s, fmin, fmax: see `welch_spectrum`.
        rate_drive: std of the population's firing rate as a fraction of its sigmoid
                   ceiling, if known (see RATE_DRIVE_MIN). nan disables the test.
        psd_null:  power spectrum of the same signal in the disconnected network, on
                   the same frequency grid. When given, the band measures are taken on
                   the network gain spectrum psd/psd_null instead of on the raw
                   spectrum - see `network_gain_spectrum` for why that matters.
        return_psd: also return (freqs, psd) alongside the feature dict.

    Returns:
        dict of features (see keys below), or (features, freqs, psd) when
        `return_psd`. All spectral features are nan when the run diverged or the
        settled window is too short.

        amplitude       std of the settled signal (absolute units)
        mean_level      mean of the settled signal (the DC operating point)
        fluctuation     std / max|x| of the settled signal
        amp_ratio       late/early segment amplitude: ~1 sustained, <<1 ringing down
        band            band holding the strongest peak ('none' if the regime is not
                        'oscillation')
        peak_freq       that peak's frequency (Hz)
        peak_prominence its height over the band's local flanks, in log10 power ratio
        peak_power      raw power at peak_freq
        raw_peak_freq   argmax of the raw spectrum (for reference; 1/f-dominated)
        {theta,alpha,beta,gamma}_prom / _freq
                        the same measure per band, so alpha can be read off directly
        scored_vs       'loopfree' if the band measures are network gains over the
                        loop-free control run, 'flat' if they had to be taken on the
                        raw spectrum (no usable control)
        regime          one of REGIMES
    """
    x = np.asarray(signal, dtype=float)
    keys = ["amplitude", "mean_level", "fluctuation", "amp_ratio", "peak_freq",
            "peak_prominence", "peak_power", "raw_peak_freq"]
    keys += [f"{name}_{suffix}" for name, _, _ in BANDS for suffix in ("prom", "freq")]
    feats = {k: np.nan for k in keys}
    feats["regime"] = "diverged"
    feats["band"] = "none"
    feats["scored_vs"] = "none"

    if x.size == 0 or not np.all(np.isfinite(x)):
        return (feats, np.array([]), np.array([])) if return_psd else feats

    settled = _settled(x, step_size, settle_s)
    feats["amplitude"] = float(settled.std())
    feats["mean_level"] = float(settled.mean())
    feats["fluctuation"] = relative_fluctuation(x, step_size, settle_s)

    freqs, psd, seg_stds = welch_spectrum(x, step_size, seg_dur, overlap, settle_s,
                                          fmin, fmax)
    feats["amp_ratio"] = amplitude_ratio(seg_stds)

    band_peak = "none"
    if len(freqs):
        gain = network_gain_spectrum(psd, psd_null)
        # not the string "null": pandas.read_csv parses that back as NaN
        feats["scored_vs"] = "flat" if gain is None else "loopfree"
        proms = band_prominences(freqs, psd if gain is None else gain)
        for name, (prom, freq) in proms.items():
            feats[f"{name}_prom"], feats[f"{name}_freq"] = prom, freq
        band_peak, peak_freq, prom = dominant_peak(proms)
        feats["peak_freq"], feats["peak_prominence"] = peak_freq, prom
        if np.isfinite(peak_freq):
            feats["peak_power"] = float(psd[int(np.argmin(np.abs(freqs - peak_freq)))])
        feats["raw_peak_freq"] = float(freqs[int(np.argmax(psd))])

    feats["regime"] = classify_regime(feats["fluctuation"], feats["amp_ratio"],
                                      feats["peak_prominence"], rate_drive)
    feats["band"] = band_peak if feats["regime"] == "oscillation" else "none"

    return (feats, freqs, psd) if return_psd else feats
