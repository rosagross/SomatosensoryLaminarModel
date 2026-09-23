"""
File: step004_1_dipole_oscillation_peakfreq.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-08-31
Description: File for analysing oscillatory behaviour of signal in dipole form.

Looking at the ongoing (stimulus-free) spectrum of each ROI dipole:
1. What is the oscillation peak frequency?
2. What is the alpha envelope?
3. Which dynamical regime / band is the ROI in at all?

Driven by step004_run_model_analysis.py, which runs the sweeps; this file only reads the
tables that sweep wrote and plots them. It never simulates.

Two things this analysis deliberately does NOT do, both of which cost a day to learn:

- It does not score a peak on the raw spectrum. A network with every loop cut already
  reports an "alpha peak" of 0.22 and a "gamma peak" of 0.30, purely from the
  Ornstein-Uhlenbeck background whose corner sits at 1/(2*pi*Ib_noise_tau) ~ 9.9 Hz -
  inside the alpha band - plus the synaptic kernel roll-off. Every peak here is scored on
  the network gain spectrum: the run divided by the same parameter set with its coupling
  gains scaled down by 1e-3 (run_parameter_sweep._null_reference). Control / itself scores
  0.000 in every band.
- It does not fit a global 1/f slope. These spectra are not power laws; one log-log line
  through the whole range leaves a broad positive bow that reads as a ~13 Hz "peak" in
  every single run whatever the parameters (see the header of
  Simulations/model/oscillation_metrics.py). The aperiodic part is removed locally
  instead, by referencing each band against flanks either side of it -
  signal_preprocessing.spectral_prominence, which oscillation_metrics.band_prominences
  already calls.

Both are done inside run_parameter_sweep, so the columns read here (`peak_freq`, `band`,
`alpha_prom`, `sig_alpha_gain`) already carry them.

Parameter combinations are defined in step004_run_model_analysis.py.
"""

import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import plotting_functions_step004 as pf                            # noqa: E402


# (column, cmap, label, kwargs) for the continuous ROI maps.
ROI_VALUES = [
    ("peak_freq", "magma",
     "peak frequency (Hz)", dict(mask_non_oscillating=True)),
    ("alpha_prom", "RdBu_r",
     "alpha prominence (log10 gain)", dict(mask_non_oscillating=True, diverging=True)),
    ("peak_prominence", "Reds",
     "peak prominence (log10 gain)", dict(mask_non_oscillating=True, diverging=True)),
    # the alpha envelope of the ROI dipole itself, referred to the loop-free control.
    # The raw envelope (sig_alpha_amp) is not plotted: with the OU corner inside the alpha
    # band roughly half of any fluctuation lands there whatever the network does, so the
    # raw number tracks the total fluctuation and not the rhythm.
    ("sig_alpha_gain", "rocket",
     "alpha envelope / loop-free control", dict(mask_non_oscillating=False,
                                                log_scale=True)),
    ("amplitude", "viridis",
     "ongoing fluctuation (SD)", dict(mask_non_oscillating=False)),
]


def run_dipole_oscillation_analysis(sweep_dir, figure_dir, name):
    """Every ROI-dipole figure of one 2-D grid sweep. Returns the paths written."""
    df, _ = pf.load_sweep(sweep_dir)
    df = df[df["signal_kind"] == "dipole"]
    if df.empty:
        print(f"    {name}: no dipole rows (sweep run with --no-dipoles?) - skipped")
        return []

    x, y = df["x_param"].iloc[0], df["y_param"].iloc[0]
    os.makedirs(figure_dir, exist_ok=True)
    made = []

    # read the regime map first: peak frequency is only meaningful where there is a peak
    for kind in ("regime", "band"):
        made.append(pf.plot_roi_map(df, x, y, kind, figure_dir, f"{name}_roi_{kind}",
                                    categorical=kind))

    for value, cmap, label, kwargs in ROI_VALUES:
        if value not in df.columns:
            continue
        made.append(pf.plot_roi_map(df, x, y, value, figure_dir, f"{name}_roi_{value}",
                                    cmap=cmap, label=label, **kwargs))

    # the line companions: a few representative traces rather than the whole dense grid
    made.append(pf.plot_grid_slices(df, x, y, "peak_freq", figure_dir,
                                    f"{name}_roi_slices_peak_freq",
                                    label="peak frequency (Hz)",
                                    mask_non_oscillating=True))
    if "sig_alpha_gain" in df.columns:
        made.append(pf.plot_grid_slices(df, x, y, "sig_alpha_gain", figure_dir,
                                        f"{name}_roi_slices_alpha_gain",
                                        label="alpha envelope / control",
                                        log_scale=True))

    # and the spectra themselves, so a peak can be checked by eye
    made.append(pf.plot_gain_spectra(pf.load_spectra(sweep_dir), df, x, y, figure_dir,
                                     f"{name}_roi_gainspectra"))
    return [m for m in made if m]
