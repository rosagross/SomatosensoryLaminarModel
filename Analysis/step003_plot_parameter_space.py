"""
File: step003_plot_parameter_space.py
Author: Rosa Grossmann
Description: Overview figures of the model's parameter space: where it oscillates,
    in which band, and at what frequency - with the alpha band called out separately.

    Input: the sweeps written by Simulations/run_parameter_sweep.py under
    $SIMDIR/parameter_space/<tag>/<sweep_name>/. Run those first, e.g.

        python Simulations/run_parameter_sweep.py grid   --n 21
        python Simulations/run_parameter_sweep.py line   --n 25
        python Simulations/run_parameter_sweep.py screen --n-samples 2000

    Output: $SIMDIR/Figures/parameter_space/<tag>/, one set of figures per sweep.

    Per 2-D grid sweep:
      *_regimemap    which dynamical regime each parameter set is in (read this first)
      *_bandmap      which band the dominant peak falls in
      *_peak_freq    peak frequency, masked to the oscillating cells
      *_alpha_prom   how far the 8-13 Hz peak rises above the loop-free control
      *_amplitude    size of the ongoing fluctuation
      *_<area>_*     the same, resolved by layer and cell type
    Per 1-D line sweep:
      *_line_peak_freq / *_line_alpha_prom   one-at-a-time sensitivity of each parameter
      *_line_regimes                         the regime along each parameter
      *_gainspectra_<param>                  spectra along the sweep, as a waterfall
    Per Latin-hypercube screen:
      *_screen_alpha_prom / *_screen_peak_freq   marginals against every parameter
      *_screen_summary                           regime and band composition
      top_alpha_<signal>.csv                     the strongest alpha parameter sets

Usage:
    python Analysis/step003_plot_parameter_space.py                # every sweep found
    python Analysis/step003_plot_parameter_space.py --tag noise0   # one tag
    python Analysis/step003_plot_parameter_space.py --sweeps grid_e1_tau_vs_e2_tau
"""

# %%
import argparse
import glob
import os
import sys

import matplotlib
import pandas as pd

if not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")     # these sweeps are usually plotted on a headless node

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import plotting_functions_parameter_space as pps  # noqa: E402

SIMDIR = os.getenv("SIMDIR", "/data/pt_02989")
SWEEP_ROOT = os.path.join(SIMDIR, "parameter_space")
FIGURE_ROOT = os.path.join(SIMDIR, "Figures", "parameter_space")

# The continuous measures mapped for every grid sweep: (column, colour map, label,
# diverging?). alpha_prom and peak_prominence are signed - a negative value is a dip
# where a peak was looked for - so they get a diverging map centred on zero.
GRID_VALUES = [
    ("peak_freq", "magma", "peak frequency (Hz)", False),
    ("alpha_prom", "Blues", "alpha prominence (log10 gain)", True),
    ("peak_prominence", "Reds", "peak prominence (log10 gain)", True),
    ("amplitude", "viridis", "fluctuation amplitude", False),
]


# %%
def plot_grid_sweep(sweep_dir, figure_dir, name):
    """Every figure of one 2-D grid sweep."""
    df, _ = pps.load_sweep(sweep_dir)
    x, y = df["x_param"].iloc[0], df["y_param"].iloc[0]
    made = []

    made.append(pps.plot_regime_map(df, x, y, figure_dir, name, kind="regime"))
    made.append(pps.plot_regime_map(df, x, y, figure_dir, name, kind="band"))
    for value, cmap, label, diverging in GRID_VALUES:
        made.append(pps.plot_value_map(df, x, y, value, figure_dir, name,
                                       cmap=cmap, label=label, diverging=diverging,
                                       # amplitude is meaningful everywhere, the
                                       # spectral measures only where there is a peak
                                       mask_non_oscillating=(value != "amplitude")))
    # laminar detail: the ROI dipoles average the layers together
    for area in ("A3b", "S1", "S2"):
        made.append(pps.plot_population_grid_map(df, x, y, "regime", area, figure_dir,
                                                 name, categorical="regime",
                                                 label="regime"))
        made.append(pps.plot_population_grid_map(df, x, y, "peak_freq", area,
                                                 figure_dir, name,
                                                 label="peak frequency (Hz)"))
        made.append(pps.plot_population_grid_map(df, x, y, "alpha_prom", area,
                                                 figure_dir, name, cmap="RdBu_r",
                                                 label="alpha prominence"))
    return made


def plot_line_sweep(sweep_dir, figure_dir, name, focus="roi_A1"):
    """Every figure of the one-at-a-time line sweep."""
    df, cfg = pps.load_sweep(sweep_dir)
    # the point the lines pass through: recorded explicitly by newer sweeps (it can be
    # moved off the defaults with --center), otherwise read back from the base params
    reference = cfg.get("reference") or {
        p: cfg["base_params"].get(p) for p in pps.PARAM_LABELS
        if p in cfg.get("base_params", {})}
    made = [
        pps.plot_line_sensitivity(df, "peak_freq", figure_dir, name,
                                  label="peak frequency (Hz)",
                                  mask_non_oscillating=True, reference=reference),
        pps.plot_line_sensitivity(df, "alpha_prom", figure_dir, name,
                                  label="alpha prominence (log10 gain)",
                                  reference=reference),
        pps.plot_line_sensitivity(df, "amplitude", figure_dir, name,
                                  label="fluctuation amplitude", reference=reference),
    ]
    for signal in pps.DEFAULT_SIGNALS:
        made.append(pps.plot_line_regimes(df, figure_dir, name, signal=signal))

    # spectral waterfalls along the parameters that move the peak the most
    try:
        spec = pps.load_spectra(sweep_dir)
    except (OSError, KeyError):
        return made
    # waterfalls for the parameters that move the focus ROI's alpha peak the most
    spread = (df[df["signal"] == focus].groupby("varied_param")["alpha_prom"]
              .agg(lambda s: s.max() - s.min()).sort_values(ascending=False))
    for param in spread.head(4).index:
        sub = df[(df["varied_param"] == param) & (df["signal"] == focus)]
        # the waterfall needs the sweep's own simulations, in parameter order
        sub_spec = dict(spec)
        idx = sub["sim_id"].to_numpy()
        sub_spec["psd"] = spec["psd"][idx]
        sub_spec["psd_null"] = spec["psd_null"][idx]
        sub_spec["theta"] = spec["theta"][idx]
        made.append(pps.plot_gain_waterfall(sub_spec, param, figure_dir, name,
                                            signal=focus))
    return made


def plot_screen_sweep(sweep_dir, figure_dir, name):
    """Every figure of the Latin-hypercube screen, plus the ranked parameter sets."""
    df, _ = pps.load_sweep(sweep_dir)
    made = [pps.plot_screen_summary(df, figure_dir, name)]
    for signal in pps.DEFAULT_SIGNALS:
        made.append(pps.plot_screen_marginals(df, "alpha_prom", figure_dir, name,
                                              signal=signal,
                                              label="alpha prominence (log10 gain)"))
        made.append(pps.plot_screen_marginals(df, "peak_freq", figure_dir, name,
                                              signal=signal,
                                              label="peak frequency (Hz)"))
        top = pps.top_oscillatory_points(df, signal=signal, value="alpha_prom", n=25)
        out = os.path.join(figure_dir, f"top_alpha_{signal}.csv")
        top.to_csv(out, index=False)
        made.append(out)
        if not top.empty:
            print(f"\n  strongest alpha parameter sets for {signal}:")
            print(top.head(5).to_string(index=False))
    return made


PLOTTERS = {"grid": plot_grid_sweep, "line": plot_line_sweep,
            "screen": plot_screen_sweep}


# %%
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sweep-root", default=SWEEP_ROOT)
    parser.add_argument("--figure-root", default=FIGURE_ROOT)
    parser.add_argument("--tag", default=None,
                        help="only this sweep tag (default: every tag found)")
    parser.add_argument("--sweeps", nargs="+", default=None,
                        help="only these sweep directory names")
    parser.add_argument("--focus", default="roi_A1", choices=pps.DEFAULT_SIGNALS,
                        help="ROI the spectral waterfalls of a line sweep follow")
    args = parser.parse_args(argv)

    pps.set_style()
    pattern = os.path.join(args.sweep_root, args.tag or "*", "*", "sweep_features.csv")
    sweeps = sorted(glob.glob(pattern))
    if not sweeps:
        raise SystemExit(f"no sweeps found under {args.sweep_root} - run "
                         f"Simulations/run_parameter_sweep.py first")

    n_figures = 0
    for features in sweeps:
        sweep_dir = os.path.dirname(features)
        name = os.path.basename(sweep_dir)
        tag = os.path.basename(os.path.dirname(sweep_dir))
        if args.sweeps and name not in args.sweeps:
            continue
        _, cfg = pps.load_sweep(sweep_dir)
        mode = cfg.get("mode", "grid")
        figure_dir = os.path.join(args.figure_root, tag)
        print(f"\n=== {tag}/{name}  ({mode}, {cfg.get('n_points')} parameter sets) ===")
        kwargs = {"focus": args.focus} if mode == "line" else {}
        made = PLOTTERS[mode](sweep_dir, figure_dir, name, **kwargs)
        n_figures += len(made)
        print(f"  {len(made)} outputs -> {figure_dir}")

    print(f"\ndone: {n_figures} figures/tables written under {args.figure_root}")


if __name__ == "__main__":
    main()
