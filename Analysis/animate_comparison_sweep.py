"""
Animate saved TF / time-course comparison runs across a parameter sweep.

For a chosen swept parameter (coupling_strength or strength_I) it groups the
saved per-run HDF5 files by holding all *other* parameters fixed, then builds one
animation per (comparison type x parameter x ROI). Each frame shows the SIMULATED
map/trace at one parameter value (value shown in the title).

Mirrors Analysis/animate_frequency_spectra.py (glob HDF5 -> sort by parameter ->
FuncAnimation + PillowWriter -> gif).
"""
# %%
import argparse
import glob
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

SIMDIR = "/data/p_02989/Modelling/output_grossmannr/"
print('simdir', SIMDIR)
ROIS = ("A3b", "A1", "S2")
SWEEP_PARAMS = ("coupling_strength", "strength_I")

# Per comparison type: input folder name, file prefix, and the error attr stored in HDF5.
COMPARISON_CONFIG = {
    "tf": {"dirname": "timefreq_comparison", "prefix": "tf_comparison", "error_attr": "tf_error"},
    "tc": {"dirname": "timecourse_comparison", "prefix": "tc_comparison", "error_attr": "tc_error"},
}

# All parameters stored on each run; used to decide which are "held fixed" for a sweep.
ALL_PARAMS = (
    "g_intercortical", "coupling_strength", "strength_I", "Ib_strength",
    "Iext_strength", "Iext_duration", "step_size", "input_onset", "input_type", "area",
)

# Display window: -200..+250 ms; the trim mask is derived per file from its time axis.
PLOT_TMIN_MS = -200.0
PLOT_TMAX_MS = 250.0
EPS = 1e-10


def load_comparison_file(filepath, data_type):
    """Load one run's simulated ROI maps/traces + time axis + parameter attrs."""
    with h5py.File(filepath, "r") as h5f:
        sim = {roi: h5f["sim"][roi][:] for roi in ROIS}
        times_ms = h5f["times_ms"][:]
        params = {}
        for key in ALL_PARAMS:
            if key in h5f.attrs:
                value = h5f.attrs[key]
                params[key] = value.decode("utf-8") if isinstance(value, bytes) else value
        error = float(h5f.attrs[COMPARISON_CONFIG[data_type]["error_attr"]])
    return {"filepath": filepath, "sim": sim, "times_ms": times_ms, "params": params, "error": error}


def select_sweep_group(records, sweep_param):
    """
    Group runs by the tuple of all parameters except `sweep_param`, then return the
    group with the most distinct sweep-parameter values, sorted by that value.
    """
    fixed_keys = [k for k in ALL_PARAMS if k != sweep_param]
    groups = {}
    for rec in records:
        if sweep_param not in rec["params"]:
            continue
        key = tuple(rec["params"].get(k) for k in fixed_keys)
        groups.setdefault(key, []).append(rec)

    best_key, best_group = None, []
    for key, group in groups.items():
        n_distinct = len({r["params"][sweep_param] for r in group})
        if n_distinct > len({r["params"][sweep_param] for r in best_group}):
            best_key, best_group = key, group

    if len({r["params"][sweep_param] for r in best_group}) < 2:
        raise ValueError(
            f"Need >=2 distinct '{sweep_param}' values with all other parameters fixed; "
            f"found at most {len({r['params'][sweep_param] for r in best_group})}. "
            f"Run a sweep over '{sweep_param}' first."
        )

    fixed = dict(zip(fixed_keys, best_key))
    print(f"  [{sweep_param}] fixed params: {fixed}")
    # de-duplicate per sweep value (keep last written), then sort
    by_value = {r["params"][sweep_param]: r for r in best_group}
    return [by_value[v] for v in sorted(by_value)]


def build_interpolated_frames(records, interpolation_steps):
    """Interpolated in-between frames between consecutive parameter values (from animate_frequency_spectra)."""
    frames = []
    if len(records) == 1:
        return [(records[0], records[0], 0.0)]
    for i in range(len(records) - 1):
        for step in range(interpolation_steps):
            alpha = step / interpolation_steps
            frames.append((records[i], records[i + 1], alpha))
    frames.append((records[-1], records[-1], 0.0))
    return frames


def animate_tf(records, roi, sweep_param, output_path, interpolation_steps, fps):
    """Animate the simulated TF map for one ROI across the swept parameter."""
    times_ms = records[0]["times_ms"]
    mask = (times_ms >= PLOT_TMIN_MS) & (times_ms <= PLOT_TMAX_MS)
    t_ms = times_ms[mask]
    freqs = np.arange(1, 41, 1).astype(float)

    maps = [rec["sim"][roi][:, mask] for rec in records]
    vmax = max(float(np.percentile(m, 95)) for m in maps)
    vmax = vmax if vmax > 0 else 1.0
    frames = build_interpolated_frames(records, interpolation_steps)

    fig, ax = plt.subplots(figsize=(7, 4))
    im = ax.imshow(
        maps[0], aspect="auto", origin="lower",
        extent=[t_ms[0], t_ms[-1], freqs[0], freqs[-1]],
        cmap="hot_r", vmin=0.0, vmax=vmax,
    )
    ax.axvline(0, color="white", lw=0.8, ls="--")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Frequency (Hz)")
    fig.colorbar(im, ax=ax, label="Power (a.u.²)", fraction=0.046, pad=0.04)

    def update(frame_data):
        left, right, alpha = frame_data
        m = (1.0 - alpha) * left["sim"][roi][:, mask] + alpha * right["sim"][roi][:, mask]
        value = (1.0 - alpha) * left["params"][sweep_param] + alpha * right["params"][sweep_param]
        im.set_data(m)
        ax.set_title(f"{roi} time-frequency - {sweep_param}={value:.3g}")
        return (im,)

    ani = FuncAnimation(fig, update, frames=frames,
                        interval=max(10, int(1000 / max(fps, 1))), blit=False, repeat=True)
    ani.save(output_path, writer=PillowWriter(fps=fps))
    plt.close(fig)
    print(f"  Saved animation: {output_path}")


def animate_tc(records, roi, sweep_param, output_path, interpolation_steps, fps):
    """Animate the simulated, peak-normalized TC trace for one ROI across the swept parameter."""
    times_ms = records[0]["times_ms"]
    mask = (times_ms >= PLOT_TMIN_MS) & (times_ms <= PLOT_TMAX_MS)
    t_ms = times_ms[mask]

    def peaknorm(rec):
        x = rec["sim"][roi][mask]
        return x / (np.max(np.abs(x)) + EPS)

    frames = build_interpolated_frames(records, interpolation_steps)

    fig, ax = plt.subplots(figsize=(7, 3.4))
    line, = ax.plot([], [], lw=1.5, color="C1")
    ax.set_xlim(t_ms[0], t_ms[-1])
    ax.set_ylim(-1.1, 1.1)
    ax.axvline(0, color="grey", lw=0.8, ls="--")
    ax.axhline(0, color="k", lw=0.5, alpha=0.3)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Peak-normalized amplitude")

    def update(frame_data):
        left, right, alpha = frame_data
        y = (1.0 - alpha) * peaknorm(left) + alpha * peaknorm(right)
        value = (1.0 - alpha) * left["params"][sweep_param] + alpha * right["params"][sweep_param]
        line.set_data(t_ms, y)
        ax.set_title(f"{roi} timecourse  {sweep_param}={value:.3g}")
        return (line,)

    ani = FuncAnimation(fig, update, frames=frames,
                        interval=max(10, int(1000 / max(fps, 1))), blit=True, repeat=True)
    ani.save(output_path, writer=PillowWriter(fps=fps))
    plt.close(fig)
    print(f"  Saved animation: {output_path}")


# %%
def main():
    parser = argparse.ArgumentParser(description="Animate TF/TC comparison runs across a parameter sweep.")
    parser.add_argument("--data-type", choices=["tf", "tc", "both"], default="both")
    parser.add_argument("--param", choices=["coupling_strength", "strength_I", "both"], default="both")
    parser.add_argument("--roi", choices=["A3b", "A1", "S2", "all"], default="all")
    parser.add_argument("--interpolation-steps", type=int, default=8)
    parser.add_argument("--fps", type=int, default=12)
    args = parser.parse_args()

    data_types = ["tf", "tc"] if args.data_type == "both" else [args.data_type]
    params = list(SWEEP_PARAMS) if args.param == "both" else [args.param]
    rois = list(ROIS) if args.roi == "all" else [args.roi]
    animate_fn = {"tf": animate_tf, "tc": animate_tc}

    for data_type in data_types:
        cfg = COMPARISON_CONFIG[data_type]
        input_dir = os.path.join(SIMDIR, cfg["dirname"])
        filepaths = sorted(glob.glob(os.path.join(input_dir, f"{cfg['prefix']}_*.hdf5")))
        if not filepaths:
            raise FileNotFoundError(f"No {data_type.upper()} runs found in {input_dir}")
        records = [load_comparison_file(p, data_type) for p in filepaths]

        out_dir = os.path.join(input_dir, "animations")
        os.makedirs(out_dir, exist_ok=True)

        for sweep_param in params:
            print(f"[{data_type.upper()}] sweeping {sweep_param} ...")
            group = select_sweep_group(records, sweep_param)
            for roi in rois:
                out_path = os.path.join(out_dir, f"{data_type}_sweep-{sweep_param}_roi-{roi}.gif")
                animate_fn[data_type](
                    group, roi, sweep_param, out_path,
                    args.interpolation_steps, args.fps,
                )


if __name__ == "__main__":
    main()

# %%
