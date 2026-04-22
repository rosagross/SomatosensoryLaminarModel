"""
Animate saved frequency spectra over increasing g_intercortical.
"""
# %%
import argparse
import glob
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

SIMDIR =  "/data/p_02989/Modelling/output_grossmannr/" #os.getenv("WDDIR")
INPUT_DIR = os.path.join(SIMDIR, "spectrum_results")


def load_spectrum_file(filepath):
    with h5py.File(filepath, "r") as h5f:
        freqs = h5f["freqs"][:]
        spectra = h5f["spectra"][:]
        labels_raw = h5f["population_labels"][:]
        population_labels = np.array([lbl.decode("utf-8") for lbl in labels_raw])
        g_intercortical = float(h5f.attrs["g_intercortical"])
        print('g_ionterc', g_intercortical)

    return {
        "filepath": filepath,
        "freqs": freqs,
        "spectra": spectra,
        "population_labels": population_labels,
        "g_intercortical": g_intercortical,
    }


def validate_compatibility(records):
    if not records:
        raise ValueError("No spectrum files found.")

    freqs_ref = records[0]["freqs"]
    labels_ref = records[0]["population_labels"]

    for rec in records[1:]:
        if not np.array_equal(rec["freqs"], freqs_ref):
            raise ValueError(f"Frequency axis mismatch: {rec['filepath']}")
        if not np.array_equal(rec["population_labels"], labels_ref):
            raise ValueError(f"Population labels mismatch: {rec['filepath']}")


def build_interpolated_frames(records, interpolation_steps):
    frames = []

    if len(records) == 1:
        frames.append((records[0], records[0], 0.0))
        return frames

    for i in range(len(records) - 1):
        left = records[i]
        right = records[i + 1]
        for step in range(interpolation_steps):
            alpha = step / interpolation_steps
            frames.append((left, right, alpha))

    frames.append((records[-1], records[-1], 0.0))
    return frames

# %%
def main():
    parser = argparse.ArgumentParser(description="Animate spectra across g_intercortical values.")
    parser.add_argument(
        "--pattern",
        type=str,
        default="spectrum_*.hdf5",
        help="Glob pattern to select spectrum files",
    )
    parser.add_argument(
        "--population",
        type=str,
        default="E1",
        help="Population label to animate",
    )
    parser.add_argument(
        "--max-freq-hz",
        type=float,
        default=100.0,
        help="Maximum displayed frequency",
    )
    parser.add_argument(
        "--interpolation-steps",
        type=int,
        default=8,
        help="Interpolated in-between frames per g_intercortical interval",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=12,
        help="Animation frames per second",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output GIF path. Default: <input-dir>/spectrum_animation_<population>.gif",
    )
    args = parser.parse_args()

    search_pattern = os.path.join(INPUT_DIR, args.pattern)
    filepaths = sorted(glob.glob(search_pattern))
    if len(filepaths) == 0:
        raise FileNotFoundError(f"No files found for pattern: {search_pattern}")

    records = [load_spectrum_file(path) for path in filepaths]
    validate_compatibility(records)
    records = sorted(records, key=lambda r: r["g_intercortical"])

    population_labels = records[0]["population_labels"]
    if args.population not in population_labels:
        raise ValueError(
            f"Population '{args.population}' not found. Available: {', '.join(population_labels)}"
        )
    pop_idx = int(np.where(population_labels == args.population)[0][0])

    freqs = records[0]["freqs"]
    freq_mask = freqs <= args.max_freq_hz
    freqs_plot = freqs[freq_mask]
    if freqs_plot.size == 0:
        raise ValueError(f"No frequencies <= {args.max_freq_hz} Hz available in input files.")

    for rec in records:
        pop_spectrum = rec["spectra"][pop_idx][freq_mask]
        peak_freq = freqs_plot[np.argmax(pop_spectrum)]
        peak_power = np.max(pop_spectrum)
        print(f"g={rec['g_intercortical']:.4f} | peak at {peak_freq:.2f} Hz, power={peak_power:.4f}")

    max_power = max(np.max(rec["spectra"][pop_idx][freq_mask]) for rec in records)
    if max_power <= 0:
        max_power = 1.0

    frames = build_interpolated_frames(records, args.interpolation_steps)

    fig, ax = plt.subplots(figsize=(10, 5))
    line, = ax.plot([], [], lw=2, color="tab:blue")
    ax.set_xlim(0, args.max_freq_hz)
    ax.set_ylim(0, max_power * 1.05)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power")
    ax.grid(True, alpha=0.3)

    def update(frame_data):
        left, right, alpha = frame_data
        y_left = left["spectra"][pop_idx][freq_mask]
        y_right = right["spectra"][pop_idx][freq_mask]
        y_interp = (1.0 - alpha) * y_left + alpha * y_right

        g_interp = (1.0 - alpha) * left["g_intercortical"] + alpha * right["g_intercortical"]
        line.set_data(freqs_plot, y_interp)
        ax.set_title(
            f"Population {args.population} spectrum | g_intercortical={g_interp:.4f}"
        )
        return (line,)

    ani = FuncAnimation(
        fig,
        update,
        frames=frames,
        interval=max(10, int(1000 / max(args.fps, 1))),
        blit=True,
        repeat=True,
    )

    output_path = args.output
    if output_path is None:
        output_path = os.path.join(INPUT_DIR, f"spectrum_animation_{args.population}.gif")
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    ani.save(output_path, writer=PillowWriter(fps=args.fps))
    print(f"Saved animation: {output_path}")


if __name__ == "__main__":
    main()

# %%
