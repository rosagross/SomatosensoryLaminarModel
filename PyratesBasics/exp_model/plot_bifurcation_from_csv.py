import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# USER SETTINGS
# ============================================================
WDDIR = os.getenv("WDDIR")
output_dir = os.path.join(WDDIR, 'PyratesBasics', 'exp_model','complete_model_continuations')
csv_name = "complete_model_bifurcation_sI026_stability_only_g0_10.csv"
csv_path = os.path.join(output_dir, csv_name)

# variable to plot on y-axis
y_var = "E3S2"

# continuation parameter
x_var = "G/g_definition/g_input"

save_path = f"bifurcation_{y_var}_lines.png"

# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(csv_path)

if y_var not in df.columns:
    raise ValueError(f"{y_var} not found in csv.")

if x_var not in df.columns:
    raise ValueError(f"{x_var} not found in csv.")

x = df[x_var].to_numpy()
y = df[y_var].to_numpy()

stability = (
    df["stability"]
    .fillna("unknown")
    .astype(str)
    .str.lower()
    .to_numpy()
)

bifurcation = (
    df["bifurcation"]
    .fillna("")
    .astype(str)
    .to_numpy()
)


# ============================================================
# PLOT
# ============================================================

fig, ax = plt.subplots(figsize=(10, 6))

# ------------------------------------------------------------
# Plot stable/unstable segments as lines
# ------------------------------------------------------------

start = 0

for i in range(1, len(x) + 1):

    end_segment = (
        i == len(x)
        or stability[i] != stability[start]
    )

    if end_segment:

        xs = x[start:i]
        ys = y[start:i]

        stab = stability[start]

        # sort points for smooth line plotting
        idx = np.argsort(xs)
        xs = xs[idx]
        ys = ys[idx]

        # skip tiny segments
        if len(xs) > 1:

            if stab == "stable":

                ax.plot(
                    xs,
                    ys,
                    linestyle="-",
                    linewidth=2,
                    color="tab:blue",
                    label="stable" if start == 0 else None
                )

            elif stab == "unstable":

                ax.plot(
                    xs,
                    ys,
                    linestyle="--",
                    linewidth=2,
                    color="tab:red",
                    label="unstable" if start == 0 else None
                )

            else:

                ax.plot(
                    xs,
                    ys,
                    linestyle=":",
                    linewidth=1.5,
                    color="0.5",
                    label="unknown" if start == 0 else None
                )

        start = i


# ============================================================
# BIFURCATION POINTS
# ============================================================

bif_mask = (
    (bifurcation != "")
    & (bifurcation != "RG")
)

if bif_mask.any():

    ax.scatter(
        x[bif_mask],
        y[bif_mask],
        s=80,
        facecolors="white",
        edgecolors="black",
        linewidths=1.5,
        zorder=5,
        label="bifurcation"
    )

    for xi, yi, lbl in zip(
        x[bif_mask],
        y[bif_mask],
        bifurcation[bif_mask]
    ):

        ax.annotate(
            lbl,
            (xi, yi),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=9
        )


# ============================================================
# FIGURE FORMATTING
# ============================================================

ax.set_xlabel(x_var)
ax.set_ylabel(f"{y_var} [mV]")

ax.set_title(
    f"Bifurcation Diagram: {y_var}"
)

ax.grid(True, alpha=0.25)

# remove duplicate legend entries
handles, labels = ax.get_legend_handles_labels()

unique = dict(zip(labels, handles))

ax.legend(
    unique.values(),
    unique.keys(),
    loc="best"
)

fig.tight_layout()

plt.savefig(save_path, dpi=300)

plt.show()

print(f"Saved figure to: {save_path}")