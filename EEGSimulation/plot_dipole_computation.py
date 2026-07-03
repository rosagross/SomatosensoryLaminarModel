"""
File: plot_dipole_computation.py
Description:
    Demonstrate how the EEG-scale current dipole is computed from the simulated
    laminar membrane potentials (see SomatoModel.compute_dipoles).

    Runs one stimulated simulation, computes the dipole for a subject forward
    model, and produces two publication-style figures:

      Figure 1 - dipole parameters per cell:
          the geometric dipole model (signed length x orientation) assigned to
          every population, one panel per area. Reuses plot_dipole_parameters().

      Figure 2 - S2 dipole computation example:
          the full pipeline for area S2 - the excitatory membrane potential
          feeding each layer (top), and each layer's resulting dipole trace plus
          the summed S2 dipole (bottom).

    Run with the project env (WDDIR/SIMDIR/DATADIR/SUBJECTS_DIR must be exported):
        python EEGSimulation/plot_dipole_computation.py
"""
# %%
import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- make the project modules importable (mirror plot_dipole_parameters.py) ---
_eeg_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.dirname(_eeg_dir)
for _p in [_eeg_dir,
           os.path.join(_project_root, 'Simulations'),
           os.path.join(_project_root, 'Simulations', 'model')]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

from somato_model import SomatoModel, read_simulation_params
from plot_dipole_parameters import plot_dipole_parameters
from plotting_style import figure_style

colors, _ = figure_style()

# output directory
SIMDIR = os.getenv("SIMDIR")
figure_dir = os.path.join(SIMDIR, "Figures", "dipole_computation")
os.makedirs(figure_dir, exist_ok=True)

# subject forward model(s) used to project the dipoles (as in simulation_main.py)
subjects = [15]


# %%
# ---------------------------------------------------------------------------
# 1) Run one stimulated simulation
# ---------------------------------------------------------------------------
params = read_simulation_params()

# override for a single stimulated run (same values as the simulation_main.py loop)
params['g_intercortical'] = 2.0
params['coupling_strength'] = 11.24
params['strength_I'] = 0.593
params['Iext_duration'] = 0.044
params['Iext_strength'] = 80
params['Ib_strength'] = 3
params['area'] = 'all'
params['resistance_factor'] = 1
params['g_thal'] = 2
params['g_thalPOm'] = 1
params['sI_thal'] = 0.5
params['delay_factor'] = 0.005
params['extI_cellcounts'] = 1000
params['bI_cellcounts'] = 100

model = SomatoModel(params)
model.simulate()

# ---------------------------------------------------------------------------
# 2) Compute the dipole from the simulated potentials
#    simDipoles shape (9, ntimes): 0 = A3b, 1-4 = A1 L{1,4,5,6}_E, 5-8 = S2 L{1,4,5,6}_E
# ---------------------------------------------------------------------------
simDipoles = model.compute_dipoles(subjects)


# %%
# ---------------------------------------------------------------------------
# Figure 1: dipole parameters per cell (publication / poster ready)
# ---------------------------------------------------------------------------
# Cell-type colours (shared with EEGSimulation/plot_dipole_parameters.py)
_CELL_COLORS = {'E': '#4477AA', 'PV': '#EE6677', 'SST': '#228833', 'VIP': '#CCBB44'}


def plot_dipole_parameters_clean(json_path, figure_dir):
    """
    Publication/poster-ready view of the dipole model parameters.

    The dipole model factorises into two independent parameters:
      - dipole LENGTH: depends on the cortical layer (and area A3b),
      - ORIENTATION:   depends on the source cell type (sign of the current).

    Two panels:
      (left)  dipole length per layer/area (mm),
      (right) orientation per cell type (inward -1 / none 0 / outward +1).
    """
    with open(json_path, 'r') as f:
        dp = json.load(f)

    # --- length per layer/area (values are uniform within each array, VIP = 0) ---
    def _repr_len(arr):
        arr = np.asarray(arr, dtype=float)
        nz = arr[arr != 0]
        return float(nz.max()) if nz.size else 0.0

    area_labels = ['A3b', 'L2/3', 'L4', 'L5', 'L6']
    lengths = [_repr_len(dp['dipole_lengths']['A3b'])] + \
              [_repr_len(a) for a in dp['dipole_lengths']['S2']]  # S1 == S2

    # --- orientation per source cell type (from the A3b block: E, PV, SST, VIP) ---
    cell_types = ['E', 'PV', 'SST', 'VIP']
    orient = [int(v) for v in dp['dipole_orientation']['A3b'][:4]]

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8),
                             gridspec_kw={'width_ratios': [1.25, 1]})

    # ---- Panel A: dipole length by layer ----
    ax = axes[0]
    depth_colors = sns.color_palette('mako_r', len(area_labels))
    bars = ax.bar(area_labels, lengths, color=depth_colors,
                  edgecolor='black', linewidth=0.6, width=0.72)
    for b, v in zip(bars, lengths):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.03, f'{v:g}',
                ha='center', va='bottom', fontsize=10)
    ax.set_ylabel('Dipole length (mm)', fontsize=12)
    ax.set_ylim(0, max(lengths) * 1.18)
    ax.set_title('Dipole length by layer', fontsize=13, fontweight='bold')
    ax.tick_params(axis='both', labelsize=11)
    ax.margins(x=0.02)

    # ---- Panel B: orientation by cell type ----
    ax = axes[1]
    bar_colors = [_CELL_COLORS[c] for c in cell_types]
    bars = ax.bar(cell_types, orient, color=bar_colors,
                  edgecolor='black', linewidth=0.6, width=0.72)
    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_ylim(-1.45, 1.45)
    ax.set_yticks([-1, 0, 1])
    ax.set_yticklabels(['inward\n(−1)', 'none\n(0)', 'outward\n(+1)'], fontsize=10)
    ax.set_title('Orientation by cell type', fontsize=13, fontweight='bold')
    ax.tick_params(axis='x', labelsize=11)
    # label the zero-height VIP bar so it is not invisible
    for b, v in zip(bars, orient):
        if v == 0:
            ax.text(b.get_x() + b.get_width() / 2, 0.06, '0',
                    ha='center', va='bottom', fontsize=10, color='0.35')

    fig.suptitle('Current-dipole model parameters', fontsize=14, fontweight='bold')
    sns.despine(fig=fig, trim=True)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    name = os.path.join(figure_dir, 'dipoleParameters_perCell')
    fig.savefig(name + '.pdf', bbox_inches='tight')
    fig.savefig(name + '.png', dpi=300, bbox_inches='tight')
    print(f'Figure 1 saved to {name}.pdf/.png')
    return fig


fig1 = plot_dipole_parameters_clean(
    os.path.join(_eeg_dir, 'dipole_parameters.json'), figure_dir)


# %%
# ---------------------------------------------------------------------------
# Figure 2: S2 dipole computation example (potentials -> dipoles -> sum)
# ---------------------------------------------------------------------------
def plot_dipole_computation_S2(model, simDipoles, figure_dir):
    """
    Demonstrate the dipole computation for area S2.

    Top panel  : the excitatory membrane potential feeding each S2 layer
                 (input to the dipole projection).
    Bottom panel: the per-layer dipole trace (potential weighted by the layer's
                 dipole model) and the summed S2 dipole (bold black).

    Parameters
    ----------
    model : SomatoModel
        A simulated model (model.potential populated).
    simDipoles : np.ndarray, shape (9, ntimes)
        Output of model.compute_dipoles(...). Indices 5-8 are the S2 layers.
    figure_dir : str
        Where to save the figure.
    """
    step_size = model.step_size
    input_onset = model.input_onset
    Iext_dur = model.Iext_duration

    # S2 excitatory populations, one per cortical layer (see get_population_mapping)
    s2 = model.get_population_mapping()['S2']
    s2_e_idx = [s2['L1_E'], s2['L4_E'], s2['L5_E'], s2['L6_E']]
    layer_labels = ['L2/3', 'L4', 'L5', 'L6']
    # qualitative palette: maximally distinct, colourblind-friendly colour per layer
    layer_colors = ['#4477AA', '#EE6677', '#228833', '#AA3377']  # blue, red, green, purple

    # S2 dipole traces live at simDipoles rows 5..8
    s2_dipoles = simDipoles[5:9]
    s2_sum = s2_dipoles.sum(axis=0)

    ntimes = simDipoles.shape[1]
    t = np.arange(ntimes) * step_size - input_onset  # time re-zeroed to stimulus onset
    window = (-0.05, 0.25)
    mask = (t >= window[0]) & (t <= window[1])

    fig, axes = plt.subplots(2, 1, figsize=(5, 5), sharex=True)

    # --- Top: excitatory membrane potentials feeding each S2 layer ---
    ax = axes[0]
    ax.axvspan(0, Iext_dur, color='0.85', linewidth=0, zorder=0)
    for idx, lab, col in zip(s2_e_idx, layer_labels, layer_colors):
        pot = model.potential[idx].sum(axis=0)  # summed synaptic (membrane) potential
        ax.plot(t[mask], pot[mask], color=col, linewidth=1.2, label=lab)
    ax.set_ylabel('E potential (mV)')
    ax.set_title('S2 excitatory membrane potential (input)')
    ax.legend(title='Layer', loc='upper right', ncol=2, fontsize=7)

    # --- Bottom: per-layer dipole traces + summed S2 dipole ---
    ax = axes[1]
    ax.axvspan(0, Iext_dur, color='0.85', linewidth=0, zorder=0)
    for row, lab, col in zip(s2_dipoles, layer_labels, layer_colors):
        ax.plot(t[mask], row[mask], color=col, linewidth=1.2, label=lab)
    ax.plot(t[mask], s2_sum[mask], color='black', linewidth=2.2, label='Sum (S2)', zorder=5)
    ax.set_ylabel('Dipole (a.u.)')
    ax.set_xlabel('Time from stimulus onset (s)')
    ax.set_title('S2 layer dipoles and their sum (output)')
    ax.legend(loc='upper right', ncol=2, fontsize=7)

    sns.despine(trim=True)
    plt.tight_layout()

    name = os.path.join(figure_dir, 'dipoleComputation_S2_example')
    fig.savefig(name + '.pdf', bbox_inches='tight')
    fig.savefig(name + '.png', dpi=300, bbox_inches='tight')
    print(f'Figure 2 saved to {name}.pdf/.png')
    return fig


fig2 = plot_dipole_computation_S2(model, simDipoles, figure_dir)

plt.show()
