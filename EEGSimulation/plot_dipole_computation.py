"""
File: plot_dipole_computation.py
Description:
    Demonstrate how the EEG-scale current dipole is computed from the simulated
    laminar membrane potentials (see SomatoModel.compute_dipoles).

    Runs one stimulated simulation, computes the dipole for a subject forward
    model, and produces seven publication-style figures:

      Figure 1 - dipole parameters per cell:
          the geometric dipole model (signed length x orientation) assigned to
          every population, one panel per area (plot_dipole_parameters_clean).

      Figures 2-4 - dipole computation example, one per area (S2, S1/A1, A3b):
          the full pipeline for that area - the excitatory membrane potential
          feeding each layer (top), and each layer's resulting dipole trace plus
          the summed area dipole (bottom). A3b is unlaminated, so it has a single
          population and a single dipole trace.

      Figures 5-7 - source contributions per layer (S2, S1/A1, A3b):
          the dipole model is defined per *input* cell, so PV and SST inhibition
          enter a layer's dipole with opposite orientations and VIP input is
          dipole-silent. The signs come from the parameter file the model loads
          (dipole_parameters_flippedPVSST.json: E -1, PV +1, SST -1, VIP 0) and are
          read from it rather than hard-coded. Structured by source: one column per layer, one row per
          source cell type showing the individual contributions of all 33 source
          populations (colour = source area, its lightness = source layer, line
          style = source cell type, as in plot_all_potentials), plus a summary row
          with the cell-type sums adding up to that layer's dipole. A3b is
          unlaminated, so its figure has a single column.

      Figures 8-10 - dipole parameters per source->target connection, one per area:
          the dipole length and orientation the parameter file assigns to every one of
          the 33 source populations, for each target layer of the area
          (plot_dipole_parameters_per_connection). Figure 1 summarises these as one
          length per layer and one orientation per cell type, which is a
          simplification: the orientation varies by source beyond its cell type (for
          the target A1 L5, S1 L2/3 E is -1 but S1 L4 E is +1, and S2 L4 E is 0).
          These three figures come from the parameter file alone - no simulation and
          no forward model is involved.

    Run with the project env (WDDIR/SIMDIR/DATADIR/SUBJECTS_DIR must be exported):
        python EEGSimulation/plot_dipole_computation.py
"""
# %%
import os
import sys
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

from somato_model import (SomatoModel, read_simulation_params, load_optimized_params,
                          read_dipole_params, DIPOLE_SOURCE_PATHS)
from plotting_style import figure_style
# same source grammar as Simulations/plotting_functions.plot_all_potentials:
# colour = source area, line style = source cell type
from plotting_functions import (SOURCE_AREA_BLOCKS, CELLTYPE_LINESTYLES,
                                _build_source_styles, _plot_connectivity_grid)

try:
    colors, _ = figure_style()
except Exception:
    # figure_style() opens a tkinter window to query the screen size, which fails on
    # headless nodes - the module must stay importable there (only _SOURCE_AREA_COLORS
    # reads `colors`, and it falls back on an empty dict)
    colors = {}


# %%
# ---------------------------------------------------------------------------
# Figure 1: dipole parameters per cell (publication / poster ready)
# ---------------------------------------------------------------------------
# Cell-type colours (shared with EEGSimulation/plot_dipole_parameters.py)
_CELL_COLORS = {'E': '#4477AA', 'PV': '#EE6677', 'SST': '#228833', 'VIP': '#CCBB44'}


def plot_dipole_parameters_clean(json_path, figure_dir=None):
    """
    Publication/poster-ready view of the dipole model parameters.

    The dipole model factorises into two independent parameters:
      - dipole LENGTH: depends on the cortical layer (and area A3b),
      - ORIENTATION:   depends on the source cell type (sign of the current).

    Two panels:
      (left)  dipole length per layer/area (mm),
      (right) orientation per cell type (inward -1 / none 0 / outward +1).
    """
    dp = read_dipole_params(json_path)

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

    if figure_dir is not None:
        name = os.path.join(figure_dir, 'dipoleParameters_perCell')
        fig.savefig(name + '.pdf', bbox_inches='tight')
        fig.savefig(name + '.png', dpi=300, bbox_inches='tight')
        print(f'Figure 1 saved to {name}.pdf/.png')
    return fig


# %%
# ---------------------------------------------------------------------------
# Figures 8-10: dipole length and orientation of every source->target connection
# ---------------------------------------------------------------------------
# Rows of these figures, as (row label, y-axis label): the signed dipole - the two
# parameters combined, which is what prepDipoles_normal multiplies the potential by -
# and the two factors it is made of.
_PER_CONNECTION_ROWS = [
    ('Signed dipole', 'length x orientation (mm)'),
    ('Length', 'Length (mm)'),
    ('Orientation', ''),   # the tick labels below already name the two directions
]


def _default_source_labels():
    """Source names from the parameter-file key paths ('S1 L4 PV', 'Thalamus ThalPOm').

    Used when no population labels are supplied, so the figure can be built without a
    SomatoModel. Only the *position* in this list matters to _build_source_styles (it
    picks the source area block), and DIPOLE_SOURCE_PATHS shares that order with
    get_population_labels() by construction.
    """
    return [' '.join(path) for path in DIPOLE_SOURCE_PATHS]


def plot_dipole_parameters_per_connection(json_path, area, figure_dir=None,
                                          pop_labels=None):
    """
    Dipole length and orientation of every source -> target connection of one area.

    The dipole model is defined per (target, source) pair: a target layer's dipole
    weights the synaptic potential arriving from each of the 33 source populations by
    that pair's length and orientation. This figure shows those parameters directly,
    laid out like Simulations/plotting_functions.plot_connectivity - one column per
    target layer, the 33 sources along x, bar colour = source area, dashed lines
    between the source blocks - with one row per quantity: the signed dipole
    (length x orientation), the length in mm and the orientation (-1 / 0 / +1).

    Only the excitatory population of a layer is a dipole target, so here the columns
    are the target layers and the rows carry the quantity, where plot_connectivity has
    one panel per target population and rows of cell types. A3b is unlaminated, so its
    figure has a single column. The dipole has no background / external / modulatory
    sources - compute_dipoles drops those synapse columns - so the x axis is the 33
    populations only.

    Parameters
    ----------
    json_path : str
        Dipole parameter file, e.g. EEGSimulation/dipole_parameters_flippedPVSST.json
        (the one SomatoModel.load_dipole_params reads).
    area : str
        Key of _AREA_SPECS: 'A3b', 'A1' (= S1) or 'S2'.
    figure_dir : str, optional
        Where to save the figure. None (the default) only builds and returns it.
    pop_labels : sequence of str, optional
        Source population names in array order, e.g.
        SomatoModel.get_population_labels(). Falls back to the parameter file's key
        paths, so no model is needed.
    """
    spec = _AREA_SPECS[area]
    display = spec['display']
    layer_labels = spec['labels']

    dp = read_dipole_params(json_path)
    # A3b is stored as a single source vector, A1/S2 as one per target layer
    lengths = np.atleast_2d(np.asarray(dp['dipole_lengths'][area], dtype=float))
    orientations = np.atleast_2d(np.asarray(dp['dipole_orientation'][area], dtype=float))
    if lengths.shape[0] != len(layer_labels):
        raise ValueError(f'{area}: parameter file has {lengths.shape[0]} target layers, '
                         f'expected {len(layer_labels)}')

    if pop_labels is None:
        pop_labels = _default_source_labels()
    pop_labels = [str(label) for label in pop_labels]
    n_sources = lengths.shape[1]

    try:
        style_colors, _ = figure_style()
    except Exception:
        # figure_style() opens a tkinter window to query the screen size, which fails
        # on headless nodes - the figure should still be produced there
        style_colors = {}
    dark2 = sns.color_palette('Dark2')
    area_colors = {'A3b': dark2[2], 'S1': dark2[0], 'S2': dark2[1],
                   'Thalamus': style_colors.get('Thal', dark2[3])}
    source_styles = _build_source_styles(pop_labels, n_sources, area_colors)

    # one vector per (column, row) panel; the grid plotter only hands these keys back
    # to `weights`, so they can be anything hashable
    values = {}
    for c in range(len(layer_labels)):
        values[(c, 0)] = lengths[c] * orientations[c]
        values[(c, 1)] = lengths[c]
        values[(c, 2)] = orientations[c]

    nrow = len(_PER_CONNECTION_ROWS)
    # the target of a column is named once, in the title of its top panel; the rows are
    # named by the shared left-hand labels. A3b is unlaminated - its single column is
    # the area itself, so repeating the layer label there would read 'A3b A3b E'.
    targets = ([f'{display} {lab} E' for lab in layer_labels] if len(layer_labels) > 1
               else [f'{display} E'])
    columns = [(lab, [f'{target} (target)'] + [''] * (nrow - 1),
                [(c, r) for r in range(nrow)])
               for c, (lab, target) in enumerate(zip(layer_labels, targets))]

    savepath = None
    if figure_dir is not None:
        os.makedirs(figure_dir, exist_ok=True)
        savepath = os.path.join(figure_dir,
                                f'dipoleParameters_perConnection_{area}.png')

    fig = _plot_connectivity_grid(
        columns, values.__getitem__, [row for row, _ in _PER_CONNECTION_ROWS],
        source_styles, pop_labels,
        f'Dipole model per source -> target connection - {display}', savepath,
        xlabel='Source population', legend_title='Source area',
        ylabels=[ylabel for _, ylabel in _PER_CONNECTION_ROWS],
        col_width=4.6 if len(layer_labels) > 1 else 7.0)

    # the orientation is a sign, not a continuous value - label it as one. Only
    # possible once the axes exist, so the figure is written again below.
    ncol = len(layer_labels)
    for ax in fig.axes[(nrow - 1) * ncol:nrow * ncol]:
        ax.set_yticks([-1, 0, 1])
        ax.set_yticklabels(['inward\n(-1)', 'none\n(0)', 'outward\n(+1)'], fontsize=8)

    if figure_dir is not None:
        name = os.path.splitext(savepath)[0]
        fig.savefig(name + '.png', dpi=300, bbox_inches='tight')
        fig.savefig(name + '.pdf', bbox_inches='tight')
        print(f'{display} per-connection dipole parameters saved to {name}.pdf/.png')
    return fig


# %%
# ---------------------------------------------------------------------------
# Figures 2-4: dipole computation example per area (potentials -> dipoles -> sum)
# ---------------------------------------------------------------------------
# Which excitatory populations and simDipoles rows belong to each area. Row order
# of simDipoles is 0 = A3b, 1-4 = A1 L{1,4,5,6}_E, 5-8 = S2 L{1,4,5,6}_E, and the
# population names come from model.get_population_mapping(). A3b is a single
# (unlaminated) population, so it has one component instead of four.
_AREA_SPECS = {
    'A3b': {'display': 'A3b', 'pops': ['E'],
            'rows': slice(0, 1), 'labels': ['A3b']},
    'A1':  {'display': 'S1 (A1)', 'pops': ['L1_E', 'L4_E', 'L5_E', 'L6_E'],
            'rows': slice(1, 5), 'labels': ['L2/3', 'L4', 'L5', 'L6']},
    'S2':  {'display': 'S2', 'pops': ['L1_E', 'L4_E', 'L5_E', 'L6_E'],
            'rows': slice(5, 9), 'labels': ['L2/3', 'L4', 'L5', 'L6']},
}
# qualitative palette: maximally distinct, colourblind-friendly colour per layer
_LAYER_COLORS = ['#4477AA', '#EE6677', '#228833', '#AA3377']  # blue, red, green, purple


def plot_dipole_computation_area(model, simDipoles, area, figure_dir=None):
    """
    Demonstrate the dipole computation for one area ('A3b', 'A1' or 'S2').

    Top panel  : the excitatory membrane potential feeding each layer of the area
                 (input to the dipole projection).
    Bottom panel: the per-layer dipole trace (potential weighted by the layer's
                 dipole model) and the summed area dipole (bold black). A3b has a
                 single population, so its dipole is shown on its own.

    Parameters
    ----------
    model : SomatoModel
        A simulated model (model.potential populated).
    simDipoles : np.ndarray, shape (9, ntimes)
        Output of model.compute_dipoles(...).
    area : str
        Key of _AREA_SPECS: 'A3b', 'A1' (= S1) or 'S2'.
    figure_dir : str, optional
        Where to save the figure. None (the default) only builds and returns it, for
        callers that just want to look at the current run (see simulation_main.py).
    """
    spec = _AREA_SPECS[area]
    display = spec['display']
    labels = spec['labels']
    colors = _LAYER_COLORS[:len(labels)]

    step_size = model.step_size
    input_onset = model.input_onset
    Iext_dur = model.Iext_duration

    # excitatory populations of this area (see get_population_mapping)
    pop_map = model.get_population_mapping()[area]
    e_idx = [pop_map[p] for p in spec['pops']]

    area_dipoles = simDipoles[spec['rows']]
    area_sum = area_dipoles.sum(axis=0)
    laminated = len(labels) > 1   # A3b is a single population -> no per-layer split

    ntimes = simDipoles.shape[1]
    t = np.arange(ntimes) * step_size - input_onset  # time re-zeroed to stimulus onset
    window = (-0.05, 0.25)
    mask = (t >= window[0]) & (t <= window[1])

    fig, axes = plt.subplots(2, 1, figsize=(5, 5), sharex=True)

    # --- Top: excitatory membrane potentials feeding each layer ---
    ax = axes[0]
    ax.axvspan(0, Iext_dur, color='0.85', linewidth=0, zorder=0)
    for idx, lab, col in zip(e_idx, labels, colors):
        pot = model.potential[idx].sum(axis=0)  # summed synaptic (membrane) potential
        ax.plot(t[mask], pot[mask], color=col, linewidth=1.2, label=lab)
    ax.set_ylabel('E potential (mV)')
    ax.set_title(f'{display} excitatory membrane potential (input)')
    if laminated:
        ax.legend(title='Layer', loc='upper right', ncol=2, fontsize=7)

    # --- Bottom: per-layer dipole traces + summed area dipole ---
    ax = axes[1]
    ax.axvspan(0, Iext_dur, color='0.85', linewidth=0, zorder=0)
    if laminated:
        for row, lab, col in zip(area_dipoles, labels, colors):
            ax.plot(t[mask], row[mask], color=col, linewidth=1.2, label=lab)
        ax.plot(t[mask], area_sum[mask], color='black', linewidth=2.2,
                label=f'Sum ({area})', zorder=5)
        ax.set_title(f'{display} layer dipoles and their sum (output)')
        ax.legend(loc='upper right', ncol=2, fontsize=7)
    else:
        ax.plot(t[mask], area_sum[mask], color='black', linewidth=2.2,
                label=f'{area} dipole', zorder=5)
        ax.set_title(f'{display} dipole (output)')
        ax.legend(loc='upper right', fontsize=7)
    ax.set_ylabel('Dipole (a.u.)')
    ax.set_xlabel('Time from stimulus onset (s)')

    sns.despine(fig=fig, trim=True)
    fig.tight_layout()

    if figure_dir is not None:
        name = os.path.join(figure_dir, f'dipoleComputation_{area}_example')
        fig.savefig(name + '.pdf', bbox_inches='tight')
        fig.savefig(name + '.png', dpi=300, bbox_inches='tight')
        print(f'{display} dipole computation figure saved to {name}.pdf/.png')
    return fig


# %%
# ---------------------------------------------------------------------------
# Figures 5-7: how each layer's interneurons contribute to its dipole
# ---------------------------------------------------------------------------
# The dipole model assigns a length and an orientation to every *presynaptic*
# (input) population, so each layer's dipole is a weighted sum over the synaptic
# potentials the layer's E population receives (see compute_dipoles):
#
#   dipole_layer(t) = sum_j w_j * potential[L*_E, j, t]
#   w_j = length_j * orientation_j * resistance_factor * cellcount_rel
#
# With the orientations of the file the model loads (E -1, PV +1, SST -1, VIP 0 in
# dipole_parameters_flippedPVSST.json), PV and SST inhibition push the dipole in
# opposite directions, and VIP input is dipole-silent even though it is a real
# synaptic input. _celltype_orientations reads those signs from the file rather than
# repeating them here, since dipole_parameters.json uses the opposite PV/SST
# convention. The dipole length is
# layer-specific (L2/3 0.2, L4 0.1, L5 1.5, L6 1.0 mm), so how much each
# interneuron class matters differs from layer to layer.

def _source_population_table(model):
    """One record per source population, in W/potential order.

    Returns a list of (label, area, layer_index, cell type), where

      label       : canonical short name from model.get_population_labels()
                    ('E1', 'SST4S2', 'ThalPOm', ...),
      area        : source-area name of Simulations/plotting_functions.SOURCE_AREA_BLOCKS
                    ('A3b', 'S1', 'S2', 'Thalamus'), which drives the trace colour,
      layer_index : 0-3 for L2/3, L4, L5, L6 (None for A3b and the thalamus, which are
                    unlaminated here), which drives the lightness of that colour,
      cell type   : 'E', 'PV', 'SST', 'VIP' or 'I' (the reticular nucleus).

    The thalamic labels have no E/PV/SST/VIP prefix, so they are mapped explicitly:
    ThalE (VPM) and ThalPOm are excitatory, ThalI is the inhibitory reticular nucleus
    (which is dipole-silent, length 0).
    """
    labels = model.get_population_labels()
    thal_types = {'ThalE': 'E', 'ThalI': 'I', 'ThalPOm': 'E'}

    table = []
    for j in range(model.nPop):
        label = str(labels[j])
        area = next(name for name, block in SOURCE_AREA_BLOCKS if j in block)

        if area in ('S1', 'S2'):
            # names are <celltype><layer index>[S2], layer index 1 = L2/3 ... 4 = L6
            layer_index = int(next(c for c in label if c.isdigit())) - 1
        else:
            layer_index = None

        if label in thal_types:
            cell_type = thal_types[label]
        else:
            cell_type = next(c for c in ('SST', 'VIP', 'PV', 'E') if label.startswith(c))

        table.append((label, area, layer_index, cell_type))

    return table


_CELLTYPE_ROWS = ('E', 'PV', 'SST', 'VIP')
# colour per source area, as in Simulations/plotting_functions.plot_all_potentials
_dark2 = sns.color_palette('Dark2')
_SOURCE_AREA_COLORS = {'A3b': _dark2[2], 'S1': _dark2[0], 'S2': _dark2[1],
                       'Thalamus': colors.get('Thal', _dark2[3])}
# fraction blended towards white per source layer: L2/3 lightest ... L6 = the area colour
_LAYER_SHADES = [0.55, 0.40, 0.22, 0.0]


# Sources that must be darkened to stay distinguishable. The lightness of a trace
# normally encodes the source layer, but the thalamus is unlaminated, so ThalE (VPM) and
# ThalPOm both fall back to the plain Thalamus colour - and both are excitatory, so they
# share the E row and its legend. Darkening POm keeps the same hue while separating the
# two. Fraction blended towards black.
_SOURCE_DARKEN = {'ThalPOm': 0.45}


def _source_color(area, layer_index, label=None):
    """Area colour, lightened according to the source layer (None = unlaminated source).

    `label` is the source population's name; sources listed in _SOURCE_DARKEN are
    additionally darkened so that two unlaminated sources of the same area and cell type
    do not end up with an identical colour.
    """
    base = np.array(_SOURCE_AREA_COLORS[area], dtype=float)
    if layer_index is not None:
        base = base + (1.0 - base) * _LAYER_SHADES[layer_index]
    darken = _SOURCE_DARKEN.get(label, 0.0)
    if darken:
        base = base * (1.0 - darken)
    return tuple(base)


def _celltype_orientations(model, area, sources):
    """Dipole orientation per source cell type, read from the JSON the model loads.

    Never hard-code these: dipole_parameters_flippedPVSST.json (the file
    load_dipole_params reads) flips PV and SST relative to dipole_parameters.json.
    Returns {cell type: -1 / 0 / +1}, or None for a cell type whose sources do not
    share one sign in this area.
    """
    _, orientation, _, _ = model.load_dipole_params()
    rows = np.atleast_2d(np.asarray(orientation[area], dtype=float))
    signs = {}
    for ctype in _CELLTYPE_ROWS:
        values = {rows[r, j] for r in range(rows.shape[0])
                  for j, rec in enumerate(sources) if rec[3] == ctype}
        signs[ctype] = int(values.pop()) if len(values) == 1 else None
    return signs


def _layer_weight(cache, area, layer_index):
    """Mean-over-subjects dipole weight vector (33,) for one target population.

    compute_dipoles averages w_k . potential over subjects and the potential does
    not depend on subject, so averaging the weight vectors reproduces it exactly.
    A3b is stored as a single vector, A1/S2 as a list of four layer vectors.
    """
    return np.mean([(s[area] if layer_index is None else s[area][layer_index])
                    for s in cache['per_subject']], axis=0)


def plot_layer_interneuron_contributions(model, simDipoles, area, subjects,
                                         figure_dir=None):
    """
    Show how every source population contributes to each layer's dipole.

    One column per layer (A3b has a single unlaminated population, so one column)
    and one row per source cell type, plus a summary row:

    Rows E / PV / SST / VIP: the individual contributions w_j * potential[L*_E, j]
           of every source population of that cell type, with the cell type's sum
           in bold black. Colour encodes the source area and its lightness the
           source layer, line style the source cell type - the same grammar as
           plot_all_potentials in Simulations/plotting_functions.py.
           PV and SST have opposite orientations, so their contributions push the
           dipole in opposite directions; VIP has dipole length 0, so its row is a
           flat zero even though VIP input is a real synaptic input.
    Summary row: the four cell-type sums, summing exactly to that layer's dipole
           (bold black).

    Sources that contribute exactly zero (no anatomical connection) are left out of
    the traces but still counted in the sums; each panel's legend title reports how
    many of the cell type's sources are drawn. Dipole-silent cell types keep their
    zero traces, since there the zero is the result.

    Parameters
    ----------
    model : SomatoModel
        A simulated model (model.potential populated).
    simDipoles : np.ndarray, shape (9, ntimes)
        Output of model.compute_dipoles(subjects).
    area : str
        'A3b', 'A1' (= S1) or 'S2'.
    subjects : list of int
        The subject list passed to compute_dipoles, used to look up the cached
        dipole projection vectors.
    figure_dir : str, optional
        Where to save the figure. None (the default) only builds and returns it.
    """
    spec = _AREA_SPECS[area]
    display = spec['display']
    pop_map = model.get_population_mapping()[area]
    laminated = len(spec['labels']) > 1

    # one entry per column: (layer title, target E index, weight-vector index,
    #                        simDipoles row)
    if laminated:
        columns = [(lab, pop_map[f'{lay}_E'], i, spec['rows'].start + i)
                   for i, (lay, lab) in enumerate(zip(['L1', 'L4', 'L5', 'L6'],
                                                      spec['labels']))]
    else:
        columns = [(spec['labels'][0], pop_map['E'], None, spec['rows'].start)]

    # dipole projections built by compute_dipoles; rebuilt only if this figure is
    # made standalone (reading every subject's forward solution is slow).
    cache = model._dipole_projection_cache.get(tuple(subjects))
    if cache is None:
        cache = model._build_dipole_projections(subjects)

    sources = _source_population_table(model)
    orientations = _celltype_orientations(model, area, sources)

    def _orient_label(name, ctype):
        sign = orientations[ctype]
        if sign is None:
            return name          # not one sign for this cell type in this area
        return f'{name} ({sign:+d})' if sign else f'{name} (0)'

    ntimes = simDipoles.shape[1]
    t = np.arange(ntimes) * model.step_size - model.input_onset
    window = (-0.05, 0.25)
    mask = (t >= window[0]) & (t <= window[1])

    ncol = len(columns)
    nrow = len(_CELLTYPE_ROWS) + 1
    fig, axes = plt.subplots(nrow, ncol, sharex=True, squeeze=False,
                             figsize=(3.6 * ncol if laminated else 5.0, 2.3 * nrow))

    # legend entries per cell-type row, collected over all columns: colour + shade
    # identify the (source area, source layer) pair uniquely within a row, so one
    # legend outside the last column labels every trace of that row.
    row_legend = {ctype: {} for ctype in _CELLTYPE_ROWS}

    for c, (lab, e_idx, li, dipole_row) in enumerate(columns):
        w = _layer_weight(cache, area, li)
        # synaptic potentials arriving at this layer's E population, one row per
        # source. Only the first nPop columns are population-to-population synapses;
        # the remaining ones are the external drives (background, external and, since
        # the mechanical input was added, a third), which compute_dipoles excludes from
        # the dipole. Indexing by nPop rather than counting input columns from the end
        # keeps this correct however many drives the model grows - a hard-coded [:-2]
        # silently became wrong when the third one appeared.
        pot = model.potential[e_idx, :model.nPop]
        contrib = w[:, None] * pot               # (33, ntimes) signed contributions

        groups = {}
        for j, (_, _, _, ctype) in enumerate(sources):
            groups.setdefault(ctype, np.zeros(pot.shape[1]))
            groups[ctype] += contrib[j]
        # 'I' (the reticular nucleus) is the only cell type outside _CELLTYPE_ROWS;
        # it has dipole length 0, so the four plotted classes are complete. The
        # total sums over every class regardless, so the assertion stays exact.
        total = sum(groups.values())
        # relative tolerance, not a fixed 1e-8: the dipole amplitude depends entirely on
        # the parameter set (order 1 for the optimized run, order 100 for some of the
        # high-gain sets simulation_main.py explores), and summing 33 float contributions
        # cannot hold an absolute bound that tight at the top of that range.
        assert np.allclose(total, simDipoles[dipole_row]), \
            f'cell-type decomposition does not reproduce the {area} {lab} dipole'

        # --- one row per source cell type: the individual source contributions ---
        for r, ctype in enumerate(_CELLTYPE_ROWS):
            ax = axes[r, c]
            ax.axvspan(0, model.Iext_duration, color='0.85', linewidth=0, zorder=0)

            members = [(j, rec) for j, rec in enumerate(sources) if rec[3] == ctype]
            # a dipole-silent cell type (dipole length or orientation 0, e.g. VIP)
            # keeps its zero traces: there the zero is the result, not a missing
            # connection. Decided from the contributions themselves rather than from
            # the orientation, so it holds whichever parameter file is loaded.
            keep_zero = all(np.allclose(contrib[j], 0) for j, _ in members)
            drawn = 0
            for j, (label, src_area, layer_index, _) in members:
                if not keep_zero and np.allclose(contrib[j], 0):
                    continue
                color = _source_color(src_area, layer_index, label)
                ax.plot(t[mask], contrib[j][mask], color=color,
                        linestyle=CELLTYPE_LINESTYLES[ctype], linewidth=1.1)
                row_legend[ctype].setdefault(label, color)
                drawn += 1
            ax.plot(t[mask], groups[ctype][mask], color='black', linewidth=2.0,
                    zorder=5)
            ax.axhline(0, color='0.6', linewidth=0.6, zorder=1)
            # how many of this cell type's sources reach this layer at all; as a
            # right-aligned title so it cannot collide with the traces
            ax.set_title(f'{drawn}/{len(members)} sources', loc='right', fontsize=6,
                         color='0.4')
            if keep_zero:
                ax.text(0.5, 0.75, 'dipole-silent', transform=ax.transAxes,
                        ha='center', va='center', fontsize=7, color='0.4')

            if r == 0 and laminated:
                ax.set_title(lab)     # A3b has one column, named by the suptitle
            if c == 0:
                ax.set_ylabel('Contribution (a.u.)')
                ax.text(-0.42, 0.5, _orient_label(f'{ctype} inputs', ctype),
                        transform=ax.transAxes, rotation=90, va='center',
                        ha='center', fontweight='bold')

        # --- summary row: the cell-type sums and the resulting dipole ---
        ax = axes[-1, c]
        ax.axvspan(0, model.Iext_duration, color='0.85', linewidth=0, zorder=0)
        for ctype in _CELLTYPE_ROWS:
            ax.plot(t[mask], groups[ctype][mask], color=_CELL_COLORS[ctype],
                    linewidth=1.4, label=_orient_label(ctype, ctype))
        if not np.allclose(groups.get('I', 0), 0):
            ax.plot(t[mask], groups['I'][mask], color='0.4', linewidth=1.4,
                    label='ThalI')
        ax.plot(t[mask], total[mask], color='black', linewidth=2.2,
                label=f'Total {"layer" if laminated else "A3b"} dipole', zorder=5)
        ax.axhline(0, color='0.6', linewidth=0.6, zorder=1)
        ax.set_xlabel('Time from stimulus onset (s)')
        if c == 0:
            ax.set_ylabel('Dipole (a.u.)')
            ax.text(-0.42, 0.5, 'All inputs', transform=ax.transAxes, rotation=90,
                    va='center', ha='center', fontweight='bold')
        if c == ncol - 1:
            ax.legend(title='Sum by input cell type', frameon=False, fontsize=6,
                      title_fontsize=6, loc='center left', bbox_to_anchor=(1.02, 0.5))

    # one legend per cell-type row, outside the last column: within a row the
    # colour (source area) and its shade (source layer) name the source uniquely.
    for r, ctype in enumerate(_CELLTYPE_ROWS):
        entries = row_legend[ctype]
        handles = [plt.Line2D([0], [0], color=col, linewidth=1.4,
                              linestyle=CELLTYPE_LINESTYLES[ctype])
                   for col in entries.values()]
        handles.append(plt.Line2D([0], [0], color='black', linewidth=2.0))
        axes[r, -1].legend(handles, list(entries) + [f'{ctype} sum'],
                           title=f'{ctype} sources', frameon=False, fontsize=6,
                           title_fontsize=6, ncol=1 if len(handles) <= 6 else 2,
                           loc='center left', bbox_to_anchor=(1.02, 0.5))

    fig.suptitle(f'{display}: source contributions to the layer dipoles'
                 if laminated else f'{display}: source contributions to the dipole')
    sns.despine(fig=fig, trim=True)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    if figure_dir is not None:
        name = os.path.join(figure_dir,
                            f'dipoleComputation_layerInterneurons_{area}_example')
        fig.savefig(name + '.pdf', bbox_inches='tight')
        fig.savefig(name + '.png', dpi=300, bbox_inches='tight')
        print(f'{display} layer interneuron figure saved to {name}.pdf/.png')
    return fig


# %%
# ---------------------------------------------------------------------------
# Standalone demo: run one simulation of its own and build all seven figures.
#
# Guarded so the functions above can be imported and applied to someone else's
# model - Simulations/simulation_main.py plots the dipole computation of the run it
# has just simulated this way. Without the guard, importing this module would run
# the extra simulation below and write its figures as a side effect.
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    # output directory
    SIMDIR = os.getenv("SIMDIR")
    figure_dir = os.path.join(SIMDIR, "Figures", "dipole_computation")
    os.makedirs(figure_dir, exist_ok=True)

    # subject forward model(s) used to project the dipoles (as in simulation_main.py)
    subjects = [15]

    # 1) Run one stimulated simulation.
    # Read the parameters from an optimization run. Same configuration as
    # Simulations/simulation_main.py: the simulation_parameter.json base params
    # updated with the run's best_params, background noise switched off.
    # Keep opt_run in sync with simulation_main.py's - the two are set independently.
    opt_run = "opt_20260806_141628_tc_roi-S2" #"opt_20260729_093613_tc_roi-S2" #"opt_20260729_114525_tc_roi-A1"
    params = load_optimized_params(opt_run, overrides={'Ib_noise_std': 0})

    model = SomatoModel(params)
    model.simulate()

    # 2) Compute the dipole from the simulated potentials.
    #    simDipoles shape (9, ntimes): 0 = A3b, 1-4 = A1 L{1,4,5,6}_E, 5-8 = S2 L{1,4,5,6}_E
    simDipoles = model.compute_dipoles(subjects)

    # Figure 1: dipole parameters per cell
    fig1 = plot_dipole_parameters_clean(
        os.path.join(_eeg_dir, 'dipole_parameters_flippedPVSST.json'), figure_dir)

    # Figures 2-4: dipole computation example per area
    figs = {area: plot_dipole_computation_area(model, simDipoles, area, figure_dir)
            for area in ('S2', 'A1', 'A3b')}

    # Figures 5-7: source contributions per layer
    figs_IN = {area: plot_layer_interneuron_contributions(
                   model, simDipoles, area, subjects, figure_dir)
               for area in ('S2', 'A1', 'A3b')}

    # Figures 8-10: dipole parameters per source->target connection
    figs_conn = {area: plot_dipole_parameters_per_connection(
                     os.path.join(_eeg_dir, 'dipole_parameters_flippedPVSST.json'),
                     area, figure_dir, model.get_population_labels())
                 for area in ('S2', 'A1', 'A3b')}

    plt.show()
