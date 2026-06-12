"""
Visualize dipole_parameters.json: dipole lengths and orientations per population.
"""
# %%
import json
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

_eeg_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.dirname(_eeg_dir)
for _p in [_eeg_dir,
           os.path.join(_project_root, 'Simulations'),
           os.path.join(_project_root, 'Simulations', 'model')]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

from sim_eeg import get_population_mapping


_CELL_COLORS = {
    'E':   '#4477AA',
    'PV':  '#EE6677',
    'SST': '#228833',
    'VIP': '#CCBB44',
    'I':   '#AA3377',
}

_SOURCE_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
_SOURCE_LABELS  = ['Source L1', 'Source L2', 'Source L3', 'Source L4']

# Background shade per region block on the x-axis
_REGION_SPANS = [
    ('A3b',      0,  4, '#dce8ff'),
    ('S1',       4, 17, '#ffe0e0'),
    ('S2',      17, 30, '#e0ffe0'),
    ('Thalamus',30, 32, '#fff8d0'),
]

# For each display area: (JSON key, is_multilayer)
_AREA_META = {
    'A3b':      ('A3b', False),
    'S1':       ('A1',  True),
    'S2':       ('S2',  True),
    'Thalamus': (None,  False),
}

N_POPS = 32


def load_dipole_json(json_path=None):
    if json_path is None:
        json_path = os.path.join(_eeg_dir, 'dipole_parameters.json')
    with open(json_path, 'r') as f:
        return json.load(f)


def _build_pop_labels(pop_mapping):
    """Return list of (x_label, cell_type) ordered by population index 0–31."""
    short = {'A3b': 'A3b', 'S1': 'S1', 'S2': 'S2', 'Thalamus': 'Thal'}
    idx_info = {}
    for region, pops in pop_mapping.items():
        for name, idx in pops.items():
            cell_type = name.split('_')[-1]
            idx_info[idx] = (f"{short[region]}-{name}", cell_type)
    return [idx_info[i] for i in range(N_POPS)]


def _signed_series(l_arr, o_arr, n=N_POPS):
    return [float(l_arr[i]) * int(o_arr[i]) if i < len(l_arr) and i < len(o_arr) else 0.0
            for i in range(n)]


def plot_dipole_parameters(json_path=None):
    """
    One subplot per area; every subplot shows all 32 populations on the x-axis.

    The y-axis is signed dipole length (length × orientation in mm):
      positive  → outward current (+1)
      negative  → inward current  (−1)
      zero      → no dipole        (0)

    Single-array areas (A3b, Thalamus): bar chart coloured by cell type.
    Multi-array areas (S1, S2): one coloured line per source layer.
    Background shading indicates which populations belong to which region.

    Args:
        json_path: path to dipole_parameters.json (default: same directory)

    Returns:
        matplotlib Figure
    """
    dipole_params = load_dipole_json(json_path)
    pop_mapping   = get_population_mapping()

    pop_info   = _build_pop_labels(pop_mapping)
    x_labels   = [p[0] for p in pop_info]
    cell_types = [p[1] for p in pop_info]
    bar_colors = [_CELL_COLORS.get(ct, '#888888') for ct in cell_types]
    x = np.arange(N_POPS)

    regions = ['A3b', 'S1', 'S2', 'Thalamus']
    fig, axes = plt.subplots(1, 4, figsize=(22, 6), sharey=False)
    fig.suptitle('Dipole Parameters per Population', fontsize=13, fontweight='bold', y=1.02)

    for ax, region in zip(axes, regions):
        json_key, multilayer = _AREA_META[region]

        # Shaded background per region block
        for _, start, end, color in _REGION_SPANS:
            ax.axvspan(start - 0.5, end - 0.5, alpha=0.35, color=color, zorder=0)

        # Thin vertical separators between region blocks
        for sep in [3.5, 16.5, 29.5]:
            ax.axvline(sep, color='#888888', linestyle=':', linewidth=0.8, zorder=1)

        if json_key is None:
            # Thalamus not in JSON → all-zero bars
            ax.bar(x, np.zeros(N_POPS), color=bar_colors,
                   edgecolor='black', linewidth=0.3, zorder=2)

        elif multilayer:
            # 4 source sub-arrays → one coloured line each
            l_raw = dipole_params['dipole_lengths'][json_key]
            o_raw = dipole_params['dipole_orientation'][json_key]
            for k in range(len(l_raw)):
                vals = _signed_series(l_raw[k], o_raw[k])
                ax.plot(x, vals, 'o-',
                        color=_SOURCE_COLORS[k], label=_SOURCE_LABELS[k],
                        linewidth=1.5, markersize=3.5, alpha=0.85, zorder=2)
            ax.legend(fontsize=7, loc='upper left', framealpha=0.85)

        else:
            # Single array → bar chart coloured by cell type
            l_arr = dipole_params['dipole_lengths'][json_key]
            o_arr = dipole_params['dipole_orientation'][json_key]
            vals  = _signed_series(l_arr, o_arr)
            ax.bar(x, vals, color=bar_colors,
                   edgecolor='black', linewidth=0.3, zorder=2)

        ax.axhline(0, color='black', linewidth=0.8, zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels(x_labels, rotation=90, fontsize=5.5)
        ax.set_title(region, fontsize=11, fontweight='bold')
        ax.set_xlabel('Population (index order)', fontsize=8)
        ax.set_ylabel('Signed dipole length (mm)', fontsize=8)
        ax.grid(axis='y', linestyle='--', alpha=0.5, linewidth=0.4, zorder=1)
        sns.despine(ax=ax, trim=True)

    # Bottom legend: cell-type colours + region background
    cell_patches = [
        mpatches.Patch(facecolor=c, label=t, edgecolor='black', linewidth=0.5)
        for t, c in _CELL_COLORS.items()
    ]
    region_patches = [
        mpatches.Patch(facecolor=color, alpha=0.6, label=name, edgecolor='#888888', linewidth=0.5)
        for name, _, _, color in _REGION_SPANS
    ]
    fig.legend(
        handles=cell_patches + region_patches,
        title='Cell type  |  Region background',
        loc='lower center',
        ncol=len(cell_patches) + len(region_patches),
        bbox_to_anchor=(0.5, -0.10),
        frameon=True,
        fontsize=8,
    )
    fig.text(
        0.5, -0.04,
        'Positive = outward orientation (+1)   |   Negative = inward orientation (−1)   |   Zero = no dipole (0)',
        ha='center', fontsize=8, color='#444444',
    )

    plt.tight_layout()
    return fig


if __name__ == '__main__':
    fig = plot_dipole_parameters()
    plt.show()
