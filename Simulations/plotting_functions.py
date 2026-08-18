import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from plotting_style import figure_style

# -----------------------------------------------------------------------------
# Population layout
# Index order mirrors SomatoModel.get_population_labels() (model/somato_model.py).
# Grids are laid out as rows = layers, columns = CELLTYPE_COLUMNS;
# None marks a cell type that does not exist in that layer (VIP only in L2/3).
# -----------------------------------------------------------------------------
CELLTYPE_COLUMNS = ['E', 'PV', 'SST', 'VIP']
LAYER_ROWS = ['Layer 2/3', 'Layer 4', 'Layer 5', 'Layer 6']

POP_GRID_A3B = [[0, 1, 2, 3]]
POP_GRID_S1 = [[4, 5, 6, 7], [8, 9, 10, None], [11, 12, 13, None], [14, 15, 16, None]]
POP_GRID_S2 = [[17, 18, 19, 20], [21, 22, 23, None], [24, 25, 26, None], [27, 28, 29, None]]
POP_GRID_THAL = [[30, 31, 32, None]]  # ThalE (VPM), ThalI (reticular), ThalPOm
THAL_COLUMNS = ['ThalE (VPM)', 'ThalI (reticular)', 'ThalPOm', '']

# presynaptic source blocks -> line colour (the two extra synapse columns of the
# 3D potential array, background and external input, are handled separately)
SOURCE_AREA_BLOCKS = [('A3b', range(0, 4)), ('S1', range(4, 17)),
                      ('S2', range(17, 30)), ('Thalamus', range(30, 33))]
# source cell type -> line style; thalamic sources reuse E/PV/SST styles in
# array order (ThalE, ThalI, ThalPOm)
CELLTYPE_LINESTYLES = {'E': '-', 'PV': '--', 'SST': '-.', 'VIP': ':'}


def plot_minmax(rates, coupling_strengths_Es):
    minRate = np.min(rates[:,:,-100:],axis=2)
    maxRate = np.max(rates[:,:,-100:],axis=2)

    # Plottign Area 3b Activity
    fig, axs = plt.subplots(1, 1, figsize=(3, 6)) 
    axs.plot(coupling_strengths_Es, minRate[:,0:3], linewidth=2)
    axs.plot(coupling_strengths_Es, maxRate[:,0:3], linewidth=0.5)
    axs.grid(True)
    axs.set_ylabel('Hz')
    axs.legend(['E', 'PV', 'SOM'])
    plt.tight_layout() 
    plt.legend()

    # plot results for S1
    fig, axs = plt.subplots(4, 1, figsize=(3, 6))  # Set figure size
    # Plot settings for all subplots
    for i, ax in enumerate(axs, start=1):
        ax.plot(coupling_strengths_Es, minRate[:,((i-1)*4)+3:i*4+3], linewidth=2)
        ax.plot(coupling_strengths_Es, maxRate[:,((i-1)*4)+3:i*4]+3, linewidth=0.5)
        ax.grid(True)
        ax.set_ylabel('Hz')
        ax.legend(['L2/3', 'L4', 'L5', 'L6'])

    # Set titles for each subplot
    axs[0].set_title('E')
    axs[1].set_title('PV')
    axs[2].set_title('SOM')
    axs[3].set_title('VIP')
    #sns.despine(trim=True, bottom=True)
    plt.tight_layout() 
    plt.legend()

    # plot results for the S2 column
    figS2, axsS2 = plt.subplots(4, 1, figsize=(3, 6))  # Set figure size

    # Plot settings for all subplots
    for i, ax in enumerate(axsS2, start=14):
        ax.plot(coupling_strengths_Es, minRate[:,((i-1)*4)+3:i*4+3], linewidth=2)
        ax.plot(coupling_strengths_Es, maxRate[:,((i-1)*4)+3:i*4+3], linewidth=0.5)
        ax.grid(True)
        ax.set_ylabel('Hz')
        ax.legend(['L2/3', 'L4', 'L5', 'L6'])

    # Set titles for each subplot
    axs[0].set_title('E')
    axs[1].set_title('PV')
    axs[2].set_title('SOM')
    axs[3].set_title('VIP')
    #sns.despine(trim=True, bottom=True)
    plt.tight_layout() 
    plt.legend()
    plt.show()

def _source_celltype(label):
    """Cell type of a population label ('SST2S2' -> 'SST'). Order matters: SST/PV before E."""
    for celltype in ('SST', 'VIP', 'PV', 'E'):
        if label.startswith(celltype):
            return celltype
    return 'E'


def _build_source_styles(pop_labels, n_sources, area_colors):
    """
    Line style per synapse column of the 3D potential array.

    Colour encodes the source area, line style the source cell type. The last two
    columns are the background (-2) and external (-1) input synapses.

    Returns a list of (colour, linestyle, group), one entry per source column, where
    group is the source area name (or 'Background input' / 'External input') and is
    used to build the shared colour legend.
    """
    n_pop = len(pop_labels)
    styles = []
    for src in range(n_sources):
        if src >= n_pop:
            if src == n_sources - 1:
                styles.append(('black', '-', 'External input'))
            else:
                styles.append(([0.55, 0.55, 0.55], '-', 'Background input'))
            continue

        area, offset = None, 0
        for name, block in SOURCE_AREA_BLOCKS:
            if src in block:
                area, offset = name, src - block.start
                break

        if area == 'Thalamus':
            # ThalE / ThalI / ThalPOm have no E/PV/SST/VIP name, style them by array order
            linestyle = ['-', '--', '-.'][offset]
        else:
            linestyle = CELLTYPE_LINESTYLES[_source_celltype(str(pop_labels[src]))]
        styles.append((area_colors[area], linestyle, area))

    return styles


def _legend_area_entries(source_styles):
    """Unique (colour, group) entries for the source-area legend, in array order."""
    entries, seen = [], set()
    for (color, _, group) in source_styles:
        if group is not None and group not in seen:
            seen.add(group)
            entries.append((color, group))
    return entries


def _grid_rows(pop_grid, row_labels, pop_labels, panel_titles=None):
    """Turn a population grid into (row_label, panel_titles, population_indices) rows."""
    rows = []
    for r, pop_idxs in enumerate(pop_grid):
        titles = [None if idx is None else
                  (panel_titles[c] if panel_titles is not None else str(pop_labels[idx]))
                  for c, idx in enumerate(pop_idxs)]
        rows.append((row_labels[r], titles, pop_idxs))
    return rows


def _plot_potential_grid(rows, all_potentials, steps, start_plot, source_styles,
                         stim_window, suptitle, savepath, stop_plot=None):
    """
    Plot one figure of per-synapse potential decompositions.

    Parameters:
    rows : list of (str, list, list)
        One entry per grid row: row label, panel titles, population indices
        (None leaves the cell empty).
    all_potentials : np.ndarray
        3D potential array (target population, source synapse, time).
    steps : np.ndarray
        Time axis in ms, same length as the time axis of all_potentials.
    start_plot : int
        First time sample to show.
    source_styles : list of (colour, linestyle, group)
        Style per source synapse, from _build_source_styles.
    stim_window : tuple or None
        (start, stop) in ms of the external input, shaded in grey.
    suptitle : str
    savepath : str
        Full path of the PDF to write.
    """
    nrows = len(rows)
    ncols = max(len(pop_idxs) for _, _, pop_idxs in rows)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6*ncols, 2.6*nrows),
                             sharex=True, sharey=False)
    axes = np.atleast_2d(axes)
    time_ms = steps[start_plot:stop_plot]

    for r, (row_label, titles, pop_idxs) in enumerate(rows):
        for c in range(ncols):
            ax = axes[r, c]
            idx = pop_idxs[c] if c < len(pop_idxs) else None
            if idx is None:
                ax.axis('off')
                continue

            if stim_window is not None:
                ax.axvspan(stim_window[0], stim_window[1], color=[0.9, 0.9, 0.9],
                           linewidth=0, zorder=0)
            for src, (color, linestyle, _) in enumerate(source_styles):
                ax.plot(time_ms, all_potentials[idx, src, start_plot:stop_plot], color=color,
                        linestyle=linestyle, linewidth=0.6)

            ax.set_title(titles[c], fontweight='bold')
            if c == 0:
                ax.set_ylabel('mV')
                if row_label:
                    ax.text(-0.32, 0.5, row_label, transform=ax.transAxes, rotation=90,
                            va='center', ha='center', fontweight='bold')

    # x-labels on the lowest *visible* axis of each column (VIP columns end early)
    for c in range(ncols):
        visible_rows = [r for r, (_, _, pop_idxs) in enumerate(rows)
                        if c < len(pop_idxs) and pop_idxs[c] is not None]
        if visible_rows:
            ax = axes[visible_rows[-1], c]
            ax.set_xlabel('Time (ms)')
            ax.tick_params(labelbottom=True)

    # shared legends: colour = source area, line style = source cell type
    area_entries = _legend_area_entries(source_styles)
    area_handles = [plt.Line2D([0], [0], color=color, linewidth=1.5)
                    for color, _ in area_entries]
    fig.legend(area_handles, [group for _, group in area_entries], title='Source area',
               frameon=False, loc='lower left', bbox_to_anchor=(1.0, 0.52))

    style_labels = ['E / ThalE', 'PV / ThalI', 'SST / ThalPOm', 'VIP']
    style_handles = [plt.Line2D([0], [0], color='0.3', linewidth=1.5,
                                linestyle=CELLTYPE_LINESTYLES[celltype])
                     for celltype in CELLTYPE_COLUMNS]
    fig.legend(style_handles, style_labels, title='Source cell type', frameon=False,
               loc='upper left', bbox_to_anchor=(1.0, 0.48))

    fig.suptitle(suptitle)
    sns.despine(trim=True)
    plt.tight_layout()
    plt.savefig(savepath, bbox_inches='tight', dpi=300)
    plt.show()


def plot_all_potentials(all_potentials, Iext, Ib, step_size, simulation_time, start_plot,
                        figdir, sI, g, d, sb, s, pop_labels=None, area='all', stop_plot=None):
    """
    Plot the per-synapse decomposition of every population potential.

    One figure per area (A3b + thalamus, S1, S2), laid out as rows = layers and
    columns = cell types. Each panel is one *target* population and shows the
    postsynaptic potential contributed by every *source* synapse: colour encodes the
    source area, line style the source cell type. Figures are saved as PDFs into
    figdir/single_simulations (nothing is shown interactively).

    Parameters:
    all_potentials : np.ndarray
        3D potential array of shape (nPop, nPop+2, n_timesteps) =
        (target population, source synapse, time). Source columns 0..nPop-1 are the
        PSPs from each presynaptic population, column -2 is the background input and
        column -1 the external input.
    Iext : np.ndarray
        external input array (used to shade the stimulus window)
    Ib : np.ndarray
        background input array (kept for signature consistency, not plotted)
    step_size : float
    simulation_time : float
    start_plot : int
        first time sample to plot (skips the initial transient)
    figdir : str
    sI : float
        EI balance
    g : float
        Coupling strength
    d : float
        External input duration
    sb : float
        Background input strength
    s : float
        External input strength
    pop_labels : sequence of str, optional
        Population names in array order, e.g. SomatoModel.get_population_labels().
        Falls back to generic 'pop{i}' names.
    area : str
        Area to plot ('all' or only 'A1')
    """
    all_potentials = np.asarray(all_potentials)
    if all_potentials.ndim != 3:
        raise ValueError(
            'plot_all_potentials needs the un-summed 3D potential array '
            '(target, source synapse, time), got shape '
            f'{all_potentials.shape}. With pyrates=True the potentials are already '
            'summed over synapses - use plot_potentials instead.')

    n_pop, n_sources, n_t = all_potentials.shape
    if pop_labels is None:
        pop_labels = [f'pop{i}' for i in range(n_pop)]
    steps = (np.arange(step_size, simulation_time+step_size, step_size)*1e3)[:n_t]

    try:
        colors, _ = figure_style()
    except Exception:
        # figure_style() opens a tkinter window to query the screen size, which
        # fails on headless nodes - the figures should still be produced there
        colors = {}

    dark2 = sns.color_palette('Dark2')
    area_colors = {'A3b': dark2[2], 'S1': dark2[0], 'S2': dark2[1],
                   'Thalamus': colors.get('Thal', dark2[3])}
    source_styles = _build_source_styles(pop_labels, n_sources, area_colors)

    # shade the external input window
    stim_window = None
    Iext = np.asarray(Iext)
    if Iext.ndim == 1 and np.any(Iext[:n_t] != 0):
        stim_samples = np.flatnonzero(Iext[:n_t] != 0)
        stim_window = (steps[stim_samples[0]], steps[stim_samples[-1]])

    figdir = os.path.join(figdir, 'single_simulations')
    os.makedirs(figdir, exist_ok=True)
    suffix = f'bEI-{sI}_g-{g}_Iextdur-{d}_Iextstr-{s}_Ibstr-{sb}'

    figures = []
    if area == 'all':
        figures.append(('A3b_thalamus',
                        'Potential decomposition - area 3b and thalamus',
                        _grid_rows(POP_GRID_A3B, ['A3b Layer 2/3'], pop_labels)
                        + _grid_rows(POP_GRID_THAL, ['Thalamus'], pop_labels,
                                     panel_titles=THAL_COLUMNS)))
    if area in ('all', 'A1'):
        figures.append(('S1', 'Potential decomposition - area 1 (S1)',
                        _grid_rows(POP_GRID_S1, LAYER_ROWS, pop_labels)))
    if area == 'all':
        figures.append(('S2', 'Potential decomposition - area S2',
                        _grid_rows(POP_GRID_S2, LAYER_ROWS, pop_labels)))

    for name, title, rows in figures:
        savepath = os.path.join(figdir, f'all_potentials_{name}_{suffix}.pdf')
        _plot_potential_grid(rows, all_potentials, steps, start_plot, source_styles,
                             stim_window, title, savepath, stop_plot)


def _conn_grid_columns(pop_grid, col_labels, pop_labels, panel_titles=None):
    """
    Turn a layer-major population grid into figure *columns*.

    The POP_GRID_* constants are laid out as rows = layers, columns = cell types;
    the connectivity figures use the flipped orientation (columns = layers, rows =
    cell types), so every grid row becomes one figure column.

    Returns a list of (column_label, panel_titles, population_indices), one entry per
    column, where the lists are indexed by figure row (cell type).
    """
    columns = []
    for c, pop_idxs in enumerate(pop_grid):
        titles = [None if idx is None else
                  (panel_titles[r] if panel_titles is not None else str(pop_labels[idx]))
                  for r, idx in enumerate(pop_idxs)]
        columns.append((col_labels[c], titles, list(pop_idxs)))
    return columns


def _source_block_bounds(source_styles):
    """
    Boundaries between consecutive source groups.

    Returns (separators, blocks) where separators are the x positions to draw a
    dividing line at and blocks is a list of (centre, group_name) for labelling.
    """
    separators, blocks = [], []
    start = 0
    for i in range(1, len(source_styles) + 1):
        if i == len(source_styles) or source_styles[i][2] != source_styles[start][2]:
            blocks.append(((start + i - 1) / 2, source_styles[start][2]))
            if i < len(source_styles):
                separators.append(i - 0.5)
            start = i
    return separators, blocks


def _plot_connectivity_grid(columns, weights, row_labels, source_styles, xtick_labels,
                            suptitle, savepath, xlabel='Source population',
                            legend_title='Source area'):
    """
    Plot one figure of per-population connection weights.

    Parameters:
    columns : list of (str, list, list)
        One entry per figure column (layer): column label, panel titles, population
        indices per row (None leaves the cell empty).
    weights : callable
        Maps a population index to the 1D weight vector shown in its panel.
    row_labels : list of str
        Left-hand labels, one per figure row (cell types).
    source_styles : list of (colour, linestyle, group)
        Style per entry of the weight vector, from _build_source_styles. Only the
        colour (source area) and the group (for separators/legend) are used here.
    xtick_labels : sequence of str
        Name per entry of the weight vector.
    suptitle : str
    savepath : str
        Full path of the PNG to write.
    xlabel : str
    legend_title : str
    """
    ncols = len(columns)
    nrows = max(len(pop_idxs) for _, _, pop_idxs in columns)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6*ncols, 2.4*nrows),
                             sharex=True, sharey='row', squeeze=False)

    bar_colors = [color for color, _, _ in source_styles]
    x = np.arange(len(source_styles))
    separators, _ = _source_block_bounds(source_styles)

    for c, (col_label, titles, pop_idxs) in enumerate(columns):
        for r in range(nrows):
            ax = axes[r, c]
            idx = pop_idxs[r] if r < len(pop_idxs) else None
            if idx is None:
                ax.axis('off')
                continue

            ax.bar(x, weights(idx), color=bar_colors, width=0.8, linewidth=0)
            ax.axhline(0, color='0.3', linewidth=0.5)
            for sep in separators:
                ax.axvline(sep, color='0.85', linewidth=0.5, linestyle='--', zorder=0)

            ax.set_title(titles[r], fontweight='bold')
            if c == 0:
                ax.set_ylabel('weight')
                ax.text(-0.18, 0.5, row_labels[r], transform=ax.transAxes, rotation=90,
                        va='center', ha='center', fontweight='bold')

    # column titles (layers) above the panel titles of the top row
    for c, (col_label, _, _) in enumerate(columns):
        axes[0, c].text(0.5, 1.16, col_label, transform=axes[0, c].transAxes,
                        ha='center', va='bottom', fontsize=13, fontweight='bold')

    # x tick labels on the lowest *visible* axis of each column (VIP rows end early).
    # sharex only propagates the tick *locations*, so rotation/size has to be set on
    # every axis that actually shows the labels.
    for c in range(ncols):
        pop_idxs = columns[c][2]
        visible_rows = [r for r in range(len(pop_idxs)) if pop_idxs[r] is not None]
        if not visible_rows:
            continue
        ax = axes[visible_rows[-1], c]
        ax.set_xticks(x)
        ax.set_xticklabels(xtick_labels, rotation=90, fontsize=5)
        ax.set_xlabel(xlabel)
        ax.tick_params(labelbottom=True)

    # shared legend: bar colour = area of the other side of the connection
    area_entries = _legend_area_entries(source_styles)
    area_handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color, _ in area_entries]
    fig.legend(area_handles, [group for _, group in area_entries], title=legend_title,
               frameon=False, loc='center left', bbox_to_anchor=(1.0, 0.5))

    fig.suptitle(suptitle)
    sns.despine(trim=False)
    plt.tight_layout()
    plt.savefig(savepath, bbox_inches='tight', dpi=300)
    plt.show()


def plot_connectivity(W, figdir, pop_labels=None, area='all', direction='in', suffix=''):
    """
    Plot the connectivity matrix as one figure per area.

    Each figure is a grid of bar panels laid out as columns = layers and rows = cell
    types (E, PV, SST, VIP); one panel per population. With direction='in' the panel
    shows the weights this population *receives* from every source (including the
    background and external input synapses), with direction='out' the weights it
    *sends* to every population. Bar colour encodes the area of the other side of the
    connection, the sign of the bar whether the connection is excitatory or
    inhibitory. The y-axis is shared within a row so cell types are comparable across
    layers without the large excitatory weights flattening the interneuron rows.

    Parameters:
    W : np.ndarray
        Connectivity matrix of shape (nPop, nPop+2) as returned by
        parameters.get_connectivity() / SomatoModel.W. Columns 0..nPop-1 are the
        presynaptic populations, column -2 the background and column -1 the external
        input.
    figdir : str
        Directory the PNGs are written to (created if missing).
    pop_labels : sequence of str, optional
        Population names in array order, e.g. SomatoModel.get_population_labels().
        Falls back to generic 'pop{i}' names.
    area : str
        Area the model was run with; only the figures of the simulated areas are
        produced ('all', 'A1', ...).
    direction : {'in', 'out'}
        Whether a panel shows the incoming or the outgoing weights of its population.
    suffix : str
        Appended to the file names, e.g. to tag a parameter set.
    """
    W = np.asarray(W)
    if W.ndim != 2:
        raise ValueError(f'plot_connectivity needs the 2D connectivity matrix, got shape {W.shape}.')
    if direction not in ('in', 'out'):
        raise ValueError(f"direction must be 'in' or 'out', got {direction!r}.")

    n_pop = W.shape[0]
    if pop_labels is None:
        pop_labels = [f'pop{i}' for i in range(n_pop)]
    pop_labels = [str(label) for label in pop_labels]

    try:
        colors, _ = figure_style()
    except Exception:
        # figure_style() opens a tkinter window to query the screen size, which
        # fails on headless nodes - the figures should still be produced there
        colors = {}

    dark2 = sns.color_palette('Dark2')
    area_colors = {'A3b': dark2[2], 'S1': dark2[0], 'S2': dark2[1],
                   'Thalamus': colors.get('Thal', dark2[3])}

    if direction == 'in':
        # full row of W: all presynaptic populations plus background and external input
        n_entries = W.shape[1]
        xtick_labels = pop_labels + ['Background', 'External']
        weights = lambda idx: W[idx, :]
        xlabel = 'Source population'
        legend_title = 'Source area'
        title_word = 'Incoming'
    else:
        # column of W: only the populations are valid targets (no input synapse rows)
        n_entries = n_pop
        xtick_labels = pop_labels
        weights = lambda idx: W[:, idx]
        xlabel = 'Target population'
        legend_title = 'Target area'
        title_word = 'Outgoing'

    source_styles = _build_source_styles(pop_labels, n_entries, area_colors)

    os.makedirs(figdir, exist_ok=True)

    figures = []
    if area == 'all':
        figures.append(('A3b_thalamus', 'area 3b and thalamus',
                        _conn_grid_columns(POP_GRID_A3B, ['A3b Layer 2/3'], pop_labels)
                        + _conn_grid_columns(POP_GRID_THAL, ['Thalamus'], pop_labels,
                                             panel_titles=THAL_COLUMNS)))
    if area in ('all', 'A1'):
        figures.append(('S1', 'area 1 (S1)',
                        _conn_grid_columns(POP_GRID_S1, LAYER_ROWS, pop_labels)))
    if area == 'all':
        figures.append(('S2', 'area S2',
                        _conn_grid_columns(POP_GRID_S2, LAYER_ROWS, pop_labels)))

    for name, title, columns in figures:
        savepath = os.path.join(figdir, f'connectivity_{direction}_{name}{suffix}.png')
        _plot_connectivity_grid(columns, weights, CELLTYPE_COLUMNS, source_styles,
                                xtick_labels, f'{title_word} connectivity - {title}',
                                savepath, xlabel=xlabel, legend_title=legend_title)


def plot_potentials(potentials, Iext, Ib, step_size, simulation_time, start_plot, figdir, sI, g, d, sb, s, area='all', stop_plot=None):
    """
    Plot population potentials for different areas and layers.
    Parameters:
    potentials : np.ndarray
    Iext : np.ndarray
        external input array
    Ib : np.ndarray
        background input array  
    step_size : float
    simulation_time : float
    start_plot : int
    figdir : str
    sI : float
        EI balance
    g : float
        Coupling strength
    d : float
        External input duration
    sb : float
        Background input strength
    s : float   
        External input strength 
    area : str
        Area to plot ('all' or only 'A1')
    """
    steps = np.arange(step_size, simulation_time+step_size, step_size)*1e3

    if area=='all':

        # Layout mirrors plot_results: 4 rows, 3 columns, with the input/thalamus/A3b
        # stack in column 0 and one row per cortical layer in columns 1 and 2.
        figure_style()
        fig, axes = plt.subplots(4, 3, figsize=(15, 15), sharex=True, sharey=False)
        axes = np.array(axes)

        # --- Column 0: external input, thalamus, POm, Area 3b ---
        axes[0, 0].plot(steps[start_plot:stop_plot], Iext[start_plot:stop_plot], label='Iext rate')
        axes[0, 0].plot(steps[start_plot:stop_plot], Ib[start_plot:stop_plot], label='Ib rate')
        axes[0, 0].legend(title='')
        axes[0, 0].set_title('External input')
        axes[0, 0].set_ylabel('Hz')

        # Thalamus block = last three populations, in array order: Thalamus E/VPM (-3),
        # Thalamus I/reticular (-2), POm (-1). Order matches get_population_labels().
        # POm gets its own panel: its potential is on a very different scale, so sharing
        # the axes with ThalE/ThalI flattens those two.
        plot_population_rates(axes[1, 0], [-3, -2], potentials, steps, start_plot,
                              ['Thalamus E', 'Thalamus I'], colors=['green', 'purple'],
                              stop_plot=stop_plot)
        axes[1, 0].set_title('Thalamus')
        axes[1, 0].set_ylabel('mV')

        plot_population_rates(axes[2, 0], [-1], potentials, steps, start_plot,
                              ['POm'], colors=['grey'], stop_plot=stop_plot)
        axes[2, 0].set_title('POm')
        axes[2, 0].set_ylabel('mV')

        plot_population_rates(axes[3, 0], POP_GRID_A3B[0], potentials, steps, start_plot,
                              CELLTYPE_COLUMNS, stop_plot=stop_plot)
        axes[3, 0].set_title('Area 3b')
        axes[3, 0].set_ylabel('mV')

        # --- Columns 1 and 2: Area 1 and Area S2 layers (None = cell type absent) ---
        for col, (pop_grid, area_name) in enumerate([(POP_GRID_S1, 'Area 1'),
                                                     (POP_GRID_S2, 'Area S2')], start=1):
            for i, layer_row in enumerate(pop_grid):
                layer_idx = [idx for idx in layer_row if idx is not None]
                plot_population_rates(axes[i, col], layer_idx, potentials, steps, start_plot,
                                      CELLTYPE_COLUMNS[:len(layer_idx)], stop_plot=stop_plot)
                axes[i, col].set_title(f'{area_name} - {LAYER_ROWS[i]}')
                axes[i, col].set_ylabel('mV')

        # set x-axis label for bottom row
        for ax in axes[3, :]:
            ax.set_xlabel('Time (ms)')

        fig.suptitle('Population Potentials')
        # all four cells of column 0 are used here, so the annotation goes below the
        # bottom row instead of into the empty cell plot_results leaves free
        annotate_fig(f'sI={np.round(sI, 4)}, g={np.round(g, 4)}, area={area}', xy=(-2, -0.35))
        sns.despine(trim=True)
        plt.tight_layout()
        figdir = os.path.join(figdir, 'single_simulations')
        if not os.path.exists(figdir):
            os.makedirs(figdir)
        plt.savefig(os.path.join(figdir, f'population_potentials_bEI-{sI}_g-{g}_area-{area}_Iextdur-{d}_Iextstr-{s}_Ibstr-{sb}.pdf'), dpi=300)
        #plt.show()


    elif area=='A1':

        # plot results for the S1 column 
        figS1, axs = plt.subplots(2, 2, figsize=(8, 5))  # Set figure size
    
        # Plot settings for all subplots
        axs_flat = axs.flatten()

        # Layer 2/3
        axs_flat[0].plot(steps[start_plot:stop_plot], potentials[:4].T[start_plot:stop_plot])
        axs_flat[0].set_title('L2/3')

        # Layer 4
        axs_flat[1].plot(steps[start_plot:stop_plot], potentials[4:4+3].T[start_plot:stop_plot])
        axs_flat[1].set_title('L4')
        # Layer 5
        axs_flat[2].plot(steps[start_plot:stop_plot], potentials[4+3:4+6].T[start_plot:stop_plot])
        # Layer 6
        axs_flat[3].plot(steps[start_plot:stop_plot], potentials[4+6:4+9].T[start_plot:stop_plot])
        
        axs_flat[0].legend(['E', 'PV', 'SST', 'VIP'])

        plt.show()

def plot_axis(axs, steps, start_plot, rates, idx_rates, color=None, stop_plot=None):
    axs.plot(steps[start_plot:stop_plot], rates[idx_rates].T[start_plot:stop_plot], linewidth=1, color=color)

def plot_population_rates(axs_op, idxs_pop, rates, steps, start_plot, labels, colors=None, stop_plot=None):
    """ Plot population rates for given labels and indices.

    colors : list, optional
        One colour per index; None (default) leaves the axes' colour cycle in charge.
    """
    legend_list = []
    for i, idx in enumerate(idxs_pop):
        plot_axis(axs_op, steps, start_plot, rates, idx,
                  color=None if colors is None else colors[i], stop_plot=stop_plot)
        legend_list.append(f'{labels[i]} {np.round(rates[idx].T[-1], 6)}')

    axs_op.legend(legend_list, loc='upper right')


def plot_results(rates, Iext, Ib, step_size, simulation_time, start_plot, sI, g, area, d, sb, s, figure_dir, stop_plot=None):
    steps = np.arange(step_size, simulation_time+step_size, step_size)*1e3
    fig, axs = plt.subplots(4, 3, figsize=(15, 15))  # Set figure size
    figure_style()

    # external input 
    axs_extI = axs[0][0]
    axs_extI.plot(steps[start_plot:stop_plot], Iext[start_plot:stop_plot], label='Iext rate')
    axs_extI.plot(steps[start_plot:stop_plot], Ib[start_plot:stop_plot], label='Ib rate')
    axs_extI.legend(title='')
    axs_extI.set_ylabel('Hz')
    # thalamus
    axs_thal = axs[1][0]
    # Thalamus block = last three populations, in array order: Thalamus E/VPM (-3,
    # the externally-driven population), Thalamus I/reticular (-2), POm (-1).
    # Order matches get_population_labels().
    axs_thal.plot(steps[start_plot:stop_plot], rates[-3].T[start_plot:stop_plot], color='green')
    axs_thal.plot(steps[start_plot:stop_plot], rates[-2].T[start_plot:stop_plot], color='purple')
    axs_thal.plot(steps[start_plot:stop_plot], rates[-1].T[start_plot:stop_plot], color='grey')
    axs_thal.legend(['Thalamus E', 'Thalamus I', 'POm'])
    axs_thal.set_ylabel('Hz')

    # area 3b
    axsA3b = axs[2][0]
    axsA3b.plot(steps[start_plot:stop_plot], rates[:4].T[start_plot:stop_plot], linewidth=1)
    axsA3b.legend([f'E {np.round(rates[0].T[-1], 6)}', f'PV {np.round(rates[1].T[-1], 6)}', f'SOM {np.round(rates[2].T[-1], 6)}', f'VIP {np.round(rates[3].T[-1], 6)}'])
    axsA3b.set_ylabel('Hz')

    # plot results for the S1 column 
    idxs_E = np.array([0+4, 4+4, 7+4, 10+4]) # indices of E populations in S1
    labels_pops = [['E1', 'E2', 'E3', 'E4'], ['PV1', 'PV2', 'PV3', 'PV4'], ['SST1', 'SST2', 'SST3', 'SST4'], ['VIP1']]
    
    # loop over populations for S1 and S2
    for i, labels in enumerate(labels_pops):
        if i<3:
            axs_pop = axs[i][1]
            plot_population_rates(axs_pop, idxs_E+i, rates, steps, start_plot, labels, stop_plot=stop_plot)
        else:
            # VIP
            axsVIPS1 = axs[i][1]
            axsVIPS1.plot(steps[start_plot:stop_plot], rates[3+4].T[start_plot:stop_plot], linewidth=1)
            axsVIPS1.legend([f'VIP1 {np.round(rates[3+4].T[-1], 6)}'])

        # plot results S2
        nr_pops = 13 # number of pops in S1
        if i<3:
            axs_pop = axs[i][2]
            plot_population_rates(axs_pop, idxs_E+i+nr_pops, rates, steps, start_plot, labels, stop_plot=stop_plot)
        else:
            # VIP
            axsVIPS2 = axs[i][2]
            axsVIPS2.plot(steps[start_plot:stop_plot], rates[3+nr_pops].T[start_plot:stop_plot], linewidth=1)
            axsVIPS2.legend([f'VIP1 {np.round(rates[3+nr_pops].T[-1], 6)}'])
    
    # Hide extra figure cell in col 0
    axs[3, 0].axis("off")
    
    # set x-axis label for bottom row
    for ax in axs[3, :]:
        ax.set_xlabel('Time (ms)')
    axs[2, 0].set_xlabel('Time (ms)')

    # set titles for each subplot
    fig.suptitle('Population Rates')
    axs[0][0].set_title('External input')
    axs[1][0].set_title('Thalamus')
    axs[2][0].set_title('Area 3b')
    axs[0][1].set_title('Area 1 (S1)')
    axs[0][2].set_title('Area S2')

    annotate_fig(f'sI={np.round(sI, 4)}, g={np.round(g, 4)}, area={area}')
    sns.despine(trim=True)
    plt.tight_layout() 
    plt.legend()
    figdir = os.path.join(figure_dir, 'single_simulations')
    if not os.path.exists(figdir):
        os.makedirs(figdir)

    plt.savefig(os.path.join(figdir, f'population_rates_bEI-{sI}_g-{g}_area-{area}_Iextdur-{d}_Iextstr-{s}_Ibstr-{sb}.pdf'), dpi=300)
    #plt.show()

def annotate_fig(dataname, xy=(-2, 0.3)):
    """ Write on the figure with which data it was generated, the date and the script name.

    xy : tuple, optional
        Position of the first line in axes fractions of the current axes; the two
        following lines are stacked 0.1 below it. The default puts the block in the
        empty bottom-left cell of the plot_results grid.
    """
    x, y = xy
    plt.annotate(dataname, xy=(x, y), xycoords='axes fraction', fontsize=8, ha='center')
    plt.annotate(datetime.datetime.now(), xy=(x, y - 0.1), xycoords='axes fraction', fontsize=8, ha='center')
    plt.annotate(f"generated in {os.path.basename(__file__)}", xy=(x, y - 0.2), xycoords='axes fraction', fontsize=8, ha='center')

