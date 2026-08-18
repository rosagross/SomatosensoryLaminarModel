import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.cm as cm 
from matplotlib.colors import ListedColormap, Normalize, BoundaryNorm
from plotting_style import figure_style
from helper_functions import load_simulation_data, compute_freq_spectrum

def sImulti_fingerprint_IextDurVsStr(data_df, gs, sIs, Ib_str, population, thalamus_source, figure_dir):
    '''
    1.3) MULTI fingerprint PLOT: Effect of input intensity and duration on dynamic behaviour (non-responsive, transfer and memory)
        - plot style: heatmap
        - y axis: input strength
        - x axis: input duration
        - columns: increasing g
        - rows: increasing sI
        - measure: dynamic function (the "finger print")

    Parameters:
    -----------
    - gs: list of global coupling strengths
    - sIs: list of E/I balance values
    '''

    cbar_ticks = ['non-responsive', 'transfer', 'memory']
    data_df = data_df[data_df['population'] == population]
    data_df = data_df[data_df['globalCoupling'].isin(gs)]
    data_df = data_df[data_df['strength_I'].isin(sIs)]
    data_df = data_df[data_df['BckgndInputStrength'] == Ib_str]
    data_df['InputDuration'] = data_df['InputDuration'].round(4)

    nrows = len(sIs)
    ncols = len(gs)
    fig, ax = plt.subplots(nrows, ncols, figsize=(3 * ncols, 2.5 * nrows), sharex=True, sharey=True)
    ax = np.atleast_2d(ax)

    # Define a fixed discrete colormap with colors for 1, 2, and 3
    n = 3  # there are three different functions of the dynamics --> make discrete colormap
    colors = sns.color_palette("Pastel2", n)
    cmap = ListedColormap(colors)

    # Define the boundaries for normalization
    bounds = [0.5, 1.5, 2.5, 3.5]
    norm = BoundaryNorm(bounds, ncolors=cmap.N)

    for r, sI in enumerate(sIs):
        for c, g in enumerate(gs):
            axis = ax[r, c]
            plot_df = data_df[(data_df['strength_I'] == sI) & (data_df['globalCoupling'] == g)]
            data_heatmap = plot_df.pivot(index='InputStrength', columns='InputDuration', values='dynamic_function_potential')
            show_cbar = (r == nrows - 1) and (c == ncols - 1)
            heat_ax = sns.heatmap(data_heatmap, cmap=cmap, norm=norm, ax=axis, cbar=show_cbar, vmin=1, vmax=3)

            if show_cbar:
                cbar = heat_ax.collections[0].colorbar
                cbar.set_ticks([1, 2, 3])
                cbar.set_ticklabels(cbar_ticks)

            if r == nrows - 1:
                axis.set_xlabel('Input Duration')
            else:
                axis.set_xlabel('')
            if c == 0:
                axis.set_ylabel('Input Strength')
            else:
                axis.set_ylabel('')
            axis.invert_yaxis()

            if r == 0:
                axis.set_title(f'g = {g}')
            if c == 0:
                axis.text(-0.3, 0.5, f'sI = {sI}', va='center', rotation=90, transform=axis.transAxes)

    # global labels
    fig.text(0.02, 0.5, 'Inhibition strength sI (increasing downwards)', ha='right', va='center', rotation='vertical', fontsize=20)
    fig.text(0.5, 0.98, 'Coupling strength g (increasing to the right)', ha='center', va='top', fontsize=20)

    plt.tight_layout(h_pad=1, w_pad=1, rect=[0.05, 0.05, 0.95, 0.93])
    figure_name = f'fingerprint_IextDurVsStr_{population}pop_tauVisual_{thalamus_source}.pdf'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
    plt.show()


def fingerprint_sI_vs_IextStr(data_df, g, Iext_dur, Ib_str, population, thalamus_source, figure_dir):
    """
    Fingerprint heatmap: sI (x) vs input strength (y) for a fixed coupling strength.
    Uses dynamic_function_potential values.
    """
    cbar_ticks = ['non-responsive', 'transfer', 'memory']
    plot_df = data_df[
        (data_df['population'] == population) &
        (data_df['globalCoupling'] == g) &
        (data_df['BckgndInputStrength'] == Ib_str) &
        (data_df['InputDuration'] == Iext_dur)
    ].copy()

    # discrete colormap
    n = 3
    colors = sns.color_palette("Pastel2", n)
    cmap = ListedColormap(colors)
    bounds = [0.5, 1.5, 2.5, 3.5]
    norm = BoundaryNorm(bounds, ncolors=cmap.N)

    data_heatmap = plot_df.pivot(index='InputStrength', columns='strength_I', values='dynamic_function_potential')
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    heat_ax = sns.heatmap(data_heatmap, cmap=cmap, norm=norm, ax=ax, cbar=True)
    cbar = heat_ax.collections[0].colorbar
    cbar.set_ticks([1, 2, 3])
    cbar.set_ticklabels(cbar_ticks)
    ax.invert_yaxis()
    ax.set_xlabel('sI')
    ax.set_ylabel('Input Strength')
    ax.set_title(f'Fingerprint: sI vs Input Strength, g={g}')

    figure_name = f'fingerprint_sI_vs_IextStr_{population}_g{g}_Iextdur{Iext_dur}_Ib{Ib_str}_tauVisual_{thalamus_source}.pdf'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
    plt.show()


def fingerprint_g_vs_IextStr(data_df, sI, Iext_dur, Ib_str, population, thalamus_source, figure_dir):
    """
    Fingerprint heatmap: coupling strength (x) vs input strength (y) for a fixed sI.
    Uses dynamic_function_potential values.
    """
    cbar_ticks = ['non-responsive', 'transfer', 'memory']
    plot_df = data_df[
        (data_df['population'] == population) &
        (data_df['strength_I'] == sI) &
        (data_df['BckgndInputStrength'] == Ib_str) &
        (data_df['InputDuration'] == Iext_dur)
    ].copy()

    # discrete colormap
    n = 3
    colors = sns.color_palette("Pastel2", n)
    cmap = ListedColormap(colors)
    bounds = [0.5, 1.5, 2.5, 3.5]
    norm = BoundaryNorm(bounds, ncolors=cmap.N)

    data_heatmap = plot_df.pivot(index='InputStrength', columns='globalCoupling', values='dynamic_function_potential')
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    heat_ax = sns.heatmap(data_heatmap, cmap=cmap, norm=norm, ax=ax, cbar=True)
    cbar = heat_ax.collections[0].colorbar
    cbar.set_ticks([1, 2, 3])
    cbar.set_ticklabels(cbar_ticks)
    ax.invert_yaxis()
    ax.set_xlabel('Global Coupling Strength (g)')
    ax.set_ylabel('Input Strength')
    ax.set_title(f'Fingerprint: g vs Input Strength, sI={sI}')

    figure_name = f'fingerprint_g_vs_IextStr_{population}_sI{sI}_Iextdur{Iext_dur}_Ib{Ib_str}_tauVisual_{thalamus_source}.pdf'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
    plt.show()


def multiPop_heatmap_IextDurVsStr(data_df, gs, sI, Ib_str, populations, rate_measure, figure_dir):
    """
    MULTI heatmap PLOT: Effect of input intensity and duration on firing rates
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    - subplot columns: populations
    - subplot rows: coupling strengths

    Parameters:
    -----------
    - g: list of global coupling strengths
    - sI: float balance E/I value
    - Ib_str: float, background input strength
    - populations: list of populations to plot
    - rate_measure: measure to plot (e.g. 'longtermVSbaseline_rate')
    """

    fig, axes = plt.subplots(len(gs), len(populations), figsize=(20,15) ,sharex=True, sharey=True)

    # Create a single colorbar axis
    #cbar_ax = fig.add_axes([1.01, 0.3, 0.02, 0.4])
    #cbar_ax.set_title(rate_measure)
    #cbar_ax.tick_params(labelsize=12) 


    for i,g in enumerate(gs):
        for j,p in enumerate(populations):

            minmax_df = data_df[data_df['globalCoupling']==g]
            minmax_df = minmax_df[minmax_df['strength_I']==sI]
            minmax_df = minmax_df[minmax_df['BckgndInputStrength']==Ib_str]
            minmax_df = minmax_df[minmax_df['population']==p]
            minmax_df['InputDuration'] = minmax_df['InputDuration'].round(4)
            minmax_df[rate_measure] = minmax_df[rate_measure].round(5)

            data_heatmap = minmax_df.pivot(index='InputStrength',columns='InputDuration', values=rate_measure)

            if minmax_df[rate_measure].isna().all():
                print(f'No measure value for {rate_measure} at g: {g}, sI: {sI}, pop: {p}')
                sns.heatmap(data_heatmap, cmap=ListedColormap(['green']), ax=axes[i,j], norm = Normalize(vmin=0, vmax=1))
            else:
                sns.heatmap(data_heatmap, cmap='magma', ax=axes[i,j], vmin=-2, vmax=2)

            cbar = axes[i, j].collections[0].colorbar
            # here set the labelsize by 20
            cbar.ax.tick_params(labelsize=12)
            axes[i, j].invert_yaxis()
            axes[i, j].set_ylabel('')
            axes[i, j].set_xlabel('')
            axes[i, j].tick_params(axis='both', labelsize=12)
            axes[len(gs)-1, j].set_xlabel('Input Duration')
            axes[0,j].set_title(f'pop: {p}')
            axes[i,0].set_ylabel(f'g: {g}, sI: {sI}', rotation=0, labelpad=60)

    fig.text(0.05, 0.2, 'Input Strength', va='center', rotation='vertical')
    fig.text(0.05, 0.5, 'Input Strength', va='center', rotation='vertical')
    fig.text(0.05, 0.83, 'Input Strength', va='center', rotation='vertical')

    plt.tight_layout(h_pad=15)
    figure_name = f'inputDurationVSinputStrength_{populations[0][0]}pop_sI{sI}_Ibstr{Ib_str}_{rate_measure}_tauVisual.png'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
    plt.show()


def multiLayer_couplingOnLongeterm_diffRate(data_df, Iext_dur, Iext_str, Ib_str, sI, thalamus_source, figure_dir):
    """
    Effect of Coupling Strengths on Longterm/steady state
    Plot difference between Mininum and Maximum Firing rates 
    """

    data_df = data_df[data_df['InputDuration']==Iext_dur]
    data_df = data_df[data_df['strength_I']==sI]
    data_df = data_df[data_df['InputStrength']==Iext_str]
    data_df = data_df[data_df['BckgndInputStrength']==Ib_str]

    # separate data in layers
    layersS1 = [['E1', 'PV1', 'SST1', 'VIP1'], ['E2', 'PV2', 'SST2'], ['E3', 'PV3', 'SST3'], ['E4', 'PV4', 'SST4']]
    layersS2 = [['E1S2', 'PV1S2', 'SST1S2', 'VIP1S2'], ['E2S2', 'PV2S2', 'SST2S2'], ['E3S2', 'PV3S2', 'SST3S2'], ['E4S2', 'PV4S2', 'SST4S2']]
    layers_all = [layersS1, layersS2]

    # plot results
    fig, axes = plt.subplots(4, 2, figsize=(10, 6), sharey=False, sharex=True)  # Set figure size
    fig.suptitle("Difference between min. and max. firing rate \nduring late longterm phase")

    for layers, axs in zip(layers_all, axes.T):
        for l, ax in zip(layers, axs):
            layer_df = data_df[data_df['population'].isin(l)]
            #sns.lineplot(layer_df, y='minRate_longterm', x='globalCoupling', hue='population', ax=ax)
            #sns.lineplot(layer_df, y='maxRate_longterm', x='globalCoupling', hue='population', ax=ax, legend=False)
            sns.lineplot(layer_df, y='diffRate_lateLongterm', x='globalCoupling', hue='population', ax=ax, legend=True)
            ax.set_ylabel('Rate (Hz)')
            ax.set_xlabel('Global Coupling Strength (g)')
            ax.legend("")

        axs[0].legend(prop={'size':8}, loc='upper right')
        axs[0].set_title(f'Layer 2/3')
        axs[1].set_title(f'Layer 4')
        axs[2].set_title(f'Layer 5')
        axs[3].set_title(f'Layer 6')
        axs[1].set_ylim([0,150])
        axs[2].set_ylim([0,150])
        axs[3].set_ylim([0,150])
        axs[0].set_ylim([0,150])

    sns.despine(trim=True)
    plt.tight_layout() 
    plt.annotate(f'min and max Rate_lateLongterm, \nInput Duration:{Iext_dur} \nInput Strength:{Iext_str} Background Input:{Ib_str} \nE-I Balance:{sI}', xy=(0, 0),
                xycoords='figure fraction', xytext=(1.14, 0.17), textcoords='figure fraction', ha='center', fontsize=10)
    figure_name = f'CouplingOnLongterm_AllLayers_Iextdur{Iext_dur}_Iextstr{Iext_str}_Ibstr{Ib_str}_sI{sI}_tauVisual_{thalamus_source}.pdf'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
    plt.show()


def heatmap_frequency_IextDurVsStr(data_df, g, sI, Ib_str, population, window_prefix, signal_type, Iext_dur, Iext_str, thalamus_source, figure_dir, vmin=None, vmax=None):
    """
    Heatmap of dominant frequency vs input duration/strength for a given population.
    signal_type: "Rate" or "Potential"
    window_prefix: "baseline", "duringInput", "lateLongterm"
    """
    col = f"freq{signal_type}_{window_prefix}"

    plot_df = data_df[
        (data_df["globalCoupling"] == g) &
        (data_df["strength_I"] == sI) &
        (data_df["BckgndInputStrength"] == Ib_str) &
        (data_df["population"] == population)
    ].copy()

    plot_df["InputDuration"] = plot_df["InputDuration"].round(4)
    data_heatmap = plot_df.pivot(index="InputStrength", columns="InputDuration", values=col)
    ax = sns.heatmap(data_heatmap, cmap="magma", vmin=vmin, vmax=vmax)
    ax.invert_yaxis()
    ax.set_xlabel("Input Duration")
    ax.set_ylabel("Input Strength")
    ax.set_title(f"{signal_type} freq ({window_prefix}) - pop {population}, g {g}, sI {sI}, Iext_dur {Iext_dur}")

    figure_name = f"freqHeatmap_{signal_type}_{window_prefix}_g{g}_sI{sI}_Ib{Ib_str}_pop{population}.pdf"
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches="tight")
    plt.show()

    # make a plot for A3b too
    cells = ['E3b', 'PV3b', 'SST3b', 'VIP3b']
    fig, axes = plt.subplots(1, 1, figsize=(10, 6), sharey=False, sharex=True)
    fig.suptitle(f"Dominant {signal_type} frequency during {window_prefix} (A3b)")
    for cell in cells:
        layer_df = data_df[data_df["population"] == cell]
        sns.lineplot(layer_df, y=col, x="globalCoupling", ax=axes, legend=True)
        axes.set_ylabel("Frequency (Hz)")
        axes.set_xlabel("Global Coupling Strength (g)")

    sns.despine(trim=True)
    figure_name = (
        f"A3b_CouplingOnFrequency_{signal_type}_{window_prefix}_Iextdur{Iext_dur}_"
        f"Iextstr{Iext_str}_Ibstr{Ib_str}_sI{sI}_tauVisual_{thalamus_source}.pdf"
    )
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches="tight")
    plt.show()


def multiLayer_couplingOnFrequency(data_df, Iext_dur, Iext_str, Ib_str, sI, window_prefix, signal_type, thalamus_source, figure_dir):
    """
    Effect of Coupling Strengths on dominant frequency.
    signal_type: "Rate" or "Potential"
    window_prefix: "baseline", "duringInput", "lateLongterm"
    """
    col = f"freq{signal_type}_{window_prefix}"

    data_df = data_df[data_df["InputDuration"] == Iext_dur]
    data_df = data_df[data_df["strength_I"] == sI]
    data_df = data_df[data_df["InputStrength"] == Iext_str]
    data_df = data_df[data_df["BckgndInputStrength"] == Ib_str]

    layersS1 = [['E1', 'PV1', 'SST1', 'VIP1'], ['E2', 'PV2', 'SST2'], ['E3', 'PV3', 'SST3'], ['E4', 'PV4', 'SST4']]
    layersS2 = [['E1S2', 'PV1S2', 'SST1S2', 'VIP1S2'], ['E2S2', 'PV2S2', 'SST2S2'], ['E3S2', 'PV3S2', 'SST3S2'], ['E4S2', 'PV4S2', 'SST4S2']]
    layers_all = [layersS1, layersS2]
    area_titles = ["S1", "S2"]
    layer_labels = ["Layer 2/3", "Layer 4", "Layer 5", "Layer 6"]

    # Fixed colour per cell type so a population's colour is stable across panels.
    dark2 = sns.color_palette("Dark2")
    celltype_colors = {"E": dark2[2], "PV": dark2[1], "SST": dark2[0], "VIP": dark2[3]}

    def _celltype(pop):
        for ct in ("SST", "VIP", "PV", "E"):
            if pop.startswith(ct):
                return ct
        return "E"

    nrows, ncols = 4, 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(7, 8), sharey=False, sharex=True)

    for c, (layers, axs) in enumerate(zip(layers_all, axes.T)):
        for r, (l, ax) in enumerate(zip(layers, axs)):
            layer_df = data_df[data_df["population"].isin(l)]
            palette = {p: celltype_colors[_celltype(p)] for p in l}
            sns.lineplot(layer_df, y=col, x="globalCoupling", hue="population",
                         palette=palette, marker="o", markersize=3, ax=ax, legend=True)
            # one compact, deduplicated legend per panel (cell types differ per layer)
            ax.legend(prop={"size": 7}, loc="best", frameon=False)
            # shared axis labels: only bottom row / left column
            ax.set_ylabel("Peak frequency (Hz)" if c == 0 else "")
            ax.set_xlabel("Global coupling $g$" if r == nrows - 1 else "")
            # layer label on the left, area title on the top
            if c == 0:
                ax.text(-0.35, 0.5, layer_labels[r], transform=ax.transAxes,
                        rotation=90, va="center", ha="center", fontweight="bold")
            if r == 0:
                ax.set_title(area_titles[c], fontweight="bold")

    fig.suptitle(f"Dominant {signal_type} frequency during {window_prefix}", fontweight="bold")

    sns.despine(trim=True)
    plt.tight_layout(rect=[0.04, 0.05, 1, 0.96])
    fig.text(
        0.5, 0.005,
        f"Input duration {Iext_dur} | Input strength {Iext_str} | "
        f"Background input {Ib_str} | E-I balance {sI}",
        ha="center", va="bottom", fontsize=8, color="0.3",
    )
    figure_name = (
        f"CouplingOnFrequency_{signal_type}_{window_prefix}_Iextdur{Iext_dur}_"
        f"Iextstr{Iext_str}_Ibstr{Ib_str}_sI{sI}_tauVisual_{thalamus_source}"
    )
    plt.savefig(os.path.join(figure_dir, figure_name + ".pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(figure_dir, figure_name + ".png"), bbox_inches="tight", dpi=300)
    plt.show()


def scatter_frequency_vs_diff(data_df, window_prefix, signal_type, population, figure_dir):
    """
    Scatter plot of dominant frequency vs min-max diff for a given population.
    """
    freq_col = f"freq{signal_type}_{window_prefix}"
    diff_col = f"diff{signal_type}_{window_prefix}"

    plot_df = data_df[data_df["population"] == population].copy()
    ax = sns.scatterplot(data=plot_df, x=diff_col, y=freq_col)
    ax.set_xlabel(f"{diff_col}")
    ax.set_ylabel(f"{freq_col} (Hz)")
    ax.set_title(f"{signal_type} freq vs amplitude ({window_prefix}) - {population}")

    figure_name = f"freqVsDiff_{signal_type}_{window_prefix}_pop{population}.pdf"
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches="tight")
    plt.show()


def heatmap_frequency_coupling_vs_sI(data_df, area, window_prefix, signal_type, Iext_dur, Iext_str, Ib_str, figure_dir, vmin=None, vmax=None):
    """
    Publication-ready grid of dominant-frequency heatmaps (global coupling x,
    E-I balance sI y) for all populations of an area, arranged as
    layer (rows) x cell type (columns).
    area: "A1", "S2", or "A3b"
    """
    col = f"freq{signal_type}_{window_prefix}"
    if col not in data_df.columns:
        raise ValueError(f"Column '{col}' not found in data frame.")

    # rows = cortical layers, columns = cell type (E, PV, SST, VIP);
    # None marks a cell type that does not exist in that layer.
    if area == "A1":
        layers = [['E1', 'PV1', 'SST1', 'VIP1'],
                  ['E2', 'PV2', 'SST2', None],
                  ['E3', 'PV3', 'SST3', None],
                  ['E4', 'PV4', 'SST4', None]]
        row_labels = ["Layer 2/3", "Layer 4", "Layer 5", "Layer 6"]
    elif area == "S2":
        layers = [['E1S2', 'PV1S2', 'SST1S2', 'VIP1S2'],
                  ['E2S2', 'PV2S2', 'SST2S2', None],
                  ['E3S2', 'PV3S2', 'SST3S2', None],
                  ['E4S2', 'PV4S2', 'SST4S2', None]]
        row_labels = ["Layer 2/3", "Layer 4", "Layer 5", "Layer 6"]
    elif area == "A3b":
        layers = [['E3b', 'PV3b', 'SST3b', 'VIP3b']]
        row_labels = ["A3b"]
    else:
        raise ValueError(f"Unsupported area: {area}. Use 'A1', 'S2', or 'A3b'.")

    col_labels = ["Excitatory", "PV", "SST", "VIP"]
    populations = [p for row in layers for p in row if p is not None]

    plot_df = data_df[
        (data_df["InputDuration"] == Iext_dur) &
        (data_df["InputStrength"] == Iext_str) &
        (data_df["BckgndInputStrength"] == Ib_str) &
        (data_df["population"].isin(populations))
    ].copy()

    nrows = len(layers)
    ncols = 4
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.2 * ncols, 2 * nrows), sharex=True, sharey=True)
    axes = np.atleast_2d(axes)

    # single shared colorbar
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap="magma", norm=norm)
    sm.set_array([])

    for r, row in enumerate(layers):
        for c, pop in enumerate(row):
            ax = axes[r, c]

            if pop is None:
                ax.axis("off")
                continue

            pop_df = plot_df[plot_df["population"] == pop]
            data_heatmap = pop_df.pivot(index="strength_I", columns="globalCoupling", values=col)

            # seaborn.heatmap fails on empty matrices ("zero-size array to reduction")
            # which can happen when no rows match the current population/filter values.
            if data_heatmap.empty or data_heatmap.shape[0] == 0 or data_heatmap.shape[1] == 0 or np.all(np.isnan(data_heatmap.to_numpy())):
                ax.axis("off")
                ax.text(0.5, 0.5, f"{pop}\n(no data)", ha="center", va="center", transform=ax.transAxes)
                continue

            # round the axis tick labels for a cleaner look
            data_heatmap.index = np.round(data_heatmap.index.astype(float), 2)
            data_heatmap.columns = np.round(data_heatmap.columns.astype(float), 1)

            sns.heatmap(data_heatmap, cmap="magma", vmin=vmin, vmax=vmax, ax=ax, cbar=False)
            ax.invert_yaxis()
            ax.tick_params(axis="both", length=0)

            # shared axis labels: bottom row / left column only
            ax.set_xlabel("Global coupling $g$" if r == nrows - 1 else "")
            ax.set_ylabel("Inhibitory coupling $s_I$" if c == 0 else "")

            # cell-type column headers (top row) and layer row labels (left column)
            if r == 0:
                ax.set_title(col_labels[c], fontweight="bold")
            if c == 0:
                ax.text(-0.45, 0.5, row_labels[r], transform=ax.transAxes,
                        rotation=90, va="center", ha="center", fontweight="bold")

    cbar = fig.colorbar(sm, ax=axes, orientation="vertical", fraction=0.02, pad=0.02)
    cbar.set_label("Peak frequency (Hz)")

    area_title = {"A1": "S1", "S2": "S2", "A3b": "A3b"}.get(area, area)
    # place the title above the column headers (higher offset needed for few rows)
    suptitle_y = 1.02 + 0.10 / nrows
    fig.suptitle(f"Dominant {signal_type} frequency ({window_prefix}) - {area_title}",
                 fontweight="bold", y=suptitle_y)
    figure_name = (
        f"freqHeatmap_gVsSI_{signal_type}_{window_prefix}_Iextdur{Iext_dur}_"
        f"Iextstr{Iext_str}_Ib{Ib_str}_area{area}"
    )
    plt.savefig(os.path.join(figure_dir, figure_name + ".pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(figure_dir, figure_name + ".png"), bbox_inches="tight", dpi=300)
    plt.show()


def heatmap_frequency_couplingVsSI_singlePop(data_df, population, window_prefix, signal_type,
                                             Iext_dur, Iext_str, Ib_str, figure_dir,
                                             vmin=None, vmax=None):
    """
    Publication-ready heatmap of the dominant frequency peak for a SINGLE chosen
    population, as a function of global coupling (x) and E-I balance sI (y).

    A focused, single-panel counterpart to heatmap_frequency_coupling_vs_sI
    (which draws every population of a whole area).

    signal_type:   "Rate" or "Potential"
    window_prefix: "baseline", "duringInput", "lateLongterm", ...
    """
    col = f"freq{signal_type}_{window_prefix}"
    if col not in data_df.columns:
        raise ValueError(f"Column '{col}' not found in data frame.")

    plot_df = data_df[
        (data_df["InputDuration"] == Iext_dur) &
        (data_df["InputStrength"] == Iext_str) &
        (data_df["BckgndInputStrength"] == Ib_str) &
        (data_df["population"] == population)
    ].copy()

    data_heatmap = plot_df.pivot(index="strength_I", columns="globalCoupling", values=col)
    # round the axis tick labels for a cleaner look
    data_heatmap.index = np.round(data_heatmap.index.astype(float), 3)
    data_heatmap.columns = np.round(data_heatmap.columns.astype(float), 1)

    fig, ax = plt.subplots(figsize=(4.5, 4))

    # seaborn.heatmap fails on empty/all-NaN matrices ("zero-size array to reduction")
    if (data_heatmap.empty or data_heatmap.shape[0] == 0 or data_heatmap.shape[1] == 0
            or np.all(np.isnan(data_heatmap.to_numpy()))):
        ax.axis("off")
        ax.text(0.5, 0.5, f"{population}\n(no data)", ha="center", va="center",
                transform=ax.transAxes)
    else:
        sns.heatmap(
            data_heatmap, cmap="magma", vmin=vmin, vmax=vmax, ax=ax,
            cbar_kws={"label": "Peak frequency (Hz)"},
        )
        ax.invert_yaxis()
        ax.set_xlabel("Global coupling $g$")
        ax.set_ylabel("Inhibitory coupling $s_I$")
        ax.tick_params(axis="both", length=0)

    ax.set_title(
        f"{signal_type} peak frequency ({window_prefix}) - {population}\n"
        f"Iext dur {Iext_dur}, Iext str {Iext_str}, Ib {Ib_str}",
        fontsize=10,
    )
    plt.tight_layout()

    figure_name = (
        f"freqHeatmap_gVsSI_singlePop_{signal_type}_{window_prefix}_Iextdur{Iext_dur}_"
        f"Iextstr{Iext_str}_Ib{Ib_str}_pop{population}"
    )
    plt.savefig(os.path.join(figure_dir, figure_name + ".pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(figure_dir, figure_name + ".png"), bbox_inches="tight", dpi=300)
    plt.show()


def baseline_spectrum_by_coupling(sweep_name, sweep_values, g, g_inter, sI, Ib_str,
                                  Iext_dur, Iext_str, input_onset, step_size, sample_dur,
                                  offset, thal_cellcounts, bI_cellcounts, extI_cellcounts,
                                  input_type, raw_dir, figure_dir, suffix='', fmax=80.0,
                                  g_thalPOm=1.0, Ib_noise_std=None):
    """
    Grid of baseline power spectra (summed potential), one subplot per excitatory
    population, one curve per swept coupling value (colour = value).

    Parameters:
    -----------
    - sweep_name: 'g' -> vary global coupling; 'g_inter' -> vary inter-cortical coupling
    - sweep_values: list of values for the swept coupling parameter
    - g, g_inter: fixed coupling values (the one matching sweep_name is ignored/overwritten)
    - sI, Ib_str, Iext_dur, Iext_str: fixed simulation parameters
    - fmax: maximum displayed frequency (Hz)
    """
    if sweep_name not in ('g', 'g_inter'):
        raise ValueError(f"sweep_name must be 'g' or 'g_inter', got {sweep_name}")

    exc_pops = ['E1', 'E2', 'E3', 'E4',
                'E1S2', 'E2S2', 'E3S2', 'E4S2', 'E3b']
    titles = {'E1': 'Area 1\nLayer 2/3', 'E2': 'Area 1\nLayer 4',
              'E3': 'Area 1\nLayer 5', 'E4': 'Area 1\nLayer 6',
              'E1S2': 'S2\nLayer 2/3', 'E2S2': 'S2\nLayer 4',
              'E3S2': 'S2\nLayer 5', 'E4S2': 'S2\nLayer 6',
              'E3b': 'A3b\nLayer 2/3'}

    baseline_start = int((input_onset - (sample_dur + offset)) / step_size)
    baseline_stop = int(baseline_start + sample_dur / step_size)

    # colour map across swept values; truncate the light (near-white) low end
    # of 'Greens' so every trace is visible on a white poster background.
    base_cmap = cm.get_cmap('Greens')
    n_vals = max(len(sweep_values), 1)
    trace_colors = [base_cmap(v) for v in np.linspace(0.35, 1.0, n_vals)]

    fig, axes = plt.subplots(3, 3, figsize=(10, 8), sharex=True)
    axes = axes.flatten()

    for k, val in enumerate(sweep_values):
        g_k = val if sweep_name == 'g' else g
        g_inter_k = val if sweep_name == 'g_inter' else g_inter
        try:
            _, potentials_df, _ = load_simulation_data(
                g_k, g_inter_k, sI, Ib_str, Iext_dur, Iext_str, input_onset,
                thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type,
                raw_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std)
        except (FileNotFoundError, OSError):
            print(f'Missing run: {sweep_name}={val} (g={g_k}, g_inter={g_inter_k}) - skipping')
            continue

        window = potentials_df.iloc[baseline_start:baseline_stop]
        for ax, pop in zip(axes, exc_pops):
            freqs, power = compute_freq_spectrum(window[pop].values, step_size)
            mask = freqs <= fmax
            ax.plot(freqs[mask], power[mask], color=trace_colors[k], lw=1.2, label=f'{val}')

    for ax, pop in zip(axes, exc_pops):
        ax.set_title(titles[pop], fontweight='bold')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Power (mV$^2$)')

    legend_title = {'g': 'Intracortical coupling', 'g_inter': 'Intercortical coupling'}.get(sweep_name, sweep_name)
    axes[0].legend(title=legend_title, fontsize=11, title_fontsize=12, ncol=2)
    sns.despine(trim=True)
    fig.suptitle(f'Baseline power spectrum vs {sweep_name} '
                 f'(sI={sI}, Ib={Ib_str}, Iextstr={Iext_str}, Iextdur={Iext_dur})')
    plt.tight_layout()

    figure_name = (f'baselineSpectrum_vs_{sweep_name}_sI{sI}_Ib{Ib_str}_'
                   f'Iextstr{Iext_str}_Iextdur{Iext_dur}.pdf')
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
    plt.show()


def backgroundInputOnPeakPower(data_df, g, sI, Iext_dur, Iext_str,
                               window_prefix, signal_type, thalamus_source, figure_dir):
    """
    Influence of background input strength on the peak-frequency power for each
    excitatory (E) population across the three ROIs (S1, S2, A3b), all other
    parameters fixed. One panel per ROI, one line per cortical layer.

    signal_type:   "Rate" or "Potential"
    window_prefix: "baseline", "duringInput", "lateLongterm", ...
    """
    col = f"fftPower{signal_type}_{window_prefix}"
    if col not in data_df.columns:
        raise ValueError(f"Column '{col}' not found in data frame.")

    plot_df = data_df[
        (data_df["globalCoupling"] == g) &
        (data_df["strength_I"] == sI) &
        (data_df["InputDuration"] == Iext_dur) &
        (data_df["InputStrength"] == Iext_str)
    ].copy()

    # ROI -> (E populations by layer, layer labels)
    rois = {
        "S1":  (['E1', 'E2', 'E3', 'E4'], ["Layer 2/3", "Layer 4", "Layer 5", "Layer 6"]),
        "S2":  (['E1S2', 'E2S2', 'E3S2', 'E4S2'], ["Layer 2/3", "Layer 4", "Layer 5", "Layer 6"]),
        "A3b": (['E3b'], ["Layer 2/3"]),
    }
    # cortical-depth colour gradient (L2/3 -> L6)
    depth_colors = sns.color_palette("viridis", 4)

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5), sharex=True)

    for ax, (roi, (pops, labels)) in zip(axes, rois.items()):
        for pop, lab, color in zip(pops, labels, depth_colors):
            pop_df = plot_df[plot_df["population"] == pop].sort_values("BckgndInputStrength")
            ax.plot(pop_df["BckgndInputStrength"], pop_df[col],
                    marker="o", color=color, label=lab)
        ax.set_title(roi, fontweight="bold")
        ax.set_xlabel("Background input strength")
        ax.legend(title="Layer", prop={"size": 8})

    axes[0].set_ylabel("Peak frequency power")

    fig.suptitle(
        f"Background input vs {signal_type} peak-frequency power ({window_prefix})   "
        f"[g {g}, sI {sI}, Iext dur {Iext_dur}, Iext str {Iext_str}]",
        fontsize=10,
    )
    sns.despine(trim=True)
    plt.tight_layout()

    figure_name = (
        f"bckgndInputOnPeakPower_{signal_type}_{window_prefix}_g{g}_sI{sI}_"
        f"Iextdur{Iext_dur}_Iextstr{Iext_str}_tauVisual_{thalamus_source}"
    )
    plt.savefig(os.path.join(figure_dir, figure_name + ".pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(figure_dir, figure_name + ".png"), bbox_inches="tight", dpi=300)
    plt.show()


def multisI_couplingOnMinmaxRate(data_df, sI, Ib_str, population, rate_measure, thalamus_source, figure_dir):
    '''
    3.1) Effect of coupling strength on the firing-rate range (min/max) of a population.
    - x-axis: global coupling strength
    - y-axis: min and max firing rate (shaded band = oscillation range)
    - one subplot per E-I balance (sI) value
    - one colored band per background-input strength (Ib_str)
    - choose rate measure ('lateLongterm', 'immediateLongterm', 'duringInput', 'baseline')
    '''
    figure_style()

    data_df = data_df[data_df['population'] == population]
    Ib_str = sorted(Ib_str)
    palette = sns.color_palette('Blues', len(Ib_str))

    fig, axes = plt.subplots(len(sI), 1, figsize=(3.5, 1.6 * len(sI)), sharex=True)
    axes = np.atleast_1d(axes)
    fig.suptitle(f'Min–max firing rate of {population}')

    for ax, b in zip(axes, sI):
        data_sI_df = data_df[data_df['strength_I'] == b]

        for i, s in enumerate(Ib_str):
            Istrength_df = data_sI_df[data_sI_df['BckgndInputStrength'] == s].sort_values('globalCoupling')
            x = Istrength_df['globalCoupling']
            y_min = Istrength_df[f'minRate_{rate_measure}']
            y_max = Istrength_df[f'maxRate_{rate_measure}']
            ax.plot(x, y_max, color=palette[i], linewidth=1)
            ax.plot(x, y_min, color=palette[i], linewidth=1)
            ax.fill_between(x, y_min, y_max, color=palette[i], alpha=0.2, linewidth=0)

        ax.set_ylabel('Rate (Hz)')
        ax.set_title(f'Inhibitory coupling = {b}')

    axes[-1].set_xlabel('Global coupling strength')

    # single shared legend for the background-input colors
    handles = [plt.Line2D([0], [0], color=palette[i], linewidth=2) for i in range(len(Ib_str))]
    fig.legend(handles, [str(s) for s in Ib_str], title='Background input',
               loc='center left', bbox_to_anchor=(1.0, 0.5))

    sns.despine(trim=True)
    plt.tight_layout()
    figure_name = (f'multisI_couplingOnMinmaxRate_{rate_measure}_pop{population}_'
                   f'Ibstr{Ib_str}_sI{sI}_tauVisual_{thalamus_source}')
    plt.savefig(os.path.join(figure_dir, figure_name + '.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(figure_dir, figure_name + '.png'), bbox_inches='tight', dpi=300)
    plt.show()


def multisI_couplingOnDiffRate(data_df, sI, Ib_str, population, rate_measure, thalamus_source, figure_dir):
    '''
    3.1b) Effect of coupling strength on the firing-rate oscillation amplitude of a population.
    Same layout as multisI_couplingOnMinmaxRate, but the y-axis is the difference between the
    maximum and minimum rate (diffRate_{rate_measure}) rather than the two absolute rates.
    - x-axis: global coupling strength
    - y-axis: max - min firing rate (oscillation amplitude)
    - one subplot per E-I balance (sI) value
    - one colored curve per background-input strength (Ib_str)
    - choose rate measure ('lateLongterm', 'immediateLongterm', 'duringInput', 'baseline')
    '''
    figure_style()

    data_df = data_df[data_df['population'] == population]
    Ib_str = sorted(Ib_str)
    palette = sns.color_palette('viridis', len(Ib_str))

    fig, axes = plt.subplots(len(sI), 1, figsize=(3.5, 1.6 * len(sI)), sharex=True)
    axes = np.atleast_1d(axes)
    fig.suptitle(f'Firing-rate range (max–min) of {population}')

    for ax, b in zip(axes, sI):
        data_sI_df = data_df[data_df['strength_I'] == b]

        for i, s in enumerate(Ib_str):
            Istrength_df = data_sI_df[data_sI_df['BckgndInputStrength'] == s].sort_values('globalCoupling')
            ax.plot(Istrength_df['globalCoupling'], Istrength_df[f'diffRate_{rate_measure}'],
                    color=palette[i], linewidth=1)

        ax.set_ylabel('Max − min rate (Hz)')
        ax.set_title(f'E–I balance = {b}')

    axes[-1].set_xlabel('Global coupling strength')

    handles = [plt.Line2D([0], [0], color=palette[i], linewidth=2) for i in range(len(Ib_str))]
    fig.legend(handles, [str(s) for s in Ib_str], title='Background input',
               loc='center left', bbox_to_anchor=(1.0, 0.5))

    sns.despine(trim=True)
    plt.tight_layout()
    figure_name = (f'multisI_couplingOnDiffRate_{rate_measure}_pop{population}_'
                   f'Ibstr{Ib_str}_sI{sI}_tauVisual_{thalamus_source}')
    plt.savefig(os.path.join(figure_dir, figure_name + '.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(figure_dir, figure_name + '.png'), bbox_inches='tight', dpi=300)
    plt.show()

def multiLayer_couplingOnMinMaxRate(data_df, Iext_dur, Iext_str, Ib_str, sI, rate_measure, thalamus_source, figure_dir):
    """
    Effect of Coupling Strengths on Min and Max Firing rates.
    Parameters:
    -----------
    - rate_measure: measure to plot (e.g. 'lateLongterm' or 'immediateLongterm', 'duringInput', 'baseline')
    """
    # separate data in layers
    layersS1 = [['E1', 'PV1', 'SST1', 'VIP1'], ['E2', 'PV2', 'SST2'], ['E3', 'PV3', 'SST3'], ['E4', 'PV4', 'SST4']]
    layersS2 = [['E1S2', 'PV1S2', 'SST1S2', 'VIP1S2'], ['E2S2', 'PV2S2', 'SST2S2'], ['E3S2', 'PV3S2', 'SST3S2'], ['E4S2', 'PV4S2', 'SST4S2']]
    layers_all = [layersS1, layersS2]

    # plot results
    fig, axes = plt.subplots(4, 3, figsize=(15, 10), sharey=False, sharex=True)
    fig.suptitle(f'Effect of Coupling strength on Min and Max firing rate {rate_measure}')

    # columns: S1, S2, A3b
    for col_idx, (layers, axs) in enumerate(zip(layers_all, axes.T[:2])):
        for l, ax in zip(layers, axs):
            layer_df = data_df[data_df['population'].isin(l)]
            sns.lineplot(layer_df, y=f'minRate_{rate_measure}', x='globalCoupling', hue='population', ax=ax)
            sns.lineplot(layer_df, y=f'maxRate_{rate_measure}', x='globalCoupling', hue='population', ax=ax, legend=False)
            ax.set_ylabel('Rate (Hz)')
            ax.set_xlabel('Global Coupling Strength')
            ax.legend(prop={'size':8})

        axs[0].set_title('Layer 2/3' if col_idx == 0 else 'Layer 2/3 (S2)')
        axs[1].set_title('Layer 4' if col_idx == 0 else 'Layer 4 (S2)')
        axs[2].set_title('Layer 5' if col_idx == 0 else 'Layer 5 (S2)')
        axs[3].set_title('Layer 6' if col_idx == 0 else 'Layer 6 (S2)')

    # A3b column
    a3b_cells = ['E3b', 'PV3b', 'SST3b', 'VIP3b']
    a3b_axs = axes[:, 2]
    for row_idx, ax in enumerate(a3b_axs):
        if row_idx == 0:
            for cell in a3b_cells:
                layer_df = data_df[data_df['population'] == cell]
                sns.lineplot(layer_df, y=f'minRate_{rate_measure}', x='globalCoupling', ax=ax, legend=True)
                sns.lineplot(layer_df, y=f'maxRate_{rate_measure}', x='globalCoupling', ax=ax, legend=False)
            ax.set_ylabel('Rate (Hz)')
            ax.set_xlabel('Global Coupling Strength')
            ax.set_title('A3b (Layer 2/3)')
        else:
            ax.axis('off')

    sns.despine(trim=True)
    plt.annotate(
        f'min and max Rate {rate_measure}, \nInput Duration:{Iext_dur} \nInput Strength:{Iext_str} '
        f'Background Input:{Ib_str} \nE-I Balance:{sI}',
        xy=(0, 0), xycoords='figure fraction', textcoords='figure fraction', xytext=(1, 0.2), ha='center', fontsize=10
    )
    plt.tight_layout()
    figure_name = f'multiLayer_couplingOnMinMaxRate_{rate_measure}_Iextstr{Iext_str}_Ibstr{Ib_str}_sI{sI}_tauVisual_{thalamus_source}.pdf'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
    plt.show()


def inputStrengthOnminMaxpotential(data_df, Iext_str, Ib_str, sI, population, potential_measure, thalamus_source, figure_dir):
    """
    Plot the effect of input strength on min and max potential 
    
    Parameters:
    -----------
    - rate_measure: measure to plot (e.g. 'lateLongterm' or 'immediateLongterm', 'duringInput', 'baseline')
    """
    fig, ax = plt.subplots(figsize=(8,6))
    fig.suptitle('Fixed point: steady state potential by input strength\npopulation: '+population)

    data_df = data_df[data_df['population']==population]
    sns.lineplot(data_df, x='InputStrength', y=f'minPotential_{potential_measure}', hue='globalCoupling', ax=ax)
    sns.lineplot(data_df, x='InputStrength', y=f'maxPotential_{potential_measure}', hue='globalCoupling', ax=ax, legend=False)
    sns.despine(trim=True)
    ax.set_ylabel('Mean Potential (mV)')
    ax.set_xlabel('External Input Strength')
    plt.savefig(os.path.join(figure_dir, f'FixedPoint_potential_{potential_measure}_Iextstr{Iext_str}_Ibstr{Ib_str}_sI{sI}_pop{population}_tauVisual_{thalamus_source}.pdf'), bbox_inches='tight')
    plt.show()

def plot_axis(axs, steps, start_plot, rates, idx_rates):
    axs.plot(steps[start_plot:], rates[idx_rates].T[start_plot:], linewidth=1)

def plot_population_rates(axs_op, idxs_pop, rates, steps, start_plot, labels):
    """ Plot population rates for given labels and indices."""
    legend_list = []
    for i, idx in enumerate(idxs_pop):
        plot_axis(axs_op, steps, start_plot, rates, idx)
        legend_list.append(f'{labels[i]} {np.round(rates[idx].T[-1], 6)}')

    axs_op.legend(legend_list, loc='upper right')

def plot_rate_results(rates, Iext, Ib, step_size, simulation_time, start_plot, sI, g, d, sb, s, figure_dir):
    steps = np.arange(step_size, simulation_time+step_size, step_size)*1e3
    fig, axs = plt.subplots(4, 3, figsize=(15, 15))  # Set figure size
    figure_style()

    # external input 
    axs_extI = axs[0][0]
    axs_extI.plot(steps[start_plot:], Iext[start_plot:], label='Iext rate')
    axs_extI.plot(steps[start_plot:], Ib[start_plot:], label='Ib rate')
    axs_extI.legend(title='')
    axs_extI.set_ylabel('Hz')
    # thalamus
    axs_thal = axs[1][0]
    axs_thal.plot(steps[start_plot:], rates[-2:-1].T[start_plot:], color='purple')
    axs_thal.plot(steps[start_plot:], rates[-1:].T[start_plot:], color='grey')
    axs_thal.legend(['Thalamus E', 'Thalamus I'])
    axs_thal.set_ylabel('Hz')

    # area 3b
    axsA3b = axs[2][0]
    axsA3b.plot(steps[start_plot:], rates[:4].T[start_plot:], linewidth=1)
    axsA3b.legend([f'E {np.round(rates[0].T[-1], 6)}', f'PV {np.round(rates[1].T[-1], 6)}', f'SOM {np.round(rates[2].T[-1], 6)}', f'VIP {np.round(rates[3].T[-1], 6)}'])
    axsA3b.set_ylabel('Hz')

    # plot results for the S1 column 
    idxs_E = np.array([0+4, 4+4, 7+4, 10+4]) # indices of E populations in S1
    labels_pops = [['E1', 'E2', 'E3', 'E4'], ['PV1', 'PV2', 'PV3', 'PV4'], ['SST1', 'SST2', 'SST3', 'SST4'], ['VIP1']]
    
    # loop over populations for S1 and S2
    for i, labels in enumerate(labels_pops):
        if i<3:
            axs_pop = axs[i][1]
            plot_population_rates(axs_pop, idxs_E+i, rates, steps, start_plot, labels)
        else:
            # VIP
            axsVIPS1 = axs[i][1]
            axsVIPS1.plot(steps[start_plot:], rates[3+4].T[start_plot:], linewidth=1)
            axsVIPS1.legend([f'VIP1 {np.round(rates[3+4].T[-1], 6)}'])

        # plot results S2
        nr_pops = 13 # number of pops in S1
        if i<3:
            axs_pop = axs[i][2]
            plot_population_rates(axs_pop, idxs_E+i+nr_pops, rates, steps, start_plot, labels)
        else:
            # VIP
            axsVIPS2 = axs[i][2]
            axsVIPS2.plot(steps[start_plot:], rates[3+nr_pops].T[start_plot:], linewidth=1)
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

    #annotate_fig(f'sI={np.round(sI, 4)}, g={np.round(g, 4)}, area={area}')
    sns.despine(trim=True)
    plt.tight_layout() 
    plt.legend()
    figdir = os.path.join(figure_dir, 'single_simulations')
    if not os.path.exists(figdir):
        os.makedirs(figdir)

    #plt.savefig(os.path.join(figdir, f'population_rates_sI-{sI}_g-{g}_area-{area}_Iextdur-{d}_Iextstr-{s}_Ibstr-{sb}.pdf'), dpi=300)
    #plt.show()

def plot_potentials(potentials, Iext, Ib, step_size, simulation_time, start_plot, figdir, sI, g, d, sb, s, area='all'):
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

        # Layout: 4 rows (max 4 layers), 3 columns (3 areas)
        fig, axes = plt.subplots(4, 3, figsize=(14, 10), sharex=True, sharey=False)
        axes = np.array(axes)

        # --- Column 1: Area 3b + Thalamus stacked ---
        axes[0, 0].plot(potentials[:4].T)
        axes[0, 0].legend(['E', 'PV', 'SOM', 'VIP'])
        axes[0, 0].legend([f'E {np.round(potentials[0, -1], 6)}', f'PV {np.round(potentials[1, -1], 6)}', f'SOM {np.round(potentials[2, -1], 6)}', f'VIP {np.round(potentials[3, -1], 6)}'], loc='upper right')
        axes[0, 0].set_title("Area 3b")
        axes[0, 0].set_ylabel('mV')

        axes[1, 0].plot(potentials[30:].T)
        axes[1, 0].set_title("Thalamus")

        # Hide extra rows in col 1 (since only 2 plots)
        for r in range(2, 4):
            axes[r, 0].axis("off")

        # --- Column 2: Area 1 layers ---
        area_1_layers = [[4,5,6,7],[8,9,10],[11,12,13],[14,15,16]]
        pop_names = ['E', 'PV', 'SOM', 'VIP']
        for i, layer_idx in enumerate(area_1_layers):
            axes[i, 1].plot(potentials[layer_idx].T)
            axes[i, 1].set_title(f"Area 1 - Layer {i+1}")
            axes[i, 1].set_ylabel('mV')
            if len(layer_idx)==4:
                axes[i, 1].legend([f'E {np.round(potentials[layer_idx[0], -1], 6)}', f'PV {np.round(potentials[layer_idx[1], -1], 6)}', f'SOM {np.round(potentials[layer_idx[2], -1], 6)}', f'VIP {np.round(potentials[layer_idx[3], -1], 6)}'], loc='upper right')
            else:
                axes[i, 1].legend([f'E {np.round(potentials[layer_idx[0], -1], 6)}', f'PV {np.round(potentials[layer_idx[1], -1], 6)}', f'SOM {np.round(potentials[layer_idx[2], -1], 6)}'], loc='upper right')


        # --- Column 3: Area S2 layers ---
        area_s2_layers = [[17,18,19,20],[21,22,23],[24,25,26],[27,28,29]]
        for i, layer_idx in enumerate(area_s2_layers):
            axes[i, 2].plot(potentials[layer_idx].T)
            axes[i, 2].set_ylabel('mV')
            axes[i, 2].set_title(f"Area S2 - Layer {i+1}")
            if len(layer_idx)==4:
                axes[i, 2].legend([f'E {np.round(potentials[layer_idx[0], -1], 6)}', f'PV {np.round(potentials[layer_idx[1], -1], 6)}', f'SOM {np.round(potentials[layer_idx[2], -1], 6)}', f'VIP {np.round(potentials[layer_idx[3], -1], 6)}'], loc='upper right')
            else:
                axes[i, 2].legend([f'E {np.round(potentials[layer_idx[0], -1], 6)}', f'PV {np.round(potentials[layer_idx[1], -1], 6)}', f'SOM {np.round(potentials[layer_idx[2], -1], 6)}'], loc='upper right')


        plt.tight_layout()
        figdir = os.path.join(figdir, 'single_simulations')
        if not os.path.exists(figdir):
            os.makedirs(figdir)
        plt.savefig(os.path.join(figdir, f'population_potentials_sI-{sI}_g-{g}_area-{area}_Iextdur-{d}_Iextstr-{s}_Ibstr-{sb}.pdf'), dpi=300)
        #plt.show()


def plot_response_type_examples(examples, fixed_params, population, raw_dir, figure_dir,
                                suffix='', plot_window=(-0.15, 0.25)):
    """
    Publication/poster-ready figure: example firing-rate trajectories showing how a single
    population reacts to the *same* brief stimulus under several parameter settings. Each
    setting is one vertically stacked panel (shared time axis, independent y-axes) so that
    qualitatively different response archetypes stay legible despite very different rate scales.

    Parameters
    ----------
    examples : list of dict
        Ordered list, one entry per panel, each with keys:
            'label' : str   - archetype name shown in bold (e.g. 'Memory')
            'g'     : float - global coupling strength of that run
            'sI'    : float - E-I balance of that run
            'color' : color - line colour for that panel
    fixed_params : dict
        Stimulus/simulation parameters shared by every run:
            'g_inter', 'Ib_str', 'Iext_dur', 'Iext_str', 'input_onset',
            'thal_cellcounts', 'bI_cellcounts', 'extI_cellcounts', 'input_type', 'step_size'
    population : str
        Population column to plot (e.g. 'E3').
    raw_dir : str
        Directory holding the per-run result folders.
    figure_dir : str
        Where to save the figure.
    plot_window : (float, float)
        Time window (s) relative to stimulus onset to display.
    """
    figure_style()

    fp = fixed_params
    step_size = fp['step_size']
    input_onset = fp['input_onset']
    Iext_dur = fp['Iext_dur']

    n = len(examples)
    fig, axes = plt.subplots(n, 1, figsize=(3.6, 1.7 * n), sharex=True, sharey=False)
    axes = np.atleast_1d(axes)

    for ax, ex in zip(axes, examples):
        g, sI = ex['g'], ex['sI']
        rates_df, _, _ = load_simulation_data(
            g, fp['g_inter'], sI, fp['Ib_str'], Iext_dur, fp['Iext_str'], input_onset,
            fp['thal_cellcounts'], fp['bI_cellcounts'], fp['extI_cellcounts'],
            fp['input_type'], raw_dir, suffix=suffix,
            g_thalPOm=fp.get('g_thalPOm', 1.0), Ib_noise_std=fp.get('Ib_noise_std'))

        rate = rates_df[population].values
        # time re-zeroed to stimulus onset
        t = np.arange(len(rate)) * step_size - input_onset
        mask = (t >= plot_window[0]) & (t <= plot_window[1])

        # stimulus window (0 .. Iext_dur) shaded in every panel
        ax.axvspan(0, Iext_dur, color='0.85', linewidth=0, zorder=0)
        ax.plot(t[mask], rate[mask], color=ex['color'], linewidth=1.2, zorder=2)

        # archetype label (bold) + parameter subtitle (grey) inside the panel
        ax.text(0.97, 0.9, ex['label'], transform=ax.transAxes, ha='right', va='top',
                fontweight='bold', fontsize=9, color=ex['color'])
        ax.text(0.97, 0.74, f'$g$ = {g},  $s_I$ = {sI}', transform=ax.transAxes, ha='right',
                va='top', fontsize=7.5, color='0.4')

    # stimulus marker + label on the top panel only
    top = axes[0]
    top.annotate('stimulus', xy=(Iext_dur / 2, 1.0), xycoords=('data', 'axes fraction'),
                 xytext=(0, 6), textcoords='offset points', ha='center', va='bottom',
                 fontsize=7.5, color='0.2',
                 arrowprops=dict(arrowstyle='-|>', color='0.2', lw=0.8))

    axes[-1].set_xlabel('Time from stimulus onset (s)')
    fig.text(0.005, 0.5, 'Firing rate (Hz)', va='center', rotation='vertical')

    sns.despine(trim=True)
    plt.tight_layout(rect=[0.03, 0, 1, 1])

    figure_name = f'responseType_examples_pop{population}{suffix}'
    plt.savefig(os.path.join(figure_dir, figure_name + '.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(figure_dir, figure_name + '.png'), bbox_inches='tight', dpi=300)
    plt.show()


def plot_gthalPOm_potential_sweep(g_thalPOms, fixed_params, populations, raw_dir, figure_dir,
                                  suffix='', plot_window=(-0.1, 0.95)):
    """
    Publication/poster-ready figure: how each population's summed potential responds to the SAME
    stimulus as the POm output-gain g_thalPOm is varied. One line per g_thalPOm value, each in a
    distinct colour with a shared legend. Near g_thalPOm ~ 1 the post-stimulus attractor is highly
    sensitive, so a few hand-picked values span very different response regimes.

    Parameters
    ----------
    g_thalPOms : list of float
        The g_thalPOm values to overlay (kept small so each colour is distinct).
    fixed_params : dict
        Operating point shared by every run:
            'g', 'g_inter', 'sI', 'Ib_str', 'Iext_dur', 'Iext_str', 'input_onset',
            'thal_cellcounts', 'bI_cellcounts', 'extI_cellcounts', 'input_type', 'step_size'
    populations : str or 2D list of str
        Potential column(s) to plot. A 2D list lays the populations out as a grid of panels
        (rows x cols), e.g. [['E2', 'E2S2'], ['E3', 'E3S2']]. A plain string plots a single panel.
    raw_dir : str
        Directory holding the per-run result folders.
    figure_dir : str
        Where to save the figure.
    plot_window : (float, float)
        Time window (s) relative to stimulus onset to display.
    """
    figure_style()

    fp = fixed_params
    step_size = fp['step_size']
    input_onset = fp['input_onset']
    Iext_dur = fp['Iext_dur']

    # normalise `populations` to a 2D grid
    if isinstance(populations, str):
        pop_grid = [[populations]]
    else:
        pop_grid = populations
    nrows = len(pop_grid)
    ncols = max(len(r) for r in pop_grid)

    # distinct, colourblind-friendly qualitative palette (one colour per value)
    base_colors = ['#4477AA', '#EE6677', '#228833', '#AA3377', '#CCBB44',
                   '#66CCEE', '#EE7733', '#000000']
    colors = (base_colors * (len(g_thalPOms) // len(base_colors) + 1))[:len(g_thalPOms)]

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.0 * nrows),
                             sharex=True, sharey=False)
    axes = np.atleast_2d(axes)

    # stimulus window shaded in every panel
    for ax in axes.flat:
        ax.axvspan(0, Iext_dur, color='0.85', linewidth=0, zorder=0)

    # load each run once, then plot its trace into every population panel
    for gp, col in zip(g_thalPOms, colors):
        try:
            _, potentials_df, _ = load_simulation_data(
                fp['g'], fp['g_inter'], fp['sI'], fp['Ib_str'], fp['Iext_dur'], fp['Iext_str'],
                input_onset, fp['thal_cellcounts'], fp['bI_cellcounts'], fp['extI_cellcounts'],
                fp['input_type'], raw_dir, suffix=suffix, g_thalPOm=np.round(gp, 3),
                Ib_noise_std=fp.get('Ib_noise_std'))
        except (FileNotFoundError, OSError, KeyError):
            print(f'Missing run for g_thalPOm={gp} - skipping')
            continue

        for r, row in enumerate(pop_grid):
            for c, pop in enumerate(row):
                pot = potentials_df[pop].values
                t = np.arange(len(pot)) * step_size - input_onset
                mask = (t >= plot_window[0]) & (t <= plot_window[1])
                axes[r, c].plot(t[mask], pot[mask], color=col, linewidth=1.3,
                                label=f'{gp:g}', zorder=2)

    # map population id -> two-line "area / layer" label (e.g. E2 -> 'Area 1\nLayer 4')
    pop_labels = {'E1': 'Area 1\nLayer 2/3', 'E2': 'Area 1\nLayer 4',
                  'E3': 'Area 1\nLayer 5', 'E4': 'Area 1\nLayer 6',
                  'E1S2': 'S2\nLayer 2/3', 'E2S2': 'S2\nLayer 4',
                  'E3S2': 'S2\nLayer 5', 'E4S2': 'S2\nLayer 6',
                  'E3b': 'A3b\nLayer 2/3'}

    # per-panel titles + shared axis labels (bottom row / left column only)
    for r, row in enumerate(pop_grid):
        for c, pop in enumerate(row):
            ax = axes[r, c]
            ax.set_title(pop_labels.get(pop, pop), fontweight='bold')
            if r == nrows - 1:
                ax.set_xlabel('Time from stimulus onset (s)')
            if c == 0:
                ax.set_ylabel('Summed potential (mV)')

    # single shared legend built from the g_thalPOm / colour pairs
    handles = [plt.Line2D([0], [0], color=col, linewidth=2) for col in colors]
    labels = [f'{gp:g}' for gp in g_thalPOms]
    fig.legend(handles, labels, title=r'$g_{\mathrm{thalPOm}}$', frameon=False,
               fontsize=8, title_fontsize=9, loc='center left', bbox_to_anchor=(1.0, 0.5))

    sns.despine(trim=True)
    plt.tight_layout()

    flat_pops = [p for row in pop_grid for p in row]
    figure_name = f'gthalPOm_sweep_potential_' + '_'.join(flat_pops) + suffix
    plt.savefig(os.path.join(figure_dir, figure_name + '.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(figure_dir, figure_name + '.png'), bbox_inches='tight', dpi=300)
    plt.show()


def plot_gthalPOm_dipole_sweep(g_thalPOms, sim_overrides, figure_dir, subjects=(15,),
                               area='S2', suffix='', plot_window=(-0.1, 0.25)):
    """
    Companion to plot_gthalPOm_potential_sweep: the SUMMED computed dipole of an area (default S2)
    for each g_thalPOm value. The dipole is a source-resolved weighting of the full 3-D potential
    (SomatoModel.compute_dipoles), which is not saved in results.hdf5, so each parameter set is
    re-simulated here and its dipole computed with a subject forward model.

    Parameters
    ----------
    g_thalPOms : list of float
        The g_thalPOm values to overlay (kept small so each colour is distinct).
    sim_overrides : dict
        SomatoModel parameter overrides for the operating point (model param names), e.g.
        {'coupling_strength', 'strength_I', 'Iext_duration', 'Iext_strength', 'Ib_strength',
         'g_intercortical', 'g_thal', 'sI_thal', 'area', 'resistance_factor', 'delay_factor',
         'extI_cellcounts', 'bI_cellcounts'}. g_thalPOm is set per value.
    figure_dir : str
        Where to save the figure.
    subjects : sequence of int
        Subject IDs whose forward models project the dipole (averaged); default (15,).
    area : str
        'A3b', 'A1' or 'S2' - which area's summed dipole to plot.
    plot_window : (float, float)
        Time window (s) relative to stimulus onset to display.
    """
    # lazy import so the analysis module isn't coupled to the model (and its env vars) at import
    import sys
    _wd = os.getenv('WDDIR')
    for _p in [os.path.join(_wd, 'Simulations'), os.path.join(_wd, 'Simulations', 'model')]:
        if _p not in sys.path:
            sys.path.insert(0, _p)
    from somato_model import SomatoModel, read_simulation_params

    figure_style()

    # dipole-row layout returned by compute_dipoles: 0=A3b, 1-4=A1 L{2/3,4,5,6}, 5-8=S2 L{2/3,4,5,6}
    area_rows = {'A3b': [0], 'A1': [1, 2, 3, 4], 'S2': [5, 6, 7, 8]}
    if area not in area_rows:
        raise ValueError(f"area must be one of {list(area_rows)}, got {area!r}")
    rows = area_rows[area]

    # distinct, colourblind-friendly qualitative palette (one colour per value)
    base_colors = ['#4477AA', '#EE6677', '#228833', '#AA3377', '#CCBB44',
                   '#66CCEE', '#EE7733', '#000000']
    colors = (base_colors * (len(g_thalPOms) // len(base_colors) + 1))[:len(g_thalPOms)]

    # one model, reused across values (compute_dipoles caches the forward projection after the
    # first call, so the forward model is read only once) - mirrors run_optimization.py
    base_params = read_simulation_params()
    base_params.update(sim_overrides)
    model = SomatoModel(base_params)
    step_size = model.step_size
    input_onset = model.input_onset
    Iext_dur = model.Iext_duration

    fig, ax = plt.subplots(figsize=(7.4, 3.6))
    ax.axvspan(0, Iext_dur, color='0.85', linewidth=0, zorder=0)

    for gp, col in zip(g_thalPOms, colors):
        model.apply_params({**sim_overrides, 'g_thalPOm': np.round(gp, 3)})
        model.initialize_state()
        model.simulate()
        simDipoles = model.compute_dipoles(list(subjects))
        dip = simDipoles[rows].sum(axis=0)

        t = np.arange(len(dip)) * step_size - input_onset
        mask = (t >= plot_window[0]) & (t <= plot_window[1])
        ax.plot(t[mask], dip[mask], color=col, linewidth=1.3, label=f'{gp:g}', zorder=2)

    ax.set_xlabel('Time from stimulus onset (s)')
    ax.set_ylabel(f'{area} summed dipole (a.u.)')
    # g_thalPOm legend outside, top-right
    ax.legend(title=r'$g_{\mathrm{thalPOm}}$', frameon=False, fontsize=8,
              title_fontsize=9, loc='upper left', bbox_to_anchor=(1.02, 1.0))

    # box listing the physically-meaningful fixed parameters, in mathematical notation
    # (key, description, math symbol, optional value transform for display)
    _box = [
        ('coupling_strength', 'Intra-cortical coupling', r'g', None),
        ('strength_I',        'Inhibitory coupling',     r's_I', None),
        ('g_intercortical',   'Inter-cortical coupling', r'g_\mathrm{inter}', None),
        ('Ib_strength',       'Background input',        r'I_\mathrm{b}', None),
        ('Iext_strength',     'External input',          r'I_\mathrm{ext}', None),
        ('Iext_duration',     'Input duration',          r'\Delta t', 'ms'),
    ]
    lines = ['Fixed parameters']
    for key, desc, sym, unit in _box:
        if key not in sim_overrides:
            continue
        val = sim_overrides[key]
        if unit == 'ms':
            lines.append(f'{desc}  ${sym} = {val * 1e3:g}$ ms')
        else:
            lines.append(f'{desc}  ${sym} = {val:g}$')
    ax.text(1.02, 0.62, '\n'.join(lines), transform=ax.transAxes, va='top', ha='left',
            fontsize=7.5,
            bbox=dict(boxstyle='round', facecolor='0.96', edgecolor='0.7', linewidth=0.6))

    sns.despine(trim=True)
    plt.tight_layout()

    figure_name = f'gthalPOm_sweep_dipole_{area}{suffix}'
    plt.savefig(os.path.join(figure_dir, figure_name + '.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(figure_dir, figure_name + '.png'), bbox_inches='tight', dpi=300)
    plt.show()
