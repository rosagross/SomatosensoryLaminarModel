import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.cm as cm 
from matplotlib.colors import ListedColormap, Normalize, BoundaryNorm
from plotting_style import figure_style

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

    fig, axes = plt.subplots(4, 2, figsize=(10, 6), sharey=False, sharex=True)
    fig.suptitle(f"Dominant {signal_type} frequency during {window_prefix}")

    for layers, axs in zip(layers_all, axes.T):
        for l, ax in zip(layers, axs):
            layer_df = data_df[data_df["population"].isin(l)]
            sns.lineplot(layer_df, y=col, x="globalCoupling", hue="population", ax=ax, legend=True)
            ax.set_ylabel("Frequency (Hz)")
            ax.set_xlabel("Global Coupling Strength (g)")
            ax.legend("")

        axs[0].legend(prop={"size":8}, loc="upper right")
        axs[0].set_title("Layer 2/3")
        axs[1].set_title("Layer 4")
        axs[2].set_title("Layer 5")
        axs[3].set_title("Layer 6")

    sns.despine(trim=True)
    plt.tight_layout()
    plt.annotate(
        f"freq {signal_type}_{window_prefix}, \nInput Duration:{Iext_dur} \nInput Strength:{Iext_str} "
        f"Background Input:{Ib_str} \nE-I Balance:{sI}",
        xy=(0, 0), xycoords="figure fraction", xytext=(1.14, 0.17),
        textcoords="figure fraction", ha="center", fontsize=10
    )
    figure_name = (
        f"CouplingOnFrequency_{signal_type}_{window_prefix}_Iextdur{Iext_dur}_"
        f"Iextstr{Iext_str}_Ibstr{Ib_str}_sI{sI}_tauVisual_{thalamus_source}.pdf"
    )
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches="tight")
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
    Heatmaps of dominant frequency vs coupling strength (x) and sI (y) for all populations in an area.
    area: "A1", "S2", or "A3b"
    """
    col = f"freq{signal_type}_{window_prefix}"
    if col not in data_df.columns:
        raise ValueError(f"Column '{col}' not found in data frame.")

    if area == "A1":
        populations = ['E1', 'PV1', 'SST1', 'VIP1',
                       'E2', 'PV2', 'SST2',
                       'E3', 'PV3', 'SST3',
                       'E4', 'PV4', 'SST4']
    elif area == "S2":
        populations = ['E1S2', 'PV1S2', 'SST1S2', 'VIP1S2',
                       'E2S2', 'PV2S2', 'SST2S2',
                       'E3S2', 'PV3S2', 'SST3S2',
                       'E4S2', 'PV4S2', 'SST4S2']
    elif area == "A3b":
        populations = ['E3b', 'PV3b', 'SST3b', 'VIP3b']
    else:
        raise ValueError(f"Unsupported area: {area}. Use 'A1', 'S2', or 'A3b'.")

    plot_df = data_df[
        (data_df["InputDuration"] == Iext_dur) &
        (data_df["InputStrength"] == Iext_str) &
        (data_df["BckgndInputStrength"] == Ib_str) &
        (data_df["population"].isin(populations))
    ].copy()

    n = len(populations)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), sharex=True, sharey=True)
    axes = np.atleast_1d(axes).flatten()

    # for the colorbar
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap="magma", norm=norm)
    sm.set_array([])

    for ax, pop in zip(axes, populations):
        pop_df = plot_df[plot_df["population"] == pop]
        data_heatmap = pop_df.pivot(index="strength_I", columns="globalCoupling", values=col)

        # seaborn.heatmap fails on empty matrices ("zero-size array to reduction")
        # which can happen when no rows match the current population/filter values.
        if data_heatmap.empty or data_heatmap.shape[0] == 0 or data_heatmap.shape[1] == 0 or np.all(np.isnan(data_heatmap.to_numpy())):
            ax.axis("off")
            ax.text(0.5, 0.5, f"{pop}\n(no data)", ha="center", va="center", transform=ax.transAxes)
            continue

        sns.heatmap(data_heatmap, cmap="magma", vmin=vmin, vmax=vmax, ax=ax, cbar=False)
        ax.invert_yaxis()
        ax.set_title(pop)
        ax.set_xlabel("g")
        ax.set_ylabel("sI")

    # turn off unused axes
    for ax in axes[n:]:
        ax.axis('off')

    cbar = fig.colorbar(sm, ax=axes, orientation="vertical", fraction=0.02, pad=0.02)
    cbar.set_label(col)

    fig.suptitle(f"{signal_type} freq ({window_prefix}) - area {area}")
    plt.tight_layout()
    figure_name = (
        f"freqHeatmap_gVsSI_{signal_type}_{window_prefix}_Iextdur{Iext_dur}_"
        f"Iextstr{Iext_str}_Ib{Ib_str}_area{area}.pdf"
    )
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches="tight")
    plt.show()

def multisI_couplingOnMinmaxRate(data_df, sI, Ib_str, population, rate_measure, thalamus_source, figure_dir):
    '''
    3.1) Background input strength 
    - x-axis: input strength
    - y-axis: PSP of E3 population
    - one subplot for each sI value
    - choose rate measure ('lateLongterm' or 'immediateLongterm', 'duringInput', 'baseline')
    '''
    data_df = data_df[data_df['population']==population]

    # plot results
    f_height = len(sI) * 2
    fig, axes = plt.subplots(len(sI), 1, figsize=(8,f_height)) 
    fig.suptitle(f'Minimum and Maximum rate of {population} Population')
    for axs, b in zip(axes, sI):
        
        data_sI_df = data_df[data_df['strength_I']==b]
        
        # Create a colormap
        cmap = cm.get_cmap('Dark2', len(Ib_str))  # Choose a colormap, e.g., 'viridis'
        cmap_max = cm.get_cmap('Dark2', len(Ib_str))  # Choose a colormap, e.g., 'viridis'

        for i, s in enumerate(Ib_str):
            Istrength_df = data_sI_df[data_sI_df['BckgndInputStrength']==s]
            color = cmap(i / (len(Ib_str)))  # Normalize i to [0, 1] for colormap
            color_max = cmap_max(i / (len(Ib_str)))  # Normalize i to [0, 1] for colormap
            axs.plot(Istrength_df['globalCoupling'], Istrength_df[f'minRate_{rate_measure}'], label=s, color=color)
            axs.plot(Istrength_df['globalCoupling'], Istrength_df[f'maxRate_{rate_measure}'], color=color_max)

        #sns.lineplot(data_sI_df, y='maxPotential_longterm', x='coupling_strength', hue='InputStrength')
        axs.set_xlabel('Coupling Strength')
        axs.set_ylabel('Rate (Hz)')
        axs.set_title(f'E-I Balance:{b}')

    #axs.set_xlim([0,100])
    sns.despine(trim=True)
    plt.legend(title=f'Effect of background input and coupling strength on MinMax Rate {rate_measure}', loc='right')
    plt.tight_layout() 
    figure_name = f'multisI_couplingOnMinmaxRate_{rate_measure}_pop{population}_measure{rate_measure}_Ibstr{Ib_str}_sI{sI}_tauVisual_{thalamus_source}.pdf'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
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

def plot_results(rates, Iext, Ib, step_size, simulation_time, start_plot, sI, g, d, sb, s, figure_dir):
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
            
