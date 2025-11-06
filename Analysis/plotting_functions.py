import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.cm as cm 
from matplotlib.colors import ListedColormap, Normalize, BoundaryNorm

def bEImulti_fingerprint_IextDurVsStr(data_df, g, bEIs, Ib_str, population, thalamus_source, figure_dir):
    '''
    1.3) MULTI fingerprint PLOT: Effect of input intensity and duration on dynamic behaviour (non-responsive, transfer and memory)
        - plot style: heatmap
        - y axis: intensity
        - x axis: duration
        - plotwise: coupling strength
        - measure: dynamic function (the "finger print")

    Parameters:
    -----------
    - g: list of global coupling strengths
    '''

    cbar_ticks = ['non-responsive', 'transfer', 'memory']
    data_df = data_df[data_df['population']==population]
    data_df = data_df[data_df['globalCoupling']==g]
    data_df = data_df[data_df['BckgndInputStrength']==Ib_str]
    data_df['InputDuration'] = data_df['InputDuration'].round(4)
    fig, ax = plt.subplots(2, int(len(bEIs)/2), figsize=(9,4), sharex=True, sharey=True)

    # Define a fixed discrete colormap with colors for 1, 2, and 3
    n = 3 # there are three different functions of the dynamics --> make discrete colormap
    colors = sns.color_palette("Pastel2", n)
    cmap = ListedColormap(colors)

    # Define the boundaries for normalization
    bounds = [0.5, 1.5, 2.5, 3.5]
    norm = BoundaryNorm(bounds, ncolors=cmap.N)

    for i, axis in enumerate(ax.flatten()):
        plot_df = data_df[data_df['balanceEI']==bEIs[i]]
        data_heatmap = plot_df.pivot(index='InputStrength',columns='InputDuration', values='dynamic_function_potential')
        if i == len(bEIs)-1:
            heat_ax = sns.heatmap(data_heatmap, cmap=cmap, norm=norm, ax=axis, cbar=True)
            # Add colorbar, make sure to specify tick locations to match desired ticklabels
            cbar = heat_ax.collections[0].colorbar
            cbar.set_ticks([1, 2, 3])
            cbar.set_ticklabels(cbar_ticks)
        else:
            sns.heatmap(data_heatmap, cmap=cmap, norm=norm, cbar=False, ax=axis, vmin=1, vmax=3, cbar_kws={'ticks': cbar_ticks})

        # Only set x-axis label for the bottom row
        if i >= len(bEIs) // 2:
            axis.set_xlabel('Input Duration')
        else:
            axis.set_xlabel('')
        axis.set_ylabel('')
        axis.invert_yaxis()
        axis.set_title(f'bEI = {bEIs[i]}, g = {g}')

    #ax[0][0].set_ylabel('Input Strength')
    #ax[1][0].set_ylabel('Input Strength')

    plt.tight_layout(h_pad=1)
    figure_name = f'_{population}pop_tauVisual_{thalamus_source}.pdf'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
    plt.show()


def multiPop_heatmap_IextDurVsStr(data_df, gs, bEI, Ib_str, populations, rate_measure, figure_dir):
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
    - bEI: float balance E/I value
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
            minmax_df = minmax_df[minmax_df['balanceEI']==bEI]
            minmax_df = minmax_df[minmax_df['BckgndInputStrength']==Ib_str]
            minmax_df = minmax_df[minmax_df['population']==p]
            minmax_df['InputDuration'] = minmax_df['InputDuration'].round(4)
            minmax_df[rate_measure] = minmax_df[rate_measure].round(5)

            data_heatmap = minmax_df.pivot(index='InputStrength',columns='InputDuration', values=rate_measure)

            if (minmax_df[rate_measure].isna() | (minmax_df[rate_measure] == 0)).all().all():
                print(f'No measure value for {rate_measure} at g: {g}, bEI: {bEI}, pop: {p}')
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
            axes[i,0].set_ylabel(f'g: {g}, bEI: {bEI}', rotation=0, labelpad=60)

    fig.text(0.05, 0.2, 'Input Strength', va='center', rotation='vertical')
    fig.text(0.05, 0.5, 'Input Strength', va='center', rotation='vertical')
    fig.text(0.05, 0.83, 'Input Strength', va='center', rotation='vertical')

    plt.tight_layout(h_pad=15)
    figure_name = f'inputDurationVSinputStrength_{populations[0][0]}pop_bEI{bEI}_Ibstr{Ib_str}_{rate_measure}_tauVisual.png'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
    plt.show()


def multiLayer_couplingOnLongeterm_diffRate(data_df, Iext_dur, Iext_str, Ib_str, bEI, thalamus_source, figure_dir):
    """
    Effect of Coupling Strengths on Longterm/steady state
    Plot difference between Mininum and Maximum Firing rates 
    """

    data_df = data_df[data_df['InputDuration']==Iext_dur]
    data_df = data_df[data_df['balanceEI']==bEI]
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
    plt.annotate(f'min and max Rate_lateLongterm, \nInput Duration:{Iext_dur} \nInput Strength:{Iext_str} Background Input:{Ib_str} \nE-I Balance:{bEI}', xy=(0, 0),
                xycoords='figure fraction', xytext=(1.14, 0.17), textcoords='figure fraction', ha='center', fontsize=10)
    figure_name = f'CouplingOnLongterm_AllLayers_Iextdur{Iext_dur}_Iextstr{Iext_str}_Ibstr{Ib_str}_bEI{bEI}_tauVisual_{thalamus_source}.pdf'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
    plt.show()

def multibEI_couplingOnMinmaxRate(data_df, bEI, Ib_str, population, rate_measure, thalamus_source, figure_dir):
    '''
    3.1) Background input strength 
    - x-axis: input strength
    - y-axis: PSP of E3 population
    - one subplot for each bEI value
    - choose rate measure ('lateLongterm' or 'immediateLongterm', 'duringInput', 'baseline')
    '''
    data_df = data_df[data_df['population']==population]

    # plot results
    fig, axes = plt.subplots(len(bEI), 1, figsize=(8,10)) 
    fig.suptitle(f'Minimum and Maximum rate of {population} Population')
    for axs, b in zip(axes, bEI):
        
        data_bEI_df = data_df[data_df['balanceEI']==b]
        
        # Create a colormap
        cmap = cm.get_cmap('Dark2', len(Ib_str))  # Choose a colormap, e.g., 'viridis'
        cmap_max = cm.get_cmap('Dark2', len(Ib_str))  # Choose a colormap, e.g., 'viridis'

        for i, s in enumerate(Ib_str):
            Istrength_df = data_bEI_df[data_bEI_df['BckgndInputStrength']==s]
            color = cmap(i / (len(Ib_str)))  # Normalize i to [0, 1] for colormap
            color_max = cmap_max(i / (len(Ib_str)))  # Normalize i to [0, 1] for colormap
            axs.plot(Istrength_df['globalCoupling'], Istrength_df[f'minRate_{rate_measure}'], label=s, color=color)
            axs.plot(Istrength_df['globalCoupling'], Istrength_df[f'maxRate_{rate_measure}'], color=color_max)

        #sns.lineplot(data_bEI_df, y='maxPotential_longterm', x='coupling_strength', hue='InputStrength')
        axs.set_xlabel('Coupling Strength')
        axs.set_ylabel('Rate (Hz)')
        axs.set_title(f'E-I Balance:{b}')

    #axs.set_xlim([0,100])
    sns.despine(trim=True)
    plt.legend(title=f'Effect of background input and coupling strength on MinMax Rate {rate_measure}', loc='right')
    plt.tight_layout() 
    figure_name = f'multibEI_couplingOnMinmaxRate_{rate_measure}_pop{population}_measure{rate_measure}_Ibstr{Ib_str}_bEI{bEI}_tauVisual_{thalamus_source}.pdf'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
    plt.show()

def multiLayer_couplingOnMinMaxRate(data_df, Iext_dur, Iext_str, Ib_str, bEI, rate_measure, thalamus_source, figure_dir):
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
    fig, axes = plt.subplots(4, 2, figsize=(10, 10), sharey=False, sharex=True)  # Set figure size
    fig.suptitle(f'Effect of Coupling strength on Min and Max firing rate {rate_measure}')
    for layers, axs in zip(layers_all, axes.T):
        for l, ax in zip(layers, axs):
            layer_df = data_df[data_df['population'].isin(l)]
            sns.lineplot(layer_df, y=f'minRate_{rate_measure}', x='globalCoupling', hue='population', ax=ax)
            sns.lineplot(layer_df, y=f'maxRate_{rate_measure}', x='globalCoupling', hue='population', ax=ax, legend=False)
            ax.set_ylabel('Rate (Hz)')
            ax.set_xlabel('Global Coupling Strength')
            ax.legend(prop={'size':8})

        axs[0].set_title(f'Layer 2/3')
        axs[1].set_title(f'Layer 4')
        axs[2].set_title(f'Layer 5')
        axs[3].set_title(f'Layer 6')

    sns.despine(trim=True)
    # set title
    plt.annotate(f'min and max Rate {rate_measure}, \nInput Duration:{Iext_dur} \nInput Strength:{Iext_str} Background Input:{Ib_str} \nE-I Balance:{bEI}', xy=(0, 0), xycoords='figure fraction', textcoords='figure fraction', xytext=(1, 0.2), ha='center', fontsize=10)
    plt.tight_layout() 
    figure_name = f'multiLayer_couplingOnMinMaxRate_{rate_measure}_Iextstr{Iext_str}_Ibstr{Ib_str}_bEI{bEI}_tauVisual_{thalamus_source}.pdf'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
    plt.show()


def inputStrengthOnminMaxpotential(data_df, Iext_str, Ib_str, bEI, population, potential_measure, thalamus_source, figure_dir):
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
    plt.savefig(os.path.join(figure_dir, f'FixedPoint_potential_{potential_measure}_Iextstr{Iext_str}_Ibstr{Ib_str}_bEI{bEI}_pop{population}_tauVisual_{thalamus_source}.pdf'), bbox_inches='tight')
    plt.show()

