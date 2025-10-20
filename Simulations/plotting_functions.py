import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from plotting_style import figure_style

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

def plot_all_potentials(all_potentials, Iext, Ib, step_size, simulation_time, start_plot, area='all'):
    steps = np.arange(step_size, simulation_time+step_size, step_size)*1e3
    
    # Area 3b
    fig, axes = plt.subplots(4, 1, figsize=(14, 10), sharex=True, sharey=True)
    for i, ax in enumerate(axes.flatten()):
        ax.plot(steps, all_potentials[i,:,:].T)
        ax.set_ylabel('mV')
    plt.show()

    # Area A1
    fig, axes = plt.subplots(4, 4, figsize=(14, 10), sharex=True, sharey=True)
    for i, ax in enumerate(axes.flatten()):
        ax.plot(steps, all_potentials[i+4,:,:].T)
        ax.set_ylabel('mV')
    plt.show()

    # Area S2 & Thal
    fig, axes = plt.subplots(4, 4, figsize=(14, 10), sharex=True, sharey=True)
    for i, ax in enumerate(axes.flatten()):
        if i==15:
            break
        ax.plot(steps, all_potentials[i+17,:,:].T)
        ax.set_ylabel('mV')
    
    plt.show()

    
# TODO: adapt this function to plot all areas
def plot_area_potentials(steps, potential):
    fig, axes = plt.subplots(4, 4, figsize=(14, 10), sharex=True, sharey=True)
    for i, ax in enumerate(axes.flatten()):
        ax.plot(steps, potential)
        ax.set_ylabel('mV')
    plt.show()


def plot_potentials(potentials, Iext, Ib, step_size, simulation_time, start_plot, figdir, bEI, g, d, sb, s, area='all'):
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
    bEI : float
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

        plt.tight_layout()
        figdir = os.path.join(figdir, 'single_simulations')
        if not os.path.exists(figdir):
            os.makedirs(figdir)
        plt.savefig(os.path.join(figdir, f'population_potentials_bEI-{bEI}_g-{g}_area-{area}_Iextdur-{d}_Iextstr-{s}_Ibstr-{sb}.pdf'), dpi=300)
        #plt.show()
            

    elif area=='A1':

        # plot results for the S1 column 
        figS1, axs = plt.subplots(2, 2, figsize=(8, 5))  # Set figure size
    
        # Plot settings for all subplots
        axs_flat = axs.flatten()

        # Layer 2/3
        axs_flat[0].plot(steps[start_plot:], potentials[:4].T[start_plot:])
        axs_flat[0].set_title('L2/3')

        # Layer 4
        axs_flat[1].plot(steps[start_plot:], potentials[4:4+3].T[start_plot:])
        axs_flat[1].set_title('L4')
        # Layer 5
        axs_flat[2].plot(steps[start_plot:], potentials[4+3:4+6].T[start_plot:])
        # Layer 6
        axs_flat[3].plot(steps[start_plot:], potentials[4+6:4+9].T[start_plot:])
        
        axs_flat[0].legend(['E', 'PV', 'SST', 'VIP'])

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


def plot_results(rates, Iext, Ib, step_size, simulation_time, start_plot, bEI, g, area, d, sb, s, figure_dir):
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

    annotate_fig(f'bEI={np.round(bEI, 4)}, g={np.round(g, 4)}, area={area}')
    sns.despine(trim=True)
    plt.tight_layout() 
    plt.legend()
    figdir = os.path.join(figure_dir, 'single_simulations')
    if not os.path.exists(figdir):
        os.makedirs(figdir)

    plt.savefig(os.path.join(figdir, f'population_rates_bEI-{bEI}_g-{g}_area-{area}_Iextdur-{d}_Iextstr-{s}_Ibstr-{sb}.pdf'), dpi=300)
    #plt.show()

def annotate_fig(dataname):
    """ Write on the figure with which data it was generated, the date and the script name."""
    plt.annotate(dataname, xy=(-2, 0.3), xycoords='axes fraction', fontsize=8, ha='center')
    plt.annotate(datetime.datetime.now(), xy=(-2, 0.2), xycoords='axes fraction', fontsize=8, ha='center')
    plt.annotate(f"generated in {os.path.basename(__file__)}", xy=(-2, 0.1), xycoords='axes fraction', fontsize=8, ha='center')

