import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='whitegrid')
from plotting_style import figure_style
figure_style()

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


def plot_potentials(potentials, Iext, Ib, step_size, simulation_time, start_plot, area='all'):

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
        plt.savefig('all_areas_potentials.png', dpi=300)
        plt.show()
            

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
    axs.plot(steps[start_plot:], rates[0, idx_rates].T[start_plot:], linewidth=1)

def plot_results(rates, Iext, Ib, step_size, simulation_time, start_plot):
    steps = np.arange(step_size, simulation_time+step_size, step_size)*1e3
    fig, axs = plt.subplots(1, 2, figsize=(8, 2))  # Set figure size

    # thalamus
    axs[0].plot(steps[start_plot:], Iext[start_plot:], label='Iext rate')
    axs[0].plot(steps[start_plot:], Ib[start_plot:], label='Ib rate')
    axs[0].legend(title='')
    axs[0].set_ylabel('Hz')
    axs[1].plot(steps[start_plot:], rates[0, -2:-1].T[start_plot:], color='purple')
    axs[1].plot(steps[start_plot:], rates[0, -1:].T[start_plot:], color='grey')
    axs[1].legend(['Thalamus E', 'Thalamus I'])
    axs[1].set_ylabel('Hz')

    plt.tight_layout() 

    # area 3b
    figA3b, axsA3b = plt.subplots(1, 1, figsize=(5, 5))
    axsA3b.plot(steps[start_plot:], rates[0, :4].T[start_plot:], linewidth=1)
    figA3b.suptitle('Area 3b')
    plt.tight_layout() 
    plt.legend([f'E {np.round(rates[0, 0].T[-1], 6)}', f'PV {np.round(rates[0, 1].T[-1], 6)}', f'SOM {np.round(rates[0, 2].T[-1], 6)}', f'VIP {np.round(rates[0, 3].T[-1], 6)}'])

    # plot results for the S1 column 
    figS1, axs = plt.subplots(2, 2, figsize=(8, 5))  # Set figure size
   
    # Plot settings for all subplots
    # E
    idxs_E = [0+4, 4+4, 7+4, 10+4]
    for idxE in idxs_E:
        plot_axis(axs[0][0], steps, start_plot, rates, idxE)

    axs[0][0].legend([f'E1 {np.round(rates[0, 0+4].T[-1], 6)}', f'E2 {np.round(rates[0, 4+4].T[-1], 6)}', f'E3 {np.round(rates[0, 7+4].T[-1], 6)}', f'E4 {np.round(rates[0, 10+4].T[-1], 6)}'], loc='upper right')

    # PV
    axs[0][1].plot(steps[start_plot:], rates[0, 1+4].T[start_plot:], linewidth=1)
    axs[0][1].plot(steps[start_plot:], rates[0, 5+4].T[start_plot:], linewidth=1)
    axs[0][1].plot(steps[start_plot:], rates[0, 8+4].T[start_plot:], linewidth=1)
    axs[0][1].plot(steps[start_plot:], rates[0, 11+4].T[start_plot:], linewidth=1)
    axs[0][1].legend([f'P1 {np.round(rates[0, 1+4].T[-1], 6)}', f'P2 {np.round(rates[0, 5+4].T[-1], 6)}', f'P3 {np.round(rates[0, 8+4].T[-1], 6)}', f'P4 {np.round(rates[0, 11+4].T[-1], 6)}'], loc='upper right')

    # SST
    axs[1][0].plot(steps[start_plot:], rates[0, 2+4].T[start_plot:], linewidth=1)
    axs[1][0].plot(steps[start_plot:], rates[0, 6+4].T[start_plot:], linewidth=1)
    axs[1][0].plot(steps[start_plot:], rates[0, 9+4].T[start_plot:], linewidth=1)
    axs[1][0].plot(steps[start_plot:], rates[0, 12+4].T[start_plot:], linewidth=1)
    axs[1][0].legend([f'SOM1 {np.round(rates[0, 2+4].T[-1], 6)}', f'SOM2 {np.round(rates[0, 6+4].T[-1], 6)}', f'SOM3 {np.round(rates[0, 9+4].T[-1], 6)}', f'SOM4 {np.round(rates[0, 12+4].T[-1], 6)}'], loc='upper right')

    # VIP
    axs[1][1].plot(steps[start_plot:], rates[0, 3+4].T[start_plot:], linewidth=1)
    axs[1][1].legend([f'VIP1 {np.round(rates[0, 3+4].T[-1], 6)}'])

    # Set titles for each subplot
    axs[0][0].set_title('E')
    axs[0][1].set_title('PV')
    axs[1][0].set_title('SOM')
    axs[1][1].set_title('VIP')
    axs[1][1].set_xlabel('time (s)')
    figS1.suptitle('S1')

    #sns.despine(trim=True, bottom=True)
    plt.tight_layout() 
    plt.legend()
    
    # plot results S2
    figS2, axsS2 = plt.subplots(2, 2, figsize=(8, 5))  # Set figure size
   
    # Plot settings for all subplots
    for i, ax in enumerate(axsS2.flatten(), start=1):
        if i == 4:
            ax.plot(steps[start_plot:], rates[0, 25].T[start_plot:], linewidth=1)
        else:
            ax.plot(steps[start_plot:], rates[0, (i-1)*4+17:i*4+17].T[start_plot:], linewidth=1)
        ax.grid(True)
        ax.set_ylabel('Hz')
    
    axsS2[1][1].legend(['L2/3', 'L4', 'L5', 'L6'])

    # Set titles for each subplot
    axsS2[0][0].set_title('E')
    axsS2[0][1].set_title('PV')
    axsS2[1][0].set_title('SOM')
    axsS2[1][1].set_title('VIP')
    axsS2[1][1].set_xlabel('time (s)')
    figS2.suptitle('S2')
    plt.tight_layout() 
    plt.legend()
    plt.show()

