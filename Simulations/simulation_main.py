"""
File: simulation_main.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2025-08-05
Description: Run this file to run the simulation! 

"""

# %%
import numpy as np
import os
import matplotlib.pyplot as plt
import jax.numpy as jnp
from jr_model_jax import JR_Model
import pandas as pd 
import time
import csv

# %% 
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
    axsA3b.plot(steps[start_plot:], rates[0, :3].T[start_plot:], linewidth=1)
    axsA3b.legend(['E', 'PV', 'SOM'])
    figA3b.suptitle('Area 3b')
    plt.tight_layout() 
    plt.legend(['E', 'PV', 'SOM'])

    # plot results for the S1 column 
    figS1, axs = plt.subplots(2, 2, figsize=(8, 5))  # Set figure size
   
    # Plot settings for all subplots
    for i, ax in enumerate(axs.flatten(), start=1):
        if i == 4:
            ax.plot(steps[start_plot:], rates[0, 12].T[start_plot:], linewidth=1)
        else:
            ax.plot(steps[start_plot:], rates[0, ((i-1)*4)+3:i*4+3].T[start_plot:], linewidth=1)
        ax.grid(True)
        ax.set_ylabel('Hz')
    
    axs[0][1].legend(['L2/3', 'L4', 'L5', 'L6'])

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
            ax.plot(steps[start_plot:], rates[0, (i-1)*4+16:i*4+16].T[start_plot:], linewidth=1)
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


def create_Iext(simulation_time, step_size, input_onset, input_duration, input_strength, input_type):
    ''' Creates external input.'''

    Iext = np.zeros(int(simulation_time/step_size))

    if input_type == "step":
        t  = int(input_duration/step_size)
        t0 = int(input_onset/step_size)
        Iext[t0:t0+t] = input_strength
    elif input_type == "background":
        # provide input for the entire simulation duration 
        Iext[:] = input_strength

    return Iext

def create_Ibackground(simulation_time, step_size, input_strength):
    ''' Create Background Input '''
    Ib = np.zeros(int(simulation_time/step_size))
    Ib[:] = input_strength
    return Ib

def save_results_csv(rates, potentials, filedir, filename, full=False):
    '''
    Safe the simulated data in a csv file
    '''    

    population_names = ['E3b','PV3b','SST3b','VIP3b', 'E1','PV1','SST1','VIP1','E2','PV2','SST2','E3','PV3','SST3','E4','PV4','SST4',
                                          'E1S2','PV1S2','SST1S2','VIP1S2','E2S2','PV2S2','SST2S2','E3S2','PV3S2','SST3S2','E4S2','PV4S2','SST4S2']
    cells = np.concatenate((population_names,['ThalE', 'ThalI']))

    filename = filename + '.h5'
    filename_rates = 'rates' + filename

    # only safe every second datapoint
    resolution_tstep = 0.01
    rates_downsampled = rates[:, ::int(1000*resolution_tstep)] 
    rates_df = pd.DataFrame(rates_downsampled.T)
    rates_df.to_hdf(os.path.join(filedir, filename_rates), index=False, key='data', mode='w')

    # sum the potentials together and save them 
    potential_sum = np.sum(potentials, axis=1)
    potential_sum_downsampled = potential_sum[:, ::int(1000*resolution_tstep)]
    potential_df = pd.DataFrame(potential_sum_downsampled.T, columns=cells)
    filename = 'potentials' + filename
    potential_df.to_hdf(os.path.join(filedir, filename), index=False, key='data', mode='w')

    if full:
        # save all potentials additionally
        psp_filename = 'full_' + filename
        write_3D_csv(os.path.join(filedir, psp_filename), potentials)


def write_3D_csv(filename, data):
    '''
    Write results in form of a 3D hdf5 file.
    '''

    with h5py.File(filename, 'w') as f:
        f.create_dataset(dataset_name, data=data, compression='gzip')
        
# %% 
def main():

    save_params = False
    save_results = True
    save_full_potentials = False # if True the potential matrix is 3D, otherwise 2D
    plot = True

    # set coupling strengths, step size and cortex type (visual or somato)
    # connectivity reverse factor is the absolute cell count divided by  
    connect_reverse_factor =  6448 # TODO: adapt this factor also to S2 cell populations!
    # to simulate:
    # thalamus I to E inhibition 
    
    coupling_strengths = [50] # , 150, 200, 250, 300]
    balance_EI = [0.6]  #[0.9 , 0.7, 0.5, 0.3, 0.1]
    step_size = 0.001
    cortex_type = 'somato'
    filedir = '/data/p_02989/Modelling/output/'

    # define input
    input_type = "step" # other options are "step", "baseline" (equals input strength 0) or "background"
    input_onset = 1.001 # in sec
    input_durations = [0.5]  #, 1, 1.5] # np.arange(0, 1, 1) # in sec 
    input_strengths = [80] #[0, 50, 300, 500] #np.arange(0, 500, 100)
    backgrndI_strengths = [5] #[0, 5, 10, 15, 20]

    # connections within the thalamus
    # in this order: tEE, tEI, tIE, tII 
    thal_connect = np.array([0, 0, 0, 0]) 

    for d in input_durations:
        simulation_time = int(input_onset) + 5
        for sb in backgrndI_strengths:
        
            for s in input_strengths:

                # arrays to store rate (for plotting)
                all_rates = []
                all_potentials = []
                all_durations = []
                all_durations_saving = []

                for g in coupling_strengths:
                    for bEI in balance_EI:

                        print('\nInput duration:', d)
                        print('Input strength:', s)
                        print('Background Input strength:', sb)
                        print('g', g)
                        print('bEI', bEI)
                        print(f'Thalamus EtoE:{thal_connect[0]} ItoE: {thal_connect[1]}')
                        print(f'Thalamus EtoI:{thal_connect[2]} ItoI: {thal_connect[3]}')

                        filename = f'_g{g}_bEI{bEI}_Ib{sb}_Iextd{d}_{input_type}Iexts{s}_Ionset{input_onset}_tauVisual_thalJiang_thalUncon_S1S2Uncon'
                        #filename = f'_gE{int(gE*connect_reverse_factor)}gI{int(gI*connect_reverse_factor)}_{cortex_type}_IbStrength{sb}_Iduration{d}_{input_type}IextStrength{s}_Ionset{input_onset}_tauVisual_thalJiang_thalEI0_S1S2'

                        # create input array 
                        Iext = create_Iext(simulation_time, step_size, input_onset, d, s, input_type)
                        Ib = create_Ibackground(simulation_time, step_size, sb)
                        gE = g * bEI /connect_reverse_factor
                        gI = g * (1 - bEI) /connect_reverse_factor
                        # for now we use the same coupling strength for the thalamus connections as for the cortical connections
                        coupling_thalE = gE
                        coupling_thalI = gI
                        print('gE', gE)
                        print('gI', gI)
                        thal_connect_scaled = thal_connect/connect_reverse_factor
                        
                        model = JR_Model(Iext, Ib, gE, gI, coupling_thalE, coupling_thalI, thal_connect_scaled, filedir, filename, step_size, simulation_time)
                        if save_params:
                            # safe connectivty parameter in yaml file 
                            model.p.save_to_yaml(os.path.join(filedir, 'params'+filename), gE, gI, coupling_thalE, coupling_thalI, thal_connect)

                        # perform simulation with current coupling strength g
                        start = time.time()
                        rate, potential = model.run_simulation()
                        rate.block_until_ready() 
                        stop = time.time()
                        duration = stop-start
                        all_durations.append(duration)
                        print('Simulation duration (in s):', duration)

                        # append results
                        all_rates.append(rate)
                        all_potentials.append(potential)

                        if save_results:
                            start = time.time()
                            save_results_csv(rate, potential, filedir, filename, save_full_potentials)
                            stop = time.time()
                            duration = stop-start
                            all_durations_saving.append(duration)
                            print('Saving duration (in s):', duration)

                if plot:
                    all_rates = np.array(all_rates)
                    all_potentials = np.squeeze(np.array(all_potentials))
                    if len(coupling_strengths) == 1:
                        # when to start plotting (in ms)
                        start_plot = 500
                        # plot only one coupling strength value with time on the x-axis
                        plot_results(all_rates, Iext, Ib, step_size, simulation_time, start_plot)
                        #plot_potentials(all_potentials, Iext, Ib, step_size, simulation_time, start_plot)
                    else: 
                        # used to plot with coupling strength on the x-axis and max/min rate on the y
                        plot_minmax(all_rates, coupling_strengths)
 
    print('Mean Simulation duration: ', np.mean(all_durations))
    print('Mean Saving duration: ', np.mean(all_durations_saving))

if __name__ == '__main__':
   main()

# %%
