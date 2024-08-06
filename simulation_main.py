import numpy as np
import os
import matplotlib.pyplot as plt
from jr_model import JR_Model
import pandas as pd 
import csv 

def plot_minmax(rates, coupling_strengths):
    minRate = np.min(rates[:,:,-100:],axis=2)
    maxRate = np.max(rates[:,:,-100:],axis=2)

    # plot results
    fig, axs = plt.subplots(4, 1, figsize=(3, 6))  # Set figure size

    # Plot settings for all subplots
    for i, ax in enumerate(axs, start=1):
        ax.plot(coupling_strengths, minRate[:,(i-1)*4:i*4], linewidth=2)
        ax.plot(coupling_strengths, maxRate[:,(i-1)*4:i*4], linewidth=0.5)
        ax.grid(True)
        ax.set_ylabel('Hz')
        ax.legend(['L2/3', 'L4', 'L5', 'L6'])

    # Set titles for each subplot
    axs[0].set_title('Baseline rate (E)')
    axs[1].set_title('Baseline rate (PV)')
    axs[2].set_title('Baseline rate (SOM)')
    axs[3].set_title('Baseline rate (VIP)')
    #sns.despine(trim=True, bottom=True)
    plt.tight_layout() 
    plt.legend()
    plt.show()

def plot_results(rates, Iext, step_size, simulation_time, start_plot):

    steps = np.arange(step_size, simulation_time+step_size, step_size)*1e3

    # plot results
    fig, axs = plt.subplots(5, 1, figsize=(3, 6))  # Set figure size

    plot_rates = rates[0, 0:4, :]
    
    axs[0].plot(steps[start_plot:], Iext[start_plot:])

    # Plot settings for all subplots
    for i, ax in enumerate(axs[1:], start=1):
        ax.plot(steps[start_plot:], rates[0, (i-1)*4:i*4].T[start_plot:], linewidth=1)
        ax.grid(True)
        ax.set_ylabel('Hz')
    
    axs[1].legend(['L2/3', 'L4', 'L5', 'L6'])

    # Set titles for each subplot
    axs[0].set_title('Input')
    axs[1].set_title('Baseline rate (E)')
    axs[2].set_title('Baseline rate (PV)')
    axs[3].set_title('Baseline rate (SOM)')
    axs[4].set_title('Baseline rate (VIP)')
    axs[4].set_xlabel('time (s)')
    
    #sns.despine(trim=True, bottom=True)
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

    return Iext

def save_results_csv(rates, potentials, cortex_type, filedir, filename, summed=True):
    '''
    Safe the simulated data in a csv file
    '''    

    cells = np.array(['E1', 'E2', 'E3', 'E4', 'P1', 'P2', 'P3', 'P4', 'S1', 'S2', 'S3', 'S4', 'V1', 'V2', 'V3', 'V4'])
    if cortex_type == 'somato':
        cells = cells[:13]

    rates_df = pd.DataFrame(rates.T, columns=cells)
    filename = 'rates' + filename + '.csv'
    rates_df.to_csv(os.path.join(filedir, filename), index=False)
    
    if summed:
        potential_sum = np.zeros((rates.shape[0], rates.shape[-1])) # (16x1000)
        for i in range(rates.shape[0]):
            potential_sum[i] = np.sum(potentials[i], axis=0)
        potential_df = pd.DataFrame(potential_sum.T, columns=cells)
        filename = 'potentials' + filename + '.csv'
        potential_df.to_csv(os.path.join(filedir, filename), index=False)
    else: 
        psp_filename = 'full_potentials_' + filename + '.csv'
        write_3D_csv(os.path.join(filedir, psp_filename), potentials)


def write_3D_csv(filename, data):
    '''
    Write results in form of a 3D numpy array into a csv file. 
    '''
    data = data.tolist()
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerows(data)

def main():

    save_results = True
    save_summed_potentials = True # if True the potential matrix is 2D, otherwise 3D
    plot = False

    # set coupling strengths, step size and cortex type (visual or somato)
    coupling_strengths = np.arange(0, 100, 20)
    step_size = 0.001 
    simulation_time = 1.5
    cortex_type = 'somato'
    filedir = 'output'

    # define input
    input_type = "step" # other options are "baseline"
    input_onset = 0.501 # in sec
    input_durations = np.arange(0.04, 0.2, 0.1) # in sec 
    input_strengths = np.arange(0, 20, 2)

    for d in input_durations:
        print('Input duration:', d)
        for s in input_strengths:

            # arrays to store rate (for plotting)
            all_rates = []
            all_potentials = []

            for g in coupling_strengths:
                filename = f'_G{g}_{cortex_type}_Iduration{d}_Istrength{s}_Ionset{input_onset}_tauVisual'

                # create input array 
                Iext = create_Iext(simulation_time, step_size, input_onset, d, s, input_type)
                model = JR_Model(Iext, cortex_type, g, filedir, filename, step_size, simulation_time)

                # perform simulation with current coupling strength g
                rate, potential = model.run_simulation()
                # append results
                all_rates.append(rate)
                all_potentials.append(potential)

                if save_results:
                    save_results_csv(rate, potential, cortex_type, filedir, filename, save_summed_potentials)

            all_rates = np.array(all_rates)
            all_potentials = np.array(all_potentials)

    if plot:
        if len(coupling_strengths) == 1:
            # when to start plotting (in ms)
            start_plot = 500
            # plot only one coupling strength value with time on the x-axis
            plot_results(all_rates, Iext, step_size, simulation_time, start_plot)
        else: 
            # used to plot with coupling strength on the x-axis and max/min rate on the y
            plot_minmax(all_rates, coupling_strengths)
    

if __name__ == '__main__':
   main()