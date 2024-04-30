import numpy as np
import os
import matplotlib.pyplot as plt
from jr_model import JR_Model
import pandas as pd 

def plot_minmax(rates, coupling_strengths):
    minRate = np.min(rates[:,:,-100:],axis=2)
    maxRate = np.max(rates[:,:,-100:],axis=2)

    # plot results
    fig, axs = plt.subplots(4, 1, figsize=(3, 6))  # Set figure size

    # Plot settings for all subplots
    for i, ax in enumerate(axs, start=1):
        ax.plot(coupling_strengths, minRate[:,(i-1)*4:i*4], linewidth=1)
        ax.plot(coupling_strengths, maxRate[:,(i-1)*4:i*4], linewidth=2)
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

def plot_results(rates, Iext, step_size, simulation_time):

    steps = np.arange(step_size, simulation_time+step_size, step_size)*1e3

    # plot results
    fig, axs = plt.subplots(5, 1, figsize=(3, 6))  # Set figure size

    plot_rates = rates[0, 0:4, :]
    
    axs[0].plot(steps, Iext)

    # Plot settings for all subplots
    for i, ax in enumerate(axs[1:], start=1):
        ax.plot(steps, rates[0, (i-1)*4:i*4].T, linewidth=1)
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

def create_Iext(simulation_time, step_size, input_type):
    ''' Creates external input.'''

    Iext = np.zeros(int(simulation_time/step_size))

    if input_type == "step":
        # after 50 ms for 50ms
        x  = 50
        Iext[501:520] = 5

    return Iext

def main():

    # directory to where to save the results
    output_dir = os.path.join('../data', 'firing_rates')
    safe_results = False
    plot = True

    # set coupling strengths, step size and cortex type (visual or somato)
    coupling_strengths = [1] # np.arange(0, 100, 5)
    step_size = 0.001 
    simulation_time = 1
    cortex_type = 'visual'

    # define input
    input_type = "baseline" # other options are "baseline"
    Iext = create_Iext(simulation_time, step_size, input_type)

    # arrays to store rate
    all_rates = []
    all_potentials = []

    model = JR_Model(Iext, cortex_type, step_size, simulation_time)

    for g in coupling_strengths:

        # perform simulation with current coupling strength g
        rate, potential = model.run_simulation(g)
        
        # append results 
        all_rates.append(rate)
        all_potentials.append(potential)
    
    all_rates = np.array(all_rates)
    all_potentials = np.array(all_potentials)

    if plot:
        if len(coupling_strengths) == 1:
            # plot only one coupling strength value with time on the x-axis
            plot_results(all_rates, Iext, step_size, simulation_time)
        else: 
            # used to plot with coupling strength on the x-axis and max/min rate on the y
            plot_minmax(all_rates, coupling_strengths)

    if safe_results:
        cells = np.array(['E1', 'E2', 'E3', 'E4', 'P1', 'P2', 'P3', 'P4', 'S1', 'S2', 'S3', 'S4', 'V1', 'V2', 'V3', 'V4'])
        for j, g in enumerate(coupling_strengths):
            rates_df = pd.DataFrame(all_rates[j].T, columns=cells)
            rates_df.to_csv(f'output/rates_G{g}.csv', index=False)
            
            potential_sum = np.zeros((all_rates.shape[1], all_rates.shape[-1])) # (16x1000)
            
            for i in range(all_rates.shape[1]):
                potential_sum[i] = np.sum(all_potentials[j][i], axis=0)

            potential_df = pd.DataFrame(potential_sum.T, columns=cells)
            potential_df.to_csv(f'output/potentials_G{g}.csv', index=False)
    
        

if __name__ == '__main__':
   main()