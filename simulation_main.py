# %%
import numpy as np
import os
import matplotlib.pyplot as plt
from jr_model import JR_Model
import pandas as pd 
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

    axs[0].plot(steps[start_plot:], Iext[start_plot:], label='Iext rate')
    axs[0].plot(steps[start_plot:], Ib[start_plot:], label='Ib rate')
    axs[0].legend(title='')
    axs[0].set_ylabel('Hz')
    axs[1].plot(steps[start_plot:], rates[0, -2:-1].T[start_plot:], color='purple')
    axs[1].plot(steps[start_plot:], rates[0, -1:].T[start_plot:], color='grey')
    axs[1].legend(['Thalamus E', 'Thalamus I'])
    axs[1].set_ylabel('Hz')

    plt.tight_layout() 

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

    cells = np.array(['E1', 'E2', 'E3', 'E4', 'P1', 'P2', 'P3', 'P4', 'S1', 'S2', 'S3', 'S4', 'V1', 'ThalE', 'ThalI']) 
    
    rates_df = pd.DataFrame(rates.T, columns=cells)
    filename = filename + '.csv'
    filename_rates = 'rates' + filename
    rates_df.to_csv(os.path.join(filedir, filename_rates), index=False)

    # sum the potentials together and save them 
    potential_sum = np.zeros((rates.shape[0], rates.shape[-1])) # (16x1000)
    for i in range(rates.shape[0]):
        potential_sum[i] = np.sum(potentials[i], axis=0)
    potential_df = pd.DataFrame(potential_sum.T, columns=cells)
    filename = 'potentials' + filename
    potential_df.to_csv(os.path.join(filedir, filename), index=False)

    if full:
        # save all potentials additionally
        psp_filename = 'full_' + filename
        write_3D_csv(os.path.join(filedir, psp_filename), potentials)


def write_3D_csv(filename, data):
    '''
    Write results in form of a 3D numpy array into a csv file. 
    '''
    data = data.tolist()
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerows(data)
# %% 
def main():

    save_results = False
    save_full_potentials = False # if True the potential matrix is 3D, otherwise 2D
    plot = True

    # set coupling strengths, step size and cortex type (visual or somato)
    coupling_strengths_E = [60] #np.arange(0, 100, 10)
    coupling_strengths_I = [40] #np.arange(0, 100, 10)
    step_size = 0.001
    cortex_type = 'somato'
    filedir = '' #'/data/p_02989/Modelling/output/'

    # define input
    input_type = "step" # other options are "step", "baseline" (equals input strength 0) or "background"
    input_onset = 1.001 # in sec
    input_durations = [1] #np.arange(0.5, 2, 0.5) # in sec 
    input_strengths = [80] #[0, 20, 40, 60, 80, 100] # np.arange(0, 80, 10)
    backgrndI_strengths = [2] #[1,2,3,4,5,6,7,8,9,10]

    for d in input_durations:
        simulation_time = int(input_onset) + d + 1
        for sb in backgrndI_strengths:
        
            for s in input_strengths:

                # arrays to store rate (for plotting)
                all_rates = []
                all_potentials = []

                for gE in coupling_strengths_E:
                    for gI in coupling_strengths_I:

                        print('\nInput duration:', d)
                        print('Input strength:', s)
                        print('Background Input strength:', sb)
                        print('gE', gE)
                        print('gI', gI)

                        filename = f'_gE{gE}gI{gI}_{cortex_type}_IbStrength{sb}_Iduration{d}_{input_type}IextStrength{s}_Ionset{input_onset}_tauVisual_thalJiang_thalEI0_S1S2'

                        # create input array 
                        Iext = create_Iext(simulation_time, step_size, input_onset, d, s, input_type)
                        Ib = create_Ibackground(simulation_time, step_size, sb)
                        model = JR_Model(Iext, Ib, gE, gI, filedir, filename, step_size, simulation_time)

                        # perform simulation with current coupling strength g
                        rate, potential = model.run_simulation()
                        # append results
                        all_rates.append(rate)
                        all_potentials.append(potential)

                        if save_results:
                            save_results_csv(rate, potential, filedir, filename, save_full_potentials)


                if plot:
                    all_rates = np.array(all_rates)
                    all_potentials = np.squeeze(np.array(all_potentials))
                    if len(coupling_strengths_E) == 1:
                        # when to start plotting (in ms)
                        start_plot = 500
                        # plot only one coupling strength value with time on the x-axis
                        plot_results(all_rates, Iext, Ib, step_size, simulation_time, start_plot)
                        #plot_potentials(all_potentials, Iext, Ib, step_size, simulation_time, start_plot)
                    else: 
                        # used to plot with coupling strength on the x-axis and max/min rate on the y
                        plot_minmax(all_rates, coupling_strengths_E)

# %%       
if __name__ == '__main__':
   main()

# %%
save_results = False
save_full_potentials = False # if True the potential matrix is 3D, otherwise 2D
plot = True

# set coupling strengths, step size
coupling_strengths_E = [60] #np.arange(0, 100, 10)
coupling_strengths_I = [40] #np.arange(0, 100, 10)
step_size = 0.001
cortex_type = 'somato'
filedir = '/data/p_02989/Modelling/output/'

# define input
input_type = "step" # other options are "step", "baseline" (equals input strength 0) or "background"
input_onset = 1.001 # in sec
input_durations = [1] #np.arange(0.5, 2, 0.5) # in sec 
input_strengths = [80] #[0, 20, 40, 60, 80, 100] # np.arange(0, 80, 10)
backgrndI_strengths = [2] #[1,2,3,4,5,6,7,8,9,10]

for d in input_durations:
    simulation_time = int(input_onset) + d + 1
    for sb in backgrndI_strengths:
    
        for s in input_strengths:

            # arrays to store rate (for plotting)
            all_rates = []
            all_potentials = []

            for gE in coupling_strengths_E:
                for gI in coupling_strengths_I:

                    print('\nInput duration:', d)
                    print('Input strength:', s)
                    print('Background Input strength:', sb)
                    print('gE', gE)
                    print('gI', gI)

                    filename = f'_gE{gE}gI{gI}_{cortex_type}_IbStrength{sb}_Iduration{d}_{input_type}IextStrength{s}_Ionset{input_onset}_tauVisual_thalJiang_thalEI0'

                    # create input array 
                    Iext = create_Iext(simulation_time, step_size, input_onset, d, s, input_type)
                    Ib = create_Ibackground(simulation_time, step_size, sb)
                    model = JR_Model(Iext, Ib, gE, gI, filedir, filename, step_size, simulation_time)

                    # perform simulation with current coupling strength g
                    rate, potential = model.run_simulation()
                    # append results
                    all_rates.append(rate)
                    all_potentials.append(potential)

                    if save_results:
                        save_results_csv(rate, potential, filedir, filename, save_full_potentials)


            if plot:
                all_rates = np.array(all_rates)
                all_potentials = np.squeeze(np.array(all_potentials))
                if len(coupling_strengths_E) == 1:
                    # when to start plotting (in ms)
                    start_plot = 500
                    # plot only one coupling strength value with time on the x-axis
                    plot_results(all_rates, Iext, Ib, step_size, simulation_time, start_plot)
                    #plot_potentials(all_potentials, Iext, Ib, step_size, simulation_time, start_plot)
                else: 
                    # used to plot with coupling strength on the x-axis and max/min rate on the y
                    plot_minmax(all_rates, coupling_strengths_E)
# %%
