'''
Simulating simple dipole
Based on Vincent Chien's Simulation
'''
# %% 
import matplotlib.pyplot as plt
import numpy as np
import h5py
import os, sys
import seaborn as sns
from ast import literal_eval
import csv
csv.field_size_limit(10**7)
import pandas as pd
# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
param_path = "/data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel/Simulations"
if param_path not in sys.path:
    sys.path.append(param_path)
from parameters import Parameter

# %%

def read_3D_potential(filename):
    """Read the potentials that are stored in a hdf5 file.

    
    Args:
        filename (str): filename of population data.
    
    Returns:
        numpy array: full potentials (3 dimensions).

    """
    dataset_name = 'full_potentials'
    with h5py.File(filename, "r") as f:
        data = f[dataset_name][:]
    return data
    

def simDipoles(dipole_setting, pop_counts, indices_E, psp_dir=None, psp_filename=None, psp=None):
    '''
    Simple dipole simulation: 
    We only consider input to E since they have the dipoles that are in the end measured with EEG
    Simplification: In this simple model the direction and strength of the dipole
    is the same for every targeted E population and from every source population. In reality this differs depending on
    the layer in which the presynaptic population is located.
    
    Read the potentials either directly from a parameter or from a file.

    Keyword arguments:
    dipole_setting -- E-> E, PV-> E, SOM-> E, VIP-> E, Th-> E
    pop_counts : nr of populations
    indices_E : on what place are the E population found in the psp matrix
    '''

    # population counts
    nP, nE, nPV, nSOM, nVIP = pop_counts

    # dipoles
    dipoleEE = np.ones(shape=(nE, nE)) * dipole_setting[0]
    dipolePVE = np.ones(shape=(nE, nPV)) * dipole_setting[1]
    dipoleSOME = np.ones(shape=(nE, nSOM)) * dipole_setting[2]
    dipoleVIPE = np.ones(shape=(nE, nVIP)) * dipole_setting[3]
    dipoleThE = np.ones(shape=(nE, 1)) * dipole_setting[4]

    # shape (nE X nP+th), so (4, 13+1)
    all_dipoles = np.concatenate((dipoleEE, dipolePVE, dipoleSOME, dipoleVIPE, dipoleThE), axis=1)

    # Weighted by E cell count 
    somato_params = Parameter()
    cellcounts = somato_params.get_cellcounts()
    # make cellcounts relative
    cellcounts_relative = cellcounts/np.sum(cellcounts)

    cells_E = [cellcounts_relative[i] for i in indices_E]
    cells_E = np.repeat([cells_E], nP+1, axis=0).T # n+1 because of thalamus

    #dipole weighted by cell count of E pops (nE x nP+1) so in this case (4 x 14)
    dipoles = all_dipoles * cells_E 

    # Read in data (all cells x all cells with IextIb x timepoints)
    if psp_filename:
        psp = read_3D_potential(os.path.join(psp_dir, psp_filename))
        psp = np.squeeze(psp)
   

    # extract populations for S1 
    psp_s1 = psp[indices_E[0]:indices_E[0]+nP, indices_E[0]:indices_E[0]+nP]
    # append thalamus 
    from_thal = psp[31:32, indices_E[0]:indices_E[3]+2]
    psp_s1 = np.concat((psp_s1, from_thal))
    to_thal = psp[indices_E[0]:indices_E[3]+2, 31:32]
    to_thal = np.concat((to_thal, psp[31:32, 31:32]))
    psp_s1 = np.concat((psp_s1, to_thal), axis=1)

    # Simulate LFP, shape (nE x timesteps)
    simDipoles = np.zeros((nE, psp_s1.shape[2]))
    for i in range(nE): 
        # for each E this should give a [1 x timesteps] vector 
        # current flow: absolute psp values that come into each E population
        simDipoles[i, :] = np.dot(dipoles[i, :],abs(psp_s1[i, :, :]))

    return simDipoles


def plot_dipoles(simDipoles, input_onset):
    
    # Average of all dipoles [1 x timepoints], it's the dipole direction * currents
    simMEG = np.sum(simDipoles, axis=0)
    time_step = 0.001
    plot_window = 0.2
    start = int(input_onset/time_step)
    stop = int(start + plot_window*1e3)
    steps = np.arange(0, plot_window, time_step)
    fig, axes = plt.subplots(1,2, sharex=True, sharey=True)

    input_names = ['E1', 'E2', 'E3', 'E4']
    for i, names in enumerate(input_names):
        axes[0].plot(steps, simDipoles[i, start:stop], label=names)
    axes[0].legend()
    axes[1].plot(steps, simMEG[start:stop], label='summed')
    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel('ECD (A*m)')
    axes[1].set_ylabel('ECD (A*m)')
    axes[1].legend()
    sns.despine(trim=True)
    plt.show()
    
# %%
