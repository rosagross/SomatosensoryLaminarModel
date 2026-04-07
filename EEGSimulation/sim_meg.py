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
location = "laptop"
if location == "laptop":
    WDDIR = r"C:\Users\gross\OneDrive - UvA\Documents\IMPRS_Leipzig\MyProject\Modelling\ChienReplication\SomatosensoryLaminarModel"
    SIMDIR = os.path.join(WDDIR, "output")
if location == "mpi":
    SIMDIR = os.getenv("SIMDIR")
    WDDIR = os.getenv("WDDIR")
    figure_dir = os.path.join(SIMDIR, "Figures")

param_path = sys.path.append(os.path.abspath(os.path.join(WDDIR, 'Simulations')))
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
        # psp is the potential directed to a single E (nP x times)
        # in dipoles we store the value for a specific E, (1, np)

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

def get_population_mapping():
    """
    Get mapping between model populations and brain regions.
    
    Returns:
        dict: Mapping of population indices to brain regions and layers
    """
    # Population order from parameters.py:
    # A3b: E, PV, SST, VIP (indices 0-3)
    # A1: E1, PV1, SST1, VIP1, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4 (indices 4-16)  
    # S2: E1, PV1, SST1, VIP1, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4 (indices 17-29)
    # Thalamus: ThalE, ThalI (indices 30-31)
    
    mapping = {
        # A3b populations
        'A3b': {
            'E': 0, 'PV': 1, 'SST': 2, 'VIP': 3
        },
        # A1 populations (layers 1-4)
        'A1': {
            'L1_E': 4, 'L1_PV': 5, 'L1_SST': 6, 'L1_VIP': 7,
            'L4_E': 8, 'L4_PV': 9, 'L4_SST': 10,
            'L5_E': 11, 'L5_PV': 12, 'L5_SST': 13,
            'L6_E': 14, 'L6_PV': 15, 'L6_SST': 16
        },
        # S2 populations (layers 1-4)
        'S2': {
            'L1_E': 17, 'L1_PV': 18, 'L1_SST': 19, 'L1_VIP': 20,
            'L4_E': 21, 'L4_PV': 22, 'L4_SST': 23,
            'L5_E': 24, 'L5_PV': 25, 'L5_SST': 26,
            'L6_E': 27, 'L6_PV': 28, 'L6_SST': 29
        },
        # Thalamic populations
        'Thalamus': {
            'E': 30, 'I': 31
        }
    }
    
    return mapping
    
# %%
