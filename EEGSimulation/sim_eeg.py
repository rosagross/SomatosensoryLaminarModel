'''
EEG Simulation using the full somatosensory model
Based on Vincent Chien's Simulation with extensions for the complete model
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
import mne
from mne.datasets import sample
import json

# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
param_path = "/data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel/Simulations"
if param_path not in sys.path:
    sys.path.append(param_path)
from parameters import Parameter
from somato_model import SomatoModel

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
    

def setup_eeg_forward_model(subject='fsaverage', subjects_dir=None):
    """
    Set up EEG forward model using fsaverage template.
    
    Args:
        subject (str): Subject name (default: 'fsaverage')
        subjects_dir (str): Path to FreeSurfer subjects directory
        
    Returns:
        dict: Forward model components
    """
    if subjects_dir is None:
        # Try to get from environment or use common paths
        subjects_dir = os.getenv('SUBJECTS_DIR')
        if subjects_dir is None:
            subjects_dir = '/usr/local/freesurfer/subjects'
    
    print(f"Setting up forward model for {subject}")
    print(f"Subjects directory: {subjects_dir}")
    
    # Create a template montage for EEG
    montage = mne.channels.make_standard_montage('standard_1020')
    
    # Create info structure for EEG
    ch_names = montage.ch_names
    ch_types = ['eeg'] * len(ch_names)
    info = mne.create_info(ch_names=ch_names, ch_types=ch_types, sfreq=1000.0)
    info.set_montage(montage)
    
    # Set up source space (surface-based)
    try:
        src = mne.setup_source_space(
            subject, spacing='oct6', add_dist='patch', subjects_dir=subjects_dir
        )
        print(f"Source space created with {len(src)} sources")
    except Exception as e:
        print(f"Error creating source space: {e}")
        print("Creating simplified source space...")
        # Create a simple sphere model if fsaverage is not available
        sphere = (0.0, 0.0, 0.0, 0.09)  # (x, y, z, radius) in meters
        src = mne.setup_volume_source_space(
            subject=None, pos=10.0, sphere=sphere, mindist=5.0
        )
    
    # Create BEM model (if available)
    try:
        conductivity = (0.3,)  # Single layer for simplicity
        model = mne.make_bem_model(
            subject=subject, ico=4, conductivity=conductivity, subjects_dir=subjects_dir
        )
        bem = mne.make_bem_solution(model)
        print("BEM model created successfully")
    except Exception as e:
        print(f"Error creating BEM model: {e}")
        print("Using sphere model for forward calculation")
        bem = None
    
    # Create forward solution
    try:
        if bem is not None:
            fwd = mne.make_forward_solution(
                info, trans=None, src=src, bem=bem, meg=False, eeg=True,
                mindist=5.0, n_jobs=1, verbose=True
            )
        else:
            # Use sphere model
            sphere = mne.make_sphere_model(r0=(0, 0, 0), head_radius=0.09)
            fwd = mne.make_forward_solution(
                info, trans=None, src=src, bem=sphere, meg=False, eeg=True,
                mindist=5.0, n_jobs=1, verbose=True
            )
        
        print(f"Forward solution created: {fwd['sol']['data'].shape}")
        return {
            'fwd': fwd,
            'info': info,
            'src': src,
            'bem': bem
        }
    except Exception as e:
        print(f"Error creating forward solution: {e}")
        return None


def get_population_mapping():
    """
    Get mapping between model populations and brain regions.
    
    Returns:
        dict: Mapping of population indices to brain regions and layers
    """
    # Population order from parameters.py:
    # A3b: E, PV, SST, VIP (indices 0-3)
    # S1: E1, PV1, SST1, VIP1, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4 (indices 4-16)  
    # S2: E1, PV1, SST1, VIP1, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4 (indices 17-29)
    # Thalamus: ThalE, ThalI (indices 30-31)
    
    mapping = {
        # A3b populations
        'A3b': {
            'E': 0, 'PV': 1, 'SST': 2, 'VIP': 3
        },
        # S1 populations (layers 1-4)
        'S1': {
            'L1_E': 4, 'L1_PV': 5, 'L1_SST': 6, 'L1_VIP': 7,
            'L2_E': 8, 'L2_PV': 9, 'L2_SST': 10,
            'L3_E': 11, 'L3_PV': 12, 'L3_SST': 13,
            'L4_E': 14, 'L4_PV': 15, 'L4_SST': 16
        },
        # S2 populations (layers 1-4)
        'S2': {
            'L1_E': 17, 'L1_PV': 18, 'L1_SST': 19, 'L1_VIP': 20,
            'L2_E': 21, 'L2_PV': 22, 'L2_SST': 23,
            'L3_E': 24, 'L3_PV': 25, 'L3_SST': 26,
            'L4_E': 27, 'L4_PV': 28, 'L4_SST': 29
        },
        # Thalamic populations
        'Thalamus': {
            'E': 30, 'I': 31
        }
    }
    
    return mapping


def simEEG(sim_model, forward_model, dipole_params=None):
    '''
    Simulate EEG data from the somatosensory model.
    
    Args:
        sim_model: SomatoModel instance with simulation results
        forward_model: Dictionary containing forward model components
        dipole_params: Dictionary with dipole parameters (optional)
        
    Returns:
        numpy array: Simulated EEG data (n_channels x n_timepoints)
    '''
    
    if dipole_params is None:
        dipole_params = load_dipole_parameters()
    
    # Get population mapping
    pop_mapping = get_population_mapping()
    
    # Extract excitatory populations (these generate the main EEG signal)
    exc_pops = []
    for area in ['A3b', 'S1', 'S2']:
        if area == 'A3b':
            exc_pops.append(pop_mapping[area]['E'])
        else:
            # For S1 and S2, include all layer E populations
            for layer in ['L1_E', 'L2_E', 'L3_E', 'L4_E']:
                exc_pops.append(pop_mapping[area][layer])
    
    n_exc = len(exc_pops)
    n_timepoints = sim_model.potential.shape[2]
    
    print(f"Simulating EEG from {n_exc} excitatory populations")
    print(f"Time points: {n_timepoints}")
    
    # Initialize EEG data array
    n_channels = len(forward_model['info']['ch_names'])
    eeg_data = np.zeros((n_channels, n_timepoints))
    
    # Get forward solution
    fwd = forward_model['fwd']
    leadfield = fwd['sol']['data']  # Shape: (n_channels, n_dipoles)
    
    # For each time point, compute EEG signal
    for t in range(n_timepoints):
        # Get potentials for excitatory populations at this time point
        # We sum the absolute values of PSPs as in the original simDipoles function
        pop_potentials = np.zeros(n_exc)
        for i, pop_idx in enumerate(exc_pops):
            # Sum absolute potentials from all presynaptic populations
            pop_potentials[i] = np.sum(np.abs(sim_model.potential[pop_idx, :, t]))
        
        # Apply dipole strengths (if provided)
        if 'dipole_strengths' in dipole_params:
            strengths = np.array(dipole_params['dipole_strengths'])
            pop_potentials = pop_potentials * strengths[:n_exc]
        
        # Project to EEG channels using forward model
        # For simplicity, we'll use a subset of dipoles that correspond to our populations
        n_dipoles_per_pop = leadfield.shape[1] // n_exc
        for i in range(n_exc):
            dipole_start = i * n_dipoles_per_pop
            dipole_end = (i + 1) * n_dipoles_per_pop
            
            # Average the leadfield for this population
            pop_leadfield = np.mean(leadfield[:, dipole_start:dipole_end], axis=1)
            
            # Add contribution to EEG
            eeg_data[:, t] += pop_potentials[i] * pop_leadfield
    
    return eeg_data, forward_model['info']


def load_dipole_parameters(filename='dipole_parameters.json'):
    """
    Load dipole parameters from JSON file.
    
    Args:
        filename (str): Path to JSON file
        
    Returns:
        dict: Dipole parameters
    """
    try:
        with open(filename, 'r') as f:
            params = json.load(f)
        # Ensure we have enough dipole strengths for our populations
        if 'dipole_strengths' not in params or len(params['dipole_strengths']) == 0:
            params['dipole_strengths'] = [1.0] * 10  # Default strengths
        return params
    except FileNotFoundError:
        print(f"Dipole parameters file {filename} not found, using defaults")
        return {
            'dipole_strengths': [1.0] * 10,  # Default strengths
            'dipole_orientation': [1.0] * 10,  # Default orientations
            'resistance_factor': 1.0  # Default resistance
        }


def plot_eeg_results(eeg_data, info, sim_model, input_onset=1.001):
    """
    Plot simulated EEG results.
    
    Args:
        eeg_data: Simulated EEG data (n_channels x n_timepoints)
        info: MNE info structure
        sim_model: SomatoModel instance
        input_onset: Time of input onset in seconds
    """
    # Create MNE Raw object for plotting
    raw = mne.io.RawArray(eeg_data, info)
    
    # Set time info
    raw.set_meas_date(None)
    raw.info['sfreq'] = 1000.0  # Assuming 1kHz sampling
    
    # Plot EEG data
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    # 1. Topographic plot at peak response
    time_step = 0.001
    onset_sample = int(input_onset / time_step)
    peak_time = input_onset + 0.05  # 50ms after onset
    peak_sample = int(peak_time / time_step)
    
    if peak_sample < eeg_data.shape[1]:
        im, _ = raw.plot_topomap(times=peak_time, ch_type='eeg', 
                               axes=axes[0,0], show=False)
        axes[0,0].set_title(f'EEG Topography at {peak_time:.3f}s')
    
    # 2. Time course at a representative electrode
    axes[0,1].plot(raw.times, eeg_data[0, :], 'b-', linewidth=1)
    axes[0,1].axvline(x=input_onset, color='r', linestyle='--', label='Stimulus onset')
    axes[0,1].set_xlabel('Time (s)')
    axes[0,1].set_ylabel('EEG Amplitude')
    axes[0,1].set_title('EEG Time Course')
    axes[0,1].legend()
    axes[0,1].grid(True)
    
    # 3. Global field power
    gfp = np.std(eeg_data, axis=0)
    axes[1,0].plot(raw.times, gfp, 'g-', linewidth=2)
    axes[1,0].axvline(x=input_onset, color='r', linestyle='--', label='Stimulus onset')
    axes[1,0].set_xlabel('Time (s)')
    axes[1,0].set_ylabel('Global Field Power')
    axes[1,0].set_title('Global Field Power')
    axes[1,0].legend()
    axes[1,0].grid(True)
    
    # 4. Source space visualization (if available)
    if hasattr(sim_model, 'rate'):
        time_window = [input_onset, input_onset + 0.1]  # 100ms window
        window_samples = [int(t / time_step) for t in time_window]
        
        # Plot average firing rates
        pop_mapping = get_population_mapping()
        exc_pops = []
        pop_names = []
        
        for area in ['A3b', 'S1', 'S2']:
            if area == 'A3b':
                exc_pops.append(pop_mapping[area]['E'])
                pop_names.append('A3b E')
            else:
                for layer in ['L1_E', 'L2_E', 'L3_E', 'L4_E']:
                    exc_pops.append(pop_mapping[area][layer])
                    pop_names.append(f'{area} {layer}')
        
        avg_rates = np.mean(sim_model.rate[exc_pops, window_samples[0]:window_samples[1]], axis=1)
        
        axes[1,1].bar(range(len(avg_rates)), avg_rates)
        axes[1,1].set_xlabel('Population')
        axes[1,1].set_ylabel('Average Firing Rate (Hz)')
        axes[1,1].set_title('Average Firing Rates')
        axes[1,1].set_xticks(range(len(pop_names)))
        axes[1,1].set_xticklabels(pop_names, rotation=45, ha='right')
        axes[1,1].grid(True)
    
    plt.tight_layout()
    plt.show()
    
    return raw


def run_eeg_simulation(simulation_params=None, forward_model_params=None):
    """
    Main function to run EEG simulation.
    
    Args:
        simulation_params: Parameters for the somatosensory model
        forward_model_params: Parameters for the forward model
        
    Returns:
        dict: Results including EEG data, forward model, and simulation model
    """
    print("Starting EEG simulation...")
    
    # 1. Set up the somatosensory model
    if simulation_params is None:
        simulation_params = {
            'simulation_dur': 2.0,
            'step_size': 0.001,
            'input_onset': 1.001,
            'thal_connect': [0, 0, 0, 0],
            'extI_cellcounts': 1000,
            'balance_EI': 0.5,
            'bI_cellcounts': 100,
            'thal_cellcounts': 500,
            'bEI_thal': 0.5,
            'g_thal': 2.0,
            'input_type': 'step',
            'area': 'all',
            'coupling_strength': 10.0,
            'Ib_strength': 7.0,
            'Iext_strength': 10.0,
            'Iext_duration': 0.5
        }
    
    print("Setting up somatosensory model...")
    model = SomatoModel(simulation_params)
    
    # 2. Run simulation
    print("Running simulation...")
    model.simulate()
    
    # 3. Set up forward model
    if forward_model_params is None:
        forward_model_params = {
            'subject': 'fsaverage',
            'subjects_dir': None
        }
    
    print("Setting up forward model...")
    forward_model = setup_eeg_forward_model(**forward_model_params)
    
    if forward_model is None:
        print("Error: Could not create forward model")
        return None
    
    # 4. Simulate EEG
    print("Simulating EEG...")
    eeg_data, info = simEEG(model, forward_model)
    
    # 5. Plot results
    print("Plotting results...")
    raw = plot_eeg_results(eeg_data, info, model, simulation_params['input_onset'])
    
    # 6. Save results
    results = {
        'eeg_data': eeg_data,
        'info': info,
        'raw': raw,
        'forward_model': forward_model,
        'simulation_model': model,
        'simulation_params': simulation_params
    }
    
    print("EEG simulation completed successfully!")
    return results


if __name__ == "__main__":
    # Example usage
    results = run_eeg_simulation()