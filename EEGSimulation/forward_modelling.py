"""
Perform forward modelling with simulated EEG data. 
Steps:
1. Compute forward operator 
2. Perform co-registration
3. 
"""

# %% 
import numpy as np
import mne
from mne.datasets import sample
from mne.datasets import eegbci
from sim_meg import read_3D_potential, get_population_mapping
import os
import json
import h5py
import sys
import matplotlib.pyplot as plt
import time

location = "laptop"
if location == "laptop":
    WDDIR = r"C:\Users\gross\OneDrive - UvA\Documents\IMPRS_Leipzig\MyProject\Modelling\ChienReplication\SomatosensoryLaminarModel"
    SIMDIR = os.path.join(WDDIR, "output")
    DATADIR = "C:\\Users\\gross\\OneDrive - UvA\\Documents\\IMPRS_Leipzig\\MyProject\\Experiment\\Analysis\\LocalCode\\data"
    RECONDIR = os.path.join(DATADIR, 'freesurfer')

if location == "mpi":
    DATADIR = os.getenv('DATADIR')
    RECONDIR = os.getenv('SUBJECTS_DIR')
    SIMDIR = os.getenv("SIMDIR")
    WDDIR = os.getenv("WDDIR")
    figure_dir = os.path.join(SIMDIR, "Figures")

param_path = sys.path.append(os.path.abspath(os.path.join(WDDIR, 'Simulations')))
if param_path not in sys.path:
    sys.path.append(param_path)
from parameters import Parameter

# %%
# setup sample data for forward modelling
data_path = sample.data_path()
subject = "fsaverage"
trans = "fsaverage"
src = data_path / "subjects" / "fsaverage" / "bem" / "fsaverage-ico-5-src.fif"
bem = data_path / "subjects" / "fsaverage" / "bem" / "fsaverage-5120-5120-5120-bem-sol.fif"
(raw_fname,) = eegbci.load_data(subjects=1, runs=[6])
raw = mne.io.read_raw_edf(raw_fname, preload=True)
# Read and set the EEG electrode locations, which are already in fsaverage's
# space (MNI space) for standard_1020:
eegbci.standardize(raw)

# TODO: this should be the actual montage which we used for the experiment
montage = mne.channels.make_standard_montage("standard_1005")
raw.set_montage(montage)

# %%
# load model potentials
sim_file = "full_potentials_fullpotential_test"
potentials = read_3D_potential(os.path.join(SIMDIR, 'simulation_results', sim_file))
print(f"Potentials size: {potentials.shape[0]} x {potentials.shape[1]}")

#%%
start = time.time()
# defining the dipole characteristics
# Read in preprocessing parameters
with open(os.path.join(WDDIR, 'EEGSimulation', 'dipole_parameters.json'), 'r') as json_file:
    params = json.load(json_file)

dipole_length = params['dipole_lengths']
dipole_orientation = params['dipole_orientation']
resistance_factor = 1 

# load cell count parameter
somato_params = Parameter()
cellcounts = somato_params.get_cellcounts(return_A3b=True)


# %%
# All pyramidal cells contribute to the EEG measurement 
# Get population mapping
pop_mapping = get_population_mapping()
# Extract excitatory populations (these generate the main EEG signal)
exc_pops = []
for area in ['A3b', 'A1', 'S2']:
    if area == 'A3b':
        exc_pops.append(pop_mapping[area]['E'])
    else:
        # For A1 and S2, include all layer E populations
        for layer in ['L1_E', 'L4_E', 'L5_E', 'L6_E']:
            exc_pops.append(pop_mapping[area][layer])

# %%
dipole_lengths_A3b = dipole_length['A3b']
dipole_orientation_A3b = dipole_orientation['A3b']
dipole_lengths_A1 = dipole_length['A1']
dipole_orientation_A1 = dipole_orientation['A1']
dipole_lengths_ES2 = dipole_length['S2']
dipole_orientation_ES2 = dipole_orientation['S2']
resistance_factor = 1
cellcounts_A3b = cellcounts[:4]
cellcounts_A3b_relative = cellcounts_A3b/np.sum(cellcounts_A3b)
cellcounts_A1 = cellcounts[4:17]
cellcounts_A1_relative = cellcounts_A1/np.sum(cellcounts_A1)
cellcounts_S2 = cellcounts[17:]
cellcounts_S2_relative = cellcounts_S2/np.sum(cellcounts_S2)
cellcounts_EA3b_relative = cellcounts_A3b_relative[0]
cellcounts_EA1_relative = cellcounts_A1_relative[np.array(exc_pops[1:5])-4]
cellcounts_ES2_relative = cellcounts_S2_relative[np.array(exc_pops[-4:])-17]

# %%

def prepDipoles(dipole_length, dipole_orientation, resistance_factor, cellcountsE_relative):

    # dipoles
    # each dipole set has a value for each source population
    dipole_matrix = []
    for i, s in enumerate(dipole_length):
        dipole = s * dipole_orientation[i] * resistance_factor 
        dipole_matrix.append(dipole)

    dipole_array = np.array(dipole_matrix)

    # Weighted by E cell count 
    dipoles_weighted = dipole_array * cellcountsE_relative

    return dipoles_weighted

# %%
# compute dipoles for A3b, A1, S2 
dipoles_A3b = prepDipoles(dipole_lengths_A3b, dipole_orientation_A3b, resistance_factor, cellcounts_EA3b_relative)
dipoles_A1 = []
dipoles_ES2 = []

for i, layer in enumerate(['L1_E', 'L4_E', 'L5_E', 'L6_E']):
    # compute dipoles for layers in A1
    dipole_layer_A1 = prepDipoles(dipole_lengths_A1[i], dipole_orientation_A1[i], resistance_factor, cellcounts_EA1_relative[i])
    dipoles_A1.append(dipole_layer_A1)

    # compute dipoles for layers in S2
    dipole_layer_S2 = prepDipoles(dipole_lengths_ES2[i], dipole_orientation_ES2[i], resistance_factor, cellcounts_ES2_relative[i])
    dipoles_ES2.append(dipole_layer_S2)

all_dipoles = np.concatenate((dipoles_A3b, *dipoles_A1, *dipoles_ES2))

# %%
# I have the dipole models computed for each area/layer
# Now it needs to be convolved with the simulated data
potentialsEA3b = potentials[exc_pops[0], :-2]
potentialsEA1 = potentials[exc_pops[1:5], :-2]
potentialsES2 = potentials[exc_pops[-4:], :-2]

# for each time point, compute the simulated dipole
nE = 9
simDipoles = np.zeros((nE, potentialsEA1.shape[2]))
simDipoles[0] = np.dot(dipoles_A3b, abs(potentialsEA3b))

for E in range(4):
    simDipoles[E+1] = np.dot(np.concatenate([dipoles_A1[E]]), abs(potentialsEA1[E]))
    simDipoles[E+5] = np.dot(np.concatenate([dipoles_ES2[E]]), abs(potentialsES2[E]))

end = time.time()
print(f"Total time for dipole computation: {end - start:.2f} seconds.")

# %%
# plot the dipoles 

times = np.arange(simDipoles.shape[1]) / raw.info["sfreq"]
area_groups = {
    "A3b": [0],
    "A1": [1, 2, 3, 4],
    "S2": [5, 6, 7, 8],
}
labels = {
    "A3b": ["E"],
    "A1": ["L1_E", "L4_E", "L5_E", "L6_E"],
    "S2": ["L1_E", "L4_E", "L5_E", "L6_E"],
}

fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
for ax, area in zip(axes, ["A3b", "A1", "S2"]):
    area_indices = area_groups[area]
    for idx, label in zip(area_indices, labels[area]):
        ax.plot(times, simDipoles[idx], label=label, linewidth=1.5)
    sum_trace = simDipoles[area_indices].sum(axis=0)
    ax.plot(times, sum_trace, color="black", linewidth=3.0, label="Sum")
    ax.set_title(f"Computed Dipoles - {area}")
    ax.set_ylabel("Dipole (Am)")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", frameon=False)

axes[-1].set_xlabel("Time (s)")
fig.suptitle("Computed Dipoles for E Populations", fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.98])

figuredir = os.path.join('.', 'Figures')
os.makedirs(figuredir, exist_ok=True)
fig.savefig(os.path.join(figuredir, 'computed_dipoles_by_area.png'), dpi=300, bbox_inches='tight')
#plt.show(fig)


# %%
####################################
# I computed the dipoles for all E populations in the model
# Now I need to compute forward solution for each area.
# 1. pick one area
# 2. compute source time course
# 3. from the source time course, compute the EEG signal
####################################
fwd = mne.make_forward_solution(
    raw.info, trans=trans, src=src, bem=bem, eeg=True, mindist=5.0, n_jobs=None
)
leadfield = fwd["sol"]["data"]
print(f"Leadfield size : {leadfield.shape[0]} sensors x {leadfield.shape[1]} dipoles")

# reduce forward solution to one orientation
fwd_fixed = mne.convert_forward_solution(
    fwd, surf_ori=True, force_fixed=True, use_cps=True
)
src_free = fwd["src"]
src_fixed = fwd_fixed["src"]

# %%
start = time.time()
n_events = 1
events = np.zeros((n_events, 3), int)
events[:, 0] = 200 + 500 * np.arange(n_events)  # Events sample.
events[:, 2] = 1  # All events have the sample id.

# %%
# Area 3b
label_file_soma_rh = os.path.join(data_path, 'subjects', 'fsaverage', 'label','rh.BA3b.label')
selected_label_soma_rh = mne.read_label(label_file_soma_rh, 'fsaverage')
#times = np.arange(0, 10, 0.001)  # Simulate for 10 seconds at 1000 Hz
raw.resample(200)
tstep = 1.0 / raw.info["sfreq"]
dipoles_downsampled = simDipoles[:, ::5]  # Downsample to match the time step
source_simulator = mne.simulation.SourceSimulator(src_free, tstep=tstep)
source_simulator.add_data(selected_label_soma_rh, dipoles_downsampled[0], events)
source_simulator.add_data(selected_label_soma_rh, np.sum(dipoles_downsampled[1:4], axis=0), events)
source_simulator.add_data(selected_label_soma_rh, np.sum(dipoles_downsampled[5:8], axis=0), events)

# %%
# Simulate and plot 
# TODO: this needs to be fixed in mne.make_forward_solution() 
# it should be possible to have a non type fwd, since our raw info also 
# is NonType
raw.info["dev_head_t"] = fwd["info"]["dev_head_t"]
raw_simulated = mne.simulation.simulate_raw(raw.info, source_simulator, forward=fwd)

# %%
epochs = mne.Epochs(raw_simulated, events, 1, baseline=None)#, tmin=-0.05, tmax=0.2)
evoked = epochs.average()
end = time.time()
print(f"Total time for forward simulation: {end - start:.2f} seconds.")

fig = evoked.plot(show=False)
figuredir = os.path.join('.', 'Figures')
os.makedirs(figuredir, exist_ok=True)
fig.savefig(os.path.join(figuredir, 'nullingbf_recon_simulated_evoked.pdf'), format='pdf', bbox_inches='tight')
plt.show(fig)

# %%

freqs = np.arange(8, 41, 2)
n_cycles = np.full_like(freqs, 2.0, dtype=float)
tfr = mne.time_frequency.tfr_morlet(
    epochs,
    freqs=freqs,
    n_cycles=n_cycles,
    return_itc=False,
    average=True,
    picks="eeg",
)
power = tfr.data.mean(axis=0)

tfr_fig, tfr_ax = plt.subplots(figsize=(12, 6))
mesh = tfr_ax.pcolormesh(tfr.times, tfr.freqs, power, shading="auto", cmap="viridis")
tfr_ax.set_title("Simulated Epochs Time-Frequency Power")
tfr_ax.set_xlabel("Time (s)")
tfr_ax.set_ylabel("Frequency (Hz)")
tfr_fig.colorbar(mesh, ax=tfr_ax, label="Power ")
tfr_fig.tight_layout()
tfr_fig.savefig(os.path.join(figuredir, 'nullingbf_recon_simulated_tfr.png'), dpi=300, bbox_inches='tight')
plt.show(tfr_fig)

topo_times = np.linspace(evoked.times[0], evoked.times[-1], 3)
topo_fig = evoked.plot_topomap(times=topo_times, show=False)
topo_fig.savefig(os.path.join(figuredir, 'nullingbf_recon_simulated_topomaps.png'), dpi=300, bbox_inches='tight')
plt.show(topo_fig)


