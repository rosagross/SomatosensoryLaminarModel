#%% 
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.api.queries.cell_types_api import CellTypesApi
from allensdk.core.cell_types_cache import ReporterStatus as RS
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
output_dir = '.'
import os

# %%
# Instantiate the CellTypesCache instance.  The manifest_file argument
# tells it where to store the manifest, which is a JSON file that tracks
# file paths.  If you supply a relative path it will go into your
# current working directory
try:
    ctc = CellTypesCache(manifest_file=Path(output_dir) / 'manifest.json')
except:
    os.remove(Path(output_dir) / 'manifest.json')
    ctc = CellTypesCache(manifest_file=Path(output_dir) / 'manifest.json')

# %%        

# Get datasets from all relevant specimen ids (layer 2-6 of SSp)
# this saves the NWB file to 'cell_types/specimen_464212183/ephys.nwb'

# download all electrophysiology features for all cells
ephys_features = ctc.get_ephys_features()
ef_df = pd.DataFrame(ephys_features)


L23 = 577813827
L4 = [585058163, 577821067]
L5 = [585018101, 584635140]
L6b = 504179737

current_specifmen_id = L4[0]
data_set = ctc.get_ephys_data(current_specifmen_id)


# %% get electrophysiology features


print("Ephys. features available for %d cells" % len(ef_df))

# filter down to a specific cell
cell_ephys_features = ef_df[ef_df['specimen_id']== current_specifmen_id]

# %% 

sweep_number = 38
sweep_data = data_set.get_sweep(sweep_number)

index_range = sweep_data["index_range"]
i = sweep_data["stimulus"][0:index_range[1]+1] # in A
v = sweep_data["response"][0:index_range[1]+1] # in V
i *= 1e12 # to pA
v *= 1e3 # to mV

sampling_rate = sweep_data["sampling_rate"] # in Hz
t = np.arange(0, len(v)) * (1.0 / sampling_rate)
t_vec = np.arange(0, len(v)) / t

extractor = EphysSweepFeatureExtractor(t=t, v=v, i=i)
extractor.process_spikes()
spike_times = extractor.spike_feature("peak_t")
firing_rate = len(spike_times) / 1
stimulus_intensity = np.max(i)
print(stimulus_intensity)

print('spike times', spike_times)
print('f rate', firing_rate)


plt.style.use('ggplot')
fig, axes = plt.subplots(2, 1, sharex=True)
axes[0].plot(t, v, color='black')
axes[1].plot(t, i, color='gray')
axes[0].set_ylabel("mV")
axes[1].set_ylabel("pA")
axes[1].set_xlabel("seconds")
plt.show()

# %% Plot the f-i curve 
# Example usage

# Fetch the sweep table to find sweeps of interest
sweep_info = ctc.get_ephys_sweeps(current_specifmen_id)
print(sweep_info)

sweep_numbers = [sweep['sweep_number'] for sweep in sweep_info if sweep['stimulus_name'] == 'Long Square']

fi_curve = []
vi_curve = []

for sweep_number in sweep_numbers:
    sweep_data = data_set.get_sweep(sweep_number)

    index_range = sweep_data["index_range"]
    i = sweep_data["stimulus"][0:index_range[1]+1] # in A
    v = sweep_data["response"][0:index_range[1]+1] # in V
    i *= 1e12 # to pA
    v *= 1e3 # to mV
 
    #print(np.unique(i)) # every sweep there are different intensities applied.
    # I need the mV response for each intensity
    #vi_curve.append((stimulus_intensity, ))

    sampling_rate = sweep_data["sampling_rate"] # in Hz
    t = np.arange(0, len(v)) * (1.0 / sampling_rate)
    t_vec = np.arange(0, len(v)) / t

    extractor = EphysSweepFeatureExtractor(t=t, v=v, i=i)
    extractor.process_spikes()
    spike_times = extractor.spike_feature("peak_t")
    if len(spike_times) > 0:
        firing_rate = len(spike_times) / 1
    else:
        firing_rate = 0
    stimulus_intensity = np.max(i)
    fi_curve.append((i, v))

# %% 

for curve in fi_curve:
    plt.scatter(curve[0], curve[1])

plt.show()

#fi_curve = calculate_fi_curve(data_set, sweep_numbers)

# %%
currents, firing_rates = zip(*fi_curve)
plt.figure(figsize=(8, 5))
plt.plot(currents, firing_rates, marker='o',linestyle='None')
plt.xlabel('Injected Current (pA)')
plt.ylabel('Firing Rate (Hz)')
plt.title('F-I Curve')
plt.grid(True)
plt.show()


# %% Get average firing rate for all cell types

print(type(list(cell_ephys_features.keys())))
hi = any('rate' in s for s in list(cell_ephys_features.keys()))
print(hi)

# %% 
# Get the electrophysiology sweeps data
sweeps = ctc.get_ephys_sweeps(current_specifmen_id)
# Print keys from the first sweep to see what data is available
if sweeps:
    print(sweeps[0].keys())

