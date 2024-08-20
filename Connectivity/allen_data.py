#%% 
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.api.queries.cell_types_api import CellTypesApi
from allensdk.core.cell_types_cache import ReporterStatus as RS
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.interpolate import BPoly
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

current_specifmen_id = L23
data_set = ctc.get_ephys_data(current_specifmen_id)


# %% get electrophysiology features


print("Ephys. features available for %d cells" % len(ef_df))

# filter down to a specific cell
cell_ephys_features = ef_df[ef_df['specimen_id']== current_specifmen_id]

# %% 
# Extract sweep numbers
# Make sure that data_set contains the sweep information, check documentation for exact methods
sweep_numbers = data_set.get_sweep_numbers()

print("Available Sweep Numbers:", sweep_numbers)

# %% EXAMPLE of one sweep
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

# %% ALL SWEEPS

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
    fi_curve.append((stimulus_intensity, firing_rate))
    
    vi_curve.append((i, v))


# %% Plot VI-curve
total_v = np.array([])
total_i = np.array([])
for curve in vi_curve:
    #plt.scatter(curve[0], curve[1])
    array_curve = np.array(curve)
    sorted_inputs = np.argsort(array_curve[0])
    sorted_curve = array_curve[:, sorted_inputs]
    total_i = np.concatenate([total_i, sorted_curve[0]])
    total_v = np.concatenate([total_v, sorted_curve[1]])

vi_array = pd.DataFrame(np.array([total_i,total_v]).T, columns=['I', 'V'])
vi_array = vi_array.groupby(['I'], as_index=False).mean()
vi_array = vi_array[vi_array['I']<150]
sns.lineplot(vi_array, x='I', y='V')
sns.scatterplot(vi_array, x='I', y='V')

#fi_curve = calculate_fi_curve(data_set, sweep_numbers)

# %% Extrapolate curve 
# the goal is to get the relationship between potential and firing rate
# to extract the steepness of the sigmoidal curve

x = vi_array['I']
y = vi_array['V']

# Define the degree of the polynomial
degree = 3  # Change this based on how complex you believe the relationship is
coefficients = np.polyfit(x,y,degree)
polynomial = np.poly1d(coefficients)

# Extrapolate up to I=300
extrapolated_I = np.linspace(min(x), 300, num=100)  # More points for smoother plot
extrapolated_V = polynomial(extrapolated_I)

plt.figure(figsize=(10, 5))
plt.plot(x, y, 'o', label='Original data')
plt.plot(extrapolated_I, extrapolated_V, '-', label='Extrapolated data')
plt.xlabel('I')
plt.ylabel('V')
plt.title('Extrapolation of Electrical Data')
plt.legend()
plt.grid(True)
plt.show()

# Sample values from the extrapolated graph
extrapolated_VIcurve = 

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

