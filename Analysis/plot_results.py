'''
Plots:
1) Effect of input intensity and duration on firing rates and potentials
    1.1) Line plot, single population - plots trajectory
    1.2) Heatmap, single population - difference between longterm and baseline 
    1.3) Heatmap, single population - plots the different dynamic functions (memory, transfer, ..) 
        - x axis: input duration  
        - y axis: input strength
        - color map: dynamic functions
    1.4) Heatmap, multiple populations - longterm versus baseline by G, input duration and strength
2) Baseline activity and coupling strength 
    2.1) Line plot, multiple populations - axis x:G, y: rate in (Hz)
3) Background activity in steady state
    3.1) Average longterm activity
        - x axis: input strength
        - y axis: Steady state potential of Layer 5 population 

4) Effect of connection strength from and to PV interneurons
    - 

'''

# %% 
import numpy as np
import datetime
import json
import os
from matplotlib.ticker import FormatStrFormatter, FuncFormatter, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm  # Import the colormap module
from matplotlib.colors import ListedColormap, Normalize, BoundaryNorm
from scipy.signal import find_peaks
import seaborn as sns
import pandas as pd
from plotting_style import figure_style
from helper_functions import *

colors, _ = figure_style() 

SIMDIR = os.getenv("SIMDIR")
WDDIR = os.getenv("WDDIR")
figure_dir = os.path.join(SIMDIR, "Figures", "global_dynamics")
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)
    

# define directories of stored data and figures
output_dir = os.path.join(SIMDIR, "simulation_results")

# %%

# Read in preprocessing parameters
with open(os.path.join(WDDIR, 'Simulations', 'simulation_parameter.json'), 'r') as json_file:
    params = json.load(json_file)

# read in data
coupling_strengths = params['coupling_strengths'] # coupling_strengths

# coupling strengths, balance and area selection
balance_EI = params['balance_EI']
g_thal = params['g_thal']
bEI_thal = params['bEI_thal']
step_size = params['step_size']
area = params['area']
filedir = params['filedir']

# inputs
input_type = params['input_type']
input_onset = params['input_onset']
simulation_dur = params['simulation_dur']
input_durations = params['input_durations']
input_strengths = params['input_strengths']
backgroundI_strengths = params['backgrndI_strengths']

# connectivity 
thal_connect = np.array(params['thal_connect'])
extI_cellcount = params['extI_cellcounts']
bI_cellcount = params['bI_cellcounts']
thalamus_cellcount = params['thal_cellcounts']

sample_delay_immediate = 0.3
sample_delay_late = 1 # when to start the long term behaviour "window"
sample_dur = 0.3 # amount of time in sec in which we look at the long term firing rate (min and max)
cortex_type = 'somato'
stimulation_type = 'step'
thalamus_cellcount = 500
extI_cellcount = 1000
bI_cellcount = 100
name_addOn = ''
load_trajectory = True
load_full_potentials = False
load_population_potential = 7 # note: order is E1, P1, S1, V1, E2, ... (idx 7 is E3)
thalamus_source = 'thalJiang'

# %% load rate and potential files 
summary_df, trajectory_df, potentials_df  = read_simulation_data(output_dir, figure_dir, input_durations, input_strengths, coupling_strengths, balance_EI, backgroundI_strengths,
                        step_size, sample_delay_immediate, sample_delay_late, input_onset, sample_dur, cortex_type, stimulation_type, thalamus_cellcount, bI_cellcount, extI_cellcount,
                        load_trajectory, load_full_potentials, load_population_potential)

# %% 
# Save the summary data frame in a separate CSV, so that it does not take that much time to load anymore ...
#summary_df.to_hdf(f'sampledelay{sample_delay_immediate}_late{sample_delay_late}_sampleduration{sample_dur}_{cortex_type}_{thalamus_source}_S1S2_{datetime.date.today()}_{name_addOn}.h5', key='data', index=True) #, index_label='pop')
#trajectory_df.to_hdf(f'trajectories_stepAndBackground_g_bEI_sampledelay{sample_delay_immediate}_late{sample_delay_late}_sampleduration{sample_dur}_{cortex_type}_{thalamus_source}_S1S2_{datetime.date.today()}_thalUncon_S1S2Uncon.h5', index=True, index_label='pop')

# %% Read in summary data frame
#summary_df = pd.read_hdf(f'sampledelay{sample_delay_immediate}_late{sample_delay_late}_sampleduration{sample_dur}_{cortex_type}_{thalamus_source}_S1S2_{datetime.date.today()}_{name_addOn}.h5', index_col=False)

# %% Make plots that demonstrate the sampling time line 

# choose example settings
g = 20
bEI = 0.8
population = 'E3'
input_duration = 1.5
input_strength = 10
backgroundI_strength = 5
line_df = trajectory_df[trajectory_df['global_coupling']==g]
line_df = line_df[line_df['balance_EI']==bEI]
line_df = line_df[line_df['population']==population]
line_df = line_df[line_df['InputDuration']==input_duration]
line_df = line_df[line_df['InputStrength']==input_strength]
line_df = line_df[line_df['BckgndInputStrength']==backgroundI_strength]
print(len(line_df))
print(len(trajectory_df))
# LONG TERM and DURING INPUT

# Add a red vertical line at time of input offset (start of sampling) and stop of sampling
simulation_time = 4*1e-3
start_sample = input_onset + input_duration + sample_delay_immediate
stop_sample = start_sample + sample_dur
input_offset = input_onset + input_duration
offset = 0.1 # time between baseline sampling and start of input 
baseline_start = input_onset - (sample_dur+offset)
baseline_stop = baseline_start + sample_dur
plotting_cut = 10

# plot the input
steps = np.arange(step_size, simulation_time+step_size, step_size)
input_line = np.zeros(len(line_df))
input_line[int((input_onset)/step_size):int(input_offset/step_size)] = input_strength
fig = plt.figure(figsize=(10,5))
#plt.plot(steps[:-plotting_cut], input_line[:-plotting_cut], color='green', linewidth=2)
sns.lineplot(data=line_df[:-plotting_cut], x='time', y='rate', hue='InputStrength', legend='', palette=['grey'])
#sns.lineplot(data=line_df[:-plotting_cut], x='time', y='potential', hue='InputStrength', legend='', palette=['grey'])
plt.axvline(x=baseline_start, color='purple', linestyle='--', linewidth=1, label='Baseline Sample')
plt.axvline(x=baseline_stop, color='purple', linestyle='--', linewidth=1, label='')
plt.axvspan(baseline_start, baseline_stop, alpha=0.2, color='purple')
plt.axvline(x=start_sample, color='red', linestyle='--', linewidth=1, label='Long Term Sample')
plt.axvline(x=stop_sample, color='red', linestyle='--', linewidth=1, label='')
plt.axvspan(start_sample, stop_sample, alpha=0.2, color='red')
plt.axvline(x=input_onset, color='blue', linestyle='--', linewidth=1, label='During Input Sample')
plt.axvline(x=input_offset, color='blue', linestyle='--', linewidth=1, label='')
plt.axvspan(input_onset, input_offset, alpha=0.2, color='blue')
plt.ylabel('Rate (Hz)')
plt.xlabel('Time (sec)')
plt.legend()
sns.despine(trim=True)

plt.savefig('C:/Users/gross/OneDrive - UvA/Documents/IMPRS_Leipzig/IMPRS SummerSchool/Poster/plotting_windows.pdf')
plt.show()

# %%

'''
1.1) SINGLE PLOT: Plot example trajectory
    - plot style: line plot
    - y axis: rate
    - x axis: time
'''

# choose example settings
g = 20
bEI = 0.8
population = 'E3' # 'E1'
input_duration = 1.5
input_strength = 10
line_df = trajectory_df[trajectory_df['global_coupling']==g]
line_df = line_df[line_df['balance_EI']==bEI]
line_df = line_df[line_df['population']==population]
line_df = line_df[line_df['InputDuration']==input_duration]
line_df = line_df[line_df['InputStrength']==input_strength]
line_df = line_df[line_df['BckgndInputStrength']==backgroundI_strength]
plotting_window = []

# plot the input
steps = np.arange(step_size, simulation_time+step_size, step_size)
fig = plt.figure(figsize=(10,5))
sns.lineplot(data=line_df[90:2000], x='time', y='rate', hue='InputStrength', legend='', palette=['grey'])
sns.despine(trim=True)
plt.ylabel('Rate (Hz)')
figure_name = f'trajectory_g{g}bEI{bEI}_pop{population}_Iduration{input_duration}_{input_strength}.pdf'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()


# %%

# choose a coupling strength, background input strength and a population
g = 20
bEI = 0.8
backgroundI_strength = 7
population = 'E3'
data_df = summary_df[summary_df['globalCoupling']==g]
data_df = data_df[data_df['balanceEI']==bEI]
data_df = data_df[data_df['population']==population]
data_df = data_df[data_df['BckgndInputStrength']==backgroundI_strength]

data_heatmap = data_df.pivot(index='InputStrength',columns='InputDuration', values='longtermVSbaseline_rate')
sns.heatmap(data_heatmap, cmap='magma')

# %%
'''
1.2) SINGLE PLOT: Effect of input intensity and duration on potentials/rates in the steady state in comparison to the baseline
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    - measure: longtermVSbaseline potential or rate
'''

# choose a coupling strength and a population
g = 20
bEI = 0.7
population = 'PV2'
data_df = summary_df[summary_df['globalCoupling']==g]
data_df = data_df[data_df['balanceEI']==bEI]
data_df = data_df[data_df['population']==population]
data_df = data_df[data_df['BckgndInputStrength']==backgroundI_strength]
data_heatmap = data_df.pivot(index='InputStrength',columns='InputDuration', values='longtermVSbaseline_rate')
sns.heatmap(data_heatmap, cmap='magma')

# %% 
'''
1.3) MULTI PLOT: Effect of input intensity and duration on dynamic behaviour (non-responsive, transfer and memory)
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    - plotwise: coupling strength
    - measure: dynamic function 
'''

# look at several difference E-I balance values
bEIs = balance_EI
# choose a global coupling strength and a population
g = 20
backgroundI_strength = 7
cbar_ticks = ['non-responsive', 'transfer', 'memory']
population = 'E2'
data_df = summary_df[summary_df['population']==population]
data_df = data_df[data_df['globalCoupling']==g]
data_df = data_df[data_df['BckgndInputStrength']==backgroundI_strength]
data_df['InputDuration'] = data_df['InputDuration'].round(4)
fig, ax = plt.subplots(2, int(len(bEIs)/2), figsize=(9,4), sharex=True, sharey=True)

# Define a fixed discrete colormap with colors for 1, 2, and 3
n = 3 # there are three different functions of the dynamics --> make discrete colormap
colors = sns.color_palette("Pastel2", n)
cmap = ListedColormap(colors)

# Define the boundaries for normalization
bounds = [0.5, 1.5, 2.5, 3.5]
norm = BoundaryNorm(bounds, ncolors=cmap.N)

for i, axis in enumerate(ax.flatten()):
    plot_df = data_df[data_df['balanceEI']==bEIs[i]]
    data_heatmap = plot_df.pivot(index='InputStrength',columns='InputDuration', values='dynamic_function_potential')
    if i == len(bEIs)-1:
        heat_ax = sns.heatmap(data_heatmap, cmap=cmap, norm=norm, ax=axis, cbar=True)
        # Add colorbar, make sure to specify tick locations to match desired ticklabels
        cbar = heat_ax.collections[0].colorbar
        cbar.set_ticks([1, 2, 3])
        cbar.set_ticklabels(cbar_ticks)
    else:
        sns.heatmap(data_heatmap, cmap=cmap, norm=norm, cbar=False, ax=axis, vmin=1, vmax=3, cbar_kws={'ticks': cbar_ticks})

    # Only set x-axis label for the bottom row
    if i >= len(bEIs) // 2:
        axis.set_xlabel('Input Duration')
    else:
        axis.set_xlabel('')
    axis.set_ylabel('')
    axis.invert_yaxis()
    axis.set_title(f'bEI = {bEIs[i]}, g = {g}')

#ax[0][0].set_ylabel('Input Strength')
#ax[1][0].set_ylabel('Input Strength')

plt.tight_layout(h_pad=1)
figure_name = f'_{population}pop_tauVisual_{thalamus_source}.pdf'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()

# %%

'''
1.4) MULTI PLOT: Effect of input intensity and duration on firing rates
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    - subplot columns: populations
    - subplot rows: coupling strengths
'''

rate_measure = 'longtermVSbaseline_rate'
bEI = 0.8
backinput = 7
populations = np.array(['E1', 'E2', 'E3', 'E4']) 
#populations = np.array(['P1', 'P2', 'P3', 'P4']) 
#populations = np.array(['S1', 'S2', 'S3', 'S4', 'V1']) 

fig, axes = plt.subplots(len(coupling_strengths), len(populations), figsize=(20,15) ,sharex=True, sharey=True)

# Create a single colorbar axis
#cbar_ax = fig.add_axes([1.01, 0.3, 0.02, 0.4])
#cbar_ax.set_title(rate_measure)
#cbar_ax.tick_params(labelsize=12) 


for i,g in enumerate(coupling_strengths):
    for j,p in enumerate(populations):

        minmax_df = summary_df[summary_df['globalCoupling']==g]
        minmax_df = minmax_df[minmax_df['balanceEI']==bEI]
        minmax_df = minmax_df[minmax_df['BckgndInputStrength']==backinput]
        minmax_df = minmax_df[minmax_df['population']==p]
        minmax_df['InputDuration'] = minmax_df['InputDuration'].round(4)
        minmax_df[rate_measure] = minmax_df[rate_measure].round(5)

        data_heatmap = minmax_df.pivot(index='InputStrength',columns='InputDuration', values=rate_measure)

        if (minmax_df[rate_measure].isna() | (minmax_df[rate_measure] == 0)).all().all():
            print('black')
            sns.heatmap(data_heatmap, cmap=ListedColormap(['black']), ax=axes[i,j], norm = Normalize(vmin=0, vmax=1))
        else:
            sns.heatmap(data_heatmap, cmap='magma', ax=axes[i,j], vmin=-2, vmax=2)

        cbar = axes[i, j].collections[0].colorbar
        # here set the labelsize by 20
        cbar.ax.tick_params(labelsize=12)
        axes[i, j].invert_yaxis()
        axes[i, j].set_ylabel('')
        axes[i, j].set_xlabel('')
        axes[i, j].tick_params(axis='both', labelsize=12)
        axes[len(coupling_strengths)-1, j].set_xlabel('Input Duration')
        axes[0,j].set_title(f'pop: {p}')
        axes[i,0].set_ylabel(f'g: {g}, g: {bEI}', rotation=0, labelpad=60)

fig.text(0.04, 0.2, 'Input Strength', va='center', rotation='vertical')
fig.text(0.04, 0.5, 'Input Strength', va='center', rotation='vertical')
fig.text(0.04, 0.83, 'Input Strength', va='center', rotation='vertical')

plt.tight_layout(h_pad=1)
figure_name = f'inputDurationVSinputStrength_{populations[0][0]}pop_{rate_measure}_tauVisual.png'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()

# %% 
"""
2.1) BASELINE vs Coupling Strengths
Plot Mininum and Maximum Firing rates 
"""

# choose settings (make sure that there is no input in the samples)
input_duration = 1.5
input_strength = 10
bEI = 0.8
backinput = 5
data_df = summary_df[summary_df['InputDuration']==input_duration]
data_df = data_df[data_df['balanceEI']==bEI]
data_df = data_df[data_df['InputStrength']==input_strength]
data_df = data_df[data_df['BckgndInputStrength']==backinput]

# separate data in layers
layers = [['E1', 'PV1', 'SST1', 'VIP1'], ['E2', 'PV2', 'SST2'], ['E3', 'PV3', 'SST3'], ['E4', 'PV4', 'SST4']]

# plot results
fig, axs = plt.subplots(4, 1, figsize=(3, 6), sharey=False, sharex=True)  # Set figure size

for l, ax in zip(layers, axs):
    layer_df = data_df.loc[l]
    layer_df['population'] = layer_df.index
    #sns.lineplot(layer_df, y='minRate_longterm', x='globalCoupling', hue='population', ax=ax)
    #sns.lineplot(layer_df, y='maxRate_longterm', x='globalCoupling', hue='population', ax=ax, legend=False)
    sns.lineplot(layer_df, y='diffRate_lateLongterm', x='globalCoupling', hue='population', ax=ax, legend=False)
    ax.set_ylabel('Rate (Hz)')
    ax.set_xlabel('Global Coupling Strength (g)')
    ax.legend(prop={'size':8})

axs[0].set_title(f'Layer 2/3')
axs[1].set_title(f'Layer 4')
axs[2].set_title(f'Layer 5')
axs[3].set_title(f'Layer 6')
axs[1].set_ylim([0,100])
axs[2].set_ylim([0,100])
axs[3].set_ylim([0,100])
axs[0].set_ylim([0,100])
#axs[3].set_xlim([30,60])
sns.despine(trim=True)
plt.tight_layout() 
figure_name = f'BaselineAllLayers_tauVisual_{thalamus_source}.pdf'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()


# %%
'''
3.1) Background input strength 
- x-axis: input strength
- y-axis: PSP of E3 population
'''

bEI = 0.8
backinputs = [5, 6, 7] 
input_duration = 1.5
input_strength = 10 #, 200, 300, 400]
population = 'E3'
summary_df['population'] = summary_df.index
data_df = summary_df[summary_df['population']==population]
data_df = data_df[data_df['balanceEI']==bEI]
data_df = data_df[data_df['InputDuration']==input_duration]
data_df = data_df[data_df['BckgndInputStrength'].isin(backinputs)]
data_df = data_df[data_df['InputStrength']==input_strength]

# plot results
fig, axs = plt.subplots(figsize=(6,3)) 

# Create a colormap
cmap = cm.get_cmap('Dark2', len(backinputs))  # Choose a colormap, e.g., 'viridis'
cmap_max = cm.get_cmap('Dark2', len(backinputs))  # Choose a colormap, e.g., 'viridis'
data_df['minPotential_lateLongterm'] = data_df['minPotential_lateLongterm'] *1e3
data_df['maxPotential_lateLongterm'] = data_df['maxPotential_lateLongterm'] *1e3
backinputs = data_df['BckgndInputStrength'].unique()

for i, s in enumerate(backinputs):
    Istrength_df = data_df[data_df['BckgndInputStrength']==s]
    color = cmap(i / (len(backinputs)))  # Normalize i to [0, 1] for colormap
    color_max = cmap_max(i / (len(backinputs)))  # Normalize i to [0, 1] for colormap
    plt.plot(Istrength_df['globalCoupling'], Istrength_df['minRate_lateLongterm'], label=s, color=color)
    plt.plot(Istrength_df['globalCoupling'], Istrength_df['maxRate_lateLongterm'], color=color_max)

#sns.lineplot(data_df, y='maxPotential_longterm', x='coupling_strength', hue='InputStrength')
axs.set_xlabel('Coupling Strength')
axs.set_ylabel('Rate (Hz)')
#axs.set_xlim([0,100])
sns.despine(trim=True)
plt.legend(title='BackgroundInput Strength', loc='right')
plt.tight_layout() 
figure_name = f'BackgroundSteadyState_pop{population}_tauVisual_{thalamus_source}.pdf'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()

# %%
'''
BACKGROUND INPUT
3.2) all layers
'''
input_duration = 1.5
input_strength = 10
data_df = summary_df[summary_df['InputDuration']==input_duration]
data_df = data_df[data_df['InputStrength']==input_strength]
bI = 7
data_df = data_df[data_df['BckgndInputStrength']==bI]
bEI = 0.8
data_df = data_df[data_df['balanceEI']==bEI]

# separate data in layers
layers = [['E1', 'PV1', 'SST1', 'VIP1'], ['E2', 'PV2', 'SST2'], ['E3', 'PV3', 'SST3'], ['E4', 'PV4', 'SST4']]

# plot results
fig, axs = plt.subplots(4, 1, figsize=(3, 6), sharey=False, sharex=True)  # Set figure size

for l, ax in zip(layers, axs):
    layer_df = data_df[data_df['population'].isin(l)]
    sns.lineplot(layer_df, y='minRate_lateLongterm', x='globalCoupling', hue='population', ax=ax)
    sns.lineplot(layer_df, y='maxRate_lateLongterm', x='globalCoupling', hue='population', ax=ax, legend=False)
    ax.set_ylabel('Rate (Hz)')
    ax.set_xlabel('Global Coupling Strength')
    ax.legend(prop={'size':8})

axs[0].set_title(f'Layer 2/3')
axs[1].set_title(f'Layer 4')
axs[2].set_title(f'Layer 5')
axs[3].set_title(f'Layer 6')

sns.despine(trim=True)
# set title
plt.annotate(f'min and max Rate_lateLongterm, \nInput Duration:{input_duration} \nInput Strength:{input_strength} Background Input:{bI} \nE-I Balance:{bEI}', xy=(0, 0), xycoords='figure fraction', xytext=(1.5, 0.2), textcoords='figure fraction', ha='center', fontsize=10)
plt.tight_layout() 
figure_name = f'BackgroundAllLayers_{input_strength}Istrength_tauVisual_{thalamus_source}.pdf'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()

# %%

'''
1.3) Line plot
    - x axis: time
    - y axis: rate
'''

# For this I need a dataframe with all timepoints
bEI = 0.8
coupling_strengths = [10, 20, 30]
population = 'E3'
input_duration = 1.5
input_strength = 10
line_df = trajectory_df.loc[trajectory_df['global_coupling'].isin(coupling_strengths)]
line_df = line_df[line_df['balance_EI']==bEI]
line_df = line_df[line_df['population']==population]
line_df = line_df[line_df['InputDuration']==input_duration]
line_df = line_df[line_df['InputStrength']==input_strength]
plt.title(f'Population:{population} Input Duration:{input_duration} Input Strength{input_strength}')
sns.lineplot(data=line_df, x='time', y='rate', hue='global_coupling', 
             palette=sns.dark_palette("#69d", reverse=False))

figure_name = f'population{population}_inputDur{input_duration}_inputStrength{input_strength}.png'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')


# %%  
'''
Fixed point plot
In this plot I use the potential instead of the firing rates.
    - x axis: input strengths 
    - y axis: response (steady-state/longterm) in mV
'''

mean_lt = potentials_df.groupby(['InputDuration', 'globalCoupling', 'InputStrength'], as_index=False).mean()

# For this I need a dataframe with all timepoints
coupling_strength = 20
population = 1
input_duration = 0.5
input_strength = np.arange(10, 50, 5)
line_df = mean_lt[mean_lt['globalCoupling']==coupling_strength]
line_df = line_df[line_df['population']==population]
line_df = line_df[line_df['InputDuration']==input_duration]
plt.title(f'Population:{population} Input Duration:{input_duration} Input Strength{input_strength}')
sns.lineplot(data=line_df, x='InputStrength', y='lt_potential') # , hue='coupling_strength')

figure_name = f'fixedPointCurve_population{population}_inputDur{input_duration}_G{coupling_strength}.png'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')

# %%
populations = np.array(['E1', 'E2', 'E3', 'E4', 'PV1', 'PV2', 'PV3', 'PV4', 'SST1', 'SST2', 'SST3', 'SST4', 'VIP1'])#, 'V2', 'V3', 'V4'])

fig, axes = plt.subplots(4, 4, sharex=True, sharey=False)
g = 100
bEI = 0.5
bI = 5
for ax, pop in zip(axes.flatten(), populations):
    line_df = trajectory_df[trajectory_df['global_coupling']==g]
    line_df = line_df[line_df['balance_EI']==bEI]
    line_df = line_df[line_df['BckgndInputStrength']==bI]
    line_df = line_df[line_df['population']==pop]
    line_df = line_df[line_df['InputDuration']==input_duration]
    sns.lineplot(data=line_df, x='time', y='rate', hue='InputStrength', ax=ax)
    ax.legend('')

handles, labels = plt.gca().get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center')
plt.tight_layout()
plt.legend()
plt.show()


# %% Line Plot 
'''
Look at the change of rate difference comparing populations
- x axis: input duration 
- y axis: firing rate diff
- hue: populations 
'''



# %% Heatmap gE versus gI

'''
Look at the interaction between global coupling and E-I balance in different populations.
'''

input_duration = 1.5
input_strength = [10]
coupling_strengths = [10, 20, 30]
backinput = 7
rate_measure = 'longtermVSbaseline_rate'

populations = np.array(['E1', 'E2', 'E3', 'E4']) 
#populations = np.array(['P1', 'P2', 'P3', 'P4']) 
#populations = np.array(['S1', 'S2', 'S3', 'S4', 'V1']) 

fig, axes = plt.subplots(len(input_strength), len(populations), figsize=(20,15) ,sharex=True, sharey=True)

# Create a single colorbar axis
#cbar_ax = fig.add_axes([1.01, 0.3, 0.02, 0.4])
#cbar_ax.set_title(rate_measure)
#cbar_ax.tick_params(labelsize=12) 


for i,input_s in enumerate(input_strength):
    for j,p in enumerate(populations):

        minmax_df = summary_df[summary_df['InputStrength']==input_s]
        minmax_df = minmax_df[minmax_df['InputDuration']==input_duration]
        minmax_df = minmax_df[minmax_df['BckgndInputStrength']==backinput]
        minmax_df = minmax_df[minmax_df['population']==p]
        #minmax_df['InputDuration'] = minmax_df['InputDuration'].round(4)

        data_heatmap = minmax_df.pivot(index='globalCoupling',columns='balanceEI', values=rate_measure)

        if (minmax_df[rate_measure].isna() | (minmax_df[rate_measure] == 0)).all().all():
            print('black')
            sns.heatmap(data_heatmap, cmap=ListedColormap(['black']), ax=axes[i,j], norm = Normalize(vmin=0, vmax=1))
        else:
            sns.heatmap(data_heatmap, cmap='magma', ax=axes[i,j], vmin=-1, vmax=2)

        cbar = axes[i, j].collections[0].colorbar
        # here set the labelsize by 20
        cbar.ax.tick_params(labelsize=12)
        axes[i, j].invert_yaxis()
        axes[i, j].set_ylabel('')
        axes[i, j].set_xlabel('')
        axes[i, j].tick_params(axis='both', labelsize=12)
        axes[len(input_strength)-1, j].set_xlabel('E-I balance')
        axes[0,j].set_title(f'pop: {p}')
        axes[i,0].set_ylabel(f'input strength: {input_s}', rotation=0, labelpad=60)

fig.text(0.04, 0.2, 'global g', va='center', rotation='vertical')
fig.text(0.04, 0.5, 'global g', va='center', rotation='vertical')
fig.text(0.04, 0.83, 'global g', va='center', rotation='vertical')

plt.tight_layout(h_pad=1)
figure_name = f'gvsbEI_{populations[0][0]}pop_{rate_measure}_tauVisual.png'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()

# %%
'''
Choose the 4 most important variables:
- global g
- E-I balance
- input strength
- input duration (or background strength)
Output value:
- difference of max. vs min. long term versus baseline rate --> oscillatory behaviour

For the output population of S1 (E3), for every global g and E-I balance parameter
plot the input strength and duration. 

The goal is to create an overview of the system behaviour: 
global x an y are global coupling strength (g) and E-I ratio
- plot 1 (lineplot): fixed points and limit cycles
    x-axis: background input 
    y-axis: long term firing rate (min and max)
- plot 2 (heatmap): long-term behaviour after input
    x-axis: stimulus duration 
    y-axis: stimulus strength
    value: difference between min and max value
- plot 3 (heatmap): systems responsiveness
    x-axis: stimulus duration 
    y-axis: stimulus strength
    value: 
        memory (long-term behaviour the same as during input)
        transfer (long-term behaviour goes back to base behaviour)
        non-responsive (no change of activity during input)

'''


# create figure with multiple subplots
fig, axes = plt.subplots(ncols=2, nrows=2)

# iterate through the global coupling strength values 
#for g in coupling_strengths:



# %% 
# PLOT 1: fixed points and limit cycles
bEIs = [0.2, 0.5, 0.8, 1]
gs = coupling_strengths

input_duration = 0.5
backgrndI_strengths = [0, 5, 10, 15, 20]
input_strengths = [100] #, 120, 140] #, 120, 140, 160, 180]
population = 'E3'
summary_df['population'] = summary_df.index
data_df = summary_df[summary_df['population']==population]
data_df = data_df[data_df['InputDuration']==input_duration]
data_df = data_df[data_df['InputStrength'].isin(input_strengths)]

# plot results (balance E-I and global coupling on the global axis)
fig, axs = plt.subplots(ncols=len(bEIs), nrows=len(gs),figsize=(20,10)) 

# Create a colormap
cmap = cm.get_cmap('Dark2', len(input_strengths))  # Choose a colormap, e.g., 'viridis'
cmap_max = cm.get_cmap('Dark2', len(input_strengths))  # Choose a colormap, e.g., 'viridis'
data_df['minPotential_longterm_mV'] = data_df['minPotential_longterm'] *1e3
data_df['maxPotential_longterm_mV'] = data_df['maxPotential_longterm'] *1e3
input_strengths = data_df['InputStrength'].unique()

for i, g in enumerate(gs):
    coupling_df = data_df[data_df['globalCoupling']==g]
    for j, bEI in enumerate(bEIs):
        balance_df = coupling_df[coupling_df['balanceEI']==bEI]
        ax = axs[i][j]
        if j ==0:
            ax.set_ylabel(f'g={g}\nRate (Hz)')
        if i == len(gs)-1:
            ax.set_xlabel(f'Background Input Strength\n E-I Balance={bEI}')
        for k, s in enumerate(input_strengths):
            Istrength_df = balance_df[balance_df['InputStrength']==s]
            color = cmap(k / (len(input_strengths)))  # Normalize i to [0, 1] for colormap
            color_max = cmap_max(k / (len(input_strengths)))  # Normalize i to [0, 1] for colormap
            ax.plot(Istrength_df['BckgndInputStrength'], Istrength_df['minRate_longterm'], label=s, color=color)
            ax.plot(Istrength_df['BckgndInputStrength'], Istrength_df['maxRate_longterm'], color=color_max)

#sns.lineplot(data_df, y='maxPotential_longterm', x='coupling_strength', hue='InputStrength')
#axs.set_xlim([0,100])
sns.despine(trim=True)
axs[0][-1].legend(title='Input Strength', loc='right')
plt.tight_layout() 
figure_name = f'BackgroundSteadyState_pop{population}_tauVisual_{thalamus_source}.pdf'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()

# %%
# PLOT 2: behaviour during input

Idur = 0.5
population = 'E3'

fig, axs = plt.subplots(ncols=len(bEIs), nrows=len(gs), figsize=(20,10)) 

data_df = summary_df[summary_df['population']==population]
data_df = data_df[data_df['InputDuration']==Idur]

for i, g in enumerate(gs):
    coupling_df = data_df[data_df['globalCoupling']==g]
    for j, bEI in enumerate(bEIs):
        balance_df = coupling_df[coupling_df['balanceEI']==bEI]
        ax = axs[i][j]

        data_heatmap = balance_df.pivot(index='BckgndInputStrength',columns='InputStrength', values='diffRate_duringInput')
        sns.heatmap(data_heatmap, cmap='magma', ax=ax, vmin=0, vmax=40)
        ax.invert_yaxis()
        if j ==0:
            ax.set_ylabel(f'g={g}\n Background Input')
        if i == len(gs)-1:
            ax.set_xlabel(f'Input Strength\n E-I Balance={bEI}')
        else:
            ax.set_xlabel('')
            ax.set_ylabel('')

# %%
# PLOT 3: long-term behaviour compared to baseline

Idur = 0.5
population = 'E3'

fig, axs = plt.subplots(ncols=len(bEIs), nrows=len(gs), figsize=(20,10)) 

data_df = summary_df[summary_df['population']==population]
data_df = data_df[data_df['InputDuration']==Idur]

for i, g in enumerate(gs):
    coupling_df = data_df[data_df['globalCoupling']==g]
    for j, bEI in enumerate(bEIs):
        balance_df = coupling_df[coupling_df['balanceEI']==bEI]
        ax = axs[i][j]

        data_heatmap = balance_df.pivot(index='BckgndInputStrength',columns='InputStrength', values='longtermVSbaseline_rate')
        sns.heatmap(data_heatmap, cmap='magma', ax=ax, vmin=0, vmax=10)
        ax.invert_yaxis()
        print(j)
        if i ==0:
            ax.set_ylabel(f'g={g}\n Background Input')
        if i == len(gs)-1:
            ax.set_xlabel(f'Input Strength\n E-I Balance={bEI}')
        else:
            ax.set_xlabel('')
            ax.set_ylabel('')
# %%
