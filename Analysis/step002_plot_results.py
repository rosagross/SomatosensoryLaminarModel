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
from plotting_functions import *

colors, _ = figure_style() 

# define directories of stored data
SIMDIR = os.getenv("SIMDIR")
WDDIR = os.getenv("WDDIR")
figure_dir = os.path.join(SIMDIR, "Figures", "global_dynamics")
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

raw_dir = os.path.join(SIMDIR, "simulation_results")
processed_dir = os.path.join(SIMDIR, "derivatives")
# read params
params = load_parameters(WDDIR)


# some global settings
sampling_params = params['sampling']
step_size = sampling_params['step_size']
sample_delay_immediate = sampling_params['sample_delay_immediate']
sample_delay_late = sampling_params['sample_delay_late']  # when to start the long term behaviour "window"
sample_dur = sampling_params['sample_dur']
offset = sampling_params['offset']
input_onset = params['input_onset']
input_type = 'step'
thalamus_source = 'Jiang'
# connectivity params
params = load_parameters(WDDIR)
extI_cellcounts = params['extI_cellcounts']
bI_cellcounts = params['bI_cellcounts']
thal_cellcounts = params['thal_cellcounts']

# %% Make plots that demonstrate the sampling time line 

# load trajectory to plot
g = 20
Iext_dur = 0.006
Iext_str = 0
Ib_str = 4
sI = 0.24
rates_df, potentials_df, filename = load_simulation_data(g, sI, Ib_str, Iext_dur, Iext_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, raw_dir) 
trajectory_df = load_trajectory(rates_df, potentials_df, g, sI, Iext_dur, Iext_str, Ib_str, step_size)

# choose population
population = 'E3b'
line_df = trajectory_df[trajectory_df['population']==population]

# LONG TERM and DURING INPUT
# Add a red vertical line at time of input offset (start of sampling) and stop of sampling
simulation_time = 4*1e-3
start_sample = input_onset + Iext_dur + sample_delay_immediate
stop_sample = start_sample + sample_dur
input_offset = input_onset + Iext_dur
offset = 0.1 # time between baseline sampling and start of input 
baseline_start = input_onset - (sample_dur+offset)
baseline_stop = baseline_start + sample_dur
plotting_cut = 10

# plot the input
steps = np.arange(step_size, simulation_time+step_size, step_size)
input_line = np.zeros(len(line_df))
input_line[int((input_onset)/step_size):int(input_offset/step_size)] = Iext_str
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

#plt.savefig('C:/Users/gross/OneDrive - UvA/Documents/IMPRS_Leipzig/IMPRS SummerSchool/Poster/plotting_windows.pdf')
plt.savefig(os.path.join(figure_dir, 'plotting_windows.pdf'), bbox_inches='tight')
plt.show()

    # %%

'''
1.1) SINGLE PLOT: Plot example trajectory
    - plot style: line plot
    - y axis: rate
    - x axis: time
'''

# choose example settings
g = 10
sI = 0.2
Iext_dur = 0.018
Iext_str = 20
Ib_str = 2

rates_df, potentials_df, filename = load_simulation_data(g, sI, Ib_str, Iext_dur, Iext_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, raw_dir) 
trajectory_df = load_trajectory(rates_df, potentials_df, g, sI, Iext_dur, Iext_str, Ib_str, step_size)

population = 'E3b' # 'E1'
line_df = trajectory_df[trajectory_df['global_coupling']==g]
line_df = line_df[line_df['population']==population]
plotting_window = []

# plot the input
steps = np.arange(step_size, simulation_time+step_size, step_size)
fig = plt.figure(figsize=(10,5))
sns.lineplot(data=line_df, x='time', y='rate', hue='InputStrength', legend='', palette=['grey'])
sns.despine(trim=True)
plt.ylabel('Rate (Hz)')
figure_name = f'trajectory_g{g}sI{sI}_pop{population}_Iduration{Iext_dur}_{Iext_str}.pdf'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()


# %%
"""
Plot a heatmap showing the effect of Input Strength versus Input Duration
"""

# choose a coupling strength, background input strength and a population
g = 20
sI = 0.3
Ib_str = 4
Iext_str = params['input_strengths']
Iext_dur = params['input_durations']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)
population = 'E3b'

data_df = data_df[data_df['globalCoupling']==g]
data_df = data_df[data_df['strength_I']==sI]
data_df = data_df[data_df['population']==population]
data_df = data_df[data_df['BckgndInputStrength']==Ib_str]

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
values = 'longtermVSbaseline_rate'

# choose a coupling strength and a population
g = 20
sI = 0.3
Iext_str = params['input_strengths']
Iext_dur = params['input_durations']
Ib_str = 4
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)

population = 'PV3b'
data_df = data_df[data_df['globalCoupling']==g]
data_df = data_df[data_df['strength_I']==sI]
data_df = data_df[data_df['population']==population]
data_df = data_df[data_df['BckgndInputStrength']==Ib_str]
data_heatmap = data_df.pivot(index='InputStrength',columns='InputDuration', values=values)
sns.heatmap(data_heatmap, cmap='magma')

# %% 

# 1.3) MULTI fingerprint PLOT:

# look at several difference E-I balance values
sIs = params['strength_I']
# choose a global coupling strength and a population
g = params['coupling_strengths']
Ib_str = 4
Iext_str = params['input_strengths']
Iext_dur = params['input_durations']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sIs, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)
population = 'E3b'
sImulti_fingerprint_IextDurVsStr(data_df, g, sIs, Ib_str, population, thalamus_source, figure_dir)

# %%

# 1.3a) Fingerprint: sI vs Input Strength (fixed g)
g = 50
Iext_dur = params['input_durations'][1]
Iext_str = params['input_strengths']
Ib_str = params['backgrndI_strengths'][3]
sIs = params['strength_I']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sIs, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)
population = 'E3b'
fingerprint_sI_vs_IextStr(data_df, g, Iext_dur, Ib_str, population, thalamus_source, figure_dir)


# 1.3b) Fingerprint: g vs Input Strength (fixed sI)
sI = 0.26 #params['strength_I'][0]
Iext_dur = params['input_durations'][1]
Iext_str = params['input_strengths']
Ib_str = params['backgrndI_strengths'][3]
g = params['coupling_strengths']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)
population = 'E3b'
fingerprint_g_vs_IextStr(data_df, sI, Iext_dur, Ib_str, population, thalamus_source, figure_dir)


# %%
'''
1.4) MULTI PLOT: Effect of input intensity and duration on firing rates
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    - subplot columns: populations
    - subplot rows: coupling strengths
'''

g = params['coupling_strengths']
rate_measure = 'longtermVSbaseline_rate'
sI = 0.26
Ib_str = 1
Iext_dur = params['input_durations']
Iext_str = params['input_strengths']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)

populations = np.array(['E1', 'E2', 'E3', 'E4']) 
#populations = np.array(['P1', 'P2', 'P3', 'P4']) 
#populations = np.array(['S1', 'S2', 'S3', 'S4', 'V1']) 

multiPop_heatmap_IextDurVsStr(data_df, g, sI, Ib_str, populations, rate_measure, figure_dir)

# %% 
"""
2.1) Effect of Coupling Strengths on Longterm/steady state
Plot difference between Mininum and Maximum Firing rates 
"""

# choose settings (make sure that there is no input in the samples)
Iext_dur = 0
Iext_str = 0
sI = 0.26
Ib_str = 5
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)
multiLayer_couplingOnLongeterm_diffRate(data_df, Iext_dur, Iext_str, Ib_str, sI, thalamus_source, figure_dir)


# %%
"""
4) Oscillation frequency analysis
"""

# 4.1) Heatmap of dominant frequency (late longterm)
g = 15
sI = 0.26
Iext_str = params['input_strengths']
Iext_dur = params['input_durations']
Ib_str = 5
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)
population = 'E2'
heatmap_frequency_IextDurVsStr(data_df, g, sI, Ib_str, population, "lateLongterm", "Rate", Iext_dur, Iext_str, thalamus_source, figure_dir)

# %%

# 4.2) Coupling strength vs frequency (late longterm)
Iext_dur = 0
Iext_str = 0
sI = 0.28
Ib_str = 5
g = params['coupling_strengths']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)
multiLayer_couplingOnFrequency(data_df, Iext_dur, Iext_str, Ib_str, sI, "lateLongterm", "Potential", thalamus_source, figure_dir)


# 4.3) Heatmap: frequency vs coupling strength and sI (single population)
Iext_dur = 0
Iext_str = 0
Ib_str = 5
g = params['coupling_strengths']
sI = params['strength_I']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)
area = 'A1'
heatmap_frequency_coupling_vs_sI(data_df, area, "lateLongterm", "Potential", Iext_dur, Iext_str, Ib_str, figure_dir)


# 4.4) Frequency vs oscillation amplitude scatter (late longterm)
population = 'E1'
scatter_frequency_vs_diff(data_df, "lateLongterm", "Potential", population, figure_dir)

# %%
'''
3.1) Background input strength 
- x-axis: input strength
- y-axis: PSP of E3 population
- one subplot for each sI value
'''

sI = params['strength_I'][0:4]
Ib_str = params['backgrndI_strengths']
Iext_dur = 0
Iext_str = 0 #, 200, 300, 400]
g = params['coupling_strengths']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)
population = 'E3b'
rate_measure = 'baseline' #('lateLongterm', 'immediateLongterm', 'duringInput', or 'baseline')
multisI_couplingOnMinmaxRate(data_df, sI, Ib_str, population, rate_measure, thalamus_source, figure_dir)

# %%
'''
BACKGROUND INPUT
3.2) all layers
'''
Iext_dur = 0
Iext_str = 0
Ib_str = 5
sI = 0.26
g = params['coupling_strengths']
rate_measure = 'baseline' #('lateLongterm', 'immediateLongterm', 'duringInput', or 'baseline')

data_df = load_all_derivatives(Iext_dur, Iext_str, g, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)
multiLayer_couplingOnMinMaxRate(data_df, Iext_dur, Iext_str, Ib_str, sI, rate_measure, thalamus_source, figure_dir)


# %%  
'''
Fixed point plot
In this plot I use the potential instead of the firing rates.
    - x axis: input strengths 
    - y axis: response (steady-state/longterm) in mV
'''

# load time series data
Iext_dur = 0.5
Iext_str = params['input_strengths']
sI = 0.2
Ib_str = 5
g = params['coupling_strengths']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)

potential_measure = 'lateLongterm'
population = 'E1'
inputStrengthOnminMaxpotential(data_df, Iext_str, Ib_str, sI, population, potential_measure, thalamus_source, figure_dir)


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

Iext_dur = 0.5
Iext_str = [10, 20, 30]
g = [10.0, 20.0, 30.0]
Ib_str = 7
sI = [0.7, 0.8, 0.9]
data_df = load_all_derivatives(Iext_dur, Iext_str, g, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)

rate_measure = 'longtermVSbaseline_rate'

populations = np.array(['E1', 'E2', 'E3', 'E4']) 
#populations = np.array(['P1', 'P2', 'P3', 'P4']) 
#populations = np.array(['S1', 'S2', 'S3', 'S4', 'V1']) 

fig, axes = plt.subplots(len(Iext_str), len(populations), figsize=(20,15) ,sharex=True, sharey=True)

# Create a single colorbar axis
#cbar_ax = fig.add_axes([1.01, 0.3, 0.02, 0.4])
#cbar_ax.set_title(rate_measure)
#cbar_ax.tick_params(labelsize=12) 
fig.suptitle(f'Effect of Input Strength, Global Coupling and E-I balance on {rate_measure} ', fontsize=16)

for i,input_s in enumerate(Iext_str):
    for j,p in enumerate(populations):

        minmax_df = data_df[data_df['InputStrength']==input_s]
        minmax_df = minmax_df[minmax_df['InputDuration']==Iext_dur]
        minmax_df = minmax_df[minmax_df['BckgndInputStrength']==Ib_str]
        minmax_df = minmax_df[minmax_df['population']==p]
        #minmax_df['InputDuration'] = minmax_df['InputDuration'].round(4)

        data_heatmap = minmax_df.pivot(index='globalCoupling', columns='strength_I', values=rate_measure)

        if (minmax_df[rate_measure].isna() | (minmax_df[rate_measure] == 0)).all().all():
            print('black')
            sns.heatmap(data_heatmap, cmap=ListedColormap(['black']), ax=axes[i,j], norm = Normalize(vmin=0, vmax=1))
        else:
            sns.heatmap(data_heatmap, cmap='vlag', ax=axes[i,j], center=0)

        cbar = axes[i, j].collections[0].colorbar
        # here set the labelsize by 20
        cbar.ax.tick_params(labelsize=12)
        axes[i, j].invert_yaxis()
        axes[i, j].set_ylabel('')
        axes[i, j].set_xlabel('')
        axes[i, j].tick_params(axis='both', labelsize=12)
        axes[len(Iext_str)-1, j].set_xlabel('E-I balance')
        axes[0,j].set_title(f'pop: {p}')
        axes[i,0].set_ylabel(f'input strength: {input_s}', rotation=0, labelpad=60)

fig.text(0, 0.2, 'coupling strength', va='center', rotation='vertical')
fig.text(0, 0.5, 'coupling strength', va='center', rotation='vertical')
fig.text(0, 0.83, 'coupling strength', va='center', rotation='vertical')

plt.tight_layout(h_pad=1)
figure_name = f'gvssI_{populations[0][0]}pop_{rate_measure}_tauVisual.png'
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
sIs = [0.2, 0.5, 0.8, 1]
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
fig, axs = plt.subplots(ncols=len(sIs), nrows=len(gs),figsize=(20,10)) 

# Create a colormap
cmap = cm.get_cmap('Dark2', len(input_strengths))  # Choose a colormap, e.g., 'viridis'
cmap_max = cm.get_cmap('Dark2', len(input_strengths))  # Choose a colormap, e.g., 'viridis'
data_df['minPotential_longterm_mV'] = data_df['minPotential_longterm'] *1e3
data_df['maxPotential_longterm_mV'] = data_df['maxPotential_longterm'] *1e3
input_strengths = data_df['InputStrength'].unique()

for i, g in enumerate(gs):
    coupling_df = data_df[data_df['globalCoupling']==g]
    for j, sI in enumerate(sIs):
        balance_df = coupling_df[coupling_df['strength_I']==sI]
        ax = axs[i][j]
        if j ==0:
            ax.set_ylabel(f'g={g}\nRate (Hz)')
        if i == len(gs)-1:
            ax.set_xlabel(f'Background Input Strength\n E-I Balance={sI}')
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

fig, axs = plt.subplots(ncols=len(sIs), nrows=len(gs), figsize=(20,10)) 

data_df = summary_df[summary_df['population']==population]
data_df = data_df[data_df['InputDuration']==Idur]

for i, g in enumerate(gs):
    coupling_df = data_df[data_df['globalCoupling']==g]
    for j, sI in enumerate(sIs):
        balance_df = coupling_df[coupling_df['strength_I']==sI]
        ax = axs[i][j]

        data_heatmap = balance_df.pivot(index='BckgndInputStrength',columns='InputStrength', values='diffRate_duringInput')
        sns.heatmap(data_heatmap, cmap='magma', ax=ax, vmin=0, vmax=40)
        ax.invert_yaxis()
        if j ==0:
            ax.set_ylabel(f'g={g}\n Background Input')
        if i == len(gs)-1:
            ax.set_xlabel(f'Input Strength\n E-I Balance={sI}')
        else:
            ax.set_xlabel('')
            ax.set_ylabel('')

# %%
# PLOT 3: long-term behaviour compared to baseline

Idur = 0.5
population = 'E3'

fig, axs = plt.subplots(ncols=len(sIs), nrows=len(gs), figsize=(20,10)) 

data_df = summary_df[summary_df['population']==population]
data_df = data_df[data_df['InputDuration']==Idur]

for i, g in enumerate(gs):
    coupling_df = data_df[data_df['globalCoupling']==g]
    for j, sI in enumerate(sIs):
        balance_df = coupling_df[coupling_df['strength_I']==sI]
        ax = axs[i][j]

        data_heatmap = balance_df.pivot(index='BckgndInputStrength',columns='InputStrength', values='longtermVSbaseline_rate')
        sns.heatmap(data_heatmap, cmap='magma', ax=ax, vmin=0, vmax=10)
        ax.invert_yaxis()
        print(j)
        if i ==0:
            ax.set_ylabel(f'g={g}\n Background Input')
        if i == len(gs)-1:
            ax.set_xlabel(f'Input Strength\n E-I Balance={sI}')
        else:
            ax.set_xlabel('')
            ax.set_ylabel('')
# %%
