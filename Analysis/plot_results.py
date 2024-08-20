'''
Plots:
1) Effect of input intensity and duration on firing rates and potentials
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    1.1) SINGLE PLOT - RATES
    1.2) SINGLE PLOT - POTENTIALS
    1.3) Line Plot: 
        - x axis: time  
        - y axis: rate (Hz)
        - line color: input strengths
2) Interaction of stimulus intensity and coupling strength 
3) Effect of connection strength from and to PV interneurons
    - 

'''
# %% 
import numpy as np
import os
from matplotlib.ticker import FormatStrFormatter, FuncFormatter, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, Normalize
from scipy.signal import find_peaks
import seaborn as sns
import pandas as pd
from plotting_style import figure_style
from helper_functions import *
import ast

colors, _ = figure_style() 

# define directories of stored data and figures
output_dir = "../output"
figure_dir = "../Figures"

# %%

# read in data
input_durations = [0.25] # np.arange(0, 1.5, 0.25) # in sec 
input_strengths = [4, 16] #np.arange(0, 20, 2)
coupling_strengths = [20, 40] # np.arange(0, 50, 10)
step_size = 0.001
sample_delay = 0.5
input_onset = 1.001
sample_dur = 0.1 # amount of time in sec in which we look at the long term firing rate (min and max)
cortex_type = 'somato'
thalamus_source = 'thalJi'
load_trajectory = True
load_full_potentials = False
load_population_potential = 3 # the 3rd population is E3, the output layer

summary_df, trajectory_df, potentials_df  = read_simulation_data(output_dir, figure_dir, input_durations, input_strengths, coupling_strengths, 
                        step_size, sample_delay, input_onset, sample_dur, cortex_type, thalamus_source, load_trajectory, load_full_potentials, load_population_potential)
    
# %% Make plots that demonstrate the sampling time line 

# choose example settings
coupling_strength = 20
population = 'E3'
input_duration = 0.25
input_strength = 16
line_df = trajectory_df[trajectory_df['coupling_strength']==coupling_strength]
line_df = line_df[line_df['population']==population]
line_df = line_df[line_df['InputDuration']==input_duration]
line_df = line_df[line_df['InputStrength']==input_strength]
# LONG TERM and DURING INPUT

# Add a red vertical line at time of input offset (start of sampling) and stop of sampling
simulation_time = len(line_df)*1e-3
start_sample = input_onset + input_duration + sample_delay
stop_sample = start_sample + sample_dur
input_offset = input_onset + input_duration
offset = 0.1 # time between baseline sampling and start of input 
baseline_start = input_onset - (sample_dur+offset)
baseline_stop = baseline_start + sample_dur
plotting_cut = 400

# plot the input
steps = np.arange(step_size, simulation_time+step_size, step_size)
input_line = np.zeros(len(line_df))
input_line[int((input_onset)/step_size):int(input_offset/step_size)] = input_strength
fig = plt.figure(figsize=(10,5))
plt.plot(steps[:-plotting_cut], input_line[:-plotting_cut], color='green', linewidth=2)
sns.lineplot(data=line_df[:-plotting_cut], x='time', y='rate', hue='InputStrength', legend='', palette=['grey'])
#sns.lineplot(data=line_df[:-plotting_cut], x='time', y='potential', hue='InputStrength', legend='', palette=['grey'])
plt.axvline(x=baseline_start, color='purple', linestyle='--', linewidth=1, label='Baseline Sample')
plt.axvline(x=baseline_stop, color='purple', linestyle='--', linewidth=1, label='')
plt.axvline(x=start_sample, color='red', linestyle='--', linewidth=1, label='Long Term Sample')
plt.axvline(x=stop_sample, color='red', linestyle='--', linewidth=1, label='')
plt.axvline(x=input_onset, color='blue', linestyle='--', linewidth=1, label='During Input Sample')
plt.axvline(x=input_offset, color='blue', linestyle='--', linewidth=1, label='')

plt.ylabel('Rate (Hz)')
plt.xlabel('Time (sec)')
plt.legend()
plt.savefig('../Figures/plotting_windows.pdf')
plt.show()

# %%

'''
1.1) SINGLE PLOT: Effect of input intensity and duration on firing rates in the steady state in comparison to the baseline
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    - measure: longtermVSbaseline rate
'''

# choose a coupling strength and a population
coupling_strength = 40
population = 'E1'
data_df = summary_df[summary_df['coupling_strength']==coupling_strength]
data_df = data_df[data_df['population']==population]
data_heatmap = data_df.pivot(index='InputStrength',columns='InputDuration', values='longtermVSbaseline_rate')
sns.heatmap(data_heatmap, cmap='magma')

# %%
'''
1.2) SINGLE PLOT: Effect of input intensity and duration on potentials in the steady state in comparison to the baseline
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    - measure: longtermVSbaseline potential
'''

# choose a coupling strength and a population
coupling_strength = 40
population = 'E1'
data_df = summary_df[summary_df['coupling_strength']==coupling_strength]
data_df = data_df[data_df['population']==population]
data_heatmap = data_df.pivot(index='InputStrength',columns='InputDuration', values='longtermVSbaseline_potential')
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

# choose a coupling strength and a population
coupling_strengths = [20,40]
population = 'E3'
data_df = summary_df[summary_df['population']==population]

fig, ax = plt.subplots(2, 1)

n = 3 # there are three different functions of the dynamics --> make discrete colormap
cmap = sns.color_palette("Pastel2", n)
for axis, G in zip(ax, coupling_strengths):
    plot_df = data_df[data_df['coupling_strength']==G]
    data_heatmap = plot_df.pivot(index='InputStrength',columns='InputDuration', values='dynamic_function_potential')
    sns.heatmap(data_heatmap, cmap=cmap, cbar=False, ax=axis)

# reconstruct color map
# add legend
box = ax[-1].get_position()
ax[-1].set_position([box.x0, box.y0, box.width * 0.7, box.height])

# add color map to legend
legend_ax = fig.add_axes([.7, .5, 1, .1])
legend_ax.axis('off')
patches = [mpatches.Patch(facecolor=c, edgecolor=c) for c in cmap]
legend = legend_ax.legend(patches,
    ['non-responsive', 'transfer', 'memory'],
    handlelength=0.8, loc='lower left')
for t in legend.get_texts():
    t.set_ha("left")

#plt.show()

# %%

'''
1.2) MULTI PLOT: Effect of input intensity and duration on firing rates
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    - subplot columns: populations
    - subplot rows: coupling strengths
'''

rate_measure = 'longtermVSbaseline'
coupling_strengths = [0, 40, 80]
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

        minmax_df = summary_df[summary_df['coupling_strength']==g]
        minmax_df = minmax_df[minmax_df['population']==p]
        minmax_df['InputDuration'] = minmax_df['InputDuration'].round(4)

        data_heatmap = minmax_df.pivot(index='InputStrength',columns='InputDuration', values=rate_measure)

        if (minmax_df[rate_measure].isna() | (minmax_df[rate_measure] == 0)).all().all():
            print('black')
            sns.heatmap(data_heatmap, cmap=ListedColormap(['black']), ax=axes[i,j], norm = Normalize(vmin=0, vmax=1))
        else:
            sns.heatmap(data_heatmap, cmap='magma', ax=axes[i,j])

        cbar = axes[i, j].collections[0].colorbar
        # here set the labelsize by 20
        cbar.ax.tick_params(labelsize=12)
        axes[i, j].invert_yaxis()
        axes[i, j].set_ylabel('')
        axes[i, j].set_xlabel('')
        axes[i, j].tick_params(axis='both', labelsize=12)
        axes[len(coupling_strengths)-1, j].set_xlabel('Input Duration')
        axes[0,j].set_title(f'pop: {p}')
        axes[i,0].set_ylabel(f'G: {g}', rotation=0, labelpad=60)

fig.text(0.04, 0.2, 'Input Strength', va='center', rotation='vertical')
fig.text(0.04, 0.5, 'Input Strength', va='center', rotation='vertical')
fig.text(0.04, 0.83, 'Input Strength', va='center', rotation='vertical')

plt.tight_layout(h_pad=1)
figure_name = f'inputDurationVSinputStrength_{populations[0][0]}pop_{rate_measure}_tauVisual.png'
#plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()

# %% investigate one specific population and see how this behaves..
# .. during input
# .. on the long term behaviour
# .. is the long term behaviour over baseline? 
# .. frequency of oscillations




# %%

'''
1.3) Line plot
    - x axis: time
    - y axis: rate
'''

# For this I need a dataframe with all timepoints
coupling_strength = [0, 40, 80]
population = 'P4'
input_duration = 1
input_strength = 15
line_df = trajectory_df.loc[trajectory_df['coupling_strength'].isin(coupling_strength)]
line_df = line_df[line_df['population']==population]
line_df = line_df[line_df['InputDuration']==input_duration]
line_df = line_df[line_df['InputStrength']==input_strength]
plt.title(f'Population:{population} Input Duration:{input_duration} Input Strength{input_strength}')
sns.lineplot(data=line_df, x='time', y='rate', hue='coupling_strength')

figure_name = f'population{population}_inputDur{input_duration}_inputStrength{input_strength}.png'
#plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')


# %%  
'''
Fixed point plot
In this plot I use the potential instead of the firing rates.
    - x axis: input strengths 
    - y axis: response (steady-state/longterm) in mV
'''

mean_lt = potentials_df.groupby(['InputDuration', 'coupling_strength', 'InputStrength'], as_index=False).mean()

# For this I need a dataframe with all timepoints
coupling_strength = 80
population = 1
input_duration = 0.5
input_strength = np.arange(10, 50, 5)
line_df = mean_lt[mean_lt['coupling_strength']==coupling_strength]
line_df = line_df[line_df['population']==population]
line_df = line_df[line_df['InputDuration']==input_duration]
plt.title(f'Population:{population} Input Duration:{input_duration} Input Strength{input_strength}')
sns.lineplot(data=line_df, x='InputStrength', y='lt_potential') # , hue='coupling_strength')

figure_name = f'fixedPointCurve_population{population}_inputDur{input_duration}_G{coupling_strength}.png'
#plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')

# %%
populations = np.array(['E1', 'E2', 'E3', 'E4', 'P1', 'P2', 'P3', 'P4', 'S1', 'S2', 'S3', 'S4', 'V1'])#, 'V2', 'V3', 'V4'])

fig, axes = plt.subplots(4, 4 ,sharex=True, sharey=False)
coupling_strength = 40
for ax, pop in zip(axes.flatten(), populations):
    line_df = trajectory_df[trajectory_df['coupling_strength']==coupling_strength]
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



# %% 

