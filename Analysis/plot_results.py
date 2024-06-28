'''
Plots:
1) Effect of input intensity and duration on firing rates
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    1.1) SINGLE PLOT
    1.2) MULTI PLOT
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
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from plotting_style import figure_style

colors, _ = figure_style() 

# define directories of stored data and figures
output_dir = "../output"
figure_dir = "../Figures"

# %%

# read in data
input_durations = np.arange(0, 0.04, 0.005) # in sec 
input_strengths = [0, 5, 10, 15, 20] # np.arange(0, 20, 2)
coupling_strengths = [20] # np.arange(0, 100, 10)
step_size = 0.001
sample_delay = 0.5
input_onset = 0.501
sample_dur = 0.1 # amount of time in sec in which we look at the long term firing rate (min and max)
cortex_type = 'somato'
load_trajectory = False

# create summary matrix: should include (max/min_value x population   
min_data = []
max_data = []
summary_df = pd.DataFrame()
time_data_df = pd.DataFrame()

for d in input_durations:
    for s in input_strengths:
        for g in coupling_strengths:
            df = pd.DataFrame()

            # read data matrix (datapoints x populations)
            filename = f"rates_G{g}_{cortex_type}_Iduration{d}_Istrength{s}_Ionset{input_onset}.csv"
            data_df = pd.read_csv(os.path.join(output_dir, filename))

            # LONG TERM (sample for 100ms starting at 0.5 sec after offset)
            start_sample = int((input_onset+d+sample_delay)/step_size)
            stop_sample = int(start_sample + sample_dur/step_size)
            df['minRate_longterm'] = data_df.iloc[start_sample:stop_sample].min()
            df['maxRate_longterm'] = data_df.iloc[start_sample:stop_sample].max()
            df['diffRate_longterm'] = df['maxRate_longterm'] - df['minRate_longterm']
            
            start_sample_during = int((input_onset)/step_size)
            stop_sample_during = int((input_onset+d)/step_size)
            df['minRate_duringInput'] = data_df.iloc[start_sample_during:stop_sample_during].min()
            df['maxRate_duringInput'] = data_df.iloc[start_sample_during:stop_sample_during].max()
            df['diffRate_duringInput'] = df['maxRate_duringInput'] - df['minRate_duringInput']
            df['population'] = df.index.values
            df['coupling_strength'] = g
            df['InputDuration'] = d
            df['InputStrength'] = s
            summary_df = pd.concat([summary_df, df])

            if load_trajectory:
                # load entire trajectory in data frame
                for pop_name, pop in data_df.items():
                    #print(pop_name)
                    pop_df = pd.DataFrame()
                    pop_df['population'] = [pop_name]*len(pop)
                    pop_df['rate'] = pop
                    pop_df['time'] = pop.index.values * 1e-3 # time in s
                    pop_df['coupling_strength'] = g
                    pop_df['InputDuration'] = d
                    pop_df['InputStrength'] = s 
                    time_data_df = pd.concat([time_data_df, pop_df])

# %% Make plots that demonstrate the sampling time line 

# choose example settings
coupling_strength = 50
population = 'S4'
input_duration = 0.1
input_strength = 8
line_df = time_data_df[time_data_df['coupling_strength']==coupling_strength]
line_df = line_df[line_df['population']==population]
line_df = line_df[line_df['InputDuration']==input_duration]
line_df = line_df[line_df['InputStrength']==input_strength]
# LONG TERM and DURING INPUT

# Add a red vertical line at time of input offset (start of sampling) and stop of sampling
simulation_time = len(line_df)*1e-3
start_sample = input_onset+input_duration+sample_delay
stop_sample = start_sample + sample_dur
input_offset = input_onset+input_duration

# plot the input
steps = np.arange(step_size, simulation_time+step_size, step_size)
input_line = np.zeros(len(line_df))
input_line[int((input_onset)/step_size):int(input_offset/step_size)] = input_strength
plt.plot(steps, input_line, color='green', linewidth=2)
sns.lineplot(data=line_df, x='time', y='rate', hue='InputStrength', legend='', palette=['grey'])
plt.axvline(x=start_sample, color='red', linestyle='--', linewidth=1, label='Long Term Sample')
plt.axvline(x=stop_sample, color='red', linestyle='--', linewidth=1, label='')
plt.axvline(x=input_onset, color='blue', linestyle='--', linewidth=1, label='Sample During Input')
plt.axvline(x=input_offset, color='blue', linestyle='--', linewidth=1, label='')
plt.legend()
#plt.savefig('../Figures/plotting_windows.pdf')
plt.show()

# DURING INPUT

# %% OSCIALLTIONS

# Find the combination in which the output is oscillatory
# For this plot the difference values 


# %%

'''
1.1) SINGLE PLOT: Effect of input intensity and duration on firing rates
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
'''

# choose a coupling strength and a population
coupling_strength = 20
population = 'E1'
data_df = summary_df[summary_df['coupling_strength']==coupling_strength]
data_df = data_df[data_df['population']==population]
data_heatmap = data_df.pivot(index='InputDuration',columns='InputStrength', values='minRate_duringInput')
sns.heatmap(data_heatmap, cmap='magma')

# %%

'''
1.2) MULTI PLOT: Effect of input intensity and duration on firing rates
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    - subplot columns: populations
    - subplot rows: coupling strengths
'''

rate_measure = 'diffRate_duringInput'
coupling_strengths = [20, 0]
populations = np.array(['E1', 'E2', 'E3', 'E4']) #, 'P1', 'P2', 'P3', 'P4', 'S1', 'S2', 'S3', 'S4', 'V1'])#, 'V2', 'V3', 'V4'])

fig, axes = plt.subplots(len(coupling_strengths), len(populations), figsize=(20,15) ,sharex=True, sharey=True)

# Create a single colorbar axis
#cbar_ax = fig.add_axes([1.01, 0.3, 0.02, 0.4])
#cbar_ax.set_title(rate_measure)
#cbar_ax.tick_params(labelsize=12) 

for i,g in enumerate(coupling_strengths):
    for j,p in enumerate(populations):

        minmax_df = summary_df[summary_df['coupling_strength']==g]
        minmax_df = minmax_df[minmax_df['population']==p]

        data_heatmap = minmax_df.pivot(index='InputStrength',columns='InputDuration', values=rate_measure)
        sns.heatmap(data_heatmap, cmap='magma', ax=axes[i,j])
                    #, cbar=(i == 0 and j == len(populations) - 1),
                    #cbar_ax=None if (i != 0 or j != len(populations) - 1) else cbar_ax)
        
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
figure_name = f'inputDurationVSinputStrength_{populations[0][0]}pop_{rate_measure}.png'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()
# %%

'''
1.3) Line plot
    - x axis: time
    - y axis: rate
'''

# For this I need a dataframe with all timepoints
coupling_strength = 10
population = 'E1'
input_duration = 0.1
line_df = time_data_df[time_data_df['coupling_strength']==coupling_strength]
line_df = line_df[line_df['population']==population]
line_df = line_df[line_df['InputDuration']==input_duration]
#line_df = line_df[line_df['InputStrength']==0]

sns.lineplot(data=line_df, x='time', y='rate', hue='InputStrength')

# %%
populations = np.array(['E1', 'E2', 'E3', 'E4', 'P1', 'P2', 'P3', 'P4', 'S1', 'S2', 'S3', 'S4', 'V1'])#, 'V2', 'V3', 'V4'])

fig, axes = plt.subplots(4, 4 ,sharex=True, sharey=False)

for ax, pop in zip(axes.flatten(), populations):
    line_df = time_data_df[time_data_df['coupling_strength']==coupling_strength]
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