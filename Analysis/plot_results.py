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
from matplotlib.ticker import FormatStrFormatter, FuncFormatter, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from scipy.signal import find_peaks
import seaborn as sns
import pandas as pd
from plotting_style import figure_style

colors, _ = figure_style() 

# define directories of stored data and figures
output_dir = "../output"
figure_dir = "../Figures"

# %%

# read in data
input_durations = np.arange(0.15, 0.6, 0.05) #  np.arange(0, 0.04, 0.001) # in sec 
input_strengths = np.arange(0, 20, 5)
coupling_strengths = np.arange(0, 100, 40)
step_size = 0.001
sample_delay = 0.5
input_onset = 0.501
sample_dur = 0.1 # amount of time in sec in which we look at the long term firing rate (min and max)
cortex_type = 'somato'
load_trajectory = True

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
            filename = f"rates_G{g}_{cortex_type}_Iduration{d}_Istrength{s}_Ionset{input_onset}_tauVisual.csv"
            data_df = pd.read_csv(os.path.join(output_dir, filename))

            # LONG TERM (sample for 100ms starting at 0.5 sec after offset)
            start_sample = int((input_onset+d+sample_delay)/step_size)
            stop_sample = int(start_sample + sample_dur/step_size)
            df['rate_longterm'] = data_df.iloc[start_sample:stop_sample].mean()
            df['minRate_longterm'] = data_df.iloc[start_sample:stop_sample].min()
            df['maxRate_longterm'] = data_df.iloc[start_sample:stop_sample].max()
            df['diffRate_longterm'] = df['maxRate_longterm'] - df['minRate_longterm']
            
            # DURING INPUT
            start_sample_during = int((input_onset)/step_size)
            stop_sample_during = int((input_onset+d)/step_size)
            df['minRate_duringInput'] = data_df.iloc[start_sample_during:stop_sample_during].min()
            df['maxRate_duringInput'] = data_df.iloc[start_sample_during:stop_sample_during].max()
            df['diffRate_duringInput'] = df['maxRate_duringInput'] - df['minRate_duringInput']
        
            # BASELINE SAMPLE
            offset = 0.1 # time between baseline sampling and start of input 
            baseline_start = int((input_onset - (sample_dur+offset))/step_size) 
            baseline_stop = int(baseline_start + sample_dur/step_size)
            df['rate_baseline'] = data_df.iloc[baseline_start:baseline_stop].mean()
            
            # long term to baseline comparison (here we don't take the diffRate because the baseline does not 
            # oscillate and we want to compare if the long term activity is still larger than baseline)
            df['longtermVSbaseline'] = df['rate_longterm'] - df['rate_baseline']
            
            # set info in data frame
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
coupling_strength = 40
population = 'E1'
input_duration = 0.12
input_strength = 4
line_df = time_data_df[time_data_df['coupling_strength']==coupling_strength]
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

# plot the input
steps = np.arange(step_size, simulation_time+step_size, step_size)
input_line = np.zeros(len(line_df))
input_line[int((input_onset)/step_size):int(input_offset/step_size)] = input_strength
plt.plot(steps, input_line, color='green', linewidth=2)
sns.lineplot(data=line_df, x='time', y='rate', hue='InputStrength', legend='', palette=['grey'])
plt.axvline(x=baseline_start, color='purple', linestyle='--', linewidth=1, label='Baseline Sample')
plt.axvline(x=baseline_stop, color='purple', linestyle='--', linewidth=1, label='')
plt.axvline(x=start_sample, color='red', linestyle='--', linewidth=1, label='Long Term Sample')
plt.axvline(x=stop_sample, color='red', linestyle='--', linewidth=1, label='')
plt.axvline(x=input_onset, color='blue', linestyle='--', linewidth=1, label='During Input Sample')
plt.axvline(x=input_offset, color='blue', linestyle='--', linewidth=1, label='')
plt.legend()
#plt.savefig('../Figures/plotting_windows.pdf')
plt.show()

# %% OSCIALLTIONS

# Find the combination in which the output is oscillatory
# For this plot the difference values 

rate_measure = 'diffRate_duringInput'
coupling_strengths = [0, 20, 40]
input_duration = 0.14
input_strengths = [10, 12, 14]
populations = np.array(['E1'])#,  'E2', 'E3', 'E4']) 
#populations = np.array(['P1', 'P2', 'P3', 'P4']) 
#populations = np.array(['S1', 'S2', 'S3', 'S4', 'V1']) 

for k, s in enumerate(input_strengths):
    for i,g in enumerate(coupling_strengths):
        for j,p in enumerate(populations):

            # get the summary data frame
            waves_df = summary_df[summary_df['coupling_strength']==g]
            waves_df_strength = minmax_df[minmax_df['input_strength']==s]

            # plot the trajctory
            plt.plot()

            # compute the amount of peaks


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
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
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
population = 'P1'
input_duration = 0.15
input_strength = 15
line_df = time_data_df.loc[time_data_df['coupling_strength'].isin(coupling_strength)]
line_df = line_df[line_df['population']==population]
line_df = line_df[line_df['InputDuration']==input_duration]
line_df = line_df[line_df['InputStrength']==input_strength]
plt.title(f'Population:{population} Input Duration:{input_duration} Input Strength{input_strength}')
sns.lineplot(data=line_df, x='time', y='rate', hue='coupling_strength')

figure_name = f'population{population}_inputDur{input_duration}_inputStrength{input_strength}.png'
#plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')

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

# %% Plot the difference to the baseline

