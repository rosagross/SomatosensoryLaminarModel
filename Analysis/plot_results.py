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
input_durations = np.arange(0, 0.2, 0.04) # in sec 
input_strengths = np.arange(0, 20, 4)
coupling_strengths = [0, 5, 10] 
step_size = 0.001
cortex_type = 'somato'

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
            filename = f"rates_G{g}_{cortex_type}_Iduration{d}_Istrength{s}.csv"
            data_df = pd.read_csv(os.path.join(output_dir, filename))

            # append series to summary df
            df['minRate'] = data_df.iloc[-100:].min()
            df['maxRate'] = data_df.iloc[-100:].max()
            df['population'] = df.index.values
            df['coupling_strength'] = g
            df['InputDuration'] = d
            df['InputStrength'] = s
            summary_df = pd.concat([summary_df, df])

            #data_df['time'] = data_df.index.values * 1e-3 # time in s
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
                
# %%

'''
1.1) SINGLE PLOT: Effect of input intensity and duration on firing rates
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
'''

# choose a coupling strength and a population
coupling_strength = 5
population = 'E1'
data_df = summary_df[summary_df['coupling_strength']==coupling_strength]
data_df = data_df[data_df['population']==population]
data_heatmap = data_df.pivot(index='InputDuration',columns='InputStrength', values='min_val')
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

rate_measure = 'maxRate'
coupling_strengths = [0, 5, 10]
populations = np.array(['P1', 'P2', 'P3', 'P4'])#, 'P1', 'P2', 'P3', 'P4', 'S1', 'S2', 'S3', 'S4', 'V1'])#, 'V2', 'V3', 'V4'])

fig, axes = plt.subplots(len(coupling_strengths), len(populations), figsize=(20,15) ,sharex=True, sharey=True)

# Create a single colorbar axis
cbar_ax = fig.add_axes([1.01, 0.3, 0.02, 0.4])
cbar_ax.set_title(rate_measure)
cbar_ax.tick_params(labelsize=12) 

for i,g in enumerate(coupling_strengths):
    for j,p in enumerate(populations):

        minmax_df = summary_df[summary_df['coupling_strength']==g]
        minmax_df = minmax_df[minmax_df['population']==p]

        data_heatmap = minmax_df.pivot(index='InputStrength',columns='InputDuration', values=rate_measure)
        sns.heatmap(data_heatmap, cmap='magma', ax=axes[i,j])
                    #, cbar=(i == 0 and j == len(populations) - 1),
                    #cbar_ax=None if (i != 0 or j != len(populations) - 1) else cbar_ax)
        
        axes[i, j].set_ylabel('')
        axes[i, j].set_xlabel('')
        axes[i, j].tick_params(axis='both', labelsize=12)
        axes[len(coupling_strengths)-1, j].set_xlabel('Input Duration')
        axes[0,j].set_title(f'pop: {p}')
        axes[i,0].set_ylabel(f'G: {g}', rotation=0, labelpad=40)

fig.text(0, 0.5, 'Input Strength', va='center', rotation='vertical')

plt.tight_layout(h_pad=1)
figure_name = f'inputDurationVSinputStrength_Epop_{rate_measure}.png'
#plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()
# %%

'''
1.3) Line plot
    - x axis: time
    - y axis: rate
'''

# For this I need a dataframe with all timepoints
coupling_strength = 5
population = 'E1'
input_duration = 0
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
# %%
