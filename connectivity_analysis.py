# %% 
from parameters import Parameter
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from plotting_style import figure_style

colors, _ = figure_style() 

def adjust_matrix(matrix):
    matrix = np.insert(matrix, 7, 0, axis=1)
    matrix = np.insert(matrix, 7, 0, axis=0)
    matrix = np.insert(matrix, 11, 0, axis=1)
    matrix = np.insert(matrix, 11, 0, axis=0)
    matrix = np.append(matrix, np.zeros((15, 1)), axis=1)
    matrix = np.append(matrix, np.zeros((1, 16)), axis=0)

    return matrix

def adjust_connectivity_matrix(conmatrix):
    conmatrix = np.append(conmatrix, np.zeros((13, 1)), axis=1)
    conmatrix = np.append(conmatrix, np.zeros((1, 14)), axis=0)
    conmatrix = np.append(conmatrix, np.zeros((14, 1)), axis=1)
    conmatrix = np.append(conmatrix, np.zeros((1, 15)), axis=0)
    conmatrix = np.append(conmatrix, np.zeros((15, 1)), axis=1)
    conmatrix = np.append(conmatrix, np.zeros((1, 16)), axis=0)
    return conmatrix

# Plot difference between Visual and Somatosensory connectivity
params_visual = Parameter('visual')
params_somato = Parameter('somato')

# Synaptic strengths
visual_S = params_visual.get_connectStrength()
somato_S = np.abs(params_somato.get_connectStrength())

# make the somatosensory matrix the same size as the visual one (where the VIP cells are insert a row/column of zeros)
somato_S = adjust_matrix(somato_S)

# %%
population_names = ['E1','PV1','SST1','VIP1','E2','PV2','SST2','VIP2','E3','PV3','SST3','VIP3','E4','PV4','SST4','VIP4']
plot_populations = ['E1','E2','E3','E4','PV1','PV2','PV3','PV4','SST1','SST2','SST3','SST4','VIP1','VIP2','VIP3','VIP4']
df_somato_S = pd.DataFrame(somato_S, index=population_names, columns=population_names)
df_somato_S['cortex_type'] = 'somato'
df_visual_S = pd.DataFrame(visual_S, index=population_names, columns=population_names)
df_visual_S['cortex_type'] = 'visual' 
df_S = pd.concat((df_somato_S, df_visual_S), ignore_index=False)

# %% Plot Synaptic strengths

fig, axs = plt.subplots(4, 4, figsize=(10,10), sharey=False)

for i, ax in enumerate(axs.flatten()):
    y = visual_S[i]
    y2 = somato_S[i]
    sns.barplot(data=df_S, ax=ax, x=population_names*2, y=plot_populations[i], hue='cortex_type', palette=[colors['somato'], colors['visual']])
    ax.set_ylabel("Synaptic Strength")
    ax.set_title(f'from {plot_populations[i]}')
    ax.tick_params(axis='x', labelsize=8, rotation=60)

sns.despine(fig, trim=True, bottom=False)
plt.tight_layout(h_pad=1)
plt.legend(title='cortex type')
#plt.savefig("Figures/Synaptic_strength.pdf")
plt.show()

# %%
# Connectivity Probabilities

visual_P = params_visual.get_connectProb()
somato_P = params_somato.get_connectProb()
somato_P = adjust_matrix(somato_P)
df_somato_P = pd.DataFrame(somato_P, index=population_names, columns=population_names)
df_somato_P['cortex_type'] = 'somato'
df_visual_P = pd.DataFrame(visual_P, index=population_names, columns=population_names)
df_visual_P['cortex_type'] = 'visual' 
df_P = pd.concat((df_somato_P, df_visual_P), ignore_index=False)

# %% PLOT
fig, axs = plt.subplots(4, 4, figsize=(10,10), sharey=False)
for i, ax in enumerate(axs.flatten()):
    sns.barplot(data=df_P, ax=ax, x=population_names*2, y=plot_populations[i], hue='cortex_type', palette=[colors['somato'], colors['visual']])
    ax.set_ylabel("Connection Probabilities")
    ax.set_title(f'from {plot_populations[i]}')
    legend = ax.legend(title='cortex type', facecolor='white', edgecolor='black', framealpha=1)
    if not i%4 == 0:
        ax.get_legend().remove()
    ax.set_xlabel('')
    ax.tick_params(axis='x', labelsize=8, rotation=60)

sns.despine(fig, trim=True, bottom=False)
plt.tight_layout(h_pad=1)
#plt.savefig("Figures/Connect_Probability.pdf")
plt.show()

# %% 
sns.heatmap(somato_S[:,:], annot=True, cmap='magma',xticklabels=population_names, yticklabels=population_names)


# %%
# Cell Counts
somato_C = params_somato.get_cellcounts()
somato_C = pd.DataFrame(np.insert(somato_C, [7, 10, 13], 0), index=population_names, columns=['cellcount'])
somato_C['cortex_type'] = 'somato'
visual_C = pd.DataFrame(params_visual.get_cellcounts(), index=population_names, columns=['cellcount'])
visual_C['cortex_type'] = 'visual'

df_C = pd.concat((somato_C, visual_C))

# %%
fig, ax = plt.subplots()
sns.barplot(df_C, x=df_C.index, y='cellcount', hue='cortex_type', ax=ax, palette=[colors['somato'], colors['visual']], dodge=True)
sns.despine(fig, trim=True, bottom=False)
ax.set_ylabel('Proportional Cellcount')
ax.set_xlabel('Populations')
#plt.legend('cortex type')
plt.show()

# %% 
# Connectivity matrix

visual_W = pd.DataFrame(np.abs(params_visual.get_connectivity(1, False)), index=plot_populations, columns=plot_populations)
somato_W_matrix = np.abs(params_somato.get_connectivity(1, False))
somato_W = pd.DataFrame(adjust_connectivity_matrix(somato_W_matrix), index=plot_populations, columns=plot_populations)
W_diff = np.array(somato_W) - np.array(visual_W)
somato_W['cortex_type'] = 'somato'
visual_W['cortex_type'] = 'visual' 
df_W = pd.concat((visual_W, somato_W), ignore_index=False)
# %%
sns.heatmap(somato_W_matrix[:,:], annot=True, cmap='magma',xticklabels=plot_populations, yticklabels=plot_populations)

# %% 

fig, axs = plt.subplots(4, 4, figsize=(10,10), sharey=False)
x = population_names
for i, ax in enumerate(axs.flatten()):
    sns.barplot(data=df_W, ax=ax, x=plot_populations*2, y=plot_populations[i], hue='cortex_type', palette=[colors['visual'], colors['somato']])

    #sns.barplot(x=x, y=y, ax=ax, color=colors['visual'])
    #sns.barplot(x=x, y=y2, ax=ax, color=colors['somato'])
    ax.set_ylabel("Connection weight")
    ax.set_title(f'from {plot_populations[i]}')
    legend = ax.legend(title='cortex type', facecolor='white', edgecolor='black', framealpha=1)
    if not i%4 == 0:
        ax.get_legend().remove()
    ax.set_xlabel('')
    ax.tick_params(axis='x', labelsize=8, rotation=60)

sns.despine(fig, trim=True, bottom=False)
plt.tight_layout()
#plt.savefig('Figures/Connection_weight.pdf')
plt.show()

# %% Input matrix (connections from the thalamus)

