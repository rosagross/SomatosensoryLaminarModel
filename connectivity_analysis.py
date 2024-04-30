# %% 
from parameters import Parameter
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from plotting import figure_style

colors, _ = figure_style() 

def adjust_matrix(matrix):
    matrix = np.insert(matrix, 7, 0, axis=1)
    matrix = np.insert(matrix, 7, 0, axis=0)
    matrix = np.insert(matrix, 11, 0, axis=1)
    matrix = np.insert(matrix, 11, 0, axis=0)
    matrix = np.append(matrix, np.zeros((15, 1)), axis=1)
    matrix = np.append(matrix, np.zeros((1, 16)), axis=0)

    return matrix

# Plot difference between Visual and Somatosensory connectivity
params = Parameter()

# Synaptic strengths
visual_S = params.get_connectStrength('visual')
somato_S = np.abs(params.get_connectStrength('somato'))

# make the somatosensory matrix the same size as the visual one (where the VIP cells are insert a row/column of zeros)
somato_S = adjust_matrix(somato_S)
# %%
population_names = ['E1','PV1','SST1','VIP1','E2','PV2','SST2','VIP2','E3','PV3','SST3','VIP3','E4','PV4','SST4','VIP4']
df_somato_S = pd.DataFrame(somato_S, index=population_names, columns=population_names)
df_somato_S['cortex_type'] = 'somato'
df_visual_S = pd.DataFrame(visual_S, index=population_names, columns=population_names)
df_visual_S['cortex_type'] = 'visual' 
df_S = pd.concat((df_somato_S, df_visual_S), ignore_index=False)

# %% 
fig, axs = plt.subplots(4, 4, figsize=(10,10), sharey=False)
x = population_names
for i, ax in enumerate(axs.flatten()):
    y = visual_S[i]
    y2 = somato_S[i]
    sns.barplot(data=df_S, ax=ax, x=df_S.index, y=population_names[i], hue='cortex_type', 
                palette=[colors['somato'], colors['visual']])
    #sns.barplot(x=x, y=y, ax=ax, color=colors['visual'])
    #sns.barplot(x=x, y=y2, ax=ax, color=colors['somato'])
    ax.set_ylabel("Synaptic Strength", fontsize=15)
    ax.set_title(f'from {population_names[i]}', size=15)

sns.despine(fig, trim=True, bottom=False)
plt.tight_layout(h_pad=1)
plt.legend('cortex type')
plt.show()

# %%
# Connectivity Probabilities

visual_P = params.get_connectProb('visual')
somato_P = params.get_connectProb('somato')
somato_P = adjust_matrix(somato_P)
df_somato_P = pd.DataFrame(somato_P, index=population_names, columns=population_names)
df_somato_P['cortex_type'] = 'somato'
df_visual_P = pd.DataFrame(visual_P, index=population_names, columns=population_names)
df_visual_P['cortex_type'] = 'visual' 
df_P = pd.concat((df_somato_P, df_visual_P), ignore_index=False)

# %% PLOT
fig, axs = plt.subplots(4, 4, figsize=(10,10), sharey=False)
for i, ax in enumerate(axs.flatten()):
    sns.barplot(data=df_P, ax=ax, x=df_P.index, y=population_names[i], hue='cortex_type', palette=[colors['somato'], colors['visual']])
    ax.set_ylabel("Connection Probabilities")
    ax.set_title(f'from {population_names[i]}')
    ax.get_legend().remove()
    ax.set_xlabel('Populations')
axs[0,0].legend(title='cortex type')
sns.despine(fig, trim=True, bottom=False)
plt.tight_layout(h_pad=1)
plt.show()

# %%
# Cell Counts
somato_C = params.get_cellcounts('somato')
somato_C = pd.DataFrame(np.insert(somato_C, [7, 10, 13], 0), index=population_names, columns=['cellcount'])
somato_C['cortex_type'] = 'somato'
visual_C = pd.DataFrame(params.get_cellcounts('visual'), index=population_names, columns=['cellcount'])
visual_C['cortex_type'] = 'visual'

df_C = pd.concat((somato_C, visual_C))

# %%
fig, ax = plt.subplots()
sns.barplot(df_C, x=df_C.index, y='cellcount', hue='cortex_type', ax=ax, palette=[colors['somato'], colors['visual']], dodge=True)
sns.despine(fig, trim=True, bottom=False)
ax.set_ylabel('Proportional Cellcount')
ax.set_xlabel('Populations')
plt.legend('cortex type')
plt.show()

# %% 
# Connectivity matrix

visual_W = pd.DataFrame(params.get_connectivity(1, 'visual', False), index=population_names, columns=population_names)
somato_W = adjust_matrix(params.get_connectivity(1, 'somato', False))
somato_W = pd.DataFrame(somato_W, index=population_names, columns=population_names)
somato_W['cortex_type'] = 'somato'
visual_W['cortex_type'] = 'visual' 
df_W = pd.concat((visual_W, somato_W), ignore_index=False)


# %% 

fig, axs = plt.subplots(4, 4, figsize=(10,10), sharey=False)
x = population_names
for i, ax in enumerate(axs.flatten()):
    sns.barplot(data=df_W, ax=ax, x=df_S.index, y=population_names[i], hue='cortex_type', palette=[colors['visual'], colors['somato']])
    #sns.barplot(x=x, y=y, ax=ax, color=colors['visual'])
    #sns.barplot(x=x, y=y2, ax=ax, color=colors['somato'])
    ax.set_ylabel("Connection weight")
    ax.set_title(f'from {population_names[i]}')

sns.despine(fig, trim=True, bottom=False)
plt.tight_layout()
plt.legend('cortex type')
plt.show()

# %%