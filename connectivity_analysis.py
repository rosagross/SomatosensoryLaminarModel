# %% 
from parameters import Parameter
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from plotting_style import figure_style
import matplotlib.colors as clrs

# %%

colors, _ = figure_style() 

def adjust_matrix(matrix):
    # remove the thalamus values (we look at those later)
    matrix = np.insert(matrix, 7, 0, axis=1)
    matrix = np.insert(matrix, 7, 0, axis=0)
    matrix = np.insert(matrix, 11, 0, axis=1)
    matrix = np.insert(matrix, 11, 0, axis=0)
    matrix = np.append(matrix, np.zeros((15, 1)), axis=1)
    matrix = np.append(matrix, np.zeros((1, 16)), axis=0)
    print(matrix.shape)

    return matrix


def adjust_connectivity_matrix(conmatrix):
    conmatrix = np.append(conmatrix, np.zeros((13, 1)), axis=1)
    conmatrix = np.append(conmatrix, np.zeros((1, 14)), axis=0)
    conmatrix = np.append(conmatrix, np.zeros((14, 1)), axis=1)
    conmatrix = np.append(conmatrix, np.zeros((1, 15)), axis=0)
    conmatrix = np.append(conmatrix, np.zeros((15, 1)), axis=1)
    conmatrix = np.append(conmatrix, np.zeros((1, 16)), axis=0)
    return conmatrix

# Set the directory where to save the figures
figure_dir = 'C:/Users/gross/OneDrive - UvA/Documents/IMPRS_Leipzig/IMPRS SummerSchool/Poster/PosterFigures'

params = Parameter()
population_names = ['E3b','PV3b','SST3b', 'E1','PV1','SST1','VIP1','E2','PV2','SST2','E3','PV3','SST3','E4','PV4','SST4',
                                          'E1S2','PV1S2','SST1S2','VIP1S2','E2S2','PV2S2','SST2S2','E3S2','PV3S2','SST3S2','E4S2','PV4S2','SST4S2']
plot_populations = ['E3b','PV3b','SST3b', 'E1','E2','E3','E4','PV1','PV2','PV3','PV4','SST1','SST2','SST3','SST4','VIP1', 
                                          'E1S2','E2S2','E3S2','E4S2','PV1S2','PV2S2','PV3S2','PV4S2','SST1S2','SST2S2','SST3S2','SST4S2','VIP1S2']
all_pops = np.concatenate((population_names,['ThalE', 'ThalI', 'B', 'Ext']))
# Synaptic strengths
S = np.abs(params.get_connectStrength())
df_S = pd.DataFrame(S, index=population_names, columns=population_names)
# Connection Probability
P = params.get_connectProb()
df_P = pd.DataFrame(P, index=population_names, columns=population_names)
# Cell counts
C = params.get_cellcounts()
df_C = pd.DataFrame(C, index=population_names)
# Total Connectivity
W = params.get_connectivity(60, 40)
df_W = pd.DataFrame(W, index=all_pops[:-2], columns=all_pops)

# %% Plot Synaptic strengths
plt.figure(figsize=(12, 10))
sns.heatmap(df_S, annot=False, cmap='coolwarm', center=0, xticklabels=True, yticklabels=True)

# Rotate x-axis labels for better visibility
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.ylabel('Target Population')
plt.xlabel('Source Population')
plt.title("Synaptic Strength")
plt.tight_layout() 
#plt.savefig("Figures/Synaptic_strength.pdf")
plt.show()

# %%
# Connectivity Probabilities
plt.figure(figsize=(12, 10))
sns.heatmap(df_W, annot=False, cmap='coolwarm', center=0, xticklabels=True, yticklabels=True)

# Rotate x-axis labels for better visibility
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.ylabel('Target Population')
plt.xlabel('Source Population')
plt.title("Connection Probability")
plt.tight_layout() 
#plt.savefig("Figures/Connect_Probability.pdf")
plt.show()

# %%
# Plot Cell Counts
fig, ax = plt.subplots()
sns.barplot(df_C.T, ax=ax)
sns.despine(fig, trim=True, bottom=False)
ax.set_ylabel('Proportional Cellcount')
ax.set_xlabel('Populations')
plt.xticks(rotation=45, ha='right')
plt.tight_layout() 
#plt.legend('cortex type')
plt.show()

# %% 
################ Connectivity matrix #############
plt.figure(figsize=(12, 10))
sns.heatmap(df_W, annot=False, cmap='vlag', center=0, xticklabels=True, yticklabels=True)

# Rotate x-axis labels for better visibility
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.ylabel('Target Population')
plt.xlabel('Source Population')
plt.title("Connectivity Matrix")
plt.tight_layout() 
#plt.savefig("Figures/Connect_Probability.pdf")
plt.show()

#%% 
################ Thlamus connectivity ##################
plt.figure(figsize=(5,12))
sns.heatmap(df_W[["ThalE", 'ThalI']], annot=False, cmap='vlag', center=0, xticklabels=True, yticklabels=True)

#%% 
################ Background & External input ##################
plt.figure(figsize=(5,12))
sns.heatmap(df_W[["B", 'Ext']], annot=False, cmap='coolwarm', xticklabels=True, yticklabels=True)

# %%
