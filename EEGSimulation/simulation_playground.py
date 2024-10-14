# %% imports
from sim_meg import *
import numpy as np
import seaborn as sns

# %%

coupling_strength = 40
input_duration = 0.019
input_strength = 20
input_onset = 1.001
filename = f'full_potentials_G{coupling_strength}_somato_Iduration{input_duration}_stepIstrength{input_strength}_Ionset{input_onset}_tauVisual_thalJiang.csv'
filename_visual = 'full_potentials_G20_visual_Iduration0.05_stepIstrength5_Ionset1.001_tauVisual_thalJiang.csv'
data_dir = '../output/'
figure_dir = 'C:/Users/gross/OneDrive - UvA/Documents/IMPRS_Leipzig/IMPRS SummerSchool/Poster/PosterFigures'

# E-> E, PV-> E, SOM-> E, VIP-> E, Th-> E
dipole_setting = [-1, -1, 1, 1, 1]
nCells = [13, 4, 3, 3, 3]
indices_E = [0, 4, 7, 11]
sim_dip = simDipoles(dipole_setting, nCells, indices_E, data_dir, filename)


# %% PLOT DIPOLE
    
# Average of all dipoles [1 x timepoints], it's the dipole direction * currents
simMEG = np.sum(sim_dip, axis=0)
time_step = 0.001
plot_window = 0.2
prestim_plot = 0.05
start = int((input_onset-prestim_plot)/time_step)
stop = int(start + plot_window*1e3)
steps = np.arange(-prestim_plot, plot_window-prestim_plot, time_step)*1e3
fig, axes = plt.subplots(1,2, sharex=True, sharey=True)

input_names = ['E layer 2/3', 'E layer 4', 'E layer 5', 'E layer 6']
for i, names in enumerate(input_names):
    axes[0].plot(steps, sim_dip[i, start:stop], label=names)
axes[0].axvline(x=0, color='grey', linestyle='--', linewidth=1)
axes[1].axvline(x=0, color='grey', linestyle='--', linewidth=1)
axes[0].legend(loc='lower right')
axes[1].plot(steps, simMEG[start:stop], label='averaged')
axes[0].set_xlabel('Time (ms)')
axes[1].set_xlabel('Time (ms)')
axes[0].set_ylabel('ECD')
axes[1].set_ylabel('')
axes[1].legend(loc='right')
sns.despine(trim=True)
figure_name = f'dipole_simulation_G{coupling_strength}_Istrength{input_strength}_Iduration{input_duration}.pdf'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()

# %%
