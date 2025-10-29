# %% imports
from sim_meg import *
import numpy as np
import seaborn as sns

# %%

coupling_strength = 10
coupling_strength_solution = 9.354
input_duration = 0.5
input_strength = 10
input_strength_solution = 10
bEI = 0.7
bEI_solution = 0.4
Ib = 7
input_onset = 1.001
filename = f'full_g{coupling_strength}_bEI{bEI}_Ib{Ib}_Iextd{input_duration}_stepIexts{input_strength}_Ionset{input_onset}_thalcells500_Ibcells100_Iextcells1000_thalUncon_S1S2Uncon.hdf5'
filename_solution = f'full_g{coupling_strength_solution}_bEI{bEI_solution}_Ib{Ib}_Iextd{input_duration}_stepIexts{input_strength_solution}_Ionset{input_onset}_thalcells500_Ibcells100_Iextcells1000_thalUncon_S1S2Uncon.hdf5'
SIMDIR = os.getenv("SIMDIR")
WDDIR = os.getenv("WDDIR")
figure_dir = os.path.join(SIMDIR, "Figures")
data_dir = os.path.join(SIMDIR, 'simulation_results')

# E-> E, PV-> E, SOM-> E, VIP-> E, Th-> E
dipole_setting = [-1, -1, 1, 1, 1]
nCells = [13, 4, 3, 3, 3]
indices_E = [4, 8, 11, 15]
# TODO: insert factor that scales the resistance for the populations
# since we use relative cell counts (arbitrary unit), we need to fit the strength of the resistance

sim_dip = simDipoles(dipole_setting, nCells, indices_E, data_dir, filename)
sim_dip_solution = simDipoles(dipole_setting, nCells, indices_E, data_dir, filename_solution)


# %% PLOT DIPOLE
    
# Average of all dipoles [1 x timepoints], it's the dipole direction * currents
simMEG = np.sum(sim_dip, axis=0)
simMEG_solution = np.sum(sim_dip_solution, axis=0)
pd.DataFrame(simMEG).to_csv(os.path.join(WDDIR, 'Optimization', 'optimization_reference_same.csv'), index=False)
simMEG_reference = pd.read_csv(os.path.join(WDDIR, 'Optimization', 'optimization_reference_same.csv')).to_numpy()
time_step = 0.001
plot_window = 0.2
prestim_plot = 0.05
start = int((input_onset-prestim_plot)/time_step)
stop = int(start + plot_window*1e3)
steps = np.arange(-prestim_plot, plot_window-prestim_plot, time_step)*1e3
fig, axes = plt.subplots(1,3, sharex=True, sharey=True)

input_names = ['E layer 2/3', 'E layer 4', 'E layer 5', 'E layer 6']
for i, names in enumerate(input_names):
    axes[0].plot(steps, sim_dip[i, start:stop], label=names)
axes[0].axvline(x=0, color='grey', linestyle='--', linewidth=1)
axes[1].axvline(x=0, color='grey', linestyle='--', linewidth=1)
axes[0].legend(loc='center')
axes[1].plot(steps, simMEG[start:stop], label='sum')
axes[2].plot(steps, simMEG_reference[start:stop], label='target')
axes[2].plot(steps, simMEG_solution[start:stop], label='opt result')
axes[0].set_xlabel('Time (ms)')
axes[1].set_xlabel('Time (ms)')
axes[2].set_xlabel('Time (ms)')
axes[0].set_ylabel('ECD')
axes[1].set_ylabel('')
axes[2].set_ylabel('')
axes[1].legend(loc='center')
axes[2].legend(loc='center')
sns.despine(trim=True)
figure_name = f'dipole_simulation_G{coupling_strength}_Istrength{input_strength}_Iduration{input_duration}.pdf'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()

# %%
