"""
File: compare_simulations.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-01-05
Description: I want to numerically compare the result of the pyrates simulation
 and the simulation with pure python. 

"""

# %%
import numpy as np 
import pandas as pd
from helper_functions import *

# Define paths
SIMDIR = os.getenv("SIMDIR")
sim_dir = os.path.join(SIMDIR, "simulation_results")
WDDIR = os.getenv("WDDIR")

# define output directory
output_dir = os.path.join(SIMDIR, "derivatives")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# %%
# Load the data 

params = load_parameters(WDDIR)
g = 10.0
bEI = 0.5
bI = 5
d = 0
s = 0
input_onset = 1.001
thal_cellcounts = 500
bI_cellcounts = 100
extI_cellcounts = 1000
input_type = 'step'

rates_df, potentials_df, _ = load_simulation_data(g, bEI, bI, d, s, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, sim_dir)
rates_df_pyrates, potentials_df_pyrates, _ = load_simulation_data(g, bEI, bI, d, s, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, sim_dir, pyrates=True)


# %%

diff_df = pd.DataFrame(
    rates_df.to_numpy() - rates_df_pyrates.to_numpy(),
    index=rates_df_pyrates.index,
    columns=rates_df_pyrates.columns
)

# %%
plt.plot(diff_class_grid0['E3b'], label='E3b')
plt.plot(diff_class_grid0['E1S1'], label='E1S1')
plt.plot(diff_class_grid0['PV2S1'], label='PV2S1')

plt.legend()