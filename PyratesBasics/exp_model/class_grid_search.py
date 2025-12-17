# %%
import os 
WDDIR = os.getenv("WDDIR")
SIMDIR = os.getenv("SIMDIR")
os.chdir(os.path.join(WDDIR,"PyratesBasics","exp_model")) 
import sys
from pyrates.frontend import OperatorTemplate, NodeTemplate, EdgeTemplate, CircuitTemplate
from pyrates import grid_search
from copy import deepcopy
from complete_model_class import SomatoModelPyrates
#from parameters import Parameter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numba import njit
#from yaml_saving import circuit_to_yaml
from pprint import pprint
## import dei parametri
param_path = os.path.join(WDDIR,"Simulations")

if param_path not in sys.path:
    sys.path.append(param_path)
from parameters import Parameter
import json
# %%  TODO: can I import this directly from the other code?
def read_simulation_params():
    """Read simulation parameters from json file."""
    # Read in preprocessing parameters
    with open(os.path.join(WDDIR, 'Simulations', 'simulation_parameter.json'), 'r') as json_file:
        params = json.load(json_file)
    
    return params
# %% model creation (class)
params = read_simulation_params()
somat_model = SomatoModelPyrates(params)
pmdl = somat_model.model # PyRates model

# %% Simulation parameters
simulation_time=somat_model.simulation_time
step_size=somat_model.step_size
sampling_step_size=somat_model.sampling_step_size

# %%
param_grid = {'background_input': np.array([0, 5, 7])} # values to "sweep" (defined as dictionary)  / the name is not related to the PyRates definition
param_map = {
    'background_input': {
        'vars': ['BACKGROUND/RPO_bI/bI'],
        'nodes': ['BACKGROUND']
    }
}

# %% 
cells = somat_model.cells
N_cells = len(cells)
rpo_names_extended = somat_model.rpo_names_extended
rpo_names = somat_model.rpo_names
outputs = {}
b_inputs = []
for i, target_cell in enumerate(cells):
    if i == (N_cells-2):
        for rpo_name_ext in rpo_names_extended[:N_cells+1]:
            outputs[f'V_{target_cell}/{rpo_name_ext}'] = f'{target_cell}/{rpo_name_ext}/v'
    else:
        for rpo_name in rpo_names[:N_cells]:
            outputs[f'V_{target_cell}/{rpo_name}'] = f'{target_cell}/{rpo_name}/v'

    # Include rpo_names_bI only if i in range(N_cells-2)
    #if i in range(N_cells-2):
        #print(f'{target_cell}/RPO_bI/v')
outputs['V_background/RPO_bI'] = 'BACKGROUND/RPO_bI/v_bI'

# %%
results, results_map = grid_search(circuit_template=pmdl,
                                   param_grid=param_grid,
                                   param_map=param_map,
                                   simulation_time=simulation_time,
                                   step_size=step_size,
                                   sampling_step_size=sampling_step_size,
                                   #inputs={'PC/RPO_e_in/u': noise},
                                   inputs = None, # the input is inside the model
                                   outputs=outputs,
                                   backend="scipy",
                                   permute_grid=True,
                                   vectorize=False,
                                   float_precision="float64"
                                   #solver="euler",
                                   #cutoff=cutoff
                                   )

# structure of "results": nested Pandas Dataframe with a multi-level column index
# 1) outputs requested in "grid_search"
# 2) simulation keys (parameters sweep runs)
# 3) node and variables of the model
# %%
"""
fig, axes = plt.subplots(nrows=results_map.shape[0], figsize=(8, 12))
results_map.sort_values('input_duration', inplace=True, axis=0)

for ax, key in zip(axes, results_map.index):
    all_potentials_param = []
    for i, cell in enumerate(cells):
    # always include rpo_names[:N_cells]
        if i == (N_cells-2):
            potential_keys = [f'V_{cell}/{rpo}' for rpo in rpo_names_extended[:N_cells+1]]
        else: 
            potential_keys = [f'V_{cell}/{rpo}' for rpo in rpo_names[:N_cells]]

    # include rpo_name only if i in range(N_cells-2)
        if i in range(N_cells-2):
            potential_keys += [f'V_{cell}/RPO_bI']

        idx = pd.IndexSlice
        sources = results.loc[:, idx[potential_keys, key]]
        all_potentials_param.append(np.sum(sources, axis=1))
"""

par_name1 = results_map.columns[0]
#par_name2 = results_map.columns[1]
layer_splits = [
    [i for i, c in enumerate(cells) if c.endswith("3b")],   # A3b
    [i for i, c in enumerate(cells) if c.endswith("S1")],   # S1
    [i for i, c in enumerate(cells) if c.endswith("S2")],   # S2
    [i for i, c in enumerate(cells) if c in ["thalE", "thalI"]]  # Thalamus
]


fig, axes = plt.subplots(nrows=results_map.shape[0], ncols=len(layer_splits), figsize=(14, 4*results_map.shape[0]), sharex=True)
#results_map.sort_values('input_duration', inplace=True, axis=0)
#results_map.sort_values('input_strength', inplace=True, axis=0)
#results_map.sort_values('background_input', inplace=True, axis=0)

if results_map.shape[0] == 1:
    axes = np.array([axes])  

idx = pd.IndexSlice
for row_axes, key in zip(axes, results_map.index):
    all_potentials_param = []
    names = []
    
    for i, cell in enumerate(cells):
        if i == (N_cells-2):
            potential_keys = [f'V_{cell}/{rpo}' for rpo in rpo_names_extended[:N_cells+1]]
        else: 
            potential_keys = [f'V_{cell}/{rpo}' for rpo in rpo_names[:N_cells]]

        if i in range(N_cells-2):
            potential_keys += [f'V_background/RPO_bI']

        sources = results.loc[:, idx[potential_keys, key]]
        summed = np.sum(sources, axis=1)
        all_potentials_param.append(summed)
        names.append(f"{cell} ({summed.iloc[-1]:.6f})")

    for ax, layer_slice in zip(row_axes, layer_splits):
        layer_series = [all_potentials_param[i] for i in layer_slice]
        layer_names  = [names[i] for i in layer_slice]

        for series, label in zip(layer_series, layer_names):
            ax.plot(series.index, series.values, label=label)

        ax.set_ylabel("[mV]")
        ax.set_xlabel("seconds")
        ax.legend(loc="best")
        
        #value = results_map.loc[key, par_name]*step_size # input amplitude
        value1 = results_map.loc[key, par_name1] 
        #value2 = results_map.loc[key, par_name2]
        ax.set_title(f"{par_name1}: {value1:.4f}")
        #ax.set_title(f"{par_name1}: {value1:.4f} / {par_name2}: {value2:.4f}")
        

#fig.suptitle(f"Potentials: {par_name1}/{par_name2}", fontsize=20, fontweight='bold')
fig.suptitle(f"Potentials: {par_name1}", fontsize=20, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

# %%

results_grid1 = results.loc[:, idx[:, 'model_0']]
results_grid2 = results.loc[:, idx[:, 'model_1']]
results_grid3 = results.loc[:, idx[:, 'model_2']]

# read in simulations from class
results_class0 = pd.read_csv(os.path.join(param_path, 'simulation_results_class_pyrates.csv'))  
results_class5 = pd.read_csv(os.path.join(param_path, 'simulation_results_class_pyrates_b5.csv'))  
results_class7 = pd.read_csv(os.path.join(param_path, 'simulation_results_class_pyrates_b7.csv'))  

# compare results from grid search with class simulations
results_grid1.columns = list(outputs.keys())
results_grid2.columns = list(outputs.keys())
results_grid3.columns = list(outputs.keys())

# %%
diff_class_grid0 = pd.DataFrame(
    results_class0.to_numpy() - results_grid1.to_numpy(),
    index=results_grid2.index,
    columns=results_grid2.columns
)

plt.plot(diff_class_grid0['V_E3b/RPO_E3b'])
plt.plot(diff_class_grid0['V_E3b/RPO_PV3b'])
plt.plot(diff_class_grid0['V_E3b/RPO_SST3b'])

# %%
diff_class_grid5 = pd.DataFrame(
    results_class5.to_numpy() - results_grid2.to_numpy(),
    index=results_grid1.index,
    columns=results_grid1.columns
)

plt.plot(diff_class_grid5['V_E3b/RPO_E3b'])
plt.plot(diff_class_grid5['V_E3b/RPO_PV3b'])
plt.plot(diff_class_grid5['V_E3b/RPO_SST3b'])

# %%
diff_class_grid7 = pd.DataFrame(
    results_class7.to_numpy() - results_grid3.to_numpy(),
    index=results_grid2.index,
    columns=results_grid2.columns
)

plt.plot(diff_class_grid7['V_E3b/RPO_E3b'])
plt.plot(diff_class_grid7['V_E3b/RPO_PV3b'])
plt.plot(diff_class_grid7['V_E3b/RPO_SST3b'])

# %%

diff_grid_grid = pd.DataFrame(
    results_grid1.to_numpy() - results_grid2.to_numpy(),
    index=results_grid1.index,
    columns=results_grid1.columns
)

plt.plot(diff_grid_grid['V_E3b/RPO_E3b'])
plt.plot(diff_grid_grid['V_E3b/RPO_PV3b'])
plt.plot(diff_grid_grid['V_E3b/RPO_SST3b'])

# %%
