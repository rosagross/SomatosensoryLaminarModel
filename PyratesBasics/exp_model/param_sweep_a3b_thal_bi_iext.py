# %%
from pyrates import grid_search
import os 
os.chdir("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/""") 
from pyrates.frontend import OperatorTemplate, NodeTemplate, CircuitTemplate
from copy import deepcopy
from parameters import Parameter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numba import njit
from pprint import pprint

# %% import of the .yaml file
# the path should end with the name of the file and then the name of the model
#/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/area_a3b_thal_bI_iext.yaml
model_path = "/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/area_a3b_thal_bI_iext.yaml/area_a3b_thal_bI_iext"
area_a3b_thal_bI_iext = CircuitTemplate.from_yaml(model_path)
"""attrs = dir(my_object)
print(attrs)
"""
# %% atrributes and parameters
cells_ext = list(area_a3b_thal_bI_iext.nodes.keys()) # getting the name of the populations directly from the model
cells = cells_ext[:-1]
N_cells = len(cells)

# names for the operators
pro_names = ["PRO_" + cell for cell in cells]
rpo_names = ["RPO_" + cell for cell in cells]
rpo_names_extended = rpo_names + ['RPO_Iext']

# sigmoid params for the rates
params = Parameter()
sigmoid_params = params.get_sigmoid() #already in the correct order
sigm = np.vstack((sigmoid_params[4:17], sigmoid_params[-2:])) # area1 columns + last two for thalE and thalI

r = []
v_thr = []
m_max = []
for row in range(N_cells):
    r.append(sigm[row,0])      
    v_thr.append(sigm[row,1]) 
    m_max.append(sigm[row,2]) 

# %% Setting the simulation:

input_type = "step" # other options are "step", "baseline" (equals input strength 0) or "background"
input_onset = 1.0 # in sec
simulation_dur = 2.0
input_duration = 0.5  #, 1, 1.5] # np.arange(0, 1, 1) # in sec 
input_strength = 500 #[0, 50, 300, 500] #np.arange(0, 500, 100)
backgrndI_strengths = 0 #[0, 5, 10, 15, 20]
step_size=1e-3
sampling_step_size=1e-3
simulation_time = float(int(input_onset) + simulation_dur)

# %%

# update of the model (might be different from the previous simulations)
#- Input duration:
area_a3b_thal_bI_iext.update_var(node_vars={'thalE/InputOp/A': input_strength}) 
area_a3b_thal_bI_iext.update_var(node_vars={'thalE/InputOp/dur': input_duration/step_size}) 
area_a3b_thal_bI_iext.update_var(node_vars={'thalE/InputOp/onset': input_onset/step_size})
area_a3b_thal_bI_iext.update_var(node_vars={'BACKGROUND/RPO_bI/bI': backgrndI_strengths})

# %%
"""
# INPUT:
#param_grid = {'input_duration': np.array([0, 0.5, 1, 1.5]) / step_size} # values to "sweep" (defined as dictionary)  / the name is not related to the PyRates definition
param_grid = {'input_strength': np.array([0, 50, 300, 500]),
            'input_duration': np.array([0, 0.5, 1, 1.5]) / step_size}
param_map = {'input_strength': {'vars': ['thalE/InputOp/A'], 'nodes': ['thalE']},
            'input_duration': {'vars': ['thalE/InputOp/dur'], 'nodes': ['thalE']}
            } # maps the parameters to the model 

# BACKGROUND INPUT:
"""
param_grid = {'background_input': np.array([0, 5, 10, 15, 20, 50, 100, 200])} # values to "sweep" (defined as dictionary)  / the name is not related to the PyRates definition
param_map = {
    'background_input': {
        'vars': ['BACKGROUND/RPO_bI/bI'],
        'nodes': ['BACKGROUND']
    }
}

# %% 
rpo_names_extended = rpo_names + ['RPO_Iext']
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

results, results_map = grid_search(circuit_template=area_a3b_thal_bI_iext,
                                   param_grid=param_grid,
                                   param_map=param_map,
                                   simulation_time=simulation_time,
                                   step_size=step_size,
                                   sampling_step_size=sampling_step_size,
                                   #inputs={'PC/RPO_e_in/u': noise},
                                   inputs = None, # the input is inside the model
                                   outputs=outputs,
                                   backend="scipy",
                                   permute_grid=True
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

layer_splits = [slice(0, 4), slice(4, None)]

fig, axes = plt.subplots(nrows=results_map.shape[0], ncols=2, figsize=(14, 4*results_map.shape[0]), sharex=True)
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
        layer_series = all_potentials_param[layer_slice]
        layer_names = names[layer_slice]

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