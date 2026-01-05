"""
File: grid_search.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-01-05
Description: Grid search based on circuit template implementation. 
"""


# %%
import time
import os 
WDDIR = os.getenv("WDDIR")
SIMDIR = os.getenv("SIMDIR")
os.chdir(os.path.join(WDDIR,"PyratesBasics","exp_model")) 
import sys
from pyrates.frontend import OperatorTemplate, NodeTemplate, EdgeTemplate, CircuitTemplate
from pyrates import grid_search
from copy import deepcopy
from complete_model_class import SomatoModelPyrates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numba import njit
from pprint import pprint
import json
param_path = os.path.join(WDDIR,"Simulations")

if param_path not in sys.path:
    sys.path.append(param_path)
from parameters import Parameter
import json

# %% 
# Set-Up the model 
# Parameters:
cells = ['E3b','PV3b','SST3b','VIP3b', # A3b
        'E1S1', 'PV1S1', 'SST1S1', 'VIPS1', 'E2S1', 'PV2S1', 'SST2S1', 'E3S1', 'PV3S1', 'SST3S1', 'E4S1', 'PV4S1', 'SST4S1', # S1
        'E1S2', 'PV1S2', 'SST1S2', 'VIPS2', 'E2S2', 'PV2S2', 'SST2S2', 'E3S2', 'PV3S2', 'SST3S2', 'E4S2', 'PV4S2', 'SST4S2', # S2
        'thalE', 'thalI'] # thalamus

N_cells = len(cells)
params = Parameter()

# %% Input definition
input_type = "step" # other options are "step", "baseline" (equals input strength 0) or "background"
input_onset = 1.0 # in sec
simulation_dur = 2
input_duration = 0  #, 1, 1.5] # np.arange(0, 1, 1) # in sec 
input_strength = 0 #[0, 50, 300, 500] #np.arange(0, 500, 100)
backgrndI_strengths = 0 #[0, 5, 10, 15, 20]
step_size=1e-3
sampling_step_size=1e-3
simulation_time = 2.0 #float(int(input_onset) + simulation_dur)

# %% sigmoid parameters and dendritic time constants
sigm = params.get_sigmoid()  #already in the correct order

r = []
v_thr = []
m_max = []
for row in range(N_cells):
    r.append(sigm[row,0])      
    v_thr.append(sigm[row,1]) 
    m_max.append(sigm[row,2]) 

tau,_ = params.get_params()           
tau = tau[0,:]

# %%
# Operator template for the PRO
pro_names = ["PRO_"+ cell for cell in cells]
rpo_names=["RPO_"+cell for cell in cells]
# no background input:
pro = OperatorTemplate(
    name='PRO', path=None,
    equations=["m_outC = m_max / (1 + exp(r*(V_thr - v)))"],
    variables={'m_outC': 'output',
               'v': 'input',
               'r': 0.1,
               'V_thr': 35.0,
               'm_max': 70.0},
    description="sigmoidal potential-to-rate operator")
# background input:
pro_bI = OperatorTemplate(
    name='PRO_bI', path=None,
    equations=["m_outC = m_max / (1 + exp(r*(V_thr - (v+v_bIn))))"],
    variables={'m_outC': 'output',
               'v': 'input',
               'v_bIn': 'input',
               'r': 0.1,
               'V_thr': 35.0,
               'm_max': 70.0},
    description="sigmoidal potential-to-rate operator")

pros = [
    (   # background input: all populations except thalamus
        deepcopy(pro_bI).update_template(name=pro_names[i], variables={'r': r[i],
                                                 'V_thr': v_thr[i],
                                                 'm_max': m_max[i]})
        if i < N_cells - 2
        # no background input: thalamus
        else deepcopy(pro).update_template(name=pro_names[i], variables={'r': r[i],
                                                 'V_thr': v_thr[i],
                                                 'm_max': m_max[i]})
    )
    for i in range(N_cells)
]

# %% Operator template for the RPO
rpo = OperatorTemplate(
    name='RPO', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * (m_in ) - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'm_in': 'input',
               'tau': 0.01,
               'H': 1.0},
    description="excitatory rate-to-potential operator")

rpos = [deepcopy(rpo).update_template(name=rpo_names[i], variables={"tau": tau[i]}) for i in range(N_cells)]

# %% Operator template for the background input --> only for the "BACKGROUND" population!
bI_cellcount = 100 

rpo_bI = OperatorTemplate(
    name='RPO_bI', path=None,
    equations=['d/dt * v_bI = i',
               'd/dt * i = H/tau * bI_cellcount * (bI) - 2 * i/tau - v_bI/tau^2'],
    variables={'v_bI': 'output',
               'i': 'variable',
               'bI': f'input({backgrndI_strengths})',  # external background input
               'bI_cellcount': bI_cellcount,
               'tau': 0.003,
               'H': 1.0},
    description="excitatory rate-to-potential operator-background input")

# %% Operator template for the external input --> only for thalE!
extI_cellcount = 1000

rpo_Iext_thalE = [OperatorTemplate(
    name='RPO_Iext', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * Iext_cellcounts * (Iext) - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'Iext': 'input',
               'Iext_cellcounts': extI_cellcount,
               'tau': tau[N_cells-2],
               'H': 1.0}
    ) 
]
# definition of the external input
create_I_ext = [OperatorTemplate(
    name="InputOp", path=None,
    equations=[
        # step input: A during [onset, onset+dur), else 0
        "Iext = (A/2)*(1+sign(t-onset)) - (A/2)*(1+sign(t-(onset+dur)))"
    ],
    variables={
        "Iext": "output",
        "t": "variable",
        "A": input_strength,
        "onset": input_onset/step_size,
        "dur": input_duration/step_size
        # TO UNCOMMENT FOR PYCOBI:
        #"A": float(input_strength),
        #"onset": float(input_onset),
        #"dur": float(input_duration)
    },
    description="External step input"
)]

# %% Weights 
# 'raw' weights:
g_thal = 2
g = 10.0
bEI_thal = 0.5
bEI = 0.5

thal_connect = (0, 0, 0, 0)  # tEE, tEI, tIE, tII

thal_cellcount = 500
W0, W_to_thal, W_from_thal, Wb, Wext = params.get_raw_connectivity(thal_connect, extI_cellcount, bI_cellcount, thal_cellcount)
#W0, W_to_thal, W_from_thal, Wb, Wext = params.get_raw_connectivity(thal_connect)
W0 = np.append(W0, W_to_thal, axis=0)
W0 = np.append(W0, W_from_thal.T, axis=1)
W = np.concatenate((W0, Wb, Wext), axis=1)
#rows_A3b_thal = np.r_[:4, 30, 31]   
#cols_A3b_thal = np.r_[:4, 30, 31]
#W_A3b = W[rows_A3b_thal[:, None], cols_A3b_thal]

# %% Nodes for the parameters
g_definition = OperatorTemplate(
    name="g_definition",
    equations="gC = g_input",
    variables={"gC": "output", "g_input": f"input({float(g)})"}
)
bEI_definition = OperatorTemplate(
    name="bEI_definition",
    equations="bEIC = bEI_input",
    variables={"bEIC": "output", "bEI_input": f"input({float(bEI)})"}
)
g_thal_definition = OperatorTemplate(
    name="g_thal_definition",
    equations="g_thalC = g_thal_input",
    variables={"g_thalC": "output", "g_thal_input": f"input({float(g_thal)})"}
)
bEI_thal_definition = OperatorTemplate(
    name="bEI_thal_definition",
    equations="bEI_thalC = bEI_thal_input",
    variables={"bEI_thalC": "output", "bEI_thal_input": f"input({float(bEI_thal)})"}
)
# %% Connectivity operators
connectivity_names = ["connectivity_"+ cell for cell in cells]
connectivityE = OperatorTemplate(name="connectivityE", 
                           equations= "m_out = (g*bEI)*m_outC",
                           variables={
                            "m_out": "output", 
                            #"g": g,
                            #"bEI":bEI,
                            "g": "input",
                            "bEI": "input",
                            "m_outC":"input"
                            #"connect_reverse_factor": f"input({float(connect_reverse_factor)})"
                            }
                            )
connectivityI = OperatorTemplate(name="connectivityI", 
                           equations= "m_out = ((-g)*(1-bEI))*m_outC",
                           variables={
                            "m_out": "output", 
                            #"g": -g,
                            #"bEI":bEI,
                            "g": "input",
                            "bEI": "input",
                            #"connect_reverse_factor":f"input({float(connect_reverse_factor)})",
                            "m_outC":"input"
                            }
                            )
connectivityE_thal = OperatorTemplate(name="connectivityE_thal", 
                           equations= "m_out = (g_thal*bEI_thal)*m_outC",
                           variables={
                            "m_out": "output", 
                            #"g_thal": g_thal,
                            #"bEI_thal":bEI_thal,
                            "g_thal": "input",
                            "bEI_thal": "input",
                            #"connect_reverse_factor_thal":f"input({float(connect_reverse_factor)})",
                            "m_outC":"input"
                            }
                            )
connectivityI_thal = OperatorTemplate(name="connectivityI_thal", 
                           equations= "m_out = ((-g_thal)*(1-bEI_thal))*m_outC",
                           variables={
                            "m_out": "output", 
                            #"g_thal": -g_thal,
                            #"bEI_thal":bEI_thal,
                            "g_thal": "input",
                            "bEI_thal": "input",
                            #"connect_reverse_factor_thal":f"input({float(connect_reverse_factor)})",
                            "m_outC":"input"
                            }
                            )

# the connectivity operators depend on the populations (if they are excitatory/inhibitory) --> parameters.py
# inhibitory
idx_I_A3b = np.array([1,2,3])
idx_I_S = [1,2,3,5,6,8,9,11,12]
idx_I_S1 = np.array(idx_I_S)+4 
idx_I_S2 = np.array(idx_I_S)+17
idx_I = np.concatenate((idx_I_A3b,idx_I_S1,idx_I_S2))
# excitatory
idx_E_A3b = np.array([0])
idx_E_S = [0,4,7,10]
idx_E_S1 = np.array(idx_E_S)+4
idx_E_S2 = np.array(idx_E_S)+17
idx_E = np.concatenate((idx_E_A3b,idx_E_S1,idx_E_S2))

connectivity = []
for i in range(N_cells):
    if i in idx_E:  # excitatory populations
        connectivity.append(deepcopy(connectivityE).update_template(name=connectivity_names[i]))
    elif i in idx_I:  # inhibitory populations
        connectivity.append(deepcopy(connectivityI).update_template(name=connectivity_names[i]))
    elif i == N_cells - 2: # thalE
        connectivity.append(deepcopy(connectivityE_thal).update_template(name=connectivity_names[i]))
    else: # thalI
        connectivity.append(deepcopy(connectivityI_thal).update_template(name=connectivity_names[i]))

# %%
# Node templates
cells_ext = cells + ['BACKGROUND'] + ['G'] + ['G_THAL'] + ['BEI'] + ['BEI_THAL'] 
nodes = [
    NodeTemplate(
        name=cells_ext[i],
        path=None,
        operators=(
            ([pros[i]] + [connectivity[i]] + rpos + 
             (rpo_Iext_thalE + create_I_ext if i == N_cells - 2 else [])
            ) if i < N_cells 
            #else [rpo_bI] if i == N_cells
            else [g_definition] if i == N_cells + 1 # G
            else [g_thal_definition] if i == N_cells + 2 # G_THAL
            else [bEI_definition] if i == N_cells + 3 # BEI
            else [bEI_thal_definition] if i == N_cells + 4 # BEI_THAL
            else [rpo_bI] # operator for the background input
        )
    )
    for i in range(len(cells_ext))
]

# %% 
# Edges
edges=[]
# i : target 
for i, cell_i in enumerate(cells):
    if cell_i not in ['thalE', 'thalI']:
        edges.append(('BACKGROUND/RPO_bI/v_bI', f'{cell_i}/{pro_names[i]}/v_bIn', None, {'weight': 1.0})) # background input
        edges.append(('G/g_definition/gC', f'{cell_i}/{connectivity_names[i]}/g', None, {'weight': 1.0})) # G
        edges.append(('BEI/bEI_definition/bEIC', f'{cell_i}/{connectivity_names[i]}/bEI', None, {'weight': 1.0})) # BEI
    else:
        edges.append(('G_THAL/g_thal_definition/g_thalC', f'{cell_i}/{connectivity_names[i]}/g_thal', None, {'weight': 1.0})) # G_THAL
        edges.append(('BEI_THAL/bEI_thal_definition/bEI_thalC', f'{cell_i}/{connectivity_names[i]}/bEI_thal', None, {'weight': 1.0})) # BEI_THAL
    for j, cell_j in enumerate(cells):
        edges.append((f'{cell_j}/{connectivity_names[j]}/m_out', f'{cell_i}/{rpo_names[j]}/m_in', None, {'weight': W[i,j]}))
    
# %%
# Set up the Model Circuit 
model = CircuitTemplate(
    name = 'model',
    nodes = {name: node for name, node in zip(cells_ext, nodes)},
    edges = edges,
    path = None)

# %%
# read pyrates model from yaml circuit file
start = time.time()
pmdl = CircuitTemplate.from_yaml("expanded_model.yaml/model")
end = time.time()
print("loading time:", start-end)

# %% Simulation parameters
with open(os.path.join(WDDIR, 'Simulations', 'simulation_parameter.json'), 'r') as json_file:
    params = json.load(json_file)
    
# %%
param_grid = {'background_input': np.array([0, 5])} # values to "sweep" (defined as dictionary)  / the name is not related to the PyRates definition
param_map = {
    'background_input': {
        'vars': ['BACKGROUND/RPO_bI/bI'],
        'nodes': ['BACKGROUND']
    }
}

# %% 
N_cells = len(cells)
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

# %%
results, results_map = grid_search(circuit_template=model,
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
#results_grid3 = results.loc[:, idx[:, 'model_2']]

# read in simulations from class
results_class0 = pd.read_csv(os.path.join(param_path, 'simulation_results_class_pyrates.csv')).iloc[:2000]
results_class5 = pd.read_csv(os.path.join(param_path, 'simulation_results_class_pyrates_b5.csv')).iloc[:2000] 
results_class7 = pd.read_csv(os.path.join(param_path, 'simulation_results_class_pyrates_b7.csv')).iloc[:2000]

# compare results from grid search with class simulations
results_grid1.columns = list(outputs.keys())
results_grid2.columns = list(outputs.keys())
#results_grid3.columns = list(outputs.keys())

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
