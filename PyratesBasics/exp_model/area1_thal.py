# %%
import os 
os.chdir("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/""") 
from pyrates.frontend import OperatorTemplate, NodeTemplate, CircuitTemplate
from copy import deepcopy
from parameters import Parameter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numba import njit
from yaml_saving import circuit_to_yaml
from pprint import pprint

#%%
# Parameters:
cells = ['E1', 'PV1', 'SST1', 'VIP', 'E2', 'PV2', 'SST2', 'E3', 'PV3', 'SST3', 'E4', 'PV4', 'SST4', 'thalE', 'thalI']

N_cells = len(cells)
params = Parameter()

# %%
sigmoid_params = params.get_sigmoid() #already in the correct order
sigm = np.vstack((sigmoid_params[4:17], sigmoid_params[-2:])) # area1 columns + last two for thalE and thalI

r = []
v_thr = []
m_max = []
for row in range(N_cells):
    r.append(sigm[row,0])      
    v_thr.append(sigm[row,1]) 
    m_max.append(sigm[row,2]) 

# %%
tau,_ = params.get_params()           
tau_a1_thal = np.hstack((tau[0, 4:17], tau[0, -2:]))

# default order of the taus: E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1
#tau_order = ['E1', 'E2', 'E3', 'E4', 'PV1', 'PV2', 'PV3', 'PV4', 'SST1', 'SST2', 'SST3', 'SST4', 'VIP']
#indices = [tau_order.index(c) for c in cells]

#tau_a1 = tau[indices]

# %%
bEI = 0.5
connect_reverse_factor =  6448 
g = 100.0 # (g)
gE = g * bEI /connect_reverse_factor
gI = g * (1 - bEI) /connect_reverse_factor
gEthal = 0
gIthal = 0
thal_connect = (0, 0, 0, 0)  # tEE, tEI, tIE, tII

W = params.get_connectivity(gE, gI, gEthal, gIthal, thal_connect, area='all') 

# selecting the region --> A1: FROM THE 5TH ELEMENT TO THE 17TH (INCLUDED) 
# in python the results include the start index but excludes the end index
rows = np.r_[4:17, 30, 31]
cols = np.r_[4:17, 30, 31]
W_A1_thal = W[np.ix_(rows, cols)]
"""
W[-2:,-4:-2] (within thalamus connectivity)
W[-2:,4:-4] (A1 to thalamus connectivity)
W[4:13,-4:-2] (thalamus to A1 connectivity)"""

# %%
# Operator template for the PRO

pro_names = ["PRO_" + cell for cell in cells]

rpo_names = ["RPO_" + cell for cell in cells]


pro = OperatorTemplate(
    name='PRO', path=None,
    equations=["m_out = m_max / (1 + exp(r*(V_thr - v)))"],
    variables={'m_out': 'output',
               'v': 'input',
               'r': 0.1,
               'V_thr': 35.0,
               'm_max': 70.0},
    description="sigmoidal potential-to-rate operator")

pros = [deepcopy(pro).update_template(name=pro_names[i],
                                      variables={'r': r[i],
                                                 'V_thr': v_thr[i],
                                                 'm_max': m_max[i]}) for i in range(N_cells)]

# Operator template for the RPO
rpo = OperatorTemplate(
    name='RPO', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * (m_in + bI + Iext) - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'm_in': 'input',
               'bI': 5.0,
               'Iext': 0.0,
               'tau': 0.01,
               'H': 1.0},
    description="excitatory rate-to-potential operator")


rpos_a1 = [deepcopy(rpo).update_template(name=rpo_names[i], variables={"tau": tau_a1_thal[i]}) for i in range(N_cells-2)]
rpos_thal = [deepcopy(rpo).update_template(name=rpo_names[i], variables={"tau": tau_a1_thal[i],"bI":0.0}) for i in [13,14]]
rpos = rpos_a1 + rpos_thal

"""
rpos = []
for i,cell in enumerate(cells):
    if cell not in {"thalE", "thalI"}:
        rpo_up = deepcopy(rpo).update_template(name=rpo_names[i], variables={"tau": tau_a1_thal[i]})
    else:
        rpo_up = deepcopy(rpo).update_template(name=rpo_names[i], variables={"tau": tau_a1_thal[i],"bI":0.0})
    rpos.append(rpo_up)
"""

# %%
# Node templates
nodes = [NodeTemplate(name=cells[i], path=None, operators= [pros[i]] + rpos) for i in range(N_cells)]

# Define edges
edges=[]
# i : target 
for i, cell_i in enumerate(cells):
    # j : source
    for j, cell_j in enumerate(cells):
        edges.append((f'{cell_j}/{pro_names[j]}/m_out', f'{cell_i}/{rpo_names[j]}/m_in', None, {'weight': W_A1_thal[i,j]}))

# Set up the Model Circuit 
area_1_thal = CircuitTemplate(
    name = 'area_1_thal',
    nodes = {name: node for name, node in zip(cells, nodes)},
    edges = edges,
    path = None)
"""
area_1_thal.update_var(node_vars={'thalE/RPO_thalE/bI': 0.0,
                    'thalI/RPO_thalI/bI': 0.0})"""

# %% save PyRates model template in a yaml file
#area_1.get_run_func(func_name='area_1', file_name='area_12', step_size=1e-4, auto=True, backend='python', solver='scipy',vectorize=False, float_precision='float64')
#circuit_to_yaml(area_1_thal, "area_1_thal.yaml")

# %% Input definition:

input_type = "step" # other options are "step", "baseline" (equals input strength 0) or "background"
input_onset = 1.001 # in sec
simulation_dur = 1 
input_duration = 0.5  #, 1, 1.5] # np.arange(0, 1, 1) # in sec 
input_strength = 50 #[0, 50, 300, 500] #np.arange(0, 500, 100)
backgrndI_strengths = 0 #[0, 5, 10, 15, 20]
step_size=1e-3
sampling_step_size=1e-3
simulation_time = int(input_onset) + simulation_dur

def create_Iext(simulation_time, step_size, input_onset, input_duration, input_strength, input_type):
    ''' Creates external input.
    Inputs in the function
    - simulation_time:
    - step_size:
    - input_onset:
    - input_duration:
    - input_strength: 
    - input_type:
    '''

    Iext = np.zeros(int(simulation_time/step_size))
    if input_type == "step":
        t  = int(input_duration/step_size)
        t0 = int(input_onset/step_size)
        Iext[t0:t0+t] = input_strength
    elif input_type == "background":
        # provide input for the entire simulation duration 
        Iext[:] = input_strength

    return Iext
Iext = create_Iext(simulation_time, step_size, input_onset, input_duration, input_strength, input_type)

# %%
# Run the simulation 
outputs = {}
for target_cell in cells:  
    for rpo_name in rpo_names[:N_cells]:
        outputs[f'V_{target_cell}/{rpo_name}'] = f'{target_cell}/{rpo_name}/v'
# %%
outputs_pro = {}
for idx, source_cell in enumerate(cells):
    pro_name = pro_names[idx]  
    outputs_pro[f'm_{source_cell}'] = f'{source_cell}/{pro_name}/m_out'

#area_1.get_variable_positions(outputs)
# %%
results = area_1_thal.run(simulation_time=simulation_time,
                  step_size=step_size,
                  sampling_step_size=sampling_step_size,
                  inputs={
                    'thalE/RPO_thalE/Iext': Iext
                    },
                  outputs=outputs,
                  backend ="scipy",
                  vectorize=True,
                  clear=False,
                  float_precision="float64",
                  decorator=njit
                  )

# %% Pandas Dataframe
all_potentials = []
for i in cells:
    sources = results[[f'V_{i}/{rpo_name}' for rpo_name in rpo_names[:N_cells]]]
    all_potentials.append(np.sum(sources, axis=1))
all_potentials = np.array(all_potentials).T # shape: samples x populations

potential_df = pd.DataFrame(all_potentials, columns=cells)

# potential_df.to_csv("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/simulations/area_1_thal.csv") 

# %% Rates 
n_timepoints, n_cells = all_potentials.shape
m_out_all = np.zeros((n_timepoints, n_cells))  # NumPy array, not list!

for i, source_cell in enumerate(cells):
    m_out_all[:, i] =  m_max[i] / (
        1 + np.exp(r[i] * (v_thr[i] - all_potentials[:, i] ))
    )

rates_df = pd.DataFrame(m_out_all, columns=cells)
#"m_out = m_max / (1 + exp(r*(V_thr - v)))"

  
# %% plot - potentials

# division into layers
folder = '/data/dpt_ticket102933_timeout20280820/simulations_results/pyrates/'
filename_potentials = "a1_thal_50_05.png"
save_path = os.path.join(folder, filename_potentials)
layers = [potential_df.columns[:4],   # first 4 
          potential_df.columns[4:7],  # next 3
          potential_df.columns[7:10], # next 3
          potential_df.columns[10:13],  # next 3
          potential_df.columns[13:]]  # thalamus

fig, axes = plt.subplots(2, 3, figsize=(14, 10))  # 2x2 grid
axes = axes.flatten()

for ax, cols in zip(axes, layers):
    # Create labels with final values
    labels_with_final = [f"{col} ({potential_df[col].iloc[-1]:.6f})" for col in cols]
    
    # Plot and override legend labels
    potential_df[cols].plot(ax=ax, legend=False)  # suppress default legend
    ax.legend(labels_with_final, loc="best")
    ax.set_ylabel("[mV]")
    ax.set_xlabel("samples")
    ax.set_title(", ".join(cols))
axes[-1].set_visible(False)
fig.suptitle("Potentials", fontsize=20, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(save_path, dpi=300, bbox_inches="tight")
plt.show()


# %% plot - rates

# division into layers
layers = [rates_df.columns[[0,4,7,10]],   # first 4 
          rates_df.columns[[1,5,8,11]],  # next 4
          rates_df.columns[[2,6,9,12]], # next 4
          rates_df.columns[[3]],
          rates_df.columns[[13,14]]
          ]  

fig, axes = plt.subplots(2, 3, figsize=(14, 10))  # 2x2 grid
axes = axes.flatten()

for ax, cols in zip(axes, layers):
    labels_with_final = [f"{col} ({rates_df[col].iloc[-1]:.6f})" for col in cols]
    
    rates_df[cols].plot(ax=ax,legend=False)
    ax.set_title(", ".join(cols))  # show which cols are in this subplot
    ax.legend(labels_with_final, loc="best")
    ax.set_ylabel("[Hz]")
    ax.set_xlabel("samples")
    
axes[-1].set_visible(False)
fig.suptitle("Rates", fontsize=20, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show() 
# %%
