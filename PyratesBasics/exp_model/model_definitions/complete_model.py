# %%
""" Complete PyRates model with the parameters from "parameters.py" and definitions compatible with PyCobi continuations. 
- background input in all non-thalamic populations
- time-varying external input in thalE
- connectivity parameters (g/bEI)
"""
# %%
import os 
os.chdir("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/""") 
from pyrates.frontend import OperatorTemplate, NodeTemplate, EdgeTemplate, CircuitTemplate
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
cells = ['E3b','PV3b','SST3b','VIP3b', # A3b
        'E1S1', 'PV1S1', 'SST1S1', 'VIPS1', 'E2S1', 'PV2S1', 'SST2S1', 'E3S1', 'PV3S1', 'SST3S1', 'E4S1', 'PV4S1', 'SST4S1', # S1
        'E1S2', 'PV1S2', 'SST1S2', 'VIPS2', 'E2S2', 'PV2S2', 'SST2S2', 'E3S2', 'PV3S2', 'SST3S2', 'E4S2', 'PV4S2', 'SST4S2', # S2
        'thalE', 'thalI'] # thalamus

N_cells = len(cells)
params = Parameter()

# %% Input definition
input_type = "step" # other options are "step", "baseline" (equals input strength 0) or "background"
input_onset = 1.0 # in sec
simulation_dur = 2.0
input_duration = 1.0  #, 1, 1.5] # np.arange(0, 1, 1) # in sec 
input_strength = 50.0 #[0, 50, 300, 500] #np.arange(0, 500, 100)
backgrndI_strengths = 0.0 #[0, 5, 10, 15, 20]
step_size=1e-3
sampling_step_size=1e-3
simulation_time = float(int(input_onset) + simulation_dur)

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
rpo_bI = OperatorTemplate(
    name='RPO_bI', path=None,
    equations=['d/dt * v_bI = i',
               'd/dt * i = H/tau * (bI) - 2 * i/tau - v_bI/tau^2'],
    variables={'v_bI': 'output',
               'i': 'variable',
               'bI': f'input({backgrndI_strengths})',  # external background input
               'tau': 0.003,
               'H': 1.0},
    description="excitatory rate-to-potential operator-background input")

# %% Operator template for the external input --> only for thalE!
rpo_Iext_thalE = [OperatorTemplate(
    name='RPO_Iext', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * (Iext) - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'Iext': 'input',
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
g_thal = 100.0
g = 100.0
connect_reverse_factor =  6448.0
bEI_thal = 0.5
bEI = 0.5

thal_connect = (0, 0, 0, 0)  # tEE, tEI, tIE, tII
W0, W_to_thal, W_from_thal, Wb, Wext = params.get_raw_connectivity(thal_connect)
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
                           equations= "m_out = (g*bEI/connect_reverse_factor)*m_outC",
                           variables={
                            "m_out": "output", 
                            #"g": g,
                            #"bEI":bEI,
                            "g": "input",
                            "bEI": "input",
                            "m_outC":"input",
                            "connect_reverse_factor": f"input({float(connect_reverse_factor)})"
                            }
                            )
connectivityI = OperatorTemplate(name="connectivityI", 
                           equations= "m_out = ((-g)*(1-bEI)/connect_reverse_factor)*m_outC",
                           variables={
                            "m_out": "output", 
                            #"g": -g,
                            #"bEI":bEI,
                            "g": "input",
                            "bEI": "input",
                            "connect_reverse_factor":f"input({float(connect_reverse_factor)})",
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
#circuit_to_yaml(model, "model.yaml")

# %%
# Run the simulation 
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

outputs['V_background/RPO_bI'] = 'BACKGROUND/RPO_bI/v_bI'

results = model.run(simulation_time=simulation_time,
                  step_size=step_size,
                  sampling_step_size=sampling_step_size,
                  outputs=outputs,
                  backend ="scipy",
                  vectorize=False,
                  clear=False,
                  float_precision="float64"
                  #decorator=njit
                  )

# %% Pandas Dataframe
all_potentials = []

for i, cell in enumerate(cells):
    # always include rpo_names[:N_cells]
    if i == (N_cells-2): # thalE
        potential_keys = [f'V_{cell}/{rpo}' for rpo in rpo_names_extended[:N_cells+1]]
    else: 
        potential_keys = [f'V_{cell}/{rpo}' for rpo in rpo_names[:N_cells]]

    # include rpo_name only if i in range(N_cells-2)
    if i in range(N_cells-2):
        potential_keys += [f'V_background/RPO_bI']

    # include rpo_names_extended[:N_cells+1] only if i == 13
    
    sources = results[potential_keys]
    all_potentials.append(np.sum(sources, axis=1))

all_potentials = np.array(all_potentials).T

potential_df = pd.DataFrame(all_potentials, columns=cells)

#potential_df.to_csv("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/simulations/a3bthal.csv") 

# %% Rates 
n_timepoints, n_cells = all_potentials.shape
m_out_all = np.zeros((n_timepoints, n_cells))  # NumPy array, not list!

for i, source_cell in enumerate(cells):
    m_out_all[:, i] =  m_max[i] / (
        1 + np.exp(r[i] * (v_thr[i] - all_potentials[:, i]))
    )

rates_df = pd.DataFrame(m_out_all, columns=cells)
# %% plot

# division into layers
layers = [potential_df.columns[:4],   # A3b 
          potential_df.columns[4:17],  # S1
          potential_df.columns[17:30],  # S2
          potential_df.columns[30:32] # thal
         ]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))  
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
#axes[-1].set_visible(False)
fig.suptitle("Potentials", fontsize=20, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

# %% plot - rates

# division into layers
layers = [rates_df.columns[:4],   # A3b 
          rates_df.columns[4:17],  # S1
          potential_df.columns[17:30],  # S2
          potential_df.columns[30:32] # thal
         ]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))  # 2x2 grid
axes = axes.flatten()

for ax, cols in zip(axes, layers):
    labels_with_final = [f"{col} ({rates_df[col].iloc[-1]:.6f})" for col in cols]
    
    rates_df[cols].plot(ax=ax,legend=False)
    ax.set_title(", ".join(cols))  # show which cols are in this subplot
    ax.legend(labels_with_final, loc="best")
    ax.set_ylabel("[Hz]")
    ax.set_xlabel("samples")
    
#axes[-1].set_visible(False)
fig.suptitle("Rates", fontsize=20, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show() 

# %%
