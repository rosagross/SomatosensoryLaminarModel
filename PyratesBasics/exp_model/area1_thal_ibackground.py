# %%
import os 
os.chdir("/data/hu_grossmannr/Desktop/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel/PyratesBasics/exp_model") 
from pyrates.frontend import OperatorTemplate, NodeTemplate, CircuitTemplate
from copy import deepcopy
from parameters import Parameter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#from numba import njit
from yaml_saving import circuit_to_yaml
from pprint import pprint

#%%
# Parameters:
cells = ['E1', 'PV1', 'SST1', 'VIP', 'E2', 'PV2', 'SST2', 'E3', 'PV3', 'SST3', 'E4', 'PV4', 'SST4', 'thalE', 'thalI']

N_cells = len(cells)
params = Parameter()

# %% Input definition

input_type = "step" # other options are "step", "baseline" (equals input strength 0) or "background"
input_onset = 1.001 # in sec
simulation_dur = 1 
input_duration = 0.5  #, 1, 1.5] # np.arange(0, 1, 1) # in sec 
input_strength = 0 #[0, 50, 300, 500] #np.arange(0, 500, 100)
backgrndI_strengths = 50 #[0, 5, 10, 15, 20]
step_size=1e-3
sampling_step_size=1e-3
simulation_time = int(input_onset) + simulation_dur

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

W = params.get_connectivity(gE, gI, gEthal, gIthal, thal_connect) 

# selecting the region --> A1: FROM THE 5TH ELEMENT TO THE 17TH (INCLUDED) 
# in python the results include the start index but excludes the end index
rows = np.r_[4:17, -2, -1]
cols = np.r_[4:17, -4, -3]
W_A1_thal = W[np.ix_(rows, cols)]

"""
W[-2:,-4:-2] (within thalamus connectivity)
W[-2:,4:-4] (A1 to thalamus connectivity)
W[4:13,-4:-2] (thalamus to A1 connectivity)"""

# %% operator names
pro_names = ["PRO_" + cell for cell in cells]
rpo_names = ["RPO_" + cell for cell in cells]

# %% Operator templates for the PRO
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

# %% Operator templates for the RPO
rpo = OperatorTemplate(
    name='RPO', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * (m_in) - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'm_in': 'input',
               'tau': 0.01,
               'H': 1.0},
    description="excitatory rate-to-potential operator")

rpos = [deepcopy(rpo).update_template(name=rpo_names[i], variables={"tau": tau_a1_thal[i]}) for i in range(N_cells)]

rpo_bI = OperatorTemplate(
    name='RPO_bI', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * (bI) - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'bI': f'input({bckgrndI_strengths})',  # external background input
               'tau': 0.003,
               'H': 1.0},
    description="excitatory rate-to-potential operator")

create_bI = [OperatorTemplate(
    name="backInputOp", path=None,
    equations=[
        "bI = background_input"
    ],
    variables={
        "bI": "output",
        "background_input": backgrndI_strengths
    },
    description="background input"
)]

# %% Operator template for the external input --> only for thalE!
rpo_Iext_thalE = [OperatorTemplate(
    name='RPO_Iext', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * (Iext) - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'Iext': 'input',
               'tau': tau_a1_thal[13],
               'H': 1.0}
    ) 
]
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
    },
    description="External step input"
)]
# %%
# Node templates
nodes = [
    NodeTemplate(
        name=cells[i],
        path=None,
        operators=(
            [pros[i]] + rpos +
            ([rpo_bI] if i in range(N_cells-2) else []) +
            (rpo_Iext_thalE + create_I_ext if i == 13 else [])
        )
    )
    for i in range(N_cells)
]

# %% Edges
edges=[]
# i : target 
for i, cell_i in enumerate(cells):
    # j : source
    for j, cell_j in enumerate(cells):
        edges.append((f'{cell_j}/{pro_names[j]}/m_out', f'{cell_i}/{rpo_names[j]}/m_in', None, {'weight': W_A1_thal[i,j]}))

# Set up the Model Circuit 
area_1_thal_bI = CircuitTemplate(
    name = 'area_1_thal_bI',
    nodes = {name: node for name, node in zip(cells, nodes)},
    edges = edges,
    path = None)
"""
area_1_thal.update_var(node_vars={'thalE/RPO_thalE/bI': 0.0,
                    'thalI/RPO_thalI/bI': 0.0})"""

"""
for node_name, node in area_1_thal_bI.nodes.items():
    print(f"\nNode: {node_name}")
    for op in node.operators:
        print(f"  Operator: {op.name}")
        if any("Iext" in eq for eq in op.equations):
            print("    ⚠ contains Iext equation!")
        elif "Iext" in op.variables:
            print("    ⚠ has Iext variable!")
"""

# %% save PyRates model template in a yaml file
#area_1.get_run_func(func_name='area_1', file_name='area_12', step_size=1e-4, auto=True, backend='python', solver='scipy',vectorize=False, float_precision='float64')

#circuit_to_yaml(area_1_thal_bI, "area_1_thal_bI.yaml")

rpo_names_extended = rpo_names + ['RPO_Iext']
# %%
# Run the simulation 
outputs = {}
b_inputs = []

for i, target_cell in enumerate(cells):
    if i == 13:
        for rpo_name_ext in rpo_names_extended[:N_cells+1]:
            outputs[f'V_{target_cell}/{rpo_name_ext}'] = f'{target_cell}/{rpo_name_ext}/v'
    else:
        for rpo_name in rpo_names[:N_cells]:
            outputs[f'V_{target_cell}/{rpo_name}'] = f'{target_cell}/{rpo_name}/v'

    # Include rpo_names_bI only if i in range(N_cells-2)
    if i in range(N_cells-2):
        print(f'{target_cell}/RPO_bI/v')
        outputs[f'V_{target_cell}/RPO_bI'] = f'{target_cell}/RPO_bI/v'

    # Include rpo_names_extended only if i == 13
    

# %%
def create_Ib(simulation_time, step_size, input_strength):
    ''' Creates external background input.
    Inputs in the function
    - simulation_time:
    - step_size:
    - input_onset:
    - input_duration:
    - input_strength: 
    '''

    Iext = np.zeros(int(simulation_time/step_size))
    
    # provide input for the entire simulation duration 
    Iext[:] = input_strength

    return Iext

bI_array = create_Ib(simulation_time, step_size, backgrndI_strengths)
 

# TODO: define input for bI 
results = area_1_thal_bI.run(simulation_time=simulation_time,
                  step_size=step_size,
                  sampling_step_size=sampling_step_size,
                  #inputs={
                  #  'thalE/RPO_thalE/Iext': Iext
                  #  },
                  inputs={
                    f'{cell}/RPO_bI/Ib': bI_array for cell in cells[:N_cells-2]
                    },
                  outputs=outputs,
                  backend ="scipy",
                  vectorize=True,
                  clear=False,
                  float_precision="float64",
                  #decorator=njit
                  )

# %% Pandas Dataframe
all_potentials = []

for i, cell in enumerate(cells):
    # always include rpo_names[:N_cells]
    if i == 13:
        potential_keys = [f'V_{cell}/{rpo}' for rpo in rpo_names_extended[:N_cells+1]]
    else: 
        potential_keys = [f'V_{cell}/{rpo}' for rpo in rpo_names[:N_cells]]

    # include rpo_name only if i in range(N_cells-2)
    if i in range(N_cells-2):
        potential_keys += [f'V_{cell}/RPO_bI']

    # include rpo_names_extended[:N_cells+1] only if i == 13
    

    sources = results[potential_keys]
    all_potentials.append(np.sum(sources, axis=1))

all_potentials = np.array(all_potentials).T  # shape: samples x populations


potential_df = pd.DataFrame(all_potentials, columns=cells)

# potential_df.to_csv("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/simulations/area_1_thal.csv") 


# %% Rates 
n_timepoints, n_cells = all_potentials.shape
m_out_all = np.zeros((n_timepoints, n_cells))  # NumPy array, not list!

for i, source_cell in enumerate(cells):
    m_out_all[:, i] =  m_max[i] / (
        1 + np.exp(r[i] * (v_thr[i] - all_potentials[:, i]))
    )

rates_df = pd.DataFrame(m_out_all, columns=cells)
#"m_out = m_max / (1 + exp(r*(V_thr - v)))"

  
# %% plot - potentials

# division into layers

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
