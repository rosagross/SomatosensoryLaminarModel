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

#%%
# Parameters:
cells = ['E3b','PV3b','SST3b','VIP3b','ThalE','ThalI'] # taken from the weights matrix

N_cells = len(cells)
params = Parameter()

# %% sigmoid parameters
sigm_a3b = params.get_sigmoid("A3b") #already in the correct order
sigm_s2 = params.get_sigmoid("S2") # thalE and thalI are the last two rows of S2
sigm_s2 = sigm_s2[-2:,:]
sigm = np.vstack((sigm_a3b, sigm_s2))

r = []
v_thr = []
m_max = []
for row in range(N_cells):
    r.append(sigm[row,0])      
    v_thr.append(sigm[row,1]) 
    m_max.append(sigm[row,2]) 

# %% dendritic time constants
tau,_ = params.get_params("A3b")           
tau = tau[0] # + QUARTULTIMO E TERZULTIMO DI S2
tau_s2,_ = params.get_params("S2") 
tau_s2 = tau_s2[0]
tau_thal = tau_s2[-4:-2]
tau_a3b = np.hstack((tau,tau_thal))

# %% weights
bEI = 0.5
connect_reverse_factor =  6448 
g = 100.0 # (g)
gE = g * bEI /connect_reverse_factor
gI = g * (1 - bEI) /connect_reverse_factor
gEthal = 1  # NOT 0 BECAUSE IN THAT CASE THERE IS NO COUPLING
gIthal = 1
thal_connect = (0, 0, 0, 0)  # tEE, tEI, tIE, tII

W = params.get_connectivity(gE, gI, gEthal, gIthal, thal_connect, include_Iext=True) 

# selecting the region: 
# in python the results include the start index but excludes the end index

W_A3b = np.block([
    [W[0:4, 0:4],       W[0:4, -4:-2]],
    [W[-2:, 0:4],              W[-2:, -4:-2]]
])
# %%
# Operator template for the PRO

pro_names = ["PRO_"+ cell for cell in cells]
rpo_names=["RPO_"+cell for cell in cells]


pro = OperatorTemplate(
    name='PRO', path=None,
    equations=["m_out = 2.*m_max / (1 + exp(r*(V_thr - (v*1000))))"],
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
               'd/dt * i = H/tau * (m_in + bI) - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'm_in': 'input',
               'bI': 0.0,
               'tau': 0.01,
               'H': 1.0},
    description="excitatory rate-to-potential operator")


rpos = [deepcopy(rpo).update_template(name=rpo_names[i], variables={"tau": tau_a3b[i]}) for i in range(N_cells)]

# %%
# Node templates
nodes = [NodeTemplate(name=cells[i], path=None, operators= [pros[i]] + rpos) for i in range(N_cells)]

# Define edges
edges=[]
# i : target 
for i, cell_i in enumerate(cells):
    # j : source
    for j, cell_j in enumerate(cells):
        edges.append((f'{cell_j}/{pro_names[j]}/m_out', f'{cell_i}/{rpo_names[j]}/m_in', None, {'weight': W_A3b[i,j]}))

# Set up the Model Circuit 
a3b = CircuitTemplate(
    name = 'a3b',
    nodes = {name: node for name, node in zip(cells, nodes)},
    edges = edges,
    path = None)

# %%
# Run the simulation 
outputs = {}
for target_cell in cells:  
    for rpo_name in rpo_names[:N_cells]:
        outputs[f'V_{target_cell}/{rpo_name}'] = f'{target_cell}/{rpo_name}/v'

results = a3b.run(simulation_time=2.0,
                  step_size=1e-4,
                  sampling_step_size=1e-4,
                  outputs=outputs,
                  backend='default',
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

potential_df.to_csv("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/simulations/a3bthal.csv") 

# %% plot

#potential_df.plot()
#plt.show()

# division into layers
layers = [potential_df.columns[:4],   # first 4 
          potential_df.columns[4:]  # next 4
         ]  # in case of extras (not needed here)

fig, axes = plt.subplots(1, 2, figsize=(14, 10))  
axes = axes.flatten()

for ax, cols in zip(axes, layers):
    potential_df[cols].plot(ax=ax)
    ax.set_title(", ".join(cols))  # show which cols are in this subplot
    ax.legend(loc="best")

plt.tight_layout()
plt.show()

# %%
