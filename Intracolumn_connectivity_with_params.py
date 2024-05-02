# %% 
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy  
from pyrates.frontend import OperatorTemplate, NodeTemplate, CircuitTemplate
from parameters import Parameter
import pandas as pd

# Load parameters
cortex_type = 'visual'
params = Parameter(cortex_type)
tau, N_cells = params.get_params()
sigm = params.get_sigmoid()
W = params.get_connectivity(1, include_Iext=False)

# Define Potential-to-Rate operators
pro_names = ['PRO_E1', 'PRO_E2', 'PRO_E3', 'PRO_E4', 'PRO_P1', 'PRO_P2', 
             'PRO_P3', 'PRO_P4', 'PRO_S1', 'PRO_S2', 'PRO_S3', 'PRO_S4',
             'PRO_V1', 'PRO_V2', 'PRO_V3', 'PRO_V4']
             
pro = OperatorTemplate(
    name = 'pro', 
    path = None, 
    equations = ["m_out = m_max/(1 + exp(r*(vth-v)))"], 
    variables = {'m_out': 'output',
                    'v': 'input',
                    'r': 0.14218422,  
                    'vth': 40.03107351,
                    'm_max': 166.82960408}, 
    description = "sigmoidal potential-to-rate operator")
pros = [deepcopy(pro).update_template(name=pro_names[i],
                	                    variables={"r": sigm[i,0], 'vth': sigm[i,1],
                                                   'm_max': sigm[i,2]}) for i in range(N_cells)]


# Define Rate-to-Potential operators
rpo_names = ['RPO_E1', 'RPO_E2', 'RPO_E3', 'RPO_E4', 'RPO_P1', 'RPO_P2',
             'RPO_P3', 'RPO_P4', 'RPO_S1', 'RPO_S2', 'RPO_S3', 'RPO_S4',
             'RPO_V1', 'RPO_V2', 'RPO_V3', 'RPO_V4']
rpo = OperatorTemplate(
        name = 'rpo', 
        path = None, 
        equations = ['d/dt * v = u',
                    'd/dt * u = H/tau * r - 2 * u/tau - v/(tau**2)'],  
        variables = {'v': 'output',
                    'u': 'variable',
                    'r': 'input',
                    'tau': 1,
                    'H': 1},
        description = "rate-to-potential operator")

# ToDo: is tau maybe from the target cell? 
rpos = [deepcopy(rpo).update_template(name=rpo_names[i], variables={"tau": 0.01}) for i in range(N_cells)]

# Define network nodes
cells = ['E1', 'E2', 'E3', 'E4', 'P1', 'P2', 'P3', 'P4', 'S1', 'S2', 'S3', 'S4', 'V1', 'V2', 'V3', 'V4']
# every population has the same operators (since they are fully connected)
node = NodeTemplate(name="node", path=None, operators=list(pros) + list(rpos)) 

# Define edges
edges=[]
# i : target 
for i, cell_i in enumerate(cells):
    # j : source
    for j, cell_j in enumerate(cells):
        edges.append((f'{cell_j}/{pros[j].name}/m_out', f'{cell_i}/{rpos[j].name}/r', None, {'weight': W[i,j]})) # to x from 

# Define circuit
cir = CircuitTemplate( 
    name="cir", 
    nodes={name: node for name in cells},
    edges=edges,
    path = None 
)

 
# %%

outputs = {}

#Membrane potential
for target_cell in cells:  
    for rpo_name in rpo_names:
        outputs[f'V_{target_cell}/{rpo_name}'] = f'{target_cell}/{rpo_name}/v'
        '''
#Firing Rate
for source_cell in cells: 
    for pro_name in pro_names:
        outputs[f'Rate_{source_cell}/{pro_name}'] = f'{source_cell}/{pro_name}/m_out'
'''

simulation_time = 2.0
step_size = 1e-3
sampling_step_size = 1e-3

results = cir.run(simulation_time = simulation_time,
                  step_size = step_size,
                  sampling_step_size=sampling_step_size,
                  outputs = outputs,
                  backend ='default',
                  solver ='scipy',
                  vectorize=False)

# %% 
test_potentials1 = []
sources = results[[f'V_E2/{rpo_name}' for rpo_name in rpo_names]]
test_potentials1.append(np.sum(sources, axis=1))

# %%
time_list = np.arange(0, simulation_time, sampling_step_size)

all_potentials = []
for i in cells:
    sources = results[[f'V_{i}/{rpo_name}' for rpo_name in rpo_names]]
    all_potentials.append(np.sum(sources, axis=1))

all_potentials = np.array(all_potentials).T
potential_df = pd.DataFrame(all_potentials, columns=cells)
potential_df.to_csv('output/pyrates_potential_G1.csv', index=False)

# take sigmoid values for population
e1_sigm = sigm[0]
#  r, v_thr, m_max
#fr_e1 = e1_sigm[2]/(1 + np.exp(e1_sigm[0]*(e1_sigm[1]-all_potentials[:, 0])))
#plt.plot(time_list, fr_e1)
#plt.show()

# %% 

cell_potential = np.zeros((N_cells, len(time_list)))
firing_rate = np.zeros((N_cells, len(time_list)))
for i,target in enumerate(cells):
    #Membrane Potential
    
    for rpo_name in rpo_names:
        plt.plot(time_list, results[f'V_{target}/{rpo_name}'])
        #summarized potential of all inputs for one cell
        #cell_potential[i] += results[f'V_{target}/{rpo_name}']
#     #plt.plot(time_list, cell_potential[i], label = target)
#     '''
#     #Firing Rate
#     for pro_name in pro_names:
#          firing_rate[i] += results[f'Rate_{target}/{pro_name}']
#     plt.plot(time_list, firing_rate[i], label = target)
# '''
plt.xlabel('Time (s)')
#plt.ylabel('Potential (mV)')
plt.ylabel('Firing Rate (Hz)')
plt.legend()
plt.show()

