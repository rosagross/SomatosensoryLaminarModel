import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy  
from pyrates.frontend import OperatorTemplate
from pyrates.frontend import NodeTemplate
from neu_parameters import Parameter
from pyrates.frontend import CircuitTemplate
from pyrates import grid_search


params = Parameter()
tau, N_cells = params.get_params()
sigm = params.get_sigmoid()

pro_names = np.array(['PRO_E1', 'PRO_E2', 'PRO_E3', 'PRO_E4', 'PRO_P1', 'PRO_P2', 'PRO_P3', 'PRO_P4', 'PRO_S1', 'PRO_S2', 'PRO_S3', 'PRO_S4', 'PRO_V1', 'PRO_V2', 'PRO_V3', 'PRO_V4'])
pro = np.array([])

for i in range(len(pro_names)):
    pro = np.append(pro, OperatorTemplate(
        name = f'{pro_names[i]}', 
        path = None, 
        equations = ["m_out = m_max/(1 + exp(r*(vth-v)))"], 
        variables = {'m_out': 'output',
                     'v': 'input',
                     'r': sigm[i,0],
                     'vth': sigm[i,1],
                     'm_max': sigm[i,2]}, 
        description = "sigmoidal potential-to-rate operator"
    ))

rpo_names = np.array(['RPO_E1', 'RPO_E2', 'RPO_E3', 'RPO_E4', 'RPO_P1', 'RPO_P2', 'RPO_P3', 'RPO_P4', 'RPO_S1', 'RPO_S2', 'RPO_S3', 'RPO_S4', 'RPO_V1', 'RPO_V2', 'RPO_V3', 'RPO_V4'])#, 'RPO_INPUT'])
rpo = np.array([])

W = params.get_connectivity(1)

for i in range(len(rpo_names)): 
        rpo = np.append(rpo, OperatorTemplate(
            name = f'{rpo_names[i]}', 
            path = None, 
            equations = ['d/dt * v = u',
                        'd/dt * u = H/tau * r - 2 * u/tau - v/(tau**2)'],  
            variables = {'v': 'output',
                        'u': 'variable',
                        'r': 'input',
                        'tau': tau[0,i],
                        'H': 1},
            description = "rate-to-potential operator"
        ))


cells = np.array(['E1', 'E2', 'E3', 'E4', 'P1', 'P2', 'P3', 'P4', 'S1', 'S2', 'S3', 'S4', 'V1', 'V2', 'V3', 'V4'])


E1 = NodeTemplate(name="E1", path=None, operators=[pro[0]] + list(rpo)) #alle Eingänge

pop = np.array([E1])

for i in range(1,N_cells):
    pop = np.append(pop, deepcopy(E1).update_template(
        name = f'{cells[i]}', operators=[pro[i]] + list(rpo)
    ))
pop = list(pop)

 
updated_nodes={}
for i in range(N_cells):
    updated_nodes[cells[i]] = pop[i]
edges=[]
# i : target 
for i, cell_i in enumerate(cells):
    # j : source
    for j, cell_j in enumerate(cells):
        edges.append((f'{cell_j}/{pro[j].name}/m_out', f'{cell_i}/{rpo[j].name}/r', None, {'weight': W[i,j]})) # to x from 
        

cir = CircuitTemplate( 
    name="cir", 
    nodes=updated_nodes,
    edges=edges,
    path = None 
)


outputs = {}

#Membrane potential
for target_cell in cells:  
    for rpo_name in rpo_names:
        outputs[f'V_{target_cell}/{rpo_name}'] = f'{target_cell}/{rpo_name}/v'


T = 2.0 #das muss irgendwie 2s sein. Nur 1s geht nicht
dt = 1e-4 #numerische Integration  1e-4
sampling_step_size = 1e-3  #Datenaufzeichnung   1e-3

######################################

g = [20, 0.8, 1]

#W_g1 = list(params.get_connectivity(g[0])[:,:-1].flatten())
W_g1 = params.get_connectivity(g[0]) #(16x17)
W_g2 = params.get_connectivity(g[1])
W_g3 = params.get_connectivity(g[2])

all_edges = [(f'{cell_j}/{pro_names[j]}/m_out', f'{cell_i}/{rpo[j].name}/r') for j, cell_j in enumerate(cells) for cell_i in cells]
#print(all_edges)

#with one edge for all three gs
results_g, params = grid_search(cir,
                                param_grid = {'g0': [W_g1[4,0], W_g2[4,0], W_g3[4,0]]},  #3 verschiedene g - Werte für ein edge
                                param_map = {'g0': {'vars': ['weight'],
                                                    'edges': all_edges} 
                                                    },
                                outputs = {'V_out': 'E1/RPO_P1/v'}, #outputs
                                step_size=dt, simulation_time=T, sampling_step_size=sampling_step_size,
                                cutoff=1.0, permute_grid=False, vectorize = False)

print(results_g.shape)



# time_list = np.linspace(0,T,len(results_g))
# cell_potential = np.zeros((N_cells, len(time_list)))
# firing_rate = np.zeros((N_cells, len(time_list)))
# plt.figure()
# for i,target in enumerate(cells):
#     for rpo_name in rpo_names:
#         #summarized potential of all inputs for one cell
#         cell_potential[i] += results_g[f'V_{target}/{rpo_name}']
#     plt.plot(time_list, cell_potential[i], label = target)
#plt.legend()
#

results_g.plot()
plt.show()
