import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy  
from pyrates.frontend import OperatorTemplate
from neu_parameters import Parameter

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
    for j in range(N_cells):
        rpo = np.append(rpo, OperatorTemplate(
            name = f'{rpo_names[i]}', 
            path = None, 
            equations = ['d/dt * v = u',
                        'd/dt * u = H/tau * r - 2 * u/tau - v/(tau**2)'],  
            variables = {'v': 'output',
                        'u': 'variable',
                        'r': 'input',
                        'tau': tau[0,i],
                        'H': 1,
                        'w': W[i,j]},
            description = "rate-to-potential operator"
        ))


cells = np.array(['E1', 'E2', 'E3', 'E4', 'P1', 'P2', 'P3', 'P4', 'S1', 'S2', 'S3', 'S4', 'V1', 'V2', 'V3', 'V4'])

from pyrates.frontend import NodeTemplate
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

        
from pyrates.frontend import CircuitTemplate
cir = CircuitTemplate( 
    name="cir", 
    nodes=updated_nodes,
    edges=edges,
    path = None 
)

outputs = {}
'''
#Membrane potential
for target_cell in cells:  
    for rpo_name in rpo_names:
        outputs[f'V_{target_cell}/{rpo_name}'] = f'{target_cell}/{rpo_name}/v'
        '''
#Firing Rate
for source_cell in cells: 
    for pro_name in pro_names:
        outputs[f'Rate_{source_cell}/{pro_name}'] = f'{source_cell}/{pro_name}/m_out'


simulation_time = 2.0
step_size = 1e-4
sampling_step_size = 1e-3

results = cir.run(simulation_time = simulation_time,
                  step_size = step_size,
                  sampling_step_size=sampling_step_size,
                  outputs = outputs,
                  backend ='default',
                  solver ='scipy',
                  vectorize=False)



time_list = np.arange(0, simulation_time, sampling_step_size)
print(np.shape(time_list))
print(np.shape(results))


cell_potential = np.zeros((N_cells, len(time_list)))
firing_rate = np.zeros((N_cells, len(time_list)))
for i,target in enumerate(cells):
    #Membrane Potential
    '''
    for rpo_name in rpo_names:
        #plt.plot(time_list, results[f'V_{target}/{rpo_name}'])
        #summarized potential of all inputs for one cell
        cell_potential[i] += results[f'V_{target}/{rpo_name}']
    plt.plot(time_list, cell_potential[i], label = target)
    '''
    #Firing Rate
    for pro_name in pro_names:
         firing_rate[i] += results[f'Rate_{target}/{pro_name}']
    plt.plot(time_list, firing_rate[i], label = target)

plt.xlabel('Time (s)')
#plt.ylabel('Potential (mV)')
plt.ylabel('Firing Rate (Hz)')
plt.legend()
plt.show()


