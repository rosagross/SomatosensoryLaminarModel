import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy  
from pyrates.frontend import OperatorTemplate
from pyrates.frontend import NodeTemplate
from parameters import Parameter
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
dt = 1e-3 #numerische Integration  1e-4
sampling_step_size = dt  #Datenaufzeichnung   1e-3


g = [0.5,1,20]


all_edges = [(f'{cell_j}/{pro_names[j]}/m_out', f'{cell_i}/{rpo[j].name}/r') for j, cell_j in enumerate(cells) for cell_i in cells]


all_results = []
all_parameters = []

for g_val in g:
     W_g = list(params.get_connectivity(g_val)[:,:-1].flatten(order = 'F')) # -1 wegen Inputspalte (16x17)

     results_g, parameter_map = grid_search(cir,
                                    param_grid = {f'g{[i]}': [W_g[i]] for i in range(len(W_g))}, 
                                    param_map = {f'g{[i]}': {'vars': ['weight'],
                                                        'edges': [all_edges[i]]} for i in range(len(W_g))
                                                        },
                                    outputs = outputs,
                                    step_size=dt, simulation_time=T, sampling_step_size=sampling_step_size,
                                    cutoff=0, permute_grid=False, vectorize = False)
     
     all_results.append(results_g)
     all_parameters.append(parameter_map)

# results in csv file speichern

time_list = np.linspace(0,T,len(results_g))
cell_potential = np.zeros((len(g), N_cells, len(results_g))) #g-Wert x Zellen x Zeitschritte
firing_rate = np.zeros((len(g), N_cells, len(results_g)))
fr_g = np.zeros((len(g), N_cells))


fig, ax = plt.subplots(3,1)
for p, res in enumerate(all_results):
    for i,target in enumerate(cells):
        for rpo_name in rpo_names:
            cell_potential[p,i] += res[f'V_{target}/{rpo_name}']

            r = sigm[i,0]
            vth = sigm[i,1]
            m_max = sigm[i,2]
            firing_rate[p,i] = m_max/(1 + np.exp(r*(vth-cell_potential[p,i])))
        
        
        
        
        #ax[p].plot(time_list, cell_potential[p,i], label = target)
        #ax[p].legend()
        ax[p].plot(time_list, firing_rate[p,i], label = target)
        ax[p].legend()
plt.show()

print(firing_rate)


fr_g = np.mean(firing_rate[:,:,-50:],axis = 2) #Durchschnitt der letzten 50 Feuerraten für jede Feuerrate
print(fr_g)

#plt.plot(g, fr_g)
#plt.show()





# alle Membranpotentiale für jede Population zusammenrechnen
# danach letztes Membranpotential umrechnen in Feuerrate und abhängig von g plotten
# Cooler wäre es, direkt auf die Feuerrate zugreifen zu können und z.B. den Durchschnitt der letzten 50 Feuerraten abhängig von g zu plotten
# bevor ich einen Input hinzufüge, würde ich lieber diesen Code besser und schneller haben
