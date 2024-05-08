# %%
from pyrates.frontend import OperatorTemplate, NodeTemplate, CircuitTemplate
from copy import deepcopy
from parameters import Parameter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#Parameter
N_cells = 16
cortex_type = 'visual'
params = Parameter(cortex_type)
sigm = params.get_sigmoid()


r = []
v_thr = []
m_max = []
for row in range(N_cells):
    r.append(sigm[row,0])      
    v_thr.append(sigm[row,1]) 
    m_max.append(sigm[row,2]) 

tau,_ = params.get_params()                  
tau = tau[0]


W = params.get_connectivity(1, include_Iext=False)
# %%
# Operator template for the PRO

pro_names = ['PRO_E1', 'PRO_E2', 'PRO_E3', 'PRO_E4',
             'PRO_P1', 'PRO_P2', 'PRO_P3', 'PRO_P4',
             'PRO_S1', 'PRO_S2', 'PRO_S3', 'PRO_S4',
             'PRO_V1', 'PRO_V2', 'PRO_V3', 'PRO_V4']

rpo_names = ['RPO_E1', 'RPO_E2', 'RPO_E3', 'RPO_E4',
             'RPO_P1', 'RPO_P2', 'RPO_P3', 'RPO_P4',
             'RPO_S1', 'RPO_S2', 'RPO_S3', 'RPO_S4',
             'RPO_V1', 'RPO_V2', 'RPO_V3', 'RPO_V4']

pro = OperatorTemplate(
    name='PRO', path=None,
    equations=["m_out = 2.*m_max / (1 + exp(r*(V_thr - (v*1000))))"],
    variables={'m_out': 'output',
               'v': 'input',
               'r': 0.1,
               'V_thr': 35,
               'm_max': 70},
    description="sigmoidal potential-to-rate operator")

pros = [deepcopy(pro).update_template(name=pro_names[i], 
                                       variables={'r': r[i], 
                                                  'V_thr': v_thr[i],
                                                  'm_max': m_max[i]}) for i in range(N_cells)] 

#print([pros[i].variables for i in range(4)]) 
# %%
# Operator template for the RPO

rpo = OperatorTemplate(
    name='RPO', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * m_in - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'm_in': 'input',
               'tau': 0.01,
               'H': 1},
    description="excitatory rate-to-potential operator")

rpos = [deepcopy(rpo).update_template(name=rpo_names[i], variables={"tau": tau[i]}) for i in range(N_cells)]



# %%
# Node templates
node = NodeTemplate(name="node", path=None, operators= pros + rpos)


cells = ['E1', 'E2', 'E3', 'E4', 'P1', 'P2', 'P3', 'P4', 'S1', 'S2', 'S3', 'S4', 'V1', 'V2', 'V3', 'V4']

# Define edges
edges=[]
# i : target 
for i, cell_i in enumerate(cells):
    # j : source
    for j, cell_j in enumerate(cells):
        edges.append((f'{cell_j}/{pro_names[j]}/m_out', f'{cell_i}/{rpo_names[j]}/m_in', None, {'weight': W[i,j]})) # to x from 

# Set up the Model Circuit 
jrc = CircuitTemplate(
    name = 'cir',
    nodes = {name: node for name in cells},
    edges = edges,
    path = None)

# Run the simulation 
outputs = {}
for target_cell in cells:  
    for rpo_name in rpo_names:
        outputs[f'V_{target_cell}/{rpo_name}'] = f'{target_cell}/{rpo_name}/v'

results = jrc.run(simulation_time=2.0,
                  step_size=1e-4,
                  sampling_step_size=1e-4,
                  outputs= outputs, #{'V_E1': 'V4/RPO_E1/v', #Eingänge nach P1
                                    # 'V_P1': 'V4/RPO_E2/v',
                                    # 'V_S1': 'V4/RPO_E3/v',
                                    # 'V_V1': 'V4/RPO_E4/v'},
                  backend='default',
                  solver='scipy',
                  vectorize = False)


all_potentials = []
for i in cells:
    sources = results[[f'V_{i}/{rpo_name}' for rpo_name in rpo_names]]
    all_potentials.append(np.sum(sources, axis=1))

all_potentials = np.array(all_potentials).T
potential_df = pd.DataFrame(all_potentials, columns=cells)
potential_df.to_csv('output/pyrates_upscaled.csv', index=False)

results.plot()
plt.show()
