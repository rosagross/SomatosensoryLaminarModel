# %%
from pyrates.frontend import OperatorTemplate, NodeTemplate, CircuitTemplate
from copy import deepcopy
from parameters import Parameter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#Parameter
N_cells = 11
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
             'PRO_S1', 'PRO_S2', 'PRO_S3']

rpo_names = ['RPO_E1', 'RPO_E2', 'RPO_E3', 'RPO_E4',
             'RPO_P1', 'RPO_P2', 'RPO_P3', 'RPO_P4',
             'RPO_S1', 'RPO_S2', 'RPO_S3']

pro = OperatorTemplate(
    name='PRO', path=None,
    equations=["m = 2.*m_max / (1 + exp(r*(V_thr - (v_in*1000))))"],
    variables={'m': 'output',
               'v_in': 'input',
               'r': 0.1,
               'V_thr': 35,
               'm_max': 70},
    description="sigmoidal potential-to-rate operator")
pros = [deepcopy(pro).update_template(name=pro_names[i],
                                      variables={'r': r[i],
                                                 'V_thr': v_thr[i],
                                                 'm_max': m_max[i]}) for i in range(N_cells)]

# Operator template for the RPO
rpo = OperatorTemplate(
    name='RPO', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * m - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'm': 'input',
               'tau': 0.01,
               'H': 1},
    description="excitatory rate-to-potential operator")
rpos = [deepcopy(rpo).update_template(name=rpo_names[i], variables={"tau": tau[i]}) for i in range(N_cells)]

# Node templates
cells = ['E1', 'E2', 'E3', 'E4', 'P1', 'P2', 'P3', 'P4', 'S1', 'S2', 'S3']
nodes = {c: NodeTemplate(name=c, path=None, operators= [pros[i], rpos[i]]) for i, c in enumerate(cells)}

# Define edges
edges=[]
# i : target 
for i, cell_i in enumerate(cells):
    # j : source
    for j, cell_j in enumerate(cells):
        edges.append((f'{cell_j}/{rpo_names[j]}/v', f'{cell_i}/{pro_names[i]}/v_in', None, {'weight': W[i, j]})) # to x from

# Set up the Model Circuit 
jrc = CircuitTemplate(
    name = 'cir',
    nodes = nodes,
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
                  outputs= outputs,
                  backend='default',
                  #solver= 'scipy',
                  #method = 'RK25',
                  vectorize = True,
                  clear=False)

# Spaltenüberschriften umbenennen
results.columns = cells
print(results)
results.to_csv('output/pyrates_11_copy.csv', index=False)

results.plot()
plt.show()

