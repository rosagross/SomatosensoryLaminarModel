# %%
from pyrates.frontend import OperatorTemplate, NodeTemplate, CircuitTemplate
from copy import deepcopy
import matplotlib.pyplot as plt
from parameters import Parameter
import numpy as np

# TODO:
# implement parameters.py!


#Parameter
N_cells = 5
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
pro_names = ['PRO_E1', 'PRO_P1', 'PRO_S1', 'PRO_V1', 'PRO_5']
rpo_names = ['RPO_E1', 'RPO_P1', 'RPO_S1', 'RPO_V1', 'RPO_5']

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


cells = ['E1', 'P1', 'S1', 'V1','5']

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
results = jrc.run(simulation_time=2.0,
                  step_size=1e-4,
                  sampling_step_size=1e-3,
                  outputs={'V_E1': 'S1/RPO_E1/v', #alle Eingänge
                           'V_P1': 'S1/RPO_P1/v',
                           'V_S1': 'S1/RPO_S1/v',
                           'V_V1': 'S1/RPO_V1/v'},
                  backend='default',
                  solver='scipy',
                  vectorize = False)


results.plot()
plt.show()
