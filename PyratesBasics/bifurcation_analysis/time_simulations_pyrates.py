# %%
from pyrates import grid_search
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
import pandas as pd
from pyrates.frontend import OperatorTemplate, NodeTemplate, CircuitTemplate
from copy import deepcopy
import math
import os 
os.chdir("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/bifurcation_analysis/""")

# %% ANALYSIS OF A PYRATES MODEL
# jrc.__dict__.keys() or here: https://pyrates.readthedocs.io/en/latest/frontend.html

# %% Model definition --> Pyrates

pro = OperatorTemplate(
    name='PRO', path=None,
    equations=["m_out = 2.*m_max / (1 + exp(r*(V_thr - v)))"],
    variables={'m_out': 'output',
               'v': 'input',
               'V_thr': 6e-3,
               'm_max': 2.5,
               'r': 560.0},
    description="sigmoidal potential-to-rate operator")

rpo_e = OperatorTemplate(
    name='RPO_e', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * m_in - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'm_in': 'input',
               'tau': 0.01,
               'H': 0.00325},
    description="excitatory rate-to-potential operator")

rpo_e_in = deepcopy(rpo_e).update_template(
    name='RPO_e_in', 
    path=None, 
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * (m_in+u) - 2 * i/tau - v/tau^2'], 
    variables={'u': 0.0},
    description="excitatory rate-to-potential operator with extrinsic input"
)

rpo_i = deepcopy(rpo_e).update_template(
    name='RPO_i', path=None, variables={'H': -0.022, 'tau': 0.02}
)

ein = NodeTemplate(name="EIN", path=None, operators=[pro, rpo_e])
iin = NodeTemplate(name="IIN", path=None, operators=[pro, rpo_e])
pc = NodeTemplate(name="PC", path=None, operators=[pro, rpo_e_in, rpo_i])

jrc = CircuitTemplate(
    name="JRC", nodes={'PC': pc, 'EIN': ein, 'IIN': iin},
    edges=[("PC/PRO/m_out", "IIN/RPO_e/m_in", None, {'weight': 33.75}),
           ("PC/PRO/m_out", "EIN/RPO_e/m_in", None, {'weight': 135.}),
           #("EIN/PRO/m_out", "PC/RPO_e/m_in", None, {'weight': 108.}),
           ("EIN/PRO/m_out", "PC/RPO_e_in/m_in", None, {'weight': 108.}),
           ("IIN/PRO/m_out", "PC/RPO_i/m_in", None, {'weight': 33.75})],
    path=None)

# jrc.update_var(node_vars={'pc/rpo_e_in/u': 0.0})

# %%
# Run the simulation 
results = jrc.run(simulation_time=2.0,
                  step_size=1e-4,
                  sampling_step_size=1e-3,
                  outputs={'V_PCE': 'PC/RPO_e_in/v',
                          'V_PCI': 'PC/RPO_i/v'},
                           # 'PC_m': 'PC/PRO/m_out'
                  backend='default',
                  solver='scipy')

# %%
import matplotlib.pyplot as plt
plt.plot(results['V_PCE']+results['V_PCI'], label='PC Membrane potential')
for V in results:
    plt.plot(results[V],label=V)
    
plt.legend()
plt.show()

# %%
