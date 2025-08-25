"""Parameters sweep: simulations in time with different values of the parameters--> see PyRates example"""

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
# %%
""" Example with the dendritic time constants (\tau_e)
param_grid = {'tau_pc': [0.05, 0.1, 0.2], 'tau_ein': [0.05, 0.1, 0.2]} # values to "sweep" (defined as dictionary)  / the name is not related to the PyRates definition
param_map = {'tau_pc': {'vars': ['PC/RPO_e_in/tau'], 'nodes': ['PC']},
            'tau_ein': {'vars': ['EIN/RPO_e/tau'], 'nodes': ['EIN']}
            } # maps the parameters to the model 
"""

param_grid = {'p': [0, 60, 120]} # values to "sweep" (defined as dictionary)  / the name is not related to the PyRates definition
param_map = {'p': {'vars': ['PC/RPO_e_in/u'], 'nodes': ['PC']}
            } # maps the parameters to the model 

# %% Simulation values
T = 10.0
dt = 1e-4
dts = 1e-3
cutoff = 1.0
# input as noise --> defined if the varied parameter is not p 
# noise = np.random.uniform(120.0, 320.0, size=(int(np.round(T/dt, decimals=0)), 1))

# %% perform parameter sweep
results, results_map = grid_search(circuit_template=jrc,
                                   param_grid=param_grid,
                                   param_map=param_map,
                                   simulation_time=T,
                                   step_size=dt,
                                   sampling_step_size=dts,
                                   #inputs={'PC/RPO_e_in/u': noise},
                                   inputs = None, # the input is  the swept parameter 
                                   outputs={'V_pce': 'PC/RPO_e_in/v', 'V_pci': 'PC/RPO_i/v'},
                                   backend="default",
                                   solver="euler",
                                   #cutoff=cutoff
                                   )
"""permute_grid (default: False)
        If true, all combinations of the provided param_grid values will be realized. If false, the param_grid values
        will be traversed pairwise.
- for the other arguments see CircuitTemplate.run() """
# SAVE IN .CSV

# %%
fig, axes = plt.subplots(nrows=results_map.shape[0], figsize=(8, 12))

# sort the results map via the values of p
results_map.sort_values('p', inplace=True, axis=0)

# plot the raw output variable for each condition
for ax, key in zip(axes, results_map.index):
    psp_e = results['V_pce'][key].iloc[:, 0]
    psp_i = results['V_pci'][key].iloc[:, 0]
    ax.plot(psp_e, label=r"$V_{pce}$ (excitatory)", color="tab:blue")
    ax.plot(psp_i, label=r"$V_{pci}$ (inhibitory)", color="tab:orange")
    ax.plot(psp_e + psp_i, label=r"$V_{pc}$", color="tab:green")
    ax.set_title(f"p = {results_map.at[key, 'p']} Hz")
    ax.set_ylabel(r'$[V]$')
    ax.legend(fontsize="small", loc="lower right")
axes[-1].set_xlabel('time [s]')
plt.tight_layout()
plt.show()



# %%
