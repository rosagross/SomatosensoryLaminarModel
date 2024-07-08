from pyrates.frontend import OperatorTemplate, NodeTemplate, CircuitTemplate
from copy import deepcopy

print('hä')
# %%
# Operator template for the PRO

pro = OperatorTemplate(
    name='PRO', path=None,
    equations=["m_out = 2.*m_max / (1 + exp(r*(V_thr - v)))"],
    variables={'m_out': 'output',
               'v': 'input',
               'V_thr': 6e-3,
               'm_max': 2.5,
               'r': 560.0},
    description="sigmoidal potential-to-rate operator")

# %%
# Operator template for the RPO

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

rpo_ein = OperatorTemplate(
    name='RPO_ein', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * (m_in+u) - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'm_in': 'input',
               'tau': 0.01,
               'H': 0.00325,
               'u': 220},
    description="excitatory rate-to-potential operator")

# %%
# Node templates
ein = NodeTemplate(name="EIN", path=None, operators=[pro, rpo_e])
iin = NodeTemplate(name="IIN", path=None, operators=[pro, rpo_e])

rpo_i = deepcopy(rpo_e).update_template(
    name='RPO_i', path=None, variables={'H': -0.022, 'tau': 0.02}
)
pc = NodeTemplate(name="PC", path=None, operators=[pro, rpo_i, rpo_ein])
# Set up the Model Circuit 
jrc = CircuitTemplate(
    name="JRC", nodes={'PC': pc, 'EIN': ein, 'IIN': iin},
    edges=[("PC/PRO/m_out", "IIN/RPO_e/m_in", None, {'weight': 168.75}),
           ("PC/PRO/m_out", "EIN/RPO_e/m_in", None, {'weight': 675.}),
           ("EIN/PRO/m_out", "PC/RPO_ein/m_in", None, {'weight': 540.}),
           ("IIN/PRO/m_out", "PC/RPO_i/m_in", None, {'weight': 168.75})],
    path=None)

# Run the simulation 
results = jrc.run(simulation_time=2.0,
                  step_size=1e-4,
                  sampling_step_size=1e-3,
                  outputs={'V_PCE': 'PC/RPO_ein/v',
                           'V_PCI': 'PC/RPO_i/v'},
                  backend='default',
                  solver='scipy')

import matplotlib.pyplot as plt
plt.plot(results['V_PCE']+results['V_PCI'], label='PC Membrane potential')
for V in results:
    plt.plot(results[V],label=V)
    
plt.legend()
plt.show()
