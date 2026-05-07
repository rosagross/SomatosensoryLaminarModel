# %% libraries import
import os 
WDDIR = os.getenv("WDDIR")
SIMDIR = os.getenv("SIMDIR")
os.chdir(os.path.join(WDDIR,"PyratesBasics","exp_model")) 
import sys
import json
from ruamel.yaml import YAML
from pyrates.frontend.template import CircuitTemplate
from pyrates.frontend.fileio.yaml import dump_to_yaml
from pyrates.frontend import OperatorTemplate, NodeTemplate, EdgeTemplate, CircuitTemplate
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import h5py
import pandas as pd
from pprint import pprint


# %% input parameters
input_type = "step" # other options are "step", "baseline" (equals input strength 0) or "background"
input_onset = 1.001 # in sec
simulation_dur = 2.0
Iext_duration  = 1.0  #, 1, 1.5] # np.arange(0, 1, 1) # in sec 
Iext_strength = 50  #[0, 50, 300, 500] #np.arange(0, 500, 100)
backgrndI_strengths = 5 #[0, 5, 10, 15, 20]
step_size=1e-3
sampling_step_size=1e-3
simulation_time = float(int(input_onset) + simulation_dur)

# %% I_ext check
"""
Creates external input.
THIS IS ONLY FOR PLOTTING REASONS. The actual input is computed via operators.
"""
Iext = np.zeros(int(simulation_dur / step_size))
t = int(Iext_duration / step_size)
t0 = int(input_onset / step_size)
Iext[t0 : t0 + t] = Iext_strength
        
# %%PyRates model
create_I_ext = OperatorTemplate(
            name="InputOp", path=None,
            equations=[
            # step input: A during [onset, onset+dur), else 0
                #"Iext = (A/2)*(1+sign(t-onset)) - (A/2)*(1+sign(t-(onset+dur)))"
                "Iext = (A/2)*(1+sign(t-onset)) - (A/2)*(1+sign(t-(onset+dur)))"
                ],
            variables={
                "Iext": "output",
                 "t": "variable",
                "A": Iext_strength,
                "onset": input_onset/step_size,
                "dur": Iext_duration/step_size
                # TO UNCOMMENT FOR PYCOBI:
                #"A": float(input_strength),
                #"onset": float(input_onset),
                #"dur": float(self.Iext_duration)
                },
            description="External step input"
            )


iin = NodeTemplate(name="in", path=None, operators=[create_I_ext])


input_example = CircuitTemplate(
    name="IN", nodes={'IIN': iin},
    edges=[ ], # ("IIN/PRO/m_out", "PC/RPO_i/m_in", None, {'weight': 33.75})
    path=None)

input_example.get_run_func(func_name='Iext_definition', step_size=step_size, file_name='Iext_definition_prova1', backend='numpy',
                 vectorize=False, float_precision='float64')

# %%
# Run the simulation 
results = input_example.run(simulation_time=2.0,
                  step_size=step_size,
                  sampling_step_size=sampling_step_size,
                  outputs={'IIN_out': 'in/InputOp/I_ext'},
                  backend='default',
                  solver='scipy')

# %% 
import matplotlib.pyplot as plt
plt.plot(results['V_PCE']+results['V_PCI'], label='PC Membrane potential')
for V in results:
    plt.plot(results[V],label=V)
    
plt.legend()
plt.show()
