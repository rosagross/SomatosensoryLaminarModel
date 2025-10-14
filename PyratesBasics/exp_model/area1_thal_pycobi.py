# %%
import os 
os.chdir("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/""") 
from pyrates.frontend import OperatorTemplate, NodeTemplate, CircuitTemplate
from copy import deepcopy
from parameters import Parameter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numba import njit
from yaml_saving import circuit_to_yaml
from pprint import pprint
from pycobi import ODESystem

# %% import of the .yaml file
# the path should end with the name of the file and then the name of the model
#/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/area_1_thal_iext.yaml
model_path = "/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/area_1_thal_iext.yaml/area_1_thal_iext"
area_1_thal = CircuitTemplate.from_yaml(model_path)
"""attrs = dir(my_object)
print(attrs)
"""
# %% atrributes and parameters
cells = list(area_1_thal.nodes.keys()) # getting the name of the populations directly from the model
N_cells = len(cells)

# names for the operators
pro_names = ["PRO_" + cell for cell in cells]
rpo_names = ["RPO_" + cell for cell in cells]

# sigmoid params for the rates
params = Parameter()
sigmoid_params = params.get_sigmoid() #already in the correct order
sigm = np.vstack((sigmoid_params[4:17], sigmoid_params[-2:])) # area1 columns + last two for thalE and thalI

r = []
v_thr = []
m_max = []
for row in range(N_cells):
    r.append(sigm[row,0])      
    v_thr.append(sigm[row,1]) 
    m_max.append(sigm[row,2]) 

"""# %% input definition
input_type = "step" # other options are "step", "baseline" (equals input strength 0) or "background"
input_onset = 1.001 # in sec
simulation_dur = 1 
input_duration = 0.5  #, 1, 1.5] # np.arange(0, 1, 1) # in sec 
input_strength = 50 #[0, 50, 300, 500] #np.arange(0, 500, 100)
backgrndI_strengths = 0 #[0, 5, 10, 15, 20]
step_size=1e-3
sampling_step_size=1e-3
simulation_time = int(input_onset) + simulation_dur

def create_Iext(simulation_time, step_size, input_onset, input_duration, input_strength, input_type):
    ''' Creates external input.
    Inputs in the function
    - simulation_time:
    - step_size:
    - input_onset:
    - input_duration:
    - input_strength: 
    - input_type:
    '''

    Iext = np.zeros(int(simulation_time/step_size))
    if input_type == "step":
        t  = int(input_duration/step_size)
        t0 = int(input_onset/step_size)
        Iext[t0:t0+t] = input_strength
    elif input_type == "background":
        # provide input for the entire simulation duration 
        Iext[:] = input_strength

    return Iext
Iext = create_Iext(simulation_time, step_size, input_onset, input_duration, input_strength, input_type)
"""

# %% PyCobi instance
area_1_thal_auto = ODESystem.from_template(area_1_thal, auto_dir = "/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p", init_cont=False)
pprint(area_1_thal_auto._var_map)
# %%
t_sols, t_cont = area_1_thal_auto.run(
    c='ivp',  name='time', DS=1e-4, DSMIN=1e-10, EPSL=1e-05, EPSU=1e-05, EPSS=1e-03,
    DSMAX=1e-3, IADS=1, NMX=50000, UZR={14: 2.0}, STOP={'UZ1'})

# %%
rpo_names_extended = rpo_names + ['RPO_Iext']
all_potentials = []
for i in cells:
    if i != 'thalE':
        sources = t_sols[[f'V_{i}/{rpo_name}' for rpo_name in rpo_names[:N_cells]]]
    else:
        sources = t_sols[[f'V_{i}/{rpo_name_ext}' for rpo_name_ext in rpo_names_extended[:N_cells+1]]]
    all_potentials.append(np.sum(sources, axis=1))
all_potentials = np.array(all_potentials).T # shape: samples x populations

# %%
layers = [potential_df.columns[:4],   # first 4 
          potential_df.columns[4:7],  # next 3
          potential_df.columns[7:10], # next 3
          potential_df.columns[10:13],  # next 3
          potential_df.columns[13:]]  # thalamus

fig, axes = plt.subplots(2, 3, figsize=(14, 10))  # 2x2 grid
axes = axes.flatten()

for ax, cols in zip(axes, layers):
    # Create labels with final values
    labels_with_final = [f"{col} ({potential_df[col].iloc[-1]:.6f})" for col in cols]
    
    # Plot and override legend labels
    potential_df[cols].plot(ax=ax, legend=False)  # suppress default legend
    ax.legend(labels_with_final, loc="best")
    ax.set_ylabel("[mV]")
    ax.set_xlabel("samples")
    ax.set_title(", ".join(cols))
axes[-1].set_visible(False)
fig.suptitle("Potentials", fontsize=20, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])

plt.show()

# %%
