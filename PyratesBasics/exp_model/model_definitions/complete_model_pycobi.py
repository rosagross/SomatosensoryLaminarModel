# %%
import os 
os.chdir("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/model_definitions""") 
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
#/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/area_a3b_thal_bI_iext.yaml
model_path = "/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/model_definitions/model.yaml/model"
model = CircuitTemplate.from_yaml(model_path)
"""attrs = dir(my_object)
print(attrs)
"""
# %% atrributes and parameters
cells_ext = list(model.nodes.keys()) # getting the name of the populations directly from the model
cells = cells_ext[:-5]
N_cells = len(cells)

# names for the operators
pro_names = ["PRO_" + cell for cell in cells]
rpo_names = ["RPO_" + cell for cell in cells]
rpo_names_extended = rpo_names + ['RPO_Iext']

# sigmoid params for the rates
params = Parameter()
sigm = params.get_sigmoid() #already in the correct order

r = []
v_thr = []
m_max = []
for row in range(N_cells):
    r.append(sigm[row,0])      
    v_thr.append(sigm[row,1]) 
    m_max.append(sigm[row,2]) 

# %% Setting the simulation:

input_type = "step" # other options are "step", "baseline" (equals input strength 0) or "background"
input_onset = 1.0 # in sec
simulation_dur = 2.0
input_duration = 1.0  #, 1, 1.5] # np.arange(0, 1, 1) # in sec 
input_strength = 0.0 #[0, 50, 300, 500] #np.arange(0, 500, 100) - 10000
backgrndI_strengths = 0.0 #[0, 5, 10, 15, 20] - 10000
step_size=1e-3
sampling_step_size=1e-3
simulation_time = float(input_onset + simulation_dur)

# update of the model (might be different from the previous simulations)
model.update_var(node_vars={'thalE/InputOp/A': input_strength}) 
model.update_var(node_vars={'thalE/InputOp/dur': input_duration}) 
model.update_var(node_vars={'thalE/InputOp/onset': input_onset})
model.update_var(node_vars={'BACKGROUND/RPO_bI/bI': backgrndI_strengths})
#- background input:
# Build and update in the same loop
"""
for c in range(N_cells-2):
    key = cells[c] + '/RPO_bI/bI'
    area_a3b_thal_bI_iext.update_var(node_vars={key: backgrndI_strengths})

params = [cells[c] + '/RPO_bI/bI' for c in range(N_cells-2)]
"""
# %% PyCobi instance
model_auto = ODESystem.from_template(model, auto_dir = "/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p", init_cont=False)
pprint(model_auto._var_map)
cont_param = 'thalE/InputOp/A'
# %% Setting up the simulation:
# - time continuation starting at a steady state
# - initial value of the parameter: ?
# %% 
# original time continuation (same result as PyRates)
t_sols, t_cont = model_auto.run(
    c='ivp',  name='time', NPR=10, DS=1e-4, DSMIN=1e-6, EPSL=1e-05, EPSU=1e-05, EPSS=1e-03,
    DSMAX=1e-2, NMX=100000000, UZR={14: simulation_time}, STOP={'UZ1'})
# %%

"""
equilibrium_values = t_sols.iloc[[-1]]  # storing the equilibrium values

for i, target_cell in enumerate(cells):
    # always include rpo_names[:N_cells]
    if i == (N_cells-2):
        potential_keys = [f'{target_cell}/{rpo_name_ext}/v' for rpo_name_ext in rpo_names_extended[:N_cells+1]]
    else: 
        potential_keys = [f'{target_cell}/{rpo_name}/v' for rpo_name in rpo_names[:N_cells]]

    # include rpo_name only if i in range(N_cells-2)
    if i in range(N_cells-2):
        potential_keys += [f'{target_cell}/RPO_bI/v']

    eqs = np.array(equilibrium_values[potential_keys])

    for v, var in enumerate(potential_keys):
        area_a3b_thal_bI_iext.update_var(node_vars={f'{var}': eqs[v]})


for i, col in enumerate(t_sols.columns):
    print(f"Column {i}: {col}")

jrc.update_var(node_vars={'pc/rpo_e_in/u': 400.0})

# set state variables close to steady-state value
jrc.update_var(node_vars={'ein/rpo_e/v': 1.762243e-02 }) # y0
jrc.update_var(node_vars={'pc/rpo_e_in/v': 3.052549e-02}) # y1
jrc.update_var(node_vars={'pc/rpo_i/v': -2.195986e-02})
jrc.update_var(node_vars={'iin/rpo_e/v': 4.405607e-03 })


t_sols, t_cont = area_a3b_thal_bI_iext_auto.run(
    c='ivp',  name='time', DS=1e-3, DSMIN=1e-10, EPSL=1e-05, EPSU=1e-05, EPSS=1e-03,
    DSMAX=5e-3, NMX=99999999, UZR={14: 2.0}, STOP={'UZ1'})
    """
# %% plots of the potentials

all_potentials = []

for i, target_cell in enumerate(cells):
    # always include rpo_names[:N_cells]
    if i == (N_cells-2):
        potential_keys = [f'{target_cell}/{rpo_name_ext}/v' for rpo_name_ext in rpo_names_extended[:N_cells+1]]
    else: 
        potential_keys = [f'{target_cell}/{rpo_name}/v' for rpo_name in rpo_names[:N_cells]]

    # include rpo_name only if i in range(N_cells-2)
    if i in range(N_cells-2):
        potential_keys += [f'BACKGROUND/RPO_bI/v_bI']
    
    sources = t_sols[potential_keys]
    all_potentials.append(np.sum(sources, axis=1))

all_potentials = np.array(all_potentials).T
potential_df = pd.DataFrame(all_potentials, columns=cells)
# %%
# division into layers
layers = [potential_df.columns[:4],   # A3b 
          potential_df.columns[4:17],  # S1
          potential_df.columns[17:30],  # S2
          potential_df.columns[30:32] # thal
         ]
fig, axes = plt.subplots(2, 2, figsize=(14, 10))  
axes = axes.flatten()
for ax, cols in zip(axes, layers):
    # Create labels with final values
    labels_with_final = [f"{col} ({potential_df[col].iloc[-1]:.6f})" for col in cols]
    
    # Plot and override legend labels
    ax.plot(t_sols['t'], potential_df[cols])  # suppress default legend
    ax.legend(labels_with_final, loc="best")
    ax.set_ylabel("[mV]")
    ax.set_xlabel("time")
    ax.set_title(", ".join(cols))
#axes[-1].set_visible(False)
fig.suptitle("Potentials", fontsize=20, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()


# %% Bifurcation analysis:
# example - amplitude of the external input 'thalE/InputOp/A'
"""
u_sols, u_cont = area_a3b_thal_bI_iext_auto.run(
    origin=t_cont, 
    starting_point='UZ1', 
    name='u', # Name of continuation parameter (variable to vary) 
    bidirectional=True, #  Continue both forward and backward from starting point
    ICP=cont_param,
    RL0=0.0, RL1=500.0, # Lower and upper bounds for continuation parameter(s)
    IPS=1, # Problem type: 1=Equilibrium, 2=Periodic orbit continuation
    ILP=1, # Detect limit points (folds): 1=Yes, 0=No 
    ISP=2, # Bifurcation detection level: 0=none, 1=some, 2=more 
    ISW=1, 
    NTST=400,
    NCOL=4, # recommended in the auto docs
    IAD=3, # recommended in the auto docs
    IPLT=0, 
    NBC=0, 
    NINT=0, 
    NMX=1000000, 
    NPR=100, # Print/report frequency (every NPR steps)
    MXBF={}, # Maximum number of bifurcations to locate
    IID=2,
    ITMX=1000, # max.num. of iterations allowed in the accurate location of special solutions
    ITNW=40, 
    NWTN=12, 
    JAC=0, 
    EPSL=1e-07, 
    EPSU=1e-07, 
    EPSS=1e-05, # approx.100/1000 times EPSL,EPSU
    DS=1e-04, # initial continuation step size (0.5)
    DSMIN=1e-8, # minimum allowed step size
    DSMAX=1e-02, # maximum allowed step size (4)
    IADS=1, 
    THL={}, 
    THU={}, 
    UZR={}, 
    STOP={}
)
# %%
all_pot_cont = []

for i, target_cell in enumerate(cells):
    # always include rpo_names[:N_cells]
    if i == (N_cells-2):
        potential_keys = [f'{target_cell}/{rpo_name_ext}/v' for rpo_name_ext in rpo_names_extended[:N_cells+1]]
    else: 
        potential_keys = [f'{target_cell}/{rpo_name}/v' for rpo_name in rpo_names[:N_cells]]

    # include rpo_name only if i in range(N_cells-2)
    if i in range(N_cells-2):
        potential_keys += [f'BACKGROUND/RPO_bI/v_bI']
    
    sources = u_sols[potential_keys]
    all_pot_cont.append(np.sum(sources, axis=1))

all_pot_cont = np.array(all_pot_cont).T
pot_cont_df = pd.DataFrame(all_pot_cont, columns=cells)
layers = [pot_cont_df.columns[:4],   # first 4 
          pot_cont_df.columns[4:]  # next 4
         ]
fig, axes = plt.subplots(1, 2, figsize=(14, 10))  
axes = axes.flatten()
for ax, cols in zip(axes, layers):
    # Create labels with final values
    labels_with_final = [f"{col} ({pot_cont_df[col].iloc[-1]:.6f})" for col in cols]
    
    # Plot and override legend labels
    ax.plot(u_sols['thalE/InputOp/A'], pot_cont_df[cols])  # suppress default legend
    ax.legend(labels_with_final, loc="best")
    ax.set_ylabel("[mV]")
    ax.set_xlabel("Hz")
    ax.set_title(", ".join(cols))
#axes[-1].set_visible(False)
fig.suptitle("Potentials", fontsize=20, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

# %% Plot

all_pot_cont = []

for i, target_cell in enumerate(cells):
    # always include rpo_names[:N_cells]
    if i == (N_cells-2):
        potential_keys = [f'{target_cell}/{rpo_name_ext}/v' for rpo_name_ext in rpo_names_extended[:N_cells+1]]
    else: 
        potential_keys = [f'{target_cell}/{rpo_name}/v' for rpo_name in rpo_names[:N_cells]]

    # include rpo_name only if i in range(N_cells-2)
    if i in range(N_cells-2):
        potential_keys += [f'BACKGROUND/RPO_bI/v_bI']
    
    sources = u_sols[potential_keys]
    all_pot_cont.append(np.sum(sources, axis=1))

all_pot_cont = np.array(all_pot_cont).T
pot_cont_df = pd.DataFrame(all_pot_cont, columns=cells)
#u_sols = pd.concat([u_sols, pot_cont_df], axis=1)
#u_sols.loc[:, cells] = pot_cont_df.values
for cell in cells:
    u_sols[(cell, "pot_cont")] = pot_cont_df[cell].values

fig, axes = plt.subplots(N_cells, 1, figsize=(6, 3 * N_cells), sharex=True)

if N_cells == 1:
    axes = [axes]

for ax, cell in zip(axes, cells):
    area_a3b_thal_bI_iext_auto.plot_continuation(
        cont_param,
        cell,
        cont='u',
        ax=ax
    )
    ax.set_title(f"Continuation for {cell}")

plt.tight_layout()
plt.show()

"""
# %%