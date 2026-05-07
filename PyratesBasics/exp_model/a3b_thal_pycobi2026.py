# %% 
"""
Bifurcation analysis of the "reduced model" (A3b)
14-02-2026
"""
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
#import h5py
import pandas as pd
from pprint import pprint
from pycobi import ODESystem

#from parameters import Parameter
param_path = os.path.join(WDDIR,"Simulations")

if param_path not in sys.path:
    sys.path.append(param_path)
from parameters import Parameter
import json

# %% import of the .yaml file
# the path should end with the name of the file and then the name of the model
#/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/area_a3b_thal_bI_iext.yaml

model_path = "/data/p_02989/Modelling/mecozzi_wd/SomatosensoryLaminarModel/PyratesBasics/exp_model/area_a3b_thal_bI_iext_2026.yaml/area_a3b_thal_bI_iext"
#model_path = "/data/p_02989/Modelling/mecozzi_wd/SomatosensoryLaminarModel/PyratesBasics/exp_model/area_a3b_thal_bI_iextneq_2026.yaml/area_a3b_thal_bI_iext"

area_a3b_thal_bI_iext = CircuitTemplate.from_yaml(model_path)
"""attrs = dir(my_object)
print(attrs)
"""
# %% attributes and parameters + equilibrium values
cells = [c for c in area_a3b_thal_bI_iext.nodes.keys() if c != "background_input"] # getting the name of the populations directly from the model + removing "background"
N_cells = len(cells)
pro_names = ["PRO_"+ cell for cell in cells]
rpo_names=["RPO_"+cell for cell in cells]
rpo_names_extended = rpo_names + ['RPO_Iext']

#equilibrium_df = pd.read_csv("equilibrium_df_A3b.csv")

# %% Setting the simulation:

input_type = "step" # other options are "step", "baseline" (equals input strength 0) or "background"
input_onset = 1.0 # in sec
simulation_dur = 2.0
input_duration = 1  #, 1, 1.5] # np.arange(0, 1, 1) # in sec 
input_strength = 0.0 #[0, 50, 300, 500] #np.arange(0, 500, 100)
backgrndI_strengths = 0 #[0, 5, 10, 15, 20]
step_size=1e-3
sampling_step_size=1e-3
#simulation_time = 1 # float(int(input_onset) + simulation_dur)

# %% Equilibrium values as starting points for the time continuation
equilibrium_df = pd.read_csv("equilibrium_df_A3b.csv")
print(equilibrium_df.columns.tolist())
#area_a3b_thal_bI_iext.update_var(node_vars=equilibrium_df.iloc[0].to_dict())
# %%
for var in equilibrium_df.columns:
    area_a3b_thal_bI_iext.update_var(node_vars={var: equilibrium_df[var].values[0]})
# %%
"""
# update of the model (might be different from the previous simulations)
area_a3b_thal_bI_iext.update_var(node_vars={'thalE/InputOp/A': input_strength}) 
area_a3b_thal_bI_iext.update_var(node_vars={'thalE/InputOp/dur': input_duration}) 
area_a3b_thal_bI_iext.update_var(node_vars={'thalE/InputOp/onset': input_onset})

#- background input:
# Build and update in the same loop

for c in range(N_cells-2):
    key = cells[c] + '/RPO_bI/bI'
    area_a3b_thal_bI_iext.update_var(node_vars={key: backgrndI_strengths})

params = [cells[c] + '/RPO_bI/bI' for c in range(N_cells-2)]
"""
# %% PyCobi instance
area_a3b_thal_bI_iext_auto = ODESystem.from_template(area_a3b_thal_bI_iext, auto_dir = "/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p", init_cont=False)
pprint(area_a3b_thal_bI_iext_auto._var_map) # "available" parameters
# %%
cont_param = 'background_input/bI_definition/backgrndI_strengths'

# %% Setting up the simulation:

simulation_dur = 2.0
# original time continuation (same result as PyRates)
t_sols, t_cont = area_a3b_thal_bI_iext_auto.run(
    c='ivp',  name='time', DS=1e-3, DSMIN=1e-5, EPSL=1e-05, EPSU=1e-05, EPSS=1e-03,
    DSMAX=1e-1, NMX=1000000, UZR={14: simulation_dur}, STOP={'UZ1'})
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
        potential_keys += [f'{target_cell}/RPO_bI/v']
        
    
    sources = t_sols[potential_keys]
    all_potentials.append(np.sum(sources, axis=1))

all_potentials = np.array(all_potentials).T
potential_df = pd.DataFrame(all_potentials, columns=cells)
# %%
# division into layers
layers = [potential_df.columns[:4],   # first 4 
          potential_df.columns[4:]  # next 4
         ]
fig, axes = plt.subplots(1, 2, figsize=(14, 10))  
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

u_sols, u_cont = area_a3b_thal_bI_iext_auto.run(
    origin=t_cont, 
    starting_point='UZ', 
    name='u', # Name of continuation parameter (variable to vary) 
    bidirectional=False, #  Continue both forward and backward from starting point
    ICP=cont_param,
    RL0=-5.0, RL1=5.0, # Lower and upper bounds for continuation parameter(s)
    IPS=1, # Problem type: 1=Equilibrium, 2=Periodic orbit continuation
    ILP=1, # Detect limit points (folds): 1=Yes, 0=No 
    ISP=2, # Bifurcation detection level: 0=none, 1=some, 2=more 
    get_eigenvals=True,
    get_stability=True,
    ISW=1, 
    NTST=400,
    NCOL=4, # recommended in the auto docs
    IAD=3, # recommended in the auto docs
    IPLT=0, 
    NBC=0, 
    NINT=0, 
    NMX=10000000, 
    NPR=100, # Print/report frequency (every NPR steps)
    MXBF={}, # Maximum number of bifurcations to locate
    IID=2,
    ITMX=10000, # max.num. of iterations allowed in the accurate location of special solutions
    ITNW=40, 
    NWTN=12, 
    JAC=0, 
    EPSL=1e-07, 
    EPSU=1e-07, 
    EPSS=1e-05, # approx.100/1000 times EPSL,EPSU
    DS=-1e-07, # initial continuation step size (0.5)
    DSMIN=1e-09, # minimum allowed step size
    DSMAX=1e-05, # maximum allowed step size (4)
    IADS=1, 
    THL={}, 
    THU={}, 
    UZR={}, 
    STOP={}
)
u_sols.to_csv("u_sols_A3b_neg.csv", index=False)
# %% Eigenvalues
"""
eigs_tot = u_sols["eigenvalues"]
eigs_list = eigs_tot[1] # <-- list of eigenvalues
"""
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
        potential_keys += [f'{target_cell}/RPO_bI/v']
    
    sources = u_sols[potential_keys]
    all_pot_cont.append(np.sum(sources, axis=1))

all_pot_cont = np.array(all_pot_cont).T
# %%
pot_cont_df = pd.DataFrame(all_pot_cont, columns=cells)
pot_cont_df.insert(0, cont_param, u_sols[f'{cont_param}'].values)
pot_cont_df['stability'] = u_sols['stability'].values
pot_cont_df['bifurcation'] = u_sols['bifurcation'].values
pot_cont_df.to_csv("pot_cont_df_A3b.csv", index=False)
# %%
stability_idx = pot_cont_df.columns.get_loc('stability')
cell_cols = pot_cont_df.columns[1:stability_idx]  # skip cont_param (index 0), stop before 'stability'

layers = [cell_cols[:4],   # first 4 cells
          cell_cols[4:]    # next 4 cells
         ]

fig, axes = plt.subplots(1, 2, figsize=(14, 10))  
axes = axes.flatten()
for ax, cols in zip(axes, layers):
    labels_with_final = [f"{col} ({pot_cont_df[col].iloc[-1]:.6f})" for col in cols]
    ax.plot(pot_cont_df[cont_param], pot_cont_df[cols])
    ax.legend(labels_with_final, loc="best")
    ax.set_ylabel("[mV]")
    ax.set_xlabel("Hz")
    ax.set_title(", ".join(cols))

fig.suptitle("Potentials", fontsize=20, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()
# %% Plot
"""
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