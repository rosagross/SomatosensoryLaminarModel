# %% Main
"""
Main code to test the model class
"""
# %% libraries import
import os 
WDDIR = os.getenv("WDDIR")
SIMDIR = os.getenv("SIMDIR")
#os.chdir(os.path.join(WDDIR,"PyratesBasics","exp_model")) 
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
from pycobi import ODESystem
#from parameters import Parameter
param_path = os.path.join(WDDIR,"Simulations")
if param_path not in sys.path:
    sys.path.append(param_path)
from parameters import Parameter
import json
# %%
#from somato_model_pyrates_no_conn_operators_complete_pycobi import SomatoModelPyrates
from somato_model_pyrates_sI_gstrength_operators_complete_pycobi import SomatoModelPyrates

# %%
def read_simulation_params():
    WDDIR = os.getenv("WDDIR")
    with open(os.path.join(WDDIR, 'Simulations', 'simulation_parameter.json'), 'r') as json_file:
        params = json.load(json_file)
    return params

# %% creation of the equilibrium .csv
"""
modello_prova.simulate()
results = modello_prova.simulation_results
equilibrium_df = results.iloc[[-1]]
N_cells = len(modello_prova.cells)
outputs = {}
b_inputs = []
for i, target_cell in enumerate(modello_prova.cells):
    if i == (N_cells-2):
        for rpo_name_ext in modello_prova.rpo_names_extended[:N_cells+1]:
                    outputs[f'V_{target_cell}/{rpo_name_ext}'] = f'{target_cell}/{rpo_name_ext}/v'
    else:
        for rpo_name in modello_prova.rpo_names[:N_cells]:
            outputs[f'V_{target_cell}/{rpo_name}'] = f'{target_cell}/{rpo_name}/v'

    if i in range(N_cells-2):
                #print(f'{target_cell}/RPO_bI/v')
        outputs[f'V_{target_cell}/RPO_bI'] = f'{target_cell}/RPO_bI/v'
        #outputs['V_background/RPO_bI'] = 'BACKGROUND/RPO_bI/v_bI'
        
equilibrium_df.rename(columns=outputs, inplace=True)
equilibrium_df.to_csv("equilibrium_df_complete.csv", index=False)
"""
# %%
# equilibrium csv: "equilibrium_df_complete.csv"
output_dir = os.path.join(WDDIR, 'PyratesBasics', 'exp_model','complete_model_continuations')
os.makedirs(output_dir, exist_ok=True) 
params = read_simulation_params()
range_par = [0.0, 4.0]
modello_prova = SomatoModelPyrates(params)
cont_param = 'G/g_definition/g_input'
auto_dir_path = "/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p"
modello_prova.pycobi_continuation(cont_param,range_par,auto_dir_path,"/data/p_02989/Modelling/mecozzi_wd/SomatosensoryLaminarModel/PyratesBasics/exp_model/equilibrium_df_complete.csv")
cont_df = modello_prova.continuation_df(cont_param)
# csv saving
filename = "complete_model_bifurcation_sI026.csv"
filepath = os.path.join(output_dir, filename)
cont_df.to_csv(filepath, index=False)
print("Simulation finished")
# %%