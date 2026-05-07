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
from pyrates import clear_frontend_caches

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
from somato_model_pyrates_no_conn_operators_A3bS1_pycobi import SomatoModelPyrates
# %%
def read_simulation_params():
    WDDIR = os.getenv("WDDIR")
    with open(os.path.join(WDDIR, 'Simulations', 'simulation_parameter.json'), 'r') as json_file:
        params = json.load(json_file)
    return params

# %%
params = read_simulation_params() # --> "default" parameters
# %%
bEI_list = np.arange(0.70, 1, 0.01) # [0.5,0.99] with 0.01 step
coupling_strength_list = np.arange(10, 100, 10) # [10,50] with 10 step
bI_list = np.arange(5, 10, 0.1)
range_par = [0.0, 50.0]
cont_param = 'background_input/bI_definition/backgrndI_strengths'
auto_dir_path = "/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p"
# %% output folder creation
output_dir = os.path.join(WDDIR, 'PyratesBasics', 'exp_model','bi_prova3')
os.makedirs(output_dir, exist_ok=True) 
# %% PyRates simulations 

for k in range(len(bI_list)):
    clear_frontend_caches()
    params_new = deepcopy(params)  # instead of params.copy()
    params_new["Ib_strength"] = float(bI_list[k])
    #for i in range(len(bEI_list)):
    #    params["balance_EI"] = float(bEI_list[i])
    #    for j in range(len(coupling_strength_list)):
    #        params["coupling_strength"] = float(coupling_strength_list[j])
    modello_prova = SomatoModelPyrates(params_new)
    modello_prova.simulate()
    potential_df = modello_prova.potential_df
    filename = f"bI_{bI_list[k]:.2f}"
    filepath = os.path.join(output_dir, filename)
    potential_df.to_csv(filepath, index=False)
    #modello_prova.pyrates_plot(filename,filepath)
    #plt.close()
            #fig.suptitle(filename, fontsize=20, fontweight='bold')
    del modello_prova 
    del params_new
    del potential_df

# %% plotting 

for k in range(len(bI_list)):
    filename = f"bI_{bI_list[k]:.2f}"
    filepath = os.path.join(output_dir, filename)
    df_temp = pd.read_csv(filepath)
        #param_vals = df_temp[f'{cont_param}']
    #df_temp = df_temp.drop(cont_param, axis=1)
    layers = [df_temp.columns[:4],   # first 4 
                df_temp.columns[4:]  # next 4
                ]
    fig, axes = plt.subplots(1, 2, figsize=(14, 10))  
    axes = axes.flatten()
    for ax, cols in zip(axes, layers):
            # Create labels with final values
            #labels_with_final = [f"{col} ({df_temp[col].iloc[-1]:.6f})" for col in cols]
            
            # Plot and override legend labels
        ax.plot(df_temp[cols])  # suppress default legend
        #ax.legend(labels_with_final, loc="best")
        ax.set_ylabel("[mV]")
        ax.set_xlabel("Hz")
        ax.set_title(", ".join(cols))
        #axes[-1].set_visible(False)
    fig.suptitle(filename, fontsize=20, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

# %% 
# definition of the model with different parameters
"""
for i in range(len(bEI_list)):
    params["balance_EI"] = bEI_list[i]
    for j in range(len(coupling_strength_list)):
        params["coupling_strength"] = coupling_strength_list[j]
        modello_prova = SomatoModelPyrates(params)
        modello_prova.pycobi_continuation(cont_param,range_par,auto_dir_path)
        cont_df = modello_prova.continuation_df(cont_param)
        # csv saving
        filename = f"cont_balanceEI_{bEI_list[i]:.2f}_cplstr_{coupling_strength_list[j]:.2f}.csv"
        filepath = os.path.join(output_dir, filename)
        cont_df.to_csv(filepath, index=False)
        #modello_prova.pycobi_plot(cont_param)
        print("Simulation finished")
        modello_prova.model_auto.close_session(clear_files=True)
        del modello_prova 
"""
# %% PyRates simulation
"""
modello_prova = SomatoModelPyrates(params)
modello_prova.simulate()
potential_df = modello_prova.potential_df
modello_prova.pyrates_plot()
"""
# %% PyCobi continuation
range_par = [-10.0, 10.0]
modello_prova = SomatoModelPyrates(params)
cont_param = 'background_input/bI_definition/backgrndI_strengths'
auto_dir_path = "/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p"
modello_prova.pycobi_continuation(cont_param,range_par,auto_dir_path,"equilibrium_df_A3b.csv")
# %%
modello_prova.u_sols.to_csv("u_sols_A3b_class.csv", index=False)
#modello_prova.pycobi_plot(cont_param)
cont_df = modello_prova.continuation_df(cont_param)
        # csv saving
filename = "continuation_prova.csv"
filepath = os.path.join(output_dir, filename)
cont_df.to_csv(filepath, index=False)
# %% plotting the continuations
"""
for i in range(len(bEI_list)):
    for j in range(len(coupling_strength_list)):
        filename = f"cont_balanceEI_{bEI_list[i]:.2f}_cplstr_{coupling_strength_list[j]:.2f}.csv"
        filepath = os.path.join(output_dir, filename)
        df_temp = pd.read_csv(filepath)
        #param_vals = df_temp[f'{cont_param}']
        df_temp = df_temp.drop(cont_param, axis=1)
        layers = [df_temp.columns[:4],   # first 4 
                df_temp.columns[4:7]  # next 4
                ]
        fig, axes = plt.subplots(1, 2, figsize=(14, 10))  
        axes = axes.flatten()
        for ax, cols in zip(axes, layers):
            # Create labels with final values
            #labels_with_final = [f"{col} ({df_temp[col].iloc[-1]:.6f})" for col in cols]
            
            # Plot and override legend labels
            ax.plot(df_temp[cols])  # suppress default legend
            ax.legend(labels_with_final, loc="best")
            ax.set_ylabel("[mV]")
            ax.set_xlabel("Hz")
            ax.set_title(", ".join(cols))
        #axes[-1].set_visible(False)
        fig.suptitle(filename, fontsize=20, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()
"""
# %%
