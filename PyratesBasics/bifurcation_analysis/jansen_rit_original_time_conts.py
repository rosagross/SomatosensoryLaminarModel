# %% Libraries import
from pprint import pprint
from pycobi import ODESystem
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyrates.frontend import CircuitTemplate
import os 
os.chdir("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/bifurcation_analysis/""")

# %% installation directory of AUTO-07p

auto_dir = "/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p"
model = "model_templates.neural_mass_models.jansenrit.JRC"

# structure of the .yaml model
jrc = CircuitTemplate.from_yaml("model_templates.neural_mass_models.jansenrit.JRC")

# %% updating the model (u=220Hz-->0.0Hz)
for p in [0.0,60,120,180]:
    jrc.update_var(node_vars={'pc/rpo_e_in/u': p})


    jrc_auto = ODESystem.from_yaml(model, auto_dir=auto_dir, init_cont = False, NPR=100, NMX=30000)

# V_PCE = pc/rpo_e/v
# V_PCI = pc/rpo_i/v

# Time continuation - PyCobi
    t_sols, t_cont = jrc_auto.run(
        c='ivp', name='time', DS=1e-4, DSMIN=1e-10, EPSL=1e-08, EPSU=1e-08, EPSS=1e-06,
        DSMAX=1e-4, NMX=100000, UZR={14: 2.0}, STOP={'UZ1'})


#  plot 
    v_pce = t_sols["pc/rpo_e_in/v"]
    v_pci = t_sols["pc/rpo_i/v"]
    pc = v_pce + v_pci
    t = t_sols["t"]
    plt.plot(t,v_pce,label='V_PCE')
    plt.plot(t,v_pci,label='V_PCI')
    plt.plot(t,pc,label='PC')
    plt.legend()
    plt.show()
    plt.close()