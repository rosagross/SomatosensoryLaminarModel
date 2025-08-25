"""Bifurcation analysis of the original Jansen-Rit Model (.yaml pre-implemented in Pyrates)"""

# %% Libraries import:
from pycobi import ODESystem
import matplotlib.pyplot as plt
from pyrates.frontend import CircuitTemplate
import os 
os.chdir("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/bifurcation_analysis/""")

# %% Initialize the Jansen-Rit model in PyRates:
jrc = CircuitTemplate.from_yaml("model_templates.neural_mass_models.jansenrit.JRC")

# update input variable to start in a high-activity state
jrc.update_var(node_vars={'pc/rpo_e_in/u': 400.0})

# set state variables close to steady-state value
jrc.update_var(node_vars={'ein/rpo_e/v': 1.762243e-02 }) # y0
jrc.update_var(node_vars={'pc/rpo_e_in/v': 3.052549e-02}) # y1
jrc.update_var(node_vars={'pc/rpo_i/v': -2.195986e-02})
jrc.update_var(node_vars={'iin/rpo_e/v': 4.405607e-03 })

# %% PyCobi model definition:
jrc_auto = ODESystem.from_template(jrc, auto_dir = "/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p", init_cont=False)

# %% Time continuation in PyCobi:
#   - `DS` defines the initial step-size of the time continuation (in ms)
#   - `DSMIN` defines the minimal step-size of the time continuation (in ms)
#   - `DSMAX` defines the maximal step-size of the time continuation (in ms)
#   - `NMX` defines the maximum number of continuation steps to perform
#   - `UZR={14: 1000.0}` tells auto-07p to create a user-specified marker when the parameter 14, which is the
#     default parameter field in auto-07p in which time is stored, reaches a value of 1000.0 (ms)
#   - `STOP={'UZ1'}` tells auto-07p to stop the continuation ones it hits the first user-specified marker

t_sols, t_cont = jrc_auto.run(
    c='ivp', name='time', DS=1e-4, DSMIN=1e-10, EPSL=1e-08, EPSU=1e-08, EPSS=1e-06,
    DSMAX=1e-3, NMX=50000, UZR={14: 2.0}, STOP={'UZ1'})

# %% Outputs of the time continuation:
#   - t_sols: pandas Dataframe with the summary of the results of the simulation
#   - t_cont: type 'branch' or 'bifDiag', an auto-07p object that can be used for subsequent parameter continuations

# %% plotting the solutions:
v_pce = t_sols["pc/rpo_e_in/v"]
v_pci = t_sols["pc/rpo_i/v"]
pc = v_pce + v_pci
t = t_sols["t"]
plt.plot(t,v_pce,label='V_PCE')
plt.plot(t,v_pci,label='V_PCI')
plt.plot(t,pc,label='PC')
plt.legend()
plt.show()

# %% BIFURCATION ANALYSIS
# 1D parameter continuation in the input parameter 'u':

u_sols, u_cont = jrc_auto.run(
    origin=t_cont, 
    starting_point='UZ1', 
    name='u', # Name of continuation parameter (variable to vary) 
    bidirectional=True, #  Continue both forward and backward from starting point
    ICP='pc/rpo_e_in/u',
    RL0=-50.0, RL1=600.0, # Lower and upper bounds for continuation parameter(s)
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
    NMX=8000, 
    NPR=20, # Print/report frequency (every NPR steps)
    MXBF={}, # Maximum number of bifurcations to locate
    IID=2,
    ITMX=1000, # max.num. of iterations allowed in the accurate location of special solutions
    ITNW=40, 
    NWTN=12, 
    JAC=0, 
    EPSL=1e-07, # recommended in the auto docs
    EPSU=1e-07, # recommended in the auto docs
    EPSS=1e-05, # recommended in the auto docs (approx.100/1000 times EPSL,EPSU)
    DS=1e-03, # initial continuation step size (0.5)
    DSMIN=1e-8, # minimum allowed step size
    DSMAX=0.5, # maximum allowed step size (4)
    IADS=1, 
    THL={}, 
    THU={}, 
    UZR={}, 
    STOP={}
)
jrc_auto.plot_continuation('pc/rpo_e_in/u', 'pc/rpo_e_in/v', cont='u')
plt.show()

# %% Continuation of the limit cycles
for i in range(1,4):
    jrc_auto.run(origin=u_cont, starting_point=f'HB{i}', name=f'u_lc{i}', IPS=2, ISP=2, ISW=-1, STOP=['BP2'])

# %% Plotting of the limit cycle continuations
fig, ax = plt.subplots()
jrc_auto.plot_continuation('pc/rpo_e_in/u', 'pc/rpo_e_in/v', cont='u', ax=ax)
for i in range(1,4):
    jrc_auto.plot_continuation('pc/rpo_e_in/u', 'pc/rpo_e_in/v', cont=f'u_lc{i}', ax=ax, ignore=["UZ", "BP"],
                               line_color_stable="green")
plt.show()

# %%
# This concludes our use example. As a final step, we clear all the temporary files we created and make sure all
# working directory changes that have happened in the background during file creation and parameter continuation are
# undone. This can be done via a single call:

jrc.clear()
jrc_auto.close_session(clear_files=True)
