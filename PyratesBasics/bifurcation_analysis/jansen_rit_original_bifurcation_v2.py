"""Bifurcation analysis of the original Jansen-Rit Model (.yaml pre-implemented in Pyrates)
CONTINUATION.PY VERSION (PyRates docs)"""

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

jrc.update_var(node_vars={'pc/rpo_e_in/u': 0.0})

# %%

jrc_auto = ODESystem.from_yaml(model, auto_dir=auto_dir, init_cont = False, NPR=100, NMX=30000)
# - the option "init_cont = False" does not perform a time integration when defining the system
pprint(jrc_auto._var_map)

# V_PCE = pc/rpo_e/v
# V_PCI = pc/rpo_i/v

# %%
"""For generating a parameter continuation, the initial point should be an equilibrium point or a periodic
orbit (= the solution to the initial value problem would be constant or periodic in time.)
"""
# %% Time continuation - first way: PyRates
# This can be simply achieved by a call to the
# :code:`CircuitTemplate.run()` method, before calling the :code:`get_run_func()` method.
# Then, PyRates will automatically use the values of the state variables from the last simulation step.
"""
results = jrc.run(simulation_time=2.0,
                  step_size=1e-4,
                  sampling_step_size=1e-3,
                  outputs={'V_PCE': 'PC/RPO_e/v',
                           'V_PCI': 'PC/RPO_i/v'},
                  backend='default',
                  solver='scipy')

import matplotlib.pyplot as plt
plt.plot(results['V_PCE']+results['V_PCI'], label='PC Membrane potential')
for V in results:
    plt.plot(results[V],label=V)
    
plt.legend()
plt.show()"""

# %% Time continuation - second way: PyCobi
t_sols, t_cont = jrc_auto.run(
    c='ivp', name='time', DS=1e-4, DSMIN=1e-10, EPSL=1e-08, EPSU=1e-08, EPSS=1e-06,
    DSMAX=1e-4, NMX=100000, UZR={14: 2.0}, STOP={'UZ1'})

#   - `DS=1e-3` defines the initial step-size of the time continuation (in ms)
#   - `DSMIN=1e-4` defines the minimal step-size of the time continuation (in ms)
#   - `DSMAX=1.0` defines the maximal step-size of the time continuation (in ms)
#   - `NMX=10000` defines the maximum number of continuation steps to perform
#   - `UZR={14: 1000.0}` tells auto-07p to create a user-specified marker when the parameter 14, which is the
#     default parameter field in auto-07p in which time is stored, reaches a value of 1000.0 (ms)
#   - `STOP={'UZ1'}` tells auto-07p to stop the continuation ones it hits the first user-specified marker

# Outputs:
# t_sols: dict with the summary of the results of the simulation
# t_cont: type 'branch' or 'bifDiag', an auto-07p object that is required for subsequent parameter continuations



# %% plot 
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
# alternative for plotting with the PyCobi functions
"""
jrc_auto.plot_continuation('PAR(14)', 'U(1)', cont='time')
jrc_auto.plot_continuation('PAR(14)', 'U(3)', cont='time')
plt.show()"""


# %% BIFURCATION ANALYSIS
# 'pc/rpo_e_in/u': {'cont': 1, 'plot': 'PAR(1)'}

"""c='ivp', name='time', DS=1e-4, DSMIN=1e-10, EPSL=1e-08, EPSU=1e-08, EPSS=1e-06,
    DSMAX=1e-4, NMX=100000, UZR={14: 2.0}, STOP={'UZ1'}"""

u_sols, u_cont = jrc_auto.run(
    origin=t_cont, 
    starting_point='EP', 
    name='u', # Name of continuation parameter (variable to vary) 
    bidirectional=True, #  Continue both forward and backward from starting point
    ICP=1, # Index of parameter(s) to continue in (AUTO parameter numbering) --> checl from the ODESystem!
    RL0=-50.0, RL1=700.0, # Lower and upper bounds for continuation parameter(s)
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
    NMX=999999999, 
    NPR=1, # Print/report frequency (every NPR steps)
    MXBF={}, # Maximum number of bifurcations to locate
    IID=2,
    ITMX=100, # max.num. of iterations allowed in the accurate location of special solutions
    ITNW=40, 
    NWTN=12, 
    JAC=0, 
    EPSL=1e-06, # recommended in the auto docs
    EPSU=1e-06, # recommended in the auto docs
    EPSS=1e-04, # recommended in the auto docs (approx.100/1000 times EPSL,EPSU)
    DS=0.5, # initial continuation step size (0.5)
    DSMIN=1e-8, # minimum allowed step size
    DSMAX=4, # maximum allowed step size (4)
    IADS=1, 
    THL={}, 
    THU={}, 
    UZR={"pc/rpo_e_in/u": 700}, 
    STOP={}
)
jrc_auto.plot_continuation('PAR(1)', 'U(1)', cont='u')
plt.show()
# %% 
# the variables can also be referred with their names in the model!
# - ode.results: gives the results of the simulation as a dict
# - ode._var_map: gives the relationship between the names of the ODESystem and the model (even if "params" is not passed when creating the system)
# - ode.get_summary(): pandas dataframe with the summaries of a continuation
# U(i): state variables
# PAR(14): time 

# plot of the steady-state solution

jrc_auto.plot_continuation('PAR(1)', 'U(5)', cont='u')
plt.show()

# addutionally, the timeseries of a variable or a trajectory can be plotted (respectively with .plot_timeseries or .plot_trajectory)
# look up how to refer direcly to the variables
"""ode.plot_continuation("t", "p/qif_sfa_op/r", cont=0)
plt.show()"""
# %% state variables and plots
v0 = 6 # mV
r = 0.560 # mV^-1
e0 = 2.5 # Hz
A = 3.25 # mV
a = 100 # Hz
tau_e = 1/a

y1 = u_sols["pc/rpo_e_in/v"].to_numpy()*1000
y2 = -u_sols["pc/rpo_i/v"].to_numpy()*1000
y = y1-y2
y0 = (A/a)*((2*e0)/(1+np.exp(r*(v0-y))))
p = u_sols["pc/rpo_e_in/u"]

# %% Plots

fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(p, y0)
axs[0, 0].set_title('y_0')
axs[0, 1].plot(p, y1, 'tab:orange')
axs[0, 1].set_title('y_1')
axs[1, 0].plot(p, y2, 'tab:green')
axs[1, 0].set_title('y_2')
axs[1, 1].plot(p, y, 'tab:red')
axs[1, 1].set_title('y')

for ax in axs.flat:
    ax.set(xlabel='p [Hz]', ylabel='[mV]')

plt.tight_layout()
plt.show()
# %%
jrc.clear()

jrc_auto.close_session(clear_files=True)
# %%
