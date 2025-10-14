"""Bifurcation analysis of the original Jansen-Rit Model (.yaml pre-implemented in Pyrates)"""

# %% Libraries import:
from pprint import pprint
from pycobi import ODESystem
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyrates.frontend import CircuitTemplate
import os 
os.chdir("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/bifurcation_analysis/""")

# %% PyRates model:

# installation directory of AUTO-07p
auto_dir = "/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p"
model = "model_templates.neural_mass_models.jansenrit.JRC"

# path for the .yaml Jansen-Rit
# /data/u_mecozzi_software/miniforge3/envs/pyrates_project2/lib/python3.13/site-packages/model_templates/neural_mass_models/

jrc = CircuitTemplate.from_yaml("model_templates.neural_mass_models.jansenrit.JRC")
# updating the model (u=220Hz-->0.0Hz)
jrc.update_var(node_vars={'pc/rpo_e_in/u': 400.0})

# %% changing the initial conditions --> influence on the bifurcation ??
"""
# y0 => ein/rpo_e/v / 135
# y1 => pc/rpo_e_in/v [mV]
# y2 => - pc/rpo_i/v [mV]

y0_init = 0.0 # [mV]
y1_init = 0.0 # [mV]
y2_init = 0.0 # [mV]
C = 135

# updating the model
jrc.update_var(node_vars={'ein/rpo_e/v': (y0_init/C) }) # y0
jrc.update_var(node_vars={'pc/rpo_e_in/v': y1_init}) # y1
jrc.update_var(node_vars={'pc/rpo_i/v': -y2_init}) # y2"""

# %% 
jrc.update_var(node_vars={'ein/rpo_e/v': 1.762243e-02 }) # y0
jrc.update_var(node_vars={'pc/rpo_e_in/v': 3.052549e-02}) # y1
jrc.update_var(node_vars={'pc/rpo_i/v': -2.195986e-02})
jrc.update_var(node_vars={'iin/rpo_e/v': 4.405607e-03 })

# %% PyCobi model definition:

jrc_auto = ODESystem.from_yaml(model, auto_dir=auto_dir, init_cont = False, NPR=100, NMX=30000)
# - the option "init_cont = False" does not perform a time integration when defining the system
pprint(jrc_auto._var_map)

# %% Time continuation - PyCobi:
t_sols, t_cont = jrc_auto.run(
    c='ivp', name='time', DS=1e-4, DSMIN=1e-10, EPSL=1e-08, EPSU=1e-08, EPSS=1e-06,
    DSMAX=1e-3, NMX=50000, UZR={14: 2.0}, STOP={'UZ1'})

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

# %% plot of the solutions:

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

# %% BIFURCATION ANALYSIS
# 'pc/rpo_e_in/u': {'cont': 1, 'plot': 'PAR(1)'}

"""c='ivp', name='time', DS=1e-4, DSMIN=1e-10, EPSL=1e-08, EPSU=1e-08, EPSS=1e-06,
    DSMAX=1e-4, NMX=100000, UZR={14: 2.0}, STOP={'UZ1'}""" 

u_sols, u_cont = jrc_auto.run(
    origin=t_cont, 
    starting_point='UZ1', 
    name='u', # Name of continuation parameter (variable to vary) 
    bidirectional=True, #  Continue both forward and backward from starting point
    ICP=1, # Index of parameter(s) to continue in (AUTO parameter numbering) --> checl from the ODESystem!
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
jrc_auto.plot_continuation('PAR(1)', 'U(1)', cont='u')
plt.show()

# %% Analysis of the state variables + plots
"""(notation related to Grimbert, F., Faugeras, O., 2006. Bifurcation analysis of Jansen's neural mass model.
Neural Comput. 18, 3052–3068.)

- state variables:

# y1 = U(1)
# y2 = -U(3)
"""
# Parameters:--> mV to make computation easier
v0 = 6 # mV
r = 0.560 # mV^-1
e0 = 2.5 # Hz
A = 3.25 # mV
a = 100 # Hz
tau_e = 1/a

y1 = u_sols["pc/rpo_e_in/v"].to_numpy()*1000
y2 = -u_sols["pc/rpo_i/v"].to_numpy()*1000
y = y1-y2
#y0 = (A/a)*((2*e0)/(1+np.exp(r*(v0-y))))
y0 = u_sols["ein/rpo_e/v"].to_numpy()*1000/C
p = u_sols["pc/rpo_e_in/u"]

fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(p, y0)
axs[0, 0].set_title('y_0')
axs[0, 1].plot(p, y1, 'tab:orange')
axs[0, 1].set_title('y_1')
axs[1, 0].plot(p, y2, 'tab:green')
axs[1, 0].set_title('y_2')
axs[1, 1].plot(p, y, 'tab:red')
axs[1, 1].set_title('y=y_1-y_2')

for ax in axs.flat:
    ax.set(xlabel='p [Hz]', ylabel='[mV]')

plt.tight_layout()
plt.show()
# %%
hopf_sols, hopf_cont = jrc_auto.run(
    origin=u_cont, starting_point='HB2', name='u_hopf',
    IPS=2, ISP=2, ISW=-1, STOP=['BP2']
)

# %%
fig, ax = plt.subplots()
jrc_auto.plot_continuation('PAR(1)', 'U(4)', cont='u', ax=ax)
jrc_auto.plot_continuation('PAR(1)', 'U(4)', cont='u_hopf',ax=ax, ignore=["UZ", "BP"])

plt.show()

# %%
# time continuations in a specified point to see the behav
"""
ode.run(origin=hopf_cont, starting_point="UZ1", c="ivp", name="lc")
ode.plot_continuation("t", "p/qif_sfa_op/r", cont="lc")
plt.show()"""


jrc.clear()
jrc_auto.close_session(clear_files=True)
# %%
