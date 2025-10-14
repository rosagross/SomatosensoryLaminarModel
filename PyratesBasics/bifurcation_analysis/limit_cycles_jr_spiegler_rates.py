"""Bifurcation analysis of the original Jansen-Rit Model (Spiegler with external rates)"""

# %%
from pprint import pprint
from pycobi import ODESystem
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyrates.frontend import OperatorTemplate, NodeTemplate, CircuitTemplate
from copy import deepcopy
import math
import os 
os.chdir("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/bifurcation_analysis/""")

# Model definition --> Pyrates

# %% Operators

pro = OperatorTemplate(
    name='pro', path=None,
    equations=["m_out = 2.*m_max / (1 + exp(r*(V_thr - v)))"],
    variables={'m_out': 'output',
               'v': 'input',
               'V_thr': 6.0,
               'm_max': 2.5,
               'r': 0.56},
    description="sigmoidal potential-to-rate operator")

rpo_e = OperatorTemplate(
    name='rpo_e', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * m_in - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'm_in': 'input',
               'tau': 0.01,
               'H': 3.25},
    description="excitatory rate-to-potential operator")

rpo_i = deepcopy(rpo_e).update_template(
    name='rpo_i', path=None, variables={'H': -22.0, 'tau': 0.02}
)

rpo_e_in_pc = OperatorTemplate(
    name='rpo_e_in_pc', 
    path=None, 
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * (c_it*m_it) - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'm_it': 'input',
               'tau': 0.01,
               'H': 3.25,
               'c_it': 1.0},
    description="excitatory rate-to-potential operator for the extrinsic input in the PC population"
)

rpo_e_in_ein = deepcopy(rpo_e_in_pc).update_template(
    name='rpo_e_in_ein', path=None, variables={'m_it': 0.0, 'c_it': 1.0}
)

rpo_e_in_iin = deepcopy(rpo_e_in_pc).update_template(
    name='rpo_e_in_iin', path=None, variables={'m_it': 0.0, 'c_it': 1.0} # {'H': -0.022, 'tau': 0.02} to make it inhibitory
)

# %% Populations and Circuit

# Nodes
ein = NodeTemplate(name="ein", path=None, operators=[pro, rpo_e, rpo_e_in_ein])
iin = NodeTemplate(name="iin", path=None, operators=[pro, rpo_e, rpo_e_in_iin])
pc = NodeTemplate(name="pc", path=None, operators=[pro, rpo_e, rpo_i, rpo_e_in_pc])

# Circuit
jrc = CircuitTemplate(
    name="jrc", nodes={'pc': pc, 'ein': ein, 'iin': iin},
    edges=[("pc/pro/m_out", "iin/rpo_e/m_in", None, {'weight': 33.75}),
           ("pc/pro/m_out", "ein/rpo_e/m_in", None, {'weight': 135.}),
           ("ein/pro/m_out", "pc/rpo_e/m_in", None, {'weight': 108.}),
           ("iin/pro/m_out", "pc/rpo_i/m_in", None, {'weight': 33.75})],
    path=None)

# %%

# update input variable to start in a high-activity state
jrc.update_var(node_vars={'pc/rpo_e_in_pc/m_it': 400}) # Hz

# Initial values of the external inputs:
#jrc.update_var(node_vars={'ein/rpo_e_in_ein/m_it': 0})
#jrc.update_var(node_vars={'iin/rpo_e_in_iin/m_it': 0})

# %%
# set state variables close to steady-state value
jrc.update_var(node_vars={'ein/rpo_e/v': 17.787791}) # y0
jrc.update_var(node_vars={'pc/rpo_e/v': 17.526185}) # y1
jrc.update_var(node_vars={'pc/rpo_i/v': -21.927124})
jrc.update_var(node_vars={'iin/rpo_e/v': 4.446948 })

#last_row_df = t_sols.iloc[[-1]]

# %% PyCobi model definition:
jrc_auto = ODESystem.from_template(jrc, auto_dir = "/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p", init_cont=False)
pprint(jrc_auto._var_map)
# %% Time continuation in PyCobi:
#   - `DS` defines the initial step-size of the time continuation (in ms)
#   - `DSMIN` defines the minimal step-size of the time continuation (in ms)
#   - `DSMAX` defines the maximal step-size of the time continuation (in ms)
#   - `NMX` defines the maximum number of continuation steps to perform
#   - `UZR={14: 1000.0}` tells auto-07p to create a user-specified marker when the parameter 14, which is the
#     default parameter field in auto-07p in which time is stored, reaches a value of 1000.0 (ms)
#   - `STOP={'UZ1'}` tells auto-07p to stop the continuation ones it hits the first user-specified marker

t_sols, t_cont = jrc_auto.run(
    c='ivp', name='time', NPR=10, DS=1e-2, DSMIN=1e-7, EPSL=1e-05, EPSU=1e-05, EPSS=1e-03,
    DSMAX=5e-1, NMX=100000, UZR={14: 2.0}, STOP={'UZ1'})

# %% Outputs of the time continuation:
#   - t_sols: pandas Dataframe with the summary of the results of the simulation
#   - t_cont: type 'branch' or 'bifDiag', an auto-07p object that can be used for subsequent parameter continuations

# %% plotting the solutions:
v_pce = t_sols["pc/rpo_e/v"]
v_pci = t_sols["pc/rpo_i/v"]
v_pc_ext = t_sols["pc/rpo_e_in_pc/v"]
pc = v_pce + v_pci + v_pc_ext
t = t_sols["t"]
plt.plot(t,v_pce,label='V_PCE')
plt.plot(t,v_pci,label='V_PCI')
plt.plot(t,v_pc_ext,label='V_PC_EXT')
plt.plot(t,pc,label='PC')
plt.legend()
plt.show()

# %% BIFURCATION ANALYSIS
# 1D parameter continuation in the input parameter 'u':

mT_sols, mT_cont = jrc_auto.run(
    origin=t_cont, 
    starting_point='UZ1', 
    name='m_T', # Name of continuation parameter (variable to vary) 
    bidirectional=True, #  Continue both forward and backward from starting point
    ICP='pc/rpo_e_in_pc/m_it',
    RL0=-50, RL1=600, # Lower and upper bounds for continuation parameter(s)
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
    NMX=100000, 
    NPR=10, # Print/report frequency (every NPR steps)
    MXBF={}, # Maximum number of bifurcations to locate
    IID=2,
    ITMX=1000, # max.num. of iterations allowed in the accurate location of special solutions
    ITNW=40, 
    NWTN=12, 
    JAC=0, 
    EPSL=1e-07, # recommended in the auto docs
    EPSU=1e-07, # recommended in the auto docs
    EPSS=1e-05, # recommended in the auto docs (approx.100/1000 times EPSL,EPSU)
    DS=1e-04, # initial continuation step size 
    DSMIN=1e-10, # minimum allowed step size
    DSMAX=0.1, # maximum allowed step size 
    IADS=1, 
    THL={}, 
    THU={}, 
    UZR={}, 
    STOP={}
)

jrc_auto.plot_continuation('pc/rpo_e_in_pc/m_it', 'pc/rpo_e/v', cont='m_T')
plt.show()

# %% plotting the y with the PyCobi function (=y1-y2)

y1 = mT_sols["pc/rpo_e/v"].iloc[:, 0]
y2 = mT_sols["pc/rpo_i/v"].iloc[:, 0]
v_ext = mT_sols["pc/rpo_e_in_pc/v"].iloc[:, 0]
y = y1+y2+v_ext
mT_sols[('y', 0)] = y # adding the y as a variable in the dataframe with the continuation
jrc_auto.plot_continuation('pc/rpo_e_in_pc/m_it', 'y', cont='m_T')
plt.show()

# %%
C = 135
y0 = mT_sols["ein/rpo_e/v"]/C
mT_sols["y0"] = y0 # adding the y as a variable in the dataframe with the continuation
jrc_auto.plot_continuation('pc/rpo_e_in_pc/m_it', 'y0', cont='m_T')
plt.show()

# %%
"""
u1_sols, u1_cont = jrc_auto.run(
    origin=u_cont, 
    starting_point='HB1', 
    name='u_1', # Name of continuation parameter (variable to vary) 
    bidirectional=True, #  Continue both forward and backward from starting point
    ICP='pc/rpo_e_in/u',
    RL0=-50.0, RL1=600.0, # Lower and upper bounds for continuation parameter(s)
    IPS=2, # Problem type: 1=Equilibrium, 2=Periodic orbit continuation
    ILP=1, # Detect limit points (folds): 1=Yes, 0=No 
    ISP=2, # Bifurcation detection level: 0=none, 1=some, 2=more 
    ISW=-1, 
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
# %%
fig, ax = plt.subplots()
y1_d = u1_sols["pc/rpo_e_in/v"].iloc[:, 0]
y2_d = u1_sols["pc/rpo_i/v"].iloc[:, 0]
y_d = y1_d+y2_d
u1_sols["y_d"] = y_d # adding the y as a variable in the dataframe with the continuation
jrc_auto.plot_continuation('pc/rpo_e_in/u', 'y_d', cont='u_1',ax=ax)

y1_u = u1_sols["pc/rpo_e_in/v"].iloc[:, 1]
y2_u = u1_sols["pc/rpo_i/v"].iloc[:, 1]
y_u = y1_u+y2_u
u1_sols["y_u"] = y_u # adding the y as a variable in the dataframe with the continuation
jrc_auto.plot_continuation('pc/rpo_e_in/u', 'y_u', cont='u_1',ax=ax)
plt.show()

# %%
# PLOTTHE VOLTAGE!!!
fig, ax = plt.subplots()
ax.plot(u_sols['pc/rpo_e_in/u'], y)  
ax.plot(u1_sols['pc/rpo_e_in/u'], y_d)
ax.plot(u1_sols['pc/rpo_e_in/u'], y_u)    
plt.show()
"""

# %% Continuation of the limit cycles - updated
limit_cycles_sols = {}
limit_cycles_conts = {}
for i in range(1,4):
    limit_cycles_sols[f'mT_lc{i}_sols'],limit_cycles_conts[f'mT_lc{i}_conts'] = jrc_auto.run(origin=mT_cont, starting_point=f'HB{i}', name=f'mT_lc{i}', IPS=2, ISP=2, ISW=-1, NPR=20, STOP=['BP2'], get_period = True)

# %% 
"""
fig, ax = plt.subplots()
jrc_auto.plot_continuation('pc/pro_ext/v_ext', 'pc/rpo_e/v', cont='v_T',ax=ax)
jrc_auto.plot_continuation('pc/pro_ext/v_ext', 'pc/rpo_e/v', cont=f'vT_lc1', ax=ax, ignore=["UZ", "BP"],
                               line_color_stable="green")
plt.show()
# %%
fig, ax = plt.subplots()
jrc_auto.plot_continuation('pc/pro_ext/v_ext', 'pc/rpo_i/v', cont='v_T',ax=ax)
jrc_auto.plot_continuation('pc/pro_ext/v_ext', 'pc/rpo_i/v', cont=f'vT_lc1', ax=ax, ignore=["UZ", "BP"],
                               line_color_stable="green")
plt.show()
"""
# %% Plotting of the limit cycle continuations - y0
fig, ax = plt.subplots()
jrc_auto.plot_continuation('pc/rpo_e_in_pc/m_it', 'y0', cont='m_T',ax=ax)

for i in range(1,4):
    # lower branch
    y0_d = limit_cycles_sols[f'mT_lc{i}_sols']["ein/rpo_e/v"].iloc[:, 0]/C
    # upper branch
    y0_u = limit_cycles_sols[f'mT_lc{i}_sols']["ein/rpo_e/v"].iloc[:, 1]/C
    
    limit_cycles_sols[f'mT_lc{i}_sols'][('y0', 0)] = y0_d
    limit_cycles_sols[f'mT_lc{i}_sols'][('y0', 1)] = y0_u
    jrc_auto.plot_continuation('pc/rpo_e_in_pc/m_it', 'y0', cont=f'mT_lc{i}', ax=ax, ignore=["UZ", "BP"],
                               line_color_stable="green")
#ax.set_ylim([0.0, 0.012]) 
#ax.set_xlim([-50, 150])
plt.show()

# %% Plotting of the limit cycle continuations - y
# parameters:
A = 3.25
a = 100
B = 22
b = 50
C = 135
e0 = 2.5
r = 0.56
v0 = 6
# %%
fig, ax = plt.subplots()
jrc_auto.plot_continuation('pc/rpo_e_in_pc/m_it', 'y', cont='m_T',ax=ax)

for i in range(1,4):
   
    # lower branch
    y0_d = limit_cycles_sols[f'mT_lc{i}_sols']["ein/rpo_e/v"].iloc[:, 0]/C
    # upper branch
    y0_u = limit_cycles_sols[f'mT_lc{i}_sols']["ein/rpo_e/v"].iloc[:, 1]/C

    y_d = np.real(v0 - (1/r) * np.log(((2*A*e0)/(a*y0_d.astype(complex))) - 1))
    y_u = np.real(v0 - (1/r) * np.log(((2*A*e0)/(a*y0_u.astype(complex))) - 1))

    limit_cycles_sols[f'mT_lc{i}_sols'][('yy', 0)] = y_d 
    limit_cycles_sols[f'mT_lc{i}_sols'][('yy', 1)] = y_u 
    jrc_auto.plot_continuation('pc/rpo_e_in_pc/m_it', 'yy', cont=f'mT_lc{i}', ax=ax, ignore=["UZ", "BP"],
                               line_color_stable="green")
#ax.set_ylim([0.0, 12]) 
#ax.set_xlim([-2, 12])
plt.show()

# %%
