"""Jansen-Rit definition through Pyrates (python) and bifurcation analysis using PyCobi"""
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

# %% Model definition --> Pyrates

pro = OperatorTemplate(
    name='PRO', path=None,
    equations=["m_out = 2.*m_max / (1 + exp(r*(V_thr - v)))"],
    variables={'m_out': 'output',
               'v': 'input',
               'V_thr': 6e-3,
               'm_max': 2.5,
               'r': 560.0},
    description="sigmoidal potential-to-rate operator")

rpo_e = OperatorTemplate(
    name='RPO_e', path=None,
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * m_in - 2 * i/tau - v/tau^2'],
    variables={'v': 'output',
               'i': 'variable',
               'm_in': 'input',
               'tau': 0.01,
               'H': 0.00325},
    description="excitatory rate-to-potential operator")

rpo_e_in = deepcopy(rpo_e).update_template(
    name='RPO_e_in', 
    path=None, 
    equations=['d/dt * v = i',
               'd/dt * i = H/tau * (m_in+u) - 2 * i/tau - v/tau^2'], 
    variables={'u': 0.0},
    description="excitatory rate-to-potential operator with extrinsic input"
)

rpo_i = deepcopy(rpo_e).update_template(
    name='RPO_i', path=None, variables={'H': -0.022, 'tau': 0.02}
)

ein = NodeTemplate(name="EIN", path=None, operators=[pro, rpo_e])
iin = NodeTemplate(name="IIN", path=None, operators=[pro, rpo_e])
pc = NodeTemplate(name="PC", path=None, operators=[pro, rpo_e_in, rpo_i])

jrc = CircuitTemplate(
    name="JRC", nodes={'PC': pc, 'EIN': ein, 'IIN': iin},
    edges=[("PC/PRO/m_out", "IIN/RPO_e/m_in", None, {'weight': 33.75}),
           ("PC/PRO/m_out", "EIN/RPO_e/m_in", None, {'weight': 135.}),
           #("EIN/PRO/m_out", "PC/RPO_e/m_in", None, {'weight': 108.}),
           ("EIN/PRO/m_out", "PC/RPO_e_in/m_in", None, {'weight': 108.}),
           ("IIN/PRO/m_out", "PC/RPO_i/m_in", None, {'weight': 33.75})],
    path=None)


# %% Definition of the system with PyCobi

auto_dir = "/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p"
jrc_auto = ODESystem.from_template(jrc, auto_dir=auto_dir, init_cont = False, NPR=100, NMX=30000)
# - the option "init_cont = False" does not perform a time integration when defining the system
pprint(jrc_auto._var_map)

# %% Time continuation

t_sols, t_cont = jrc_auto.run(
    c='ivp', name='time', DS=1e-4, DSMIN=1e-10, EPSL=1e-08, EPSU=1e-08, EPSS=1e-06,
    DSMAX=1e-4, NMX=100000, UZR={14: 2.0}, STOP={'UZ1'})

# %% plot of the solution

v_pce = t_sols["PC/RPO_e_in/v"]
v_pci = t_sols["PC/RPO_i/v"]
pc = v_pce + v_pci
t = t_sols["t"]
plt.plot(t,v_pce,label='V_PCE')
plt.plot(t,v_pci,label='V_PCI')
plt.plot(t,pc,label='PC')
plt.ylabel('[V]')
plt.xlabel('[s]')
plt.legend()
plt.show()

# %% Bifurcation analysis on the parameter "u" (starting from the initial equilibrium point [0,0,0,0,0,0])

u_sols, u_cont = jrc_auto.run(
    origin=t_cont, 
    starting_point='EP', 
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
    NMX=10000, 
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
    UZR={"PC/RPO_e_in/u": 600}, 
    STOP={}
)
#%%
jrc_auto.plot_continuation('PAR(1)', 'U(1)', cont='u')
plt.show()

# %% Analysis of the state variables 
"""(notation related to Grimbert, F., Faugeras, O., 2006. Bifurcation analysis of Jansen's neural mass model.
Neural Comput. 18, 3052–3068.)

- state variables:

# y1 = U(1)
# y2 = -U(3)
# y0 = U(5)-U(7) ???
"""
# Parameters:--> mV to make computation easier

v0 = 6 # mV
r = 0.560 # mV^-1
e0 = 2.5 # Hz
A = 3.25 # mV
a = 100 # Hz
tau_e = 1/a

y1 = u_sols["PC/RPO_e_in/v"].to_numpy()*1000
y2 = -u_sols["PC/RPO_i/v"].to_numpy()*1000
y = y1-y2
y0 = (A/a)*((2*e0)/(1+np.exp(r*(v0-y))))
p = u_sols["PC/RPO_e_in/u"]

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
