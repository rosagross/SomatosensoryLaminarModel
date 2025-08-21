"""Bifurcation analysis of the original Jansen-Rit Model (.yaml pre-implemented in Pyrates)"""

# %% Libraries import
from pycobi import ODESystem
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyrates.frontend import CircuitTemplate
import os 
os.chdir("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/bifurcation_analysis/""")

# %% Model definition
# path to YAML model definition
model = "model_templates.neural_mass_models.jansenrit.JRC"

# structure of the .yaml model
jrc = CircuitTemplate.from_yaml("model_templates.neural_mass_models.jansenrit.JRC")

# %% 
# installation directory of auto-07p
auto_dir = "/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p"

#ode = ODESystem.from_yaml(model, auto_dir=auto_dir, NPR=100, NMX=30000)

ode = ODESystem.from_yaml(model, auto_dir=auto_dir, init_cont = False, NPR=100, NMX=30000)

# init_cont: does not integrate in time when defining the model (when a system is defined, an initial time integration is performed in order to make it converge to a steady-state solution--> needed for the bifurcation analysis. 
# http://www.scholarpedia.org/article/Equilibrium)
# NPR: print/report frequency
# NMX: maximum number of continuation steps 

"""this is another way to define the system while changing some defined parameters of the model:
ode = ODESystem.from_yaml(
    "model_templates.neural_mass_models.qif.qif_sfa", auto_dir="/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p",
    node_vars={'p/qif_sfa_op/Delta': 2.0, 'p/qif_sfa_op/alpha': 1.0, 'p/qif_sfa_op/eta': 3.0},
    edge_vars=[('p/qif_sfa_op/r', 'p/qif_sfa_op/r_in', {'weight': 15.0*np.sqrt(2.0)})],
    NPR=100, NMX=30000
)"""


print(ode._var_map)
# %% Time integration of a system: structure of the results

# the variables can also be referred with their names in the model!
# - ode.results: gives the results of the simulation as a dict
# - ode._var_map: gives the relationship between the names of the ODESystem and the model (even if "params" is not passed when creating the system)
# - ode.get_summary(): pandas dataframe with the summaries of a continuation
# U(i): state variables
# PAR(14): time 

# plot of the steady-state solution

ode.plot_continuation("t", "U(1)", cont=0)
plt.show()

# addutionally, the timeseries of a variable or a trajectory can be plotted (respectively with .plot_timeseries or .plot_trajectory)
# look up how to refer direcly to the variables
"""ode.plot_continuation("t", "p/qif_sfa_op/r", cont=0)
plt.show()"""

# %% Time continuation:

"""FROM "continuation.py" in pycobi_prova (tutorial in the Pyrates docs)
# In parameter continuations, it is required that you start continuing the parameters from an
# `equilibrium <http://www.scholarpedia.org/article/Equilibrium>`_ or
# `periodic orbit <http://www.scholarpedia.org/article/Periodic_orbit>`_, i.e. that the solution to the
# `initial value problem <http://www.scholarpedia.org/article/Initial_value_problems>`_ would be constant or periodic
# in time. To achieve this, you can either set the initial values and parameters of your system to a known solution in
# the model definition (e.g. the YAML template), or choose a model parameterization for which a finite solution exists
# and then calculate the solution of the model in time until it converges to an equilibrium (or periodic orbit).
"""
t_sols, t_cont = ode.run(
    c='ivp', name='time', DS=1e-4, DSMIN=1e-10, EPSL=1e-08, EPSU=1e-08, EPSS=1e-06,
    DSMAX=1e-2, NMX=1000, UZR={14: 4.0}, STOP={'UZ1'})

ode.plot_continuation('PAR(14)', 'U(1)', cont='time')
plt.show()


# %% 1D parameter continuation:
# origin: indicates the key of the solution branch ??
# starting_point: starting point of the parameter continuation
# name: parameter of interest
# bidirectional: automatic continuation of the solution branch into both directions of the continuation parameter
# others--> specific of Auto07p (<https://github.com/auto-07p/auto-07p/doc>)

# in .run (<https://pycobi.readthedocs.io/en/latest/pycobi.html>) there is a list of all the possible parameters, 
# for example how to get the eigenvalues or the period of the solutions
u_sols, u_cont = ode.run(
    origin=t_cont, starting_point='UZ1', name='u', bidirectional=True,
    ICP=1, RL0=-50.0, RL1=150.0, IPS=1, ILP=1, ISP=2, ISW=1, NTST=400,
    NCOL=4, IAD=3, IPLT=0, NBC=0, NINT=0, NMX=2000, NPR=10, MXBF=5, IID=2,
    ITMX=40, ITNW=40, NWTN=12, JAC=0, EPSL=1e-06, EPSU=1e-06, EPSS=1e-04,
    DS=1e-4, DSMIN=1e-8, DSMAX=5e-2, IADS=1, THL={}, THU={}, UZR={}, STOP={}
)
# see email 12/08/2025 for the meaning of the parameters

# u_sols is a dataframe with the results of the continuation

# %% plot of the continuation
ode.plot_continuation('PAR(4)', 'U(1)', cont='u')
plt.show()

"""ode.plot_continuation("p/qif_sfa_op/eta", "p/qif_sfa_op/r", cont="eta")
plt.show()"""

# solid line: stable / dotted: unstable
# AUTO detects the bifurcations:
# - gray triangles: fold ("saddle-node") bifurcations (the critical eigenvalue of the vector field defined by the right-hand sides of the model´s ODE
#  cross the imaginary axis = the real part becomes positive/negative) --> CHANGE OF STABILITY IN THE SOLUTION
#  <http://www.scholarpedia.org/article/Saddle-node_bifurcation>
# - green circle: Hopf bifurcations <http://www.scholarpedia.org/article/Andronov-Hopf_bifurcation>


# %% CONTINUE "QUIF_SFA.PY" FOR THE LIMIT CYCLES
