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
u_pce = t_sols["PC/RPO_e_in/i"]
v_pci = t_sols["PC/RPO_i/v"]
u_pci = t_sols["PC/RPO_i/i"]
pc = v_pce + v_pci
t = t_sols["t"]
plt.plot(t,v_pce,label='V_PCE')
plt.plot(t,v_pci,label='V_PCI')
plt.plot(t,pc,label='PC')
plt.ylabel('[V]')
plt.xlabel('[s]')
plt.legend()
plt.show()


# %% Phase portraits

""" usare questo in una griglia
qif.get_run_func(func_name='qif_rhs', file_name='qif', step_size=1e-4, auto=True, backend='fortran', solver='scipy',
                 vectorize=False, float_precision='float64')"""

# next, define a grid of points at which we will show arrows
x0=np.linspace(-2,2,20)
x1=np.linspace(-2,3,20)

# create a grid
X0,X1=np.meshgrid(v_pce,v_pci)

# %%
# projections of the trajectory tangent vector 
dX0=np.zeros(X0.shape)
dX1=np.zeros(X1.shape)

shape1,shape2=X1.shape

for indexShape1 in range(shape1):
    for indexShape2 in range(shape2):
        dxdtAtX=dynamicsStateSpace([X0[indexShape1,indexShape2],X1[indexShape1,indexShape2]],0)
        dX0[indexShape1,indexShape2]=dxdtAtX[0]
        dX1[indexShape1,indexShape2]=dxdtAtX[1]
      
# %%        
#adjust the figure size
plt.figure(figsize=(8, 8))
# plot the phase portrait
plt.quiver(X0,X1,u_pce,u_pci,color='b')
# adjust the axis limits
#plt.xlim(-2,2)
#plt.ylim(-2,2)
# insert the title
plt.title('Phase Portrait', fontsize=14)
# set the axis labels
plt.xlabel('$x_{1}$',fontsize=14)
plt.ylabel('$x_{2}$',fontsize=14)
# adjust the font size of x and y axes
plt.tick_params(axis='both', which='major', labelsize=14)
# insert legend
# plt.legend(fontsize=14)
# save figure
# plt.savefig('phasePortrait.png',dpi=600)
plt.show()

# %%
# now simulate the dynamics for a certain starting trajectory and 
# plot the dynamics on the same graph
initialState=np.array([-1,-1])
simulationTime=np.linspace(0,2,200)
# generate the state-space trajectory
solutionState=odeint(dynamicsStateSpace,initialState,simulationTime)


#adjust the figure size
plt.figure(figsize=(8, 8))
# plot the phase portrait
plt.quiver(X0,X1,dX0,dX1,color='b')
# adjust the axis limits
plt.xlim(-2,2)
plt.ylim(-2,2)
# add the state trajectory plot
plt.plot(solutionState[:,0], solutionState[:,1], color='r',linewidth=3)
# insert the title
plt.title('Phase Portrait', fontsize=14)
# set the axis labels
plt.xlabel('$x_{1}$',fontsize=14)
plt.ylabel('$x_{2}$',fontsize=14)
# adjust the font size of x and y axes
plt.tick_params(axis='both', which='major', labelsize=14)
# insert legend
# plt.legend(fontsize=14)
# save figure
# plt.savefig('phasePortraitStateTrajectory.png',dpi=600)
plt.show() 
 