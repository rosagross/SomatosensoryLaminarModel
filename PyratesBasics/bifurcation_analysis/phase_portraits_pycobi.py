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

# %%
"""
jrc.get_run_func(func_name='jrc_rhs', file_name='jrc_pprova', step_size=1e-4, auto=True, backend='python', solver='scipy',
                 vectorize=False, float_precision='float64')
best version:                 
jrc.get_run_func(func_name='jrc_v1', step_size=1e-4, file_name='jrc_pprova1', backend='numpy',
                 vectorize=False, float_precision='float64')

Generate a function that evaluates the vector field of the dynamical system represented by this
        `CircuitTemplate` instance.

        Parameters
        ----------
        func_name
            Name of the vector field evaluation function.
        step_size
            Integration step-size. Required for the implementation of the extrinsic inputs.
        inputs
            Dictionary providing extrinsic, time-dependent inputs to the system. Keys are the names of system variables,
            following the `*circuit/node/op/var` notation. Values are 1D numpy arrays that represent the input over time
            with time steps of size `step_size`.
        backend
            Name of the backend that should be used for implementing the system equations. Possible choices are:
                - 'default' or 'numpy': A backend based on `numpy` functions, representing all system variables as
                    `np.ndarray`.
                - 'tensorflow': A backend that represents the system equations as a `tensorflow` graph, in which all
                    system variables are stored as `tf.constant` or `tf.Variable`.
                - 'torch': A backend based on `pytorch` which represents all variables as `torch.tensor`.
                - 'fortran': Translates all system variables and equations into Fortran90 equivalents and uses
                    `numpy.f2py` to make them available via Python. Requires `vectorize` to be set to `False`.
                - 'julia': Translates all system variables and equations into Julia equivalents and uses `PyJulia` to
                    make them available via Python. Requires `vectorize` to be set to `False`. Also requires that
                    the path to the julia executable that should be used for the simulation is provided via the
                    keyword argument `julia_path`.
        vectorize
            If true, nodes that are governed by the same equation sets will be grouped and the respective equations will
            be vectorized. If false, all equations will be scalar in nature.
        verbose
            If true updates regarding the status of the `run` procedure will be displayed.
        clear
            If true, all cached variables will be freed and all temporary files will be deleted after the `run`
            procedure. To inspect the vector field evaluation function, `clear` should be set to `False`.
        in_place
        kwargs
            Additional keyword arguments.

        Returns
        -------
        Tuple[Callable, tuple, tuple, dict]
            The vector field evaluation function, all its positional arguments, the argument keys, and the indices of
            the different state variables in the state vector.

        """

# result of this function:

# see how to call it from the file
# %%
# -------------------------------------------------------------------------------------------------------------------------------------------------------
from numpy import pi, sqrt
from numpy import exp

def jrc_rhs(y,u,tau,H,tau_v1,H_v1,tau_v2,H_v2,tau_v3,H_v3,V_thr,r,m_max,weight_v2,V_thr_v1,r_v1,m_max_v1,weight_v3,V_thr_v2,r_v2,m_max_v2,weight,weight_v1):
    """{'PC/RPO_e_in/v': 0, # y1
    'PC/RPO_e_in/i': 1,
    'PC/RPO_i/v': 2, # y2
    'PC/RPO_i/i': 3,
    'EIN/RPO_e/v': 4, # y0_ein
    'EIN/RPO_e/i': 5,
    'IIN/RPO_e/v': 6, # y0_iin
    'IIN/RPO_e/i': 7}
    """
    v = y[0] # y1 
    i = y[1]
    v_v1 = y[2] # y2
    i_v1 = y[3]
    v_v2 = y[4] # y0_ein
    i_v2 = y[5]
    v_v3 = y[6] # y0_iin
    i_v3 = y[7]

    dy = np.zeros_like(y) 

    m_out_v2 = 2.0*m_max/(exp(r*(V_thr - v - v_v1)) + 1)
    m_in_v2 = m_out_v2*weight_v2
    m_out = 2.0*m_max_v1/(exp(r_v1*(V_thr_v1 - v_v2)) + 1)
    m_in_v3 = m_out_v2*weight_v3
    m_out_v1 = 2.0*m_max_v2/(exp(r_v2*(V_thr_v2 - v_v3)) + 1)
    m_in = m_out*weight
    m_in_v1 = m_out_v1*weight_v1
	
    dy[0] = i 
    dy[1] = H*(m_in + u)/tau - 2*i/tau - v/tau**2
    dy[2] = i_v1
    dy[3] = H_v1*m_in_v1/tau_v1 - 2*i_v1/tau_v1 - v_v1/tau_v1**2
    dy[4] = i_v2
    dy[5] = H_v2*m_in_v2/tau_v2 - 2*i_v2/tau_v2 - v_v2/tau_v2**2
    dy[6] = i_v3
    dy[7] = H_v3*m_in_v3/tau_v3 - 2*i_v3/tau_v3 - v_v3/tau_v3**2

    return dy
# -------------------------------------------------------------------------------------------------------------------------------------------------------

# %%
C= 135
params = dict(
    tau=0.01, H=0.00325,
    tau_v1=0.01, H_v1=0.00325,
    tau_v2=0.01, H_v2=0.00325,
    tau_v3=0.02, H_v3=-0.022,
    V_thr=6e-3, r=560.0, m_max=2.5, weight_v2= C,
    V_thr_v1=6e-3, r_v1=560.0, m_max_v1=2.5, weight_v3= 0.25*C,
    V_thr_v2=6e-3, r_v2=560.0, m_max_v2=2.5, weight= 0.8*C,
    weight_v1= 0.25*C
)

u = 0.0

# %% grid construction
y0=np.linspace(-0.2,0.2,20)
y1=np.linspace(-0.2,0.2,20)
y2=np.linspace(-0.2,0.2,20)
y3=np.linspace(-0.2,0.2,20)
y4=np.linspace(-0.2,0.2,20)
y5=np.linspace(-0.2,0.2,20)
y6=np.linspace(-0.2,0.2,20)
y7=np.linspace(-0.2,0.2,20)
# grid for the phase portrait
Y0,Y1,Y2,Y3,Y4,Y5,Y6,Y7 =np.meshgrid(y0,y1,y2,y3,y4,y5,y6,y7)

dY0=np.zeros(Y0.shape)
dY1=np.zeros(Y1.shape)
dY2=np.zeros(Y2.shape)
dY3=np.zeros(Y3.shape)
dY4=np.zeros(Y4.shape)
dY5=np.zeros(Y5.shape)
dY6=np.zeros(Y6.shape)
dY7=np.zeros(Y7.shape)

#%%
n_iterations = Y0.size
dy_matrix = np.zeros((8, n_iterations))
for idx, (yy0, yy1, yy2, yy3) in enumerate(zip(Y0.ravel(), Y1.ravel(), Y2.ravel(), Y3.ravel())):
    y = np.array([yy0, 0, yy1, 0, yy2, 0, yy3, 0])
    dy = jrc_rhs(
        y, u,
        params["tau"], params["H"],
        params["tau_v1"], params["H_v1"],
        params["tau_v2"], params["H_v2"],
        params["tau_v3"], params["H_v3"],
        params["V_thr"], params["r"], params["m_max"], params["weight_v2"],
        params["V_thr_v1"], params["r_v1"], params["m_max_v1"], params["weight_v3"],
        params["V_thr_v2"], params["r_v2"], params["m_max_v2"], params["weight"],
        params["weight_v1"]
    )
    dy_matrix[:, idx] = dy
    #dV.ravel()[idx] = dy[0]
    #dI.ravel()[idx] = dy[1]

"""import numpy as np
import itertools

# %% grid construction
y0 = np.linspace(-0.2, 0.2, 20)
y1 = np.linspace(-0.2, 0.2, 20)
y2 = np.linspace(-0.2, 0.2, 20)
y3 = np.linspace(-0.2, 0.2, 20)
y4 = np.linspace(-0.2, 0.2, 20)
y5 = np.linspace(-0.2, 0.2, 20)
y6 = np.linspace(-0.2, 0.2, 20)
y7 = np.linspace(-0.2, 0.2, 20)

# number of total points
n_iterations = len(y0) * len(y1) * len(y2) * len(y3) * len(y4) * len(y5) * len(y6) * len(y7)

# container for results (8 x N matrix)
dy_matrix = np.zeros((8, n_iterations))

# %% iterate through the grid without building a huge array
for idx, point in enumerate(itertools.product(y0, y1, y2, y3, y4, y5, y6, y7)):
    # point is a tuple (yy0, yy1, ..., yy7)
    y = np.array(point)

    dy = jrc_rhs(
        y, u,
        params["tau"], params["H"],
        params["tau_v1"], params["H_v1"],
        params["tau_v2"], params["H_v2"],
        params["tau_v3"], params["H_v3"],
        params["V_thr"], params["r"], params["m_max"], params["weight_v2"],
        params["V_thr_v1"], params["r_v1"], params["m_max_v1"], params["weight_v3"],
        params["V_thr_v2"], params["r_v2"], params["m_max_v2"], params["weight"],
        params["weight_v1"]
    )

    dy_matrix[:, idx] = dy
"""
# %%
#   RELAZIONE TRA LE VARIABILI

# PLOT

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
 