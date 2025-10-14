from area1_prova import area_1prova
from scipy.integrate import solve_ivp
import os 
os.chdir("/data/hu_mecozzi/Documents/SomatosensoryLaminarModel/PyratesBasics/exp_model/""") 
from pyrates.frontend import OperatorTemplate, NodeTemplate, CircuitTemplate
from copy import deepcopy
from parameters import Parameter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numba import njit
from yaml_saving import circuit_to_yaml
from pprint import pprint

cells = ['E1', 'PV1', 'SST1', 'VIP', 'E2', 'PV2', 'SST2', 'E3', 'PV3', 'SST3', 'E4', 'PV4', 'SST4']

N_cells = len(cells)
params = Parameter()

# %%
sigm = params.get_sigmoid("S1") #already in the correct order

r = []
v_thr = []
m_max = []
for row in range(N_cells):
    r.append(sigm[row,0])      
    v_thr.append(sigm[row,1]) 
    m_max.append(sigm[row,2]) 

# %%
tau,_ = params.get_params("S1")           
tau_a1 = tau[0] 
# default order of the taus: E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1
#tau_order = ['E1', 'E2', 'E3', 'E4', 'PV1', 'PV2', 'PV3', 'PV4', 'SST1', 'SST2', 'SST3', 'SST4', 'VIP']
#indices = [tau_order.index(c) for c in cells]

#tau_a1 = tau[indices]

# %%
bEI = 0.5
connect_reverse_factor =  6448 
g = 100.0 # (g)
gE = g * bEI /connect_reverse_factor
gI = g * (1 - bEI) /connect_reverse_factor
gEthal = 0
gIthal = 0
thal_connect = (0, 0, 0, 0)  # tEE, tEI, tIE, tII

W = params.get_connectivity(gE, gI, gEthal, gIthal, thal_connect, include_Iext=True) 

# selecting the region --> A1: FROM THE 5TH ELEMENT TO THE 17TH (INCLUDED) 
# in python the results include the start index but excludes the end index
W_A1 = W[4:17,4:17]
results = solve_ivp(fun=area_1prova, t_span=(t0, T), y0=y, first_step=dt, args=(W_A1, tau_a1, params))