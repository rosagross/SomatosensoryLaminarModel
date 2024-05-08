# %% 

import numpy as np
import matplotlib.pyplot as plt
from parameters import Parameter
import pandas as pd

safe_results = True

# Simulation parameters
simulation_time = 2 # in s
step_size = 1e-4 # in s

# Model parameters
cortex_type = 'visual'
params = Parameter(cortex_type)
sigm = params.get_sigmoid()

nPop = 4
r = []
v_thr = []
m_max = []
for row in range(nPop):
    r.append(sigm[row,0])      
    v_thr.append(sigm[row,1]) 
    m_max.append(sigm[row,2]) 

r = np.array(r)
v_thr = np.array(v_thr)
m_max = np.array(m_max)


H_E = 1 # in V (e.g. the 3.25mV from the paper are 0.00325V)
H_I = 1
tau_E = 0.01 # in s^-1
tau_I = 0.01


# Output matrices to store computed values

# define time steps (starting with the second step since the first one is initialized already)
steps = np.arange(2*step_size, simulation_time+step_size, step_size)

# rates & potentials (E, IIN , EIN) 
rate = np.zeros((nPop, len(steps)))
potential = np.zeros((nPop, nPop, len(steps))) 

# first derivatives and time constants (E, IIN , EIN), here we don't need to save the entire time line
u_t = np.zeros((nPop, nPop)) # the initial first-order derivative: v'(t) = u(t)
u_dot = np.zeros((nPop, nPop)) # we don't need an initial value for this one actually but a place holder array
tau,_ = params.get_params()                  #[tau_E, tau_I, tau_E, tau_E]
tau = tau[0]
H = [H_E, H_I, H_E, H_E]
W = params.get_connectivity(1, include_Iext=False)
# Simulation loop 

# Initialize first values with 0 or randomly
rate[:, 0] = 0
potential[:, 0] = 0 # np.random.normal(0,1,(3,round(stimulation_time/step_size))) # noise


for iter, time in enumerate(steps):

    # store current rate and potential for now
    rate_current = rate[:, iter-1] 
    v_current = potential[:, :, iter-1] 

    # DERIVATES
    # the current 
    # RPO: update the population's membrane potential (calculated using the previous rate): u'(t) = v''(t)
    for i in range(nPop):
        for j in range(nPop):
            u_dot[i,j] = H[j]/tau[j] * (W[i,j]*rate_current[j]) - 2 * u_t[i, j]/tau[j] - v_current[i,j]/(tau[j]**2)
    
    v_current = v_current + u_t * step_size
    u_t = u_t + u_dot * step_size

    # PRO: update the population's RATE (calculated from the calculated **updated** potential)
    rate_current = 2*m_max / (1 + np.exp(r*(v_thr - np.sum(v_current*1000, axis=1)))) # All populations have the same PRO

    rate[:,iter] = rate_current
    potential[:, :, iter] = v_current


if safe_results:
    cells = np.array(['E1', 'E2', 'E3', 'E4'])
    
    potential_sum = np.zeros((nPop, len(steps))) # (16x1000)
    for i in range(nPop):  #alle Eingänge für jede Zielpopulation 
            potential_sum[i] = np.sum(potential[i], axis=0)

    potential_df = pd.DataFrame(potential_sum.T, columns=cells)
    potential_df.to_csv(f'output/python_4.csv', index=False)
# %%
plt.plot(potential[2, 0], label='E1') #alle Eingänge (target, source)
plt.plot(potential[2, 1], label='P1')
plt.plot(potential[2, 2], label='S1')
plt.plot(potential[2, 3], label='V1')
plt.legend()
plt.show()


    






