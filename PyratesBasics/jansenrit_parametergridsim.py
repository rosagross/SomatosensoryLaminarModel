# %% 
'''
This code runs simulations of the basic Jansen-Rit model with systematically changing the parameters!
It will run a loop with one simulations for every parameter setting.
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Simulation parameters
simulation_time = 2 # in s
step_size = 1e-4 # in s

# Model parameters
m_max = 2.5 # in s^-1! (so 1 means 1 per second)
r = 560 # sigmoidal steepness in V^-1 (560 here means 0.56 mV^-1)
v_thr = 0.006 # firing threshold in V
H_E = 0.00325 # in V (e.g. the 3.25mV from the paper are 0.00325V)
H_I = -0.022 
tau_Es = [0.0005, 0.001, 0.01, 0.1, 0.2] # in s^-1 
tau_I = 0.02

# Output matrices to store computed values

# define time steps (starting with the second step since the first one is initialized already)
steps = np.arange(2*step_size, simulation_time+step_size, step_size)

# rates & potentials (E, IIN , EIN) 
nPop = 3
rate = np.zeros((nPop, len(steps)))
potential = np.zeros((nPop, nPop, len(steps))) 

# first derivatives and time constants (E, IIN , EIN), here we don't need to save the entire time line
u_t = np.zeros((nPop, nPop)) # the initial first-order derivative: v'(t) = u(t)
u_dot = np.zeros((nPop, nPop)) # we don't need an initial value for this one actually but a place holder array
H = [H_E, H_I, H_E]
W = np.array([[0, 33.75, 108],[33.75, 0, 0],[135, 0, 0]]) #np.array([[0, 33.75, 108],[33.75, 0, 0],[135., 0, 0]])

# Simulation loop 

# Initialize first values with 0 or randomly
rate[:, 0] = 0
potential[:, 0] = 0 # np.random.normal(0,1,(3,round(stimulation_time/step_size))) # noise
for tau_E in tau_Es:

    tau = [tau_E, tau_I, tau_E]

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
        rate_current = 2*m_max / (1 + np.exp(r*(v_thr - np.sum(v_current, axis=1)))) # All populations have the same PRO

        rate[:,iter] = rate_current
        potential[:, :, iter] = v_current
    
    # safe the results
    # sum together the potential for all
    potential_sum = np.sum(potential, axis=1)
    potential_df = pd.DataFrame(potential_sum.T, columns=['PC', 'IIN', 'EIN'])
    potential_df.to_csv(f'Outputs/simpleJR_potentials_tauE{tau_E}_tauI{tau_I}.csv', index=False)

    # save the firing rate 
    rate_df = pd.DataFrame(rate.T, columns=['PC', 'IIN', 'EIN'])
    rate_df.to_csv(f'Outputs/simpleJR_rates_tauE{tau_E}_tauI{tau_I}.csv', index=False)

# %%
plt.plot(potential[0, 1]+potential[0, 2], label='E')
plt.plot(potential[0, 1], label='I')
plt.plot(potential[0, 2], label='EI')
plt.legend()
plt.show()

    
    







# %%
