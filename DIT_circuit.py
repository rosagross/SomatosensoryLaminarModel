'''
In this script I implement a mini DIT circuit, with 3 interneuron populations and 1 pyramidal cell population.

Connectivity:

Integration function:

External input:



'''
# %%
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns

# %%

# number of populations
nPop = 4
# coupling strength
g = 1

# Connectivity Parameter
S = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
P = np.array([[0, 0, 0, 0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
C = np.array([0, 0, 0, 0])
PS = P * S
W = PS * np.tile(C, (nPop,1))

# input weights 
Wext_basal = np.array([1, 0, 0, 0])
Wext_apical = np.array([1, 0, 0, 0])
W = np.concatenate((W*g,np.expand_dims(Wext_basal, 1)), axis=1)
W = np.concatenate((W*g,np.expand_dims(Wext_apical, 1)), axis=1)
print('W',W)

# Time constants
tau = [6,3,20,15, 3, 3] 
tau = np.tile(tau, (nPop+2,1))

# sigmoid function (parameters: r, v_thr, m_max)
sigm = [[  0.12782346,  32.10540543,  31.39696397], # E
           [  0.14218422,  40.03107351, 166.82960408], # PV
           [  0.07937015,  42.01276379,  56.95305832], # SST
           [  0.0704119 ,  37.86409387,  38.52689646]] # VIP

# Simulation paraneters
simulation_time = 3 # in s
step_size = 0.001 # in s

# Synaptic kernel 
H = np.ones((nPop, nPop+1))

# define time steps 
steps = np.arange(step_size, simulation_time+step_size, step_size)

# external input matrix APICAL
input_duration_apical = 0.5
input_strength_apical = 100
input_onset_apical = 1
Iext_apical = np.zeros(int(simulation_time/step_size))
t_apical  = int(input_duration_apical/step_size)
t0_apical = int(input_onset_apical/step_size)
Iext_apical[t0_apical:t0_apical+t_apical] = input_strength_apical
Iext_apical = np.tile(Iext_apical, (nPop,1))

# external input matrix BASAL
input_duration_basal = 0.3
input_strength_basal = 100
input_onset_basal = 1
Iext_basal = np.zeros(int(simulation_time/step_size))
t_basal  = int(input_duration_basal/step_size)
t0_basal = int(input_onset_basal/step_size)
Iext_basal[t0_basal:t0_basal+t_basal] = input_strength_basal
Iext_basal = np.tile(Iext_basal, (nPop,1))

# %% arrays to store simulated potentials and rates
rate_current = np.zeros(nPop)
rate = np.zeros((nPop, len(steps)))
potential = np.zeros((nPop, nPop+2, len(steps))) 
v_current = np.zeros((nPop, nPop+2))
u_t = np.zeros((nPop, nPop+2)) # the initial first-order derivative: v'(t) = u(t)


# Simulation loop 
for timestep, time in enumerate(steps):
    
    # Potential-to-Rate
    for i in range(nPop):
        
        #print('population', i)
        # integrate the input from basal and apical dendrites
        integrated_v = v_current[i,-1] + v_current[i,-2]
        sigmoid_v = 1 / (1 + np.exp(10*(1-integrated_v*1e3)))
        

        # after integrating it we store it in the previous last position of the array (for summation in the next step)  
        v_current[i,-2] = integrated_v * sigmoid_v

        if (i == 0) and (timestep*step_size <1.1):
            print('timestep', timestep*step_size)
            print('summed Inputs', v_current[i,-1]+v_current[i,-2])
            print(sigmoid_v)


        
        # TODO: role of the second order thalamus: the second order thalamus controls the detection threshold. 
        # It could do this by shifting the threshold parameter in the sigmoid function. Strong input from 
        # the POm would lower the threshold

        # the incoming potential has to be defined in mV because of how the sigmoid parameter are defined (also in mV!)
        rate_current[i] = sigm[i][2] / (1 + np.exp(sigm[i][0]*(sigm[i][1] - np.sum(v_current[i,:-1])*1e3)))  
        if (timestep/step_size < 1.5) and (i == 0):
            print('sum', np.sum(v_current[i,:-1]))
            print(rate_current[i])
        
            #print(v_current[i,:-1])
            #print('rate', rate_current[i])

    # Save the new values
    rate[:, timestep] = rate_current 
    potential[:, :, timestep] = v_current

    # Rate-to-Potential
    # We go through every population and evaluate the new membrane potential based on the connectivity
    for i in range(nPop):
        for j in range(nPop):

            # 1. Calculate new POTENTIAL (calculated using the current first derivative) - but it's only updated/saved in the next step
            v_dot = u_t[i, j]
            v_current[i, j] = v_current[i, j] + v_dot * step_size
            
            # 2. Update the first derivative based on the current potential (NOT the just updated one!) current first derivative
            u_dot = (H[i, j]/tau[i,j]) * (W[i, j]*rate_current[j]) - 2 * u_t[i, j]/tau[i,j] - potential[i, j, timestep]/(tau[i,j]**2)
            u_t[i, j] = u_t[i,j] + u_dot * step_size

        # BASAL INPUT
        v_dot_basal = u_t[i, -2]
        v_current[i, -2] = v_current[i, -2] + v_dot_basal * step_size
        u_dot_basal = (H[i,-2]/tau[i,-2]) * (W[i, -2] * Iext_basal[i, timestep]) - 2 * u_t[i, -2]/tau[i,-2] - potential[i, -2, timestep]/(tau[i,-2]**2)
        u_t[i, -2] = u_t[i,-2] + u_dot_basal * step_size

        # DISTAL APICAL INPUT
        v_dot_apical = u_t[i, -1]
        v_current[i, -1] = v_current[i, -1] + v_dot_apical * step_size
        u_dot_apical = (H[i,-1]/tau[i,-1]) * (W[i, -1] * Iext_apical[i, timestep]) - 2 * u_t[i, -1]/tau[i,-1] - potential[i, -1, timestep]/(tau[i,-1]**2)
        u_t[i, -1] = u_t[i,-1] + u_dot_apical * step_size
        

print(rate.shape)
print('final', rate)
labels = ['E', 'PV', 'SST', 'VIP']

fig, axs = plt.subplots(2,1)
print(Iext_apical.shape)
# plot the distal apical and basal input and the potential in the two compartments
axs[0].plot(steps[:], Iext_apical[0][:], label='apical', color=sns.color_palette('Dark2')[3])
axs[0].plot(steps[:], Iext_basal[0][:], label='basal', color=sns.color_palette('Dark2')[2])

#axs[1].plot(steps, rate[0], label='E')


for target, label in zip(rate, labels):
    axs[1].plot(steps[:], target[:], label=label)

plt.legend()
plt.show()