import numpy as np
import matplotlib.pyplot as plt
import parameters

params = parameters.Parameter()


# Simulation parameters
tau = params.tau
nPop = params.nPop
simulation_time = params.sim_dur # in s
step_size = params.step_size # in s


# sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
sigm = params.sigmoid_params

# Synaptic kernel 
H = np.ones((nPop, nPop+1))

# external input matrix
Iext = params.Iext

# define time steps 
steps = np.arange(step_size, simulation_time+step_size, step_size)

# Output matrices to store computed values for rates & potentials (E, IIN , EIN) 
rate = np.zeros((nPop, len(steps)))
potential = np.zeros((nPop, nPop+1, len(steps))) 

# Simulation loop
# Initialize first values for the potential, rate and first order derivative with 0 or randomly
v_current = np.zeros((nPop, nPop+1))
rate_current = np.zeros(nPop)
u_t = np.zeros((nPop, nPop+1)) # the initial first-order derivative: v'(t) = u(t)

coupling_strengths = [100]# np.arange(0, 100, 5)

# arrays to store min and max rate
minRate = np.zeros((nPop, len(coupling_strengths)))
maxRate = np.zeros((nPop, len(coupling_strengths)))

for sim_iter, g in enumerate(coupling_strengths):

    # Weight matrix [to x from]
    W = params.get_connectivity(params.P, params.S, params.C, g) 

    for timestep, time in enumerate(steps):
        
        # Update RATE (calculated from the current potential)
        # technically this could also go to the end of the for-loop but we need a rate value calculated from the initial potential value
        for i in range(nPop):
            rate_current[i] = sigm[i][2] / (1 + np.exp(sigm[i][0]*(sigm[i][1] - np.sum(v_current[i,:])))) 

        # Save the new values
        rate[:, timestep] = rate_current 
        potential[:, :, timestep] = v_current

        # We go through every population and evaluate the new membrane potential based on the connectivity
        for i in range(nPop):
            for j in range(nPop):

                # 1. Calculate new POTENTIAL (calculated using the current first derivative) - but it's only updated/saved in the next step
                v_dot = u_t[i, j]
                v_current[i, j] = v_current[i, j] + v_dot * step_size
                
                # 2. Update the first derivative based on the current potential (NOT the just updated one!) current first derivative
                u_dot = (H[i, j]/tau[j]) * (W[i, j]*rate_current[j]) - 2 * u_t[i, j]/tau[j] - potential[i, j, timestep]/(tau[j]**2)
                u_t[i, j] = u_t[i,j] + u_dot * step_size

            # Add external input 
            u_dot = (H[i,-1]/tau[-1]) * (W[i, -1]) * Iext[i, timestep] - 2 * u_t[i, j]/tau[-1] - potential[i, -1, timestep]/(tau[-1]**2)
            u_t[i, -1] = u_t[i,-1] + u_dot * step_size
    
    minRate[:,sim_iter] = np.min(rate[:,-100:],axis=1)
    maxRate[:,sim_iter] = np.max(rate[:,-100:],axis=1)


# plot potentials
fig, axs = plt.subplots(4, 1, figsize=(3, 6))  # Set figure size

# Plot settings for all subplots
for i, ax in enumerate(axs, start=1):
    ax.plot(steps, rate[(i-1)*4:i*4, :].T, linewidth=2)
    ax.grid(True)
    ax.set_ylabel('Hz')
    ax.legend(['L2/3', 'L4', 'L5', 'L6'])

# Set titles for each subplot
axs[0].set_title('Baseline rate (E)')
axs[1].set_title('Baseline rate (PV)')
axs[2].set_title('Baseline rate (SOM)')
axs[3].set_title('Baseline rate (VIP)')
plt.tight_layout() 
plt.legend()
plt.show()




    
    






