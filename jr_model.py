import numpy as np
import matplotlib.pyplot as plt
from parameters import Parameter

class JR_Model():

    def __init__(self, Iext, cortex_type, step_size=0.001, simulation_time=1):
        
        # load in all parameters
        self.p = Parameter(cortex_type)

        # Simulation parameters
        self.tau = self.p.tau
        self.nPop = self.p.nPop
        self.simulation_time = simulation_time# in s
        self.step_size = step_size # in s

        # sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
        self.sigm = self.p.sigmoid_params

        # Synaptic kernel 
        self.H = np.ones((self.nPop, self.nPop+1))

        # define time steps 
        self.steps = np.arange(self.step_size, self.simulation_time+self.step_size, self.step_size)

        # external input matrix
        self.Iext = np.tile(Iext, (self.nPop,1))

    def run_simulation(self, g):
        '''
        g: coupling strength
        '''
        # Output matrices to store computed values for rates & potentials (E, IIN , EIN) 
        self.rate = np.zeros((self.nPop, len(self.steps)))
        self.potential = np.zeros((self.nPop, self.nPop+1, len(self.steps))) 

        # Simulation loop
        # Initialize first values for the potential, rate and first order derivative with 0 or randomly
        self.v_current = np.zeros((self.nPop, self.nPop+1))
        self.rate_current = np.zeros(self.nPop)
        self.u_t = np.zeros((self.nPop, self.nPop+1)) # the initial first-order derivative: v'(t) = u(t)

        # Weight matrix [to x from]
        W = self.p.get_connectivity(g) 

        for timestep, time in enumerate(self.steps):
            
            # Update RATE (calculated from the current potential)
            # technically this could also go to the end of the for-loop but we need a rate value calculated from the initial potential value
            for i in range(self.nPop):
                # the incoming potential has to be defined in mV because of how the sigmoid parameter are defined (also in mV!)
                self.rate_current[i] = self.sigm[i][2] / (1 + np.exp(self.sigm[i][0]*(self.sigm[i][1] - np.sum(self.v_current[i,:])*1e3)))  

            # Save the new values
            self.rate[:, timestep] = self.rate_current 
            self.potential[:, :, timestep] = self.v_current

            # We go through every population and evaluate the new membrane potential based on the connectivity
            for i in range(self.nPop):
                for j in range(self.nPop):

                    # 1. Calculate new POTENTIAL (calculated using the current first derivative) - but it's only updated/saved in the next step
                    v_dot = self.u_t[i, j]
                    self.v_current[i, j] = self.v_current[i, j] + v_dot * self.step_size
                    
                    # 2. Update the first derivative based on the current potential (NOT the just updated one!) current first derivative
                    u_dot = (self.H[i, j]/self.tau[i,j]) * (W[i, j]*self.rate_current[j]) - 2 * self.u_t[i, j]/self.tau[i,j] - self.potential[i, j, timestep]/(self.tau[i,j]**2)
                    self.u_t[i, j] = self.u_t[i,j] + u_dot * self.step_size

                # Add external input 
                v_dot = self.u_t[i, -1]
                self.v_current[i, -1] = self.v_current[i, -1] + v_dot * self.step_size
                u_dot = (self.H[i,-1]/self.tau[i,-1]) * ((W[i, -1]) * self.Iext[i, timestep]) - 2 * self.u_t[i, -1]/self.tau[i,-1] - self.potential[i, -1, timestep]/(self.tau[i,-1]**2)
                self.u_t[i, -1] = self.u_t[i,-1] + u_dot * self.step_size
        
        
        return self.rate, self.potential
    


    





    
    






