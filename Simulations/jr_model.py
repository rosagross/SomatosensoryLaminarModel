import numpy as np
import matplotlib.pyplot as plt
from parameters import Parameter
import os

class JR_Model():

    def __init__(self, Iext, Ib, gE, gI, gEthal, gIthal, thal_connect, extI_cellcount, bI_cellcount, thal_cellcount, step_size=0.001, simulation_time=1, area='all'):
        
        # load in all parameters
        self.p = Parameter()

        # Simulation parameters
        self.tau = self.p.tau
        self.nPop = self.p.nPop
        self.simulation_time = simulation_time# in s
        self.step_size = step_size # in s
        self.gE = gE
        self.gI = gI
        self.gEthal = gEthal
        self.gIthal = gIthal
        self.thal_connect = thal_connect
        self.extI_cellcount = extI_cellcount
        self.bI_cellcount = bI_cellcount
        self.thal_cellcount = thal_cellcount

        # sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
        self.sigm = self.p.sigmoid_params

        # Synaptic kernel 
        self.H = np.ones((self.nPop, self.nPop+1))

        # define time steps 
        self.steps = np.arange(self.step_size, self.simulation_time+self.step_size, self.step_size)

        # external input matrix and background input matrix
        # update parameters based on file
        #self.__dict__.update(params)
        # extend input arrays
        self.Iext = np.tile(Iext, (self.nPop,1))
        self.Ib = np.tile(Ib, (self.nPop,1))

        # the area to isolate (default is just the entire model)
        self.area = area
        

    def run_simulation(self):
        '''
        Simulation loop
        '''
        # Output matrices to store computed values for rates & potentials (E, IIN , EIN) 
        self.rate = np.zeros((self.nPop, len(self.steps)))
        self.potential = np.zeros((self.nPop, self.nPop+2, len(self.steps))) 

        # Simulation loop
        # Initialize first values for the potential, rate and first order derivative with 0 or randomly
        self.v_current = np.zeros((self.nPop, self.nPop+2)) # +2 because 1 for background input and one for external input 
        self.rate_current = np.zeros(self.nPop)
        self.u_t = np.zeros((self.nPop, self.nPop+2)) # the initial first-order derivative: v'(t) = u(t)

        # Weight matrix [to x from]
        W = self.p.get_connectivity(self.gE, self.gI, self.gEthal, self.gIthal, self.thal_connect, self.extI_cellcount, self.bI_cellcount, self.thal_cellcount, area=self.area) 

        last_step = self.steps[-1] 

        for timestep, time in enumerate(self.steps):
            
            # Update RATE (calculated from the current potential)
            # technically this could also go to the end of the for-loop but we need a rate value calculated from the initial potential value
            for i in range(self.nPop):

                # the incoming potential has to be defined in mV because of how the sigmoid parameter are defined (also in mV!)
                self.rate_current[i] = self.sigm[i][2] / (1 + np.exp(self.sigm[i][0]*(self.sigm[i][1] - np.sum(self.v_current[i,:]))))  
                
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

                # Add external input (goes to thalamus only, so it's kind of a brain stem input) 
                v_dot = self.u_t[i, -1] # the -1 is the external input 
                self.v_current[i, -1] = self.v_current[i, -1] + v_dot * self.step_size
                u_dot = (self.H[i,-1]/self.tau[i,-1]) * (W[i, -1] * self.Iext[i, timestep]) - 2 * self.u_t[i, -1]/self.tau[i,-1] - self.potential[i, -1, timestep]/(self.tau[i,-1]**2)
                self.u_t[i, -1] = self.u_t[i,-1] + u_dot * self.step_size

                # Add background input (to all populations)
                v_dot = self.u_t[i, -2]
                self.v_current[i, -2] = self.v_current[i, -2] + v_dot * self.step_size
                u_dot = (self.H[i,-2]/self.tau[i,-2]) * (W[i, -2] * self.Ib[i, timestep]) - 2 * self.u_t[i, -2]/self.tau[i,-2] - self.potential[i, -2, timestep]/(self.tau[i,-2]**2)
                self.u_t[i, -2] = self.u_t[i,-2] + u_dot * self.step_size

        print('finished loop...')

        return self.rate, self.potential
    


    





    
    






