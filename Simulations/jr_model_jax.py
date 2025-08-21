# %%
import numpy as np
import matplotlib.pyplot as plt
from parameters import Parameter
import os
import jax.numpy as jnp
from jax import jit, lax, vmap
from jax_simulation import run_jax_simulation  # <- contains the JAX version
from parameters import Parameter
import os

# %%
class JR_Model():

    def __init__(self, Iext, Ib, gE, gI, gEthal, gIthal, thal_connect, filedir, filename, step_size=0.001, simulation_time=1):
        
        # load in all parameters
        self.p = Parameter()
        # Simulation parameters
        self.tau = jnp.array(self.p.tau)
        self.nPop = self.p.nPop
        self.simulation_time = simulation_time # in s
        self.step_size = step_size # in s
        self.gE = gE
        self.gI = gI
        self.gEthal = gEthal
        self.gIthal = gIthal
        self.thal_connect = jnp.array(thal_connect)

        # sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
        self.sigm = jnp.array(self.p.sigmoid_params)

        # Synaptic kernel 
        self.H = jnp.ones((self.nPop, self.nPop+1))

        # define time steps 
        self.steps = jnp.arange(self.step_size, self.simulation_time+self.step_size, self.step_size)

        # external input matrix and background input matrix
        self.Iext = jnp.tile(Iext, (self.nPop,1))
        self.Ib = jnp.tile(Ib, (self.nPop,1))

        # Output matrices to store computed values for rates & potentials (E, IIN , EIN) 
        self.rate = jnp.zeros((self.nPop, len(self.steps)))
        self.potential = jnp.zeros((self.nPop, self.nPop+2, len(self.steps))) 

        # Simulation loop
        # Initialize first values for the potential, rate and first order derivative with 0 or randomly
        self.v_current = jnp.zeros((self.nPop, self.nPop+2)) # +2 because 1 for background input and one for external input 
        self.rate_current = jnp.zeros(self.nPop)
        self.u_t = jnp.zeros((self.nPop, self.nPop+2)) # the initial first-order derivative: v'(t) = u(t)

        # Weight matrix [to x from]
        self.W = jnp.array(self.p.get_connectivity(self.gE, self.gI, self.gEthal, self.gIthal, self.thal_connect))

        
    def run_simulation(self):
        '''
        Simulation loop
        '''

        # parameter dictionary
        params = {
        'nPop': self.nPop,
        'sigm': self.sigm
        }

        # call JAX simulation function
        rate, potential = run_jax_simulation(
            nPop=self.nPop,
            sigm=self.sigm,
            W=self.W,
            H=self.H,
            tau=self.tau,
            Iext=self.Iext,
            Ib=self.Ib,
            step_size=self.step_size,
            steps=self.steps
        )

        return rate, potential
    


    





    
    






