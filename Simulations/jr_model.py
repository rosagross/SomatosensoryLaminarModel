import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import yaml
import json
from parameters import Parameter
import os

WDDIR = os.getenv("WDDIR")

class JR_Model():

    def __init__(self, g, sb, s, d):
        
        # load in all connectivity parameters
        self.p = Parameter()

        # Simulation parameters
        self.tau = self.p.tau
        self.nPop = self.p.nPop
        # sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
        self.sigm = self.p.sigmoid_params

        # parameters that will be updated from the json file 
        # (first initialized with default values) 
        self.simulation_time = 2 # in s
        self.step_size = 0.001 # in s
        self.input_onset = 1.001
        self.thal_connect = [0,0,0,0]
        self.extI_cellcount = 1000
        self.bI_cellcount = 100
        self.thal_cellcount = 500
        self.bEI_thal = 0.5
        self.g_thal = 2
        self.input_type = 'step'
        self.area = 'all' 

        # parameters that come as function parameters since we might loop over it outside of the class
        self.g = g
        self.sb = sb
        self.s = s
        self.d = d

        # external input matrix and background input matrix
        # update parameters based on file
        params = self.read_simulation_params()
        self.__dict__.update(params)

        # create input array
        Iext = self.create_Iext(self.simulation_time, self.step_size, self.input_onset, self.d, self.s, self.input_type)
        Ib = self.create_Ibackground(simulation_time, self.step_size, self.sb)
        self.gE = self.g * self.bEI 
        self.gI = self.g * (1 - self.bEI)
        self.gE_thal = self.g_thal * self.bEI_thal
        self.gI_thal = self.g_thal * (1 - self.bEI_thal)

        # Synaptic kernel 
        self.H = np.ones((self.nPop, self.nPop+1))

        # define time steps 
        self.steps = np.arange(self.step_size, self.simulation_time+self.step_size, self.step_size)

        # extend input arrays
        self.Iext = np.tile(Iext, (self.nPop,1))
        self.Ib = np.tile(Ib, (self.nPop,1))


    def create_Iext(self, imulation_time, step_size, input_onset, input_duration, input_strength, input_type):
        """Creates external input."""

        Iext = np.zeros(int(simulation_time / step_size))

        if input_type == "step":
            t = int(input_duration / step_size)
            t0 = int(input_onset / step_size)
            Iext[t0 : t0 + t] = input_strength
        elif input_type == "background":
            # provide input for the entire simulation duration
            Iext[:] = input_strength

        return Iext


    def create_Ibackground(self, simulation_time, step_size, input_strength):
        """Create Background Input"""
        Ib = np.zeros(int(simulation_time / step_size))
        Ib[:] = input_strength
        return Ib


    def read_simulation_params(self):
        """Read simulation parameters from json file."""
        # Read in preprocessing parameters
        with open(os.path.join(WDDIR, 'Simulations', 'simulation_parameter.json'), 'r') as json_file:
            params = json.load(json_file)
        
        return params

        
    def save_to_yaml(self, filename):
        
        S = self.p.get_connectStrength()
        P = self.p.get_connectProb()
        C = self.p.get_cellcounts()
        W = self.p.get_connectivity(self.gE, self.gI, self.gEthal, self.gIthal, self.thal_connect)

        # Convert numpy arrays to lists
        parameters = {
            'gE': self.gE,
            'gI': self.gI,
            'gEthal': self.gEthal, 
            'gIthal': self.gIthal,
            'S': S.tolist(),
            'P': P.tolist(),
            'C': C.tolist(),
            'W': W.tolist()
        }

        # Save parameters to a YAML file
        with open(filename + '.yaml', 'w') as file:
            yaml.dump(parameters, file)


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


    def save_results_csv(self, filedir, filename, full=False):
        """
        Safe the simulated data in a csv file
        """

        population_names = [
            "E3b",
            "PV3b",
            "SST3b",
            "VIP3b",
            "E1",
            "PV1",
            "SST1",
            "VIP1",
            "E2",
            "PV2",
            "SST2",
            "E3",
            "PV3",
            "SST3",
            "E4",
            "PV4",
            "SST4",
            "E1S2",
            "PV1S2",
            "SST1S2",
            "VIP1S2",
            "E2S2",
            "PV2S2",
            "SST2S2",
            "E3S2",
            "PV3S2",
            "SST3S2",
            "E4S2",
            "PV4S2",
            "SST4S2",
        ]
        cells = np.concatenate((population_names, ["ThalE", "ThalI"]))

        filename = filename + ".hdf5"

        # only safe every second datapoint
        resolution_tstep = 0.01
        print("tstep resolution", resolution_tstep)
        rates_downsampled = self.rates[:, :: int(1000 * resolution_tstep)]
        rates_df = pd.DataFrame(rates_downsampled.T, columns=cells)
        rates_df.to_hdf(
            os.path.join(filedir, filename), index=False, key="rates", mode="a"
        )

        # sum the potentials together and save them
        potential_sum = np.sum(self.potentials, axis=1)
        potential_sum_downsampled = potential_sum[:, :: int(1000 * resolution_tstep)]
        potential_df = pd.DataFrame(potential_sum_downsampled.T, columns=cells)
        potential_df.to_hdf(
            os.path.join(filedir, filename), index=False, key="summed_potential", mode="a"
        )

        if full:
            # save all potentials additionally
            psp_filename = "full_" + filename
            print('full potential file:', psp_filename)
            write_3D_csv(os.path.join(filedir, psp_filename), self.potentials)



    def write_3D_csv(filename, data):
        """
        Write results in form of a 3D hdf5 file.
        """
        dataset_name = 'full_potentials'

        with h5py.File(filename, "w") as f:
            f.create_dataset(dataset_name, data=data, compression="gzip")

        





    
    






