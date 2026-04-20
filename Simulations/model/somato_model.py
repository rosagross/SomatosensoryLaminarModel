import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import h5py
import yaml
import json
from parameters import Parameter
import os

WDDIR = os.getenv("WDDIR")
SIMDIR = os.getenv("SIMDIR")

def read_simulation_params():
    """Read simulation parameters from json file."""
    # Read in preprocessing parameters
    with open(os.path.join(WDDIR, 'Simulations', 'simulation_parameter.json'), 'r') as json_file:
        params = json.load(json_file)
    
    return params

class SomatoModel():

    def __init__(self, params={}):
        
        # load in all connectivity parameters, time constants, etc.
        self.p = Parameter()
        self.tau = self.p.tau
        self.nPop = self.p.nPop
        # sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
        self.sigm = self.p.sigmoid_params

        # parameters that will be updated from the json file 
        # (first initialized with default values) 
        self.simulation_dur = 2 # in s
        self.step_size = 0.001 # in s
        self.resolution_tstep = 0.001 # in s
        self.input_onset = 1.001
        self.thal_connect = [0,0,0,0]
        self.extI_cellcounts = 1000
        self.strength_I = 0 #0.7
        self.bI_cellcounts = 100
        self.thal_cellcounts = 500
        self.sI_thal = 0.5
        self.g_thal = 2
        self.input_type = 'step'
        self.area = 'all' 
        self.coupling_strength = 10
        self.Ib_strength = 7
        self.Iext_strength = 10
        self.Iext_duration = 0.5

        # update parameters based on params dicts
        self.__dict__.update(params)

        # create input array
        Iext = self.create_Iext()
        Ib = self.create_Ibackground()
        self.gE = self.coupling_strength
        self.gI = self.coupling_strength * self.strength_I
        self.gEthal = self.g_thal 
        self.gIthal = self.g_thal * self.sI_thal

        # Synaptic kernel 
        self.H = np.ones((self.nPop, self.nPop+1))

        # define time steps 
        self.steps = np.arange(self.step_size, self.simulation_dur+self.step_size, self.step_size)

        # extend input arrays
        self.Iext = np.tile(Iext, (self.nPop,1))
        self.Ib = np.tile(Ib, (self.nPop,1))

        self.filename = (
            f"gthal{self.g_thal}_sIthal{self.sI_thal}_g{self.coupling_strength}_sI{self.strength_I}_Ib{self.Ib_strength}_Iextd{self.Iext_duration}_"
            f"{self.input_type}Iexts{self.Iext_strength}_Ionset{self.input_onset}_thalcells{self.thal_cellcounts}_"
            f"Ibcells{self.bI_cellcounts}_Iextcells{self.extI_cellcounts}_thalUncon_S1S2Connected"
        )

        # Output matrices to store computed values for rates & potentials (E, IIN , EIN) 
        self.rate = np.zeros((self.nPop, len(self.steps)))
        self.potential = np.zeros((self.nPop, self.nPop+2, len(self.steps))) 

        # Simulation loop
        # Initialize first values for the potential, rate and first order derivative with 0 or randomly
        self.v_current = np.zeros((self.nPop, self.nPop+2)) # +2 because 1 for background input and one for external input 
        self.rate_current = np.zeros(self.nPop)
        self.u_t = np.zeros((self.nPop, self.nPop+2)) # the initial first-order derivative: v'(t) = u(t)
        self.t = 0.0

        # Weight matrix [to x from]
        self.W = self.p.get_connectivity(self.gE, self.gI, self.gEthal, self.gIthal, self.thal_connect, self.extI_cellcounts, self.bI_cellcounts, self.thal_cellcounts, area=self.area) 


    def initialize_state(self):
        """
        Reset the dynamic state for interactive simulations.
        """
        self.v_current = np.zeros((self.nPop, self.nPop + 2))
        self.rate_current = np.zeros(self.nPop)
        self.u_t = np.zeros((self.nPop, self.nPop + 2))
        self.t = 0.0
        self.rate = np.zeros((self.nPop, len(self.steps)))
        self.potential = np.zeros((self.nPop, self.nPop + 2, len(self.steps)))


    def apply_params(self, params: dict):
        """
        Update parameters and recompute derived state.
        """
        self.__dict__.update(params)

        # recompute inputs and gains
        Iext = self.create_Iext()
        Ib = self.create_Ibackground()
        self.Iext = np.tile(Iext, (self.nPop, 1))
        self.Ib = np.tile(Ib, (self.nPop, 1))

        self.gE = self.coupling_strength 
        self.gI = self.coupling_strength * self.strength_I
        self.gEthal = self.g_thal 
        self.gIthal = self.g_thal * self.sI_thal

        # update connectivity with new gains and counts
        self.W = self.p.get_connectivity(
            self.gE,
            self.gI,
            self.gEthal,
            self.gIthal,
            self.thal_connect,
            self.extI_cellcounts,
            self.bI_cellcounts,
            self.thal_cellcounts,
            area=self.area
        )


    def create_Iext(self):
        """Creates external input."""

        Iext = np.zeros(int(self.simulation_dur / self.step_size))

        if self.input_type == "step":
            t = int(self.Iext_duration / self.step_size)
            t0 = int(self.input_onset / self.step_size)
            Iext[t0 : t0 + t] = self.Iext_strength
        elif self.input_type == "background":
            # provide input for the entire simulation duration
            Iext[:] = self.Iext_strength

        return Iext


    def create_Ibackground(self):
        """Create Background Input"""
        Ib = np.zeros(int(self.simulation_dur / self.step_size))
        Ib[:] = self.Ib_strength
        return Ib

        
    def save_to_yaml(self, filename):
        
        S = self.p.get_connectStrength()
        P = self.p.get_connectProb()
        C = self.p.get_cellcounts()
        W = self.p.get_connectivity(self.gE, self.gI, self.gEthal, self.gIthal, self.thal_connect, self.extI_cellcounts, self.bI_cellcounts, self.thal_cellcounts)

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

    def plot_W_heatmap(self):
        """
        Plot connectivity matrix as heatmap.
        """
        
        W = self.p.get_connectivity(self.gE, self.gI, self.gEthal, self.gIthal, self.thal_connect, self.extI_cellcounts, self.bI_cellcounts, self.thal_cellcounts, area=self.area) 
        sns.heatmap(W, annot=False, cmap='coolwarm', center=0, xticklabels=True, yticklabels=True)


    def simulate(self):
        '''
        Simulation loop
        '''

        last_step = self.steps[-1] 

        for timestep, _ in enumerate(self.steps):

            # compute simulation step         
            self.potential[:, :, timestep] = self.v_current.copy()
            self.rate_current = self.compute_rates()
            self.rate[:, timestep] = self.rate_current
            self.v_current = self.compute_potentials(timestep)

        print('finished loop...')

    def compute_rates(self):
        """
        Compute the firing rates of all populations.
        """
        v_sum = np.sum(self.v_current, axis=1)

        self.rate_current = self.sigm[:,2] / (
            1 + np.exp(self.sigm[:,0] * (self.sigm[:,1] - v_sum))
        )

        return self.rate_current.copy()

    def compute_potentials(self, timestep):
        """
        Compute the potentials of all populations.
        Also take into account the background input and external input.
        """

        pop_slice = slice(0, self.nPop)

        # store previous potentials (needed for correct Euler update)
        v_prev = self.v_current.copy()

        # -----------------------------
        # POPULATION INTERACTIONS
        # -----------------------------

        # update potentials
        self.v_current[:, pop_slice] += (
            self.u_t[:, pop_slice] * self.step_size
        )

        # synaptic drive
        drive = self.W[:, pop_slice] * self.rate_current

        u_dot = (
            (self.H[:, pop_slice] / self.tau[:, pop_slice]) * drive
            - 2 * self.u_t[:, pop_slice] / self.tau[:, pop_slice]
            - v_prev[:, pop_slice] / (self.tau[:, pop_slice] ** 2)
        )

        self.u_t[:, pop_slice] += u_dot * self.step_size


        # -----------------------------
        # EXTERNAL INPUT
        # -----------------------------

        # update potentials
        v_dot = self.u_t[:, -1]
        self.v_current[:, -1] += v_dot * self.step_size

        u_dot = (
            (self.H[:, -1] / self.tau[:, -1])
            * (self.W[:, -1] * self.Iext[:, timestep])
            - 2 * self.u_t[:, -1] / self.tau[:, -1]
            - v_prev[:, -1] / (self.tau[:, -1] ** 2)
        )

        self.u_t[:, -1] += u_dot * self.step_size


        # -----------------------------
        # BACKGROUND INPUT
        # -----------------------------

        # update potentials
        v_dot = self.u_t[:, -2]
        self.v_current[:, -2] += v_dot * self.step_size

        u_dot = (
            (self.H[:, -2] / self.tau[:, -2])
            * (self.W[:, -2] * self.Ib[:, timestep])
            - 2 * self.u_t[:, -2] / self.tau[:, -2]
            - v_prev[:, -2] / (self.tau[:, -2] ** 2)
        )

        self.u_t[:, -2] += u_dot * self.step_size

        return self.v_current.copy()


    def simulate_step(self, timestep):
        """
        Only for interactive computation
        """
        self.rate_current = self.compute_rates()
        self.v_current = self.compute_potentials(timestep)

        self.t += self.step_size

        return self.rate_current.copy()


    def compute_ecds():
        raise NotImplementedError

    def save_results_csv(self, filedir, filename, full=False, save_params=False):
        """
        Safe the simulated data in a csv file
        """

        cells = [
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
            "ThalE", 
            "ThalI"
        ]
        cells = np.array(cells)

        filename = filename + ".hdf5"

        # only safe every X datapoint
        print("tstep resolution", self.resolution_tstep)
        rates_downsampled = self.rate[:, :: int(1000 * self.resolution_tstep)]
        rates_df = pd.DataFrame(rates_downsampled.T, columns=cells)
        rates_df.to_hdf(
            os.path.join(filedir, filename), index=False, key="rates", mode="a"
        )

        # sum the potentials together and save them
        potential_sum = np.sum(self.potential, axis=1)
        potential_sum_downsampled = potential_sum[:, :: int(1000 * self.resolution_tstep)]
        potential_df = pd.DataFrame(potential_sum_downsampled.T, columns=cells)
        potential_df.to_hdf(
            os.path.join(filedir, filename), index=False, key="summed_potential", mode="a"
        )

        if full:
            # save all potentials additionally
            psp_filename = "full_" + filename
            print('full potential file:', psp_filename)
            self.write_3D_csv(os.path.join(filedir, psp_filename))

        if save_params:
            # safe connectivty parameter in yaml file
            self.save_to_yaml(os.path.join(filedir, "params" + self.filename))


    def write_3D_csv(self, filename):
        """
        Write results in form of a 3D hdf5 file.
        """
        dataset_name = 'full_potentials'

        with h5py.File(filename, "w") as f:
            f.create_dataset(dataset_name, data=self.potential, compression="gzip")
        





    
    



