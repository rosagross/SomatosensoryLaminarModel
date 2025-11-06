# %%
import numpy as np
import h5py
import os
import sys
import json
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import time
import csv
from jr_model import JR_Model

optimizer_path = "/data/p_02989/Modelling/neuronaldynamics"
if optimizer_path not in sys.path:
    sys.path.append(optimizer_path)
from Optimizers.Optimizer import Hierarchical_Random, GA

# Add the parent directory to the sys.path
sim_path = "/data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel/EEGSimulation"
if sim_path not in sys.path:
    sys.path.append(sim_path)
from sim_meg import *

# %%
SIMDIR = os.getenv("SIMDIR")
WDDIR = os.getenv("WDDIR")
figure_dir = os.path.join(SIMDIR, "Figures")
# %%




class SimulationSomato():
    """
    Simulation instance for running somatosensory model simulation.
    """
    def __init__(self, parameters={}):
        self.coupling_strengths = np.array([[0]])
        self.balance_EI = np.array([[0]])
        self.step_size = 0.001
        self.simulation_dur = 2
        self.input_onset = 1.001
        self.t = None
        self.__dict__.update(parameters)
        self.simulation_time = int(self.input_onset) + self.simulation_dur
        self.t = np.arange(self.step_size, self.simulation_time+self.step_size, self.step_size)
    
    def simulate(self):

        # PARAMETER
        params = read_simulation_params()

        # coupling strengths, balance and area selection
        g_thal = params['g_thal']
        bEI_thal = params['bEI_thal']
        step_size = params['step_size']
        area = params['area']
        filedir = params['filedir']

        # inputs
        input_type = params['input_type']
        input_onset = params['input_onset']
        d = params['input_durations'][0]
        s = params['input_strengths'][0]
        sb = params['backgrndI_strengths'][0]

        # connectivity 
        thal_connect = np.array(params['thal_connect'])
        extI_cellcounts = params['extI_cellcounts']
        bI_cellcounts = params['bI_cellcounts']
        thal_cellcounts = params['thal_cellcounts']

        filename = f"g{self.coupling_strengths}_bEI{self.balance_EI}_Ib{sb}_Iextd{d}_{input_type}Iexts{s}_Ionset{input_onset}_thalcells{thal_cellcounts}_Ibcells{bI_cellcounts}_Iextcells{extI_cellcounts}_thalUncon_S1S2Uncon"

                
        Iext = create_Iext(
            self.simulation_time, step_size, input_onset, d, s, input_type
        )
        Ib = create_Ibackground(self.simulation_time, step_size, sb)
        gE = self.coupling_strengths * self.balance_EI 
        gI = self.coupling_strengths * (1 - self.balance_EI)
        gE_thal = g_thal * bEI_thal
        gI_thal = g_thal * (1 - bEI_thal)
        # for now we use the same coupling strength for the thalamus connections as for the cortical connections
        coupling_thalE = gE_thal
        coupling_thalI = gI_thal
        print("gE", gE)
        print("gI", gI)
        thal_connect_scaled = thal_connect 

        # init model 
        model = JR_Model(
                            Iext,
                            Ib,
                            gE,
                            gI,
                            coupling_thalE,
                            coupling_thalI,
                            thal_connect_scaled,
                            extI_cellcounts,
                            bI_cellcounts,
                            thal_cellcounts,
                            step_size,
                            self.simulation_time,
                            area=area,
                        )

        # simulation result (timecourse)
        start_sim = time.time()
        self.rate, self.potentials = model.run_simulation()
        end_sim = time.time()
        print('simulation duration', end_sim-start_sim)

        # observer model 
        dipole_setting = [-1, -1, 1, 1, 1]
        nCells = [13, 4, 3, 3, 3]
        indices_E = [4, 8, 11, 15]

        ecd_time = time.time()
        dipoles = simDipoles(dipole_setting, nCells, indices_E, psp=self.potentials)
        self.ecd = np.sum(dipoles, axis=0)
        ecd_stop_time = time.time()
        print('ECD computation time', ecd_stop_time-ecd_time)

# %%

optim_params = {'coupling_strengths':10, 'balance_EI':0.7, 'input_strengths':10}
simEEG = SimulationSomato(optim_params)
simEEG.simulate()
# %%
# optmizer setup
sim_object = SimulationSomato()
model_parameters = ['coupling_strengths', 'balance_EI', 'input_strengths']
opt_params = {}
opt_params['model_parameters'] = model_parameters
opt_params['bounds'] = [[5, 50], [0.4, 1], [0, 50]]
opt_params['n_iter'] = 5
opt_params['N1'] = 5 # population size (N2 is the population of the mutation crossover populations)
opt_params['verbose'] = 1
opt_params['x_out'] = 'ecd'
opt_params['simulation_class'] = sim_object
opt_params['simulate'] = sim_object.simulate
opt_params['reference'] = np.array(pd.read_csv(os.path.join(WDDIR, 'Optimization', 'optimization_reference_same.csv')))

optimizer = GA(parameters=opt_params)
optimizer.run()
optimal_param = optimizer.optimum

# %%
opt_somato = SimulationSomato({'coupling_strengths': optimal_param[0], 'balance_EI': optimal_param[1], 'input_strengths': optimal_param[2]})
opt_somato.simulate()
y_opt = opt_somato.ecd
print(f'y_opt = {y_opt}')

# %% 
optimizer.plot_fit()
# plotting the optimization results
# x-axis : coupling strength
# y-axis : balance EI 
# value : error 

