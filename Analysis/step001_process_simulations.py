
"""
File: step001_process_simulations.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2025-10-21
Description: Compute the charasteristics of the simulated time series. 
- baseline activity (before stimulation)
- steady state/ "long term" activity (after stimulation)
    - immdediately after stimulation
    - after a longer period
- min/max firing rates and potentials (indicating oscillations)
"""
# %%
import os
import json
import pandas as pd
import numpy as np
from helper_functions import *

# Define paths
SIMDIR = os.getenv("SIMDIR")
sim_dir = os.path.join(SIMDIR, "simulation_results")
WDDIR = os.getenv("WDDIR")
figure_dir = os.path.join(SIMDIR, "Figures", "global_dynamics")
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

# define output directory
output_dir = os.path.join(SIMDIR, "derivatives")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# %%

# read params
params = load_parameters(WDDIR)
coupling_strengths = params['coupling_strengths'] # parse_coupling() coupling_strengths

# connectivity 
extI_cellcounts = params['extI_cellcounts']
bI_cellcounts = params['bI_cellcounts']
thal_cellcounts = params['thal_cellcounts']

# coupling strengths, balance and area selection
balance_EI = params['balance_EI']
g_thal = params['g_thal']
bEI_thal = params['bEI_thal']
step_size = params['step_size']
area = params['area']
filedir = params['filedir']

# inputs
input_type = params['input_type']
input_onset = params['input_onset']
simulation_dur = params['simulation_dur']
input_durations = params['input_durations']
input_strengths = params['input_strengths']
backgroundI_strengths = params['backgrndI_strengths']

# sampling parameters
step_size = 0.01 # of saving, not of simulation!
sample_delay_immediate = 0.3
sample_delay_late = 0.5 # when to start the long term behaviour "window"
input_onset = 1.001
sample_dur = 0.3
offset=0.1 # offset : time in s between baseline sampling and start of input 

# %%


# loop over all simulations
for d in input_durations:
    for s in input_strengths:
        for g in coupling_strengths:
            for bEI in balance_EI:
                for bI in backgroundI_strengths:

                    df = pd.DataFrame()
                    
                    # read data 
                    rates_df, potentials_df, filename = load_simulation_data(g, bEI, bI, d, s, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, sim_dir)

                    # compute characteristics
                    compute_longeterm_immediate(df, rates_df, potentials_df, input_onset, d, step_size, sample_delay_immediate, sample_dur)
                    compute_longeterm_late(df, rates_df, potentials_df, input_onset, d, step_size, sample_delay_late, sample_dur)
                    compute_input_response(df, rates_df, potentials_df, input_onset, d, step_size, sample_dur, offset)
                    compute_baseline(df, rates_df, potentials_df, input_onset, step_size, sample_dur, offset)
                    compare_longterm_baseline(df)             
                    classify_response(df)
                    set_sim_info(df, potentials_df, g, bEI, d, s, bI)

                    # save to csv
                    outfilename = filename.replace('.hdf5', '_processed.csv')
                    outpath = os.path.join(output_dir, outfilename)
                    df.to_csv(outpath, index=False)
                    

# %%
