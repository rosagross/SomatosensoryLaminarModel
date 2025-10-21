
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
import pandas as pd
import numpy as np
from helper_functions import *

# Define paths
SIMDIR = os.getenv("SIMDIR")
WDDIR = os.getenv("WDDIR")
figure_dir = os.path.join(SIMDIR, "Figures", "global_dynamics")
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

# define output directory
output_dir = os.path.join(SIMDIR, "derivatives")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# %%
# define parameters for looping over simulations
input_durations = [0.5, 1, 1.5] 
input_strengths = [10, 20, 30] 
coupling_strengths = [10, 20, 30]  # global coupling strengths
balance_EI = [0.7, 0.8, 0.9]  # balance of excitation and inhibition
backgroundI_strengths = [5, 6, 7]  # background input to inhibitory populations
input_type = 'step'
input_onset = 1.001
thal_cellcounts = 500
extI_cellcounts = 1000
bI_cellcounts = 100
step_size = 0.01
sample_delay_immediate = 0.3
sample_delay_late = 1 # when to start the long term behaviour "window"
input_onset = 1.001
sample_dur = 0.3
offset=0.1 # offset : time in s between baseline sampling and start of input 

# loop over all simulations
for d in input_durations:
    for s in input_strengths:
        for g in coupling_strengths:
            for bEI in balance_EI:
                for bI in backgroundI_strengths:

                    df = pd.DataFrame()
                    
                    # read in firing rates in data matrix (datapoints x populations)
                    filename = f"g{g}_bEI{bEI}_Ib{bI}_Iextd{d}_{input_type}Iexts{s}_Ionset{input_onset}_thalcells{thal_cellcounts}_Ibcells{bI_cellcounts}_Iextcells{extI_cellcounts}_thalUncon_S1S2Uncon.hdf5"
                    rates_df = pd.read_hdf(os.path.join(output_dir, filename), key='rates')
                    
                    # read in potentials in data matrix (datapoints x populations)
                    potentials_df = pd.read_hdf(os.path.join(output_dir, filename), key='summed_potential')

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
                    
