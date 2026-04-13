
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
strength_I = params['strength_I']
g_thal = params['g_thal']
sI_thal = params['sI_thal']
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

#coupling_strengths = [10] #np.arange(0,55,5) #[100, 120, 140, 160]
#backgroundI_strengths = [0] #np.arange(0,8,2) #[40, 60, 80] #,6,7]
#input_durations = [0.004] #np.arange(0, 0.02, 0.002) # [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
#input_strengths = [10] #np.arange(0,50,10)
#strength_I = [0.2]


# sampling parameters
sampling_params = params['sampling']
step_size = sampling_params['step_size']  # of saving, not of simulation!
sample_delay_immediate = sampling_params['sample_delay_immediate']
sample_delay_late = sampling_params['sample_delay_late']  # when to start the long term behaviour "window"
sample_dur = sampling_params['sample_dur']
offset = sampling_params['offset']  # time in s between baseline sampling and start of input
rate_osc_threshold = sampling_params['rate_osc_threshold']
potential_osc_threshold = sampling_params['potential_osc_threshold']
response_window = sampling_params['response_window']

# %%

# loop over all simulations
for d in input_durations:
    print('Iextd', d)
    for s in input_strengths:
        print('Iexts', s)
        for g in coupling_strengths:
            #print('coupling', g)
            for sI in strength_I:
                #print('sI', sI)
                for bI in backgroundI_strengths:
                    #print('bI', bI)

                    df = pd.DataFrame()
                    
                    # read data 
                    rates_df, potentials_df, filename = load_simulation_data(g, sI, bI, d, s, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, sim_dir)

                    # compute characteristics
                    compute_longeterm_immediate(df, rates_df, potentials_df, input_onset, d, step_size, sample_delay_immediate, sample_dur)
                    compute_longeterm_late(df, rates_df, potentials_df, input_onset, d, step_size, sample_delay_late, sample_dur)
                    compute_input_response(df, rates_df, potentials_df, input_onset, d, step_size, response_window)
                    compute_baseline(df, rates_df, potentials_df, input_onset, step_size, sample_dur, offset)
                    compare_longterm_baseline(df)             
                    classify_response(df)
                    
                    # compute oscillation frequency in different windows
                    baseline_start = int((input_onset - (sample_dur+offset))/step_size) 
                    baseline_stop = int(baseline_start + sample_dur/step_size)
                    compute_window_frequency(
                        df, rates_df, potentials_df,
                        baseline_start, baseline_stop,
                        "baseline", step_size,
                        rate_osc_threshold, potential_osc_threshold
                    )

                    start_sample_during = int((input_onset)/step_size)
                    stop_sample_during = int((input_onset+d)/step_size)
                    compute_window_frequency(
                        df, rates_df, potentials_df,
                        start_sample_during, stop_sample_during,
                        "duringInput", step_size,
                        rate_osc_threshold, potential_osc_threshold
                    )

                    start_sample_late = int((input_onset+d+sample_delay_late)/step_size)
                    stop_sample_late = int(start_sample_late + sample_dur/step_size)
                    compute_window_frequency(
                        df, rates_df, potentials_df,
                        start_sample_late, stop_sample_late,
                        "lateLongterm", step_size,
                        rate_osc_threshold, potential_osc_threshold
                    )

                    set_sim_info(df, potentials_df, g, sI, d, s, bI)

                    # save to csv
                    outfilename = filename.replace('.hdf5', '_processed.csv')
                    outpath = os.path.join(output_dir, outfilename)
                    df.to_csv(outpath, index=False)
                    
