
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
import glob
import pandas as pd
import numpy as np
from helper_functions import *

# Define paths
SIMDIR = "/data/pt_02989"
sim_dir = os.path.join(SIMDIR, "output_grossmannr")
WDDIR = os.getenv("WDDIR")
figure_dir = os.path.join(SIMDIR, "Figures", "global_dynamics")
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

# %%

# read params
params = load_parameters(WDDIR)

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

# Which runs to process: every per-run folder in sim_dir that holds a results.hdf5 and whose
# name does not end in _thalUncon (the older, thalamus-unconnected batch). Each folder's own
# params.json supplies its parameters, so no filename stem has to be reconstructed here and
# batches with different parameter grids are handled in one pass.
run_dirs = sorted(
    d for d in glob.glob(os.path.join(sim_dir, "*"))
    if os.path.isdir(d)
    and not os.path.basename(d).endswith("_thalUncon")
    and os.path.exists(os.path.join(d, "results.hdf5"))
)
print(f"found {len(run_dirs)} runs to process in {sim_dir}")

n_done, failed = 0, []
for i, run_dir in enumerate(run_dirs, start=1):
    try:
        with open(os.path.join(run_dir, "params.json"), 'r') as f:
            run_params = json.load(f)

        # simulation parameters of this run (not from analysis_parameter.json, which
        # describes a different setup)
        g = run_params['coupling_strength']
        sI = run_params['strength_I']
        d = run_params['Iext_duration']
        s = run_params['Iext_strength']
        bI = run_params['Ib_strength']
        g_thalPOm = run_params['g_thalPOm']
        g_inter = run_params['g_intercortical']
        Ib_noise_std = run_params['Ib_noise_std']
        input_onset = run_params['input_onset']
        step_size = run_params['resolution_tstep']  # of saving, not of simulation!

        df = pd.DataFrame()

        # read data
        results_path = os.path.join(run_dir, "results.hdf5")
        rates_df = pd.read_hdf(results_path, key='rates')
        potentials_df = load_hdf_safe(results_path)

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

        set_sim_info(df, potentials_df, g, sI, d, s, bI,
                     g_thalPOm=g_thalPOm, g_intercortical=g_inter, Ib_noise_std=Ib_noise_std)

        # save to csv inside the per-run folder, next to results.hdf5
        outpath = os.path.join(run_dir, "processed.csv")
        df.to_csv(outpath, index=False)
        n_done += 1

    except Exception as err:
        # one unreadable run should not abort the batch
        failed.append((os.path.basename(run_dir), err))
        print(f"  FAILED {os.path.basename(run_dir)}: {err}")

    if i % 25 == 0 or i == len(run_dirs):
        print(f"  {i}/{len(run_dirs)} processed")

print(f"done: {n_done} processed, {len(failed)} failed")
for name, err in failed:
    print(f"  {name}: {err}")


# %%
