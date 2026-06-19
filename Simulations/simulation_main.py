"""
File: simulation_main.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2025-08-05
Description: Run this file to run the simulation! 
"""

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
import mne
from mne.datasets import sample

location = "mpi"
if location == "laptop":
    WDDIR = r"C:\Users\gross\OneDrive - UvA\Documents\IMPRS_Leipzig\MyProject\Modelling\ChienReplication\SomatosensoryLaminarModel"
    SIMDIR = os.path.join(WDDIR, "output")
    DATADIR = "C:\\Users\\gross\\OneDrive - UvA\\Documents\\IMPRS_Leipzig\\MyProject\\Experiment\\Analysis\\LocalCode\\data"
    RECONDIR = os.path.join(DATADIR, 'freesurfer')

if location == "mpi":
    DATADIR = os.getenv('DATADIR')
    RECONDIR = os.getenv('SUBJECTS_DIR')
    SIMDIR = os.getenv("SIMDIR")
    WDDIR = os.getenv("WDDIR")
    RESDIR = os.getenv("RESDIR")
    
figure_dir = os.path.join(SIMDIR, "Figures")
tf_target_path = os.path.join(
    RESDIR, "Figures", "Main", "eeg_results", "source_reconstruction",
    "group", "_preprestim_corrected", "roi_epochswise",
    "group_roi_tf_morlet_ses-elec_preprestim_corrected.csv"
)
tc_target_path = os.path.join(
    RESDIR, "Figures", "Main", "eeg_results", "source_reconstruction",
    "group", "_preprestim_corrected", "roi_epochswise",
    "group_roi_timecourse_pooled_ses-elec_preprestim_corrected.csv"
)


# add model to datapath
sys.path.append(os.path.join(WDDIR, 'Simulations', 'model'))
from somato_model import SomatoModel, read_simulation_params
#from somato_model_pyrates_no_conn_operators_complete_pycobi import SomatoModelPyrates, read_simulation_params
import plotting_functions as pf

# %%   

# we parallelize over different coupling strengths (in srun HPC script)   
def parse_params():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--g",
        type=float,
        nargs="+",
        help="coupling strengths",
        required=False,
    )
    g = parser.g

    return g

# load EEG data and forward solution from real subject derivatives
subID    = 29           # example subject; edit to change
modality = 'elec'      # 'mecha' or 'elec'
l_freq_orig = 0.1
l_freq      = 1
h_freq      = 40
suffix      = '_preprestim_corrected'

sub_dir = os.path.join(DATADIR, 'derivatives', 'eeg-preproc',
                       f'sub-0{subID}', f'ses-{modality}')
epoch_fif = os.path.join(sub_dir,
    f"sub-0{subID}_ses-{modality}_task-NT_"
    f"lfreqori-{l_freq_orig}_lfreq-{l_freq}_hfreq-{h_freq}_epochs{suffix}.fif")
fwd_file  = os.path.join(sub_dir, f'sub-0{subID}_ico-5_ses-{modality}_fwd.fif')

epochs = mne.read_epochs(epoch_fif, preload=True)

fwd = mne.read_forward_solution(fwd_file)
fwd_fixed = mne.convert_forward_solution(
    fwd, surf_ori=True, force_fixed=True, use_cps=True
)

leadfield = fwd["sol"]["data"]
print(f"Leadfield size : {leadfield.shape[0]} sensors x {leadfield.shape[1]} dipoles")

src_free  = fwd["src"]
src_fixed = fwd_fixed["src"]

# needed for rh.BA3b.label lookup inside simulate_eeg()
data_path_labels = sample.data_path()

#%%
# Assign variables from loaded parameters
params = read_simulation_params()
input_onset = params['input_onset']
simulation_dur = params['simulation_dur']
save_params = params['save_params']
save_results = params['save_results']
save_full_potentials = params['save_full_potentials']
plot_rates = params['plot_rates']
plot_potentials = params['plot_potentials']
plot_all_potentials = params['plot_all_potentials']
jax_mode = params['jax_mode']

# specify output directory
filedir = params['filedir']
if not os.path.exists(filedir):
    os.makedirs(filedir)

# set parameters to loop over 
coupling_strengths = np.arange(0,55,5) #[100, 120, 140, 160]
backgrndI_strengths = np.arange(0,8,2) #[40, 60, 80] #,6,7]
input_durations = np.arange(0, 0.02, 0.004) # [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
input_strengths = np.arange(0,50,10)
strength_I = np.arange(0.2,0.9,0.02) #, 0.25, 0.26, 0.36]
ginters = np.arange(0,2,0.2)
resistance_factor = 1
area = 'all'
pyrates = False

for ginter in ginters:
    for d in input_durations:
        for sb in backgrndI_strengths:
            for s in input_strengths:
                # arrays to store simulation duration
                all_durations = []
                all_durations_saving = []

                for g in coupling_strengths:
                    for sI in strength_I:
                        
                        params['g_intercortical'] = ginter
                        params['coupling_strength'] = g 
                        params['strength_I'] = np.round(sI, 3)
                        params['Iext_duration'] = np.round(d, 3)
                        params['Iext_strength'] = s
                        params['Ib_strength'] = sb
                        params['area'] = area
                        params['resistance_factor'] = resistance_factor

                        # additional parameters (that are usually fixed)
                        params['g_thal'] = 2
                        params['g_thalPOm'] = 1
                        params['sI_thal'] = 0.5
                        params['delay_factor'] = 0.005
                        params['extI_cellcounts'] = 1000
                        params['bI_cellcounts'] = 100
                        # thalE_cellcounts / thalI_cellcounts / pom_cellcounts come from the JSON

                        if pyrates:
                            model = SomatoModelPyrates(params)
                        else:
                            model = SomatoModel(params)
                            model.plot_W_heatmap()
                        
                        # simulate rates and potentials
                        start = time.time()
                        model.simulate()
                        stop = time.time()
                        duration = stop - start
                        all_durations.append(duration)
                        print("Simulation duration (in s):", duration)

                        # analyse signal (frequency spectra)
                        #model.analyse_signal(save_spectrum=True)

                        """
                        # time-frequency error against measured data
                        tf_error, tf_sim, tf_target = model.compute_error_timefreq(tf_target_path)
                        print("TF error (log-MSE):", tf_error)
                        model.plot_timefreq_comparison(tf_sim, tf_target)

                        # time-course error against measured data
                        tc_error, tc_sim, tc_target = model.compute_error_timecourse(tc_target_path)
                        print("Time-course error (peak-norm MSE):", tc_error)
                        model.plot_timecourse_comparison(tc_sim, tc_target)

                        # persist per-run comparison maps/traces + values for later animation
                        run_dir = model.prepare_run_dir(filedir)
                        model.save_timefreq_comparison(
                            run_dir, tf_sim, tf_target, tf_error, filename="tf_comparison")
                        model.save_timecourse_comparison(
                            run_dir, tc_sim, tc_target, tc_error, filename="tc_comparison")
                        model.append_comparison_summary(tf_error=tf_error, tc_error=tc_error)

                        # compute dipoles
                        #sim_dip = model.compute_dipoles()
                        #model.plot_dipoles(sim_dip, epochs.info)
                        """

                        #error = compute_error_timefreq()

                        #stc, evoked, epochs_sim = model.simulate_eeg(epochs, data_path_labels, sim_dip, fwd, src_fixed)
                        #model.plot_eeg(evoked, epochs_sim)

                        # print important parameters
                        """
                        print('simulation_dur', model.simulation_dur)
                        print('step_size', model.step_size)
                        print('input_onset', model.input_onset) 
                        print('thal_connect', model.thal_connect) 
                        print('extI_cellcounts', model.extI_cellcounts) 
                        print('strength_I', model.strength_I) 
                        print('bI_cellcounts', model.bI_cellcounts) 
                        print('thalE_cellcounts', model.thalE_cellcounts)
                        print('thalI_cellcounts', model.thalI_cellcounts)
                        print('pom_cellcounts', model.pom_cellcounts)
                        print('sI_thal', model.sI_thal) 
                        print('g_thal', model.g_thal) 
                        print('input_type', model.input_type) 
                        print('area', model.area) 
                        """
                        print('coupling strength', model.coupling_strength) 
                        print('b input', model.Ib_strength) 
                        print('Iext strength', model.Iext_strength) 
                        print('Iext dur', model.Iext_duration) 

                        if save_results:
                            # create per-run output folder (holds params.json + all HDF5s for this run)
                            run_dir = model.prepare_run_dir(filedir)
                            start = time.time()
                            model.save_results_csv(run_dir, "results", save_full_potentials)
                            stop = time.time()
                            duration = stop - start
                            all_durations_saving.append(duration)
                            print("Saving duration (in s):", duration)

                        if plot_rates:
                            start_plot = 0
                            pf.plot_results(
                                model.rate,
                                model.Iext[-2],
                                model.Ib[0],
                                model.step_size,
                                simulation_dur,
                                start_plot,
                                sI,
                                g,
                                model.area,
                                d,
                                sb,
                                s,
                                figure_dir
                            )

                        if plot_potentials:
                            resolution_tstep = 1e-2
                            if pyrates:
                                potential_sum = model.potential
                            else:
                                potential_sum = np.sum(model.potential, axis=1)


                            pf.plot_potentials(
                                potential_sum,
                                model.Iext[-2],
                                model.Ib[0],
                                model.step_size,
                                simulation_dur,
                                start_plot,
                                figure_dir,
                                sI,
                                g,
                                d,
                                sb,
                                s
                            )

                        if plot_all_potentials:
                            pf.plot_all_potentials(
                                model.potential,
                                model.Iext[-2],
                                model.Ib[0],
                                model.step_size,
                                simulation_dur,
                                start_plot,
                            )

            # if (len(coupling_strengths) > 1):
            #    # used to plot with coupling strength on the x-axis and max/min rate on the y
            #    pf.plot_minmax(all_rates, coupling_strengths)

print("Mean Simulation duration: ", np.mean(all_durations))
print("Mean Saving duration: ", np.mean(all_durations_saving))


# %%
