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
ps_target_path = os.path.join(
    RESDIR, "Figures", "Main", "eeg_results", "source_reconstruction",
    "group", "_preprestim_corrected", "roi_epochswise",
    "group_roi_prestim_spectrum_ses-elec_preprestim_corrected.csv"
)


# add model to datapath
sys.path.append(os.path.join(WDDIR, 'Simulations', 'model'))
from somato_model import SomatoModel, read_simulation_params, load_optimized_params
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

# electrical-modality subjects; compute_dipoles reads each one's forward model
# (same list as Analysis/SourceReconstruction/step002_inverse_solution_multisub_epochswise.py)
subID_elec = [15] #, 16, 17, 18, 23, 24, 25, 26, 27, 28, 29, 34, 35, 36, 37, 38, 39, 40,
              #42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52]

# needed for rh.BA3b.label lookup inside simulate_eeg()
data_path_labels = sample.data_path()

#%%
# Assign variables from loaded parameters
params = read_simulation_params()
input_onset = params['input_onset']
simulation_dur = params['simulation_dur']
save_params = params['save_params']
save_results = params['save_results']
save_connectivity = params['save_connectivity']
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
# Values below are the best fit of the GA run opt_20260803_124557_tc_roi-A1
# (error mode "tc", ROI A1, best combined error 0.008113): best_params from its
# optimization_summary.json plus the fixed base_params from its run_config.json.
# They are overrides on top of the simulation_parameter.json defaults loaded above,
# and are only applied in the `if not use_opt:` branch below - with use_opt = True
# the whole block is unused and the parameters come from load_optimized_params().
coupling_strengths = [3.293441306598355]  # np.arange(0,20,1)#[100, 120, 140, 160]
backgrndI_strengths = [9.994769237606103] # np.arange(4,10,1)#[40, 60, 80] #,6,7]
input_durations = [0.0291] # np.arange(0, 0.02, 0.004)# [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
input_strengths = [20] # np.arange(0,50,10)
strength_Is = [0.5950270376909454] # np.arange(0.68,0.76,0.005)#, 0.25, 0.26, 0.36]
ginters = [2.0] #np.arange(0,2,0.2)
g_thalPOms = [1.258205044608075] #np.arange(0.1,1.01,0.1) # scales the POm population's output connectivity
delay_factor = 0.0017424379080450052
receptor_thalamus_delay = 0.01127 # from periphery to thal
thal_delay_factor = 0.004355885296978261 # from thal to cortex delay
e3b_tau = 9.98520373331793
e1_tau = 8.680392487165907
e2_tau = 4.290026558869845
p_2PVE = 16.351739686255572 # L4 PV <- E connection probability (S1 and S2); original value 37.4
p_4PVE = 15.198587424769887 # L6 PV <- E connection probability (S1 and S2); original value 39.6
resistance_factor = 1
area = 'all'
pyrates = False
use_opt = True
save_connectivity = True

Ib_noise_std = 0


if use_opt:
    # option to read file from optimization run: simulation_parameter.json base params
    # updated with the run's best_params (keeps input_onset, cell counts, filedir, ...)
    # EEGSimulation/plot_dipole_computation.py plots this same configuration - update
    # its opt_run too when changing it here.
    opt_run = "opt_20260804_090235_tc_roi-S2" #"opt_20260803_124557_tc_roi-A1"
    params = load_optimized_params(opt_run, overrides={'Ib_noise_std': Ib_noise_std})


    
for ginter in ginters:
    print("ginter", np.round(ginter, 3))
    for g_thalPOm in g_thalPOms:
      print("g_thalPOm", np.round(g_thalPOm, 3))
      for d in input_durations:
        for sb in backgrndI_strengths:
            print("background", sb)
            for s in input_strengths:
                # arrays to store simulation duration
                all_durations = []
                all_durations_saving = []

                for g in coupling_strengths:
                    for sI in strength_Is:
                        
                        
                        if not use_opt:
                            params['g_intercortical'] = np.round(ginter, 4)
                            params['coupling_strength'] = g 
                            params['strength_I'] = np.round(sI, 4)
                            params['Iext_duration'] = np.round(d, 4)
                            params['Iext_strength'] = s
                            params['Ib_strength'] = sb
                            params['Ib_noise_std'] = Ib_noise_std
                            params['area'] = area
                            params['resistance_factor'] = resistance_factor

                            # additional parameters (that are usually fixed)
                            params['g_thalPOm'] = np.round(g_thalPOm, 4)
                            params['delay_factor'] = delay_factor
                            params['receptor_thalamus_delay'] = receptor_thalamus_delay
                            params['thal_delay_factor'] = thal_delay_factor
                            params['e3b_tau'] = e3b_tau
                            params['e1_tau'] = e1_tau
                            params['e2_tau'] = e2_tau
                            params['p_2PVE'] = p_2PVE
                            params['p_4PVE'] = p_4PVE


                        # everything not set above (g_thal, sI_thal, the cell counts,
                        # thal_connect, ...) comes from simulation_parameter.json

                        if pyrates:
                            model = SomatoModelPyrates(params)
                        else:
                            model = SomatoModel(params)

                        # save the connectivity heatmap into this run's directory
                        if save_connectivity and not pyrates:
                            run_dir = model.prepare_run_dir(filedir)
                            model.plot_W_heatmap(save_dir=run_dir)
                            # per-area connectivity grids (columns = layers, rows = cell types)
                            pop_labels = model.get_population_labels()
                            pf.plot_connectivity(model.W, run_dir, pop_labels=pop_labels,
                                                    area=model.area, direction='in')
                            pf.plot_connectivity(model.W, run_dir, pop_labels=pop_labels,
                                                    area=model.area, direction='out')


                        # simulate rates and potentials
                        start = time.time()
                        model.simulate()
                        stop = time.time()
                        duration = stop - start
                        all_durations.append(duration)
                        #print("Simulation duration (in s):", duration)

                        # analyse signal (frequency spectra)
                        #model.analyse_signal(save_spectrum=True)
                        # time-frequency error against measured data
                        sim_dip = model.compute_dipoles(subID_elec)
                        tf_error, tf_sim, tf_target = model.compute_error_timefreq(tf_target_path, sim_dip)
                        print("TF error (log-MSE):", tf_error)
                        model.plot_timefreq_comparison(tf_sim, tf_target)

                        # time-course error against measured data
                        tc_error, tc_sim, tc_target = model.compute_error_timecourse(tc_target_path, sim_dip)
                        print("Time-course error (peak-norm MSE):", tc_error)
                        model.plot_timecourse_comparison(tc_sim, tc_target)

                        # pre-stimulus spectrum error against measured data
                        # ps_target_path is the raw measured CSV, so flatten both sides
                        ps_error, ps_sim, ps_target = model.compute_error_prestim_spectrum(
                            ps_target_path, sim_dip, flatten_sim=True, flatten_target=True)
                        print("Pre-stim spectrum error (rel-power MSE):", ps_error)
                        model.plot_prestim_spectrum_comparison(ps_target_path, sim_dip)  # saves to PRESTIM_SPECTRUM_DIR

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
                        #print('coupling strength', model.coupling_strength) 
                        #print('b input', model.Ib_strength) 
                        #print('Iext strength', model.Iext_strength) 
                        #print('Iext dur', model.Iext_duration) 

                        if save_results:
                            # create per-run output folder (holds params.json + all HDF5s for this run)
                            run_dir = model.prepare_run_dir(filedir)
                            start = time.time()
                            model.save_results_csv(run_dir, "results", save_full_potentials)
                            stop = time.time()
                            duration = stop - start
                            all_durations_saving.append(duration)
                            #print("Saving duration (in s):", duration)

                        start_plot = 1980
                        if plot_rates:
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
                                figure_dir,
                                sI,
                                g,
                                d,
                                sb,
                                s,
                                pop_labels=model.get_population_labels(),
                            )

            # if (len(coupling_strengths) > 1):
            #    # used to plot with coupling strength on the x-axis and max/min rate on the y
            #    pf.plot_minmax(all_rates, coupling_strengths)

print("Mean Simulation duration: ", np.mean(all_durations))
print("Mean Saving duration: ", np.mean(all_durations_saving))


# %%
