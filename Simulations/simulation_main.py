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
from mne.datasets import eegbci

location = "laptop"
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
    
figure_dir = os.path.join(SIMDIR, "Figures")


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

# load EEG data and forward model computation
# setup sample data for forward modelling
data_path_labels = sample.data_path()
subject = "fsaverage"
trans = "fsaverage"
src = data_path_labels / "subjects" / "fsaverage" / "bem" / "fsaverage-ico-5-src.fif"
bem = data_path_labels / "subjects" / "fsaverage" / "bem" / "fsaverage-5120-5120-5120-bem-sol.fif"
(raw_fname,) = eegbci.load_data(subjects=1, runs=[6])
raw = mne.io.read_raw_edf(raw_fname, preload=True)
# Read and set the EEG electrode locations, which are already in fsaverage's
# space (MNI space) for standard_1020:
eegbci.standardize(raw)
montage = mne.channels.make_standard_montage("standard_1005")
raw.set_montage(montage)

# %%
# TODO: read the forward solution instead of computing
# forward solution and leadfield computation
fwd = mne.make_forward_solution(
    raw.info, trans=trans, src=src, bem=bem, eeg=True, mindist=5.0, n_jobs=None
)
leadfield = fwd["sol"]["data"]
print(f"Leadfield size : {leadfield.shape[0]} sensors x {leadfield.shape[1]} dipoles")

# reduce forward solution to one orientation
fwd_fixed = mne.convert_forward_solution(
    fwd, surf_ori=True, force_fixed=True, use_cps=True
)
src_free = fwd["src"]
src_fixed = fwd_fixed["src"]

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
filedir = os.path.join(SIMDIR, 'simulation_results')
if not os.path.exists(filedir):
    os.makedirs(filedir)

# set parameters to loop over 
coupling_strengths = [15] #np.arange(0,55,5) #[100, 120, 140, 160]
backgrndI_strengths = [1] #np.arange(0,8,2) #[40, 60, 80] #,6,7]
input_durations = [0.3] #np.arange(0, 0.02, 0.002) # [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
input_strengths = [15] #np.arange(0,50,10)
strength_I = [0.26] #np.arange(0.2,0.44,0.02) #, 0.25, 0.26, 0.36]
ginters = [0.5] #np.arange(0,2,0.01)
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

                        # additional parameters (that are usually fixed)
                        params['g_thal'] = 2
                        params['sI_thal'] = 0.5
                        params['extI_cellcounts'] = 1000
                        params['bI_cellcounts'] = 100
                        params['thal_cellcounts'] = 500

                        if pyrates:
                            model = SomatoModelPyrates(params)
                        else:
                            model = SomatoModel(params)
                            #model.plot_W_heatmap()
                        
                        # simulate rates and potentials
                        start = time.time()
                        model.simulate()
                        stop = time.time()
                        duration = stop - start
                        all_durations.append(duration)
                        print("Simulation duration (in s):", duration)

                        # analyse signal (frequency spectra)
                        #   model.analyse_signal(save_spectrum=True)

                        # compute dipoles
                        sim_dip = model.compute_dipoles()
                        model.plot_dipoles(sim_dip, raw.info)
                        evoked, epochs = model.simulate_eeg(raw, data_path_labels, sim_dip, fwd, src_fixed)
                        model.plot_eeg(evoked, epochs)

                        # print important parameters
                        """
                        print('simulation_dur', model.simulation_dur)
                        print('step_size', model.step_size)
                        print('input_onset', model.input_onset) 
                        print('thal_connect', model.thal_connect) 
                        print('extI_cellcounts', model.extI_cellcounts) 
                        print('strength_I', model.strength_I) 
                        print('bI_cellcounts', model.bI_cellcounts) 
                        print('thal_cellcounts', model.thal_cellcounts) 
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
                            start = time.time()
                            model.save_results_csv(filedir, model.filename, save_full_potentials)
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
