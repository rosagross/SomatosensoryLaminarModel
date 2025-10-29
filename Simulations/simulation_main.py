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
from jr_model import JR_Model
import plotting_functions as pf

# %%
SIMDIR = os.getenv("SIMDIR")
WDDIR = os.getenv("WDDIR")
figure_dir = os.path.join(SIMDIR, "Figures")

# %%    

def main():
    # read simulation params
    params = read_simulation_params()
    
    # we parallelize over different coupling strengths (in srun HPC script)   
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--g",
        type=float,
        nargs="+",
        help="coupling strengths",
        required=False,
    )
    #coupling_strengths = parser.parse_args().g
    coupling_strengths = params['coupling_strengths'] # coupling_strengths

    # Assign variables from loaded parameters
    save_params = params['save_params']
    save_results = params['save_results']
    save_full_potentials = params['save_full_potentials']
    plot_rates = params['plot_rates']
    plot_potentials = params['plot_potentials']
    plot_all_potentials = params['plot_all_potentials']
    jax_mode = params['jax_mode']

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
    backgrndI_strengths = params['backgrndI_strengths']

    # connectivity 
    thal_connect = np.array(params['thal_connect'])
    extI_cellcounts = params['extI_cellcounts']
    bI_cellcounts = params['bI_cellcounts']
    thal_cellcounts = params['thal_cellcounts']

    # implement option to choose only one part of the cortical circuit
    # this can be done by putting the connectivity for those parts to zero
    # options:
    # 1. only thalamus & Area 3b
    # 2. A3b and Area 1
    # 3. A3b and Area 1 and thalamus
    # 4. only Area 1
    # 5. only Area 1 and S2 and thalamus
    # 6. only Area 1 and S2
    # 7. only S2
    # 8. default: all

    filedir = os.path.join(SIMDIR, 'simulation_results')
    if not os.path.exists(filedir):
        os.makedirs(filedir)


    for d in input_durations:
        simulation_time = int(input_onset) + simulation_dur
        for sb in backgrndI_strengths:
            for s in input_strengths:
                # arrays to store rate (for plotting)
                all_rates = []
                all_potentials = []
                all_durations = []
                all_durations_saving = []

                for g in coupling_strengths:
                    for bEI in balance_EI:
                        print("\nInput duration:", d)
                        print("Input strength:", s)
                        print("Background Input strength:", sb)
                        print("g", g)
                        print("bEI", bEI)
                        print(
                            f"Thalamus EtoE:{thal_connect[0]} ItoE: {thal_connect[1]}"
                        )
                        print(
                            f"Thalamus EtoI:{thal_connect[2]} ItoI: {thal_connect[3]}"
                        )

                        filename = f"g{g}_bEI{bEI}_Ib{sb}_Iextd{d}_{input_type}Iexts{s}_Ionset{input_onset}_thalcells{thal_cellcounts}_Ibcells{bI_cellcounts}_Iextcells{extI_cellcounts}_thalUncon_S1S2Uncon"

                        model = JR_Model(
                            sb,
                            s,
                            d,
                            
                        )

                        if save_params:
                            # safe connectivty parameter in yaml file
                            model.p.save_to_yaml(
                                os.path.join(filedir, "params" + filename),
                                sb,
                                s,
                                d,
                                g,
                                g_thal,
                                bEI,
                                bEI_thal,
                                thal_connect,
                            )

                        # perform simulation with current coupling strength g
                        start = time.time()
                        rate, potential = model.run_simulation()

                        if jax_mode:
                            rate.block_until_ready()

                        stop = time.time()
                        duration = stop - start
                        all_durations.append(duration)
                        print("Simulation duration (in s):", duration)

                        # append results
                        all_rates.append(rate)
                        all_potentials.append(potential)

                        if save_results:
                            start = time.time()
                            save_results_csv(
                                rate, potential, filedir, filename, save_full_potentials
                            )
                            stop = time.time()
                            duration = stop - start
                            all_durations_saving.append(duration)
                            print("Saving duration (in s):", duration)

                        if plot_rates:
                            start_plot = 0
                            pf.plot_results(
                                rate,
                                Iext,
                                Ib,
                                step_size,
                                simulation_time,
                                start_plot,
                                bEI,
                                g,
                                area,
                                d,
                                sb,
                                s,
                                figure_dir
                            )

                        if plot_potentials:
                            potential_sum = np.sum(potential, axis=1)
                            resolution_tstep = 1e-2
                            potential_sum_downsampled = potential_sum[
                                :, :: int(1000 * resolution_tstep)
                            ]
                            pf.plot_potentials(
                                potential_sum,
                                Iext,
                                Ib,
                                step_size,
                                simulation_time,
                                start_plot,
                                figure_dir,
                                bEI,
                                g,
                                d,
                                sb,
                                s
                            )

                        if plot_all_potentials:
                            pf.plot_all_potentials(
                                potential,
                                Iext,
                                Ib,
                                step_size,
                                simulation_time,
                                start_plot,
                            )

                # if (len(coupling_strengths) > 1):
                #    # used to plot with coupling strength on the x-axis and max/min rate on the y
                #    pf.plot_minmax(all_rates, coupling_strengths)

    print("Mean Simulation duration: ", np.mean(all_durations))
    print("Mean Saving duration: ", np.mean(all_durations_saving))

    return potential, rate

if __name__ == "__main__":
    potential, rate = main()

