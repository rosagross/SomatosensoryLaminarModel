"""
File: simulation_main.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2025-08-05
Description: Run this file to run the simulation! 

"""

# %%
import numpy as np
import os
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

def create_Iext(
    simulation_time, step_size, input_onset, input_duration, input_strength, input_type
):
    """Creates external input."""

    Iext = np.zeros(int(simulation_time / step_size))

    if input_type == "step":
        t = int(input_duration / step_size)
        t0 = int(input_onset / step_size)
        Iext[t0 : t0 + t] = input_strength
    elif input_type == "background":
        # provide input for the entire simulation duration
        Iext[:] = input_strength

    return Iext


def create_Ibackground(simulation_time, step_size, input_strength):
    """Create Background Input"""
    Ib = np.zeros(int(simulation_time / step_size))
    Ib[:] = input_strength
    return Ib


def save_results_csv(rates, potentials, filedir, filename, full=False):
    """
    Safe the simulated data in a csv file
    """

    population_names = [
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
    ]
    cells = np.concatenate((population_names, ["ThalE", "ThalI"]))

    filename = filename + ".hdf5"

    # only safe every second datapoint
    resolution_tstep = 0.01
    print("tstep resolution", resolution_tstep)
    rates_downsampled = rates[:, :: int(1000 * resolution_tstep)]
    rates_df = pd.DataFrame(rates_downsampled.T, columns=cells)
    rates_df.to_hdf(
        os.path.join(filedir, filename), index=False, key="rates", mode="a"
    )

    # sum the potentials together and save them
    potential_sum = np.sum(potentials, axis=1)
    potential_sum_downsampled = potential_sum[:, :: int(1000 * resolution_tstep)]
    potential_df = pd.DataFrame(potential_sum_downsampled.T, columns=cells)
    potential_df.to_hdf(
        os.path.join(filedir, filename), index=False, key="summed_potential", mode="a"
    )

    if full:
        # save all potentials additionally
        psp_filename = "full_" + filename
        write_3D_csv(os.path.join(filedir, psp_filename), potentials)

# TODO: implement saving in hdf5 format
def write_3D_csv(filename, data):
    """
    Write results in form of a 3D hdf5 file.
    """

    with h5py.File(filename, "w") as f:
        f.create_dataset(dataset_name, data=data, compression="gzip")


# %%

def read_simulation_params():
    """Read simulation parameters from json file."""
    # Read in preprocessing parameters
    with open(os.path.join(WDDIR, 'Simulations', 'simulation_parameter.json'), 'r') as json_file:
        params = json.load(json_file)
    
    return params

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
    coupling_strengths = parser.parse_args().g
    #coupling_strengths = params['coupling_strengths'] # coupling_strengths

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

                        # create input array
                        Iext = create_Iext(
                            simulation_time, step_size, input_onset, d, s, input_type
                        )
                        Ib = create_Ibackground(simulation_time, step_size, sb)
                        gE = g * bEI 
                        gI = g * (1 - bEI)
                        gE_thal = g_thal * bEI_thal
                        gI_thal = g_thal * (1 - bEI_thal)
                        # for now we use the same coupling strength for the thalamus connections as for the cortical connections
                        coupling_thalE = gE_thal
                        coupling_thalI = gI_thal
                        print("gE", gE)
                        print("gI", gI)
                        thal_connect_scaled = thal_connect 

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
                            simulation_time,
                            area=area,
                        )
                        if save_params:
                            # safe connectivty parameter in yaml file
                            model.p.save_to_yaml(
                                os.path.join(filedir, "params" + filename),
                                gE,
                                gI,
                                coupling_thalE,
                                coupling_thalI,
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

# %%
