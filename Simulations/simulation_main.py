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
import matplotlib.pyplot as plt
import jax.numpy as jnp
import pandas as pd
import time
import csv
from jr_model import JR_Model
import plotting_functions as pf

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

    filename = filename + ".h5"
    filename_rates = "rates" + filename

    # only safe every second datapoint
    resolution_tstep = 0.01
    print("tstep resolution", resolution_tstep)
    rates_downsampled = rates[:, :: int(1000 * resolution_tstep)]
    rates_df = pd.DataFrame(rates_downsampled.T)
    rates_df.to_hdf(
        os.path.join(filedir, filename_rates), index=False, key="data", mode="w"
    )

    # sum the potentials together and save them
    potential_sum = np.sum(potentials, axis=1)
    potential_sum_downsampled = potential_sum[:, :: int(1000 * resolution_tstep)]
    potential_df = pd.DataFrame(potential_sum_downsampled.T, columns=cells)
    filename = "potentials" + filename
    potential_df.to_hdf(
        os.path.join(filedir, filename), index=False, key="data", mode="w"
    )

    if full:
        # save all potentials additionally
        psp_filename = "full_" + filename
        write_3D_csv(os.path.join(filedir, psp_filename), potentials)


def write_3D_csv(filename, data):
    """
    Write results in form of a 3D hdf5 file.
    """

    with h5py.File(filename, "w") as f:
        f.create_dataset(dataset_name, data=data, compression="gzip")


# %%
def main():
    save_params = False
    save_results = False
    save_full_potentials = False  # if True the potential matrix is 3D, otherwise 2D
    plot_rates = True
    plot_potentials = False
    plot_all_potentials = False
    jax_mode = False

    # set coupling strengths, step size and cortex type (visual or somato)
    # connectivity reverse factor is the absolute cell count divided by
    connect_reverse_factor = (
        1  # 6448 # TODO: adapt this factor also to S2 cell populations!
    )
    # to simulate:
    # thalamus I to E inhibition

    coupling_strengths = [20, 40, 60]  # , 150, 200, 250, 300]
    balance_EI = [
        0.7, 0.8, 0.9
    ]  # np.arange(0, 1.1, 0.1) # excitation-inhibition balance (between 0 and 1)
    g_thal = 2
    bEI_thal = 0.5  # if g_thal is 0, this does not matter
    step_size = 1e-3
    area = "all"
    filedir = "/data/p_02989/Modelling/output/"

    # define input
    input_type = "step"  # other options are "step", "baseline" (equals input strength 0) or "background"
    input_onset = 1.001  # in sec
    simulation_dur = 2
    input_durations = [1.5]  # , 1, 1.5] # np.arange(0, 1, 1) # in sec
    input_strengths = [500]  # np.arange(0, 500, 100)
    backgrndI_strengths = [7]

    # connections within the thalamus
    # in this order: tEE, tEI, tIE, tII
    thal_connect = np.array([0, 0, 0, 0])
    extI_cellcounts = 1000
    bI_cellcounts = 100

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

                        filename = f"_g{g}_bEI{bEI}_Ib{sb}_Iextd{d}_{input_type}Iexts{s}_Ionset{input_onset}_tauVisual_thalJiang_thalUncon_S1S2Uncon"
                        # filename = f'_gE{int(gE*connect_reverse_factor)}gI{int(gI*connect_reverse_factor)}_{cortex_type}_IbStrength{sb}_Iduration{d}_{input_type}IextStrength{s}_Ionset{input_onset}_tauVisual_thalJiang_thalEI0_S1S2'

                        # create input array
                        Iext = create_Iext(
                            simulation_time, step_size, input_onset, d, s, input_type
                        )
                        Ib = create_Ibackground(simulation_time, step_size, sb)
                        gE = g * bEI / connect_reverse_factor
                        gI = g * (1 - bEI) / connect_reverse_factor
                        gE_thal = g_thal * bEI_thal
                        gI_thal = g_thal * (1 - bEI_thal)
                        # for now we use the same coupling strength for the thalamus connections as for the cortical connections
                        coupling_thalE = gE_thal
                        coupling_thalI = gI_thal
                        print("gE", gE)
                        print("gI", gI)
                        thal_connect_scaled = thal_connect / connect_reverse_factor

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
