import numpy as np
import os
import json
import argparse
from matplotlib.ticker import FormatStrFormatter, FuncFormatter, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from scipy.signal import find_peaks
import seaborn as sns
import pandas as pd
import ast
import h5py
from pandas.errors import InvalidIndexError

def compute_fft_peak_frequency(signal, step_size, fmin=0.0, fmax=None, window="hann"):
    """
    Compute dominant frequency of a 1D signal using FFT peak.
    Returns (frequency_hz, peak_power).
    """
    if signal is None:
        return np.nan, np.nan
    signal = np.asarray(signal)
    n = signal.size
    if n < 2:
        return np.nan, np.nan

    # detrend
    x = signal - np.mean(signal)

    # windowing
    if window == "hann":
        win = np.hanning(n)
        x = x * win

    # FFT
    fft_vals = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(n, d=step_size)
    power = np.abs(fft_vals) ** 2

    # exclude DC
    valid = freqs > 0.0
    if fmin is not None:
        valid &= freqs >= fmin
    if fmax is not None:
        valid &= freqs <= fmax

    if not np.any(valid):
        return np.nan, np.nan

    idx = np.argmax(power[valid])
    freq = freqs[valid][idx]
    peak_power = power[valid][idx]
    return freq, peak_power


def compute_window_frequency(df, rates_df, potentials_df, start_idx, stop_idx, prefix, step_size,
                             osc_gate_rate, osc_gate_potential, fmin=0.0, fmax=None):
    """
    Compute dominant frequency for each population within a time window.
    Gating is done via diffRate_{prefix} and diffPotential_{prefix}.
    Adds freqRate_{prefix}, freqPotential_{prefix}, fftPowerRate_{prefix}, fftPowerPotential_{prefix}.
    """
    rate_col = f"freqRate_{prefix}"
    pot_col = f"freqPotential_{prefix}"
    rate_pow_col = f"fftPowerRate_{prefix}"
    pot_pow_col = f"fftPowerPotential_{prefix}"

    # initialize columns with NaN
    for col in [rate_col, pot_col, rate_pow_col, pot_pow_col]:
        if col not in df.columns:
            df[col] = np.nan

    if start_idx is None or stop_idx is None or stop_idx <= start_idx:
        return

    window_rates = rates_df.iloc[start_idx:stop_idx]
    window_pots = potentials_df.iloc[start_idx:stop_idx]

    for pop in rates_df.columns:
        # gating by min-max diff
        gate_rate = df.loc[pop, f"diffRate_{prefix}"] if f"diffRate_{prefix}" in df.columns else np.nan
        gate_pot = df.loc[pop, f"diffPotential_{prefix}"] if f"diffPotential_{prefix}" in df.columns else np.nan

        if np.isfinite(gate_rate) and gate_rate > osc_gate_rate:
            freq, pwr = compute_fft_peak_frequency(window_rates[pop].values, step_size, fmin=fmin, fmax=fmax)
            df.loc[pop, rate_col] = freq
            df.loc[pop, rate_pow_col] = pwr

        if np.isfinite(gate_pot) and gate_pot > osc_gate_potential:
            freq, pwr = compute_fft_peak_frequency(window_pots[pop].values, step_size, fmin=fmin, fmax=fmax)
            df.loc[pop, pot_col] = freq
            df.loc[pop, pot_pow_col] = pwr

def read_simulation_data(output_dir, figure_dir, input_durations, input_strengths, coupling_strengths, strength_I,  
                        backgroundI_strengths, step_size, sample_delay_immediate, sample_delay_late, input_onset, sample_dur, cortex_type, input_type,
                        thal_cellcounts, bI_cellcounts, extI_cellcounts, load_trajectory, load_full_potentials, load_population_potential = 3, offset=0.1):
    '''
    Parameter:
    sample_dur: duration of sampling baseline and longterm activity in s 
    load_population_potential : the index of the cell population of which the potential (all incoming input to the population) will be loaded - (if it was saved after the simulation)
    offset : time in s between baseline sampling and start of input 
    '''
    # create summary matrix: should include (max/min_value x population   
    min_data = []
    max_data = []
    summary_df = pd.DataFrame()
    time_data_df = pd.DataFrame()
    potential_data_df = pd.DataFrame()

    for d in input_durations:
        for s in input_strengths:
            for g in coupling_strengths:
                for sI in strength_I:
                    for bI in backgroundI_strengths:

                        print('\nInput duration:', d)
                        print('Input strength:', s)
                        print('Background Input strength:', bI)
                        print('g', g)
                        print('sI', sI)

                        df = pd.DataFrame()

                        # read in firing rates in data matrix (datapoints x populations)
                        filename = f"g{g}_sI{sI}_Ib{bI}_Iextd{d}_{input_type}Iexts{s}_Ionset{input_onset}_thalcells{thal_cellcounts}_Ibcells{bI_cellcounts}_Iextcells{extI_cellcounts}_thalUncon_S1S2Uncon.hdf5"
                        rates_df = pd.read_hdf(os.path.join(output_dir, filename), key='rates')
                        
                        # read in potentials in data matrix (datapoints x populations)
                        potentials_df = pd.read_hdf(os.path.join(output_dir, filename), key='summed_potential')

                        # LONG TERM immediate (sample for x ms starting x secs after offset)
                        start_sample_immediate = int((input_onset+d+sample_delay_immediate)/step_size)
                        stop_sample_immediate = int(start_sample_immediate + sample_dur/step_size)
                        df['rate_immediateLongterm'] = rates_df.iloc[start_sample_immediate:stop_sample_immediate].mean()
                        df['minRate_immediateLongterm'] = rates_df.iloc[start_sample_immediate:stop_sample_immediate].min()
                        df['maxRate_immediateLongterm'] = rates_df.iloc[start_sample_immediate:stop_sample_immediate].max()
                        df['diffRate_immediateLongterm'] = df['maxRate_immediateLongterm'] - df['minRate_immediateLongterm']
                        df['potential_immediateLongterm'] = potentials_df.iloc[start_sample_immediate:stop_sample_immediate].mean()
                        df['minPotential_immediateLongterm'] = potentials_df.iloc[start_sample_immediate:stop_sample_immediate].min()
                        df['maxPotential_immediateLongterm'] = potentials_df.iloc[start_sample_immediate:stop_sample_immediate].max()
                        df['diffPotential_immediateLongterm'] = df['maxPotential_immediateLongterm'] - df['minPotential_immediateLongterm']

                        # LONG TERM late (sample for x ms starting at x sec after offset)
                        start_sample_late = int((input_onset+d+sample_delay_late)/step_size)
                        stop_sample_late = int(start_sample_late + sample_dur/step_size)
                        df['rate_lateLongterm'] = rates_df.iloc[start_sample_late:stop_sample_late].mean()
                        df['minRate_lateLongterm'] = rates_df.iloc[start_sample_late:stop_sample_late].min()
                        df['maxRate_lateLongterm'] = rates_df.iloc[start_sample_late:stop_sample_late].max()
                        df['diffRate_lateLongterm'] = df['maxRate_lateLongterm'] - df['minRate_lateLongterm']
                        df['potential_lateLongterm'] = potentials_df.iloc[start_sample_late:stop_sample_late].mean()
                        df['minPotential_lateLongterm'] = potentials_df.iloc[start_sample_late:stop_sample_late].min()
                        df['maxPotential_lateLongterm'] = potentials_df.iloc[start_sample_late:stop_sample_late].max()
                        df['diffPotential_lateLongterm'] = df['maxPotential_lateLongterm'] - df['minPotential_lateLongterm']
                        
                        # DURING INPUT
                        start_sample_during = int((input_onset)/step_size)
                        stop_sample_during = int((input_onset+d)/step_size)
                        df['rate_duringInput'] = rates_df.iloc[start_sample_during:stop_sample_during].mean()
                        df['minRate_duringInput'] = rates_df.iloc[start_sample_during:stop_sample_during].min()
                        df['maxRate_duringInput'] = rates_df.iloc[start_sample_during:stop_sample_during].max()
                        df['diffRate_duringInput'] = df['maxRate_duringInput'] - df['minRate_duringInput']
                        df['potential_duringInput'] = potentials_df.iloc[start_sample_during:stop_sample_during].mean()
                        df['minPotential_duringInput'] = potentials_df.iloc[start_sample_during:stop_sample_during].min()
                        df['maxPotential_duringInput'] = potentials_df.iloc[start_sample_during:stop_sample_during].max()
                        df['diffPotential_duringInput'] = df['maxPotential_duringInput'] - df['minPotential_duringInput']
                    
                        # BASELINE SAMPLE
                        baseline_start = int((input_onset - (sample_dur+offset))/step_size) 
                        baseline_stop = int(baseline_start + sample_dur/step_size)
                        df['rate_baseline'] = rates_df.iloc[baseline_start:baseline_stop].mean()
                        df['minRate_baseline'] = rates_df.iloc[baseline_start:baseline_stop].min()
                        df['maxRate_baseline'] = rates_df.iloc[baseline_start:baseline_stop].max()
                        df['diffRate_baseline'] = df['maxRate_baseline'] - df['minRate_baseline']
                        df['potential_baseline'] = potentials_df.iloc[baseline_start:baseline_stop].mean()
                        df['minPotential_baseline'] = potentials_df.iloc[baseline_start:baseline_stop].min()
                        df['maxPotential_baseline'] = potentials_df.iloc[baseline_start:baseline_stop].max()
                        df['diffPotential_baseline'] = df['maxPotential_baseline'] - df['minPotential_baseline']
                        
                        # long term to baseline comparison (here we don't take the diffRate because the baseline does not 
                        # oscillate and we want to compare if the long term activity is still larger than baseline)
                        df['longtermVSbaseline_rate'] = df['rate_lateLongterm'] - df['rate_baseline']
                        df['duringInputVSbaseline_rate'] = df['rate_duringInput'] - df['rate_baseline']
                        df['longtermVSbaseline_potential'] = df['potential_lateLongterm'] - df['potential_baseline']
                        df['duringInputVSbaseline_potential'] = df['potential_duringInput'] - df['potential_baseline']
                        
                        # name the behaviour switch through, memory and non reponsive
                        #! TODO: adjust the threshold values for dynamic functions 
                        non_responsive = np.abs(df['duringInputVSbaseline_rate'])<0.1
                        transfer = ((np.abs(df['duringInputVSbaseline_rate'])>0.1) & (np.abs(df['longtermVSbaseline_rate'])<0.1))
                        memory = (np.abs(df['longtermVSbaseline_rate'])>0.1)
                        df.loc[non_responsive,'dynamic_function_rate'] = 1 # 'non-responsive'
                        df.loc[transfer,'dynamic_function_rate'] = 2 # 'transfer'
                        df.loc[memory,'dynamic_function_rate'] = 3 # 'memory'
                        
                        non_responsive = np.abs(df['duringInputVSbaseline_potential'])<0.001
                        transfer = ((np.abs(df['duringInputVSbaseline_potential'])>0.001) & (np.abs(df['longtermVSbaseline_potential'])<0.001))
                        memory = (np.abs(df['longtermVSbaseline_potential'])>0.001)

                        """if d == 1.5:
                            if (not non_responsive.iloc[2]) & (not transfer.iloc[2]) & (not memory.iloc[2]): 
                                print('no function')
                                print(np.abs(df['duringInputVSbaseline_potential'][2]))
                                print(np.abs(df['longtermVSbaseline_potential'][2]))"""

                        df.loc[non_responsive,'dynamic_function_potential'] = 1 #'non-responsive'
                        df.loc[transfer,'dynamic_function_potential'] = 2 #'transfer'
                        df.loc[memory,'dynamic_function_potential'] = 3 #'memory'

                        # TODO: also include where there are oscillations or not (maybe in the non-responsive behaviour)

                        # set info in data frame
                        df['population'] = potentials_df.columns
                        df['globalCoupling'] = g
                        df['strength_I'] = sI
                        df['InputDuration'] = d
                        df['InputStrength'] = s
                        df['BckgndInputStrength'] = bI
                        summary_df = pd.concat([summary_df, df])

                        if load_trajectory:
                            # load entire trajectory in data frame
                            for pop_rate, pop_potential in zip(rates_df.items(), potentials_df.items()):
                                pop_name = pop_potential[0]
                                rate_trajectory = pop_rate[1]
                                potential_trajectory = pop_potential[1]
                                pop_df = pd.DataFrame()
                                pop_df['population'] = [pop_name]*len(rate_trajectory)
                                pop_df['rate'] = rate_trajectory
                                pop_df['potential'] = potential_trajectory
                                pop_df['time'] = rate_trajectory.index.values * step_size # time in s
                                pop_df['global_coupling'] = g
                                pop_df['strength_I'] = sI
                                pop_df['InputDuration'] = d
                                pop_df['InputStrength'] = s 
                                pop_df['BckgndInputStrength'] = bI
                                time_data_df = pd.concat([time_data_df, pop_df])

                        if load_full_potentials:
                            # load full 3D (target x source x timestep) potential of selected population 
                            filename = f"full_potentials_g{g}_sI{sI}_{cortex_type}_IbStrength{bI}_Iduration{d}_{input_type}IextStrength{s}_Ionset{input_onset}_tauVisual_{thalcells}_thalEI0.csv"
                            potential_df = pd.read_csv(os.path.join(output_dir, filename), sep=',', header=None)
                            potential_df = potential_df.applymap(lambda x: ast.literal_eval(x) if isinstance(x, str) and x.startswith('[') and x.endswith(']') else x)
                            potential_df = potential_df.iloc[load_population_potential]
                            longterm_pot_df = pd.DataFrame()

                            longterm_input = []
                            for sourcepop in potential_df:
                                # extract the longterm behaviour
                                longterm_input.append(sourcepop[start_sample_late:stop_sample_late])
                            
                            longterm_input = np.sum(np.array(longterm_input), axis=0)
                            longterm_mean = np.sum(longterm_input)
                            print(longterm_input.shape)

                            longterm_pot_df['lt_potential'] = longterm_input
                            longterm_pot_df['population'] = load_population_potential
                            longterm_pot_df['global_coupling'] = g
                            longterm_pot_df['strength_I'] = sI
                            longterm_pot_df['InputDuration'] = d
                            longterm_pot_df['InputStrength'] = s 
                            potential_data_df = pd.concat([potential_data_df, longterm_pot_df]) 

    
    return summary_df, time_data_df, potential_data_df


def compute_longeterm_immediate(df, rates_df, potentials_df, input_onset, d, step_size, sample_delay_immediate, sample_dur):
    # LONG TERM immediate (sample for x ms starting x secs after offset)
    start_sample_immediate = int((input_onset+d+sample_delay_immediate)/step_size)
    stop_sample_immediate = int(start_sample_immediate + sample_dur/step_size)
    df['rate_immediateLongterm'] = rates_df.iloc[start_sample_immediate:stop_sample_immediate].mean()
    df['minRate_immediateLongterm'] = rates_df.iloc[start_sample_immediate:stop_sample_immediate].min()
    df['maxRate_immediateLongterm'] = rates_df.iloc[start_sample_immediate:stop_sample_immediate].max()
    df['diffRate_immediateLongterm'] = df['maxRate_immediateLongterm'] - df['minRate_immediateLongterm']
    df['potential_immediateLongterm'] = potentials_df.iloc[start_sample_immediate:stop_sample_immediate].mean()
    df['minPotential_immediateLongterm'] = potentials_df.iloc[start_sample_immediate:stop_sample_immediate].min()
    df['maxPotential_immediateLongterm'] = potentials_df.iloc[start_sample_immediate:stop_sample_immediate].max()
    df['diffPotential_immediateLongterm'] = df['maxPotential_immediateLongterm'] - df['minPotential_immediateLongterm']

def compute_longeterm_late(df, rates_df, potentials_df, input_onset, d, step_size, sample_delay_late, sample_dur):
    # LONG TERM late (sample for x ms starting at x sec after offset)
    start_sample_late = int((input_onset+d+sample_delay_late)/step_size)
    stop_sample_late = int(start_sample_late + sample_dur/step_size)
    df['rate_lateLongterm'] = rates_df.iloc[start_sample_late:stop_sample_late].mean()
    df['minRate_lateLongterm'] = rates_df.iloc[start_sample_late:stop_sample_late].min()
    df['maxRate_lateLongterm'] = rates_df.iloc[start_sample_late:stop_sample_late].max()
    df['diffRate_lateLongterm'] = df['maxRate_lateLongterm'] - df['minRate_lateLongterm']
    df['potential_lateLongterm'] = potentials_df.iloc[start_sample_late:stop_sample_late].mean()
    df['minPotential_lateLongterm'] = potentials_df.iloc[start_sample_late:stop_sample_late].min()
    df['maxPotential_lateLongterm'] = potentials_df.iloc[start_sample_late:stop_sample_late].max()
    df['diffPotential_lateLongterm'] = df['maxPotential_lateLongterm'] - df['minPotential_lateLongterm']

def compute_input_response(df, rates_df, potentials_df, input_onset, d, step_size, response_window):
    # DURING INPUT
    start_sample_during = int((input_onset)/step_size)
    stop_sample_during = int((input_onset+d+response_window)/step_size)
    df['rate_duringInput'] = rates_df.iloc[start_sample_during:stop_sample_during].mean()
    df['minRate_duringInput'] = rates_df.iloc[start_sample_during:stop_sample_during].min()
    df['maxRate_duringInput'] = rates_df.iloc[start_sample_during:stop_sample_during].max()
    df['diffRate_duringInput'] = df['maxRate_duringInput'] - df['minRate_duringInput']
    df['potential_duringInput'] = potentials_df.iloc[start_sample_during:stop_sample_during].mean()
    df['minPotential_duringInput'] = potentials_df.iloc[start_sample_during:stop_sample_during].min()
    df['maxPotential_duringInput'] = potentials_df.iloc[start_sample_during:stop_sample_during].max()
    df['diffPotential_duringInput'] = df['maxPotential_duringInput'] - df['minPotential_duringInput']

def compute_baseline(df, rates_df, potentials_df, input_onset, step_size, sample_dur, offset):       
    # BASELINE SAMPLE (bofore input onset)
    baseline_start = int((input_onset - (sample_dur+offset))/step_size) 
    baseline_stop = int(baseline_start + sample_dur/step_size)
    df['rate_baseline'] = rates_df.iloc[baseline_start:baseline_stop].mean()
    df['minRate_baseline'] = rates_df.iloc[baseline_start:baseline_stop].min()
    df['maxRate_baseline'] = rates_df.iloc[baseline_start:baseline_stop].max()
    df['diffRate_baseline'] = df['maxRate_baseline'] - df['minRate_baseline']
    
    df['potential_baseline'] = potentials_df.iloc[baseline_start:baseline_stop].mean()
    df['minPotential_baseline'] = potentials_df.iloc[baseline_start:baseline_stop].min()
    df['maxPotential_baseline'] = potentials_df.iloc[baseline_start:baseline_stop].max()
    df['diffPotential_baseline'] = df['maxPotential_baseline'] - df['minPotential_baseline']

def compare_longterm_baseline(df):                    
    # long term to baseline comparison (here we don't take the diffRate because the baseline does not 
    # oscillate and we want to compare if the long term activity is still larger than baseline)
    df['longtermVSbaseline_rate'] = df['rate_lateLongterm'] - df['rate_baseline']
    df['duringInputVSbaseline_rate'] = df['rate_duringInput'] - df['rate_baseline']
    df['longtermVSbaseline_potential'] = df['potential_lateLongterm'] - df['potential_baseline']
    df['duringInputVSbaseline_potential'] = df['potential_duringInput'] - df['potential_baseline']

def classify_response(df):
    """name the behaviour after stimulation: transfer, memory or non reponsive"""
    
    # rates
    non_responsive = np.abs(df['duringInputVSbaseline_rate'])<0.1
    transfer = ((np.abs(df['duringInputVSbaseline_rate'])>0.1) & (np.abs(df['longtermVSbaseline_rate'])<0.1))
    memory = (np.abs(df['longtermVSbaseline_rate'])>0.1)
    df.loc[non_responsive,'dynamic_function_rate'] = 1 # 'non-responsive'
    df.loc[transfer,'dynamic_function_rate'] = 2 # 'transfer'
    df.loc[memory,'dynamic_function_rate'] = 3 # 'memory'
    
    # potentials
    non_responsive = np.abs(df['duringInputVSbaseline_potential'])<0.1
    transfer = ((np.abs(df['duringInputVSbaseline_potential'])>0.1) & (np.abs(df['longtermVSbaseline_potential'])<0.1))
    memory = (np.abs(df['longtermVSbaseline_potential'])>0.1)
    df.loc[non_responsive,'dynamic_function_potential'] = 1 #'non-responsive'
    df.loc[transfer,'dynamic_function_potential'] = 2 #'transfer'
    df.loc[memory,'dynamic_function_potential'] = 3 #'memory'
    

def set_sim_info(df, potentials_df, g, sI, d, s, bI):
    """ set info in data frame """
    df['population'] = potentials_df.columns
    df['globalCoupling'] = g
    df['strength_I'] = sI
    df['InputDuration'] = d
    df['InputStrength'] = s
    df['BckgndInputStrength'] = bI


def load_trajectory(rates_df, potentials_df, g, sI, d, s, bI, step_size):
    """ load entire trajectory in data frame """
    time_data_df = pd.DataFrame()
    for pop_rate, pop_potential in zip(rates_df.items(), potentials_df.items()):
        pop_name = pop_potential[0]
        rate_trajectory = pop_rate[1]
        potential_trajectory = pop_potential[1]
        pop_df = pd.DataFrame()
        pop_df['population'] = [pop_name]*len(rate_trajectory)
        pop_df['rate'] = rate_trajectory
        pop_df['potential'] = potential_trajectory
        pop_df['time'] = rate_trajectory.index.values * step_size # time in s
        pop_df['global_coupling'] = g
        pop_df['strength_I'] = sI
        pop_df['InputDuration'] = d
        pop_df['InputStrength'] = s 
        pop_df['BckgndInputStrength'] = bI
        time_data_df = pd.concat([time_data_df, pop_df])

    return time_data_df

def parse_coupling():
    """ we parallelize over different coupling strengths (in srun HPC script), so we need to read in those """   
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--g",
        type=float,
        nargs="+",
        help="coupling strengths",
        required=False,
    )
    coupling_strengths = parser.parse_args().g
    return coupling_strengths


def load_parameters(WDDIR): 
    """ Read in preprocessing parameters """
    with open(os.path.join(WDDIR, 'Analysis', 'analysis_parameter.json'), 'r') as json_file:
        params = json.load(json_file)
    return params

def load_simulation_data(g, sI, bI, d, s, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, data_dir, pyrates=False):
    """ Read simulation data """
    # read in firing rates in data matrix (datapoints x populations)
    if pyrates:
        filename = f"gthal2_sIthal0.5_g{g}_sI{sI}_Ib{bI}_Iextd{d}_{input_type}Iexts{s}_Ionset{input_onset}_thalcells{thal_cellcounts}_Ibcells{bI_cellcounts}_Iextcells{extI_cellcounts}_PYRATES.hdf5"
    else:    
        filename = f"gthal2_sIthal0.5_g{g}_sI{sI}_Ib{bI}_Iextd{d}_{input_type}Iexts{s}_Ionset{input_onset}_thalcells{thal_cellcounts}_Ibcells{bI_cellcounts}_Iextcells{extI_cellcounts}_thalUncon_S1S2Uncon.hdf5"

    rates_df = pd.read_hdf(os.path.join(data_dir, filename), key='rates')
    
    # read in potentials in data matrix (datapoints x populations)
    potentials_df = load_hdf_safe(os.path.join(data_dir, filename))

    return rates_df, potentials_df, filename

def load_hdf_safe(fname):
    key = "summed_potential"
    with h5py.File(fname, 'r') as f:
        grp    = f[key]
        values = grp['block0_values'][:]
        cols   = [c.decode() for c in grp['block0_items'][:]]  # use block0_items, NOT axis0
        idx    = grp['axis1'][:]
    potentials_safe = pd.DataFrame(values, columns=cols, index=idx)
    return potentials_safe

def load_derivative(g, sI, bI, d, s, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, deriv_dir):
    """ load the characteristics/processed data of one simulation """
    deriv_file = f"gthal2_sIthal0.5_g{g}_sI{sI}_Ib{bI}_Iextd{d}_{input_type}Iexts{s}_Ionset{input_onset}_thalcells{thal_cellcounts}_Ibcells{bI_cellcounts}_Iextcells{extI_cellcounts}_thalUncon_S1S2Uncon_processed.csv"
    deriv_df = pd.read_csv(os.path.join(deriv_dir, deriv_file))
    return deriv_df

def check_list(sim_param):
    if not isinstance(sim_param, list):
        sim_param = [sim_param]

    return sim_param

def load_all_derivatives(Iext_dur, Iext_str, gs, sIs, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir):
    """ load all selected combinations of processed simulation (with their characteristics) into one dataframe
    All parameters can be either lists of single values. If they are single values, 
    convert them into lists with one entry.

    Parameters:
    -----------
    Iext_dur : float or list of input durations
    Iext_str : float or list of input strengths
    g : float or list of coupling strengths
    sI : float or list of E/I balance values
    Ib_str : float or list of background input strengths
    """

    # check if parameters are lists, if not convert to list
    Iext_dur = check_list(Iext_dur)
    Iext_str = check_list(Iext_str)
    gs = check_list(gs)
    sIs = check_list(sIs)
    Ib_str = check_list(Ib_str)

    data_df = pd.DataFrame()
    for d in Iext_dur:
            for s in Iext_str:
                for g in gs:
                    for sI in sIs:
                        for bI in Ib_str:
                            data_single_df = load_derivative(g, sI, bI, d, s, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir)
                            data_df = pd.concat([data_df, data_single_df], ignore_index=True)

    return data_df
