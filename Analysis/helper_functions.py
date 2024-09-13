import numpy as np
import os
from matplotlib.ticker import FormatStrFormatter, FuncFormatter, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from scipy.signal import find_peaks
import seaborn as sns
import pandas as pd
import ast

def read_simulation_data(output_dir, figure_dir, input_durations, input_strengths, coupling_strengths, 
                        step_size, sample_delay, input_onset, sample_dur, cortex_type, input_type,
                        thalamus_source, load_trajectory, load_full_potentials, load_population_potential = 3):
    
    # create summary matrix: should include (max/min_value x population   
    min_data = []
    max_data = []
    summary_df = pd.DataFrame()
    time_data_df = pd.DataFrame()
    potential_data_df = pd.DataFrame()

    for d in input_durations:
        for s in input_strengths:
            for g in coupling_strengths:
                df = pd.DataFrame()

                # read in firing rates in data matrix (datapoints x populations)
                filename_rates = f"rates_G{g}_{cortex_type}_Iduration{d}_{input_type}Istrength{s}_Ionset{input_onset}_tauVisual_{thalamus_source}.csv"
                rates_df = pd.read_csv(os.path.join(output_dir, filename_rates))
                
                # read in potentials in data matrix (datapoints x populations)
                filename_potentials = f"potentials_G{g}_{cortex_type}_Iduration{d}_{input_type}Istrength{s}_Ionset{input_onset}_tauVisual_{thalamus_source}.csv"
                potentials_df = pd.read_csv(os.path.join(output_dir, filename_potentials))

                # LONG TERM (sample for 100ms starting at 0.5 sec after offset)
                start_sample = int((input_onset+d+sample_delay)/step_size)
                stop_sample = int(start_sample + sample_dur/step_size)
                df['rate_longterm'] = rates_df.iloc[start_sample:stop_sample].mean()
                df['minRate_longterm'] = rates_df.iloc[start_sample:stop_sample].min()
                df['maxRate_longterm'] = rates_df.iloc[start_sample:stop_sample].max()
                df['diffRate_longterm'] = df['maxRate_longterm'] - df['minRate_longterm']
                df['potential_longterm'] = potentials_df.iloc[start_sample:stop_sample].mean()
                df['minPotential_longterm'] = potentials_df.iloc[start_sample:stop_sample].min()
                df['maxPotential_longterm'] = potentials_df.iloc[start_sample:stop_sample].max()
                df['diffPotential_longterm'] = df['maxPotential_longterm'] - df['minPotential_longterm']
                
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
                offset = 0.1 # time between baseline sampling and start of input 
                baseline_start = int((input_onset - (sample_dur+offset))/step_size) 
                baseline_stop = int(baseline_start + sample_dur/step_size)
                df['rate_baseline'] = rates_df.iloc[baseline_start:baseline_stop].mean()
                df['potential_baseline'] = potentials_df.iloc[baseline_start:baseline_stop].mean()
                
                # long term to baseline comparison (here we don't take the diffRate because the baseline does not 
                # oscillate and we want to compare if the long term activity is still larger than baseline)
                df['longtermVSbaseline_rate'] = df['rate_longterm'] - df['rate_baseline']
                df['duringInputVSbaseline_rate'] = df['rate_duringInput'] - df['rate_baseline']
                df['longtermVSbaseline_potential'] = df['potential_longterm'] - df['potential_baseline']
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

                if d == 1.5:
                    if (not non_responsive.iloc[2]) & (not transfer.iloc[2]) & (not memory.iloc[2]): 
                        print('no function')
                        print(np.abs(df['duringInputVSbaseline_potential'][2]))
                        print(np.abs(df['longtermVSbaseline_potential'][2]))

                df.loc[non_responsive,'dynamic_function_potential'] = 1 #'non-responsive'
                df.loc[transfer,'dynamic_function_potential'] = 2 #'transfer'
                df.loc[memory,'dynamic_function_potential'] = 3 #'memory'

                # TODO: also include where there are oscillations or not (maybe in the non-responsive behaviour)

                # set info in data frame
                df['population'] = df.index.values
                df['coupling_strength'] = g
                df['InputDuration'] = d
                df['InputStrength'] = s
                summary_df = pd.concat([summary_df, df])

                if load_trajectory:
                    # load entire trajectory in data frame
                    for pop_rate, pop_potential in zip(rates_df.items(), potentials_df.items()):
                        pop_name = pop_rate[0]
                        rate_trajectory = pop_rate[1]
                        potential_trajectory = pop_potential[1]
                        pop_df = pd.DataFrame()
                        pop_df['population'] = [pop_name]*len(rate_trajectory)
                        pop_df['rate'] = rate_trajectory
                        pop_df['potential'] = potential_trajectory
                        pop_df['time'] = rate_trajectory.index.values * 1e-3 # time in s
                        pop_df['coupling_strength'] = g
                        pop_df['InputDuration'] = d
                        pop_df['InputStrength'] = s 
                        time_data_df = pd.concat([time_data_df, pop_df])

                if load_full_potentials:
                    # load full 3D (target x source x timestep) potential of selected population 
                    filename = f"full_potentials_G{g}_{cortex_type}_Iduration{d}_{input_type}Istrength{s}_Ionset{input_onset}_tauVisual_{thalamus_source}.csv"
                    potential_df = pd.read_csv(os.path.join(output_dir, filename), sep=',', header=None)
                    potential_df = potential_df.applymap(lambda x: ast.literal_eval(x) if isinstance(x, str) and x.startswith('[') and x.endswith(']') else x)
                    potential_df = potential_df.iloc[load_population_potential]
                    longterm_pot_df = pd.DataFrame()

                    longterm_input = []
                    for sourcepop in potential_df:
                        # extract the longterm behaviour
                        longterm_input.append(sourcepop[start_sample:stop_sample])
                    
                    longterm_input = np.sum(np.array(longterm_input), axis=0)
                    longterm_mean = np.sum(longterm_input)
                    print(longterm_input.shape)

                    longterm_pot_df['lt_potential'] = longterm_input
                    longterm_pot_df['population'] = load_population_potential
                    longterm_pot_df['coupling_strength'] = g
                    longterm_pot_df['InputDuration'] = d
                    longterm_pot_df['InputStrength'] = s 
                    potential_data_df = pd.concat([potential_data_df, longterm_pot_df]) 

    
    return summary_df, time_data_df, potential_data_df