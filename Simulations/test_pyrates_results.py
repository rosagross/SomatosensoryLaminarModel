"""
File: test_pyrates_a3b.py
Author: Rosa Grossmann
Contact: grossmannr@cbs.mpg.de
Date: 2026-01-08
Description: Write tests if the current model is corresponding to the results simulated with pyrates. 
- extracted area 3b from big python model =? small area 3b pyrates model
- complete model test 
    - try different cases here (different input strengths, g, bEI, etc.)

"""

import csv
import os
import numpy as np
from numpy import testing
import matplotlib.pyplot as plt
import pandas as pd
from somato_model import SomatoModel, read_simulation_params
import plotting_functions as pf


SIMDIR = os.getenv("SIMDIR")
WDDIR = os.getenv("WDDIR")
param_path = os.path.join(WDDIR, 'Simulations')
figure_dir = os.path.join(SIMDIR, "Figures")


def test_pyrates_ThalA3b():

    params = read_simulation_params()

    # set parameters 
    params['coupling_strength'] = 10
    params['balance_EI'] = 0.5
    params['Iext_duration'] = 0
    params['Iext_strength'] = 0
    params['Ib_strength'] = 0
    params['area'] = 'ThalA3b'
    
    # additional parameters (that are usually fixed)
    params['g_thal'] = 2
    params['bEI_thal'] = 0
    params['extI_cellcounts'] = 1000
    params['bI_cellcounts'] = 100
    params['thal_cellcounts'] = 500

    # setup model
    model = SomatoModel(params)
    
    # simulate rates and potentials
    model.simulate()

    # extract results for area A3b and thal
    rate = model.rate
    potential_sum = np.sum(model.potential, axis=1)
    potential_ThalA3b = np.concat((potential_sum[:4], potential_sum[-2:])).T

    # load pyrates results
    potential_sum_pyrates = pd.read_csv(os.path.join(param_path, 'test_data', 'pyrates_ThalA3b_bEI0.5_g10_gthal2_bEIthal0.5_thalcellcount500_test.csv'), index_col=False)
    potential_pyrates_array = np.array(potential_sum_pyrates)
    print(potential_pyrates_array.shape)

    plt.plot(potential_ThalA3b[:,0])
    plt.plot(potential_ThalA3b[:,1])
    plt.plot(potential_ThalA3b[:,2])
    plt.plot(potential_ThalA3b[:,3])
    legend_list = []
    for i in range(4):
        legend_list.append(f'{i} {np.round(potential_ThalA3b[:,i][-1],10)}')
    plt.legend(legend_list)
    plt.title('Python - Thal A3b')
    plt.show()

    plt.plot(potential_pyrates_array[:,0]-potential_ThalA3b[:,0])
    plt.plot(potential_pyrates_array[:,1]-potential_ThalA3b[:,1])
    plt.plot(potential_pyrates_array[:,2]-potential_ThalA3b[:,2])
    plt.plot(potential_pyrates_array[:,3]-potential_ThalA3b[:,3])
    plt.title('Thal & Area 3b  - Difference Pyrates vs. Pure Python')
    plt.show()

    # check if they are (almost) equal
    testing.assert_almost_equal(potential_ThalA3b, potential_pyrates_array)

    #assert potential_sum.equals(potential_sum_pyrates), "Results don't match!!"



def test_pyrates_A1():

    params = read_simulation_params()

    # set parameters 
    params['coupling_strength'] = 10
    params['balance_EI'] = 0.5
    params['Iext_duration'] = 0
    params['Iext_strength'] = 0
    params['Ib_strength'] = 0
    params['area'] = 'A1'

    # setup model
    model = SomatoModel(params)
    
    # simulate rates and potentials
    model.simulate()

    # extract results for area A3b and thal
    rate = model.rate
    potential_sum = np.sum(model.potential, axis=1)
    potential_A1 = potential_sum[4:17].T

    # load pyrates results
    potential_sum_pyrates = pd.read_csv(os.path.join(param_path, 'test_data', 'pyrates_A1_test.csv'), index_col=False)
    potential_pyrates_array = np.array(potential_sum_pyrates)
    print(potential_pyrates_array.shape)
    plt.plot(potential_pyrates_array[:,0]-potential_A1[:,0])
    plt.plot(potential_pyrates_array[:,1]-potential_A1[:,1])
    plt.plot(potential_pyrates_array[:,2]-potential_A1[:,2])
    plt.plot(potential_pyrates_array[:,3]-potential_A1[:,3])
    plt.plot(potential_pyrates_array[:,4]-potential_A1[:,4])
    plt.plot(potential_pyrates_array[:,5]-potential_A1[:,5])
    plt.plot(potential_pyrates_array[:,6]-potential_A1[:,6])
    plt.plot(potential_pyrates_array[:,7]-potential_A1[:,7])
    plt.plot(potential_pyrates_array[:,8]-potential_A1[:,8])
    plt.plot(potential_pyrates_array[:,9]-potential_A1[:,9])
    plt.plot(potential_pyrates_array[:,10]-potential_A1[:,10])
    plt.plot(potential_pyrates_array[:,11]-potential_A1[:,11])
    plt.title('Area 1 - Difference Pyrates vs. Pure Python')
    plt.show()

    # check if they are (almost) equal
    testing.assert_almost_equal(potential_A1, potential_pyrates_array, 5)



def test_pyrates_ThalA1():
    """
    Testing different parameter settings.
    Pyrates data is computed beforehand and saved. Python model is running live.
    """

    params = read_simulation_params()

    test_cases = [[10, 0.5, 0, 0, 0, 50, 1, 1000, 100, 500],
                  [10, 0.5, 0, 0, 0, 2, 0.5, 1000, 100, 500],
                  [10, 0.5, 0.5, 50, 5, 50, 1, 1000, 100, 500]]

    case = 1
    
    # set parameters 
    params['coupling_strength'] = test_cases[case][0]
    params['balance_EI'] = test_cases[case][1]
    params['Iext_duration'] = test_cases[case][2]
    params['Iext_strength'] = test_cases[case][3]
    params['Ib_strength'] = test_cases[case][4]
    params['area'] = 'ThalA1'
    
    # additional parameters (that are usually fixed)
    params['g_thal'] = test_cases[case][5]
    params['bEI_thal'] = test_cases[case][6]
    params['extI_cellcounts'] = test_cases[case][7]
    params['bI_cellcounts'] = test_cases[case][8]
    params['thal_cellcounts'] = test_cases[case][9]

    test_filename = f'pyrates_ThalA1_bEI{test_cases[case][1]}_g{test_cases[case][0]}_gthal{test_cases[case][5]}_bEIthal{test_cases[case][6]}_thalcellcount500_Iext{test_cases[case][3]}_Idur{test_cases[case][2]}_bI{test_cases[case][4]}_test.csv'
    # setup model
    model = SomatoModel(params)
    
    # simulate rates and potentials
    model.simulate()

    # extract results for area A3b and thal
    rate = model.rate
    potential_sum = np.sum(model.potential, axis=1)
    potential_A1 = np.concat((potential_sum[4:17], potential_sum[-2:])).T

    # load pyrates results
    potential_sum_pyrates = pd.read_csv(os.path.join(param_path, 'test_data', test_filename), index_col=False)
    cells = potential_sum_pyrates.columns
    potential_pyrates_array = np.array(potential_sum_pyrates)
    print(potential_pyrates_array.shape)
    plt.plot(potential_pyrates_array[:,0]-potential_A1[:,0])
    plt.plot(potential_pyrates_array[:,1]-potential_A1[:,1])
    plt.plot(potential_pyrates_array[:,2]-potential_A1[:,2])
    plt.plot(potential_pyrates_array[:,3]-potential_A1[:,3])
    plt.plot(potential_pyrates_array[:,4]-potential_A1[:,4])
    plt.plot(potential_pyrates_array[:,5]-potential_A1[:,5])
    plt.plot(potential_pyrates_array[:,6]-potential_A1[:,6])
    plt.plot(potential_pyrates_array[:,7]-potential_A1[:,7])
    plt.plot(potential_pyrates_array[:,8]-potential_A1[:,8])
    plt.plot(potential_pyrates_array[:,9]-potential_A1[:,9])
    plt.plot(potential_pyrates_array[:,10]-potential_A1[:,10])
    plt.plot(potential_pyrates_array[:,11]-potential_A1[:,11])
    plt.plot(potential_pyrates_array[:,12]-potential_A1[:,12])
    plt.plot(potential_pyrates_array[:,13]-potential_A1[:,13])
    plt.title('Thal & Area 1 - Difference Pyrates vs. Pure Python')
    plt.legend(cells)
    plt.show()

    # check if they are (almost) equal
    testing.assert_almost_equal(potential_A1, potential_pyrates_array, 5)

def test_complete_python_model():

    input_onset = 1.001
    input_type = 'step'
    step_size = 0.01
    simulation_time = 2


    params = read_simulation_params()

    test_cases = [[1, 1, 0, 0, 0, 0, 0, 1000, 100, 500]] # tens!
                  #[0, 0, 0, 0, 0, 0, 0, 1000, 100, 500], # fine
                  #[10, 0, 0, 0, 0, 0, 0, 1000, 100, 500], # fine
                  #[10, 1, 0, 0, 0, 0, 0, 1000, 100, 500], # error
                  #[10, 0.5, 0, 0, 0, 0, 0, 1000, 100, 500], # error!
                  #[10, 0, 0, 0, 0, 50, 0, 1000, 100, 500], # fine
                  #[10, 0, 0, 0, 0, 50, 0.5, 1000, 100, 500], # fine
                  #[10, 0, 0, 0, 0, 50, 1, 1000, 100, 500]] # fine!

    
    for case in test_cases:
        # set parameters 
        params['coupling_strength'] = case[0]
        params['balance_EI'] = case[1]
        params['Iext_duration'] = case[2]
        params['Iext_strength'] = case[3]
        params['Ib_strength'] = case[4]
        params['area'] = 'all'
        
        # additional parameters (that are usually fixed)
        params['g_thal'] = case[5]
        params['bEI_thal'] = case[6]
        params['extI_cellcounts'] = case[7]
        params['bI_cellcounts'] = case[8]
        params['thal_cellcounts'] = case[9]

        # setup model
        model = SomatoModel(params)
        
        # simulate rates and potentials
        model.simulate()

        rate = model.rate
        potential_sum_python = np.sum(model.potential, axis=1).T

        # load pyrates results
        test_filename = f'TENSA3b_gthal{case[5]}_bEIthal{case[6]}_g{case[0]}_bEI{case[1]}_Ib{case[4]}_Iextd{case[2]}_stepIexts{case[3]}_Ionset1.001_thalcells500_Ibcells100_Iextcells1000_PYRATES.hdf5'
        potential_sum_pyrates = pd.read_hdf(os.path.join(param_path, 'test_data', test_filename), key='summed_potential')

        # difference 
        diff_df = pd.DataFrame(
        potential_sum_python - potential_sum_pyrates,
        index=potential_sum_pyrates.index,
        columns=potential_sum_pyrates.columns
        )

        # plot difference 
        cells = potential_sum_pyrates.columns
        for cell in cells[:4]:
            plt.plot(diff_df[cell], label=cell)
        plt.legend()
        plt.show()

        for cell in cells[4:17]:
            plt.plot(diff_df[cell], label=cell)
        plt.legend()
        plt.show()

        for cell in cells[17:-2]:
            plt.plot(diff_df[cell], label=cell)
        plt.legend()
        plt.show()

        for cell in cells[-2:]:
            plt.plot(diff_df[cell], label=cell)
        plt.legend()
        plt.show()

        python_df = pd.DataFrame(potential_sum_python,
        index=potential_sum_pyrates.index,
        columns=potential_sum_pyrates.columns
        )
        pyrates_df = pd.DataFrame(potential_sum_pyrates,
        index=potential_sum_pyrates.index,
        columns=potential_sum_pyrates.columns
        )

        for cell in cells[-2:]:
            plt.plot(pyrates_df[cell], label=f'{cell} - pyrates')
        plt.legend()
        plt.show()

        for cell in cells[-2:]:
            plt.plot(python_df[cell], label=cell)
        plt.legend()
        plt.show()

        start_plot = 0

        # assert if different
        #testing.assert_almost_equal(potential_sum_python, potential_sum_pyrates, 8)


    



def test_complete_pyrates_model():
    """
    Compare to previously computed python results to make sure they sill match.
    """


#test_pyrates_ThalA3b()
#test_pyrates_A1()
#test_pyrates_ThalA1()
test_complete_python_model()