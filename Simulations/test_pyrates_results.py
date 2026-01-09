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

SIMDIR = os.getenv("SIMDIR")
WDDIR = os.getenv("WDDIR")
param_path = os.path.join(WDDIR, 'Simulations')


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
    Until now I only test one case. More should be implemented
    """


    params = read_simulation_params()

    # set parameters 
    params['coupling_strength'] = 10
    params['balance_EI'] = 0.5
    params['Iext_duration'] = 0
    params['Iext_strength'] = 0
    params['Ib_strength'] = 0
    params['area'] = 'ThalA1'
    
    # additional parameters (that are usually fixed)
    params['g_thal'] = 2
    params['bEI_thal'] = 1
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
    potential_A1 = np.concat((potential_sum[4:17], potential_sum[-2:])).T

    # load pyrates results
    potential_sum_pyrates = pd.read_csv(os.path.join(param_path, 'test_data', 'pyrates_ThalA1_bEI0.5_g10_gthal2_bEIthal1_thalcellcount500_test.csv'), index_col=False)
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

def test_complete_model():
    params = read_simulation_params()

    # set parameters 
    params['coupling_strength'] = 10
    params['balance_EI'] = 0.5
    params['Iext_duration'] = 0
    params['Iext_strength'] = 0
    params['Ib_strength'] = 0
    params['area'] = 'all'

    # setup model
    model = SomatoModel(params)
    
    # simulate rates and potentials
    model.simulate()

    # extract results for area A3b and thal
    rate = model.rate
    potential_sum = np.sum(model.potential, axis=1)
    potential_A1 = np.concat((potential_sum[4:17], potential_sum[-2:])).T

    # load pyrates results
    potential_sum_pyrates = pd.read_csv(os.path.join(param_path, 'test_data', 'pyrates_complete_test.csv'), index_col=False)
    potential_pyrates_array = np.array(potential_sum_pyrates)



#test_pyrates_ThalA3b()
#test_pyrates_A1()
test_pyrates_ThalA1()
