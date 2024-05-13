# %% 

'''
In this file I am experimenting with a bifurcation analysis of the Jansen-Rit model.
Structure:
- reading in simulation files 
'''
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from plotting import figure_style


colors, _ = figure_style()

# %% 

test = pd.read_csv(f'../PyratesBasics/Outputs/simpleJR_rates_tauE{0.1}_tauI{0.02}.csv')

# %% Reading in simulation files 
tau_Es = [0.0005, 0.001, 0.01, 0.1, 0.2]
tau_I = 0.02
jr_simulations = []

# plot them in a subplot
fig, axes = plt.subplots(5)

for tau_E, ax in zip(tau_Es, axes):

    simulation_data = pd.read_csv(f'../PyratesBasics/Outputs/simpleJR_rates_tauE{tau_E}_tauI{tau_I}.csv')
    # ToDo: get the max and min membrane potential value for every population

    ax.plot(simulation_data)
    ax.set_ylabel('Hz')
    ax.legend(['PC', 'IIN', 'EIN'])

    # Set titles for each subplot
    ax.set_title(f'tau E = {tau_E}')
    
    #jr_simulations.append(, sep=',', header=None))

sns.despine(fig, trim=True, bottom=False)
plt.tight_layout(h_pad=1)

# %% Plotting phase diagram 

# x-axis: 