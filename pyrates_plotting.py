# %%
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd 


potential_G1 = pd.read_csv('output/potentials_G1_visual.csv')
pyrates_potential_G1 = pd.read_csv('output/pyrates_potential_G1.csv')
potential_G1_W1 = pd.read_csv('output/potentials_G1_visual_W1.csv')
pyrates_potential_G1_W1 = pd.read_csv('output/pyrates_potential_G1_W1.csv')

# %%
plt.plot(potential_G1['E2'][:500], label='python')
plt.plot(pyrates_potential_G1['E2'][:500], label='pyrates')
#plt.plot(potential_G1_W1['E2'][:500], label='python W1')
#plt.plot(pyrates_potential_G1_W1['E2'][:500], label='pyrates W1')
plt.legend()
plt.show()

# %%
