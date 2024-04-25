# %%
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd 


potential_G1 = pd.read_csv('output/potentials_G1.csv')
pyrates_potential_G1 = pd.read_csv('output/pyrates_potential_G1.csv')


# %%
plt.plot(potential_G1['P1'])
plt.plot(pyrates_potential_G1['P1'][:1000])

# %%
