# %%
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd 


potential_G1 = pd.read_csv('output/potentials_G1_W1.csv')
pyrates_potential_G1 = pd.read_csv('output/pyrates_potential_G1_W1.csv')


# %%
plt.plot(potential_G1['P1'], label='python')
plt.plot(pyrates_potential_G1['P1'][:1000], label='pyrates')
plt.legend()
plt.show()

# %%
