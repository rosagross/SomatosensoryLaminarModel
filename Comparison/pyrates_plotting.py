# %%
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd 


python = pd.read_csv('output/python_upscaled.csv')
pyrates = pd.read_csv('output/pyrates_upscaled.csv')

python_15 = pd.read_csv('output/python_15.csv')
pyrates_15 = pd.read_csv('output/pyrates_15.csv')

python_11 = pd.read_csv('output/python_11.csv')
pyrates_11 = pd.read_csv('output/pyrates_11.csv')

python_10 = pd.read_csv('output/python_10.csv') 
pyrates_10 = pd.read_csv('output/pyrates_10.csv')

python_4 = pd.read_csv('output/python_4.csv')
pyrates_4 = pd.read_csv('output/pyrates_4.csv')
# %%
fig, ax = plt.subplots(4,1)
ax[0].set_title('16')
ax[0].plot(python['S1'][:], label='python')
ax[0].plot(pyrates['S1'][:], label='pyrates')

ax[1].set_title('15')
ax[1].plot(python_15['S1'][:], label='python')
ax[1].plot(pyrates_15['S1'][:], label='pyrates')

ax[2].set_title('11')
ax[2].plot(python_11['S1'][:], label='python')
ax[2].plot(pyrates_11['S1'][:], label='pyrates')
  
ax[3].set_title('10')
ax[3].plot(python_10['S1'][:], label='python')
ax[3].plot(pyrates_10['S1'][:], label='pyrates')

# ax[4].set_title('4')
# ax[4].plot(python_4['E3'][:], label='python')
# ax[4].plot(pyrates_4['E3'][:], label='pyrates')
plt.legend()
plt.show()

#Unterschied zwischen 10 & 11 bei E1, E3, E3, P1, P2, P3, S1, S2
#Unterschied erst zwischen 11 & 16 bei E4, P4

