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

# %%
import matplotlib.pyplot as plt

cells = ['E1', 'E2', 'E3', 'E4', 'P1', 'P2', 'P3', 'P4', 'S1', 'S2']

fig, ax = plt.subplots(5, 2, figsize=(10, 20))
fig.suptitle('PyRates(--) vs pure Python(-) with 10, 11 and 16 populations')

time_r = np.linspace(0,2,len(python_10)+1)
time = np.linspace(0,2,len(python_10))

for i, cell in enumerate(cells):
    row = i // 2
    col = i % 2
    ax[row, col].set_title(cell)

    ax[row, col].plot(time,python_10[cell][:], label='Python 10', color='darkblue', alpha=0.6)
    ax[row, col].plot(time_r,pyrates_10[cell][:], label='Pyrates 10', color='darkblue', ls='--', alpha=0.6)
    ax[row, col].plot(time,python_11[cell][:], label='Python 11', color='darkgreen', alpha=0.8)
    ax[row, col].plot(time_r,pyrates_11[cell][:], label='Pyrates 11', color='darkgreen', ls='--', alpha=0.8)
    ax[row, col].plot(time,python[cell][:], label='Python 16', color='darkorange', alpha=0.4)
    ax[row, col].plot(time_r,pyrates[cell][:], label='Pyrates 16', color='darkorange', ls='--', alpha=0.4)

    ax[row,col].set_xlabel('time (s)')
    ax[row,col].set_ylabel('potential (mV)') 

    ax[0, 0].legend()

plt.tight_layout(h_pad=7)
plt.savefig('output/compare.pdf') 
plt.show()