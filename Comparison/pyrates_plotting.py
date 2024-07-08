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
pyrates_11_copy = pd.read_csv('output/pyrates_11_copy.csv')
pyrates_11_copy_ij = pd.read_csv('output/pyrates_11_copy_ij.csv')
pyrates_11_copy_ji = pd.read_csv('output/pyrates_11_copy_ji.csv')

python_10 = pd.read_csv('output/python_10.csv') 
pyrates_10 = pd.read_csv('output/pyrates_10.csv')

python_4 = pd.read_csv('output/python_4.csv')
pyrates_4 = pd.read_csv('output/pyrates_4.csv')


# %%
fig, ax = plt.subplots(1,1)
time_r = np.linspace(0,2,len(python_10)+1)
time = np.linspace(0,2,len(python_10))
ax.set_title('11 cells')
ax.plot(time, python_11['S1'][:], label='python', color = 'b', alpha = 0.5)
ax.plot(time_r, pyrates_11['S1'][:], label='pyrates', color = 'r', alpha= 0.5)
#ax.plot(time_r, pyrates_11_copy['S1'][:], label= 'pyrates_copy', color = 'g', alpha = 0.9)
#ax.plot(time_r, pyrates_11_copy_ij['S1'][:], label= 'pyrates_copy_ij', color = 'magenta', alpha = 0.9)
#ax.plot(time_r, pyrates_11_copy_ji['S1'][:], label= 'pyrates_copy_ji', color = 'orange', alpha = 0.9)
ax.legend()
  

# %%
time_r = np.linspace(0,2,len(python_10)+1)
time = np.linspace(0,2,len(python_10))

fig, ax = plt.subplots(4,1)
fig.suptitle('Difference between PyRates and pure Python for example of S1')
ax[0].set_title('16 cells')
ax[0].plot(time, python['S1'][:], label='python')
ax[0].plot(time_r, pyrates['S1'][:], label='pyrates')

ax[1].set_title('15 cells')
ax[1].plot(time, python_15['S1'][:], label='python')
ax[1].plot(time_r, pyrates_15['S1'][:], label='pyrates')

ax[2].set_title('11 cells')
ax[2].plot(time, python_11['S1'][:], label='python', color = 'b', alpha = 0.5)
ax[2].plot(time_r, pyrates_11['S1'][:], label='pyrates', color = 'r', alpha= 0.5)
ax[2].plot(time_r, pyrates_11_copy['S1'][:], label= 'pyrates_copy', color = 'g', alpha = 0.9)
ax[2].legend()
  
ax[3].set_title('10 cells')
ax[3].plot(time, python_10['S1'][:], label='python')
ax[3].plot(time_r, pyrates_10['S1'][:], label='pyrates')

# ax[4].set_title('4')
# ax[4].plot(python_4['E3'][:], label='python')
# ax[4].plot(pyrates_4['E3'][:], label='pyrates')
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.legend()
plt.show()

#Unterschied zwischen 10 & 11 bei E1, E3, E3, P1, P2, P3, S1, S2
#Unterschied erst zwischen 11 & 16 bei E4, P4

# %%
cells = ['E1', 'E2', 'E3', 'E4', 'P1', 'P2', 'P3', 'P4', 'S1', 'S2']

fig, ax = plt.subplots(5, 2, figsize=(10, 20))
fig.suptitle('Total cell potentials compared between PyRates(--) and pure Python(-) with 10, 11, and 16 populations')

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

    ax[0, 0].legend(loc='upper right')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.subplots_adjust(hspace=0.8)
plt.savefig('output/compare.pdf') 
plt.show()