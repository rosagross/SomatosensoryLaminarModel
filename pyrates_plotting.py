# %%
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd 


python = pd.read_csv('output/python_up_syn.csv')
pyrates = pd.read_csv('output/pyrates_up_syn.csv')
# %%
plt.plot(python[:][:500], label='python')
plt.plot(pyrates[:][:500], label='pyrates')

plt.legend()
plt.show()

# %%
