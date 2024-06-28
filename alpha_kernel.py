
# %% 
import matplotlib.pyplot as plt
import numpy as np


# %% Alpha kernel visualizations

tau = 0.2
H = 2
t = np.arange(0, 1, 0.01)
alpha_kernel = H * t/tau * np.exp(-t/tau)

plt.plot(t, alpha_kernel)
plt.show()
