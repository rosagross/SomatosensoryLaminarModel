
# %% 
import matplotlib.pyplot as plt
import numpy as np


# %% Alpha kernel visualizations

tau1 = 0.2
tau2 = 0.3
tau3 = 0.4
H = 2
t = np.arange(0, 2, 0.01)
alpha_kernel1 = H * t/tau1 * np.exp(-t/tau1)
alpha_kernel2 = H * t/tau2 * np.exp(-t/tau2)
alpha_kernel3 = H * t/tau3 * np.exp(-t/tau3)

plt.plot(t, alpha_kernel1, label="1")
plt.plot(t, alpha_kernel2, label="2")
plt.plot(t, alpha_kernel3, label="3")
plt.legend()
plt.show()


# %% Sigmoid curve visualization 


r = 0.62
v0 = 2
max_firing = 1
v = np.arange(20, 40, 0.01)
sigm = max_firing / (1 + r*np.exp((v0-v)))

plt.plot(v, sigm)
plt.show()

# %%
