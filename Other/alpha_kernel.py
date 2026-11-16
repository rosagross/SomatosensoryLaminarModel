
# %% 
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# %% Alpha kernel visualizations

tau1 = 0.2
tau2 = 0.3
tau3 = 0.4
tau4 = 0.2
H = 1
t = np.arange(0, 2, 0.01)
alpha_kernel1 = H * t/tau1 * np.exp(-t/tau1)
alpha_kernel2 = 2 * t/tau2 * np.exp(-t/tau2)
alpha_kernel3 = H * t/tau3 * np.exp(-t/tau3)
alpha_kernel4 = H * t/tau4 * np.exp(-t/tau4)

plt.plot(t, alpha_kernel1, label="1")
plt.plot(t, alpha_kernel2, label="2")
plt.plot(t, alpha_kernel3, label="3")
plt.plot(t, alpha_kernel4, label="4")
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
#################### NMDA and AMPA vizualisation ####################################

# Time constants (milliseconds)
tau_ampa = 6.0
tau_nmda = 100.0

# Relative contributions
w_ampa = 1.0
w_nmda = 0.36

H = 1.0
t = np.arange(0, 500, 0.5)

# Individual alpha kernels
def alpha_kernel(t, H, tau):
    return H * (t / tau) * np.exp(-t / tau)

ampa = w_ampa * alpha_kernel(t, H, tau_ampa)
nmda = w_nmda * alpha_kernel(t, H, tau_nmda)

# Combined AMPA + NMDA response
combined = ampa + nmda

# Fit a single alpha kernel to the combined response
def effective_kernel(t, H_eff, tau_eff):
    return H_eff * (t / tau_eff) * np.exp(-t / tau_eff)

# Fit amplitude and time constant
popt, pcov = curve_fit(
    effective_kernel,
    t,
    combined
)

H_eff, tau_eff = popt
approximation = effective_kernel(t, H_eff, tau_eff)

# Calculate goodness of fit
r_squared = 1 - np.sum((combined - approximation)**2) / np.sum(
    (combined - np.mean(combined))**2
)

# Visualization
plt.figure(figsize=(10, 5))

plt.plot(t, ampa, '--', label=f'AMPA (tau={tau_ampa} ms)')
plt.plot(t, nmda, '--', label=f'NMDA (tau={tau_nmda} ms)')
plt.plot(t, combined, linewidth=2.5, label='AMPA + NMDA')
plt.plot(
    t, approximation, ':', linewidth=2.5,
    label=f'Single alpha approximation (tau={tau_eff:.2f} ms)'
)

plt.xlabel('Time (ms)')
plt.ylabel('Synaptic response')
plt.title(f'AMPA/NMDA vs. single alpha kernel (R²={r_squared:.4f})')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

print(f"Effective amplitude: {H_eff:.4f}")
print(f"Effective time constant: {tau_eff:.2f} ms")
print(f"R²: {r_squared:.4f}")

# %%
