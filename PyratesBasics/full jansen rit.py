# %% 
import numpy as np
import matplotlib.pyplot as plt
import time

start_time = time.time()

def sigmoid(v):
    m_max = 2.5 # in s^-1! (so 1 means 1 per second)
    r = 0.56 # sigmoidal steepness in V^-1 (560 here means 0.56 mV^-1)
    v_thr = 6 # v0

    return 2*m_max / (1 + np.exp(r*(v_thr - v)))


def jrc(C): 
    # Simulation parameters
    simulation_time = 2 # in s
    dt = 1e-4 # in s

    # Model parameters
    H_e = 3.25 # A
    H_i = 22 # B  # Woher kommt das Inhibitorische wenn H nicht negativ ist??
    tau_e = 0.01 
    tau_i = 0.02 

    # Output matrices to store computed values

    # define time steps (starting with the second step since the first one is initialized already)
    steps = np.arange(0, simulation_time, dt)

    C1 = C
    C2 = 0.8 * C
    C3 = C4 = 0.25 * C
    C2i = 0.2 * C

    # Initialize first values with 0 or randomly
    v_p = 0
    u_p = 0
    v_e = 0
    u_e = 0
    v_i = 0
    u_i= 0

    potential = np.zeros((len(steps)))

    Iext = np.random.uniform(120, 320, size= len(steps)) 

    for t in range(len(steps)):
        # DERIVATES
        v_p = v_p + u_p * dt
        v_e = v_e + u_e * dt
        v_i = v_i + u_i * dt

        u_p = u_p + (H_e/tau_e * sigmoid(v_e - v_i)  - 2 * 1/tau_e * u_p - 1/(tau_e**2) * v_p) *dt
        u_e = u_e + (H_e * 1/tau_e * (Iext[t] + C2 * sigmoid(C1 * v_p)) - 2*1/tau_e*u_e - 1/tau_e**2 * v_e) *dt
        u_i= u_i + (H_i * 1/tau_i * (C4 * sigmoid(C3 * v_p)) - 2*1/tau_i*u_i - 1/tau_i**2 * v_i) *dt

        

        potential[t] = v_p # u_p sieht irgendwie besser aus

    plt.plot(steps, potential)
    plt.show()

    return potential


def bifurcation_diagram(C_values):
    bifurcation_points = []
    for c in C_values:
        pot = jrc(c)
        for i in range(1,100):
            val = pot[-i]
            bifurcation_points.append((c, val))
    return np.array(bifurcation_points)


# C variieren
C_values =  [68, 128, 135, 270, 675, 1350] 
bifurcation_data = bifurcation_diagram(C_values)

# # Diagramm erstellen
# plt.scatter(bifurcation_data[:,0], bifurcation_data[:,1], s=0.5, color='black')
# plt.xlabel('C')
# plt.ylabel('Potential')
# plt.show()


#Computation time
end_time = time.time()
execution_time = end_time - start_time
print("Die Ausführungszeit beträgt: ", execution_time, " Sekunden")
