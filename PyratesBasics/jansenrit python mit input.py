#steady state des Potentials in Abhängigkeit vom Input?

# %% 
import numpy as np
import matplotlib.pyplot as plt
import time

start_time = time.time()

def jrc(C): 
    # Simulation parameters
    simulation_time = 2 # in s
    step_size = 1e-4 # in s

    # Model parameters
    m_max = 2.5 # in s^-1! (so 1 means 1 per second)
    r = 0.56 # sigmoidal steepness in V^-1 (560 here means 0.56 mV^-1)
    v_thr = 6 # v0
    H_E = 3.25 # A
    H_I = -22 # B  # Woher kommt das Inhibitorische wenn H nicht negativ ist??
    tau_E = 0.01 
    tau_I = 0.02 

    # Output matrices to store computed values

    # define time steps (starting with the second step since the first one is initialized already)
    steps = np.arange(2*step_size, simulation_time+step_size, step_size)

    # rates & potentials (E, IIN , EIN) 
    nPop = 3
    rate = np.zeros((nPop, len(steps)))
    potential = np.zeros((nPop, nPop+1, len(steps))) 

    u_t = np.zeros((nPop, nPop+1)) # the initial first-order derivative: v'(t) = u(t)
    u_dot = np.zeros((nPop, nPop+1)) # we don't need an initial value for this one actually but a place holder array
    tau = np.array([tau_E, tau_I, tau_E, tau_E])
    H = np.array([H_E, H_I, H_E, H_E])

    C1 = C
    C2 = 0.8 * C
    C3 = C4 = 0.25 * C
    W = np.array([[0, C4, C2, 0],[C3, 0, 0, 0],[C1, 0, 0, 0]]) # PYR, I, E, Input (welches Weight?)

    # Simulation loop 

    # Initialize first values with 0 or randomly
    rate[:, 0] = 0
    potential[:, 0] = 0 #np.random.normal(0,1,(3,round(simulation_time/step_size)-1)) # noise

    Iext = np.random.uniform(120, 320, size=(nPop, len(steps))) # Und welches weight?

    for iter, time in enumerate(steps):

        # store current rate and potential for now
        rate_current = rate[:, iter-1] 
        v_current = potential[:, :, iter-1] 

        # DERIVATES
        # the current 
        # RPO: update the population's membrane potential (calculated using the previous rate): u'(t) = v''(t)
        for i in range(nPop):
            for j in range(nPop):

                if i == 0 and j == 2:
                    u_dot[i,j] = H[j]/tau[j] * (W[i,j]*(rate_current[j]+Iext[i,iter]))  - 2 * u_t[i, j]/tau[j] - v_current[i,j]/(tau[j]**2)
                else:
                    u_dot[i,j] = H[j]/tau[j] * (W[i,j]*rate_current[j])  - 2 * u_t[i, j]/tau[j] - v_current[i,j]/(tau[j]**2)


            #external Input durch RPO:
            #v_dot = u_t[i, -1]
            #v_current[i, -1] = v_current[i, -1] + v_dot * step_size
            #u_dot[i,-1] = H[-1]/tau[-1] * (W[i, -1] * Iext[i, iter]) - 2 * u_t[i, -1]/tau[-1] - v_current[i, -1]/(tau[-1]**2)
            #u_t[i, -1] = u_t[i,-1] + u_dot[i,-1] * step_size

        v_current = v_current + u_t * step_size
        u_t = u_t + u_dot * step_size

        # PRO: update the population's RATE (calculated from the calculated **updated** potential)
        rate_current = 2*m_max / (1 + np.exp(r*(v_thr - np.sum(v_current, axis=1)))) # All populations have the same PRO

        rate[:,iter] = rate_current
        potential[:, :, iter] = v_current

    plt.plot(potential[0, 1]+potential[0, 2], label='P')
    plt.plot(potential[0, 1], label='I')
    plt.plot(potential[0, 2], label='E')
    plt.legend()
    plt.show()

    PYR = potential[0, 1]+potential[0, 2]
    return PYR


def bifurcation_diagram(C_values):
    bifurcation_points = []
    for c in C_values:
        pot = jrc(c)
        for i in range(1,100):
            val = pot[-i]
            bifurcation_points.append((c, val))
    return np.array(bifurcation_points)


# C variieren
C_values =  [68, 128, 135, 270, 675, 1350] # np.linspace(20, 1200, 100)
bifurcation_data = bifurcation_diagram(C_values)

# Diagramm erstellen
plt.figure(figsize=(10, 6))
plt.scatter(bifurcation_data[:,0], bifurcation_data[:,1], s=0.5, color='black')
plt.xlabel('C')
plt.ylabel('Potential')
plt.show()


#Computation time
end_time = time.time()
execution_time = end_time - start_time
print("Die Ausführungszeit beträgt: ", execution_time, " Sekunden")