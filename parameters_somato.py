import numpy as np

class Parameter():

    def __init__(self):
        self.tau, self.nPop = self.get_params()
        self.S = self.get_connectStrength()
        self.P = self.get_connectProb()
        self.C = self.get_cellcounts()
        self.sigmoid_params = self.get_sigmoid()

    def get_params(self):

        # time constants
        # E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1, VIP2, VIP3, VIP4, Iext
        tau = np.tile(np.array([])*1e-3, (17,1)) # sec

        # nr. of populations
        nPop = 16

        # synaptic kernel efficacy
        return tau, nPop

    def get_connectProb(self):

        # Connection probabilies
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        P = np.array([[]])
        return P


    def get_connectStrength(self):

        # The postsynaptic currents are computed using the following formula

        # synaptical time constants (all E and all I have the same)
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        tau_syn = np.array([2, 4, 4, 4, 2, 4, 4, 2, 4, 4, 2, 4, 4])
        # timeconstant for the membrane (time until the potential reaches equilibrium)
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        tau_m = np.array([5.2, 3.0, 11.2, 10.4, 5.2, 3.0, 11.2, 5.9, 3.8, 11.1, 5.9, 3.8, 11.1])

        a = tau_syn/tau_m
        Cm = np.array([229.8, 93.9, 123.3, 86.5, 229.8, 93.9, 123.3, 269.2, 81.0, 146.8, 269.2, 81.0, 146.8])
        psp = np.array([[0.25, -2, -2, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0], # E1 
                        [0.75, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0], # PV1
                        [0.5, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0], # SST1
                        [0.25, -2.0, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0], # VIP
                        [0.5, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0], # E2
                        [0.75, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0], # PV2
                        [0.5, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0], # SST2
                        [0.5, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0], # E3
                        [0.25, -2.0, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0], # PV3
                        [0.25, -2.0, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0] # SST3
                        [0.5, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -3.0, -2.0], # E4
                        [0.25, -2.0, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0], # PV4
                        [0.25, -2.0, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -3.0, -2.0], # SST4                        
                        ])
        
        # Thalamic input to E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        thalamic_psp = np.array([])
        
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        # psc is a (13x13) matrix 
        psc = (Cm*(a-1)*psp) / (tau_syn (a^(1/(1-a))) - a^(a/(1-a)))



        # Synaptic strength
        S = np.array([])	
        return S

    def get_cellcounts(self):
        # Cell counts (only VIP in first layer!)
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        C = np.array([1691, 90, 74, 85, 1656, 85, 48, 1095, 109, 105, 1288, 56, 66])
        # TODO: Scale the values so that they are between 0 and 1!

        return C 

    def get_connectivity(self, g):
        # g is a scaling factor scaling the general coupling strength

        # Final connectivity matrix 
        # PS = P * S (16 x 16) --> Connectivity probability if all cells were equally distributed
        # W = PS * C (16 x 16) --> make sure that cell counts are accounted for! C is shaped (16 x 1) so tile it before multiplying  
        PS = self.P * self.S
        W = PS * np.tile(self.C, (16,1))

        # sort so that this order yields:  E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1, VIP2, VIP3, VIP4, Iext
        # define index for E, PV, SOM, and VIP
        iE = np.arange(0, 16, 4)  # E1, E2, E3, E4
        iP = np.arange(1, 16, 4)  # PV1, PV2, PV3, PV4
        iS = np.arange(2, 16, 4)  # SOM1, SOM2, SOM3, SOM4
        iV = np.arange(3, 16, 4)  # VIP1, VIP2, VIP3, VIP4

        # Extracting submatrices based on the defined index sets
        Wee = W[np.ix_(iE, iE)]
        Wpe = W[np.ix_(iP, iE)]
        Wse = W[np.ix_(iS, iE)]
        Wve = W[np.ix_(iV, iE)]

        Wep = W[np.ix_(iE, iP)]
        Wpp = W[np.ix_(iP, iP)]
        Wsp = W[np.ix_(iS, iP)]
        Wvp = W[np.ix_(iV, iP)]

        Wes = W[np.ix_(iE, iS)]
        Wps = W[np.ix_(iP, iS)]
        Wss = W[np.ix_(iS, iS)]
        Wvs = W[np.ix_(iV, iS)]

        Wev = W[np.ix_(iE, iV)]
        Wpv = W[np.ix_(iP, iV)]
        Wsv = W[np.ix_(iS, iV)]
        Wvv = W[np.ix_(iV, iV)]

        # put them back together
        # Reconstructing W0 from the submatrices, negating where indicated
        row1 = np.hstack([Wee, -Wep, -Wes, -Wev])
        row2 = np.hstack([Wpe, -Wpp, -Wps, -Wpv])
        row3 = np.hstack([Wse, -Wsp, -Wss, -Wsv])
        row4 = np.hstack([Wve, -Wvp, -Wvs, -Wvv])

        # Vertically stack the rows to form W0
        W0 = np.vstack([row1, row2, row3, row4])

        # now append the external input matrix 
        Wext = np.zeros((16,1))
        Wext[1] = 1
        W = np.concatenate((W0*g, Wext), axis=1)

        return W


    def get_sigmoid(self):
        # sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
        sigmoid_params = np.array([[]])

        return sigmoid_params
    


