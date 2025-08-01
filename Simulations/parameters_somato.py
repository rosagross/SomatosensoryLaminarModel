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
        # nr. of populations
        nPop = 13
        tau = np.tile(np.array([6, 3, 10, 15, 6, 3, 20, 6, 3, 20, 6, 3, 20,3])*1e-3, (nPop+1,1)) # sec

        # synaptic kernel efficacy
        return tau, nPop

    def get_connectProb(self):

        # Connection probabilies
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        P = np.array([[ 6.7, 27.1, 28.,  4.3, 11.,  2.,  4.5,  2.4,  0.1,  1.5,  0., 0, 0,  6.2 ],
                      [28.5, 31.8, 11.8,  5.5,  0.4,  1.7,  3.2,  0.05,  0.1,  1.4,  0.01, 0.01, 0.6,  6.2 ],
                      [24.3, 29.1,  0., 27.1,  0.4, 1.2, 4.4, 0.03, 0.1, 2.1, 0., 0., 0, 0],
                      [16., 20.6, 46.3,  6.3,  0.2,  1.,  1.7,  0.02,  0.1,  0.6, 0, 0.01,  0.5, 0],
                      [ 1.3,  3.5,  5.9,  7.,  9.9, 37.4, 19.8,  0.6,  2.8,  5.2, 0, 0.2, 2.4, 40.],
                      [ 3.6,  1.6,  2.3,  5.3, 37.4, 28.8, 36.3,  0.9,  2.8,  4.7,  0.1, 0.3, 1.8, 40.],
                      [ 3.6,  1.8,  2.4,  4.9, 19., 33.2,  1.2,  0.9,  3.2,  4.7,  0.2, 0.2, 2.3, 20.],
                      [ 8.6,  4.2,  9.4,  6.2,  8.7,  7.,  7.6,  9.8, 39.6, 15.1,  1., 1.8, 5.1, 25.9],
                      [ 1.7,  0.5,  1.1,  1.9,  3.5,  3.,  2.3, 39.6, 24.6, 24.1,  0.4, 1.3, 2.6, 25.9],
                      [ 1.6,  0.4,  1.,  1.9,  3.8,  3.0,  2.5, 20.5, 31.1,  0.6,  0.6, 1-.5,  3.2, 0],
                      [0, 0.2, 0.3, 0.8, 3.0, 1.5, 1.0, 4.0, 1.9, 1.7, 2.1, 39.6, 15.1, 9.0],
                      [0.3, 0.01, 0.02, 0.2, 0.5, 0.2, 0.05, 1.6, 0, 0.3, 39.6, 24.6, 24.1, 9.0],
                      [0.2, 0, 0.02, 0.1, 0.5, 0.2, 0.04, 1.5, 0, 0.2, 20.5, 31.1, 0.6, 0]
                      ])

        return P


    def get_connectStrength(self):

        # The postsynaptic currents are computed using synaptical time constants and timeconstant for the membrane
        # synaptical time constants (all incoming EPSPs and all IPSPS have the same)
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        tau_syn = np.tile(np.array([2, 4, 4, 4, 2, 4, 4, 2, 4, 4, 2, 4, 4]),(self.nPop, 1))

        # (13x13)
        # timeconstant for the membrane (time until the potential reaches equilibrium), TARGETS
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4, Thalamus
        tau_m = np.tile(np.array([5.2, 3.0, 11.2, 10.4, 5.2, 3.0, 11.2, 5.9, 3.8, 11.1, 5.9, 3.8, 11.1]),(self.nPop, 1)).T

        # (13x13)
        a = (tau_syn/tau_m)
        Cm = np.tile(np.array([229.8, 93.9, 123.3, 86.5, 229.8, 93.9, 123.3, 269.2, 81.0, 146.8, 269.2, 81.0, 146.8]),(self.nPop, 1)).T

        # Postsynaptic potential (13x14), input from: E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4, Thalamus
        psp = np.array([[0.25, -2, -2, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.49], # E1 
                        [0.75, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.49], # PV1
                        [0.5, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.245], # SST1
                        [0.25, -2.0, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0], # VIP
                        [0.5, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.49], # E2
                        [0.75, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.49], # PV2
                        [0.5, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.245], # SST2
                        [0.5, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.49], # E3
                        [0.25, -2.0, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.49], # PV3
                        [0.25, -2.0, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.245], # SST3
                        [0.5, -2.0, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -2.0, -2.0, 0.5, -3.0, -2.0, 0.49], # E4
                        [0.25, -2.0, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.49], # PV4
                        [0.25, -2.0, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -2.0, -2.0, 0.25, -3.0, -2.0, 0.245], # SST4                        
                        ])
    
        
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        # psc is a (13x14) matrix
        psc = (Cm*(a-1)) * psp[:, :-1] / (tau_syn * ((a**(1/(1-a))) - a**(a/(1-a))))

        # set connection where probability is zero to 0 
        psc[0, 10:] = 0
        psc[2, 10:] = 0
        psc[2, 2] = 0
        psc[3, 10] = 0
        psc[4, 10] = 0
        psc[10, 0] = 0
        psc[11, 8] = 0
        psc[12, 8] = 0
        psc[12, 1] = 0

        return psc

    def get_cellcounts(self):
        # Cell counts (only VIP in first layer!)
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        C_absolute = np.array([1691, 90, 74, 85, 1656, 85, 48, 1095, 109, 105, 1288, 56, 66])
        
        # translate in proportion
        C = C_absolute/np.sum(C_absolute)

        return C 

    def get_connectivity(self, g):
        # g is a scaling factor scaling the general coupling strength

        # Final connectivity matrix 
        # PS = P * S (16 x 16) --> Connectivity probability if all cells were equally distributed
        # W = PS * C (16 x 16) --> make sure that cell counts are accounted for! C is shaped (16 x 1) so tile it before multiplying  
        PS = self.P[:, :13] * self.S
        W = PS * np.tile(self.C, (self.nPop,1))

        # sort so that this order yields:  E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1, VIP2, VIP3, VIP4, Iext
        # define index for E, PV, SOM, and VIP
        '''
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
        '''
        # now append the external input matrix 
        Wext = np.zeros((self.nPop, 1))
        Wext[1] = 1

        W = np.concatenate((W*g, Wext), axis=1)

        return W


    def get_sigmoid(self):
        # sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
        sigmoid_params = np.array([[  0.12782346,  32.10540543,  31.39696397],
                                    [  0.14218422,  40.03107351, 166.82960408],
                                    [  0.07937015,  42.01276379,  56.95305832],
                                    [  0.0704119 ,  37.86409387,  38.52689646],
                                    [  0.12782346,  32.10540543,  31.39696397],
                                    [  0.14218422,  40.03107351, 166.82960408],
                                    [  0.07937015,  42.01276379,  56.95305832],
                                    [  0.12782346,  32.10540543,  31.39696397],
                                    [  0.14218422,  40.03107351, 166.82960408],
                                    [  0.07937015,  42.01276379,  56.95305832],
                                    [  0.12782346,  32.10540543,  31.39696397],
                                    [  0.14218422,  40.03107351, 166.82960408],
                                    [  0.07937015,  42.01276379,  56.95305832],])

        return sigmoid_params
    


