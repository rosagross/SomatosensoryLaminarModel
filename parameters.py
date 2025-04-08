# %%
import numpy as np
import yaml

# %%
class Parameter():
    
    def __init__(self):
        self.tau, self.nPop = self.get_params()
        self.S = self.get_connectStrength()
        self.P = self.get_connectProb()
        self.C = self.get_cellcounts()
        self.sigmoid_params = self.get_sigmoid()

    def get_params(self):

        # nr. of populations
        nPopA3b = 3 # E2, PV2, SOM2 
        nPopS1 = 13 #  
        nPopS2 = 13 # 
        nPopThal = 2 # one forward excitatory and one inhibitory feedback neuron
        nPopTotal = nPopS1+nPopS2+nPopA3b+nPopThal
        # SYNAPTIC DECAY (depends on the connection type excitatory/inhibitory)
        # E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1
        #tau = np.tile(np.array([2,2,2,2,4,4,4,4,4,4,4,4,4,3])*1e-3, (nPop+1,1)) 
        #tau = np.tile(np.array([5.2, 5.2, 5.9, 5.9, 3, 3, 3.8, 3.8, 11.2, 11.2, 11.1, 11.1, 10.4])*1e-3, (nPop+1,1)) 
        # Visual cortex values
        # the last two values are used for the external input and background input
        # with nPopTotal = 28 the shape of tau should be (28, 32) 
        
        tauA3b = np.tile(np.array([6,3,20])*1e-3, (nPopTotal,1)) # sec
        tauS1 = np.tile(np.array([6,3,20,15,6,3,20,6,3,20,6,3,20])*1e-3, (nPopTotal,1)) # sec
        tauS2 = np.tile(np.array([6,3,20,15,6,3,20,6,3,20,6,3,20,3,3,3,3])*1e-3, (nPopTotal,1)) # sec
        tau = np.hstack((tauA3b,tauS1,tauS2))

        # TODO: MEMBRANE CONSTANT (NOTE: this value depends on post synaptic neuron whereas the synaptic decay depends on the presynapse)
        
        return tau, nPopTotal

    def get_connectProb(self):

        # Connection probabilies
        # Targets in rows, sources in columns 
        # Area 3b: E2, PV2, SST2
        
        #P_A3b = np.array([[9.9, 37.4, 19.8],[37.4, 28.8, 36.3],[19., 33.2,  1.2]])*1e-2 
        
        # E1, PV1, SST1, VIP1, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4 (Thalamus is added later!)
        P_S1 = np.array([[ 6.7, 27.1, 28.,  4.3, 11.,  2.,  4.5,  2.4,  0.1,  1.5,  0., 0, 0],
                    [28.5, 31.8, 11.8,  5.5,  0.4,  1.7,  3.2,  0.05,  0.1,  1.4,  0.01, 0.01, 0.6],
                    [24.3, 29.1,  0., 27.1,  0.4, 1.2, 4.4, 0.03, 0.1, 2.1, 0., 0., 0],
                    [16., 20.6, 46.3,  6.3,  0.2,  1.,  1.7,  0.02,  0.1,  0.6, 0, 0.01,  0.5],
                    [ 1.3,  3.5,  5.9,  7.,  9.9, 37.4, 19.8,  0.6,  2.8,  5.2, 0, 0.2, 2.4],
                    [ 3.6,  1.6,  2.3,  5.3, 37.4, 28.8, 36.3,  0.9,  2.8,  4.7,  0.1, 0.3, 1.8],
                    [ 3.6,  1.8,  2.4,  4.9, 19., 33.2,  1.2,  0.9,  3.2,  4.7,  0.2, 0.2, 2.3],
                    [ 8.6,  4.2,  9.4,  6.2,  8.7,  7.,  7.6,  9.8, 39.6, 15.1,  1., 1.8, 5.1],
                    [ 1.7,  0.5,  1.1,  1.9,  3.5,  3.,  2.3, 39.6, 24.6, 24.1,  0.4, 1.3, 2.6],
                    [ 1.6,  0.4,  1.,  1.9,  3.8,  3.0,  2.5, 20.5, 31.1,  0.6,  0.6, 1-.5,  3.2],
                    [0, 0.2, 0.3, 0.8, 3.0, 1.5, 1.0, 4.0, 1.9, 1.7, 2.1, 39.6, 15.1],
                    [0.3, 0.01, 0.02, 0.2, 0.5, 0.2, 0.05, 1.6, 0, 0.3, 39.6, 24.6, 24.1],
                    [0.2, 0, 0.02, 0.1, 0.5, 0.2, 0.04, 1.5, 0, 0.2, 20.5, 31.1, 0.6]])*1e-2      

        P_S2 = np.array([[ 6.7, 27.1, 28.,  4.3, 11.,  2.,  4.5,  2.4,  0.1,  1.5,  0., 0, 0],
                    [28.5, 31.8, 11.8,  5.5,  0.4,  1.7,  3.2,  0.05,  0.1,  1.4,  0.01, 0.01, 0.6],
                    [24.3, 29.1,  0., 27.1,  0.4, 1.2, 4.4, 0.03, 0.1, 2.1, 0., 0., 0],
                    [16., 20.6, 46.3,  6.3,  0.2,  1.,  1.7,  0.02,  0.1,  0.6, 0, 0.01,  0.5],
                    [ 1.3,  3.5,  5.9,  7.,  9.9, 37.4, 19.8,  0.6,  2.8,  5.2, 0, 0.2, 2.4],
                    [ 3.6,  1.6,  2.3,  5.3, 37.4, 28.8, 36.3,  0.9,  2.8,  4.7,  0.1, 0.3, 1.8],
                    [ 3.6,  1.8,  2.4,  4.9, 19., 33.2,  1.2,  0.9,  3.2,  4.7,  0.2, 0.2, 2.3],
                    [ 8.6,  4.2,  9.4,  6.2,  8.7,  7.,  7.6,  9.8, 39.6, 15.1,  1., 1.8, 5.1],
                    [ 1.7,  0.5,  1.1,  1.9,  3.5,  3.,  2.3, 39.6, 24.6, 24.1,  0.4, 1.3, 2.6],
                    [ 1.6,  0.4,  1.,  1.9,  3.8,  3.0,  2.5, 20.5, 31.1,  0.6,  0.6, 1-.5,  3.2],
                    [0, 0.2, 0.3, 0.8, 3.0, 1.5, 1.0, 4.0, 1.9, 1.7, 2.1, 39.6, 15.1],
                    [0.3, 0.01, 0.02, 0.2, 0.5, 0.2, 0.05, 1.6, 0, 0.3, 39.6, 24.6, 24.1],
                    [0.2, 0, 0.02, 0.1, 0.5, 0.2, 0.04, 1.5, 0, 0.2, 20.5, 31.1, 0.6]])*1e-2      
        
        # E population projects to S1 layer 4 E and PV 
        #P_A3btoS1 = np.array([[0,0,0],[0,0,0],[0,0,0],[0,0,0],
        #                        [100,0,0],[50,0,0],[0,0,0],
        #                        [0,0,0],[0,0,0],[0,0,0],
        #                        [0,0,0],[0,0,0],[0,0,0]])*1e-2
        
        #P_A3btoS2 = np.zeros((13,3)) 
        P_S1toS2 = np.zeros((13,13)) # TODO: implement S1 to S2 connection! 
        P_S1toA3b = np.zeros((3,13)) # TODO: implement feedback S1 to A3b! 
        P_S2toS1 = np.zeros((13,13)) # TODO: implement feedback from S2 to S1
        P_S2toA3b = np.zeros((3,13)) 
        #P_toA3b = np.hstack((P_A3b,P_S1toA3b,P_S2toA3b))
        P_toS1 = np.hstack((P_S1, P_S2toS1))
        P_toS2 = np.hstack((P_S1toS2, P_S2))
        P = np.vstack((P_toS1, P_toS2)) #  26x26 (with nPop=13 per S1/S2)
    
        return P


    def get_connectStrength(self):

        # Based on values below 
        #psp_A3b = [[1.59, 1.0, 1.0],[1.5, 1.0, 1.0],[0.5, 1.0, 0.19]]
                    
        # Postsynaptic potential (13x13) averages from Isbister, Jiang, and more (see excel file FinalConnectivity_PSP.csv for mor info)
        # order: E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        # Targets in rows, sources in columns
        psp_S1 = np.array([
            [0.75, 1.0, 1.0, 1.0, 0.5, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 1.0, 1.0],
            [0.75, 1.0, 1.0, 1.0, 1.2, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.25, 1.0, 1.0, 1.0, 0.3, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 2.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 1.59, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 1.5, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 0.5, 1.0, 0.19, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 0.6, 1.1, 0.38, 0.8, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.25, 1.0, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 1.0, 1.0],
            [0.25, 1.0, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 1.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 0.29, 1.0, 1.0, 0.5, 1.0, 1.0, 0.75, 2.0, 1.0],
            [0.25, 1.0, 1.0, 1.0, 0.25, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [0.25, 1.0, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 2.0, 1.0]])

        psp_S2 = [
            [0.75, 1.0, 1.0, 1.0, 0.5, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 1.0, 1.0],
            [0.75, 1.0, 1.0, 1.0, 1.2, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.25, 1.0, 1.0, 1.0, 0.3, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 2.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 1.59, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 1.5, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 0.5, 1.0, 0.19, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 0.6, 1.1, 0.38, 0.8, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.25, 1.0, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 1.0, 1.0],
            [0.25, 1.0, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 1.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 0.29, 1.0, 1.0, 0.5, 1.0, 1.0, 0.75, 2.0, 1.0],
            [0.25, 1.0, 1.0, 1.0, 0.25, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [0.25, 1.0, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 1.0, 1.0, 0.25, 2.0, 1.0]]

        # to Area 3b
        #psp_S1toA3b = [[0.5, 1.0, 1.0, 1.0, 1.59, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
        #    [0.5, 1.0, 1.0, 1.0, 1.5, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
        #    [0.5, 1.0, 1.0, 1.0, 0.5, 1.0, 0.19, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
        #    ] # same as S1 to layer 4, since A3b is only modelled as the input layer!
        #psp_S2toA3b = np.zeros((3,13))
        #psp_toA3b = np.hstack((psp_A3b, psp_S1toA3b, psp_S2toA3b))
        
        # to S1
        psp_S2toS1 = np.zeros((13,13))
        #psp_A3btoS1 = psp_S1[:,4:7] # same as layer 4 of S1
        psp_toS1 = np.hstack((psp_S1, psp_S2toS1))
        
        # to S2
        #psp_A3btoS2 = psp_S1[:,4:7] # same as layer 4 of S1, but 3b to S2 connection is unlikely
        # TODO: insert connectivity between area S1 and S2 
        psp_S1toS2 = np.zeros((13,13))
        psp_toS2 = np.hstack((psp_S1toS2, psp_S2))
        psp = np.vstack((psp_toS1, psp_toS2)) #  26x26 (with nPop=13 per S1/S2)

        return psp

    def get_cellcounts(self):
        
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        C_S1 = np.array([1691, 90, 74, 85, 1656, 85, 48, 1095, 109, 105, 1288, 56, 66])
        CS1_total = np.sum(C_S1)
        C_relativeS1 = C_S1/CS1_total
        C_S2 = np.array([1691, 90, 74, 85, 1656*(3/4), 85, 48, 1095*(5/4), 109, 105, 1288, 56, 66]) # for S2 adapt this to more neurons in layer 5 than in layer 4
        C_relativeS2 = C_S2/np.sum(C_S2)

        # For Area 3b we take the cell counts summed for each population
        # So we end up with a "collapsed" area 
        C_A3b_E = np.sum(C_relativeS1[[0,4,7,10]])
        C_A3b_PV = np.sum(C_relativeS1[[1,5,8,11]])
        C_A3b_SST = np.sum(C_relativeS1[[2,6,9,12]])
        C_A3b_VIP = np.sum(C_relativeS1[[3]])
        C = np.hstack((C_relativeS1, C_relativeS2))

        return C 

    def get_connectivity(self, gE, gI, include_Iext=True):
        # g is a scaling factor scaling the general coupling strength

        S = self.get_connectStrength()
        P = self.get_connectProb()
        C = self.get_cellcounts()

        # Final connectivity matrix 
        # PS = P * S (16 x 16) --> Connectivity probability if all cells were equally distributed
        # W0 = PS * C (16 x 16) --> make sure that cell counts are accounted for! C is shaped (16 x 1) so tile it before multiplying  
        PS = P * S
        nPop = len(PS)
        W0 = PS * np.tile(C, (nPop,1))

        idx_E = [0,4,7,10]
        idx_P = [1,5,8,11]
        idx_S = [2,6,9,12]
        idx_V = [2,6,9,12]

        # set up the connectivity for area 3b as a weighted average of the connectivity from S1
        # ... within A3b
        # ... merge input from all E to all other populations
        W_A3bA3b_E = np.sum(W0[:,idx_E], axis=1)  
        W_A3bA3b_EE = np.dot(W_A3bA3b_E[idx_E], C[idx_E]) # .. now we have to compute the weighted average for the all receiving cells
        W_A3bA3b_PE = np.dot(W_A3bA3b_E[idx_P], C[idx_P])
        W_A3bA3b_SE = np.dot(W_A3bA3b_E[idx_S], C[idx_S])
        W_A3bA3b_VE = np.dot(W_A3bA3b_E[idx_V], C[idx_V])
        # ... merge input from all P targeting all other populations
        W_A3bA3b_P = np.sum(W0[:,idx_P], axis=1)  
        W_A3bA3b_EP = np.dot(W_A3bA3b_P[idx_E], C[idx_E])
        W_A3bA3b_PP = np.dot(W_A3bA3b_P[idx_P], C[idx_P])
        W_A3bA3b_SP = np.dot(W_A3bA3b_P[idx_S], C[idx_S])
        W_A3bA3b_VP = np.dot(W_A3bA3b_P[idx_V], C[idx_V])
        # ... merge input from all S targeting all other populations
        W_A3bA3b_S = np.sum(W0[:,idx_S], axis=1)  
        W_A3bA3b_ES = np.dot(W_A3bA3b_S[idx_E], C[idx_E])
        W_A3bA3b_PS = np.dot(W_A3bA3b_S[idx_P], C[idx_P])
        W_A3bA3b_SS = np.dot(W_A3bA3b_S[idx_S], C[idx_S])
        W_A3bA3b_VS = np.dot(W_A3bA3b_S[idx_V], C[idx_V])
        # ... merge input from all V targeting all other populations
        W_A3bA3b_V = np.sum(W0[:,idx_V], axis=1)  
        W_A3bA3b_EV = np.dot(W_A3bA3b_V[idx_E], C[idx_E])
        W_A3bA3b_PV = np.dot(W_A3bA3b_V[idx_P], C[idx_P])
        W_A3bA3b_SV = np.dot(W_A3bA3b_V[idx_S], C[idx_S])
        W_A3bA3b_VV = np.dot(W_A3bA3b_V[idx_V], C[idx_V])

        # stack them all together to create the A3b to A3b connection matrix
        W_A3bA3b = np.hstack((np.vstack((W_A3bA3b_EE, W_A3bA3b_PE, W_A3bA3b_SE, W_A3bA3b_VE)),
                            np.vstack((W_A3bA3b_EP, W_A3bA3b_PP, W_A3bA3b_SP, W_A3bA3b_VP)),
                            np.vstack((W_A3bA3b_ES, W_A3bA3b_PS, W_A3bA3b_SS, W_A3bA3b_VS)),
                            np.vstack((W_A3bA3b_EV, W_A3bA3b_PV, W_A3bA3b_SV, W_A3bA3b_VV))))

        # ... from  area 3b to S1



        # ... from S1 to area 3b

        
   
        # create the external input matrix, based on thalamus connectivity (average of findings, see FinalConnectivity_PSPs.ods)
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        # S_thal = np.array([0.49, 0.49, 0.245, 0, 1.45, 2.3, 0.245, 0.5, 0.49,0.245, 0.85, 2.2, 0.245]) # Thal to S1 based on Isbister & jiang 

        # E2, PV2, SST2
        # input to area 3b is collapsed S1 (PS weighted average by cell count of each population)
        S_thalToA3b_E = PS
        #S_thalToA3b_PV = 
        #S_thalToA3b_SST = 

        S_thalToA3b = np.array([[0.49, 0.49, 0.245], # copied from Thal to S1 (Jiang et al 2023 only)
                                [0, 0, 0]])
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        S_thalToS1 = np.array([[0.49, 0.49, 0.245, 0, 0.49, 0.49, 0.245, 0.49, 0.49, 0.245, 0.49, 0.49, 0.245], # Thal to S1: Based on Jiang et al. 2023 only! 
                              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) # no thal E/I target (integrated in last two entries of below S2 array) 
        S_thalToS2 = np.zeros((2,15)) # TODO: check if there should still be some connection from thalamus to S1 (area 1)!
        S_thalToS2[0, -1] = 0.5 # ThalE to ThalI
        S_thalToS2[1, -2] = 0.5 # ThalI to ThalE
                    #np.array([[0.49, 0.49, 0.245, 0, 0.49, 0.49, 0.245, 0.49, 0.49, 0.245, 0.49, 0.49, 0.245, 0, 0.5], # last two are Thal E and Thal I 
                    #           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0]]) # Reticular inhibition: just an assumption
        
        #print('S', S_thalToS1.shape, S_thalToS2.shape)
        S_from_thal = np.hstack((S_thalToA3b, S_thalToS1, S_thalToS2))
        
        S_A3bto_thal = np.array([[0, 0, 0], [0, 0, 0]])
        S_S1to_thal = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
        S_S2to_thal = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
        S_to_thal = np.hstack((S_A3bto_thal, S_S1to_thal,S_S2to_thal))

        #print('S', S_from_thal.shape)
        
        # cell count of the thalamus 
        C_thal = [1, 1] # Thal E, Thal I
        
        # connection probabilities
        # E2, PV2, SOM2
        P_thalToA3b = np.array([[40, 40, 20], # thalamus excitatory
                                [0, 0, 0]])*1e-2 # reticlar nucleus inhibitory (estimation)
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        P_thalToS1 = np.array([[6.2, 6.2, 0, 0, 40, 40, 20, 25.9, 25.9, 0, 9, 9, 0], # thalamus excitatory
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])*1e-2 # reticlar nucleus inhibitory (estimation)

        # TODO: the last two values are ThalE and ThalI and the value 10 is intuitively assigned. Check if this makes sense!
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        P_thalToS2 = np.array([[6.2, 6.2, 0, 0, 40, 40, 20, 25.9, 25.9, 0, 9, 9, 0, 0, 10], # thalamus excitatory
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0]])*1e-2 # reticlar nucleus inhibitory (estimation)
        
        
        # TODO: Check how much S2 receives from the thalamus! In reality S2 receives way less than S1    
        P_thalToS2 = P_thalToS2 * 0.2 # S2 receives just a fith of what S1 receives from the thalamus
        #print('P from thal', P_thalToS1.shape, P_thalToS2.shape)
        P_from_thal = np.hstack((P_thalToA3b, P_thalToS1, P_thalToS2))
        
        # order: E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        P_A3bto_thal = np.array([[0, 0, 0], [0, 0, 0]]) # this is only from the cortex!
        P_S1to_thal = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) # this is only from the cortex!
        P_S2to_thal = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) # this is only from the cortex!
        P_to_thal = np.hstack((P_A3bto_thal, P_S1to_thal,P_S2to_thal))
        
        #print('P to thal', P_to_thal.shape)

        # calculate final thalamus connectivity
        PS_to_thal = np.multiply(P_to_thal, S_to_thal)
        W_to_thal = np.multiply(PS_to_thal, C)
        PS_from_thal = np.multiply(P_from_thal, S_from_thal)
        W_from_thal = np.multiply(PS_from_thal, np.expand_dims(C_thal,1))

        # weight the inhibitory popultaion (from the reticular nucleus in the thalamus) as negative
        W_from_thal[1] = W_from_thal[1]*-1

        # Create the external input matrix
        # Only the thalamus E population receives the external input (from the brain stem)
        Wext = np.zeros((W_from_thal.shape[1],1))
        Wext[29] = 1 # thalamus E population
        Wext[30] = 0 # reticular I population

        # .. and also for the background input (all cells receive input except the thalamus)
        Wb = np.zeros((W_from_thal.shape[1],1))
        Wb[:-2] = 1

        # make inhibitory connections negative and apply weights gI and gE respectively
        idx_I_A3b = np.array([1,2])
        idx_I_S = [1,2,3,5,6,8,9,11,12]
        idx_I_S1 = np.array(idx_I_S)+3 
        idx_I_S2 = np.array(idx_I_S)+16
        idx_I = np.concatenate((idx_I_A3b,idx_I_S1,idx_I_S2))
        idx_E_A3b = np.array([0])
        idx_E_S = [0,4,7,10]
        idx_E_S1 = np.array(idx_E_S)+3 
        idx_E_S2 = np.array(idx_E_S)+16
        idx_E = np.concatenate((idx_E_A3b,idx_E_S1,idx_E_S2))

        W0[:,idx_I] = W0[:,idx_I] * -gI # negative weight for inhibitory connections
        W0[:,idx_E] = W0[:,idx_E] * gE

        # include the external input to the matrix 
        if include_Iext:
            print('append input connectivity')
            # append the thalamus population(s) values to the matrix
            W0 = np.append(W0, W_to_thal, axis=0)
            print('W from thalamus 2\n', W_from_thal.shape)
            print(W_from_thal)

            W0 = np.append(W0, W_from_thal.T, axis=1)

            W = np.concatenate((W0, Wb, Wext), axis=1)
            print('\nW matrix shape', W.shape) 
            print('W matrix\n', W[:,-4:]) 

        return W

    def get_sigmoid(self):
        # sigmoid function (nPop x 3) --> 3 stands for parameters: r(1/mV), v_thr(mV), m_max (1/s)
        sigmoid_params_A3b = np.array([[  0.12782346,  32.10540543,  31.39696397],
                                    [  0.14218422,  40.03107351, 166.82960408],
                                    [  0.07937015,  42.01276379,  56.95305832]])
        
        # order:
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        sigmoid_params_S1 = np.array([[  0.12782346,  32.10540543,  31.39696397],
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
                                    [  0.07937015,  42.01276379,  56.95305832]]) 

        sigmoid_params_S2 = np.array([[  0.12782346,  32.10540543,  31.39696397],
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
                                    [  0.07937015,  42.01276379,  56.95305832],
                                    [  0.1,  40,  30], # Thalamus E
                                    [  0.1,  40,  30]]) # Thalamus I

        sigmoid_params = np.vstack((sigmoid_params_A3b, sigmoid_params_S1, sigmoid_params_S2))
                                       
        
            
        return sigmoid_params
    
    def save_to_yaml(self, filename, gE, gI):
        
        S = self.get_connectStrength()
        P = self.get_connectProb()
        C = self.get_cellcounts()
        W = self.get_connectivity(gE, gI)

        # Convert numpy arrays to lists
        parameters = {
            'S': S.tolist(),
            'P': P.tolist(),
            'C': C.tolist(),
            'W': W.tolist()
        }

        # Save parameters to a YAML file
        with open(filename + '.yaml', 'w') as file:
            yaml.dump(parameters, file)


# %%
