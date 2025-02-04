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
        nPopS1 = 13 #  E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1
        nPopS2 = 13 #  E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1
        nPopThal = 2 # one forward excitatory and one inhibitory feedback neuron
        nPopTotal = nPopS1+nPopS2+nPopA3b+nPopThal
        # SYNAPTIC DECAY (depends on the connection type excitatory/inhibitory)
        #tau = np.tile(np.array([2,2,2,2,4,4,4,4,4,4,4,4,4,3])*1e-3, (nPop+1,1)) 
        #tau = np.tile(np.array([5.2, 5.2, 5.9, 5.9, 3, 3, 3.8, 3.8, 11.2, 11.2, 11.1, 11.1, 10.4])*1e-3, (nPop+1,1)) 
        # Visual cortex values
        # the last two values are used for the external input and background input
        # with nPopTotal = 28 the shape of tau should be (28, 32) 
        tauA3b = np.tile(np.array([6,3,20])*1e-3, (nPopTotal,1)) # sec
        tauS1 = np.tile(np.array([6,6,6,6,3,3,3,3,20,20,20,20,15])*1e-3, (nPopTotal,1)) # sec
        tauS2 = np.tile(np.array([6,6,6,6,3,3,3,3,20,20,20,20,15,3,3,3,3])*1e-3, (nPopTotal,1)) # sec
        tau = np.hstack((tauA3b,tauS1,tauS2))

        # TODO: MEMBRANE CONSTANT (NOTE: this value depends on post synaptic neuron whereas the synaptic decay depends on the presynapse)
        
        
        return tau, nPopTotal

    def get_connectProb(self):

        # Connection probabilies
        # Targets in rows, sources in columns 
        # Area 3b: E2, PV2, SST2
        P_A3b = np.array([[9.9, 37.4, 19.8],[37.4, 28.8, 36.3],[19., 33.2,  1.2]])
        
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
        P_A3btoS1 = np.array([[0,0,0],[0,0,0],[0,0,0],[0,0,0],
                                [100,0,0],[50,0,0],[0,0,0],
                                [0,0,0],[0,0,0],[0,0,0],
                                [0,0,0],[0,0,0],[0,0,0]])*1e-2
        
        P_A3btoS2 = np.zeros((13,3)) 
        P_S1toS2 = np.zeros((13,13)) # TODO: implement S1 to S2 connection! 
        P_S1toA3b = np.zeros((3,13)) # TODO: implement feedback S1 to A3b! 
        P_S2toS1 = np.zeros((13,13)) # TODO: implement feedback from S2 to S1
        P_S2toA3b = np.zeros((3,13)) 
        P_toA3b = np.hstack((P_A3b,P_S1toA3b,P_S2toA3b))
        P_toS1 = np.hstack((P_A3btoS1, P_S1, P_S2toS1))
        P_toS2 = np.hstack((P_A3btoS2, P_S1toS2, P_S2))
        P = np.vstack((P_toA3b, P_toS1, P_toS2)) #  26x26 (with nPop=13 per S1/S2)

    
        return P


    def get_connectStrength(self):

        # Based on values below 
        psp_A3b = [[1.59, 1.0, 1.0],[1.5, 1.0, 1.0],[0.5, 1.0, 0.19]]
                    
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
        psp_S1toA3b = [[0.5, 1.0, 1.0, 1.0, 1.59, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 1.5, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            [0.5, 1.0, 1.0, 1.0, 0.5, 1.0, 0.19, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0],
            ] # same as S1 to layer 4, since A3b is only modelled as the input layer!
        psp_S2toA3b = np.zeros((3,13))
        psp_toA3b = np.hstack((psp_A3b, psp_S1toA3b, psp_S2toA3b))
        # to S1
        psp_S2toS1 = np.zeros((13,13))
        psp_A3btoS1 = psp_S1[:,4:7] # same as layer 4 of S1
        psp_toS1 = np.hstack((psp_A3btoS1, psp_S1, psp_S2toS1))
        # to S2
        psp_A3btoS2 = psp_S1[:,4:7] # same as layer 4 of S1, but 3b to S2 connection is unlikely
        psp_S1toS2 = np.zeros((13,13))
        psp_toS2 = np.hstack((psp_A3btoS2, psp_S1toS2, psp_S2))
        psp = np.vstack((psp_toA3b, psp_toS1, psp_toS2)) #  26x26 (with nPop=13 per S1/S2)

        
        return psp

    def get_cellcounts(self):
        
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        C_S1 = np.array([1691, 90, 74, 85, 1656, 85, 48, 1095, 109, 105, 1288, 56, 66])
        C_absoluteS1 = C_S1/np.sum(C_S1)
        C_S2 = np.array([1691, 90, 74, 85, 1656*(3/4), 85, 48, 1095*(5/4), 109, 105, 1288, 56, 66]) # for S2 adapt this to more neurons in layer 5 than in layer 4
        C_absoluteS2 = C_S2/np.sum(C_S2)

        # For Area 3b we take the same cell counts as for S1 
        C_absoluteA3b = C_absoluteS1[4:7]
        C = np.hstack((C_absoluteA3b,C_absoluteS1,C_absoluteS2))

        return C 

    def get_connectivity(self, gE, gI, include_Iext=True):
        # g is a scaling factor scaling the general coupling strength

        S = self.get_connectStrength()
        P = self.get_connectProb()
        C = self.get_cellcounts()

        # Final connectivity matrix 
        # PS = P * S (16 x 16) --> Connectivity probability if all cells were equally distributed
        # W = PS * C (16 x 16) --> make sure that cell counts are accounted for! C is shaped (16 x 1) so tile it before multiplying  
        PS = P * S
        nPop = len(PS)
        W = PS * np.tile(C, (nPop,1))

        # create the external input matrix, based on thalamus connectivity (average of findings, see FinalConnectivity_PSPs.ods)
        # order: 'E1','E2','E3','E4','PV1','PV2','PV3','PV4','SST1','SST2','SST3','SST4','VIP1'
        #S_thal = np.array([0.49, 1.45, 0.5, 0.85, 0.49, 2.3, 0.49, 2.2, 0.245, 0.245, 0.245, 0.245, 0]) # Thal to S1 based on Isbister & jiang 
        S_thalToA3b = np.array([[0.49, 0.49, 0.245], # Thal to S1: Based on Jiang et al. 2023 only! 
                                [0, 0, 0]])
        S_thalToS1 = np.array([[0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.245, 0.245, 0.245, 0.245, 0], # Thal to S1: Based on Jiang et al. 2023 only! 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) # no thal E/I target (integrated in last two entries of below S2 array) 
        S_thalToS2 = np.array([[0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.245, 0.245, 0.245, 0.245, 0, 0, 0.5], # last to are Thal E and Thal I 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0]]) # Reticular inhibition: just an assumption
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
        # order: 'E1','E2','E3','E4','PV1','PV2','PV3','PV4','SST1','SST2','SST3','SST4','VIP1'
        P_thalToS1 = np.array([[6.2, 40, 25.9, 9, 6.2, 40, 25.9, 9, 0, 20, 0, 0, 0], # thalamus excitatory
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])*1e-2 # reticlar nucleus inhibitory (estimation)
    
        # TODO: change this!! In reality S2 receives way less than S1    
        P_thalToS2 = np.array([[6.2, 40, 25.9, 9, 6.2, 40, 25.9, 9, 0, 20, 0, 0, 0, 0, 0], # thalamus excitatory
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])*1e-2 # reticlar nucleus inhibitory (estimation)
        
        P_thalToS2 = P_thalToS2 * 0.2 # S2 receives just a fith of what S1 receives from the thalamus
        #print('P from thal', P_thalToS1.shape, P_thalToS2.shape)
        P_from_thal = np.hstack((P_thalToA3b, P_thalToS1, P_thalToS2))
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

        # S1
        # indices to reorder the matrix to E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1
        iE_S1 = np.array([3, 7, 10, 13])  # E1, E2, E3, E4 of S1 
        iP_S1 = iE_S1+1  # PV1, PV2, PV3, PV4
        iS_S1 = iE_S1+2  # SOM1, SOM2, SOM3, SOM4
        iV_S1 = [6]  # VIP1

        # S2
        iE_S2 = np.array([16, 20, 23, 26])  # E1, E2, E3, E4 of S2 
        iP_S2 = iE_S2+1  # PV1, PV2, PV3, PV4
        iS_S2 = iE_S2+2  # SOM1, SOM2, SOM3, SOM4
        iV_S2 = [19]  # VIP1

        # Extracting submatrices based on the defined index sets
        Wee_A3b = W[np.ix_(0, iE_S1)]*gE
        Wpe_A3b = W[np.ix_(0, iE_S1)]*gI
        Wse_A3b = W[np.ix_(0, iE_S1)]*gI
        
        # within S1 and S2
        Wee_S1 = W[np.ix_(iE_S1, iE_S1)]*gE
        Wpe_S1 = W[np.ix_(iP_S1, iE_S1)]*gE
        Wse_S1 = W[np.ix_(iS_S1, iE_S1)]*gE
        Wve_S1 = W[np.ix_(iV_S1, iE_S1)]*gE
        Wee_S2 = W[np.ix_(iE_S2, iE_S2)]*gE
        Wpe_S2 = W[np.ix_(iP_S2, iE_S2)]*gE
        Wse_S2 = W[np.ix_(iS_S2, iE_S2)]*gE
        Wve_S2 = W[np.ix_(iV_S2, iE_S2)]*gE
        # connectivity between S1 and S2
        Wee_S1S2 = W[np.ix_(iE_S1, iE_S2)]*gE
        Wpe_S1S2 = W[np.ix_(iP_S1, iE_S2)]*gE
        Wse_S1S2 = W[np.ix_(iS_S1, iE_S2)]*gE
        Wve_S1S2 = W[np.ix_(iV_S1, iE_S2)]*gE
        Wee_S2S1 = W[np.ix_(iE_S2, iE_S1)]*gE
        Wpe_S2S1 = W[np.ix_(iP_S2, iE_S1)]*gE
        Wse_S2S1 = W[np.ix_(iS_S2, iE_S1)]*gE
        Wve_S2S1 = W[np.ix_(iV_S2, iE_S1)]*gE

        # within S1 and S2
        Wep_S1 = W[np.ix_(iE_S1, iP_S1)]*gI
        Wpp_S1 = W[np.ix_(iP_S1, iP_S1)]*gI
        Wsp_S1 = W[np.ix_(iS_S1, iP_S1)]*gI
        Wvp_S1 = W[np.ix_(iV_S1, iP_S1)]*gI
        Wep_S2 = W[np.ix_(iE_S2, iP_S2)]*gI
        Wpp_S2 = W[np.ix_(iP_S2, iP_S2)]*gI
        Wsp_S2 = W[np.ix_(iS_S2, iP_S2)]*gI
        Wvp_S2 = W[np.ix_(iV_S2, iP_S2)]*gI
        # connectivity between S1 and S2
        Wep_S1S2 = W[np.ix_(iE_S1, iP_S2)]*gI
        Wpp_S1S2 = W[np.ix_(iP_S1, iP_S2)]*gI
        Wsp_S1S2 = W[np.ix_(iS_S1, iP_S2)]*gI
        Wvp_S1S2 = W[np.ix_(iV_S1, iP_S2)]*gI
        Wep_S2S1 = W[np.ix_(iE_S2, iP_S1)]*gI
        Wpp_S2S1 = W[np.ix_(iP_S2, iP_S1)]*gI
        Wsp_S2S1 = W[np.ix_(iS_S2, iP_S1)]*gI
        Wvp_S2S1 = W[np.ix_(iV_S2, iP_S1)]*gI

        # within S1 and S2
        Wes_S1 = W[np.ix_(iE_S1, iS_S1)]*gI
        Wps_S1 = W[np.ix_(iP_S1, iS_S1)]*gI
        Wss_S1 = W[np.ix_(iS_S1, iS_S1)]*gI
        Wvs_S1 = W[np.ix_(iV_S1, iS_S1)]*gI
        Wes_S2 = W[np.ix_(iE_S2, iS_S2)]*gI
        Wps_S2 = W[np.ix_(iP_S2, iS_S2)]*gI
        Wss_S2 = W[np.ix_(iS_S2, iS_S2)]*gI
        Wvs_S2 = W[np.ix_(iV_S2, iS_S2)]*gI
        # connectivity between S1 and S2
        Wes_S1S2 = W[np.ix_(iE_S1, iS_S2)]*gI
        Wps_S1S2 = W[np.ix_(iP_S1, iS_S2)]*gI
        Wss_S1S2 = W[np.ix_(iS_S1, iS_S2)]*gI
        Wvs_S1S2 = W[np.ix_(iV_S1, iS_S2)]*gI
        Wes_S2S1 = W[np.ix_(iE_S2, iS_S1)]*gI
        Wps_S2S1 = W[np.ix_(iP_S2, iS_S1)]*gI
        Wss_S2S1 = W[np.ix_(iS_S2, iS_S1)]*gI
        Wvs_S2S1 = W[np.ix_(iV_S2, iS_S1)]*gI

        # within S1 and S2
        Wev_S1 = W[np.ix_(iE_S1, iV_S1)]*gI
        Wpv_S1 = W[np.ix_(iP_S1, iV_S1)]*gI
        Wsv_S1 = W[np.ix_(iS_S1, iV_S1)]*gI
        Wvv_S1 = W[np.ix_(iV_S1, iV_S1)]*gI
        Wev_S2 = W[np.ix_(iE_S2, iV_S2)]*gI
        Wpv_S2 = W[np.ix_(iP_S2, iV_S2)]*gI
        Wsv_S2 = W[np.ix_(iS_S2, iV_S2)]*gI
        Wvv_S2 = W[np.ix_(iV_S2, iV_S2)]*gI
        # connectivity between S1 and S2
        Wev_S1S2 = W[np.ix_(iE_S1, iV_S2)]*gI
        Wpv_S1S2 = W[np.ix_(iP_S1, iV_S2)]*gI
        Wsv_S1S2 = W[np.ix_(iS_S1, iV_S2)]*gI
        Wvv_S1S2 = W[np.ix_(iV_S1, iV_S2)]*gI
        Wev_S2S1 = W[np.ix_(iE_S2, iV_S1)]*gI
        Wpv_S2S1 = W[np.ix_(iP_S2, iV_S1)]*gI
        Wsv_S2S1 = W[np.ix_(iS_S2, iV_S1)]*gI
        Wvv_S2S1 = W[np.ix_(iV_S2, iV_S1)]*gI

        # put them back together
        # Reconstructing W0 from the submatrices, negating where indicated
        # TODO: either leave the values in this order or sort them in the parameter file!
        #rowA3b1 = np.hstack([
        #rowA3b2 = np.hstack([
        #rowA3b3 = np.hstack([
        row1 = np.hstack([Wee_S1, -Wep_S1, -Wes_S1, -Wev_S1, Wee_S1S2, -Wep_S1S2, -Wes_S1S2, -Wev_S1S2])
        row2 = np.hstack([Wpe_S1, -Wpp_S1, -Wps_S1, -Wpv_S1, Wpe_S1S2, -Wpp_S1S2, -Wps_S1S2, -Wpv_S1S2])
        row3 = np.hstack([Wse_S1, -Wsp_S1, -Wss_S1, -Wsv_S1, Wse_S1S2, -Wsp_S1S2, -Wss_S1S2, -Wsv_S1S2])
        row4 = np.hstack([Wve_S1, -Wvp_S1, -Wvs_S1, -Wvv_S1, Wve_S1S2, -Wvp_S1S2, -Wvs_S1S2, -Wvv_S1S2])
        row5 = np.hstack([Wee_S2S1, -Wep_S2S1, -Wes_S2S1, -Wev_S2S1, Wee_S2, -Wep_S2, -Wes_S2, -Wev_S2])
        row6 = np.hstack([Wpe_S2S1, -Wpp_S2S1, -Wps_S2S1, -Wpv_S2S1, Wpe_S2, -Wpp_S2, -Wps_S2, -Wpv_S2])
        row7 = np.hstack([Wse_S2S1, -Wsp_S2S1, -Wss_S2S1, -Wsv_S2S1, Wse_S2, -Wsp_S2, -Wss_S2, -Wsv_S2])
        row8 = np.hstack([Wve_S2S1, -Wvp_S2S1, -Wvs_S2S1, -Wvv_S2S1, Wve_S2, -Wvp_S2, -Wvs_S2, -Wvv_S2])
        
        # Vertically stack the rows to form W0
        W0_S1S2 = np.vstack([row1, row2, row3, row4, row5, row6, row7, row8])
        # add area 3b
        to_area3b = W[0:3, 3:]
        #from_area3b = W[:, :3# TODO: this has to be sorted so that S1 receives input from 3b in the correct order!!!
        W0 = np.vstack((to_area3b, W0_S1S2))
        W0 = np.hstack((from_area3b, W0))
        #print(Wee_S1.shape)
        #print(W.shape, W0.shape)

        # include the external input to the matrix 
        if include_Iext:
            #print('W0', W0.shape) # this should be 13x13
            
            # append the thalamus population(s) values to the matrix
            W0 = np.append(W0, W_to_thal, axis=0)
            #print('W from thalamus 2\n', W_from_thal.shape)
            #print(W_from_thal)

            W0 = np.append(W0, W_from_thal.T, axis=1)

            W = np.concatenate((W0, Wb, Wext), axis=1)
            #print('\nW matrix shape', W.shape) 
            #print('W matrix', W[:,-4:]) 

        return W


    def get_sigmoid(self):
        # sigmoid function (nPop x 3) --> 3 stands for parameters: r(1/mV), v_thr(mV), m_max (1/s)
        
        sigmoid_params_A3b = np.array([[  0.12782346,  32.10540543,  31.39696397],
                                    [  0.14218422,  40.03107351, 166.82960408],
                                    [  0.07937015,  42.01276379,  56.95305832]])
        
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
