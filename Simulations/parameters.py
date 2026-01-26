# %%
import numpy as np
import pandas as pd
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
        nPopA3b = 4 # E, PV, SST, VIP
        nPopS1 = 13 # 4*E, 4*PV, 4*SST, 1*VIP
        nPopS2 = 13 # 4*E, 4*PV, 4*SST, 1*VIP
        nPopThal = 2 # one forward excitatory and one inhibitory feedback neuron
        nPopTotal = nPopS1+nPopS2+nPopA3b+nPopThal

        # SYNAPTIC DECAY (depends on the connection type excitatory/inhibitory)
        # Based on visual cortex values
        # the last two values are used for the external input and background input
        # with nPopTotal = 32 the shape of tau should be (32, 32) 
        # order: E1, PV1, SST1, VIP1, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4 
        tauA3b = np.tile(np.array([6,3,20,15])*1e-3, (nPopTotal,1)) # sec
        tauS1 = np.tile(np.array([6,3,20,15,6,3,20,6,3,20,6,3,20])*1e-3, (nPopTotal,1)) # sec
        tauS2 = np.tile(np.array([6,3,20,15,6,3,20,6,3,20,6,3,20,3,3])*1e-3, (nPopTotal,1)) # sec
        tau = np.hstack((tauA3b,tauS1,tauS2))

        # Alternative, old values for tau
        # E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1
        #tau = np.tile(np.array([2,2,2,2,4,4,4,4,4,4,4,4,4,3])*1e-3, (nPop+1,1)) 
        #tau = np.tile(np.array([5.2, 5.2, 5.9, 5.9, 3, 3, 3.8, 3.8, 11.2, 11.2, 11.1, 11.1, 10.4])*1e-3, (nPop+1,1)) 

        # TODO: MEMBRANE CONSTANT (NOTE: this value depends on post-synaptic neuron whereas the synaptic decay depends on the pre-synapse)
        
        return tau, nPopTotal

    def get_connectProb(self):

        # Connection probabilies
        # Targets in rows, sources in columns 
        # Values for area 3b are merged values from S1 (after already computed W = P*S*C)
        
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
        
        
        # forward connections: S1 Layer 5, E3 --> S2 layer 4 E2 
        P_S1toS2 = np.zeros((13,13)) # TODO: implement S1 to S2 connection! 
        
        # feedback connections: S2 Layer 5, E3 --> S2 layer 4 E2 
        P_S2toS1 = np.zeros((13,13)) # TODO: implement feedback from S2 to S1

        P_toS1 = np.hstack((P_S1, P_S2toS1))
        P_toS2 = np.hstack((P_S1toS2, P_S2))
        P = np.vstack((P_toS1, P_toS2)) #  26x26 (with nPop=13 per S1/S2)
    
        return P


    def get_connectStrength(self):
        """ 
        Postsynaptic potential (13x13) averages from Isbister, Jiang, and more (see excel file FinalConnectivity_PSP.csv for mor info).
        Again, values for area 3b are merged from S1 later on. 
        """
                    
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
        
        # to S1
        psp_S2toS1 = np.zeros((13,13))
        psp_toS1 = np.hstack((psp_S1, psp_S2toS1))
        
        # to S2
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
        C = np.hstack((C_S1, C_S2))

        return C 

    def get_raw_connectivity(self, thal_connect, extI_cellcount, bI_cellcount, thal_cellcount):
        """Put together the connevtivity matrix (not yet scaled by coupling strength and EI-balance parameter)
        
        Args:
            thal_connect: connecivity between thalamic neurons (E and I)
            extI_cellcount: number of external input neurons (in the thalamus)
            bI_cellcount: number of background input neurons (from other cortical areas)
            thal_cellcount: number of neurons in the thalamus connecting to the somatosensory area

        Returns: 
            W0 (2D numpy array): connectivity matrix of only cortical populations (no thalamus)! 
            W_to_thal (2D numpy array): 
            W_from_thal (2D numpy array):
            Wb (1D numpy array): background input array
            Wext (1D numpy array): external input array
        """

        S = self.get_connectStrength()
        P = self.get_connectProb()
        C = self.get_cellcounts()

        # Final connectivity matrix s
        # PS = P * S (16 x 16) --> Connectivity probability if all cells were equally distributed
        # W0 = PS * C (16 x 16) --> make sure that cell counts are accounted for! C is shaped (16 x 1) so tile it before multiplying  
        PS = P * S
        nPop = len(PS)
        W0 = PS * np.tile(C, (nPop,1))

        idx_S1_E = [0,4,7,10]
        idx_S1_P = [1,5,8,11]
        idx_S1_S = [2,6,9,12]
        idx_S1_V = [3]

        # set up the connectivity for area 3b as a weighted average of the connectivity from S1
        # ... within A3b
        # ... merge input from all E to all other populations
        WS1_collapse_sources_E = np.dot(W0[:int(len(W0)/2),idx_S1_E], C[idx_S1_E])/np.sum(C[idx_S1_E])  
        W_A3bA3b_EE = np.dot(WS1_collapse_sources_E[idx_S1_E], C[idx_S1_E])/np.sum(C[idx_S1_E]) # .. now we have to compute the weighted average for the all receiving cells
        W_A3bA3b_PE = np.dot(WS1_collapse_sources_E[idx_S1_P], C[idx_S1_P])/np.sum(C[idx_S1_P])
        W_A3bA3b_SE = np.dot(WS1_collapse_sources_E[idx_S1_S], C[idx_S1_S])/np.sum(C[idx_S1_S])
        W_A3bA3b_VE = np.dot(WS1_collapse_sources_E[idx_S1_V], C[idx_S1_V])/np.sum(C[idx_S1_V])
        # ... merge input from all P targeting all other populations
        WS1_collapse_sources_P = np.dot(W0[:int(len(W0)/2),idx_S1_P], C[idx_S1_P])/np.sum(C[idx_S1_P])  
        W_A3bA3b_EP = np.dot(WS1_collapse_sources_P[idx_S1_E], C[idx_S1_E])/np.sum(C[idx_S1_E])
        W_A3bA3b_PP = np.dot(WS1_collapse_sources_P[idx_S1_P], C[idx_S1_P])/np.sum(C[idx_S1_P])
        W_A3bA3b_SP = np.dot(WS1_collapse_sources_P[idx_S1_S], C[idx_S1_S])/np.sum(C[idx_S1_S])
        W_A3bA3b_VP = np.dot(WS1_collapse_sources_P[idx_S1_V], C[idx_S1_V])/np.sum(C[idx_S1_V])
        # ... merge input from all S targeting all other populations
        WS1_collapse_sources_S = np.dot(W0[:int(len(W0)/2),idx_S1_S], C[idx_S1_S])/np.sum(C[idx_S1_S])  
        W_A3bA3b_ES = np.dot(WS1_collapse_sources_S[idx_S1_E], C[idx_S1_E])/np.sum(C[idx_S1_E])
        W_A3bA3b_PS = np.dot(WS1_collapse_sources_S[idx_S1_P], C[idx_S1_P])/np.sum(C[idx_S1_P])
        W_A3bA3b_SS = np.dot(WS1_collapse_sources_S[idx_S1_S], C[idx_S1_S])/np.sum(C[idx_S1_S])
        W_A3bA3b_VS = np.dot(WS1_collapse_sources_S[idx_S1_V], C[idx_S1_V])/np.sum(C[idx_S1_V])
        # ... merge input from all V targeting all other populations
        WS1_collapse_sources_V = np.dot(W0[:int(len(W0)/2),idx_S1_V], C[idx_S1_V])/np.sum(C[idx_S1_V])  
        W_A3bA3b_EV = np.dot(WS1_collapse_sources_V[idx_S1_E], C[idx_S1_E])/np.sum(C[idx_S1_E])
        W_A3bA3b_PV = np.dot(WS1_collapse_sources_V[idx_S1_P], C[idx_S1_P])/np.sum(C[idx_S1_P])
        W_A3bA3b_SV = np.dot(WS1_collapse_sources_V[idx_S1_S], C[idx_S1_S])/np.sum(C[idx_S1_S])
        W_A3bA3b_VV = np.dot(WS1_collapse_sources_V[idx_S1_V], C[idx_S1_V])/np.sum(C[idx_S1_V])

        # stack them all together to create the A3b to A3b connection matrix
        W_A3bA3b = np.hstack((np.vstack((W_A3bA3b_EE, W_A3bA3b_PE, W_A3bA3b_SE, W_A3bA3b_VE)),
                            np.vstack((W_A3bA3b_EP, W_A3bA3b_PP, W_A3bA3b_SP, W_A3bA3b_VP)),
                            np.vstack((W_A3bA3b_ES, W_A3bA3b_PS, W_A3bA3b_SS, W_A3bA3b_VS)),
                            np.vstack((W_A3bA3b_EV, W_A3bA3b_PV, W_A3bA3b_SV, W_A3bA3b_VV))))
        
        #W_A3bA3b = np.zeros((4, 4))
        
        # ... from  area 3b to S1
        # this means we have to fuse the source populations together
        W_S1A3b = np.column_stack((WS1_collapse_sources_E, WS1_collapse_sources_P, WS1_collapse_sources_S, WS1_collapse_sources_V))
        W_S1A3b = np.zeros((13, 4))
 
        # ... from S1 to area 3b
        # we need to compress the target populations 
        WS1_collapse_targets_E = np.dot(W0[idx_S1_E,:int(len(W0)/2)].T, C[idx_S1_E])/np.sum(C[idx_S1_E])
        WS1_collapse_targets_P = np.dot(W0[idx_S1_P,:int(len(W0)/2)].T, C[idx_S1_P])/np.sum(C[idx_S1_P])
        WS1_collapse_targets_S = np.dot(W0[idx_S1_S,:int(len(W0)/2)].T, C[idx_S1_S])/np.sum(C[idx_S1_S])
        WS1_collapse_targets_V = np.dot(W0[idx_S1_V,:int(len(W0)/2)].T, C[idx_S1_V])/np.sum(C[idx_S1_V])
        W_A3bS1 = np.vstack((WS1_collapse_targets_E, WS1_collapse_targets_P, WS1_collapse_targets_S, WS1_collapse_targets_V))

        # TODO: connectivity between S2 and A3b has to be defined manually
        W_S2A3b = np.zeros(W_S1A3b.shape) 
        #W_S2A3b[4,0] = 100
        #W_S2A3b[5,0] = 100
        W_A3bS2 = np.zeros(W_A3bS1.shape) 
        #W_A3bS2[0,4] = 50
        #W_A3bS2[1,4] = 50

        # TODO: connectivity between S2 and S1 also has to be defined manually
        # feedforward
        #W0[idx_S1_E[1]+13,idx_S1_E[2]] = 100 # S1 layer 5 E to S2 layer 4 E
        #W0[idx_S1_P[1]+13,idx_S1_E[2]] = 100
        # feedback
        #W0[idx_S1_E[0],idx_S1_E[2]+13] = 50 # S2 layer 5 E to S1 layer 1 E 
        #W0[idx_S1_P[0],idx_S1_E[2]+13] = 50 # S2 layer 5 E to S1 layer 1 PV

        # stack the matrices together 
        W0 = np.vstack((np.hstack((W_A3bS1, W_A3bS2)), W0))
        W0 = np.hstack((np.vstack((W_A3bA3b, W_S1A3b, W_S2A3b)), W0))

        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        S_thalToS1 = np.array([[0.49, 0.49, 0.245, 0, 0.49, 0.49, 0.245, 0.49, 0.49, 0.245, 0.49, 0.49, 0.245], # Thal to S1: Based on Jiang et al. 2023 only! 
                              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) # no thal E/I target (integrated in last two entries of below S2 array) 
        
        S_A3bto_thal = np.array([[0, 0, 0, 0], [0, 0, 0, 0]])
        S_S1to_thal = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
        S_S2to_thal = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
        S_to_thal = np.hstack((S_A3bto_thal, S_S1to_thal,S_S2to_thal))

        # connection probabilities
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        P_thalToS1 = np.array([[6.2, 6.2, 0, 0, 40, 40, 20, 25.9, 25.9, 0, 9, 9, 0], # thalamus excitatory
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])*1e-2 # reticlar nucleus inhibitory (estimation)

        # order: E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        P_A3bto_thal = np.array([[0, 0, 0, 0], [0, 0, 0, 0]])
        P_S1to_thal = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) # this is only from the cortex!
        P_S2to_thal = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) # this is only from the cortex!
        P_to_thal = np.hstack((P_A3bto_thal, P_S1to_thal,P_S2to_thal))
        
        # calculate final thalamus connectivity
        PS_to_thal = np.multiply(P_to_thal, S_to_thal)
        a3b_C = np.array([np.sum(C[idx_S1_E]), np.sum(C[idx_S1_P]), np.sum(C[idx_S1_S]), np.sum(C[idx_S1_V])])
        C_all = np.concatenate((a3b_C, C))
        W_to_thal = np.multiply(PS_to_thal, C_all)
        # thalamus input to S1 connectivity  
        PS_from_thal_S1 = np.multiply(P_thalToS1, S_thalToS1)
        cell_count_thal = thal_cellcount # it is 230 in Jiang et al. 2023 but for us it might be different!
        W_from_thal_S1 = np.multiply(PS_from_thal_S1, cell_count_thal)

        # now that we have the connectivity from- and to thalamus and S1/S2,
        # we can collapse it to get input array to area 3b 
        W_A3bThal_E = np.sum(W_from_thal_S1[:,idx_S1_E], axis=1)
        W_A3bThal_P = np.sum(W_from_thal_S1[:,idx_S1_P], axis=1)
        W_A3bThal_S = np.sum(W_from_thal_S1[:,idx_S1_S], axis=1)
        W_A3bThal_V = np.sum(W_from_thal_S1[:,idx_S1_V], axis=1)
        W_A3bThal = np.vstack((W_A3bThal_E, W_A3bThal_P, W_A3bThal_S, W_A3bThal_V))

        # adjust input from thalamus to S1 (should be slightly weaker compared to A3b)
        W_from_thal_S1 = W_from_thal_S1 * 0.8
        # thalamus input to S2 is just a fifth of what A3b receives 
        W_from_thal_S2 = W_from_thal_S1 * 0.2
        # join area 3b, S1 and S2 input arrays
        W_from_thal = np.hstack((W_A3bThal.T, W_from_thal_S1, W_from_thal_S2))
        # add within thalamus (E and I) connectivity 
        tEE, tEI, tIE, tII = thal_connect
        W_from_thal = np.hstack((W_from_thal, [[tEE,tEI],[tIE,tII]]))

        # Create the external input matrix
        # Only the thalamus E population receives the external input (from the brain stem)
        Wext = np.zeros((W_from_thal.shape[1],1))
        Wext[-2] = 1 * extI_cellcount # thalamus E population
        Wext[-1] = 0 # reticular I population

        # .. and also for the background input (all cells receive input except the thalamus)
        Wb = np.zeros((W_from_thal.shape[1],1))
        Wb[:-2] = 1 * bI_cellcount # cellcount from background input

        return W0, W_to_thal, W_from_thal, Wb, Wext


    def get_connectivity(self, gE, gI, gEthal, gEthal, thal_connect, extI_cellcount, bI_cellcount, thal_cellcount, area='all'):
        """Apply coupling strength parameter and compute the final connectivity matrix.  
        
        Args:
            gE, gI (float): Coupling strength parameter for E and I populations.   
                            It should be the result of the equation: g * bEI (where g is the general coupling strength and bEI is the EI-balance).
            gEthal, gEthal (float): Same as above, just for thalamic neurons.
            thal_connect (int): connecivity between thalamic neurons (E and I)
            extI_cellcount (int): number of external input neurons (in the thalamus)
            bI_cellcount (int): number of background input neurons (from other cortical areas)
            thal_cellcount (int): number of neurons in the thalamus connecting to the somatosensory area
            area (string): area to chose if you want to look at an isolated area (sets all other connections to zero).
                           Options: 'ThalA3b', 'S1', 'A3b', 'S1','ThalS1', 'A1', 'ThalA1', 'ThalA1S2', 'S2', 'A1S2'
        
        Returns:
            2D numpy array: complete connectivity matrix (32x34)
        """
        

        W0, W_to_thal, W_from_thal, Wb, Wext = self.get_raw_connectivity(thal_connect, extI_cellcount, bI_cellcount, thal_cellcount)

        # make inhibitory connections negative and apply weights gI and gE respectively
        idx_I_A3b = np.array([1,2,3])
        idx_I_S = [1,2,3,5,6,8,9,11,12]
        idx_I_S1 = np.array(idx_I_S)+4 
        idx_I_S2 = np.array(idx_I_S)+17
        idx_I = np.concatenate((idx_I_A3b,idx_I_S1,idx_I_S2))
        idx_E_A3b = np.array([0])
        idx_E_S = [0,4,7,10]
        idx_E_S1 = np.array(idx_E_S)+4
        idx_E_S2 = np.array(idx_E_S)+17
        idx_E = np.concatenate((idx_E_A3b,idx_E_S1,idx_E_S2))

        # scale the coupling strength
        W0[:,idx_I] = W0[:,idx_I] * -gI # negative weight for inhibitory connections
        W0[:,idx_E] = W0[:,idx_E] * gE
        W_to_thal[:,idx_I] = W_to_thal[:,idx_I] * -gI
        W_to_thal[:,idx_E] = W_to_thal[:,idx_E] * gE
        W_from_thal[0,:] = W_from_thal[0,:] * gEthal
        # weight the inhibitory popultaion (from the reticular nucleus in the thalamus) as negative
        W_from_thal[1,:] = W_from_thal[1,:] * -gIthal
        
        # include the external input to the matrix 
        # append the thalamus population(s) values to the matrix
        W0 = np.append(W0, W_to_thal, axis=0)
        W0 = np.append(W0, W_from_thal.T, axis=1)
        W = np.concatenate((W0, Wb, Wext), axis=1)

        
        # if we don't want to use the entire model, just use one area 
        # based on what area was chosen, set the respective connectivity values to 0 
        if area=='ThalA3b':
            W[4:-2,:] = 0 # to A1, S2 from everywhere
            W[:,4:-4] = 0 # to everywhere from A1, S2
        elif area=='A3b':
            W[4:] = 0 # to A1, S2, Thal from everywhere
            W[:,4:-2] = 0 
        elif area=='S1':
            W[4+13:,4+13:] = 0
        elif area=='ThalS1':
            # this still includes A3b
            W[4+13:-2,:] = 0  # cut out S2 
            W[:,4+13:-4] = 0  # cut out S2 
        elif area=='A1':
            W[:4,:] = 0 # cut out Area3b
            W[:,:4] = 0 # cut out Area3b
            W[4+13:,:] = 0 
            W[:,4+13:-2] = 0 
            W[-2:,:] = 0 # cut out Thalamus
            W[:,-4:-2] = 0 # cut out input from Thalamus
        elif area=='ThalA1':
            W[:4,:] = 0 # cut out Area3b
            W[:,:4] = 0 # cut out Area3b
            W[4+13:-2,:] = 0  # cut out S2 
            W[:,4+13:-4] = 0  # cut out S2
        elif area=='ThalA1S2':
            W[:4,:] = 0 # cut out Area3b
            W[:,:4] = 0 # cut out Area3b
        elif area=='A1S2':
            W[:4,:] = 0 # cut out Area3b
            W[:,:4] = 0 # cut out Area3b  
            W[-2:,:] = 0 # cut out Thalamus
            W[:,-4:-2] = 0 # cut out input from Thalamus
        elif area=='S2':
            W[:4,:] = 0 # cut out Area3b
            W[:,:4] = 0 # cut out Area3b 
            W[4:4+13,:] = 0 
            W[:,4:4+13] = 0 
            W[-2:,:] = 0 # cut out Thalamus
            W[:,-4:-2] = 0 # cut out input from Thalamus

        return W

    def get_sigmoid(self):
        """
        Sigmoid function (nPop x 3) --> 3 stands for parameters: r(1/mV), v_thr(mV), m_max (1/s)
        
        Returns:
            numpy array: sigmoid values.
        """
        
        # sigmoid function (nPop x 3) --> 3 stands for parameters: r(1/mV), v_thr(mV), m_max (1/s)
        
        
        sigmoid_params_A3b = [[  0.12782346,  32.10540543,  31.39696397],
                                    [  0.14218422,  40.03107351, 166.82960408],
                                    [  0.07937015,  42.01276379,  56.95305832],
                                    [  0.0704119 ,  37.86409387,  38.52689646]]
        
        # order:
        # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
        sigmoid_params_S1 = [[  0.12782346,  32.10540543,  31.39696397],
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
                                    [  0.07937015,  42.01276379,  56.95305832]]

        sigmoid_params_S2 = [[  0.12782346,  32.10540543,  31.39696397],
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
                                    [  0.1,  40,  30]] # Thalamus I

        sigmoid_params = np.vstack((sigmoid_params_A3b, sigmoid_params_S1, sigmoid_params_S2))
                                       
        return sigmoid_params
