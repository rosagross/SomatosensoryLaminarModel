# %%
import numpy as np
import yaml
# %%
class Parameter():

    def __init__(self, cortex_type='somato'):
        self.cortex_type = cortex_type
        self.tau, self.nPop = self.get_params()
        self.S = self.get_connectStrength()
        self.P = self.get_connectProb()
        self.C = self.get_cellcounts()
        self.sigmoid_params = self.get_sigmoid()

    def get_params(self):

        if self.cortex_type == 'somato':
            # nr. of populations
            nPopS1 = 13 #  E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1
            nPopS2 = 13 #  E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1
            nPopThal = 2 # one forward excitatory and one inhibitory feedback neuron
            nPopTotal = nPopS1+nPopS2+nPopThal
            # SYNAPTIC DECAY (depends on the connection type excitatory/inhibitory)
            #tau = np.tile(np.array([2,2,2,2,4,4,4,4,4,4,4,4,4,3])*1e-3, (nPop+1,1)) 
            #tau = np.tile(np.array([5.2, 5.2, 5.9, 5.9, 3, 3, 3.8, 3.8, 11.2, 11.2, 11.1, 11.1, 10.4])*1e-3, (nPop+1,1)) 
            # Visual cortex values
            # the last two values are used for the external input and background input
            # with nPopTotal = 28 the shape of tau should be (28, 32) 
            tauS1 = np.tile(np.array([6,6,6,6,3,3,3,3,20,20,20,20,15])*1e-3, (nPopTotal,1)) # sec
            tauS2 = np.tile(np.array([6,6,6,6,3,3,3,3,20,20,20,20,15,3,3,3,3])*1e-3, (nPopTotal,1)) # sec
            tau = np.hstack((tauS1,tauS2))

            # TODO: MEMBRANE CONSTANT (NOTE: this value depends on post synaptic neuron whereas the synaptic decay depends on the presynapse)
        
        
        elif self.cortex_type == 'visual':
            # nr. of populations
            nPopTotal = 16
            # E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1, VIP2, VIP3, VIP4, Iext
            tau = np.tile(np.array([6,6,6,6,3,3,3,3,20,20,20,20,15,15,15,15,3])*1e-3, (nPop+1,1)) # sec
            # synaptic kernel efficacy
        
        return tau, nPopTotal

    def get_connectProb(self):

        # Connection probabilies
        if self.cortex_type == 'somato':
            # Connection probabilies
            # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4 (Thalamus is added later!)
            # Targets in rows, sources in columns 
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

            P_S1toS2 = np.zeros((13,13))
            P_S2toS1 = np.zeros((13,13))
            P_toS1 = np.hstack((P_S1,P_S2toS1))
            P_toS2 = np.hstack((P_S1toS2,P_S2))
            P = np.vstack((P_toS1, P_toS2)) #  26x26 (with nPop=13 per S1/S2)

        elif self.cortex_type == 'visual':
            # E1, PV1, SST1, VIP1, E2, PV2, SST2, VIP2, E3, PV3, SST3, VIP3, E4, PV4, SST4, VIP4
            P = np.array([[0.110000000000000,	0.494409448818898,	0.319464566929134,	0.0836692913385827,	0.170000000000000,	0.532440944881890,	0.463984251968504,	0,	0, 0.380314960629921, 0.220582677165354, 0, 0,	0,	0,	0],
                        [0.621051968503937, 0.788569291338583, 0.439203149606299, 0.0399275590551181, 0.416655118110236, 0.529040157480315, 0, 0, 0.0776893006021306, 0.449185039370079, 0, 0, 0, 0.459166929133858, 0, 0],
                        [0.597467716535433, 0.279492913385827, 0.169692125984252, 0.489112598425197, 0.0707527559055118, 0.159710236220472, 0.179674015748031, 0.998188976377953, 0, 0, 0.189655905511811, 0, 0, 0, 0, 0],
                        [0.298733858267717, 0, 0.918333858267717, 0.0499094488188976, 0.0393070866141732, 0.269511023622047, 0.558985826771654, 0, 0, 0, 0.339384251968504, 0, 0, 0, 0.798551181102362, 0],
                        [0.120000000000000, 0.304251968503937, 0.0456377952755906, 0, 0.220000000000000, 0.395527559055118, 0.410740157480315, 0, 0, 0.684566929133858, 0.266220472440945, 0, 0, 0, 0, 0],
                        [0.479546456692913, 0.419239370078740, 0.658804724409449, 0, 0.283011023622047, 0.968243307086614, 0.668786614173228, 0.0399275590551181, 0.0493265400648448, 0.658804724409449, 0.628859055118110, 0, 0, 0.998188976377953, 0, 0],
                        [0.314456692913386, 0, 0, 0.439203149606299, 0.0786141732283465, 0.139746456692913, 0.0399275590551181, 0.499094488188976, 0.106641834988018, 0, 0.109800787401575, 0.479130708661417, 0, 0, 0, 0],
                        [0, 0, 0.429221259842520, 0.249547244094488, 0, 0.279492913385827, 0.818514960629921, 0.129764566929134, 0.0362392164394085, 0.189655905511811, 0.658804724409449, 0, 0, 0, 0, 0],
                        [0.222222222222222, 0.249106299212598, 0, 0, 0, 0.230292846372927, 0.132764495347173, 0, 0.0950048875855327, 0.302338727597701, 0.169197900262467, 0, 0, 0.273826771653543, 0.0624803149606299, 0],
                        [0, 0.558985826771654, 0.249547244094488, 0.109800787401575, 0.542437795275591, 0.409257480314961, 0, 0.179674015748031, 0.162590469254927, 0.499094488188976, 0.289474803149606, 0.0199637795275591, 0, 0.159710236220472, 0.219601574803150, 0],
                        [0.0786141732283465, 0, 0, 0, 0.0628913385826772, 0, 0.119782677165354, 0.389293700787402, 0.199036793128132, 0.139746456692913, 0.0998188976377953, 0.189655905511811, 0, 0.0598913385826772, 0.0399275590551181, 0],
                        [0, 0.568967716535433, 0, 0.269511023622047, 0.235842519685039, 0, 0.259529133858268, 0.229583464566929, 0.0434834645669291, 0.0998188976377953, 0.329402362204724, 0.149728346456693, 0, 0.249547244094488, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.144519685039370, 0.0100000000000000, 0.456377952755906, 0.114094488188976, 0],
                        [0, 0, 0, 0, 0, 0.998188976377953, 0, 0, 0.212258267716535, 0.269511023622047, 0.119782677165354, 0, 0.338040944881890, 0.748641732283465, 0.339384251968504, 0.139746456692913],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0.0650600054303557, 0.109800787401575, 0.0698732283464567, 0.259529133858268, 0.243703937007874, 0.0598913385826772, 0.159710236220472, 0.189655905511811],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.129764566929134, 0, 0.0598913385826772, 0.708714173228346, 0]])

        return P


    def get_connectStrength(self):

        # Synaptic strength
        if self.cortex_type == 'somato':

            # Postsynaptic potential (13x14) averages from Isbister, Jiang, and more (see excel file FinalConnectivity_PSP.csv for mor info)
            # order: E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
            # Targets in rows, sources in columns
            psp_S1 = [
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

            psp_S1toS2 = np.zeros((13,13))
            psp_S2toS1 = np.zeros((13,13))
            psp_toS1 = np.hstack((psp_S1,psp_S2toS1))
            psp_toS2 = np.hstack((psp_S1toS2,psp_S2))
            psp = np.vstack((psp_toS1, psp_toS2)) #  26x26 (with nPop=13 per S1/S2)

        elif self.cortex_type == 'visual':
            # E1, PV1, SST1, VIP1, E2, PV2, SST2, VIP2, E3, PV3, SST3, VIP3, E4, PV4, SST4, VIP4
            psp = np.array([[0.36,	1.49,	0.86,	1.31,	0.34,	1.39,	0.69,	0.91,	0.74,	1.32,	0.53,	0,	0,	0,	0,	0],
                            [0.48,	0.68,	0.42,	0.41,	0.56,	0.68,	0.42,	0.41,	0.2,	0.79,	0,	0,	0,	0,	0,	0],
                            [0.31,	0.5,	0.15,	0.52,	0.3,	0.5,	0.15,	0.52,	0.22,	0,	0,	0,	0,	0,	0,	0],
                            [0.28,	0.18,	0.32,	0.37,	0.29,	0.18,	0.32,	0.37,	0,	0,	0,	0,	0,	0,	0,	0],
                            [0.78,	1.39,	0.69,	0.91,	0.83,	1.29,	0.51,	0.51,	0.63,	1.25,	0.52,	0.91,	0.96,	0,	0,	0],
                            [0.56,	0.68,	0.42,	0.41,	0.64,	0.68,	0.42,	0.41,	0.73,	0.94,	0.42,	0.41,	0,	0,	0,	0],
                            [0.3,	0.5,	0.15,	0.52,	0.29,	0.5,	0.15,	0.52,	0.28,	0.45,	0.28,	0.52,	0,	0,	0,	0],
                            [0.29,	0.18,	0.32,	0.37,	0.29,	0.18,	0.32,	0.37,	0,	0.18,	0.33,	0.37,	0,	0,	0,	0],
                            [0.47,	1.25,	0.52,	0.91,	0.38,	1.25,	0.52,	0.91,	0.75,	1.2,	0.52,	1.31,	0.4,	2.5,	0.52,	1.31],
                            [0,	0.51,	0,	0,	0,	0.94,	0.42,	0.41,	0.81,	1.19,	0.41,	0.41,	0.81,	1.19,	0.41,	0.41],
                            [0.25,	0,	0.39,	0,	0.28,	0.45,	0.28,	0.52,	0.27,	0.4,	0.4,	0.52,	0.27,	0.4,	0.4,	0.52],
                            [0,	0,	0,	0,	0.29,	0.18,	0.33,	0.37,	0.28,	0.18,	0.33,	0.37,	0.28,	0.18,	0.33,	0.37],
                            [0,	0,	0,	0,	0,	0,	0,	0,	0.23,	2.5,	0.52,	1.31,	0.94,	3.8,	0.52,	1.31],
                            [0.81,	0,	0,	0,	0.81,	0,	0,	0,	0.81,	1.19,	0.41,	0.41,	0.81,	1.19,	0.41,	0.41],
                            [0,	0,	0,	0,	0,	0,	0,	0,	0.27,	0.4,	0.4,	0.52,	0.27,	0.4,	0.4,	0.52],
                            [0,	0,	0,	0,	0,	0,	0,	0,	0.28,	0.18,	0.33,	0.37,	0.28,	0.18,	0.33,	0.37]]).T
            
        return psp

    def get_cellcounts(self):
        
        # Cell counts
        if self.cortex_type == 'somato':
            # E1, PV1, SST1, VIP, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4
            C_absoluteS1 = np.array([1691, 90, 74, 85, 1656, 85, 48, 1095, 109, 105, 1288, 56, 66])
            C_absoluteS2 = np.array([1691, 90, 74, 85, 1656, 85, 48, 1095, 109, 105, 1288, 56, 66])
            C_absolute = np.hstack((C_absoluteS1,C_absoluteS2))

            # translate in proportion
            C = C_absolute/np.sum(C_absolute)
            
        elif self.cortex_type == 'visual':
            # E1, PV1, SST1, VIP1, E2, PV2, SST2, VIP2, E3, PV3, SST3, VIP3, E4, PV4, SST4, VIP4
            C = np.array([0.317584074016321, 0.0173259178404976, 0.0131590345998579, 0.0138260824367171, 0.175424918135038, 0.00974582878527990, 
                        0.00594278981929068, 0.00315331704697056, 0.166432766775820, 0.0216573973006220, 0.0127345496127657, 
                        0.00266819134743663, 0.195020531212641, 0.0119115685153421, 0.00848969974184382, 0.00291075419720360])

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

        # sort the values like this: 
        if self.cortex_type == 'visual':
            iE = np.arange(0, 16, 4)  # E1, E2, E3, E4
            iP = np.arange(1, 16, 4)  # PV1, PV2, PV3, PV4
            iS = np.arange(2, 16, 4)  # SOM1, SOM2, SOM3, SOM4
            iV = np.arange(3, 16, 4)  # VIP1, VIP2, VIP3, VIP4

            # now create the external input matrix 
            Wext = np.zeros((nPop,1))
            Wext[1] = 1

        elif self.cortex_type == 'somato':
        
            # create the external input matrix, based on thalamus connectivity (average of findings, see FinalConnectivity_PSPs.ods)
            # order: 'E1','E2','E3','E4','PV1','PV2','PV3','PV4','SST1','SST2','SST3','SST4','VIP1', 'Thal E'
            #S_thal = np.array([0.49, 1.45, 0.5, 0.85, 0.49, 2.3, 0.49, 2.2, 0.245, 0.245, 0.245, 0.245, 0]) # Thal to S1 based on Isbister & jiang 
            S_thalToS1 = np.array([[0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.245, 0.245, 0.245, 0.245, 0], # Thal to S1: Based on Jiang et al. 2023 only! 
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) # no thal E/I target (integrated in last two entries of below S2 array) 
            S_thalToS2 = np.array([[0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.245, 0.245, 0.245, 0.245, 0, 0, 0.5], # last to are Thal E and Thal I 
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0]]) # Reticular inhibition: just an assumption
            #print('S', S_thalToS1.shape, S_thalToS2.shape)
            
            S_from_thal = np.hstack((S_thalToS1, S_thalToS2))
            S_S1to_thal = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
            S_S2to_thal = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
            S_to_thal = np.hstack((S_S1to_thal,S_S2to_thal))


            #print('S', S_from_thal.shape)
            
            # cell count of the thalamus 
            C_thal = [1, 1] # Thal E, Thal I
            
            # connection probabilities
            #P_from_thal = np.array([[0,0,0,0,0,0,0,0,0,0, 0, 0, 0, 0, 0], # thalamus excitatory
            P_thalToS1 = np.array([[6.2, 40, 25.9, 9, 6.2, 40, 25.9, 9, 0, 20, 0, 0, 0], # thalamus excitatory
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])*1e-2 # reticlar nucleus inhibitory (estimation)
            P_thalToS2 = np.array([[6.2, 40, 25.9, 9, 6.2, 40, 25.9, 9, 0, 20, 0, 0, 0, 0, 0], # thalamus excitatory
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])*1e-2 # reticlar nucleus inhibitory (estimation)
            #P_thalToS2 = np.array([[6.2, 40, 25.9, 9, 6.2, 40, 25.9, 9, 0, 20, 0, 0, 0, 0, 0], # thalamus excitatory
            #                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])*1e-2 # reticlar nucleus inhibitory (estimation)
            # no projections from thalamus to S2 
            #P_thalToS2 = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], # thalamus excitatory
            #                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])*1e-2 # reticlar nucleus inhibitory (estimation)


            #print('P from thal', P_thalToS1.shape, P_thalToS2.shape)
            P_from_thal = np.hstack((P_thalToS1, P_thalToS2))
            P_S1to_thal = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) # this is only from the cortex!
            P_S2to_thal = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) # this is only from the cortex!
            P_to_thal = np.hstack((P_S1to_thal,P_S2to_thal))
            
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
            Wext[26] = 1 # thalamus E population
            Wext[27] = 0 # reticular I population

            # .. and also for the background input (all cells receive input except the thalamus)
            Wb = np.zeros((W_from_thal.shape[1],1))
            Wb[:-2] = 1

            # S1
            # indices to reorder the matrix to E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1
            iE_S1 = np.array([0, 4, 7, 10])  # E1, E2, E3, E4 of S1 
            iP_S1 = iE_S1+1  # PV1, PV2, PV3, PV4
            iS_S1 = iE_S1+2  # SOM1, SOM2, SOM3, SOM4
            iV_S1 = [3]  # VIP1

            # S2
            iE_S2 = np.array([13, 17, 20, 23])  # E1, E2, E3, E4 of S2 
            iP_S2 = iE_S2+1  # PV1, PV2, PV3, PV4
            iS_S2 = iE_S2+2  # SOM1, SOM2, SOM3, SOM4
            iV_S2 = [16]  # VIP1

        # Extracting submatrices based on the defined index sets
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
        row1 = np.hstack([Wee_S1, -Wep_S1, -Wes_S1, -Wev_S1, Wee_S1S2, -Wep_S1S2, -Wes_S1S2, -Wev_S1S2])
        row2 = np.hstack([Wpe_S1, -Wpp_S1, -Wps_S1, -Wpv_S1, Wpe_S1S2, -Wpp_S1S2, -Wps_S1S2, -Wpv_S1S2])
        row3 = np.hstack([Wse_S1, -Wsp_S1, -Wss_S1, -Wsv_S1, Wse_S1S2, -Wsp_S1S2, -Wss_S1S2, -Wsv_S1S2])
        row4 = np.hstack([Wve_S1, -Wvp_S1, -Wvs_S1, -Wvv_S1, Wve_S1S2, -Wvp_S1S2, -Wvs_S1S2, -Wvv_S1S2])
        row5 = np.hstack([Wee_S2S1, -Wep_S2S1, -Wes_S2S1, -Wev_S2S1, Wee_S2, -Wep_S2, -Wes_S2, -Wev_S2])
        row6 = np.hstack([Wpe_S2S1, -Wpp_S2S1, -Wps_S2S1, -Wpv_S2S1, Wpe_S2, -Wpp_S2, -Wps_S2, -Wpv_S2])
        row7 = np.hstack([Wse_S2S1, -Wsp_S2S1, -Wss_S2S1, -Wsv_S2S1, Wse_S2, -Wsp_S2, -Wss_S2, -Wsv_S2])
        row8 = np.hstack([Wve_S2S1, -Wvp_S2S1, -Wvs_S2S1, -Wvv_S2S1, Wve_S2, -Wvp_S2, -Wvs_S2, -Wvv_S2])
        
        # Vertically stack the rows to form W0
        W0 = np.vstack([row1, row2, row3, row4, row5, row6, row7, row8])
        #print(Wee_S1.shape)
        #print(W.shape, W0.shape)

        # include the external input to the matrix 
        if include_Iext:
            #print('W0', W0.shape) # this should be 13x13
            if self.cortex_type == 'somato':
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
        if self.cortex_type == 'somato':
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

            sigmoid_params = np.vstack((sigmoid_params_S1, sigmoid_params_S2))
                                       
        elif self.cortex_type == 'visual':
            sigmoid_params = np.array([[  0.12782346,  32.10540543,  31.39696397],
                                    [  0.12782346,  32.10540543,  31.39696397],
                                    [  0.12782346,  32.10540543,  31.39696397],
                                    [  0.12782346,  32.10540543,  31.39696397],
                                    [  0.14218422,  40.03107351, 166.82960408],
                                    [  0.14218422,  40.03107351, 166.82960408],
                                    [  0.14218422,  40.03107351, 166.82960408],
                                    [  0.14218422,  40.03107351, 166.82960408],
                                    [  0.07937015,  42.01276379,  56.95305832],
                                    [  0.07937015,  42.01276379,  56.95305832],
                                    [  0.07937015,  42.01276379,  56.95305832],
                                    [  0.07937015,  42.01276379,  56.95305832],
                                    [  0.0704119 ,  37.86409387,  38.52689646],
                                    [  0.0704119 ,  37.86409387,  38.52689646],
                                    [  0.0704119 ,  37.86409387,  38.52689646],
                                    [  0.0704119 ,  37.86409387,  38.52689646]])
            
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

