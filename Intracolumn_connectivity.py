import numpy as np
from pyrates.frontend import OperatorTemplate

N_cells = 16
pro_variables = np.tile(np.array(['pro_E', 'pro_PV', 'pro_SOM', 'pro_VIP']),4)
pros = np.array(['PRO_E', 'PRO_PV', 'PRO_SOM', 'PRO_VIP'])
pros = np.tile(pros,4)

# sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
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

#values for E from Chien
pro_E = OperatorTemplate(
    name='PRO_E', path=None,
    equations=["m_out = 1/(1 + exp(r*(vth-v))) - 1/(1 + exp(r*vth))"], #welche Formel mit m_max?
    variables={'m_out': 'output',
               'v': 'input',
               'r': 0.12,
               'vth': 32.10540543,
               'm_max': 31.39696397}, 
    description="sigmoidal potential-to-rate operator")


from copy import deepcopy  
for i in range(1,N_cells):
    new_template = pro_variables[i]
    new_template = deepcopy(pro_E).update_template(
        name = f'{pros[i]}', path=None, variables={'r': sigmoid_params[i,0], 'vth': sigmoid_params[i,1], 'm_max': sigmoid_params[i,2]}
    )
    

# E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1, VIP2, VIP3, VIP4, Iext
tau_L = np.array([6,6,6,6,3,3,3,3,20,20,20,20,15,15,15,15,3])*1e-3 # sec

rpos = np.array([['RPO_EE','RPO_E_PV','RPO_E_SOM','RPO_E_VIP'],
                 ['RPO_PV_E','RPO_PV_PV','RPO_PV_SOM','RPO_PV_VIP'],
                 ['RPO_SOM_E','RPO_SOM_PV','RPO_SOM_SOM','RPO_SOM_VIP'],
                 ['RPO_VIP_E','RPO_VIP_PV','RPO_VIP_SOM','RPO_VIP_VIP']])
rpos_up = np.tile(rpos, (4, 4))

# values for E -> E AMPA
rpo_E_E = OperatorTemplate(
    name='RPO_EE', path=None,
    equations=['d/dt * v = u',
               'd/dt * u = H * r - (tau1+tau2)/(tau1*tau2) * u - 1/(tau1*tau2) * v'], #was ist mit w(t) in der Gleichung im Paper?
    variables={'v': 'output',
               'u': 'variable',
               'r': 'input',
               'tau': 0.006,
               'H': 1},
    description="rate-to-potential operator")

for i in range(4):
    for j in range(4):
        new_template = rpos[i,j]
        new_template = deepcopy(rpo_E_E).update_template(
            name = f'{rpos[i,j]}', path=None, variables={'tau': tau_L[i]} #stimmt die Reihenfolge der tau-values?
        )


from pyrates.frontend import NodeTemplate
pop_E = NodeTemplate(name="E1", path=None, operators=[pro_E, rpo_E_E, rpo_PV_E, rpo_SOM_E, rpo_VIP_E]) #alle Eingänge
pop_PV = NodeTemplate(name="P1", path=None, operators=[pro_PV, rpo_E_PV, rpo_PV_PV, rpo_SOM_PV, rpo_VIP_PV])
pop_SOM = NodeTemplate(name='S1', path=None, operators=[pro_SOM, rpo_E_SOM, rpo_PV_SOM, rpo_SOM_SOM, rpo_VIP_SOM])
pop_VIP = NodeTemplate(name='V1', path=None, operators=[pro_VIP, rpo_E_VIP, rpo_PV_VIP, rpo_SOM_VIP, rpo_VIP_VIP])



cells = np.array(['E1','P1','S1','V1','E2','P2','S2','V2','E3','P3','S3','V3','E4','P4','S4','V4'])


# Connection probabilies
# E1, PV1, SST1, VIP1, E2, PV2, SST2, VIP2, E3, PV3, SST3, VIP3, E4, PV4, SST4, VIP4
W = np.array([[0.110000000000000,	0.494409448818898,	0.319464566929134,	0.0836692913385827,	0.170000000000000,	0.532440944881890,	0.463984251968504,	0,	0, 0.380314960629921, 0.220582677165354, 0, 0,	0,	0,	0],
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

from pyrates.frontend import CircuitTemplate
cir = CircuitTemplate(
    name="cir", nodes={'E1': pop_E, 'E2': pop_E, 'E3': pop_E, 'E4': pop_E, 'P1': pop_PV, 'P2': pop_PV, 'P3': pop_PV, 'P4': pop_PV, 'S1': pop_SOM, 'S2': pop_SOM, 'S3': pop_SOM, 'S4': pop_SOM, 'V1': pop_VIP, 'V2': pop_VIP, 'V3': pop_VIP, 'V4': pop_VIP},
    edges=[],
            path = None 
)

for i in range(len(cells)):
    for j in range(len(cells)):
        cir = deepcopy(cir).update_template(
        edges=[(f'{cells[i]}/{pros[i]}/m_out', f'{cells[j]}/{rpos_up[i,j]}/r', None, {'weight': W[i,j]})]
) #Wie kann ich kontrollieren ob das stimmt?



results = cir.run(simulation_time=2.0,
                  step_size=1e-4,
                  sampling_step_size=1e-3,
                  outputs={'V_E1': 'E1/RPO_E_E/v'},
                  backend='default',
                  solver='scipy')



import matplotlib.pyplot as plt
results.plot()
plt.show()
