import numpy as np
from pyrates.frontend import OperatorTemplate

N_cells = 16
pro_names = np.array(['PRO_E1', 'PRO_P1', 'PRO_S1', 'PRO_V1', 'PRO_E2', 'PRO_P2', 'PRO_S2', 'PRO_V2', 'PRO_E3', 'PRO_P3', 'PRO_S3', 'PRO_V3', 'PRO_E4', 'PRO_P4', 'PRO_S4', 'PRO_V4'])


# sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
sigmoid_params = np.array([ [  0.12782346,  32.10540543,  31.39696397], #Zellen in welcher Reihenfolge?
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

PRO_E1 = OperatorTemplate(
    name='PRO_E1', path=None,
    equations=["m_out = 1/(1 + exp(r*(vth-v))) - 1/(1 + exp(r*vth))"], #welche Formel mit m_max?
    variables={'m_out': 'output',
               'v': 'input',
               'r': 0.12,
               'vth': 32.10540543,
               'm_max': 31.39696397}, 
    description="sigmoidal potential-to-rate operator")

pro = np.array([PRO_E1])


from copy import deepcopy  
for i in range(1,len(pro_names)):
    pro = np.append(pro, deepcopy(PRO_E1).update_template(
        name = f'{pro_names[i]}', path=None, variables={'r': sigmoid_params[i,0], 'vth': sigmoid_params[i,1], 'm_max': sigmoid_params[i,2]}
    ))

#print(pro)

# E1, E2, E3, E4, PV1, PV2, PV3, PV4, SOM1, SOM2, SOM3, SOM4, VIP1, VIP2, VIP3, VIP4, Iext 
tau_L = np.array([6,3,20,15,6,3,20,15,6,3,20,15,6,3,20,15,3])*1e-3 # sec  #NEUE REIHENFOLGE NACH OPERATOREN

rpo_names = np.array(['RPO_E1', 'RPO_P1', 'RPO_S1', 'RPO_V1', 'RPO_E2', 'RPO_P2', 'RPO_S2', 'RPO_V2', 'RPO_E3', 'RPO_P3', 'RPO_S3', 'RPO_V3', 'RPO_E4', 'RPO_P4', 'RPO_S4', 'RPO_V4', 'RPO_INPUT'])
#rpos_up = np.tile(rpo_names, (4, 4))

RPO_E1 = OperatorTemplate(
    name='RPO_E1', path=None,
    equations=['d/dt * v = u',
               'd/dt * u = H * r - tau/tau * u - 1/tau * v'],  #Woher kommt die Gleichung?
    variables={'v': 'output',
               'u': 'variable',
               'r': 'input',
               'tau': 0.006,
               'H': 1},
    description="rate-to-potential operator")

rpo = np.array([RPO_E1])

for i in range(1, len(rpo_names)):
    rpo = np.append(rpo, deepcopy(RPO_E1).update_template(
    name = f'{rpo_names[i]}', variables={'tau': tau_L[i]} 
    ))

#print(rpo)


cells = np.array(['E1','P1','S1','V1','E2','P2','S2','V2','E3','P3','S3','V3','E4','P4','S4','V4'])

from pyrates.frontend import NodeTemplate
E1 = NodeTemplate(name="E1", path=None, operators=[PRO_E1] + list(rpo)) #alle Eingänge

pop = np.array([E1])

for i in range(1,N_cells):
    pop = np.append(pop, deepcopy(E1).update_template(
        name = f'{cells[i]}', operators=[pro[i]] + list(rpo)
    ))

#print("Pop:",pop)
#print("Operatoren von P1:", pop[1].operators)



#Spalte für Input-connections?

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
cir = CircuitTemplate(  #Schleife
    name="cir", nodes={},
    edges=[],
            path = None 
)

updated_nodes={}
for i in range(N_cells):
    updated_nodes[cells[i]] = pop[i]

cir = cir.update_template(nodes=updated_nodes)


for i, cell_i in enumerate(cells):
    for j, cell_j in enumerate(cells):
        cir = deepcopy(cir).update_template(
        edges=[(f'{cell_i}/{pro[i].name}/m_out', f'{cell_j}/{rpo[j].name}/r', None, {'weight': W[i,j]})]
)
print(cir.edges)


results = cir.run(simulation_time=2.0,
                  step_size=1e-4,
                  sampling_step_size=1e-3,
                  outputs={'V_E1': 'E1/RPO_E/v'},
                  backend='default',
                  solver='scipy')


import matplotlib.pyplot as plt
results.plot()
plt.show()
