# %%
""" Complete PyRates model with the parameters from "parameters.py" and definitions compatible with PyCobi continuations. 
- background input in all non-thalamic populations
- time-varying external input in thalE
- connectivity parameters (g/bEI)
"""
# %% libraries import
import os 
import sys
import json
from ruamel.yaml import YAML
from pyrates.frontend.template import CircuitTemplate
from pyrates.frontend.fileio.yaml import dump_to_yaml
from pyrates.frontend import OperatorTemplate, NodeTemplate, EdgeTemplate, CircuitTemplate
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import h5py
import pandas as pd
from pprint import pprint
from parameters import Parameter

# %% TODO: ???
WDDIR = os.getenv("WDDIR")

def read_simulation_params():
    """Read simulation parameters from json file."""
    # Read in preprocessing parameters
    with open(os.path.join(WDDIR, 'Simulations', 'simulation_parameter.json'), 'r') as json_file:
        params = json.load(json_file)
    
    return params

# %%
class SomatoModelPyrates():
    def __init__(self, params={}):
        # load in all connectivity parameters, time constants, etc.
        self.p = Parameter()
        self.tau = self.p.tau
        self.nPop = self.p.nPop
        # sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
        self.sigm = self.p.sigmoid_params
        self.cells = ['E3b','PV3b','SST3b','VIP3b', # A3b
        'E1S1', 'PV1S1', 'SST1S1', 'VIPS1', 'E2S1', 'PV2S1', 'SST2S1', 'E3S1', 'PV3S1', 'SST3S1', 'E4S1', 'PV4S1', 'SST4S1', # S1
        'E1S2', 'PV1S2', 'SST1S2', 'VIPS2', 'E2S2', 'PV2S2', 'SST2S2', 'E3S2', 'PV3S2', 'SST3S2', 'E4S2', 'PV4S2', 'SST4S2', # S2
        'ThalE', 'ThalI']
        self.cells_ext = self.cells + ['BACKGROUND'] + ['G'] + ['G_THAL'] + ['BEI'] + ['BEI_THAL'] 
        # model parameters
        self.pro_names = ["PRO_"+ cell for cell in self.cells]
        self.rpo_names=["RPO_"+cell for cell in self.cells]
        self.rpo_names_extended = self.rpo_names + ['RPO_Iext']
        self.connectivity_names = ["connectivity_"+ cell for cell in self.cells]

        # SIMULATION PARAMETERS - parameters that will be updated from the json file 
        # (first initialized with default values) 
        self.simulation_dur = 2 # in s
        self.step_size = 0.001 # in s
        self.sampling_step_size = 0.001 # in s
        self.input_onset = 1.001
        self.thal_connect = [0,0,0,0]
        self.simulation_dur = 5.0
        self.extI_cellcounts = 1000
        self.balance_EI = 0.7
        self.bI_cellcounts = 100
        self.thal_cellcounts = 500
        self.bEI_thal = 0.5
        self.g_thal = 2
        self.input_type = 'step'
        self.area = 'all' 
        self.coupling_strength = 10
        self.Ib_strength = 7
        self.Iext_strength = 10
        self.Iext_duration = 0.5

        # update parameters based on params dicts
        self.__dict__.update(params)

        # sigmoid parameters 
        self._extract_sigmoid_params()
        # tau parameters 
        self._load_tau_params()
        # connectivity parameters
        self.connectivity_weights()
        self.connectivity_populations()

        # create the PyRates model
        self.create_pyrates_model()

        # define filename for saving results
        self.filename = (
            f"g{self.coupling_strength}_bEI{self.balance_EI}_Ib{self.Ib_strength}_Iextd{self.Iext_duration}_"
            f"{self.input_type}Iexts{self.Iext_strength}_Ionset{self.input_onset}_thalcells{self.thal_cellcounts}_"
            f"Ibcells{self.bI_cellcounts}_Iextcells{self.extI_cellcounts}_PYRATES"
        )

    def _extract_sigmoid_params(self):
        """Extract r, v_thr, m_max from sigmoid parameter matrix."""
        N_cells = len(self.cells)
        self.r = [self.sigm[row, 0] for row in range(N_cells)]
        self.v_thr = [self.sigm[row, 1] for row in range(N_cells)]
        self.m_max = [self.sigm[row, 2] for row in range(N_cells)]

    def _load_tau_params(self):
        """Load time constants from Parameter object."""
        tau, _ = self.p.get_params()
        self.tau = tau[0, :]
    
    def connectivity_weights(self):
        self.W0, self.W_to_thal, self.W_from_thal, self.Wb, self.Wext = self.p.get_raw_connectivity(self.thal_connect, self.extI_cellcounts, self.bI_cellcounts, self.thal_cellcounts)
        self.W0 = np.append(self.W0, self.W_to_thal, axis=0)
        self.W0 = np.append(self.W0, self.W_from_thal.T, axis=1)
        self.W = np.concatenate((self.W0, self.Wb, self.Wext), axis=1)
    
    def connectivity_populations(self):
        # inhibitory
        idx_I_A3b = np.array([1,2,3])
        idx_I_S = [1,2,3,5,6,8,9,11,12]
        idx_I_S1 = np.array(idx_I_S)+4 
        idx_I_S2 = np.array(idx_I_S)+17
        self.idx_I = np.concatenate((idx_I_A3b,idx_I_S1,idx_I_S2))
        # excitatory
        idx_E_A3b = np.array([0])
        idx_E_S = [0,4,7,10]
        idx_E_S1 = np.array(idx_E_S)+4
        idx_E_S2 = np.array(idx_E_S)+17
        self.idx_E = np.concatenate((idx_E_A3b,idx_E_S1,idx_E_S2))


    def create_Iext(self):
        """
        Creates external input.
        THIS IS ONLY FOR PLOTTING REASONS. The actual input is computed via operators.
        """

        Iext = np.zeros(int(self.simulation_dur / self.step_size))

        if self.input_type == "step":
            t = int(self.Iext_duration / self.step_size)
            t0 = int(self.input_onset / self.step_size)
            Iext[t0 : t0 + t] = self.Iext_strength
        elif self.input_type == "background":
            # provide input for the entire simulation duration
            Iext[:] = self.Iext_strength

        return Iext

    def create_Ibackground(self):
        """
        Create Background Input
        THIS IS ONLY FOR PLOTTING REASONS. The actual input is computed via operators.
        """
        
        Ib = np.zeros(int(self.simulation_dur / self.step_size))
        Ib[:] = self.Ib_strength
        return Ib

    def create_pyrates_model(self):
        "creates the complete PyRates model"
        # Operator template for the PRO
        # no background input:
        N_cells = len(self.cells)
        pro = OperatorTemplate(
            name='PRO', path=None,
            equations=["m_outC = m_max / (1 + exp(r*(V_thr - v)))"],
            variables={'m_outC': 'output',
               'v': 'input',
               'r': 0.1,
               'V_thr': 35.0,
               'm_max': 70.0},
            description="sigmoidal potential-to-rate operator")
        # background input:
        pro_bI = OperatorTemplate(
        name='PRO_bI', path=None,
        equations=["m_outC = m_max / (1 + exp(r*(V_thr - (v+v_bIn))))"],
        variables={'m_outC': 'output',
               'v': 'input',
               'v_bIn': 'input',
               'r': 0.1,
               'V_thr': 35.0,
               'm_max': 70.0},
        description="sigmoidal potential-to-rate operator")

        pros = [
            (   # background input: all populations except thalamus
                deepcopy(pro_bI).update_template(name=self.pro_names[i], variables={'r': self.r[i],
                                                 'V_thr': self.v_thr[i],
                                                 'm_max': self.m_max[i]})
                if i < N_cells - 2
                # no background input: thalamus
                else deepcopy(pro).update_template(name=self.pro_names[i], variables={'r': self.r[i],
                                                 'V_thr': self.v_thr[i],
                                                 'm_max': self.m_max[i]})
            )
            for i in range(N_cells)
            ]   

        # Operator template for the RPO
        rpo = OperatorTemplate(
        name='RPO', path=None,
        equations=['d/dt * v = i',
               'd/dt * i = H/tau * (m_in ) - 2 * i/tau - v/tau^2'],
        variables={'v': 'output',
               'i': 'variable',
               'm_in': 'input',
               'tau': 0.01,
               'H': 1.0},
        description="excitatory rate-to-potential operator")

        rpos = [deepcopy(rpo).update_template(name=self.rpo_names[i], variables={"tau": self.tau[i]}) for i in range(N_cells)]

        # Operator template for the background input --> only for the "BACKGROUND" population!
        rpo_bI = OperatorTemplate(
            name='RPO_bI', path=None,
            equations=['d/dt * v_bI = i',
               'd/dt * i = H/tau * bI_cellcount * (bI) - 2 * i/tau - v_bI/tau^2'],
            variables={'v_bI': 'output',
               'i': 'variable',
               'bI': f'input({self.Ib_strength})',  # external background input
               'bI_cellcount': self.bI_cellcounts,
               'tau': 0.003,
               'H': 1.0},
            description="excitatory rate-to-potential operator-background input")

        # also create an Ib array for plotting later
        Ib = self.create_Ibackground()
        self.Ib = np.tile(Ib, (self.nPop,1))

        # Operator template for the external input --> only for thalE!
        rpo_Iext_thalE = [OperatorTemplate(
        name='RPO_Iext', path=None,
        equations=['d/dt * v = i',
               'd/dt * i = H/tau * Iext_cellcounts * (Iext) - 2 * i/tau - v/tau^2'],
        variables={'v': 'output',
               'i': 'variable',
               'Iext': 'input',
               'Iext_cellcounts': self.extI_cellcounts,
               'tau': self.tau[N_cells-2],
               'H': 1.0}
            ) 
        ]
        
        # definition of the external input
        create_I_ext = [OperatorTemplate(
            name="InputOp", path=None,
            equations=[
            # step input: A during [onset, onset+dur), else 0
                "Iext = (A/2)*(1+sign(t-onset)) - (A/2)*(1+sign(t-(onset+dur)))"
                ],
            variables={
                "Iext": "output",
                "t": "variable",
                "A": self.Iext_strength,
                "onset": self.input_onset/self.step_size,
                "dur": self.Iext_duration/self.step_size
                # TO UNCOMMENT FOR PYCOBI:
                #"A": float(input_strength),
                #"onset": float(input_onset),
                #"dur": float(self.Iext_duration)
                },
            description="External step input"
            )]

        # also create Iext array for later plotting 
        Iext = self.create_Iext()
        self.Iext = np.tile(Iext, (self.nPop,1))
        
        # Operators and nodes for the connectivity parameters
        g_definition = OperatorTemplate(
            name="g_definition",
            equations="gC = g_input",
            variables={"gC": "output", "g_input": f"input({float(self.coupling_strength)})"}
        )
        bEI_definition = OperatorTemplate(
            name="bEI_definition",
            equations="bEIC = bEI_input",
            variables={"bEIC": "output", "bEI_input": f"input({float(self.balance_EI)})"}
        )       
        g_thal_definition = OperatorTemplate(
            name="g_thal_definition",
            equations="g_thalC = g_thal_input",
            variables={"g_thalC": "output", "g_thal_input": f"input({float(self.g_thal)})"}
        )   
        bEI_thal_definition = OperatorTemplate(
            name="bEI_thal_definition",
            equations="bEI_thalC = bEI_thal_input",
            variables={"bEI_thalC": "output", "bEI_thal_input": f"input({float(self.bEI_thal)})"}
        )

        # Connectivity operators
        connectivityE = OperatorTemplate(name="connectivityE", 
                           equations= "m_out = (g*bEI)*m_outC",
                           variables={
                            "m_out": "output", 
                            #"g": g,
                            #"bEI":bEI,
                            "g": "input",
                            "bEI": "input",
                            "m_outC":"input"
                            #"connect_reverse_factor": f"input({float(connect_reverse_factor)})"
                            }
                            )
        connectivityI = OperatorTemplate(name="connectivityI", 
                           equations= "m_out = ((-g)*(1-bEI))*m_outC",
                           variables={
                            "m_out": "output", 
                            #"g": -g,
                            #"bEI":bEI,
                            "g": "input",
                            "bEI": "input",
                            #"connect_reverse_factor":f"input({float(connect_reverse_factor)})",
                            "m_outC":"input"
                            }
                            )
        connectivityE_thal = OperatorTemplate(name="connectivityE_thal", 
                           equations= "m_out = (g_thal*bEI_thal)*m_outC",
                           variables={
                            "m_out": "output", 
                            #"g_thal": g_thal,
                            #"bEI_thal":bEI_thal,
                            "g_thal": "input",
                            "bEI_thal": "input",
                            #"connect_reverse_factor_thal":f"input({float(connect_reverse_factor)})",
                            "m_outC":"input"
                            }
                            )
        connectivityI_thal = OperatorTemplate(name="connectivityI_thal", 
                           equations= "m_out = ((-g_thal)*(1-bEI_thal))*m_outC",
                           variables={
                            "m_out": "output", 
                            #"g_thal": -g_thal,
                            #"bEI_thal":bEI_thal,
                            "g_thal": "input",
                            "bEI_thal": "input",
                            #"connect_reverse_factor_thal":f"input({float(connect_reverse_factor)})",
                            "m_outC":"input"
                            }
                            )
        connectivity = []
        for i in range(N_cells):
            if i in self.idx_E:  # excitatory populations
                connectivity.append(deepcopy(connectivityE).update_template(name=self.connectivity_names[i]))
            elif i in self.idx_I:  # inhibitory populations
                connectivity.append(deepcopy(connectivityI).update_template(name=self.connectivity_names[i]))
            elif i == N_cells - 2: # thalE
                connectivity.append(deepcopy(connectivityE_thal).update_template(name=self.connectivity_names[i]))
            else: # thalI
                connectivity.append(deepcopy(connectivityI_thal).update_template(name=self.connectivity_names[i]))

        # Node templates
        nodes = [
            NodeTemplate(
                name=self.cells_ext[i],
                path=None,
                operators=(
                    ([pros[i]] + [connectivity[i]] + rpos + 
                    (rpo_Iext_thalE + create_I_ext if i == N_cells - 2 else [])
                    ) if i < N_cells 
                    #else [rpo_bI] if i == N_cells
                    else [g_definition] if i == N_cells + 1 # G
                    else [g_thal_definition] if i == N_cells + 2 # G_THAL
                    else [bEI_definition] if i == N_cells + 3 # BEI
                    else [bEI_thal_definition] if i == N_cells + 4 # BEI_THAL
                    else [rpo_bI] # operator for the background input
                )
            )
            for i in range(len(self.cells_ext))
        ]
        # Edges
        edges=[]
        # i : target 
        for i, cell_i in enumerate(self.cells):
            if cell_i not in ['ThalE', 'ThalI']:
                edges.append(('BACKGROUND/RPO_bI/v_bI', f'{cell_i}/{self.pro_names[i]}/v_bIn', None, {'weight': 1.0})) # background input
                edges.append(('G/g_definition/gC', f'{cell_i}/{self.connectivity_names[i]}/g', None, {'weight': 1.0})) # G
                edges.append(('BEI/bEI_definition/bEIC', f'{cell_i}/{self.connectivity_names[i]}/bEI', None, {'weight': 1.0})) # BEI
            else:
                edges.append(('G_THAL/g_thal_definition/g_thalC', f'{cell_i}/{self.connectivity_names[i]}/g_thal', None, {'weight': 1.0})) # G_THAL
                edges.append(('BEI_THAL/bEI_thal_definition/bEI_thalC', f'{cell_i}/{self.connectivity_names[i]}/bEI_thal', None, {'weight': 1.0})) # BEI_THAL
            for j, cell_j in enumerate(self.cells):
                edges.append((f'{cell_j}/{self.connectivity_names[j]}/m_out', f'{cell_i}/{self.rpo_names[j]}/m_in', None, {'weight': self.W[i,j]}))
                    
        # Set up the Model Circuit 
        self.model = CircuitTemplate(
            name = 'model',
            nodes = {name: node for name, node in zip(self.cells_ext, nodes)},
            edges = edges,
            path = None)

        # TODO: add the optional saving to a .yaml file 
        #circuit_to_yaml(model, "model.yaml")


    def simulate(self):
        """run the simulation"""
        N_cells = len(self.cells)
        outputs = {}
        b_inputs = []
        for i, target_cell in enumerate(self.cells):
            if i == (N_cells-2):
                for rpo_name_ext in self.rpo_names_extended[:N_cells+1]:
                    outputs[f'V_{target_cell}/{rpo_name_ext}'] = f'{target_cell}/{rpo_name_ext}/v'
            else:
                for rpo_name in self.rpo_names[:N_cells]:
                    outputs[f'V_{target_cell}/{rpo_name}'] = f'{target_cell}/{rpo_name}/v'

        outputs['V_background/RPO_bI'] = 'BACKGROUND/RPO_bI/v_bI'
        results = self.model.run(simulation_time=self.simulation_dur,
                  step_size=self.step_size,
                  sampling_step_size=self.sampling_step_size,
                  outputs=outputs,
                  backend ="scipy",
                  vectorize=False,
                  clear=False,
                  float_precision="float64"
                  )
        all_potentials = []
        for i, cell in enumerate(self.cells):
                
            # ThalE and ThalI do not receive background input
            if cell in ['ThalE', 'ThalI']:   
                # thalE receives external input from RPO_Iext
                if cell == 'ThalE': 
                    potential_keys = [f'V_{cell}/{rpo}' for rpo in self.rpo_names_extended]
                    sources = results[potential_keys]
                    
                    # append background input potential as zero array
                    sources = np.array(sources)
                    zero_array = np.zeros((sources.shape[0], 1))
                    sources = np.hstack((sources, zero_array))
                else:
                    potential_keys = [f'V_{cell}/{rpo}' for rpo in self.rpo_names]
                    sources = results[potential_keys]
                    zero_array = np.zeros((sources.shape[0], 1))
                    # append background input potential and external input as zero array
                    sources = np.hstack((sources, zero_array))
                    sources = np.hstack((sources, zero_array))  
            
            else: 
                potential_keys = [f'V_{cell}/{rpo}' for rpo in self.rpo_names]
                potential_keys += [f'V_background/RPO_bI']
                sources = results[potential_keys]

                # append external input potential as zero array
                sources = np.array(sources)
                zero_array = np.zeros(((sources.shape[0]),1))
                sources = np.hstack((sources, zero_array))

            all_potentials.append(sources)

        all_potentials = np.array(all_potentials)
        all_potentials = np.rollaxis(all_potentials, 2, 1) # shape: target cells x source cells x timepoints

        # rates:
        summed_potentials = np.sum(all_potentials, axis=1)

        n_cells, n_timepoints  = summed_potentials.shape
        m_out_all = np.zeros((n_cells, n_timepoints))  # NumPy array, not list!

        for i, _ in enumerate(self.cells):
            m_out_all[i, :] =  self.m_max[i] / (
                1 + np.exp(self.r[i] * (self.v_thr[i] - summed_potentials[i, :]))
            )

        self.potential = all_potentials
        self.rate = m_out_all


    def circuit_to_yaml(self, circuit: CircuitTemplate, path: str):
        """
        Save a Pyrates CircuitTemplate to YAML safely, converting NumPy types to Python types.
        """
        # Patch ruamel.yaml to handle numpy scalars
        yaml = YAML()
        def represent_numpy_scalar(dumper, data):
            return dumper.represent_float(float(data))
        yaml.representer.add_representer(np.float64, represent_numpy_scalar)
        yaml.representer.add_representer(np.float32, represent_numpy_scalar)
        yaml.representer.add_representer(np.int32, lambda d, x: d.represent_int(int(x)))
        yaml.representer.add_representer(np.int64, lambda d, x: d.represent_int(int(x)))

        dump_to_yaml(circuit, path=path)
        print(f"CircuitTemplate successfully saved to {path}")


    def save_results_csv(self, filedir, filename, full=False):
        """
        Safe the simulated data in a csv file
        """

        filename = filename + ".hdf5"

        # only safe every second datapoint
        resolution_tstep = 0.01
        print("tstep resolution", resolution_tstep)
        rates_downsampled = self.rate[:, :: int(1000 * resolution_tstep)]
        # TODO: double check if the potential array needs to be transposed 
        rates_df = pd.DataFrame(rates_downsampled.T, columns=self.cells)
        rates_df.to_hdf(
            os.path.join(filedir, filename), index=False, key="rates", mode="a"
        )

        # sum the potentials together and save them
        potential_sum = np.sum(self.potential, axis=1)
        potential_sum_downsampled = potential_sum[:, :: int(1000 * resolution_tstep)]
        potential_df = pd.DataFrame(potential_sum_downsampled.T, columns=self.cells)
        potential_df.to_hdf(
            os.path.join(filedir, filename), index=False, key="summed_potential", mode="a"
        )

        if full:
            # save all potentials additionally
            psp_filename = "full_" + filename
            print('full potential file:', psp_filename)
            self.write_3D_csv(os.path.join(filedir, psp_filename))


    def write_3D_csv(self, filename):
        """
        Write results in form of a 3D hdf5 file.
        """
        dataset_name = 'full_potentials'

        with h5py.File(filename, "w") as f:
            f.create_dataset(dataset_name, data=self.potential, compression="gzip")
        
# TODOs:
# - after simulation we assign the attributes rate and potential to the class. The attribute "model.potential" should be a numpy array of 3 dimensions 
# - in "save_results" also implement to save the connectivity parameter