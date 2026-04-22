import numpy as np
import os
import h5py
import yaml
import json
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from parameters import Parameter



WDDIR = os.getenv("WDDIR")
SIMDIR = os.getenv("SIMDIR")

helper_path = os.path.join(WDDIR, 'Analysis')
sys.path.insert(0, helper_path)
import helper_functions as hf

def read_simulation_params():
    """Read simulation parameters from json file."""
    # Read in preprocessing parameters
    with open(os.path.join(WDDIR, 'Simulations', 'simulation_parameter.json'), 'r') as json_file:
        params = json.load(json_file)
    
    return params

def read_analysis_params():
    analysis_params = hf.load_parameters(WDDIR)
    return analysis_params

class SomatoModel():

    def __init__(self, params={}):
        
        # load in all connectivity parameters, time constants, etc.
        self.p = Parameter()
        self.tau = self.p.tau
        self.nPop = self.p.nPop
        # sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
        self.sigm = self.p.sigmoid_params

        # parameters that will be updated from the json file 
        # (first initialized with default values) 
        self.simulation_dur = 2 # in s
        self.step_size = 0.001 # in s
        self.resolution_tstep = 0.001 # in s
        self.input_onset = 1.001
        self.thal_connect = [0,0,0,0]
        self.extI_cellcounts = 1000
        self.strength_I = 0 #0.7
        self.bI_cellcounts = 100
        self.thal_cellcounts = 500
        self.sI_thal = 0.5
        self.g_thal = 2
        self.input_type = 'step'
        self.area = 'all' 
        self.coupling_strength = 10
        self.Ib_strength = 7
        self.Iext_strength = 10
        self.Iext_duration = 0.5

        # scaling the coupling strength between the cortical areas
        self.g_intercortical = 1

        # update parameters based on params dicts
        self.__dict__.update(params)

        # create input array
        Iext = self.create_Iext()
        Ib = self.create_Ibackground()
        self.gE = self.coupling_strength
        self.gI = self.coupling_strength * self.strength_I
        self.gEthal = self.g_thal 
        self.gIthal = self.g_thal * self.sI_thal

        # Synaptic kernel 
        self.H = np.ones((self.nPop, self.nPop+1))

        # define time steps 
        self.steps = np.arange(self.step_size, self.simulation_dur+self.step_size, self.step_size)

        # extend input arrays
        self.Iext = np.tile(Iext, (self.nPop,1))
        self.Ib = np.tile(Ib, (self.nPop,1))

        self.filename = (
            f"gthal{self.g_thal}_sIthal{self.sI_thal}_g{self.coupling_strength}_sI{self.strength_I}_Ib{self.Ib_strength}_Iextd{self.Iext_duration}_"
            f"{self.input_type}Iexts{self.Iext_strength}_Ionset{self.input_onset}_thalcells{self.thal_cellcounts}_"
            f"Ibcells{self.bI_cellcounts}_Iextcells{self.extI_cellcounts}_gInter{self.g_intercortical}_thalUncon"
        )

        # Output matrices to store computed values for rates & potentials (E, IIN , EIN) 
        self.rate = np.zeros((self.nPop, len(self.steps)))
        self.potential = np.zeros((self.nPop, self.nPop+2, len(self.steps))) 

        # Simulation loop
        # Initialize first values for the potential, rate and first order derivative with 0 or randomly
        self.v_current = np.zeros((self.nPop, self.nPop+2)) # +2 because 1 for background input and one for external input 
        self.rate_current = np.zeros(self.nPop)
        self.u_t = np.zeros((self.nPop, self.nPop+2)) # the initial first-order derivative: v'(t) = u(t)
        self.t = 0.0

        # Weight matrix [to x from]
        self.W = self.p.get_connectivity(self.g_intercortical, self.gE, self.gI, self.gEthal, self.gIthal, self.thal_connect, self.extI_cellcounts, self.bI_cellcounts, self.thal_cellcounts, area=self.area) 


    def initialize_state(self):
        """
        Reset the dynamic state for interactive simulations.
        """
        self.v_current = np.zeros((self.nPop, self.nPop + 2))
        self.rate_current = np.zeros(self.nPop)
        self.u_t = np.zeros((self.nPop, self.nPop + 2))
        self.t = 0.0
        self.rate = np.zeros((self.nPop, len(self.steps)))
        self.potential = np.zeros((self.nPop, self.nPop + 2, len(self.steps)))


    def apply_params(self, params: dict):
        """
        Update parameters and recompute derived state.
        """
        self.__dict__.update(params)

        # recompute inputs and gains
        Iext = self.create_Iext()
        Ib = self.create_Ibackground()
        self.Iext = np.tile(Iext, (self.nPop, 1))
        self.Ib = np.tile(Ib, (self.nPop, 1))

        self.gE = self.coupling_strength 
        self.gI = self.coupling_strength * self.strength_I
        self.gEthal = self.g_thal 
        self.gIthal = self.g_thal * self.sI_thal

        # update connectivity with new gains and counts
        self.W = self.p.get_connectivity(
            self.g_intercortical,
            self.gE,
            self.gI,
            self.gEthal,
            self.gIthal,
            self.thal_connect,
            self.extI_cellcounts,
            self.bI_cellcounts,
            self.thal_cellcounts,
            area=self.area
        )


    def create_Iext(self):
        """Creates external input."""

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
        """Create Background Input"""
        Ib = np.zeros(int(self.simulation_dur / self.step_size))
        Ib[:] = self.Ib_strength
        return Ib

        
    def save_to_yaml(self, filename):
        
        S = self.p.get_connectStrength()
        P = self.p.get_connectProb()
        C = self.p.get_cellcounts()
        W = self.p.get_connectivity(self.g_intercortical,self.gE, self.gI, self.gEthal, self.gIthal, self.thal_connect, self.extI_cellcounts, self.bI_cellcounts, self.thal_cellcounts)

        # Convert numpy arrays to lists
        parameters = {
            'gE': self.gE,
            'gI': self.gI,
            'gEthal': self.gEthal, 
            'gIthal': self.gIthal,
            'S': S.tolist(),
            'P': P.tolist(),
            'C': C.tolist(),
            'W': W.tolist()
        }

        # Save parameters to a YAML file
        with open(filename + '.yaml', 'w') as file:
            yaml.dump(parameters, file)

    def plot_W_heatmap(self):
        """
        Plot connectivity matrix as heatmap.
        """
        
        W = self.p.get_connectivity(self.g_intercortical, self.gE, self.gI, self.gEthal, self.gIthal, self.thal_connect, self.extI_cellcounts, self.bI_cellcounts, self.thal_cellcounts, area=self.area) 
        sns.heatmap(W, annot=False, cmap='coolwarm', center=0, xticklabels=True, yticklabels=True)


    def simulate(self):
        '''
        Simulation loop
        '''

        last_step = self.steps[-1] 

        for timestep, _ in enumerate(self.steps):

            # compute simulation step         
            self.potential[:, :, timestep] = self.v_current.copy()
            self.rate_current = self.compute_rates()
            self.rate[:, timestep] = self.rate_current
            self.v_current = self.compute_potentials(timestep)

        print('finished loop...')

    def compute_rates(self):
        """
        Compute the firing rates of all populations.
        """
        v_sum = np.sum(self.v_current, axis=1)

        self.rate_current = self.sigm[:,2] / (
            1 + np.exp(self.sigm[:,0] * (self.sigm[:,1] - v_sum))
        )

        return self.rate_current.copy()

    def compute_potentials(self, timestep):
        """
        Compute the potentials of all populations.
        Also take into account the background input and external input.
        """

        pop_slice = slice(0, self.nPop)

        # store previous potentials (needed for correct Euler update)
        v_prev = self.v_current.copy()

        # -----------------------------
        # POPULATION INTERACTIONS
        # -----------------------------

        # update potentials
        self.v_current[:, pop_slice] += (
            self.u_t[:, pop_slice] * self.step_size
        )

        # synaptic drive
        drive = self.W[:, pop_slice] * self.rate_current

        u_dot = (
            (self.H[:, pop_slice] / self.tau[:, pop_slice]) * drive
            - 2 * self.u_t[:, pop_slice] / self.tau[:, pop_slice]
            - v_prev[:, pop_slice] / (self.tau[:, pop_slice] ** 2)
        )

        self.u_t[:, pop_slice] += u_dot * self.step_size


        # -----------------------------
        # EXTERNAL INPUT
        # -----------------------------

        # update potentials
        v_dot = self.u_t[:, -1]
        self.v_current[:, -1] += v_dot * self.step_size

        u_dot = (
            (self.H[:, -1] / self.tau[:, -1])
            * (self.W[:, -1] * self.Iext[:, timestep])
            - 2 * self.u_t[:, -1] / self.tau[:, -1]
            - v_prev[:, -1] / (self.tau[:, -1] ** 2)
        )

        self.u_t[:, -1] += u_dot * self.step_size


        # -----------------------------
        # BACKGROUND INPUT
        # -----------------------------

        # update potentials
        v_dot = self.u_t[:, -2]
        self.v_current[:, -2] += v_dot * self.step_size

        u_dot = (
            (self.H[:, -2] / self.tau[:, -2])
            * (self.W[:, -2] * self.Ib[:, timestep])
            - 2 * self.u_t[:, -2] / self.tau[:, -2]
            - v_prev[:, -2] / (self.tau[:, -2] ** 2)
        )

        self.u_t[:, -2] += u_dot * self.step_size

        return self.v_current.copy()


    def simulate_step(self, timestep):
        """
        Only for interactive computation
        """
        self.rate_current = self.compute_rates()
        self.v_current = self.compute_potentials(timestep)

        self.t += self.step_size

        return self.rate_current.copy()


    def compute_ecds():
        raise NotImplementedError

    def get_population_labels(self):
        return np.array([
            "E3b",
            "PV3b",
            "SST3b",
            "VIP3b",
            "E1",
            "PV1",
            "SST1",
            "VIP1",
            "E2",
            "PV2",
            "SST2",
            "E3",
            "PV3",
            "SST3",
            "E4",
            "PV4",
            "SST4",
            "E1S2",
            "PV1S2",
            "SST1S2",
            "VIP1S2",
            "E2S2",
            "PV2S2",
            "SST2S2",
            "E3S2",
            "PV3S2",
            "SST3S2",
            "E4S2",
            "PV4S2",
            "SST4S2",
            "ThalE",
            "ThalI",
        ])

    def get_population_spectrum_groups(self):
        return [
            ("A3b Layer 2/3", ["E3b", "PV3b", "SST3b", "VIP3b"]),
            ("S1 Layer 2/3", ["E1", "PV1", "SST1", "VIP1"]),
            ("S1 Layer 4", ["E2", "PV2", "SST2"]),
            ("S1 Layer 5", ["E3", "PV3", "SST3"]),
            ("S1 Layer 6", ["E4", "PV4", "SST4"]),
            ("S2 Layer 2/3", ["E1S2", "PV1S2", "SST1S2", "VIP1S2"]),
            ("S2 Layer 4", ["E2S2", "PV2S2", "SST2S2"]),
            ("S2 Layer 5", ["E3S2", "PV3S2", "SST3S2"]),
            ("S2 Layer 6", ["E4S2", "PV4S2", "SST4S2"]),
            ("Thalamus", ["ThalE", "ThalI"]),
        ]

    def prepare_dataframes(self):

        cells = self.get_population_labels()

        # only safe every X datapoint
        print("tstep resolution", self.resolution_tstep)
        rates_downsampled = self.rate[:, :: int(1000 * self.resolution_tstep)]
        rates_df = pd.DataFrame(rates_downsampled.T, columns=cells)

        # sum the potentials together and save them
        potential_sum = np.sum(self.potential, axis=1)
        potential_sum_downsampled = potential_sum[:, :: int(1000 * self.resolution_tstep)]
        potential_df = pd.DataFrame(potential_sum_downsampled.T, columns=cells)

        return rates_df, potential_df

    def save_results_csv(self, filedir, filename, full=False, save_params=False):
        """
        Safe the simulated data in a csv file
        """

        rates_df, potential_df = self.prepare_dataframes()

        filename = filename + ".hdf5"
        
        rates_df.to_hdf(
            os.path.join(filedir, filename), index=False, key="rates", mode="a"
        )

        potential_df.to_hdf(
            os.path.join(filedir, filename), index=False, key="summed_potential", mode="a"
        )

        if full:
            # save all potentials additionally
            psp_filename = "full_" + filename
            print('full potential file:', psp_filename)
            self.write_3D_csv(os.path.join(filedir, psp_filename))

        if save_params:
            # safe connectivty parameter in yaml file
            self.save_to_yaml(os.path.join(filedir, "params" + self.filename))


    def write_3D_csv(self, filename):
        """
        Write results in form of a 3D hdf5 file.
        """
        dataset_name = 'full_potentials'

        with h5py.File(filename, "w") as f:
            f.create_dataset(dataset_name, data=self.potential, compression="gzip")
    
    def compute_late_longterm_spectrum(self):
        """
        Compute late-longterm frequency spectra for all populations.
        """
        analysis_params = read_analysis_params()
        sampling_params = analysis_params['sampling']

        rates_df, potentials_df = self.prepare_dataframes()
        df = pd.DataFrame()

        hf.compute_longeterm_late(
            df,
            rates_df,
            potentials_df,
            self.input_onset,
            self.Iext_duration,
            self.step_size,
            sampling_params['sample_delay_late'],
            sampling_params['sample_dur']
        )

        start_sample_late = int(
            (self.input_onset + self.Iext_duration + sampling_params['sample_delay_late']) / self.step_size
        )
        stop_sample_late = int(start_sample_late + sampling_params['sample_dur'] / self.step_size)

        spectra, freqs = hf.compute_window_frequency(
            df,
            rates_df,
            potentials_df,
            start_sample_late,
            stop_sample_late,
            "lateLongterm",
            self.step_size,
            sampling_params['rate_osc_threshold'],
            sampling_params['potential_osc_threshold'],
            compute_spectrum=True
        )

        return spectra, freqs

    def save_frequency_spectra(self, filedir, filename=None, spectra=None, freqs=None, window_prefix="lateLongterm"):
        """
        Save frequency spectra to an HDF5 file with simulation metadata.
        """
        os.makedirs(filedir, exist_ok=True)

        if spectra is None or freqs is None:
            spectra, freqs = self.compute_late_longterm_spectrum()

        if spectra.ndim != 2:
            raise ValueError("spectra must be 2D (n_populations x n_frequencies)")

        population_labels = self.get_population_labels()
        if spectra.shape[0] != len(population_labels):
            raise ValueError(
                f"spectra rows ({spectra.shape[0]}) do not match number of populations ({len(population_labels)})"
            )

        if filename is None:
            filename = (
                f"spectrum_ginter{self.g_intercortical}_"
                f"g{self.coupling_strength}_sI{self.strength_I}_Ib{self.Ib_strength}_"
                f"Iextd{self.Iext_duration}_{self.input_type}Iexts{self.Iext_strength}.hdf5"
            )
        elif not filename.endswith(".hdf5"):
            filename = filename + ".hdf5"

        filepath = os.path.join(filedir, filename)
        with h5py.File(filepath, "w") as h5f:
            h5f.create_dataset("freqs", data=freqs)
            h5f.create_dataset("spectra", data=spectra)
            h5f.create_dataset("population_labels", data=np.asarray(population_labels, dtype="S32"))

            h5f.attrs["window_prefix"] = window_prefix
            h5f.attrs["g_intercortical"] = self.g_intercortical
            h5f.attrs["coupling_strength"] = self.coupling_strength
            h5f.attrs["strength_I"] = self.strength_I
            h5f.attrs["Ib_strength"] = self.Ib_strength
            h5f.attrs["Iext_strength"] = self.Iext_strength
            h5f.attrs["Iext_duration"] = self.Iext_duration
            h5f.attrs["step_size"] = self.step_size
            h5f.attrs["input_onset"] = self.input_onset
            h5f.attrs["area"] = self.area
            h5f.attrs["input_type"] = self.input_type

        return filepath


    def analyse_signal(self, save_spectrum=False):
        """
        Get the following aspects from the signal:
        - oscillation frequency peak
        - oscillation yes/no
        """
        spectra, freqs = self.compute_late_longterm_spectrum()


        self.plot_freq_spectrum(spectra, freqs)
        self.plot_freq_spectrum_all_populations(spectra, freqs)

        if save_spectrum:
            spectrum_dir = os.path.join(SIMDIR, "spectrum_results")
            path = self.save_frequency_spectra(spectrum_dir, spectra=spectra, freqs=freqs)
            print(f"saved spectra: {path}")


    def plot_freq_spectrum(self, spectra, freqs, pop_idx=0, pop_name=None, max_freq_hz=100):
        """
        Plot frequency spectrum for a single population.

        Parameters
        ----------
        spectra : np.ndarray
            Array of shape (n_populations, n_frequencies)
        freqs : np.ndarray
            Frequency vector (Hz)
        pop_idx : int
            Index of population to plot
        pop_name : str (optional)
            Name of the population (for title)
        """

        if spectra.ndim != 2:
            raise ValueError("spectra must be 2D (n_populations x n_frequencies)")

        if pop_idx < 0 or pop_idx >= spectra.shape[0]:
            raise IndexError(f"pop_idx {pop_idx} out of range")

        population_labels = self.get_population_labels()
        resolved_name = pop_name
        if resolved_name is None:
            if pop_idx < len(population_labels):
                resolved_name = population_labels[pop_idx]
            else:
                resolved_name = f"Population {pop_idx}"

        power = spectra[pop_idx]

        plt.figure()
        plt.plot(freqs, power, label=resolved_name)

        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power Spectrum Density")
        plt.title(f"Frequency Spectrum - {resolved_name}")
        plt.legend(title="Population")
        plt.xlim(0, max_freq_hz)

        plt.grid(True)
        plt.tight_layout()
        plt.show()


    def plot_freq_spectrum_all_populations(self, spectra, freqs, log_scale=False, max_freq_hz=100):
        """
        Plot frequency spectra for all populations with one subplot per layer/area group.
        """

        if spectra.ndim != 2:
            raise ValueError("spectra must be 2D (n_populations x n_frequencies)")

        population_labels = self.get_population_labels()
        if spectra.shape[0] != len(population_labels):
            raise ValueError(
                f"spectra rows ({spectra.shape[0]}) do not match number of populations ({len(population_labels)})"
            )

        groups = self.get_population_spectrum_groups()
        n_groups = len(groups)
        n_cols = 2
        n_rows = int(np.ceil(n_groups / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 3 * n_rows), sharex=True)
        axes = np.atleast_1d(axes).flatten()
        label_to_idx = {label: idx for idx, label in enumerate(population_labels)}

        for ax, (group_title, group_labels) in zip(axes, groups):
            for label in group_labels:
                idx = label_to_idx.get(label)
                if idx is None:
                    continue
                ax.plot(freqs, spectra[idx], label=label)

            if log_scale:
                ax.set_yscale("log")
            ax.set_title(group_title)
            ax.set_xlabel("Frequency (Hz)")
            ax.set_ylabel("Power Spectrum Density")
            ax.set_xlim(0, max_freq_hz)
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)

        for ax in axes[n_groups:]:
            ax.axis("off")

        fig.suptitle("Frequency spectra by area and layer")
        plt.tight_layout(rect=[0, 0, 1, 0.98])
        plt.show()


        





    
    
