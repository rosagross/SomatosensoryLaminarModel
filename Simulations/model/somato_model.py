import numpy as np
import os
import h5py
import yaml
import json
import sys
import time
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import mne
from parameters import Parameter

location = "mpi"
if location == "laptop":
    WDDIR = r"C:\Users\gross\OneDrive - UvA\Documents\IMPRS_Leipzig\MyProject\Modelling\ChienReplication\SomatosensoryLaminarModel"
    SIMDIR = os.path.join(WDDIR, "output")
    DATADIR = "C:\\Users\\gross\\OneDrive - UvA\\Documents\\IMPRS_Leipzig\\MyProject\\Experiment\\Analysis\\LocalCode\\data"
    RECONDIR = os.path.join(DATADIR, 'freesurfer')

if location == "mpi":
    DATADIR = os.getenv('DATADIR')
    RECONDIR = os.getenv('SUBJECTS_DIR')
    SIMDIR = os.getenv("SIMDIR")
    WDDIR = os.getenv("WDDIR")
    
figure_dir = os.path.join(SIMDIR, "Figures")

# Structured per-comparison output folders (per-run PNG + HDF5 + summary CSV).
TIMEFREQ_DIR = os.path.join(SIMDIR, "timefreq_comparison")
TIMECOURSE_DIR = os.path.join(SIMDIR, "timecourse_comparison")
PRESTIM_SPECTRUM_DIR = os.path.join(SIMDIR, "prestim_spectrum_comparison")

# Travel time from the fingertip receptor to the thalamus. The real data is
# time-locked to the fingertip pulse (t=0), but the model injects its external
# input directly at the thalamus, so the model's stimulus onset corresponds to
# real-data t = +RECEPTOR_THALAMUS_DELAY_S (≈ the cortical N20 latency). The
# simulated error window is shifted earlier by this amount so the model response
# aligns with the measured response instead of leading it.
RECEPTOR_THALAMUS_DELAY_S = 0.020  # 20 ms


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

    def __init__(self, params={}, WDDIR=None):
        
        # load in all connectivity parameters, time constants, etc.
        # these are read here (before Parameter is built) since they shape the tau matrix
        self.delay_factor      = params.get('delay_factor', 5e-3)
        self.thal_delay_factor = params.get('thal_delay_factor', 3e-3)
        self.e3b_tau           = params.get('e3b_tau', 6)
        self.e1_tau            = params.get('e1_tau', 6)
        self.e2_tau            = params.get('e2_tau', 6)
        self.p = Parameter(delay_factor=self.delay_factor, thal_delay_factor=self.thal_delay_factor,
                           e3b_tau=self.e3b_tau, e1_tau=self.e1_tau, e2_tau=self.e2_tau)
        self.tau = self.p.tau
        self.nPop = self.p.nPop
        # sigmoid function (16 x 3) --> 3 stands for parameters: r, v_thr, m_max
        self.sigm = self.p.sigmoid_params

        # parameters that will be updated from the json file 
        # (first initialized with default values) 
        self.simulation_dur = 2 # in s
        self.step_size = 0.001 # in s
        self.resolution_tstep = 0.001 # in s
        self.sfreq_saved = 1 / self.resolution_tstep 
        self.input_onset = 1.001
        # periphery→thalamus alignment delay (s); overridable via params / optimized by the GA
        self.receptor_thalamus_delay = RECEPTOR_THALAMUS_DELAY_S
        self.thal_connect = [0,0,0,0]
        self.extI_cellcounts = 1000
        self.strength_I = 0 #0.7
        self.bI_cellcounts = 100
        self.thalE_cellcounts = 500
        self.thalI_cellcounts = 500
        self.pom_cellcounts = 500
        self.sI_thal = 0.5
        self.g_thal = 2
        self.g_thalPOm = 1
        self.input_type = 'step'
        self.area = 'all' 
        self.coupling_strength = 10
        self.Ib_strength = 7
        self.Ib_noise_std = 1.0     # stationary std of OU background noise (units of Ib_strength); 0 = off
        self.Ib_noise_tau = 0.016   # OU correlation time constant (s)
        self.Ib_noise_seed = None   # seed for the isolated background-noise RNG (None = fresh draws)
        self.Iext_strength = 10
        self.Iext_duration = 0.5
        self.resistance_factor = 1

        # scaling the coupling strength between the cortical areas
        self.g_intercortical = 1

        # update parameters based on params dicts
        self.__dict__.update(params)

        # keep a copy of the params used for this run (written to the run folder)
        self.params = dict(params)

        # create input array
        Iext = self.create_Iext()
        Ib = self.create_Ibackground()
        self.gE = self.coupling_strength
        self.gI = self.coupling_strength * self.strength_I
        self.gEthal = self.g_thal
        self.gIthal = self.g_thal * self.sI_thal
        self.gPOmthal = self.g_thalPOm

        # Synaptic kernel
        self.H = np.ones((self.nPop, self.nPop+1))

        # define time steps 
        self.steps = np.arange(self.step_size, self.simulation_dur+self.step_size, self.step_size)

        # extend input arrays
        self.Iext = np.tile(Iext, (self.nPop,1))
        self.Ib = np.tile(Ib, (self.nPop,1))
        self.Ib = self.add_background_noise(self.Ib)

        self.filename = (
            f"gthal{self.g_thal}_gthalPOm{self.g_thalPOm}_sIthal{self.sI_thal}_g{self.coupling_strength}_sI{self.strength_I}_Ib{self.Ib_strength}_Ibnoise{self.Ib_noise_std}_Iextd{self.Iext_duration}_"
            f"{self.input_type}Iexts{self.Iext_strength}_Ionset{self.input_onset}_thalcells{self.thalE_cellcounts}_"
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
        self.W = self.p.get_connectivity(self.g_intercortical, self.gE, self.gI, self.gEthal, self.gIthal, self.gPOmthal, self.thal_connect, self.extI_cellcounts, self.bI_cellcounts, self.thalE_cellcounts, self.thalI_cellcounts, self.pom_cellcounts, area=self.area)

        # per-subject dipole projection vectors (forward model + labels + geometry), invariant
        # across simulation runs → built once per subjects list and reused (see compute_dipoles).
        self._dipole_projection_cache = {}


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

        # rebuild the tau matrix so tau/delay-shaping params (delay_factor,
        # thal_delay_factor, e3b_tau, e1_tau, e2_tau) take effect on every update
        self.p = Parameter(delay_factor=self.delay_factor, thal_delay_factor=self.thal_delay_factor,
                           e3b_tau=self.e3b_tau, e1_tau=self.e1_tau, e2_tau=self.e2_tau)
        self.tau = self.p.tau

        # recompute inputs and gains
        Iext = self.create_Iext()
        Ib = self.create_Ibackground()
        self.Iext = np.tile(Iext, (self.nPop, 1))
        self.Ib = np.tile(Ib, (self.nPop, 1))
        self.Ib = self.add_background_noise(self.Ib)

        self.gE = self.coupling_strength
        self.gI = self.coupling_strength * self.strength_I
        self.gEthal = self.g_thal
        self.gIthal = self.g_thal * self.sI_thal
        self.gPOmthal = self.g_thalPOm

        # update connectivity with new gains and counts
        self.W = self.p.get_connectivity(
            self.g_intercortical,
            self.gE,
            self.gI,
            self.gEthal,
            self.gIthal,
            self.gPOmthal,
            self.thal_connect,
            self.extI_cellcounts,
            self.bI_cellcounts,
            self.thalE_cellcounts,
            self.thalI_cellcounts,
            self.pom_cellcounts,
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


    def add_background_noise(self, Ib_matrix):
        """Add independent Ornstein-Uhlenbeck (colored) noise to the tiled background input.

        Each population gets its own zero-mean OU trace: an exponentially
        autocorrelated process with stationary std `Ib_noise_std` (units of
        `Ib_strength`) and correlation time `Ib_noise_tau` (s). The exact discrete
        OU update is used, so the noise statistics do not depend on `step_size`.
        A std of 0 returns the input unchanged (identical to the constant-background
        model).

        Note: with a fixed `Ib_noise_seed` every call regenerates the same noise;
        leaving `Ib_noise_seed=None` gives fresh independent noise each call (e.g.
        per parameter set when `apply_params` is called repeatedly during SBI).
        """
        if self.Ib_noise_std <= 0:
            return Ib_matrix
        rng = np.random.default_rng(self.Ib_noise_seed)
        n_pop, n_steps = Ib_matrix.shape
        a = np.exp(-self.step_size / self.Ib_noise_tau)       # decay per step
        b = self.Ib_noise_std * np.sqrt(1.0 - a**2)           # scale for stationary std
        xi = rng.standard_normal((n_pop, n_steps))
        noise = np.empty((n_pop, n_steps))
        noise[:, 0] = self.Ib_noise_std * xi[:, 0]            # start at stationary distribution
        for t in range(1, n_steps):
            noise[:, t] = a * noise[:, t - 1] + b * xi[:, t]
        return Ib_matrix + noise


    def save_to_yaml(self, filename):
        
        S = self.p.get_connectStrength()
        P = self.p.get_connectProb()
        C = self.p.get_cellcounts()
        W = self.p.get_connectivity(self.g_intercortical,self.gE, self.gI, self.gEthal, self.gIthal, self.gPOmthal, self.thal_connect, self.extI_cellcounts, self.bI_cellcounts, self.thalE_cellcounts, self.thalI_cellcounts, self.pom_cellcounts)

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

    def plot_W_heatmap(self, save_dir=None):
        """
        Plot the connectivity matrix as a heatmap.

        Args:
            save_dir (str, optional): if given, save the heatmap (PNG) and the
                connectivity matrix (CSV) into this run directory. Defaults to
                self.run_dir when it exists; if neither is available the figure
                is only drawn, not saved.
        """
        pop_names = self.get_population_labels()
        W = self.p.get_connectivity(self.g_intercortical, self.gE, self.gI, self.gEthal, self.gIthal, self.gPOmthal, self.thal_connect, self.extI_cellcounts, self.bI_cellcounts, self.thalE_cellcounts, self.thalI_cellcounts, self.pom_cellcounts, area=self.area)
        # drop the background (B) and external (Ext) input columns so the matrix is square
        W_df = pd.DataFrame(W[:, :-2], index=pop_names, columns=pop_names)

        fig, ax = plt.subplots(figsize=(14, 12))
        sns.heatmap(W_df, annot=False, cmap='coolwarm', center=0, xticklabels=pop_names, yticklabels=pop_names, ax=ax)
        ax.set_xlabel("Source population")
        ax.set_ylabel("Target population")
        ax.set_title("Connectivity matrix W")
        fig.tight_layout()

        if save_dir is None:
            save_dir = getattr(self, "run_dir", None)
        if save_dir is not None:
            # save the connectivity heatmap and matrix in the run directory
            os.makedirs(save_dir, exist_ok=True)
            fig.savefig(os.path.join(save_dir, "connectivity_heatmap.png"), dpi=300, bbox_inches="tight")
            W_df.to_csv(os.path.join(save_dir, "connectivity_matrix.csv"))
            plt.close(fig)
        return fig

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
        # with added noise
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
            "ThalPOm"
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
        rates_downsampled = self.rate[:, :: int(1000 * self.resolution_tstep)]
        rates_df = pd.DataFrame(rates_downsampled.T, columns=cells)

        # sum the potentials together and save them
        potential_sum = np.sum(self.potential, axis=1)
        potential_sum_downsampled = potential_sum[:, :: int(1000 * self.resolution_tstep)]
        potential_df = pd.DataFrame(potential_sum_downsampled.T, columns=cells)

        return rates_df, potential_df

    def prepare_run_dir(self, base_dir):
        """Create base_dir/<self.filename>/ for this run, dump params.json, return the path."""
        run_dir = os.path.join(base_dir, self.filename)
        os.makedirs(run_dir, exist_ok=True)
        def _to_jsonable(o):
            if isinstance(o, np.generic):
                return o.item()
            if isinstance(o, np.ndarray):
                return o.tolist()
            raise TypeError(f"Object of type {o.__class__.__name__} is not JSON serializable")
        with open(os.path.join(run_dir, "params.json"), "w") as f:
            json.dump(self.params, f, indent=2, default=_to_jsonable)
        self.run_dir = run_dir
        return run_dir

    def save_results_csv(self, filedir, filename, full=False, save_params=False):
        """
        Safe the simulated data in a csv file
        """
        rates_df, potential_df = self.prepare_dataframes()
        #print('saving rates', len(rates_df))

        filename = filename + ".hdf5"
        
        rates_df.to_hdf(
            os.path.join(filedir, filename), index=False, key="rates", mode="w"
        )

        potential_df.to_hdf(
            os.path.join(filedir, filename), index=False, key="summed_potential", mode="a"
        )

        if full:
            # save all potentials additionally
            psp_filename = "full_" + filename
            #print('full potential file:', psp_filename)
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
    

    def compute_full_spectrum(self):
        """
        Compute full signal frequency spectra for all populations.
        """
        analysis_params = read_analysis_params()
        sampling_params = analysis_params['sampling']

        rates_df, potentials_df = self.prepare_dataframes()
        spectra, freqs = hf.compute_spectra_full(rates_df, potentials_df, self.step_size)
        


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
        #spectra, freqs = self.compute_late_longterm_spectrum()

        #self.plot_freq_spectrum(spectra, freqs, min_freq_hz=3, max_freq_hz=50)
        #self.plot_freq_spectrum_all_populations(spectra, freqs)

        spectra_full, freqs_full = self.compute_full_spectrum()
        self.plot_freq_spectrum(spectra_full, freqs_full, min_freq_hz=3, max_freq_hz=50)
        self.plot_freq_spectrum_all_populations(spectra_full, freqs_full)
        
        if save_spectrum:
            spectrum_dir = os.path.join(SIMDIR, "spectrum_results")
            path = self.save_frequency_spectra(spectrum_dir, spectra=spectra_full, freqs=freqs_full)
            #print(f"saved spectra: {path}")


    def plot_freq_spectrum(self, spectra, freqs, pop_idx=0, pop_name=None, min_freq_hz=0, max_freq_hz=100):
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

    def load_dipole_params(self):
        # Read in preprocessing parameters
        with open(os.path.join(WDDIR, 'EEGSimulation', 'dipole_parameters.json'), 'r') as json_file:
            dipole_params = json.load(json_file)

        dipole_length = dipole_params['dipole_lengths']
        dipole_orientation = dipole_params['dipole_orientation']
        resistance_factor = self.resistance_factor

        # load cell count parameter
        cellcounts = self.p.get_cellcounts(return_A3b=True)

        return dipole_length, dipole_orientation, resistance_factor, cellcounts

    def get_population_mapping(self):
        """
        Get mapping between model populations and brain regions.
        
        Returns:
            dict: Mapping of population indices to brain regions and layers
        """
        # Population order from parameters.py:
        # A3b: E, PV, SST, VIP (indices 0-3)
        # A1: E1, PV1, SST1, VIP1, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4 (indices 4-16)  
        # S2: E1, PV1, SST1, VIP1, E2, PV2, SST2, E3, PV3, SST3, E4, PV4, SST4 (indices 17-29)
        # Thalamus: ThalE, ThalI (indices 30-31)
        
        mapping = {
            # A3b populations
            'A3b': {
                'E': 0, 'PV': 1, 'SST': 2, 'VIP': 3
            },
            # A1 populations (layers 1-4)
            'A1': {
                'L1_E': 4, 'L1_PV': 5, 'L1_SST': 6, 'L1_VIP': 7,
                'L4_E': 8, 'L4_PV': 9, 'L4_SST': 10,
                'L5_E': 11, 'L5_PV': 12, 'L5_SST': 13,
                'L6_E': 14, 'L6_PV': 15, 'L6_SST': 16
            },
            # S2 populations (layers 1-4)
            'S2': {
                'L1_E': 17, 'L1_PV': 18, 'L1_SST': 19, 'L1_VIP': 20,
                'L4_E': 21, 'L4_PV': 22, 'L4_SST': 23,
                'L5_E': 24, 'L5_PV': 25, 'L5_SST': 26,
                'L6_E': 27, 'L6_PV': 28, 'L6_SST': 29
            },
            # Thalamic populations
            'Thalamus': {
                'E': 30, 'I': 31
            }
        }
        
        return mapping

    
    def prepDipoles(self, dipole_length, dipole_orientation, resistance_factor, cellcountsE_relative):

        # dipoles
        # each dipole set has a value for each source population
        dipole_matrix = []
        for i, s in enumerate(dipole_length):
            dipole = s * dipole_orientation[i] * resistance_factor 
            dipole_matrix.append(dipole)

        dipole_array = np.array(dipole_matrix)

        # Weighted by E cell count 
        dipoles_weighted = dipole_array * cellcountsE_relative

        return dipoles_weighted


    def prepDipoles_normal(self, label, src_fixed, dipole_length, dipole_orientation, resistance_factor, cellcountsE_relative):
        """
        Differently to the prepDipoles() function, here we consider the vertex normal orientation.
        Only like this we can actually compare it to the source reconstructed values.
        """

        # get the vertex normal orientation

        # Find which source space vertices are in the label
        src_rh = src_fixed[1]  # right hemisphere
        label_verts = np.intersect1d(label.vertices, src_rh['vertno'])
        idx = np.searchsorted(src_rh["vertno"], label_verts)
        normals = src_rh["nn"][idx]
        # compute average norm from all vertices
        mean_normal = normals.mean(axis=0)
        mean_normal /= np.linalg.norm(mean_normal)
        #print(label.name)
        #print(mean_normal)
        # TODO: instead of using the mean vertex orientation compute the dipole for each single vertex. 


        # dipoles
        # each dipole set has a value for each source population
        dipole_matrix = []
        for i, s in enumerate(dipole_length):
            #dipole = s * np.dot(dipole_orientation[i], mean_normal) * resistance_factor 
            dipole = s * dipole_orientation[i] * resistance_factor 
            #print("dipole before vertex normal transformation:", s * dipole_orientation[i] * resistance_factor)
            #print(".. after vertex normal transformation:", s * np.dot((dipole_orientation)[i], mean_normal) * resistance_factor)
            mean_dipole = np.mean(dipole)
            #print('mean dipole:', mean_dipole)
            dipole_matrix.append(mean_dipole)

        dipole_array = np.array(dipole_matrix)

        # Weighted by E cell count 
        dipoles_weighted = dipole_array * cellcountsE_relative

        return dipoles_weighted


    def _build_dipole_projections(self, subjects):
        """Precompute the per-subject dipole projection vectors that are invariant across
        simulation runs.

        The forward models, subject labels, dipole geometry (dipole_parameters.json),
        resistance factor and cell counts do not depend on the optimized parameters
        (apply_params only changes inputs/gains/connectivity), so this expensive work —
        reading + converting each subject's forward solution and projecting it through
        prepDipoles_normal — is done once per subjects list and cached (see compute_dipoles).
        The Forward objects are read here and discarded; only the small projection vectors
        are kept.

        Returns:
            dict with 'exc_pops' (potential indices) and 'per_subject' (list, in `subjects`
            order, of {'A3b', 'A1', 'S2'} projection vectors).
        """
        # load dipole parameters
        dipole_length, dipole_orientation, resistance_factor, cellcounts = self.load_dipole_params()

        # get population mapping
        pop_mapping = self.get_population_mapping()

        # Extract excitatory populations (these generate the main EEG signal)
        exc_pops = []
        for area in ['A3b', 'A1', 'S2']:
            if area == 'A3b':
                exc_pops.append(pop_mapping[area]['E'])
            else:
                # For A1 and S2, include all layer E populations
                for layer in ['L1_E', 'L4_E', 'L5_E', 'L6_E']:
                    exc_pops.append(pop_mapping[area][layer])

        # define parameters and cellcounts
        dipole_lengths_A3b = dipole_length['A3b']
        dipole_orientation_A3b = dipole_orientation['A3b']
        dipole_lengths_A1 = dipole_length['A1']
        dipole_orientation_A1 = dipole_orientation['A1']
        dipole_lengths_ES2 = dipole_length['S2']
        dipole_orientation_ES2 = dipole_orientation['S2']
        cellcounts_A3b = cellcounts[:4]
        cellcounts_A3b_relative = cellcounts_A3b/np.sum(cellcounts_A3b)
        cellcounts_A1 = cellcounts[4:17]
        cellcounts_A1_relative = cellcounts_A1/np.sum(cellcounts_A1)
        cellcounts_S2 = cellcounts[17:]
        cellcounts_S2_relative = cellcounts_S2/np.sum(cellcounts_S2)
        cellcounts_EA3b_relative = cellcounts_A3b_relative[0]
        cellcounts_EA1_relative = cellcounts_A1_relative[np.array(exc_pops[1:5])-4]
        cellcounts_ES2_relative = cellcounts_S2_relative[np.array(exc_pops[-4:])-17]

        # for each subject, compute the (parameter-independent) dipole projections
        per_subject = []
        for subID in subjects:
            soma_labels, _ = self.read_labels(subID)
            label_A3b = soma_labels[7]
            label_A1 = soma_labels[4]
            label_A2 = soma_labels[9]   # RH S2 (Lat_Fis-post-rh)

            # read the per-subject elec forward model and take its fixed-orientation source
            # space, mirroring step002_inverse_solution_multisub_epochswise.py
            fwd_file = os.path.join(DATADIR, 'derivatives', 'eeg-preproc',
                                    f'sub-0{subID}', 'ses-elec',
                                    f'sub-0{subID}_ico-5_ses-elec_fwd.fif')
            fwd_vector = mne.read_forward_solution(fwd_file, verbose=0)
            fwd_fixed = mne.convert_forward_solution(
                fwd_vector, surf_ori=True, force_fixed=True, use_cps=True)
            src_fixed_sub = fwd_fixed['src']

            dipoles_A3b_sub = self.prepDipoles_normal(label_A3b, src_fixed_sub, dipole_lengths_A3b, dipole_orientation_A3b, resistance_factor, cellcounts_EA3b_relative)

            dipoles_A1_layers = []
            dipoles_ES2_layers = []

            for i, layer in enumerate(['L1_E', 'L4_E', 'L5_E', 'L6_E']):
                # compute dipoles for layers in A1
                dipole_layer_A1 = self.prepDipoles_normal(label_A1, src_fixed_sub, dipole_lengths_A1[i], dipole_orientation_A1[i], resistance_factor, cellcounts_EA1_relative[i])
                dipoles_A1_layers.append(dipole_layer_A1)

                # compute dipoles for layers in S2
                dipole_layer_S2 = self.prepDipoles_normal(label_A2, src_fixed_sub, dipole_lengths_ES2[i], dipole_orientation_ES2[i], resistance_factor, cellcounts_ES2_relative[i])
                dipoles_ES2_layers.append(dipole_layer_S2)

            per_subject.append({"A3b": dipoles_A3b_sub, "A1": dipoles_A1_layers, "S2": dipoles_ES2_layers})

        return {"exc_pops": exc_pops, "per_subject": per_subject}


    def compute_dipoles(self, subjects):
        # The forward-model projections are invariant across simulation runs, so build them
        # once per subjects list and reuse; only the dot products with self.potential (which
        # changes every simulation) are recomputed on each call.
        key = tuple(subjects)
        cache = self._dipole_projection_cache.get(key)
        if cache is None:
            cache = self._build_dipole_projections(subjects)
            self._dipole_projection_cache[key] = cache
        exc_pops = cache["exc_pops"]

        # for each subject compute dipoles for A3b, A1, S2
        simDipoles_all = []
        for proj in cache["per_subject"]:
            dipoles_A3b_sub = proj["A3b"]
            dipoles_A1_layers = proj["A1"]
            dipoles_ES2_layers = proj["S2"]

            # convolve the precomputed dipole models with the simulated data
            potentialsEA3b = self.potential[exc_pops[0], :-2]
            potentialsEA1 = self.potential[exc_pops[1:5], :-2]
            potentialsES2 = self.potential[exc_pops[-4:], :-2]

            # for each time point, compute the simulated dipole
            nE = 9
            simDipoles = np.zeros((nE, potentialsEA1.shape[2]))
            simDipoles[0] = np.dot(dipoles_A3b_sub, potentialsEA3b)

            for E in range(4):
                simDipoles[E+1] = np.dot(np.concatenate([dipoles_A1_layers[E]]), potentialsEA1[E])
                simDipoles[E+5] = np.dot(np.concatenate([dipoles_ES2_layers[E]]), potentialsES2[E])

            simDipoles_all.append(simDipoles)

        # average over subjects
        simDipoles_avg = np.mean(np.stack(simDipoles_all, axis=0), axis=0)


        return simDipoles_avg


    def plot_dipoles(self, simDipoles, raw_info):
        time = np.arange(simDipoles.shape[1]) / self.sfreq_saved
        area_groups = {
            "A3b": [0],
            "A1": [1, 2, 3, 4],
            "S2": [5, 6, 7, 8],
        }
        labels = {
            "A3b": ["E"],
            "A1": ["L1_E", "L4_E", "L5_E", "L6_E"],
            "S2": ["L1_E", "L4_E", "L5_E", "L6_E"],
        }

        fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
        for ax, area in zip(axes, ["A3b", "A1", "S2"]):
            area_indices = area_groups[area]
            for idx, label in zip(area_indices, labels[area]):
                ax.plot(time, simDipoles[idx], label=label, linewidth=1.5)
            sum_trace = simDipoles[area_indices].sum(axis=0)
            ax.plot(time, sum_trace, color="black", linewidth=3.0, label="Sum")
            ax.set_title(f"Computed Dipoles - {area}")
            ax.set_ylabel("Dipole (Am)")
            ax.grid(alpha=0.3)
            ax.legend(loc="upper right", frameon=False)

        axes[-1].set_xlabel("Time (s)")
        fig.suptitle("Computed Dipoles for E Populations", fontsize=14)
        fig.tight_layout(rect=[0, 0, 1, 0.98])

        figuredir = os.path.join('.', 'Figures')
        os.makedirs(figuredir, exist_ok=True)
        fig.savefig(os.path.join(figuredir, 'computed_dipoles_by_area.png'), dpi=300, bbox_inches='tight')
        plt.show()


    def simulate_eeg(self, raw, data_path_labels, simDipoles, fwd, src_fixed):
        """
        Simulate EEG directly from dipoles using the forward model.
        No raw data, no SourceSimulator, no event bookkeeping needed.
        """
        # --- 1. Build source estimate (stc) directly ---
        # simDipoles shape: (n_dipole_groups, n_times)
        # Combine your dipole groups as before
        combined_dipoles = (
            simDipoles[0]
            + np.sum(simDipoles[1:4], axis=0)
            + np.sum(simDipoles[5:8], axis=0)
        )  # shape: (n_times,)

        # Get vertices for the label
        label_file = os.path.join(
            data_path_labels, 'subjects', 'fsaverage', 'label', 'rh.BA3b.label'
        )
        label = mne.read_label(label_file, 'fsaverage')

        # Find which source space vertices are in the label
        src_rh = src_fixed[1]  # right hemisphere
        label_verts = np.intersect1d(label.vertices, src_rh['vertno'])

        # Build source data: broadcast the dipole signal to all label vertices
        # shape: (n_label_verts, n_times)
        n_times = combined_dipoles.shape[0]
        source_data_rh = np.outer(
            np.ones(len(label_verts)),  # all verts get the same signal
            combined_dipoles
        )

        # Construct a SourceEstimate
        # lh vertices empty, rh vertices = label verts
        stc = mne.SourceEstimate(
            data=np.vstack([
                np.zeros((len(src_fixed[0]['vertno']), n_times)),  # lh: zeros
                source_data_rh
            ]),
            vertices=[src_fixed[0]['vertno'], label_verts],
            tmin=0.0,
            tstep=1.0 / self.sfreq_saved,
            subject='fsaverage'
        )

        # --- 2. Apply forward model: stc → EEG sensor data ---
        # result shape: (n_channels, n_times)
        eeg_data = mne.apply_forward(fwd, stc, raw.info).data

        # --- 3. Wrap in EpochsArray (1 epoch, no noise) ---
        info = fwd['info'].copy()
        # Pick EEG channels from raw.info to match eeg_data
        raw.resample(self.sfreq_saved)
        eeg_picks = mne.pick_types(raw.info, meg=False, eeg=True)
        eeg_info = mne.pick_info(raw.info, eeg_picks)
        
        # EpochsArray expects shape: (n_epochs, n_channels, n_times)
        epochs = mne.EpochsArray(
            eeg_data[np.newaxis, :, :],
            info=eeg_info,
            tmin=0.0,
            baseline=None
        )
        evoked = epochs.average()

        return stc, evoked, epochs
        

    def simulate_eeg_mnesimulator(self, raw, data_path_labels, simDipoles, fwd, src_fixed):
        n_events = 1
        events = np.zeros((n_events, 3), int)
        events[:, 0] = 200 + 500 * np.arange(n_events)  # Events sample.
        events[:, 2] = 1  # All events have the sample id.

        # Area 3b
        label_file_soma_rh = os.path.join(data_path_labels, 'subjects', 'fsaverage', 'label','rh.BA3b.label')
        selected_label_soma_rh = mne.read_label(label_file_soma_rh, 'fsaverage')
        #times = np.arange(0, 10, 0.001)  # Simulate for 10 seconds at 1000 Hz
        tstep = 1.0 / self.sfreq_saved
        dipoles_downsampled = simDipoles  #[:, ::5]  # Downsample to match the time step
        #print("shape of dipoles:", dipoles_downsampled.shape)
        source_simulator = mne.simulation.SourceSimulator(src_fixed, tstep=tstep)
        source_simulator.add_data(selected_label_soma_rh, dipoles_downsampled[0], events)
        source_simulator.add_data(selected_label_soma_rh, np.sum(dipoles_downsampled[1:4], axis=0), events)
        source_simulator.add_data(selected_label_soma_rh, np.sum(dipoles_downsampled[5:8], axis=0), events)

        # TODO: this needs to be fixed in mne.make_forward_solution() 
        # it should be possible to have a non type fwd, since our raw info also 
        # is NonType
        raw.info["dev_head_t"] = fwd["info"]["dev_head_t"]
        raw_simulated = mne.simulation.simulate_raw(raw.info, source_simulator, forward=fwd)
        
        # create 
        epochs = mne.Epochs(raw_simulated, events, 1, baseline=None)#, tmin=-0.05, tmax=0.2)
        evoked = epochs.average()

        return evoked, epochs
    

    def plot_eeg(self, evoked, epochs):
        """Plot simulated EEG results."""
        fig = evoked.plot(show=False)
        figuredir = os.path.join('.', 'Figures')
        os.makedirs(figuredir, exist_ok=True)
        #fig.savefig(os.path.join(figuredir, 'simulated_evoked.pdf'))
        plt.show()

        freqs = np.arange(8, 50, 2)
        n_cycles = np.full_like(freqs, 2.0, dtype=float)
        tfr = mne.time_frequency.tfr_morlet(
            epochs,
            freqs=freqs,
            n_cycles=n_cycles,
            return_itc=False,
            average=True,
            picks="eeg",
        )
        power = tfr.data.mean(axis=0)

        tfr_fig, tfr_ax = plt.subplots(figsize=(12, 6))
        mesh = tfr_ax.pcolormesh(tfr.times, tfr.freqs, power, shading="auto", cmap="viridis")
        tfr_ax.set_title("Simulated Epochs Time-Frequency Power")
        tfr_ax.set_xlabel("Time (s)")
        tfr_ax.set_ylabel("Frequency (Hz)")
        tfr_fig.colorbar(mesh, ax=tfr_ax, label="Power ")
        tfr_fig.tight_layout()
        tfr_fig.savefig(os.path.join(figuredir, 'simulated_tfr_epochs.png'), dpi=300, bbox_inches='tight')
        plt.show()

        topo_times = np.linspace(evoked.times[0], evoked.times[-1], 3)
        topo_fig = evoked.plot_topomap(times=topo_times, show=False)
        topo_fig.savefig(os.path.join(figuredir, 'simulated_topomaps_epochs.png'), dpi=300, bbox_inches='tight')
        plt.show()



    def compute_prestim_spectrum(self, sim_dip, fmin=1.0, fmax=40.0):
        """Power spectrum of the 400 ms pre-stimulus window of the simulated ROI dipoles.

        Mirrors helper_functions.compute_freq_spectrum (detrend + Hann + |rfft|^2 / n^2).
        The window is the 400 ms ending at stimulus onset; captures ongoing oscillatory
        activity before stimulation.

        Returns:
            (freqs, dict roi -> power) restricted to [fmin, fmax] Hz.
        """
        rois = {"A3b": sim_dip[0], "A1": np.sum(sim_dip[1:5], axis=0), "S2": np.sum(sim_dip[5:9], axis=0)}
        stim_idx = round(self.input_onset / self.step_size) - 1
        n_pre = round(0.4 / self.step_size)                  # 400 ms -> 400 samples @ 1 kHz
        win = np.hanning(n_pre)
        freqs = np.fft.rfftfreq(n_pre, d=self.step_size)
        fmask = (freqs >= fmin) & (freqs <= fmax)
        spectra = {}
        for roi, sig in rois.items():
            seg = sig[stim_idx - n_pre:stim_idx]
            x = (seg - seg.mean()) * win
            spectra[roi] = ((np.abs(np.fft.rfft(x)) ** 2) / n_pre**2)[fmask]
        return freqs[fmask], spectra


    def load_target_prestim_spectrum(self, data_path):
        """Load the measured group pre-stimulus spectrum CSV written by step009.

        Returns:
            (freqs, dict roi -> power) for the electrical modality.
        """
        df = pd.read_csv(data_path)
        df = df[df["modality"] == "elec"]
        roi_map = {"BA3b": "A3b", "BA1": "A1", "S2": "S2"}
        freqs, spectra = None, {}
        for src, dst in roi_map.items():
            roi_df = df[df["roi"] == src].sort_values("freq_hz")
            if freqs is None:
                freqs = roi_df["freq_hz"].to_numpy()
            spectra[dst] = roi_df["power"].to_numpy()
        return freqs, spectra


    def compute_error_prestim_spectrum(self, data_path, sim_dip, fmin=1.0, fmax=40.0, target_dip=None,
                                       flatten_sim=False, rois=None):
        """Error between simulated and measured pre-stimulus spectra.

        Each spectrum is normalised to unit sum (relative power) over the compared
        frequencies — removing the sim/data amplitude-scale mismatch — then compared by
        MSE, averaged over ROIs. Both signals fall on a shared 5 Hz frequency grid.

        `rois` selects which ROIs the returned error averages over (default: all three).
        The returned spectra dicts always hold every ROI, so plotting/saving is unaffected.

        If `target_dip` is given (a saved dipole trace), the target spectrum is computed
        from it instead of loading `data_path` — used for parameter-recovery tests against
        a synthetic ground-truth trace (sim and target then share an identical freq grid).

        If `flatten_sim` is True, the FOOOF aperiodic (1/f) component is removed from the
        simulated spectrum before comparison — use this when the measured target loaded
        from `data_path` has already been 1/f-removed (see run_optimization.py), so both
        sides are compared on the same aperiodic-removed footing.

        Returns:
            (float mean MSE over ROIs, sim spectra dict, target spectra dict).
        """
        f_sim, spec_sim = self.compute_prestim_spectrum(sim_dip, fmin, fmax)
        f_tgt, spec_tgt = (self.compute_prestim_spectrum(target_dip, fmin, fmax)
                           if target_dip is not None else self.load_target_prestim_spectrum(data_path))
        if flatten_sim:
            from signal_preprocessing import remove_aperiodic
            for roi in spec_sim:
                f_flat, spec_sim[roi] = remove_aperiodic(f_sim, spec_sim[roi], fmin, fmax)
            f_sim = f_flat
        tmask = (f_tgt >= fmin) & (f_tgt <= fmax)
        # align on the shared frequency bins (both are multiples of 5 Hz)
        common = np.intersect1d(np.round(f_sim, 6), np.round(f_tgt[tmask], 6))
        sim_sel = np.isin(np.round(f_sim, 6), common)
        tgt_sel = np.isin(np.round(f_tgt, 6), common)
        eps, errors = 1e-10, []
        for roi in (tuple(rois) if rois else ("A3b", "A1", "S2")):
            s = spec_sim[roi][sim_sel]
            s = s / (s.sum() + eps)
            t = spec_tgt[roi][tgt_sel]
            t = t / (t.sum() + eps)
            errors.append(np.mean((s - t) ** 2))


        print('prestim error', float(np.mean(errors)))

        return float(np.mean(errors)), spec_sim, spec_tgt


    def compute_timefreq(self, simulated_dip):
        """
        Compute Morlet TF power for simulated dipoles per ROI.

        The simulated window is shifted earlier by RECEPTOR_THALAMUS_DELAY_S so the
        model's thalamic stimulus onset maps to +20 ms in the real-data frame (the
        target is time-locked to the fingertip pulse, the model to thalamic arrival).

        Returns:
            dict: roi label → (n_freqs=40, n_times=376) array, time axis -500..250 ms at 2 ms steps.
        """
        sfreq    = self.sfreq_saved                              # 1000 Hz
        stim_idx = round(self.input_onset / self.step_size) - 1  # 0-based stimulus index

        a3b_dip = simulated_dip[0]
        a1_dip  = np.sum(simulated_dip[1:5], axis=0)
        s2_dip  = np.sum(simulated_dip[5:9], axis=0)

        # Window of interest: -500 ms to +250 ms relative to stimulus onset (901 samples
        # at 1 kHz), shifted 20 ms earlier for the receptor→thalamus travel time.
        delay_samples = round(self.receptor_thalamus_delay / self.step_size)  # 20 samples @ 1 kHz (default)
        window_start = stim_idx - 500 - delay_samples
        window_end   = stim_idx + 251 - delay_samples

        tf_freqs    = np.arange(1, 41, 1).astype(float)
        tf_n_cycles = tf_freqs / 2
        tf_out = {}
        for roi_label, dip in [("A3b", a3b_dip), ("A1", a1_dip), ("S2", s2_dip)]:
            segment = dip

            power = mne.time_frequency.tfr_array_morlet(
                segment[np.newaxis, np.newaxis, :],
                sfreq=sfreq, freqs=tf_freqs,
                n_cycles=tf_n_cycles, output="power",
                decim=1, n_jobs=1,
            )[0, 0]  # (n_freqs=40, n_samples)

            # Crop to the stimulus-aligned window → (n_freqs=40, 901), downsample to 2 ms → (40, 451)
            power = power[:, window_start:window_end]
            tf_out[roi_label] = power[:,::2]

        return tf_out


    def load_target_timefreq(self, data_path):
        """
        Load raw Morlet TF power for electrical stimulation from the CSV produced by
        step009_plot_roi_epochswise_response_intensity.py.

        Returns:
            dict: model roi label → (n_freqs=40, n_times=376) array
        """
        df = pd.read_csv(data_path)
        df = df[(df["modality"] == "elec") & (df["norm_mode"] == "raw")]
        # cap the post-stimulus window at +250 ms (in-model crop; CSV spans -500..+400 ms)
        df = df[df["time_ms"] <= 250 + 1e-6]

        roi_map = {"BA3b": "A3b", "BA1": "A1", "S2": "S2"}
        tf_target = {}

        for csv_roi, model_roi in roi_map.items():
            roi_df = df[df["roi"] == csv_roi].sort_values(["freq_hz", "time_ms"])
            freqs  = np.sort(roi_df["freq_hz"].unique())
            times  = np.sort(roi_df["time_ms"].unique())
            power  = roi_df["power"].to_numpy().reshape(len(freqs), len(times))
            tf_target[model_roi] = power  # (n_freqs=40, n_times=451)

        return tf_target


    def compute_error_timefreq(self, data_path, sim_dip, target_dip=None, rois=None):
        """
        Compute the dimensionless TF error between simulation and measured data.

        Each area is normalized by its own baseline (-500 to -430 ms per frequency),
        then compared in log10 space via MSE. This makes the error scale-invariant —
        robust to the absolute amplitude difference between simulated and measured dipoles.

        If `target_dip` is given (a saved dipole trace), the target is derived from it
        via compute_timefreq instead of loading `data_path` — used for parameter-recovery
        tests against a synthetic ground-truth trace.

        `rois` selects which ROIs the returned error averages over (default: all three).
        The returned tf dicts always hold every ROI, so plotting/saving is unaffected.

        Returns:
            tuple: (float error, tf_sim dict, tf_target dict)
                - float: Mean log10-normalized MSE across `rois` (default A3b, A1, S2).
                - tf_sim: roi label → (n_freqs=40, n_times=451) simulated power array
                - tf_target: roi label → (n_freqs=40, n_times=451) measured power array
        """
        tf_sim    = self.compute_timefreq(simulated_dip=sim_dip)
        tf_target = (self.compute_timefreq(simulated_dip=target_dip)
                     if target_dip is not None else self.load_target_timefreq(data_path))

        # Baseline: -500 to -430 ms → indices 0:36 at 2 ms steps (35 * 2 = 70 ms / 2 = 35+1=36 pts)
        baseline_slice = slice(150, 200)
        # Analysis: -200 ms onward → index 150 on the 2 ms / -500..400 axis.
        # Baseline is still normalized against -500..-430 ms; only the error window is trimmed.
        analysis_slice = slice(200, 400)
        eps = 1e-10

        errors = []
        for roi in (tuple(rois) if rois else ("A3b", "A1", "S2")):
            if roi not in tf_sim or roi not in tf_target:
                raise RuntimeError(f"Missing TF data for ROI '{roi}'.")
            P_sim = tf_sim[roi]    # (n_freqs, n_times)
            P_tgt = tf_target[roi]

            if P_sim.shape != P_tgt.shape:
                # cut simulated signal to same size
                P_sim = P_sim[:, :P_tgt.shape[1]]

            # each area normalized by its own per-frequency baseline
            bl_sim = P_sim[:, baseline_slice].mean(axis=1, keepdims=True)
            bl_tgt = P_tgt[:, baseline_slice].mean(axis=1, keepdims=True)

            log_sim = np.log10(P_sim[:, analysis_slice] / (bl_sim + eps) + eps)
            log_tgt = np.log10(P_tgt[:, analysis_slice] / (bl_tgt + eps) + eps)

            errors.append(float(np.mean((log_sim - log_tgt) ** 2)))

        return float(np.mean(errors)), tf_sim, tf_target


    def plot_timefreq_comparison(self, tf_sim, tf_target):
        """
        Plot simulated vs. measured Morlet TF power side-by-side for each ROI.

        Layout: 2 rows (Simulated / Measured) × 3 cols (A3b, A1, S2).
        All panels share one colour scale (95th-percentile of simulated power across
        all areas) so the power ratios between areas are visible.
        Saves to ./Figures/tf_comparison_g-<g>_sI-<sI>_area-<area>.png.
        """
        rois      = ("A3b", "A1", "S2")
        row_labels = ("Simulated", "Measured")
        tf_freqs  = np.arange(1, 41, 1).astype(float)

        # Time axis: -500 ms to +250 ms at 2 ms steps → 376 points.
        # Display window trimmed to index 200 onward; baseline at -500..-430 kept upstream.
        start_idx = 200
        stop_idx = None
        t_ms = np.linspace(-500, 250, 376)[start_idx:stop_idx]

        figuredir = TIMEFREQ_DIR
        os.makedirs(figuredir, exist_ok=True)

        fig, axes = plt.subplots(2, 3, figsize=(13, 6), sharex=True, sharey=True)
        fig.suptitle(
            f"TF power — g={self.coupling_strength}, sI={self.strength_I}, area={self.area}",
            fontsize=11,
        )

        # one shared colour scale across all areas -> cross-area power ratios stay visible
        vmin = 0.0
        vmax = max(
            float(np.percentile(tf_sim[roi][:, start_idx:stop_idx], 95)) for roi in rois
        )

        for col, roi in enumerate(rois):
            P_sim = tf_sim[roi][:, start_idx:stop_idx]    # (n_freqs=40, n_times_trimmed)
            P_tgt = tf_target[roi][:, start_idx:stop_idx]

            for row, (data, label) in enumerate(
                zip((P_sim, P_tgt), row_labels)
            ):
                ax = axes[row, col]
                im = ax.imshow(
                    data,
                    aspect="auto",
                    origin="lower",
                    extent=[t_ms[0], t_ms[-1], tf_freqs[0], tf_freqs[-1]],
                    cmap="hot_r",
                    vmin=vmin,
                    vmax=vmax,
                )
                ax.axvline(0, color="white", lw=0.8, ls="--")
                if row == 0:
                    ax.set_title(roi)
                if col == 0:
                    ax.set_ylabel(f"{label}\nFrequency (Hz)")
                if row == 1:
                    ax.set_xlabel("Time (ms)")
                plt.colorbar(im, ax=ax, label="Power (a.u.²)", fraction=0.046, pad=0.04)

        fig.tight_layout()
        fname = (
            f"tf_comparison_g-{self.coupling_strength}"
            f"_sI-{self.strength_I}"
            f"_area-{self.area}.png"
        )
        fig.savefig(os.path.join(figuredir, fname), dpi=300, bbox_inches="tight")
        plt.show()


    def compute_timecourse(self, simulated_dip):
        """
        Extract the simulated dipole time course per ROI (no Morlet).

        Mirrors compute_timefreq's windowing/grouping. Returns the full
        -500..+250 ms window (baseline kept) so the analysis/plot window can be
        trimmed to -200 ms later while still normalizing against the
        -500..-430 ms baseline.

        Returns:
            dict: roi label -> (n_times=376,) baseline-corrected trace,
                  time axis -500..+250 ms at 2 ms steps.
        """
        stim_idx = round(self.input_onset / self.step_size) - 1  # 0-based stimulus index

        a3b_dip = simulated_dip[0]
        a1_dip  = np.sum(simulated_dip[1:5], axis=0)
        s2_dip  = np.sum(simulated_dip[5:9], axis=0)

        # Window of interest: -500 ms to +250 ms relative to stimulus onset (751 samples at 1 kHz),
        # shifted 20 ms earlier so the model's thalamic onset aligns with the real cortical
        # response (which lags the fingertip pulse by the receptor→thalamus travel time).
        delay_samples = round(self.receptor_thalamus_delay / self.step_size)  # 20 samples @ 1 kHz (default)
        win_start = stim_idx - 500 - delay_samples
        win_end   = stim_idx + 251 - delay_samples

        # Baseline window: -500 to -430 ms -> indices 0:36 on the 2 ms axis.
        baseline_slice = slice(0, 36)

        tc_out = {}
        for roi_label, dip in [("A3b", a3b_dip), ("A1", a1_dip), ("S2", s2_dip)]:
            segment = dip[win_start:win_end][::2]                  # (376,) at 2 ms steps
            segment = segment - segment[baseline_slice].mean()     # baseline-correct to -500..-430 ms
            tc_out[roi_label] = segment

        return tc_out


    def load_target_timecourse(self, data_path):
        """
        Load the measured group-average ROI time course for electrical stimulation
        from the CSV produced by step009_plot_roi_epochswise_response_intensity.py.

        Returns:
            dict: model roi label -> (n_times=376,) baseline-corrected trace.
        """
        df = pd.read_csv(data_path)
        df = df[df["modality"] == "elec"]
        # cap the post-stimulus window at +250 ms (in-model crop; CSV spans -500..+400 ms)
        df = df[df["time_s"] <= 0.250 + 1e-9]

        roi_map = {"BA3b": "A3b", "BA1": "A1", "S2": "S2"}
        baseline_slice = slice(0, 36)  # -500..-430 ms
        tc_target = {}

        for csv_roi, model_roi in roi_map.items():
            roi_df = df[df["roi"] == csv_roi].sort_values("time_s")
            trace  = roi_df["amplitude"].to_numpy()
            # Already baseline-corrected upstream; re-subtract for symmetry with the simulation.
            trace  = trace - trace[baseline_slice].mean()
            tc_target[model_roi] = trace

        return tc_target


    def compute_error_timecourse(self, data_path, sim_dip, scaling_factor=None, target_dip=None, rois=None):
        """
        Compute the dimensionless time-course error between simulation and measured data.

        Normalization depends on `scaling_factor`:
          - if given, every area's simulated trace is divided by that single scalar
            (measured traces by 1), keeping the cross-area amplitude ratios.
          - if not given, each area is normalized by *its own* peak (a separate
            scaling factor per area), computed separately for sim and measured. This
            compares the per-area shape and does NOT preserve cross-area ratios.

        If `target_dip` is given (a saved dipole trace), the target is derived from it
        via compute_timecourse instead of loading `data_path` — used for parameter-recovery
        tests against a synthetic ground-truth trace.

        `rois` selects which ROIs the returned error averages over (default: all three).
        The returned tc dicts always hold every ROI, so plotting/saving is unaffected.
        With a single ROI the `scaling_factor` distinction collapses — there are no
        cross-area ratios left to preserve, only that ROI's own amplitude scale.

        Returns:
            tuple: (float error, tc_sim dict, tc_target dict)
                - float: Mean MSE across `rois` (default A3b, A1, S2).
                - tc_sim: roi label -> (n_times=451,) simulated trace
                - tc_target: roi label -> (n_times=451,) measured trace
        """
        tc_sim    = self.compute_timecourse(simulated_dip=sim_dip)
        tc_target = (self.compute_timecourse(simulated_dip=target_dip)
                     if target_dip is not None else self.load_target_timecourse(data_path))

        # Analysis window: stimulus onset (0 ms) onward -> index 250 on the 2 ms / -500..400 axis.
        analysis_slice = slice(250, None)
        eps = 1e-10
        rois = tuple(rois) if rois else ("A3b", "A1", "S2")

        # collect the sliced, shape-matched traces per ROI
        sims, tgts = {}, {}
        for roi in rois:
            if roi not in tc_sim or roi not in tc_target:
                raise RuntimeError(f"Missing time-course data for ROI '{roi}'.")
            x_sim = tc_sim[roi][analysis_slice]
            x_tgt = tc_target[roi][analysis_slice]
            if x_sim.shape != x_tgt.shape:
                # cut simulated signal to same size
                x_sim = x_sim[:x_tgt.shape[0]]
            sims[roi], tgts[roi] = x_sim, x_tgt


        if scaling_factor:
            # one shared scalar across all areas -> preserves cross-area ratios
            sim_peak = {roi: scaling_factor for roi in rois}
            tgt_peak = {roi: 1 for roi in rois}
        else:
            # separate peak per area -> each ROI normalized to its own amplitude
            sim_peak = {roi: np.max(np.abs(sims[roi])) + eps for roi in rois}
            tgt_peak = {roi: np.max(np.abs(tgts[roi])) + eps for roi in rois}

        # the error of the ERP time course should be computed starting from stimulation onset 0 ms
        errors = [
            float(np.mean((sims[roi] / sim_peak[roi] - tgts[roi] / tgt_peak[roi]) ** 2))
            for roi in rois
        ]

        return float(np.mean(errors)), tc_sim, tc_target


    def plot_timecourse_comparison(self, tc_sim, tc_target, scaling_factor=None):
        """
        Plot simulated vs. measured ROI time courses (peak-normalized) per ROI.

        Layout: 1 row x 3 cols (A3b, A1, S2). Normalization mirrors the error
        (compute_error_timecourse): with a `scaling_factor` all areas share that one
        scalar (cross-area ratios visible); without it each area is normalized by its
        own peak over the -200..+250 ms window (per-area shape comparison).
        Saves to ./Figures/tc_comparison_g-<g>_sI-<sI>_area-<area>.png.
        """
        rois = ("A3b", "A1", "S2")
        eps  = 1e-10

        # Time axis: -500..+250 ms at 2 ms steps, trimmed to -200 ms onward.
        t_ms = np.linspace(-500, 250, 376)[150:]

        figuredir = TIMECOURSE_DIR
        os.makedirs(figuredir, exist_ok=True)

        fig, axes = plt.subplots(1, 3, figsize=(13, 3.4), sharex=True, sharey=True)
        fig.suptitle(
            f"Time course - g={self.coupling_strength}, sI={self.strength_I}, area={self.area}",
            fontsize=11,
        )


        if scaling_factor:
            # one shared scalar across all areas -> preserves cross-area ratios
            sim_peak = {roi: scaling_factor for roi in rois}
            tgt_peak = {roi: 1 for roi in rois}
        else:
            # separate peak per area -> each ROI normalized to its own amplitude
            sim_peak = {roi: np.max(np.abs(tc_sim[roi][150:])) + eps for roi in rois}
            tgt_peak = {roi: np.max(np.abs(tc_target[roi][150:])) + eps for roi in rois}

        for ax, roi in zip(axes, rois):
            x_sim = tc_sim[roi][150:] / sim_peak[roi]
            x_tgt = tc_target[roi][150:] / tgt_peak[roi]

            ax.plot(t_ms, x_sim, color="C1", lw=1.5, label="Simulated")
            ax.plot(t_ms, x_tgt, color="black", lw=1.5, label="Measured")
            ax.axvline(0, color="grey", lw=0.8, ls="--")
            ax.axhline(0, color="k", lw=0.5, alpha=0.3)
            ax.set_title(roi)
            ax.set_xlabel("Time (ms)")
            if ax is axes[0]:
                ax.set_ylabel("Amplitude (shared-scalar normalized)" if scaling_factor
                              else "Amplitude (per-area peak normalized)")
            ax.legend(frameon=False, fontsize=8)

        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fname = (
            f"tc_comparison_g-{self.coupling_strength}"
            f"_sI-{self.strength_I}"
            f"_area-{self.area}.png"
        )
        fig.savefig(os.path.join(figuredir, fname), dpi=300, bbox_inches="tight")
        plt.show()


    def plot_prestim_spectrum_comparison(self, data_path, sim_dip, fmin=1.0, fmax=40.0, target_dip=None):
        """Plot simulated vs measured pre-stim power spectra (unit-sum relative power) per ROI.

        Layout 1x3 (A3b, A1, S2) on the shared 5 Hz grid — the same quantity the error uses.
        Self-contained (recomputes the sim spectrum and reloads the target) so the target's own
        frequency grid is available for alignment. When `target_dip` is given, the target
        spectrum comes from that saved trace instead of `data_path` (parameter-recovery mode).
        Saves to PRESTIM_SPECTRUM_DIR/prestim_spectrum_comparison_g-<g>_sI-<sI>_area-<area>.png.
        """
        rois, eps = ("A3b", "A1", "S2"), 1e-10
        f_sim, spec_sim = self.compute_prestim_spectrum(sim_dip, fmin, fmax)
        f_tgt, spec_tgt = (self.compute_prestim_spectrum(target_dip, fmin, fmax)
                           if target_dip is not None else self.load_target_prestim_spectrum(data_path))
        tmask  = (f_tgt >= fmin) & (f_tgt <= fmax)
        f_tgt  = f_tgt[tmask]
        common = np.intersect1d(np.round(f_sim, 6), np.round(f_tgt, 6))
        sim_sel = np.isin(np.round(f_sim, 6), common)
        tgt_sel = np.isin(np.round(f_tgt, 6), common)

        figuredir = PRESTIM_SPECTRUM_DIR
        os.makedirs(figuredir, exist_ok=True)
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.4), sharex=True, sharey=True)
        fig.suptitle(
            f"Pre-stim spectrum - g={self.coupling_strength}, sI={self.strength_I}, area={self.area}",
            fontsize=11,
        )
        for ax, roi in zip(axes, rois):
            s = spec_sim[roi][sim_sel];        s = s / (s.sum() + eps)
            t = spec_tgt[roi][tmask][tgt_sel]; t = t / (t.sum() + eps)
            ax.plot(common, s, color="C1", lw=1.5, marker="o", ms=3, label="Simulated")
            ax.plot(common, t, color="black", lw=1.5, marker="o", ms=3, label="Measured")
            ax.set_title(roi)
            ax.set_xlabel("Frequency (Hz)")
            if ax is axes[0]:
                ax.set_ylabel("Relative power (unit-sum normalized)")
            ax.legend(frameon=False, fontsize=8)

        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fname = (
            f"prestim_spectrum_comparison_g-{self.coupling_strength}"
            f"_sI-{self.strength_I}"
            f"_area-{self.area}.png"
        )
        fig.savefig(os.path.join(figuredir, fname), dpi=300, bbox_inches="tight")
        plt.show()


    def _comparison_param_attrs(self):
        """Parameter metadata stored on every saved comparison file (mirrors save_frequency_spectra)."""
        return {
            "g_intercortical": self.g_intercortical,
            "coupling_strength": self.coupling_strength,
            "strength_I": self.strength_I,
            "Ib_strength": self.Ib_strength,
            "Iext_strength": self.Iext_strength,
            "Iext_duration": self.Iext_duration,
            "step_size": self.step_size,
            "input_onset": self.input_onset,
            "area": self.area,
            "input_type": self.input_type,
        }

    def _comparison_stem(self, prefix):
        """Per-run filename stem encoding the full parameter set (uniqueness across the sweep)."""
        return (
            f"{prefix}_ginter{self.g_intercortical}_g{self.coupling_strength}_sI{self.strength_I}"
            f"_Ib{self.Ib_strength}_Iextd{self.Iext_duration}_{self.input_type}Iexts{self.Iext_strength}"
            f"_area{self.area}"
        )

    def append_comparison_summary(self, tf_error=None, tc_error=None):
        """Append one summary row (params + scalar errors) to SIMDIR/comparison_summary.csv (call once per run)."""
        summary_path = os.path.join(SIMDIR, "comparison_summary.csv")
        row = self._comparison_param_attrs()
        row["tf_error"] = tf_error
        row["tc_error"] = tc_error
        pd.DataFrame([row]).to_csv(
            summary_path,
            mode="a",
            header=not os.path.exists(summary_path),
            index=False,
        )


    def save_timefreq_comparison(self, filedir, tf_sim, tf_target, tf_error, filename=None):
        """
        Save the simulated/measured TF maps + scalar error for one run to HDF5.

        Layout: groups 'sim' and 'target', each with datasets A3b/A1/S2 (40 x 376);
        plus 'freqs' (1-40 Hz) and 'times_ms' (-500..+250 ms). Parameter metadata is
        stored as attrs so the animation can filter on them.
        """
        os.makedirs(filedir, exist_ok=True)
        if filename is None:
            filename = self._comparison_stem("tf_comparison") + ".hdf5"
        elif not filename.endswith(".hdf5"):
            filename = filename + ".hdf5"

        freqs = np.arange(1, 41, 1).astype(float)
        # Derive the time axis from the actual map width so the stored axis always
        # matches the data (the TF window spans -500..+400 ms by design).
        n_times  = np.asarray(tf_sim["A3b"]).shape[1]
        times_ms = np.linspace(-500, 400, n_times)

        filepath = os.path.join(filedir, filename)
        with h5py.File(filepath, "w") as h5f:
            sim_grp = h5f.create_group("sim")
            tgt_grp = h5f.create_group("target")
            for roi in ("A3b", "A1", "S2"):
                sim_grp.create_dataset(roi, data=np.asarray(tf_sim[roi]))
                tgt_grp.create_dataset(roi, data=np.asarray(tf_target[roi]))
            h5f.create_dataset("freqs", data=freqs)
            h5f.create_dataset("times_ms", data=times_ms)
            h5f.attrs["tf_error"] = float(tf_error)
            for key, value in self._comparison_param_attrs().items():
                h5f.attrs[key] = value

        return filepath

    def save_timecourse_comparison(self, filedir, tc_sim, tc_target, tc_error, filename=None):
        """
        Save the simulated/measured TC traces + scalar error for one run to HDF5.

        Layout: groups 'sim' and 'target', each with datasets A3b/A1/S2 (376,);
        plus 'times_ms' (-500..+250 ms). Parameter metadata stored as attrs.
        """
        os.makedirs(filedir, exist_ok=True)
        if filename is None:
            filename = self._comparison_stem("tc_comparison") + ".hdf5"
        elif not filename.endswith(".hdf5"):
            filename = filename + ".hdf5"

        # Derive the time axis from the trace length (TC window spans -500..+250 ms).
        n_times  = np.asarray(tc_sim["A3b"]).shape[0]
        times_ms = np.linspace(-500, 250, n_times)

        filepath = os.path.join(filedir, filename)
        with h5py.File(filepath, "w") as h5f:
            sim_grp = h5f.create_group("sim")
            tgt_grp = h5f.create_group("target")
            for roi in ("A3b", "A1", "S2"):
                sim_grp.create_dataset(roi, data=np.asarray(tc_sim[roi]))
                tgt_grp.create_dataset(roi, data=np.asarray(tc_target[roi]))
            h5f.create_dataset("times_ms", data=times_ms)
            h5f.attrs["tc_error"] = float(tc_error)
            for key, value in self._comparison_param_attrs().items():
                h5f.attrs[key] = value

        return filepath


    def save_dipole_trace(self, sim_dip, filedir, filename=None):
        """
        Save the full simulated dipole trace + its pre-stimulus spectrum to HDF5.

        Used to build a synthetic ground-truth target for parameter-recovery testing
        (see run_optimization.py): the saved trace can be reloaded and fed to the same
        error functions as the measured data. The full (9, n_times) trace is stored so
        the tf/tc/ps errors can all be computed against it later.

        Layout: dataset 'dipole' (9 x n_times); group 'prestim_spectrum' with 'freqs'
        plus A3b/A1/S2 power datasets. Parameter metadata + sampling info stored as attrs.
        """
        os.makedirs(filedir, exist_ok=True)
        if filename is None:
            filename = self._comparison_stem("dipole_trace") + ".hdf5"
        elif not filename.endswith(".hdf5"):
            filename = filename + ".hdf5"

        freqs, spectra = self.compute_prestim_spectrum(sim_dip)

        filepath = os.path.join(filedir, filename)
        with h5py.File(filepath, "w") as h5f:
            h5f.create_dataset("dipole", data=np.asarray(sim_dip))
            ps_grp = h5f.create_group("prestim_spectrum")
            ps_grp.create_dataset("freqs", data=np.asarray(freqs))
            for roi in ("A3b", "A1", "S2"):
                ps_grp.create_dataset(roi, data=np.asarray(spectra[roi]))
            for key, value in self._comparison_param_attrs().items():
                h5f.attrs[key] = value
            h5f.attrs["sfreq_saved"] = self.sfreq_saved

        return filepath

    def load_dipole_trace(self, filepath):
        """
        Load a dipole trace saved by save_dipole_trace.

        Returns the (9, n_times) dipole array, ready to plug into the error functions
        in place of a freshly computed sim_dip (as the `target_dip` argument).
        """
        with h5py.File(filepath, "r") as h5f:
            return h5f["dipole"][()]


    def read_labels(self, subID):

        # read labels
        somatosensory_label_files = [
            os.path.join(RECONDIR, f'sub-0{subID}', 'label', 'lh.BA1_exvivo.label'),
            os.path.join(RECONDIR, f'sub-0{subID}', 'label', 'lh.BA2_exvivo.label'),
            os.path.join(RECONDIR, f'sub-0{subID}', 'label', 'lh.BA3a_exvivo.label'),
            os.path.join(RECONDIR, f'sub-0{subID}', 'label', 'lh.BA3b_exvivo.label'),
            os.path.join(RECONDIR, f'sub-0{subID}', 'label', 'rh.BA1_exvivo.label'),
            os.path.join(RECONDIR, f'sub-0{subID}', 'label', 'rh.BA2_exvivo.label'),
            os.path.join(RECONDIR, f'sub-0{subID}', 'label', 'rh.BA3a_exvivo.label'),
            os.path.join(RECONDIR, f'sub-0{subID}', 'label', 'rh.BA3b_exvivo.label'),
        ]


        # Read labels
        somatosensory_labels = [mne.read_label(label_file) for label_file in somatosensory_label_files]
        somatosensory_label_names = [
            'LH BA1',
            'LH BA2',
            'LH BA3a',
            'LH 3b',
            'RH BA1',
            'RH BA2',
            'RH BA3a',
            'RH 3b'
        ]

        # S2 label
        sub_name = f"sub-0{subID}"
        labels = mne.read_labels_from_annot(sub_name, parc='aparc.a2009s', subjects_dir=RECONDIR, verbose=0)
        # Relevant labels:
        # 'S_circular_insula_sup' - superior circular sulcus of insula
        # 'G_temp_sup-Plan_polar' 
        # 'Lat_Fis-post' - posterior lateral fissure (parietal operculum = S2)
        # 'G_and_S_transv_frontopol'
        s2_lh_label = next(label for label in labels if label.name == 'Lat_Fis-post-lh')
        s2_rh_label = next(label for label in labels if label.name == 'Lat_Fis-post-rh')

        # append S2 to the label list (new indices 8 = LH S2, 9 = RH S2)
        somatosensory_labels += [s2_lh_label, s2_rh_label]
        somatosensory_label_names += ['LH S2', 'RH S2']

        return somatosensory_labels, somatosensory_label_names














        
    
