'''
Plots:
1) Effect of input intensity and duration on firing rates and potentials
    1.1) Line plot, single population - plots trajectory
    1.2) Heatmap, single population - difference between longterm and baseline 
    1.3) Heatmap, single population - plots the different dynamic functions (memory, transfer, ..) 
        - x axis: input duration  
        - y axis: input strength
        - color map: dynamic functions
    1.4) Heatmap, multiple populations - longterm versus baseline by G, input duration and strength
2) Baseline activity and coupling strength 
    2.1) Line plot, multiple populations - axis x:G, y: rate in (Hz)
3) Background activity in steady state
    3.1) Average longterm activity
        - x axis: input strength
        - y axis: Steady state potential of Layer 5 population 

4) Effect of connection strength from and to PV interneurons
'''

# %% 
import numpy as np
import datetime
import json
import os
from matplotlib.ticker import FormatStrFormatter, FuncFormatter, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm  # Import the colormap module
from matplotlib.colors import ListedColormap, Normalize, BoundaryNorm
from scipy.signal import find_peaks
import seaborn as sns
import pandas as pd
from plotting_style import figure_style
from helper_functions import *
from plotting_functions_analysis import *

colors, _ = figure_style() 

# define directories of stored data
SIMDIR = "/data/pt_02989"
WDDIR = os.getenv("WDDIR")

raw_dir = os.path.join(SIMDIR, "output_grossmannr")
# processed characteristics now live inside each per-run folder under raw_dir
processed_dir = raw_dir
# read params
params = load_parameters(WDDIR)


# some global settings
sampling_params = params['sampling']
step_size = sampling_params['step_size']
sample_delay_immediate = sampling_params['sample_delay_immediate']
sample_delay_late = sampling_params['sample_delay_late']  # when to start the long term behaviour "window"
sample_dur = sampling_params['sample_dur']
offset = sampling_params['offset']
input_onset = params['input_onset']
input_type = 'step'
thalamus_source = 'Jiang'
# Fixed tokens of the runs analysed here. All of these are part of the run-folder name
# (see SomatoModel.filename), so every load_* call has to pass them or it builds a stem
# that matches nothing. Change here to look at a different slice of the sweep; cells that
# sweep g_thalPOm or Im themselves pass their own value.
# They come from analysis_parameter.json, which must describe the batch in filedir - if a
# load_* call raises FileNotFoundError, check that file against the folder names on disk
# before touching anything else.
g_thalPOm = params['g_thalPOms'][0]
Ib_noise_std = params['Ib_noise_std']
g_thal = params['g_thal']
sI_thal = params['sI_thal']
Im_strength = params['Im_strengths'][0]
mI_cellcounts = params['mI_cellcounts']

# Representative operating point for the cells below that only need "one run from this
# batch" rather than a specific one. Derived from the swept grid instead of hardcoded, so
# they keep working when the sweep moves - hardcoded literals here are what silently broke
# every one of these cells when the batch changed. A cell that needs a *different* point
# should still pick it out of params rather than writing a number.
g_example = params['coupling_strengths'][len(params['coupling_strengths']) // 2]
sI_example = params['strength_I'][0]
Ib_example = params['backgrndI_strengths'][0]
g_inter_example = params['g_inter']
# connectivity params
params = load_parameters(WDDIR)
extI_cellcounts = params['extI_cellcounts']
bI_cellcounts = params['bI_cellcounts']
thal_cellcounts = params['thal_cellcounts']

# load connected and unconnected network
suffix = ''

figure_dir = os.path.join(SIMDIR, "Figures", "global_dynamics", suffix)
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)


# %% Make plots that demonstrate the sampling time line 

# load trajectory to plot
g = 1.0
Iext_dur = 0.005
Iext_str = 0
Ib_str = 40
sI = 0.6
g_inter = 0.4

rates_df, potentials_df, filename = load_simulation_data(g, g_inter, sI, Ib_str, Iext_dur, Iext_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, raw_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal) 
trajectory_df = load_trajectory(rates_df, potentials_df, g, sI, Iext_dur, Iext_str, Ib_str, step_size)

# choose population
population = 'E3S2'
line_df = trajectory_df[trajectory_df['population']==population]

# LONG TERM and DURING INPUT
# Add a red vertical line at time of input offset (start of sampling) and stop of sampling
simulation_time = 4*1e-3
start_sample = input_onset + Iext_dur + sample_delay_immediate
stop_sample = start_sample + sample_dur
input_offset = input_onset + Iext_dur
offset = 0.1 # time between baseline sampling and start of input 
baseline_start = input_onset - (sample_dur+offset)
baseline_stop = baseline_start + sample_dur
plotting_cut = 10

# plot the input
steps = np.arange(step_size, simulation_time+step_size, step_size)
input_line = np.zeros(len(line_df))
input_line[int((input_onset)/step_size):int(input_offset/step_size)] = Iext_str
fig = plt.figure(figsize=(10,5))
#plt.plot(steps[:-plotting_cut], input_line[:-plotting_cut], color='green', linewidth=2)
sns.lineplot(data=line_df[:-plotting_cut], x='time', y='rate', hue='InputStrength', legend='', palette=['grey'])
#sns.lineplot(data=line_df[:-plotting_cut], x='time', y='potential', hue='InputStrength', legend='', palette=['grey'])
plt.axvline(x=baseline_start, color='purple', linestyle='--', linewidth=1, label='Baseline Sample')
plt.axvline(x=baseline_stop, color='purple', linestyle='--', linewidth=1, label='')
plt.axvspan(baseline_start, baseline_stop, alpha=0.2, color='purple')
plt.axvline(x=start_sample, color='red', linestyle='--', linewidth=1, label='Long Term Sample')
plt.axvline(x=stop_sample, color='red', linestyle='--', linewidth=1, label='')
plt.axvspan(start_sample, stop_sample, alpha=0.2, color='red')
plt.axvline(x=input_onset, color='blue', linestyle='--', linewidth=1, label='During Input Sample')
plt.axvline(x=input_offset, color='blue', linestyle='--', linewidth=1, label='')
plt.axvspan(input_onset, input_offset, alpha=0.2, color='blue')
plt.ylabel('Rate (Hz)')
plt.xlabel('Time (sec)')
plt.legend()
sns.despine(trim=True)

#plt.savefig('C:/Users/gross/OneDrive - UvA/Documents/IMPRS_Leipzig/IMPRS SummerSchool/Poster/plotting_windows.pdf')
plt.savefig(os.path.join(figure_dir, 'plotting_windows.pdf'), bbox_inches='tight')
plt.show()

# %%

'''
1.1) SINGLE PLOT: Plot example trajectory
    - plot style: line plot
    - y axis: rate
    - x axis: time
'''

# choose example settings: the grid point closest to the A1 best fit (g=3.293, sI=0.595)
g = 1.0
sI = 0.6
Iext_dur = 0.005 #params['input_durations'][0]
Iext_str = 0 #params['input_strengths'][0]
Ib_str = 40 #params['backgrndI_strengths'][0]
g_inter = 0.4 #params['g_inter']

rates_df, potentials_df, filename = load_simulation_data(g, g_inter, sI, Ib_str, Iext_dur, Iext_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, raw_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal) 
trajectory_df = load_trajectory(rates_df, potentials_df, g, sI, Iext_dur, Iext_str, Ib_str, step_size)

population = 'E3' # 'E1'
line_df = trajectory_df[trajectory_df['global_coupling']==g]
line_df = line_df[line_df['population']==population]
plotting_window = []

# plot the input
steps = np.arange(step_size, simulation_time+step_size, step_size)
fig = plt.figure(figsize=(10,5))
sns.lineplot(data=line_df, x='time', y='rate', hue='InputStrength', legend='', palette=['grey'])
sns.despine(trim=True)
plt.ylabel('Rate (Hz)')
figure_name = f'trajectory_g{g}sI{sI}_pop{population}_Iduration{Iext_dur}_{Iext_str}.pdf'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.show()


# %%
'''
1.1b) RESPONSE-TYPE EXAMPLES (publication / poster figure)
Show how one population (E3S2, S2 L5) reacts to the SAME brief stimulus under four
different parameter settings, one stacked panel each:
    - non-responsive : barely reacts, returns to baseline
    - transfer       : transient response, returns to baseline
    - memory         : switches to a persistent active state
    - oscillation onset / continued : post-stimulus rhythm that either appears with the
      pulse or was already running in the baseline
Only g / sI differ between panels; the stimulus (duration/strength/background) is fixed.

CAVEAT for the current batch: the four labels above describe regimes this batch does not
actually contain. Of its 253 stimulated runs (Iexts=50, Im=0), all 253 classify as
'transfer' - every one responds (|duringInputVSbaseline_rate| > 0.1) and not one retains
(max |longtermVSbaseline_rate| = 0.025, an order of magnitude below the 0.1 gate). So there
is no non-responsive and no memory example to show. Oscillation is present (157 runs ring
already in the baseline) but is not a separate regime here: response strength and baseline
oscillation correlate at r = 0.98, i.e. they are one axis, driven by low g.

The four (g, sI) pairs below are therefore quartiles along that single response/oscillation
axis, read off this batch's processed characteristics, and are labelled for what they show
rather than for the four regimes the figure was designed around. Restore the original
four-regime version once a batch spans them - a wider g/sI range, or a stronger stimulus.

The population is E3S2 rather than E3 because in this batch the S1 populations never
oscillate - every rhythm lives in S2.
'''

dark2 = sns.color_palette('Dark2')
# quartiles of duringInputVSbaseline_rate for E3S2 across the stimulated runs; the
# diffRate_baseline of each is quoted because it tracks the response almost exactly
response_examples = [
    {'label': 'Weak transfer',      'g': 1.58, 'sI': 1.1, 'color': dark2[1]},   # osc 0.07
    {'label': 'Moderate transfer',  'g': 0.98, 'sI': 1.0, 'color': dark2[2]},   # osc 0.10
    {'label': 'Strong transfer',    'g': 1.46, 'sI': 0.7, 'color': dark2[3]},   # osc 0.12
    {'label': 'Transfer + ringing', 'g': 0.84, 'sI': 0.6, 'color': [0.5, 0.5, 0.5]},  # osc 0.31
]

fixed_params = {
    'g_inter': params['g_inter'],
    'Ib_str': params['backgrndI_strengths'][0],
    'Iext_dur': params['input_durations'][0],
    'Iext_str': params['input_strengths'][-1],   # the stimulated runs, not the 0-strength control
    'input_onset': input_onset,
    'thal_cellcounts': thal_cellcounts,
    'bI_cellcounts': bI_cellcounts,
    'extI_cellcounts': extI_cellcounts,
    'input_type': input_type,
    'step_size': step_size,
    'g_thalPOm': g_thalPOm,
    'Ib_noise_std': Ib_noise_std,
    'g_thal': g_thal,
    'sI_thal': sI_thal,
}

plot_response_type_examples(response_examples, fixed_params, 'E3S2', raw_dir, figure_dir, suffix=suffix)


# %%
'''
1.1c) POm-GAIN SENSITIVITY (publication figure)
Show how the E3S2 summed potential responds to the SAME stimulus as the POm output gain
g_thalPOm is varied. Near g_thalPOm ~ 1 the post-stimulus attractor is highly sensitive, so a
few hand-picked values each land in a different regime (oscillation / elevated / giant transient
/ quenched / extreme high). All other parameters are fixed at the processed operating point.
'''

# the g_thalPOm values simulated in this batch
gthalPOm_values = params['g_thalPOms']

gthalPOm_fixed = {
    'g': 3.25,
    'g_inter': params['g_inter'],
    'sI': 0.6,
    'Ib_str': params['backgrndI_strengths'][0],
    'Iext_dur': params['input_durations'][0],
    'Iext_str': params['input_strengths'][0],
    'input_onset': input_onset,
    'thal_cellcounts': thal_cellcounts,
    'bI_cellcounts': bI_cellcounts,
    'extI_cellcounts': extI_cellcounts,
    'input_type': input_type,
    'step_size': step_size,
    'g_thalPOm': g_thalPOm,
    'Ib_noise_std': Ib_noise_std,
    'g_thal': g_thal,
    'sI_thal': sI_thal,
}

plot_gthalPOm_potential_sweep(gthalPOm_values, gthalPOm_fixed, [['E2', 'E2S2'], ['E3', 'E3S2']], raw_dir, figure_dir, suffix=suffix, plot_window=(-0.1,0.3))


# %%
'''
1.1d) POm-GAIN SENSITIVITY - summed S2 dipole (publication figure)
Companion to 1.1c: the summed computed dipole of area S2 for each g_thalPOm value. The dipole
needs the full source-resolved potential (SomatoModel.compute_dipoles), which is not saved, so
each parameter set is re-simulated here and projected through a subject forward model.
Requires DATADIR / SUBJECTS_DIR env vars (the model reads them at import).
'''
gthalPOm_values = params['g_thalPOms']

# model-side overrides (model parameter names) for the same operating point as 1.1c
gthalPOm_sim = {
    'coupling_strength': 3.25,
    'strength_I': 0.6,
    'Iext_duration': params['input_durations'][0],
    'Iext_strength': params['input_strengths'][0],
    'Ib_strength': params['backgrndI_strengths'][0],
    'g_intercortical': params['g_inter'],
    'g_thal': 2,
    'sI_thal': 0.5,
    'area': 'all',
    'resistance_factor': 1,
    'delay_factor': 0.005,
    'extI_cellcounts': 1000,
    'bI_cellcounts': 100,
}

plot_gthalPOm_dipole_sweep(gthalPOm_values, gthalPOm_sim, figure_dir, subjects=[15], area='S2', suffix=suffix, plot_window=(-0.1,0.3))


# %%
"""
Plot a heatmap showing the effect of Input Strength versus Input Duration
"""

# choose a coupling strength, background input strength and a population
g = g_example
sI = sI_example
Ib_str = Ib_example
g_inters = [g_inter_example]
Iext_str = params['input_strengths']
Iext_dur = params['input_durations']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
population = 'E3b'

data_df = data_df[data_df['globalCoupling']==g]
data_df = data_df[data_df['strength_I']==sI]
data_df = data_df[data_df['population']==population]
data_df = data_df[data_df['BckgndInputStrength']==Ib_str]

data_heatmap = data_df.pivot(index='InputStrength',columns='InputDuration', values='longtermVSbaseline_rate')
sns.heatmap(data_heatmap, cmap='magma')

# %%
'''
1.2) SINGLE PLOT: Effect of input intensity and duration on potentials/rates in the steady state in comparison to the baseline
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    - measure: longtermVSbaseline potential or rate
'''
values = 'longtermVSbaseline_rate'
g_inters = [g_inter_example]

# choose a coupling strength and a population
g = g_example
sI = sI_example
Iext_str = params['input_strengths']
Iext_dur = params['input_durations']
Ib_str = Ib_example
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)

population = 'E3b'
data_df = data_df[data_df['globalCoupling']==g]
data_df = data_df[data_df['strength_I']==sI]
data_df = data_df[data_df['population']==population]
data_df = data_df[data_df['BckgndInputStrength']==Ib_str]
data_heatmap = data_df.pivot(index='InputStrength',columns='InputDuration', values=values)
sns.heatmap(data_heatmap, cmap='magma')

# %% 

# 1.3) MULTI fingerprint PLOT:
print('MULTI fingerprint plot')
#populations = np.array(['E1', 'E2', 'E3', 'E4','P1', 'P2', 'P3', 'P4','S1', 'S2', 'S3', 'S4', 'V1']) 
# look at several difference E-I balance values
sIs = params['strength_I'][-10:-3]
# choose a global coupling strength and a population
g = params['coupling_strengths'][5:15]
Ib_strs = [Ib_example]
g_inters = [g_inter_example]
Iext_str = params['input_strengths']
Iext_dur = params['input_durations']
populations = ['E3b']
for ib in Ib_strs:
    for p in populations:
        data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sIs, ib, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
        sImulti_fingerprint_IextDurVsStr(data_df, g, sIs, Ib_str, p, thalamus_source, figure_dir)

# %%

# 1.3a) Fingerprint: sI vs Input Strength (fixed g)
g = g_example
# scalars, not lists: fingerprint_sI_vs_IextStr filters InputDuration == Iext_dur and
# BckgndInputStrength == Ib_str, and a list makes pandas compare elementwise and raise
Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths']
Ib_str = Ib_example
sIs = params['strength_I']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sIs, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
population = 'E3b'
fingerprint_sI_vs_IextStr(data_df, g, Iext_dur, Ib_str, population, thalamus_source, figure_dir)

# %%
# 1.3b) Fingerprint: g vs Input Strength (fixed sI)
sI = sI_example
Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths']
Ib_str = params['backgrndI_strengths'][0]
g = params['coupling_strengths']
g_inters = [g_inter_example]
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
population = 'E3b'
fingerprint_g_vs_IextStr(data_df, sI, Iext_dur, Ib_str, population, thalamus_source, figure_dir)


# %%
'''
1.4) MULTI PLOT: Effect of input intensity and duration on firing rates
    - plot style: heatmap
    - y axis: intensity
    - x axis: duration
    - subplot columns: populations
    - subplot rows: coupling strengths
'''

g = params['coupling_strengths']
rate_measure = 'diffRate_lateLongterm'
sI = sI_example
Ib_str = Ib_example
Iext_dur = params['input_durations']
Iext_str = params['input_strengths']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)

populations = np.array(['E1', 'E2', 'E3', 'E4']) 
#populations = np.array(['P1', 'P2', 'P3', 'P4']) 
#populations = np.array(['S1', 'S2', 'S3', 'S4', 'V1']) 

multiPop_heatmap_IextDurVsStr(data_df, g, sI, Ib_str, populations, rate_measure, figure_dir)

# %% 
"""
2.1) Effect of Coupling Strengths on Longterm/steady state
Plot difference between Mininum and Maximum Firing rates 
"""

# choose settings. This batch has no stimulus-free runs, so the single stimulus it does have
# is used; the late-longterm window sampled here still sits well after the input offset.
g = params['coupling_strengths']
g_inters = [params['g_inter']]
Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths'][0]
sI = 0.6
Ib_str = params['backgrndI_strengths'][0]
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
multiLayer_couplingOnLongeterm_diffRate(data_df, Iext_dur, Iext_str, Ib_str, sI, thalamus_source, figure_dir)



# %%
"""
4) Oscillation frequency analysis
"""

# 4.1) Heatmap of dominant frequency (late longterm)
g = g_example
sI = sI_example
Iext_str = params['input_strengths']
Iext_dur = params['input_durations']
Ib_str = Ib_example
g_inters = [g_inter_example]
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
population = 'E1S2'
heatmap_frequency_IextDurVsStr(data_df, g, sI, Ib_str, population, "lateLongterm", "Rate", Iext_dur, Iext_str, thalamus_source, figure_dir)

# %%

# 4.2) Coupling strength vs frequency (late longterm)
# S2 populations: in this batch the S1 side never crosses the oscillation gate.
Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths'][0]
sI = 0.6
population = 'E1S2'
Ib_str = params['backgrndI_strengths'][0]
g_inters = [params['g_inter']]
g = params['coupling_strengths']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
multiLayer_couplingOnFrequency(data_df, Iext_dur, Iext_str, Ib_str, sI, "baseline", "Potential", thalamus_source, figure_dir)

#%%
# 4.3) Heatmap: frequency vs coupling strength and sI (single population)
Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths'][0]
Ib_strs = [params['backgrndI_strengths'][0]]
g_inters = [params['g_inter']]
vmax=40

areas = ['A3b', 'A1', 'S2']
g = params['coupling_strengths']
sI = params['strength_I']
for a in areas:
    for ib in Ib_strs:
        data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, ib, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
        heatmap_frequency_coupling_vs_sI(data_df, a, "baseline", "Potential", Iext_dur, Iext_str, Ib_str, figure_dir, vmin=0, vmax=vmax)

#%%
# 4.3b) Heatmap: frequency vs coupling strength and sI (single chosen population)
Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths'][0]
Ib_str = params['backgrndI_strengths'][0]
population = 'E3S2'
g = params['coupling_strengths']
sI = params['strength_I']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
heatmap_frequency_couplingVsSI_singlePop(data_df, population, "baseline", "Potential", Iext_dur, Iext_str, Ib_str, figure_dir, vmin=0, vmax=vmax)

# %%

# 4.4) Frequency vs oscillation amplitude scatter (late longterm)
population = 'E3S2'  # S2: the S1 populations never cross the oscillation gate in this batch
scatter_frequency_vs_diff(data_df, "lateLongterm", "Potential", population, figure_dir)

# %%
# 4.5) Background input strength vs peak-frequency power (E populations, 3 ROIs)
g = g_example
sI = sI_example
Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths'][0]
g_inters = [g_inter_example]
Ib_str = params['backgrndI_strengths']  # list -> background input sweep
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
backgroundInputOnPeakPower(data_df, g, sI, Iext_dur, Iext_str, "baseline", "Potential", thalamus_source, figure_dir)

# %%
"""
5) Baseline power spectrum of the excitatory populations vs coupling strength
   (summed potential, computed over the baseline window only). One subplot per
   excitatory population, one curve per swept coupling value.
"""

# 5.1) Baseline spectrum vs global coupling g
sI = 0.6
Ib_str = params['backgrndI_strengths'][0]
Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths'][0]
g_inter = params['g_inter']
g_values = params['coupling_strengths']
baseline_spectrum_by_coupling('g', g_values, None, g_inter, sI, Ib_str,
    Iext_dur, Iext_str, input_onset, step_size, sample_dur, offset,
    thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type,
    raw_dir, figure_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)

# %%

# 5.2) Baseline spectrum vs inter-cortical coupling g_inter
sI = sI_example
Ib_str = Ib_example
Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths'][0]
g = g_example
g_inter_values = [0.4, 0.6, 0.8, 1.0]  # set to the g_inter values on disk
baseline_spectrum_by_coupling('g_inter', g_inter_values, g, None, sI, Ib_str,
    Iext_dur, Iext_str, input_onset, step_size, sample_dur, offset,
    thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type,
    raw_dir, figure_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)

# %%
'''
3.1) Background input strength
- x-axis: input strength
- y-axis: PSP of E3 population
- one subplot for each sI value
'''

g_inters = [g_inter_example]
sI = params['strength_I']
Ib_str = params['backgrndI_strengths']
Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths'][-1]
g = params['coupling_strengths']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
population = 'E3S2'
rate_measure = 'baseline' #('lateLongterm', 'immediateLongterm', 'duringInput', or 'baseline')
multisI_couplingOnMinmaxRate(data_df, sI, Ib_str, population, rate_measure, thalamus_source, figure_dir)
multisI_couplingOnDiffRate(data_df, sI, Ib_str, population, rate_measure, thalamus_source, figure_dir)

# %%
'''
BACKGROUND INPUT
3.2) all layers
'''
Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths'][0]
Ib_str = params['backgrndI_strengths'][0]
sI = 0.6
g = params['coupling_strengths']
g_inters = [params['g_inter']]
rate_measure = 'baseline' #('lateLongterm', 'immediateLongterm', 'duringInput', or 'baseline')

data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
multiLayer_couplingOnMinMaxRate(data_df, Iext_dur, Iext_str, Ib_str, sI, rate_measure, thalamus_source, figure_dir)


# %%  
'''
Fixed point plot
In this plot I use the potential instead of the firing rates.
    - x axis: input strengths 
    - y axis: response (steady-state/longterm) in mV
'''

# load time series data
Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths']
sI = sI_example
Ib_str = Ib_example
g = params['coupling_strengths']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)

potential_measure = 'lateLongterm'
population = 'E1S2'
inputStrengthOnminMaxpotential(data_df, Iext_str, Ib_str, sI, population, potential_measure, thalamus_source, figure_dir)


# %% Line Plot 
'''
Look at the change of rate difference comparing populations
- x axis: input duration 
- y axis: firing rate diff
- hue: populations 
'''



# %% Heatmap gE versus gI

'''
Look at the interaction between global coupling and E-I balance in different populations.
'''

Iext_dur = params['input_durations'][0]
Iext_str = params['input_strengths']  # single value in this batch -> one row of panels
g = params['coupling_strengths']
Ib_str = params['backgrndI_strengths'][0]
sI = params['strength_I']
data_df = load_all_derivatives(Iext_dur, Iext_str, g, g_inters, sI, Ib_str, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts, input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm, Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)

rate_measure = 'longtermVSbaseline_rate'

populations = np.array(['E1', 'E2', 'E3', 'E4']) 
#populations = np.array(['P1', 'P2', 'P3', 'P4']) 
#populations = np.array(['S1', 'S2', 'S3', 'S4', 'V1']) 

fig, axes = plt.subplots(len(Iext_str), len(populations), figsize=(20,15) ,sharex=True, sharey=True)
# with a single input strength plt.subplots returns a 1-D array; the loop below indexes [row, col]
axes = np.atleast_2d(axes)

# Create a single colorbar axis
#cbar_ax = fig.add_axes([1.01, 0.3, 0.02, 0.4])
#cbar_ax.set_title(rate_measure)
#cbar_ax.tick_params(labelsize=12) 
fig.suptitle(f'Effect of Input Strength, Global Coupling and E-I balance on {rate_measure} ', fontsize=16)

for i,input_s in enumerate(Iext_str):
    for j,p in enumerate(populations):

        minmax_df = data_df[data_df['InputStrength']==input_s]
        minmax_df = minmax_df[minmax_df['InputDuration']==Iext_dur]
        minmax_df = minmax_df[minmax_df['BckgndInputStrength']==Ib_str]
        minmax_df = minmax_df[minmax_df['population']==p]
        #minmax_df['InputDuration'] = minmax_df['InputDuration'].round(4)

        data_heatmap = minmax_df.pivot(index='globalCoupling', columns='strength_I', values=rate_measure)

        if (minmax_df[rate_measure].isna() | (minmax_df[rate_measure] == 0)).all().all():
            print('black')
            sns.heatmap(data_heatmap, cmap=ListedColormap(['black']), ax=axes[i,j], norm = Normalize(vmin=0, vmax=1))
        else:
            sns.heatmap(data_heatmap, cmap='vlag', ax=axes[i,j], center=0)

        cbar = axes[i, j].collections[0].colorbar
        # here set the labelsize by 20
        cbar.ax.tick_params(labelsize=12)
        axes[i, j].invert_yaxis()
        axes[i, j].set_ylabel('')
        axes[i, j].set_xlabel('')
        axes[i, j].tick_params(axis='both', labelsize=12)
        axes[len(Iext_str)-1, j].set_xlabel('E-I balance')
        axes[0,j].set_title(f'pop: {p}')
        axes[i,0].set_ylabel(f'input strength: {input_s}', rotation=0, labelpad=60)

fig.text(0, 0.2, 'coupling strength', va='center', rotation='vertical')
fig.text(0, 0.5, 'coupling strength', va='center', rotation='vertical')
fig.text(0, 0.83, 'coupling strength', va='center', rotation='vertical')

plt.tight_layout(h_pad=1)
figure_name = f'gvssI_{populations[0][0]}pop_{rate_measure}_tauVisual.png'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.close()

# %%
'''
Choose the 4 most important variables:
- global g
- E-I balance
- input strength
- input duration (or background strength)
Output value:
- difference of max. vs min. long term versus baseline rate --> oscillatory behaviour

For the output population of S1 (E3), for every global g and E-I balance parameter
plot the input strength and duration. 

The goal is to create an overview of the system behaviour: 
global x an y are global coupling strength (g) and E-I ratio
- plot 1 (lineplot): fixed points and limit cycles
    x-axis: background input 
    y-axis: long term firing rate (min and max)
- plot 2 (heatmap): long-term behaviour after input
    x-axis: stimulus duration 
    y-axis: stimulus strength
    value: difference between min and max value
- plot 3 (heatmap): systems responsiveness
    x-axis: stimulus duration 
    y-axis: stimulus strength
    value: 
        memory (long-term behaviour the same as during input)
        transfer (long-term behaviour goes back to base behaviour)
        non-responsive (no change of activity during input)

'''


# create figure with multiple subplots
fig, axes = plt.subplots(ncols=2, nrows=2)

# iterate through the global coupling strength values 
#for g in coupling_strengths:



# %% 
# PLOT 1: fixed points and limit cycles
# These last three cells share summary_df / gs / sIs / population. summary_df used to be a
# leftover global that no longer exists, so load it here the way every other cell does.
# gs is subsampled because this batch sweeps 40 coupling values and one panel per value is
# unreadable; sIs comes from the grid rather than a hardcoded list that predates it.
sIs = params['strength_I']
gs = params['coupling_strengths'][::8]

input_duration = params['input_durations'][0]
backgrndI_strengths = params['backgrndI_strengths']
input_strengths = params['input_strengths']
population = 'E3'
summary_df = load_all_derivatives(
    input_duration, input_strengths, params['coupling_strengths'], [g_inter_example], sIs,
    backgrndI_strengths, input_onset, thal_cellcounts, bI_cellcounts, extI_cellcounts,
    input_type, processed_dir, suffix=suffix, g_thalPOm=g_thalPOm,
    Ib_noise_std=Ib_noise_std, g_thal=g_thal, sI_thal=sI_thal)
data_df = summary_df[summary_df['population']==population]
data_df = data_df[data_df['InputDuration']==input_duration]
data_df = data_df[data_df['InputStrength'].isin(input_strengths)]

# plot results (balance E-I and global coupling on the global axis)
fig, axs = plt.subplots(ncols=len(sIs), nrows=len(gs), figsize=(20,10), squeeze=False)

# Create a colormap
cmap = plt.get_cmap('Dark2', len(input_strengths))
cmap_max = plt.get_cmap('Dark2', len(input_strengths))
# step001 names the late-longterm window *_lateLongterm, not *_longterm
data_df = data_df.copy()
data_df['minPotential_lateLongterm_mV'] = data_df['minPotential_lateLongterm'] * 1e3
data_df['maxPotential_lateLongterm_mV'] = data_df['maxPotential_lateLongterm'] * 1e3
input_strengths = data_df['InputStrength'].unique()

for i, g in enumerate(gs):
    coupling_df = data_df[data_df['globalCoupling']==g]
    for j, sI in enumerate(sIs):
        balance_df = coupling_df[coupling_df['strength_I']==sI]
        ax = axs[i][j]
        if j ==0:
            ax.set_ylabel(f'g={g}\nRate (Hz)')
        if i == len(gs)-1:
            ax.set_xlabel(f'Background Input Strength\n E-I Balance={sI}')
        for k, s in enumerate(input_strengths):
            Istrength_df = balance_df[balance_df['InputStrength']==s]
            color = cmap(k / (len(input_strengths)))  # Normalize i to [0, 1] for colormap
            color_max = cmap_max(k / (len(input_strengths)))  # Normalize i to [0, 1] for colormap
            ax.plot(Istrength_df['BckgndInputStrength'], Istrength_df['minRate_lateLongterm'], label=s, color=color)
            ax.plot(Istrength_df['BckgndInputStrength'], Istrength_df['maxRate_lateLongterm'], color=color_max)

#sns.lineplot(data_df, y='maxPotential_longterm', x='coupling_strength', hue='InputStrength')
#axs.set_xlim([0,100])
sns.despine(trim=True)
axs[0][-1].legend(title='Input Strength', loc='right')
plt.tight_layout() 
figure_name = f'BackgroundSteadyState_pop{population}_tauVisual_{thalamus_source}.pdf'
plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches='tight')
plt.close()

# %%
# PLOT 2: behaviour during input

Idur = params['input_durations'][0]
population = 'E3'

fig, axs = plt.subplots(ncols=len(sIs), nrows=len(gs), figsize=(20,10), squeeze=False)

data_df = summary_df[summary_df['population']==population]
data_df = data_df[data_df['InputDuration']==Idur]

for i, g in enumerate(gs):
    coupling_df = data_df[data_df['globalCoupling']==g]
    for j, sI in enumerate(sIs):
        balance_df = coupling_df[coupling_df['strength_I']==sI]
        ax = axs[i][j]

        data_heatmap = balance_df.pivot(index='BckgndInputStrength',columns='InputStrength', values='diffRate_duringInput')
        sns.heatmap(data_heatmap, cmap='magma', ax=ax, vmin=0, vmax=40)
        ax.invert_yaxis()
        if j ==0:
            ax.set_ylabel(f'g={g}\n Background Input')
        if i == len(gs)-1:
            ax.set_xlabel(f'Input Strength\n E-I Balance={sI}')
        else:
            ax.set_xlabel('')
            ax.set_ylabel('')

# %%
# PLOT 3: long-term behaviour compared to baseline

Idur = params['input_durations'][0]
population = 'E3'

fig, axs = plt.subplots(ncols=len(sIs), nrows=len(gs), figsize=(20,10), squeeze=False)

data_df = summary_df[summary_df['population']==population]
data_df = data_df[data_df['InputDuration']==Idur]

for i, g in enumerate(gs):
    coupling_df = data_df[data_df['globalCoupling']==g]
    for j, sI in enumerate(sIs):
        balance_df = coupling_df[coupling_df['strength_I']==sI]
        ax = axs[i][j]

        data_heatmap = balance_df.pivot(index='BckgndInputStrength',columns='InputStrength', values='longtermVSbaseline_rate')
        sns.heatmap(data_heatmap, cmap='magma', ax=ax, vmin=0, vmax=10)
        ax.invert_yaxis()
        print(j)
        if i ==0:
            ax.set_ylabel(f'g={g}\n Background Input')
        if i == len(gs)-1:
            ax.set_xlabel(f'Input Strength\n E-I Balance={sI}')
        else:
            ax.set_xlabel('')
            ax.set_ylabel('')
# %%

# %%
# The modulatory (VIP) input figures used to live here as section 1.1e. They re-simulate
# the model at a hand-set operating point instead of reading this batch, so they moved to
# their own script: Analysis/step005_plot_modulatory_input.py
