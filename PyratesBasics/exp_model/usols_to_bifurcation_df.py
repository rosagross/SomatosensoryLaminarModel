# %% [1] User parameters — edit this cell before running
import os
WDDIR = os.getenv("WDDIR", "/data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel")

input_path  = os.path.join(WDDIR, "PyratesBasics", "exp_model",
                           "complete_model_continuations",
                           "complete_model_bifurcation_sI026_sara.csv")
output_path = os.path.join(WDDIR, "PyratesBasics", "exp_model",
                           "complete_model_continuations",
                           "complete_model_bifurcation_sI026_sara_processed.csv")
cont_param  = "G/g_definition/g_input"

# Optional: set paths to fort.7 and fort.9 for stability parsing.
# Leave as None to skip (output will have NaN stability / empty bifurcation).
fort7_path  = None
fort9_path  = None

# %% [2] Model constants — must match SomatoModelPyrates.__init__
import numpy as np
import pandas as pd
import re

cells = [
    'E3b', 'PV3b', 'SST3b', 'VIP3b',                          # A3b
    'E1S1', 'PV1S1', 'SST1S1', 'VIPS1',                       # S1 L1
    'E2S1', 'PV2S1', 'SST2S1',                                 # S1 L2
    'E3S1', 'PV3S1', 'SST3S1',                                 # S1 L3
    'E4S1', 'PV4S1', 'SST4S1',                                 # S1 L4
    'E1S2', 'PV1S2', 'SST1S2', 'VIPS2',                       # S2 L1
    'E2S2', 'PV2S2', 'SST2S2',                                 # S2 L2
    'E3S2', 'PV3S2', 'SST3S2',                                 # S2 L3
    'E4S2', 'PV4S2', 'SST4S2',                                 # S2 L4
    'ThalE', 'ThalI',                                          # Thalamus
]
N_cells            = len(cells)                                # 32
rpo_names          = ['RPO_' + c for c in cells]
rpo_names_extended = rpo_names + ['RPO_Iext']

# %% [3] Load raw u_sols CSV
u_sols = pd.read_csv(input_path)
print(f"Loaded {input_path}")
print(f"  shape: {u_sols.shape}")
print(f"  cont_param present: {cont_param in u_sols.columns}")
print(f"  stability column present: {'stability' in u_sols.columns}")
print(f"  bifurcation column present: {'bifurcation' in u_sols.columns}")

# %% [4] Replicate continuation_df() logic
# For each cell, sum all incoming /v potential columns.
all_pot_cont = []

for i, target_cell in enumerate(cells):
    if i == N_cells - 2:   # ThalE: include RPO_Iext
        potential_keys = [f'{target_cell}/{rpo}/v' for rpo in rpo_names_extended[:N_cells + 1]]
    else:
        potential_keys = [f'{target_cell}/{rpo}/v' for rpo in rpo_names[:N_cells]]

    if i in range(N_cells - 2):   # non-thalamic: include background RPO
        potential_keys += [f'{target_cell}/RPO_bI/v']

    missing = [k for k in potential_keys if k not in u_sols.columns]
    if missing:
        print(f"WARNING: {len(missing)} expected columns missing for {target_cell}. First few: {missing[:3]}")

    present_keys = [k for k in potential_keys if k in u_sols.columns]
    all_pot_cont.append(u_sols[present_keys].sum(axis=1).values)

all_pot_cont = np.array(all_pot_cont).T  # shape: (n_points, n_cells)
cont_df = pd.DataFrame(all_pot_cont, columns=cells)
cont_df.insert(0, cont_param, u_sols[cont_param].values)

# %% [5] Stability and bifurcation columns

def inject_stability(u_sols_df, cont_param, fort7_path, fort9_path, ndim):
    """Replicate _inject_stability_from_auto(): parse fort.7/fort.9 and label each point."""
    pt_to_stable_by_branch = {}
    pt_to_stable_abs = {}
    try:
        with open(fort9_path, "r") as f9:
            for line in f9:
                m = re.match(r"^\s*(\d+)\s+(-?\d+)\s+Eigenvalues\s+:\s+Stable:(\d+)", line)
                if m:
                    br        = int(m.group(1))
                    pt_abs    = abs(int(m.group(2)))
                    stable_d  = int(m.group(3))
                    pt_to_stable_by_branch[(br, pt_abs)] = stable_d
                    pt_to_stable_abs[pt_abs]              = stable_d
    except OSError:
        pass

    fort_points = []
    try:
        with open(fort7_path, "r") as f7:
            for line in f7:
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    br  = int(parts[0])
                    pt  = int(parts[1])
                    par = float(parts[4])
                    fort_points.append((br, pt, par))
                except ValueError:
                    continue
    except OSError:
        pass

    stable_points = []
    for br, pt, par in fort_points:
        s = pt_to_stable_by_branch.get((br, abs(pt)))
        if s is None:
            s = pt_to_stable_abs.get(abs(pt))
        if s is not None:
            stable_points.append((par, s))

    stable_dims    = []
    matched        = 0
    for g in u_sols_df[cont_param].values:
        if not fort_points:
            stable_dims.append(np.nan)
            continue
        br, pt, par = min(fort_points, key=lambda x: abs(x[2] - g))
        tol = max(1e-6, 1e-4 * max(1.0, abs(g)))
        s   = None
        if abs(par - g) <= tol:
            s = pt_to_stable_by_branch.get((br, abs(pt)))
            if s is None:
                s = pt_to_stable_abs.get(abs(pt))
        if s is None and stable_points:
            par_s, s_near = min(stable_points, key=lambda x: abs(x[0] - g))
            if abs(par_s - g) <= max(1e-3, 5.0 * tol):
                s = s_near
        if s is None:
            stable_dims.append(np.nan)
        else:
            stable_dims.append(float(s))
            matched += 1

    stable_dims = np.array(stable_dims, dtype=float)
    labels = []
    for s in stable_dims:
        if np.isnan(s):
            labels.append(np.nan)
        elif s >= float(ndim):
            labels.append("stable")
        else:
            labels.append("unstable")

    print(f"Stability mapping: matched={matched}, unmatched={len(stable_dims)-matched}")
    return labels

if 'stability' in u_sols.columns:
    cont_df['stability']   = u_sols['stability'].values
    cont_df['bifurcation'] = u_sols['bifurcation'].values if 'bifurcation' in u_sols.columns else ''
    print("Stability copied from u_sols.")
elif fort7_path is not None and fort9_path is not None:
    # ndim = number of unique state variable columns (each has /v and /i → 2 per RPO-cell pair)
    state_cols = [c for c in u_sols.columns if c != cont_param
                  and 'stability' not in c and 'bifurcation' not in c]
    ndim = len(state_cols)
    print(f"Running stability parsing (ndim estimate = {ndim}) ...")
    stability_labels = inject_stability(u_sols, cont_param, fort7_path, fort9_path, ndim)
    cont_df['stability']   = stability_labels
    cont_df['bifurcation'] = u_sols['bifurcation'].values if 'bifurcation' in u_sols.columns else ''
else:
    cont_df['stability']   = np.nan
    cont_df['bifurcation'] = ''
    print("No stability info available. Columns set to NaN/empty.")

# %% [6] Save and summarise
cont_df.to_csv(output_path, index=False)
print(f"\nSaved to {output_path}")
print(f"Output shape: {cont_df.shape}")
print(f"Columns: {list(cont_df.columns)}")
print(f"\ng_input range: {cont_df[cont_param].min():.4f} – {cont_df[cont_param].max():.4f}")
if cont_df['stability'].notna().any():
    print(f"Stability counts: {cont_df['stability'].value_counts().to_dict()}")
