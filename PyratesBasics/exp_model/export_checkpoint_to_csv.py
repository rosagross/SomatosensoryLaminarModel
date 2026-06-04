# %%

import os
import re
import glob
import numpy as np
import pandas as pd


# ============================================================
# USER SETTINGS
# ============================================================
WDDIR = os.getenv("WDDIR")
if not WDDIR:
    raise EnvironmentError("WDDIR is not set. Please export WDDIR to the repository root.")

output_dir = os.path.join(WDDIR, "PyratesBasics", "exp_model", "complete_model_continuations")
checkpoint_root = os.path.join(output_dir, "checkpoints")

# Use a specific checkpoint timestamp folder (e.g., "20260521T084433")
# Leave as None to use the latest checkpoint
checkpoint_timestamp = None

# Output CSV path (leave None to write inside latest checkpoint folder)
output_csv = None

# Mapping CSV (u_sols) to map U(i) -> state names for population-level export
state_map_csv = "/data/hu_grossmannr/Desktop/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel/PyratesBasics/exp_model/complete_model_continuations/stability_write_test_g0_1p2.csv"

# Continuation parameter name (PAR(14) for this model)
cont_param = "G/g_definition/g_input"

# %%
# ============================================================
# HELPERS
# ============================================================
def resolve_checkpoint_dir(root_dir, timestamp=None):
    if timestamp:
        candidate = os.path.join(root_dir, timestamp)
        if not os.path.isdir(candidate):
            raise FileNotFoundError(f"Checkpoint folder not found: {candidate}")
        return candidate

    if not os.path.isdir(root_dir):
        raise FileNotFoundError(f"Checkpoint root not found: {root_dir}")

    entries = [d for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))]
    if not entries:
        raise FileNotFoundError(f"No checkpoint folders found in {root_dir}")

    entries.sort()
    return os.path.join(root_dir, entries[-1])


def read_ndim(checkpoint_dir):
    c_ivp_path = os.path.join(checkpoint_dir, "c.ivp")
    if not os.path.isfile(c_ivp_path):
        c_matches = glob.glob(os.path.join(checkpoint_dir, "c.*"))
        c_ivp_path = c_matches[0] if c_matches else None

    if not c_ivp_path or not os.path.isfile(c_ivp_path):
        return None

    try:
        with open(c_ivp_path, "r") as c_file:
            for line in c_file:
                match = re.search(r"NDIM\s*=\s*(\d+)", line)
                if match:
                    return int(match.group(1))
    except OSError:
        return None
    return None


def parse_fort7(path):
    rows = []
    max_u = 0
    with open(path, "r") as f7:
        for line in f7:
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                br = int(parts[0])
                pt = int(parts[1])
            except ValueError:
                continue

            ty = parts[2]
            try:
                lab = int(parts[3])
            except ValueError:
                lab = np.nan

            try:
                par = float(parts[4])
                l2_norm = float(parts[5])
            except ValueError:
                continue

            u_vals = []
            for token in parts[6:]:
                try:
                    u_vals.append(float(token))
                except ValueError:
                    break
            max_u = max(max_u, len(u_vals))
            rows.append((br, pt, ty, lab, par, l2_norm, u_vals))

    if not rows:
        raise ValueError(f"No continuation data lines found in {path}")

    data = {
        "BR": [],
        "PT": [],
        "TY": [],
        "LAB": [],
        "PAR(14)": [],
        "L2-NORM": [],
    }
    for br, pt, ty, lab, par, l2_norm, _ in rows:
        data["BR"].append(br)
        data["PT"].append(pt)
        data["TY"].append(str(ty))
        data["LAB"].append(lab)
        data["PAR(14)"].append(par)
        data["L2-NORM"].append(l2_norm)

    for i in range(1, max_u + 1):
        data[f"U({i})"] = []

    for _, _, _, _, _, _, u_vals in rows:
        for i in range(1, max_u + 1):
            value = u_vals[i - 1] if i - 1 < len(u_vals) else np.nan
            data[f"U({i})"].append(value)

    return pd.DataFrame(data)


def parse_fort9(path):
    pt_to_stable_by_branch = {}
    pt_to_stable_abs = {}
    if not path or not os.path.isfile(path):
        return pt_to_stable_by_branch, pt_to_stable_abs

    try:
        with open(path, "r") as f9:
            for line in f9:
                match = re.match(r"^\s*(\d+)\s+(-?\d+)\s+Eigenvalues\s+:\s+Stable:(\d+)", line)
                if match:
                    br = int(match.group(1))
                    pt_abs = abs(int(match.group(2)))
                    stable_dims = int(match.group(3))
                    pt_to_stable_by_branch[(br, pt_abs)] = stable_dims
                    pt_to_stable_abs[pt_abs] = stable_dims
    except OSError:
        return {}, {}

    return pt_to_stable_by_branch, pt_to_stable_abs


def apply_stability(df, ndim, pt_to_stable_by_branch, pt_to_stable_abs):
    if ndim is None:
        df["stability"] = "unknown"
        return df
    if not pt_to_stable_by_branch and not pt_to_stable_abs:
        df["stability"] = "unknown"
        return df

    stable_points = []
    for _, row in df.iterrows():
        br = int(row["BR"])
        pt_abs = abs(int(row["PT"]))
        par = float(row["PAR(14)"])
        s = pt_to_stable_by_branch.get((br, pt_abs))
        if s is None:
            s = pt_to_stable_abs.get(pt_abs)
        if s is not None:
            stable_points.append((par, s))

    stability_labels = []
    for _, row in df.iterrows():
        br = int(row["BR"])
        pt_abs = abs(int(row["PT"]))
        par = float(row["PAR(14)"])
        s = pt_to_stable_by_branch.get((br, pt_abs))
        if s is None:
            s = pt_to_stable_abs.get(pt_abs)

        if s is None and stable_points:
            par_s, s_near = min(stable_points, key=lambda x: abs(x[0] - par))
            tol = max(1e-6, 1e-4 * max(1.0, abs(par)))
            fallback_tol = max(1e-3, 5.0 * tol)
            if abs(par_s - par) <= fallback_tol:
                s = s_near

        if s is None:
            stability_labels.append("unknown")
        elif float(s) >= float(ndim):
            stability_labels.append("stable")
        else:
            stability_labels.append("unstable")

    df["stability"] = stability_labels
    return df


def add_bifurcation_labels(df):
    labels = []
    for ty in df["TY"].astype(str).values:
        ty_clean = ty.strip()
        if ty_clean in {"0", "", "RG"}:
            labels.append("")
        else:
            labels.append(ty_clean)
    df["bifurcation"] = labels
    return df


def load_state_names(csv_path, ndim):
    if not csv_path or not os.path.isfile(csv_path):
        raise FileNotFoundError(f"State map CSV not found: {csv_path}")
    if ndim is None:
        raise ValueError("NDIM not available; cannot map U(i) to state names.")
    headers = pd.read_csv(csv_path, nrows=0).columns.tolist()
    if ndim > len(headers):
        raise ValueError(f"NDIM={ndim} exceeds header count in {csv_path} ({len(headers)}).")
    return headers[:ndim]


def build_state_df(df, state_names):
    data = {}
    for i, name in enumerate(state_names, start=1):
        col = f"U({i})"
        if col in df.columns:
            data[name] = df[col]
    if not data:
        raise ValueError("No U(i) columns found to map into state names.")
    return pd.DataFrame(data)


def compute_population_df(state_df, cells):
    n_cells = len(cells)
    rpo_names = [f"RPO_{cell}" for cell in cells]
    rpo_names_extended = rpo_names + ["RPO_Iext"]
    all_pot_cont = []
    for i, target_cell in enumerate(cells):
        if i == (n_cells - 2):
            potential_keys = [f"{target_cell}/{rpo_name_ext}/v" for rpo_name_ext in rpo_names_extended[: n_cells + 1]]
        else:
            potential_keys = [f"{target_cell}/{rpo_name}/v" for rpo_name in rpo_names[:n_cells]]
        if i in range(n_cells - 2):
            potential_keys += [f"{target_cell}/RPO_bI/v"]
        missing = [key for key in potential_keys if key not in state_df.columns]
        if missing:
            raise ValueError(f"Missing state variables for {target_cell}: {missing[:5]}...")
        sources = state_df[potential_keys]
        all_pot_cont.append(np.sum(sources, axis=1))
    all_pot_cont = np.array(all_pot_cont).T
    return pd.DataFrame(all_pot_cont, columns=cells)


# ============================================================
# MAIN
# ============================================================
checkpoint_dir = resolve_checkpoint_dir(checkpoint_root, checkpoint_timestamp)
fort7_path = os.path.join(checkpoint_dir, "fort.7")
if not os.path.isfile(fort7_path):
    raise FileNotFoundError(f"fort.7 not found in {checkpoint_dir}")

fort9_path = os.path.join(checkpoint_dir, "fort.9")

df = parse_fort7(fort7_path)
df = add_bifurcation_labels(df)

ndim = read_ndim(checkpoint_dir)
pt_to_stable_by_branch, pt_to_stable_abs = parse_fort9(fort9_path)
df = apply_stability(df, ndim, pt_to_stable_by_branch, pt_to_stable_abs)

cells = [
    "E3b", "PV3b", "SST3b", "VIP3b",
    "E1S1", "PV1S1", "SST1S1", "VIPS1", "E2S1", "PV2S1", "SST2S1", "E3S1", "PV3S1", "SST3S1", "E4S1", "PV4S1", "SST4S1",
    "E1S2", "PV1S2", "SST1S2", "VIPS2", "E2S2", "PV2S2", "SST2S2", "E3S2", "PV3S2", "SST3S2", "E4S2", "PV4S2", "SST4S2",
    "ThalE", "ThalI",
]

state_names = load_state_names(state_map_csv, ndim)
state_df = build_state_df(df, state_names)
pop_df = compute_population_df(state_df, cells)

cont_series = df["PAR(14)"].rename(cont_param)
out_df = pd.concat(
    [cont_series, pop_df, df["stability"], df["bifurcation"]],
    axis=1,
)

if output_csv is None:
    output_csv = os.path.join(checkpoint_dir, "checkpoint_population_export.csv")

out_df.to_csv(output_csv, index=False)

print(f"Checkpoint: {checkpoint_dir}")
print(f"Rows: {len(out_df)}")
print(f"Wrote CSV: {output_csv}")
