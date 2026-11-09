"""
File: dipole_parameterization.py
Author: Rosa Grossmann
Description:
    Layer-resolved parameterization of the dipole lengths with 10 parameters, and a check
    that prints what the parameters produce. No simulation is run and run_optimization_dipole.py
    does not use this file yet.

    The length of a dipole depends on the target (A3b, or layer L2/3, L4, L5, L6 of A1 and S2),
    on the source population, and on the source cell type. The 10 parameters (variables below):

     1  L_BASE        mm      cortical sources -> L2/3 of A1 and S2 (the base length)
     2  F_L5          factor  the same projections -> L5 of A1 and S2 (must be >= F_L5_MIN = 2)
     3  F_L6          factor  the same projections -> L6 of A1 and S2
     4  F_L4_A1       factor  x L_BASE, everything (thalamus included) -> L4 of A1
     5  F_L4_S2       factor  x L_BASE, everything (thalamus included) -> L4 of S2
     6  L_THALE_CTX   mm      ThalE    -> L2/3 of A1 and S2 (L5 and L6: x F_L5, x F_L6)
     7  L_THALPOM_CTX mm      ThalPOm  -> L2/3 of A1 and S2 (L5 and L6: x F_L5, x F_L6)
     8  L_THALE_A3B   mm      ThalE    -> A3b
     9  L_THALPOM_A3B mm      ThalPOm  -> A3b
    10  L_A3B_CORT    mm      all cortical sources -> A3b

    Source ratios (cortical sources, parameters 1-3, 6-7 and 10): E = SST = 1, PV = 0.5, VIP = 0;
    sources in L4 (S1 L4, S2 L4): 0.5 for E, SST and PV. Targets A1 L4 and S2 L4 (parameters 4, 5)
    have no ratios: every non-VIP source gets the same length. ThalI, ThalPOmI and VIP are always 0.
    Orientation and the E cell count weighting are not touched: only dipole_lengths is built.

    Usage:
        python dipole_parameterization.py                         # print the check
        python dipole_parameterization.py --plot lengths.png      # + heatmaps (file, parameterized, difference)
        python dipole_parameterization.py --write-json out.json   # + the parameter file with these lengths
    The warm values below were read off dipole_parameters_flippedPVSST.json.
"""

import argparse
import copy
import json
import os
import sys
import numpy as np

WDDIR = os.getenv("WDDIR", "/data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel")
SIM_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, SIM_DIR)
from somato_model import read_dipole_params, DIPOLE_SOURCE_PATHS, DIPOLE_TARGET_LAYERS   # noqa: E402

REFERENCE_JSON = os.path.join(WDDIR, "EEGSimulation", "dipole_parameters_flippedPVSST.json")

# ── the 10 parameters (warm values = current file) ─────────────────────────────
L_BASE        = 2.0     # 1  mm
F_L5          = 5.0     # 2  factor on L_BASE, target L5 (>= F_L5_MIN)
F_L6          = 3.5     # 3  factor on L_BASE, target L6
F_L4_A1       = 0.05    # 4  factor on L_BASE, target A1 L4 (0.1 mm in the file)
F_L4_S2       = 0.05    # 5  factor on L_BASE, target S2 L4
L_THALE_CTX   = 1.0     # 6  mm
L_THALPOM_CTX = 2.0     # 7  mm
L_THALE_A3B   = 0.1     # 8  mm
L_THALPOM_A3B = 2.0     # 9  mm
L_A3B_CORT    = 2.0     # 10 mm

PARAM_NAMES = ["L_BASE", "F_L5", "F_L6", "F_L4_A1", "F_L4_S2",
               "L_THALE_CTX", "L_THALPOM_CTX", "L_THALE_A3B", "L_THALPOM_A3B", "L_A3B_CORT"]
PARAM_UNITS = ["mm", "x base", "x base", "x base", "x base", "mm", "mm", "mm", "mm", "mm"]

# ── fixed rules ────────────────────────────────────────────────────────────────
TYPE_RATIO = {"E": 1.0, "SST": 1.0, "PV": 0.5, "VIP": 0.0}   # cortical sources outside L4
L4_SOURCE_RATIO = 0.5                                         # E, SST, PV in S1 L4 / S2 L4 (VIP does not exist there)
F_L5_MIN = 2.0

ROWS = [("A3b", None)] + [("A1", l) for l in DIPOLE_TARGET_LAYERS] + [("S2", l) for l in DIPOLE_TARGET_LAYERS]
ROW_NAMES = ["A3b"] + [f"{a} {l}" for a in ("A1", "S2") for l in DIPOLE_TARGET_LAYERS]


def current_params():
    """The 10 module-level variables as a dict."""
    g = globals()
    return {k: float(g[k]) for k in PARAM_NAMES}


def _validate(p):
    bad = [k for k in PARAM_NAMES if not np.isfinite(p[k]) or p[k] < 0]
    if bad:
        raise ValueError(f"parameters must be finite and >= 0: {bad}")
    if p["F_L5"] < F_L5_MIN:
        raise ValueError(f"F_L5 = {p['F_L5']} < {F_L5_MIN}: the L5 length must be at least twice the L2/3 length")


def source_ratio(path):
    """Ratio of a cortical source (path from DIPOLE_SOURCE_PATHS); 0 for VIP, 0.5 for L4 sources."""
    celltype = path[-1]
    if celltype == "VIP":
        return 0.0
    if len(path) == 3 and path[1] == "L4":
        return L4_SOURCE_RATIO
    return TYPE_RATIO[celltype]


def length(p, target, layer, path):
    """Dipole length (mm) of source `path` onto target ('A3b' or 'A1'/'S2' with `layer`)."""
    celltype = path[-1]
    thal = path[0] == "Thalamus"
    if target == "A3b":
        if thal:
            return {"ThalE": p["L_THALE_A3B"], "ThalPOm": p["L_THALPOM_A3B"]}.get(celltype, 0.0)
        return p["L_A3B_CORT"] * source_ratio(path)
    if layer == "L4":
        f4 = p["F_L4_A1"] if target == "A1" else p["F_L4_S2"]
        if thal:
            return p["L_BASE"] * f4 if celltype in ("ThalE", "ThalPOm") else 0.0
        return 0.0 if path[-1] == "VIP" else p["L_BASE"] * f4     # no type or L4-source ratios in L4 targets
    f = {"L2/3": 1.0, "L5": p["F_L5"], "L6": p["F_L6"]}[layer]
    if thal:
        return {"ThalE": p["L_THALE_CTX"], "ThalPOm": p["L_THALPOM_CTX"]}.get(celltype, 0.0) * f
    return p["L_BASE"] * f * source_ratio(path)


def _set(tree, path, value):
    for k in path[:-1]:
        tree = tree.setdefault(k, {})
    tree[path[-1]] = value


def _block(p, target, layer):
    """Nested {source region: ...: cell type: length} dict of one target, in the file's schema."""
    block = {}
    for path in DIPOLE_SOURCE_PATHS:
        _set(block, path, length(p, target, layer, path))
    return block


def build_lengths(p=None):
    """`dipole_lengths` in the schema of the parameter JSON (A3b / A1 / S2 -> layer -> source -> ...)."""
    p = current_params() if p is None else {**current_params(), **p}
    _validate(p)
    out = {"A3b": _block(p, "A3b", None)}
    for target in ("A1", "S2"):
        out[target] = {layer: _block(p, target, layer) for layer in DIPOLE_TARGET_LAYERS}
    return out


def build_params_json(out_path, p=None, reference_path=REFERENCE_JSON):
    """Copy of the reference parameter file with dipole_lengths replaced (orientation etc. unchanged)."""
    with open(reference_path) as f:
        ref = json.load(f)
    new = copy.deepcopy(ref)
    new["dipole_lengths"] = build_lengths(p)
    with open(out_path, "w") as f:
        json.dump(new, f, indent=2)
    read_dipole_params(out_path)          # raises if the model could not read it
    return out_path


def lengths_matrix(lengths):
    """(9, 34): rows A3b, A1 L2/3-L6, S2 L2/3-L6; columns in get_population_labels() order."""
    tmp = {"dipole_lengths": lengths, "dipole_orientation": lengths}
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_tmp_dipole_lengths.json")
    with open(path, "w") as f:
        json.dump(tmp, f)
    try:
        flat = read_dipole_params(path)["dipole_lengths"]
    finally:
        os.remove(path)
    return np.vstack([flat["A3b"]] + list(flat["A1"]) + list(flat["S2"]))


def reference_matrix(path=REFERENCE_JSON):
    flat = read_dipole_params(path)["dipole_lengths"]
    return np.vstack([flat["A3b"]] + list(flat["A1"]) + list(flat["S2"]))


# ── linear form for the optimizer ──────────────────────────────────────────────
# For fixed F_L5 and F_L6 the lengths are linear in 8 coefficients (mm). L_L4_A1 = L_BASE * F_L4_A1 and
# L_L4_S2 = L_BASE * F_L4_S2 are the lengths of the L4 targets themselves.
COEFF_NAMES = ["L_BASE", "L_L4_A1", "L_L4_S2", "L_THALE_CTX", "L_THALPOM_CTX",
               "L_THALE_A3B", "L_THALPOM_A3B", "L_A3B_CORT"]
_ROW_L4 = {"A1": 2, "S2": 6}                  # rows of the L4 targets in the (9, 34) matrix
_ROWS_L5, _ROWS_L6 = [3, 7], [4, 8]           # rows of the L5 / L6 targets of A1 and S2
_unit_cache = {}


def length_matrix(p):
    """(9, 34) length matrix straight from the 10 parameters (no file access); rows as in ROWS."""
    return np.array([[length(p, t, l, path) for path in DIPOLE_SOURCE_PATHS] for t, l in ROWS])


def _unit_matrices():
    """Length matrices per unit coefficient at F_L5 = F_L6 = 1 (8, 9, 34)."""
    if "u" not in _unit_cache:
        zero = {k: 0.0 for k in PARAM_NAMES}
        unit = []
        for name in COEFF_NAMES:
            p = {**zero, "F_L5": 1.0, "F_L6": 1.0}
            if name == "L_BASE":
                p["L_BASE"] = 1.0
                M = length_matrix(p)                               # F_L4 = 0: the L4 rows are 0
            elif name in ("L_L4_A1", "L_L4_S2"):
                area = name[-2:]
                p.update({"L_BASE": 1.0, f"F_L4_{area}": 1.0})
                M = np.zeros((9, len(DIPOLE_SOURCE_PATHS)))
                M[_ROW_L4[area]] = length_matrix(p)[_ROW_L4[area]]
            else:
                p[name] = 1.0
                M = length_matrix(p)
            unit.append(M)
        _unit_cache["u"] = np.array(unit)
    return _unit_cache["u"]


def coefficient_basis(F_L5=None, F_L6=None):
    """(8, 9, 34): length matrix per unit of each coefficient in COEFF_NAMES.

    The L5 / L6 target rows of L_BASE, L_THALE_CTX and L_THALPOM_CTX are multiplied by F_L5 / F_L6.
    """
    F_L5 = globals()["F_L5"] if F_L5 is None else F_L5
    F_L6 = globals()["F_L6"] if F_L6 is None else F_L6
    if F_L5 < F_L5_MIN:
        raise ValueError(f"F_L5 = {F_L5} < {F_L5_MIN}")
    rf = np.ones(9)
    rf[_ROWS_L5], rf[_ROWS_L6] = F_L5, F_L6
    B = _unit_matrices().copy()
    for g in (0, 3, 4):
        B[g] = B[g] * rf[:, None]
    return B


def lengths_from_coefficients(c, F_L5=None, F_L6=None):
    """(9, 34) length matrix = sum_g c_g B_g."""
    return np.tensordot(np.asarray(c, dtype=float), coefficient_basis(F_L5, F_L6), axes=1)


def coefficients_from_params(p=None):
    """The 8 coefficients of the 10 parameters in `p` (default: module values)."""
    p = current_params() if p is None else {**current_params(), **p}
    return np.array([p["L_BASE"], p["L_BASE"] * p["F_L4_A1"], p["L_BASE"] * p["F_L4_S2"], p["L_THALE_CTX"],
                     p["L_THALPOM_CTX"], p["L_THALE_A3B"], p["L_THALPOM_A3B"], p["L_A3B_CORT"]])


def params_from_coefficients(c, F_L5, F_L6):
    """The 10 parameters (dict) from the 8 coefficients and the two layer factors."""
    c = np.asarray(c, dtype=float)
    base = c[0] if c[0] > 0 else np.nan
    return {"L_BASE": float(c[0]), "F_L5": float(F_L5), "F_L6": float(F_L6),
            "F_L4_A1": float(c[1] / base), "F_L4_S2": float(c[2] / base),
            "L_THALE_CTX": float(c[3]), "L_THALPOM_CTX": float(c[4]),
            "L_THALE_A3B": float(c[5]), "L_THALPOM_A3B": float(c[6]), "L_A3B_CORT": float(c[7])}


# ── check output ───────────────────────────────────────────────────────────────
def _source_groups():
    """Source columns (indices in DIPOLE_SOURCE_PATHS order) by group, for the printed table."""
    groups = {"cort E/SST": [], "cort PV": [], "L4 src": [], "VIP": [], "ThalE": [], "ThalPOm": [], "ThalI/POmI": []}
    for i, path in enumerate(DIPOLE_SOURCE_PATHS):
        ct = path[-1]
        if path[0] == "Thalamus":
            groups["ThalE" if ct == "ThalE" else "ThalPOm" if ct == "ThalPOm" else "ThalI/POmI"].append(i)
        elif ct == "VIP":
            groups["VIP"].append(i)
        elif len(path) == 3 and path[1] == "L4":
            groups["L4 src"].append(i)
        elif ct == "PV":
            groups["cort PV"].append(i)
        else:
            groups["cort E/SST"].append(i)
    return groups


def _fmt(vals):
    lo, hi = float(np.min(vals)), float(np.max(vals))
    return f"{lo:g}" if np.isclose(lo, hi) else f"{lo:g}..{hi:g}"


def print_check(p=None):
    p = current_params() if p is None else {**current_params(), **p}
    new_m = lengths_matrix(build_lengths(p))
    ref_m = reference_matrix()
    groups = _source_groups()

    print("parameters")
    for i, (k, u) in enumerate(zip(PARAM_NAMES, PARAM_UNITS), start=1):
        print(f"  {i:2d}  {k:<14s} {p[k]:<8g} {u}")
    print(f"  L4_SOURCE_RATIO = {L4_SOURCE_RATIO}, "
          f"TYPE_RATIO = {TYPE_RATIO}")

    print("\nlength (mm) per target and source group:  parameterized  |  file")
    head = f"  {'target':<9s}" + "".join(f"{g:>20s}" for g in groups)
    print(head)
    for r, name in enumerate(ROW_NAMES):
        cells = []
        for cols in groups.values():
            cells.append(f"{_fmt(new_m[r, cols])} | {_fmt(ref_m[r, cols])}")
        print(f"  {name:<9s}" + "".join(f"{c:>20s}" for c in cells))
    n_diff = int(np.sum(~np.isclose(new_m, ref_m)))
    print(f"\n{n_diff} of {new_m.size} entries differ from the file")

    # checks
    names = [str(x) for x in __import__("parameters").POPULATION_LABELS]
    zero_cols = [i for i, n in enumerate(names) if n.startswith("VIP") or n in ("ThalI", "ThalPOmI")]
    assert np.all(new_m[:, zero_cols] == 0), "VIP / ThalI / ThalPOmI must be 0"
    assert new_m.shape == (9, 34)
    for a, base in (("A1", 1), ("S2", 5)):
        l23, l5 = new_m[base], new_m[base + 2]
        assert np.all(l5 >= F_L5_MIN * l23 - 1e-12), f"{a}: L5 must be >= {F_L5_MIN} x L2/3"
    print("checks passed: VIP/ThalI/ThalPOmI = 0, L5 >= 2 x L2/3, shape (9, 34)")
    check_basis()
    return new_m, ref_m


def check_basis(n_random=5, seed=0):
    """lengths_from_coefficients == lengths_matrix(build_lengths(p)) for the defaults and random valid p."""
    rng = np.random.default_rng(seed)
    cases = [current_params()]
    for _ in range(n_random):
        cases.append({"L_BASE": rng.uniform(.5, 4), "F_L5": rng.uniform(2, 8), "F_L6": rng.uniform(.5, 8),
                      "F_L4_A1": rng.uniform(0, .3), "F_L4_S2": rng.uniform(0, .3),
                      **{k: rng.uniform(0, 4) for k in ("L_THALE_CTX", "L_THALPOM_CTX", "L_THALE_A3B",
                                                        "L_THALPOM_A3B", "L_A3B_CORT")}})
    worst = 0.0
    for q in cases:
        q = {**current_params(), **q}
        ref = lengths_matrix(build_lengths(q))
        lin = lengths_from_coefficients(coefficients_from_params(q), q["F_L5"], q["F_L6"])
        worst = max(worst, float(np.abs(ref - lin).max()))
        back = params_from_coefficients(coefficients_from_params(q), q["F_L5"], q["F_L6"])
        assert all(np.isclose(back[k], q[k]) for k in PARAM_NAMES), "coefficients <-> parameters round trip"
    assert worst < 1e-12, f"coefficient basis differs from build_lengths by {worst:g}"
    print(f"coefficient basis == build_lengths for {len(cases)} parameter sets (max diff {worst:.1e})")


def plot_matrices(new_m, ref_m, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = [str(x) for x in __import__("parameters").POPULATION_LABELS]
    fig, axes = plt.subplots(3, 1, figsize=(13, 10), sharex=True)
    vmax = max(new_m.max(), ref_m.max())
    for ax, m, title, kw in ((axes[0], ref_m, "file (dipole_parameters_flippedPVSST.json)", dict(cmap="viridis", vmin=0, vmax=vmax)),
                             (axes[1], new_m, "parameterized (10 parameters)", dict(cmap="viridis", vmin=0, vmax=vmax)),
                             (axes[2], new_m - ref_m, "parameterized - file",
                              dict(cmap="RdBu_r", vmin=-np.abs(new_m - ref_m).max() or -1, vmax=np.abs(new_m - ref_m).max() or 1))):
        im = ax.imshow(m, aspect="auto", **kw)
        ax.set_yticks(range(9)); ax.set_yticklabels(ROW_NAMES)
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, label="mm")
    axes[2].set_xticks(range(34)); axes[2].set_xticklabels(labels, rotation=90, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", default=None, help="save heatmaps (file, parameterized, difference) to this png")
    ap.add_argument("--write-json", default=None, help="write the parameter file with these lengths to this path")
    args = ap.parse_args()
    new_m, ref_m = print_check()
    if args.plot:
        plot_matrices(new_m, ref_m, args.plot)
        print(f"heatmaps: {args.plot}")
    if args.write_json:
        build_params_json(args.write_json)
        print(f"parameter file: {args.write_json}")
