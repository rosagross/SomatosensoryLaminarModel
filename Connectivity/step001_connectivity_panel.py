# %% Connectivity overview panel (publication figure)
# A) A1 connectivity matrix   B) A3b - A1 - S2 inter-area connections
# C) VPM projection matrix    D) POm connections
# The weights are the final connectivity matrix W of Parameter.get_connectivity(), built
# from the values in Simulations/simulation_parameter.json exactly as SomatoModel does.
import os
import sys
import json
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

WDDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) if '__file__' in globals() \
    else os.path.abspath('..')
sys.path.insert(0, os.path.join(WDDIR, 'Simulations'))
from parameters import Parameter, POPULATION_LABELS

try:
    from plotting_style import figure_style
    figure_style()      # sets the seaborn theme before it fails on a headless node (tkinter)
except Exception:
    pass

# %% Weights
A3B, A1, S2, THAL = slice(0, 4), slice(4, 17), slice(17, 30), slice(30, 34)
I_VPM, I_TRN, I_POM, I_POMI = 30, 31, 32, 33

# A1 / S2 population order: E, PV, SST, VIP (L2/3) | E, PV, SST (L4) | E, PV, SST (L5) | E, PV, SST (L6)
LAYERS = ['L2/3', 'L4', 'L5', 'L6']
CELL_OF = ['E', 'PV', 'SST', 'VIP', 'E', 'PV', 'SST', 'E', 'PV', 'SST', 'E', 'PV', 'SST']
LAYER_OF = [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]
CELL_COLOR = {'E': '#4477AA', 'PV': '#EE6677', 'SST': '#228833', 'VIP': '#AA3377'}
AREA_COLOR = {'A3b': '#CCBB44', 'A1': '#4477AA', 'S2': '#228833', 'Thal': '#888888'}
AREA_TINT = {'A3b': '#FBF6DC', 'A1': '#E6EDF5', 'S2': '#E3F1E6', 'Thal': '#EEEEEE'}


def load_weights():
    c = json.load(open(os.path.join(WDDIR, 'Simulations', 'simulation_parameter.json')))
    p = Parameter(c['delay_factor_short'], c['delay_factor'], c['thal_delay_factor'],
                  c['e3b_tau'], c['e1_tau'], c['e2_tau'], c['p_2PVE'], c['p_4PVE'],
                  c['f_intra'], c['f_intra_S2'])
    gE = c['coupling_strength']
    gI = gE * c['strength_I']
    W = p.get_connectivity(c['g_intercortical'], gE, gI,
                           c['g_thal'], c['g_thal'] * c['sI_thal'],
                           c['g_thalPOm'], c['g_thal'] * c['sI_thalPOm'],
                           c['thal_connect'], c['extI_cellcounts'], c['bI_cellcounts'],
                           c['thalE_cellcounts'], c['thalI_cellcounts'], c['pom_cellcounts'],
                           c['pomI_cellcounts'], c['mI_cellcounts'], area=c['area'])
    return W[:, :W.shape[0]]    # drop the background / external / modulatory input columns


def pretty(idx):
    """Short label like 'E L2/3' for a cortical population index."""
    if idx < 4:
        return f"{POPULATION_LABELS[idx][:-2]} BA3b"
    j = (idx - 4) % 13 if idx < 17 else (idx - 17)
    return f"{CELL_OF[j]} {LAYERS[LAYER_OF[j]]}"


# %% Matrix panels
def make_norm(*mats):
    vmax = max(np.nanmax(np.abs(m)) for m in mats)
    return mcolors.SymLogNorm(linthresh=5, linscale=0.5, vmin=-vmax, vmax=vmax)


CMAP = 'RdBu_r'


def matrix_panel(ax, M, xlabels, ylabels, norm, xlabel='Source', ylabel='Target', fs=6, cmap=CMAP):
    im = ax.imshow(M, cmap=cmap, norm=norm, aspect='equal' if M.shape[0] == M.shape[1] else 'auto',
                   interpolation='nearest')
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels(xlabels, rotation=90, fontsize=fs)
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels(ylabels, fontsize=fs)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(length=1.5, pad=1)
    for s in ax.spines.values():
        s.set_linewidth(0.4)
    return im


def add_cbar(fig, im, ax, label='Weight', **kw):
    cb = fig.colorbar(im, ax=ax, **kw)
    vmax = im.norm.vmax
    t = [v for v in (-1000, -100, -10, 0, 10, 100, 1000) if abs(v) <= vmax]
    cb.set_ticks(t); cb.set_ticklabels([str(v) for v in t]); cb.ax.minorticks_off()
    cb.set_label(label, fontsize=7)
    cb.ax.tick_params(labelsize=6, length=2)
    cb.outline.set_linewidth(0.4)
    return cb


# %% Graph panels
class Node:
    def __init__(self, key, xy, cell, label=None):
        self.key, self.xy, self.cell, self.label = key, np.array(xy, float), cell, label


def draw_area_box(ax, name, x0, y0, w, h):
    ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle='round,pad=0.0,rounding_size=0.12',
                                fc=AREA_TINT[name], ec=AREA_COLOR[name], lw=1.0, zorder=0))


def draw_nodes(ax, nodes, r, fs=5):
    for n in nodes.values():
        ax.add_patch(Circle(n.xy, r, fc=CELL_COLOR.get(n.cell, '#BBBBBB'), ec='k', lw=0.4, zorder=3))
        if n.label:
            ax.text(*n.xy, n.label, ha='center', va='center', fontsize=fs, color='w', zorder=4)


def draw_edge(ax, a, b, w, r, rad, color, wmax):
    lw = 0.3 + 2.2 * np.sqrt(abs(w) / wmax)
    style = '-|>,head_length=3.2,head_width=1.5' if w > 0 else '-[,widthB=0.9,lengthB=0.0'
    ax.add_patch(FancyArrowPatch(a, b, connectionstyle=f'arc3,rad={rad}', arrowstyle=style,
                                 lw=lw, color=color, alpha=0.85, shrinkA=r * 30, shrinkB=r * 30,
                                 mutation_scale=1, zorder=2))


def area_grid(nodes, area, x0, y0, dx, dy, pop0, npop, only_E=False):
    """Place the nodes of A1/S2 in a layer (rows, L2/3 on top) x cell class (columns) grid."""
    cols = {'E': 0, 'PV': 1, 'SST': 2, 'VIP': 3}
    for j in range(npop):
        if only_E and CELL_OF[j] != 'E':
            continue
        c = 0 if only_E else cols[CELL_OF[j]]
        nodes[pop0 + j] = Node(pop0 + j, (x0 + c * dx, y0 - LAYER_OF[j] * dy), CELL_OF[j])


def inter_area_edges(W):
    """Nonzero connections between A3b, A1 and S2 (target, source, weight)."""
    areas = {'A3b': range(0, 4), 'A1': range(4, 17), 'S2': range(17, 30)}
    out = []
    for ta, tr in areas.items():
        for sa, sr in areas.items():
            if ta == sa:
                continue
            for t in tr:
                for s in sr:
                    if W[t, s] != 0:
                        out.append((t, s, W[t, s], sa, ta))
    return out


def panel_B(ax, W):
    ax.set_aspect('equal', anchor='N'); ax.axis('off')
    nodes = {}
    r = 0.17
    dx, dy = 0.5, 0.5
    # A3b: collapsed, one row
    for j in range(4):
        nodes[j] = Node(j, (-0.75 + j * dx, 3.3), ['E', 'PV', 'SST', 'VIP'][j])
    area_grid(nodes, 'A1', -3.2, 1.0, dx, dy, 4, 13)
    area_grid(nodes, 'S2', 1.7, 1.0, dx, dy, 17, 13)
    for name, (x0, y0, w, h) in {'A3b': (-1.15, 2.95, 2.3, 0.7), 'A1': (-3.6, -0.75, 2.4, 2.5),
                                 'S2': (1.3, -0.75, 2.4, 2.5)}.items():
        draw_area_box(ax, name, x0, y0, w, h)
    ax.text(0, 3.8, 'BA3b', ha='center', fontsize=8, color=AREA_COLOR['A3b'], weight='bold')
    ax.text(-2.4, 1.95, 'BA1', ha='center', fontsize=8, color=AREA_COLOR['A1'], weight='bold')
    ax.text(2.5, 1.95, 'S2', ha='center', fontsize=8, color=AREA_COLOR['S2'], weight='bold')
    for area, x in (('A1', -3.4), ('S2', 1.5)):
        for l, name in enumerate(LAYERS):
            ax.text(x, 1.0 - l * dy, name, ha='right', fontsize=5, va='center')
    draw_nodes(ax, nodes, r)
    edges = inter_area_edges(W)
    wmax = max(abs(e[2]) for e in edges)
    for t, s, w, sa, ta in edges:
        rad = 0.2
        draw_edge(ax, nodes[s].xy, nodes[t].xy, w, r, rad, AREA_COLOR[sa], wmax)
    ax.set_xlim(-4.0, 4.0); ax.set_ylim(-1.9, 4.0)
    legend_nodes(ax, loc='lower center', wmax=wmax, ncol=3)


def legend_nodes(ax, loc, wmax, ncol=2, extra=None):
    h = [Line2D([], [], marker='o', ls='', mfc=c, mec='k', mew=0.4, ms=4, label=k)
         for k, c in CELL_COLOR.items()]
    for w in (10, 100, 500):
        if w <= wmax:
            h.append(Line2D([], [], color='0.3', lw=0.3 + 2.2 * np.sqrt(w / wmax), label=f'{w}'))
    if extra:
        h += extra
    ax.legend(handles=h, loc=loc, ncol=ncol, fontsize=5.5, frameon=False, handlelength=1.4,
              columnspacing=0.8, borderpad=0, labelspacing=0.3, bbox_to_anchor=(0.5, -0.1))


def panel_D(ax, W):
    ax.set_aspect('equal', anchor='N'); ax.axis('off')
    nodes = {}
    r = 0.17
    dy = 0.55
    # rows: L2/3, L4, L5, L6
    ys = {0: 2.8, 1: 2.25, 2: 1.7, 3: 1.15}
    E_idx = {0: 0, 1: 4, 2: 7, 3: 10}      # offset of E within A1/S2 per layer
    nodes[0] = Node(0, (-3.0, 2.0), 'E')
    for area, x0, pop0 in (('A1', -1.2, 4), ('S2', 1.5, 17)):
        for l in range(4):
            nodes[pop0 + E_idx[l]] = Node(pop0 + E_idx[l], (x0, ys[l]), 'E')
    for name, (x0, y0, w, h) in {'A3b': (-3.4, 1.65, 0.8, 0.7), 'A1': (-1.55, 0.85, 0.7, 2.3),
                                 'S2': (1.15, 0.85, 0.7, 2.3)}.items():
        draw_area_box(ax, name, x0, y0, w, h)
    ax.text(-3.0, 2.5, 'BA3b', ha='center', fontsize=8, color=AREA_COLOR['A3b'], weight='bold')
    ax.text(-1.2, 3.3, 'BA1', ha='center', fontsize=8, color=AREA_COLOR['A1'], weight='bold')
    ax.text(1.5, 3.3, 'S2', ha='center', fontsize=8, color=AREA_COLOR['S2'], weight='bold')
    for x in (-1.55, 1.15):
        for l, name in enumerate(LAYERS):
            ax.text(x - 0.05, ys[l], name, ha='right', va='center', fontsize=5)
    # thalamic nodes
    nodes[I_POM] = Node(I_POM, (0.15, -0.4), 'E', 'POm')
    nodes[I_POMI] = Node(I_POMI, (0.15, -1.55), 'PV', 'TRN')
    rt = 0.3
    ax.add_patch(FancyBboxPatch((-0.45, -2.0), 1.2, 2.1, boxstyle='round,pad=0.0,rounding_size=0.12',
                                fc=AREA_TINT['Thal'], ec=AREA_COLOR['Thal'], lw=1.0, zorder=0))
    ax.text(0.9, -1.0, 'Thalamus', ha='left', rotation=90, va='center', fontsize=8, color='0.35', weight='bold')
    draw_nodes(ax, {k: nodes[k] for k in nodes if k < 4 or 4 <= k < 34 and k not in (I_POM, I_POMI)}, r)
    for k in (I_POM, I_POMI):
        n = nodes[k]
        ax.add_patch(Circle(n.xy, rt, fc=CELL_COLOR[n.cell], ec='k', lw=0.4, zorder=3))
        ax.text(*n.xy, n.label, ha='center', va='center', fontsize=5.5, color='w', zorder=4)
    # edges
    cort = [k for k in nodes if k < 30]
    wmax = max(abs(W[k, I_POM]) for k in cort)
    wmax = max(wmax, max(abs(W[I_POM, k]) for k in cort), abs(W[I_POMI, I_POM]), abs(W[I_POM, I_POMI]))
    area_of = lambda k: 'A3b' if k < 4 else ('A1' if k < 17 else 'S2')
    for k in cort:
        for pre, post, rad in ((I_POM, k, 0.12), (k, I_POM, 0.12)):
            w = W[post, pre]
            if w == 0:
                continue
            col = AREA_COLOR['Thal'] if pre == I_POM else AREA_COLOR[area_of(k)]
            rp = rt if pre == I_POM else r
            rq = r if pre == I_POM else rt
            lw = 0.3 + 2.2 * np.sqrt(abs(w) / wmax)
            ax.add_patch(FancyArrowPatch(nodes[pre].xy, nodes[post].xy, connectionstyle=f'arc3,rad={rad}',
                                         arrowstyle='-|>,head_length=3.2,head_width=1.5', lw=lw, color=col,
                                         alpha=0.85, shrinkA=rp * 30, shrinkB=rq * 30,
                                         mutation_scale=1, zorder=2))
    for pre, post, rad in ((I_POM, I_POMI, -0.3), (I_POMI, I_POM, -0.3)):
        w = W[post, pre]
        lw = 0.3 + 2.2 * np.sqrt(abs(w) / wmax)
        style = '-|>,head_length=3.2,head_width=1.5' if w > 0 else '-[,widthB=0.9,lengthB=0.0'
        ax.add_patch(FancyArrowPatch(nodes[pre].xy, nodes[post].xy, connectionstyle=f'arc3,rad={rad}',
                                     arrowstyle=style, lw=lw, color='0.35' if w > 0 else CELL_COLOR['PV'],
                                     shrinkA=rt * 30, shrinkB=rt * 30, mutation_scale=1, zorder=2))
    ax.set_xlim(-4.0, 2.4); ax.set_ylim(-2.6, 3.6)
    legend_nodes(ax, loc='lower center', wmax=wmax, ncol=3,
                 extra=[Line2D([], [], color=CELL_COLOR['PV'], lw=1.2, label='inhibitory')])


# %% Figure
def area_labels(ax, y=1.0, dy=0.015, fs=7, left=False):
    """Area names above (and optionally left of) a 30 x 30 cortical matrix."""
    for x, name in ((1.5, 'A3b'), (10, 'A1'), (23, 'S2')):
        txt = {'A3b': 'BA3b', 'A1': 'BA1'}.get(name, name)
        ax.text(x, -0.9, txt, ha='center', va='bottom', fontsize=fs, weight='bold',
                color=AREA_COLOR[name], transform=ax.transData, clip_on=False)
        if left:
            ax.text(30.2, x, txt, ha='left', va='center', fontsize=fs, weight='bold', rotation=-90,
                    color=AREA_COLOR[name], transform=ax.transData, clip_on=False)


def make_figure(W, outdir):
    cort = list(range(30))
    clab = [pretty(i) for i in cort]
    # all matrices are built first so that they share one colour scale
    A = W[A1, A1]
    M_C = np.vstack((W[:30, I_VPM], W[I_VPM, :30], W[I_TRN, :30]))
    L_C = W[I_VPM:I_TRN + 1, I_VPM:I_TRN + 1]
    M_B = W[:30, :30].astype(float).copy()          # inter-area connections only
    for sl in (A3B, A1, S2):
        M_B[sl, sl] = np.nan
    M_D = np.column_stack((W[:30, I_POM], W[I_POM, :30], W[I_POMI, :30]))     # rows: ctx
    L_D = W[I_POM:I_POMI + 1, I_POM:I_POMI + 1]
    norm = make_norm(A, M_C, L_C, M_B, M_D, L_D)

    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 3, left=0.05, right=0.93, top=0.95, bottom=0.07,
                          width_ratios=[1, 1.1, 1.05], wspace=0.22, hspace=0.30)

    # A) A1 connectivity matrix
    axA = fig.add_subplot(gs[0, 0])
    imA = matrix_panel(axA, A, clab[4:17], clab[4:17], norm)
    for b in (3.5, 6.5, 9.5):
        axA.axhline(b, color='k', lw=0.3); axA.axvline(b, color='k', lw=0.3)
    axA.set_title('BA1', fontsize=8, pad=3)

    # C) VPM projection matrix
    gsC = gs[1, 0].subgridspec(2, 2, height_ratios=[3, 2], width_ratios=[1, 1], hspace=1.2, wspace=0.6)
    axC = fig.add_subplot(gsC[0, :])
    matrix_panel(axC, M_C, clab, ['VPM → ctx', 'ctx → VPM', 'ctx → TRN'], norm, xlabel='', ylabel='', fs=4)
    for b in (3.5, 16.5):
        axC.axvline(b, color='k', lw=0.6)
    area_labels(axC, fs=6)
    axL = fig.add_subplot(gsC[1, 0])
    matrix_panel(axL, L_C, ['VPM', 'TRN'], ['VPM', 'TRN'], norm, fs=6)
    annotate(axL, L_C)
    axL.set_title('Thalamic loop', fontsize=7, pad=2)

    # B) inter-area graph + matrix
    axB = fig.add_subplot(gs[0, 1])
    panel_B(axB, W)
    axBm = fig.add_subplot(gs[0, 2])
    cmap = plt.get_cmap(CMAP).copy(); cmap.set_bad('0.92')
    matrix_panel(axBm, M_B, clab, clab, norm, fs=4, cmap=cmap)
    for b in (3.5, 16.5):
        axBm.axhline(b, color='k', lw=0.4); axBm.axvline(b, color='k', lw=0.4)
    area_labels(axBm, fs=6, left=True)

    # D) POm graph + matrix
    axD = fig.add_subplot(gs[1, 1])
    gsD = gs[1, 2].subgridspec(1, 2, width_ratios=[0.8, 1], wspace=0.9)
    panel_D(axD, W)
    axDm = fig.add_subplot(gsD[0, 0])
    matrix_panel(axDm, M_D, ['POm → ctx', 'ctx → POm', 'ctx → TRN'], clab, norm, xlabel='', ylabel='', fs=4)
    axDm.set_aspect('auto')
    for b in (3.5, 16.5):
        axDm.axhline(b, color='k', lw=0.4)
    axDm.tick_params(axis='x', labelsize=5)
    gsDl = gsD[0, 1].subgridspec(3, 1, height_ratios=[1, 1.2, 1], hspace=0.3)
    axDl = fig.add_subplot(gsDl[1])
    matrix_panel(axDl, L_D, ['POm', 'TRN'], ['POm', 'TRN'], norm, fs=6)
    annotate(axDl, L_D)
    axDl.set_title('POm - TRN loop', fontsize=7, pad=2)

    # one shared colour bar
    cax = fig.add_axes([0.955, 0.07, 0.01, 0.88])
    cb = add_cbar(fig, imA, None, cax=cax)

    # panel letters in figure coordinates
    for letter, cell in (('A', gs[0, 0]), ('B', gs[0, 1]), ('C', gs[1, 0]), ('D', gs[1, 1])):
        pos = cell.get_position(fig)
        fig.text(pos.x0 - (0.035 if letter in 'AC' else 0.01), pos.y1 + 0.005, letter, fontsize=12, weight='bold', va='center')
    os.makedirs(outdir, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(outdir, f'connectivity_panel.{ext}'), dpi=300)
    # copy for the manuscript
    img_dir = os.path.join(WDDIR, 'Writing', 'img')
    os.makedirs(img_dir, exist_ok=True)
    fig.savefig(os.path.join(img_dir, 'Fig01_connectivity.png'), dpi=300)
    return fig


def annotate(ax, M, fs=6):
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f'{M[i, j] + 0:.0f}'.replace('-0', '0'), ha='center', va='center', fontsize=fs,
                    color='w' if abs(M[i, j]) > 400 else 'k')


# %%
if __name__ == '__main__' or 'get_ipython' in globals():
    W = load_weights()
    print('W (populations only):', W.shape)
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    outdir = os.path.join(WDDIR, 'Connectivity', 'Figures', stamp)
    make_figure(W, outdir)
    print('saved to', outdir)
