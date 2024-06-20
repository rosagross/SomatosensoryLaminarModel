import seaborn as sns
import matplotlib as plt
import numpy as np
import tkinter as tk

def figure_style():
    """
    Set style for plotting figures.
    Modified from Guido Meijer (@guidomeijer in 2024)
    """
    sns.set_theme(style="ticks", context="paper",
            font="Arial",
            rc={"font.size": 7,
                "figure.titlesize": 7,
                "axes.titlesize": 7,
                "axes.labelsize": 7,
                "axes.linewidth": 0.5,
                "lines.linewidth": 1,
                "lines.markersize": 3,
                "xtick.labelsize": 7,
                "ytick.labelsize": 7,
                "savefig.transparent": True,
                "xtick.major.size": 2.5,
                "ytick.major.size": 2.5,
                "xtick.major.width": 0.5,
                "ytick.major.width": 0.5,
                "xtick.minor.size": 2,
                "ytick.minor.size": 2,
                "xtick.minor.width": 0.5,
                "ytick.minor.width": 0.5,
                'legend.fontsize': 7,
                'legend.title_fontsize': 7,
                'legend.frameon': False,
                 })
    
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42


    subject_pal = sns.color_palette(
        np.concatenate((sns.color_palette('tab20'),
                        [plt.colors.to_rgb('maroon'), np.array([0, 0, 0])])))
    
    colors = {'visual' : sns.color_palette('Dark2')[3],
              'somato' : sns.color_palette('Dark2')[2],
              'subject_palette': subject_pal,
              'grey': [0.7, 0.7, 0.7],
              'sert': sns.color_palette('Dark2')[0],
              'wt': [0.7, 0.7, 0.7],
              'awake': sns.color_palette('Dark2')[2],
              'anesthesia': sns.color_palette('Dark2')[3],
              'enhanced': sns.color_palette('colorblind')[3],
              'suppressed': sns.color_palette('colorblind')[0],
              'down-state': sns.color_palette('colorblind')[3],
              'up-state': [1, 1, 1],
              'stim': [0, 0, 0],
              'no-stim': [0.7, 0.7, 0.7],
              'NS': sns.color_palette('Set2')[0],
              'WS': sns.color_palette('Set2')[1],
              'WS1': sns.color_palette('Set2')[1],
              'WS2': sns.color_palette('Set2')[2],
              'states': 'Dark2',
              'states_light': 'Set2',
              'main_states': sns.diverging_palette(20, 210, l=55, center='dark'),
              'OFC': sns.color_palette('Dark2')[0],
              'mPFC': sns.color_palette('Dark2')[1],
              'M2': sns.color_palette('Dark2')[2],
              'Amyg': sns.color_palette('Dark2')[3],
              'Hipp': sns.color_palette('Dark2')[4],
              'VIS': sns.color_palette('Dark2')[5],
              'Pir': sns.color_palette('Dark2')[6],
              'SC': sns.color_palette('Dark2')[7],
              'Thal': sns.color_palette('tab10')[9],
              'PAG': sns.color_palette('Set1')[7],
              'BC': sns.color_palette('Accent')[0],
              'Str': sns.color_palette('Accent')[1],
              'MRN': sns.color_palette('Accent')[2],
              'OLF': sns.color_palette('tab10')[8],
              'RSP': 'r',
              'SNr': [0.75, 0.75, 0.75]}
    screen_width = tk.Tk().winfo_screenwidth()
    dpi = screen_width / 10
    return colors, dpi
