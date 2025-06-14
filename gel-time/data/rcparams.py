import matplotlib as mpl

# inspired by magma
named_colors = {
    "black": "k",
    "orange": "#f89540",
    "pink": "#cc4778",
    "purple": "#7e03a8",
    "blue": "#0d0887",
    "yellow": "#f0f921",
    "green": "#00b35a",
    "red": "#e63946",
    "cyan": "#00cfc1",
    "gray": "#7f7f7f",
}

custom_rcparams = {
    'lines.linewidth': 2,
    'axes.labelsize': 18,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 18,
    'figure.figsize': (8, 8),
    'figure.titlesize': 18,
    'savefig.bbox': 'tight',
    'axes.grid': True,
    'legend.frameon': True,
    'axes.facecolor': 'white',      
    'axes.labelweight': 'bold',       
    'xtick.direction': 'inout',      
    'ytick.direction': 'inout',
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'xtick.major.width': 2,
    'ytick.major.width': 2,
    'xtick.minor.width': 1.5,
    'ytick.minor.width': 1.5,
}

# Apply custom rcParams on import
mpl.rcParams.update(custom_rcparams)