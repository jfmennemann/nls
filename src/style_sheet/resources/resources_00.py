"""
for tick in ax_00.xaxis.get_major_ticks():
    tick.tick1line.set_visible(False)
    
for tick in ax_00.xaxis.get_minor_ticks():
    tick.tick1line.set_visible(False)
    
for tick in ax_00.yaxis.get_major_ticks():
    tick.tick1line.set_visible(False)
    
for tick in ax_00.yaxis.get_minor_ticks():
    tick.tick1line.set_visible(False)
"""

# ax_30.xaxis.labelpad = 20


"""
# Say, "the default sans-serif font is COMIC SANS"
# matplotlib.rcParams['font.sans-serif'] = "Comic Sans MS"
# Then, "ALWAYS use sans-serif fonts"
# matplotlib.rcParams['font.family'] = "sans-serif"

plt.rcParams.update({
    # "font.family": "serif",             # use serif/main font for text elements
    # "font.family": "sans-serif",
    # "font.serif": [],                   # use latex default serif font
    # "font.sans-serif": ["DejaVu Sans"],   # use a specific sans-serif font
    "text.usetex": True,                  # use inline math for ticks
    "pgf.rcfonts": False,                 # don't setup fonts from rc parameters
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "axes.linewidth": 0.85,
    "xtick.bottom": False,
    "ytick.left": False,
    "xtick.major.pad": 2,
    "ytick.major.pad": 2,
    "axes.labelpad": 4,
    "legend.frameon": True,
    "legend.framealpha": 1.0,
    "legend.edgecolor": mycolors.wet_asphalt,
    "legend.fontsize": 12,
    "legend.fancybox": False,
    "legend.shadow": True,
    "legend.handlelength": 3.0, 
    "pgf.preamble": "\n".join([           # load additional packages
        "\\usepackage{amsmath}",
        "\\usepackage{bbm}",                    
        ])
    })
"""



