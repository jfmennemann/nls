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


"""
spacing_x = 0.0
spacing_y = 0.1

gridspec = fig.add_gridspec(ncols=1, nrows=4, left=0.075, right=0.95, bottom=0.1, top=0.925, wspace=spacing_x, hspace=spacing_y, width_ratios=[1], height_ratios=[1, 1, 1, 1])
"""

# plt.tight_layout(0.5)




