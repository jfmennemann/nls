from style_sheet import mycolors

width_figure_n_colums_1 = 6



# color_u_ref = mycolors.peter_river
# color_u_ref = mycolors.alizarin
# color_u_ref = mycolors.turquoise
# color_u_ref = mycolors.emerald
color_u_ref = mycolors.fluorescent_red

color_u     = mycolors.black


linewidth_u_ref = 1.75
linewidth_u     = 0.75

linestyle_u_ref = '-'
linestyle_u     = '-'




color_gridlines_major = '#d3d3d3'
color_gridlines_minor = '#e5e5e5'

linestyle_gridlines_major = '-'
linestyle_gridlines_minor = '-'

linewidth_gridlines_major = 0.5
linewidth_gridlines_minor = 0.5




color_stability_regions = '#d3d3d3'


rc_parameters = {
    "text.usetex": True,
    "pgf.rcfonts": False,
    #
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "axes.linewidth": 0.85,
    "axes.labelpad": 4,
    "axes.grid": True,
    "axes.grid.axis": "both",
    "axes.grid.which": "both",
    #
    # 'grid.color': 'b0b0b0',
    # 'grid.linestyle': '--',
    # 'grid.linewidth': 0.8,
    # 'grid.alpha': 0.5,
    #
    "xtick.bottom": False,
    "ytick.left": False,
    "xtick.major.pad": 1,
    "ytick.major.pad": 1,
    #
    "legend.frameon": True,
    "legend.framealpha": 1.0,
    "legend.edgecolor": mycolors.wet_asphalt,
    "legend.fontsize": 12,
    "legend.fancybox": False,
    "legend.shadow": True,
    "legend.handlelength": 3.0,
    # 
    "pgf.preamble": ("\\usepackage{amsmath}\n\\usepackage{amsfonts}\n\\usepackage{amssymb}\n\\usepackage{bbm}")
    }





