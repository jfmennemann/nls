import matplotlib as mpl
import matplotlib.pyplot as plt

export_pdf = True

if export_pdf == True:

    mpl.use("pgf")
    
    # Say, "the default sans-serif font is COMIC SANS"
    # matplotlib.rcParams['font.sans-serif'] = "Comic Sans MS"
    # Then, "ALWAYS use sans-serif fonts"
    # matplotlib.rcParams['font.family'] = "sans-serif"
    
    plt.rcParams.update({
    # "font.family": "serif",             # use serif/main font for text elements
    "font.family": "sans-serif",
    # "font.serif": [],                   # use latex default serif font
    "font.sans-serif": ["DejaVu Sans"],   # use a specific sans-serif font
    "text.usetex": True,                  # use inline math for ticks
    "pgf.rcfonts": False,                 # don't setup fonts from rc parameters
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "pgf.preamble": "\n".join([           # load additional packages
        "\\usepackage{amsmath}",
        "\\usepackage{bbm}",                    
        ])
    })
    
    
color_gridlines_major = '#666666'
color_gridlines_minor = '#999999'

color_stable = 'lightgray'



import numpy as np


N = 100 + 1

x = np.linspace(-6.0, 6.0, N+1, endpoint=True)
y = np.linspace(-4.0, 4.0, N+1, endpoint=True)

x_min = x[0]
y_min = y[0]

x_max = x[-1]
y_max = y[-1]

x_ticks_major = [-6, -4, -2, 0, 2, 4, 6]
x_ticks_minor = [-5, -3, -1, 1, 3, 5]

y_ticks_major = [-4, -2, 0, 2, 4]
y_ticks_minor = [-3, -1, 1, 3]


X, Y = np.meshgrid(x, y)

Z = X + 1j * Y

z_explicit_Euler = np.abs(1+Z)
z_implicit_Euler = np.abs(1/(1-Z))

z_Crank_Nicolson = np.abs((1+0.5*Z)/(1-0.5*Z))
z_RK4            = np.abs(1 + Z + Z**2/2 + Z**3/6+ Z**4/24)




width  = 8.00
height = 5.65

fig = plt.figure("figure_main", figsize=(width, height), facecolor="white", constrained_layout=False)


gridspec = fig.add_gridspec(ncols=2, nrows=2, left=0.075, right=0.95, bottom=0.1, top=0.925, wspace=0.5, hspace=0.5, width_ratios=[1, 1], height_ratios=[1, 1])
# gridspec = fig.add_gridspec(ncols=2, nrows=2, wspace=0.5, hspace=0.5, width_ratios=[1, 1], height_ratios=[1, 1])



ax_explicit_Euler = fig.add_subplot(gridspec[0, 0])
ax_implicit_Euler = fig.add_subplot(gridspec[0, 1])
ax_Crank_Nicolson = fig.add_subplot(gridspec[1, 0])
ax_RK4            = fig.add_subplot(gridspec[1, 1])

"""
x_min = -4
x_max = +4

y_min = -4
y_max = +4
"""

#------------------------------------------------------------------------------
# ax_explicit_Euler = plt.subplot(2,2,1)

levels = [0, 1]

ax_explicit_Euler.set_title("Explicit Euler")

ax_explicit_Euler.contourf(X, Y, z_explicit_Euler, levels, colors=(color_stable,))

ax_explicit_Euler.contour(X, Y, z_explicit_Euler, levels, colors=('black',), linewidths=1.0)

ax_explicit_Euler.axis('equal')

ax_explicit_Euler.set_xlim([x_min, x_max])
ax_explicit_Euler.set_ylim([y_min, y_max])

ax_explicit_Euler.set_xticks(x_ticks_major, minor=False)
ax_explicit_Euler.set_xticks(x_ticks_minor, minor=True)

ax_explicit_Euler.set_yticks(y_ticks_major, minor=False)
ax_explicit_Euler.set_yticks(y_ticks_minor, minor=True)

ax_explicit_Euler.grid(b=True, which='major', color=color_gridlines_major, linestyle='-', linewidth=0.5, alpha=0.5)
ax_explicit_Euler.grid(b=True, which='minor', color=color_gridlines_minor, linestyle='-', linewidth=0.5, alpha=0.5)

ax_explicit_Euler.set_xlabel(r'$\operatorname{Re}(z)$')
ax_explicit_Euler.set_ylabel(r'$\operatorname{Im}(z)$')
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# ax_implicit_Euler = plt.subplot(2,2,2)

levels = [0, 1]

ax_implicit_Euler.set_title("Implicit Euler")

ax_implicit_Euler.contourf(X, Y, z_implicit_Euler, levels, colors=(color_stable ,))

ax_implicit_Euler.contour(X, Y, z_implicit_Euler, levels, colors=('black',), linewidths=1.0)

ax_implicit_Euler.axis('equal')

ax_implicit_Euler.set_xlim([x_min, x_max])
ax_implicit_Euler.set_ylim([y_min, y_max])

ax_implicit_Euler.set_xticks(x_ticks_major, minor=False)
ax_implicit_Euler.set_xticks(x_ticks_minor, minor=True)

ax_implicit_Euler.set_yticks(y_ticks_major, minor=False)
ax_implicit_Euler.set_yticks(y_ticks_minor, minor=True)

ax_implicit_Euler.grid(b=True, which='major', color=color_gridlines_major, linestyle='-', linewidth=0.5, alpha=0.5)
ax_implicit_Euler.grid(b=True, which='minor', color=color_gridlines_minor, linestyle='-', linewidth=0.5, alpha=0.5)

ax_implicit_Euler.set_xlabel(r'$\operatorname{Re}(z)$')
ax_implicit_Euler.set_ylabel(r'$\operatorname{Im}(z)$')
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# ax_Crank_Nicolson = plt.subplot(2,2,3)

levels = [0, 1]

ax_Crank_Nicolson.set_title("Crank-Nicolson")

ax_Crank_Nicolson.contourf(X, Y, z_Crank_Nicolson, levels, colors=(color_stable,))

ax_Crank_Nicolson.contour(X, Y, z_Crank_Nicolson, levels, colors=('black',), linewidths=1.0)

ax_Crank_Nicolson.axis('equal')

ax_Crank_Nicolson.set_xlim([x_min, x_max])
ax_Crank_Nicolson.set_ylim([y_min, y_max])

ax_Crank_Nicolson.set_xticks(x_ticks_major, minor=False)
ax_Crank_Nicolson.set_xticks(x_ticks_minor, minor=True)

ax_Crank_Nicolson.set_yticks(y_ticks_major, minor=False)
ax_Crank_Nicolson.set_yticks(y_ticks_minor, minor=True)

ax_Crank_Nicolson.grid(b=True, which='major', color=color_gridlines_major, linestyle='-', linewidth=0.5, alpha=0.5)
ax_Crank_Nicolson.grid(b=True, which='minor', color=color_gridlines_minor, linestyle='-', linewidth=0.5, alpha=0.5)

ax_Crank_Nicolson.set_xlabel(r'$\operatorname{Re}(z)$')
ax_Crank_Nicolson.set_ylabel(r'$\operatorname{Im}(z)$')
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# ax_RK4 = plt.subplot(2,2,4)

levels = [0, 1]

ax_RK4.set_title("RK4")

ax_RK4.contourf(X, Y, z_RK4, levels, colors=(color_stable,))

ax_RK4.contour(X, Y, z_RK4, levels, colors=('black',), linewidths=1.0)

ax_RK4.axis('equal')

ax_RK4.set_xlim([x_min, x_max])
ax_RK4.set_ylim([y_min, y_max])

ax_RK4.set_xticks(x_ticks_major, minor=False)
ax_RK4.set_xticks(x_ticks_minor, minor=True)

ax_RK4.set_yticks(y_ticks_major, minor=False)
ax_RK4.set_yticks(y_ticks_minor, minor=True)

ax_RK4.grid(b=True, which='major', color=color_gridlines_major, linestyle='-', linewidth=0.5, alpha=0.5)
ax_RK4.grid(b=True, which='minor', color=color_gridlines_minor, linestyle='-', linewidth=0.5, alpha=0.5)

ax_RK4.set_xlabel(r'$\operatorname{Re}(z)$')
ax_RK4.set_ylabel(r'$\operatorname{Im}(z)$')
#------------------------------------------------------------------------------


# plt.tight_layout(0.5)

if export_pdf == True:

    plt.savefig('figure_stability_regions.pdf', backend='pgf')

else:

    plt.show()


