from scipy.sparse import diags, eye

from scipy.sparse.linalg import spsolve

import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)

from linear_schroedinger.wave_packet.reference_solutions import gaussian

from linear_schroedinger.wave_packet.figure_1 import Figure1



order_spatial_discretization = 2


x0 = 0.0

sigma_0 = 0.5

k0 = 4



x_min = -8
x_max = +8

Jx = 400

x = np.linspace(x_min, x_max, Jx+1, endpoint=True)

dx = x[1] - x[0]



T = 2

dt = 0.0025

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)




if order_spatial_discretization == 2:
    
    D_xx = diags([1, -2, 1], [-1, 0, 1], shape=(Jx-1, Jx-1))
    
    D_xx = D_xx / dx**2

""" 
if order_spatial_discretization == 4:
    
    D_xx = diags([-1, 16, -30, 16, -1], [-2, -1, 0, 1, 2], shape=(Jx-1, Jx-1))

    D_xx = D_xx / (12 * dx**2)
    
if order_spatial_discretization == 6:
    
    D_xx = diags([2, -27, 270, -490, 270, -27, 2], [-3, -2, -1, 0, 1, 2, 3], shape=(Jx-1, Jx-1))
    
    D_xx = D_xx / (180 * dx**2)
    
if order_spatial_discretization == 8:
    
    D_xx = diags([-9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9], [-4, -3, -2, -1, 0, 1, 2, 3, 4], shape=(Jx-1, Jx-1))
    
    D_xx = D_xx / (5040 * dx**2)
"""


E = eye(Jx-1)


A = E - 0.25 * 1j * dt * D_xx
B = E + 0.25 * 1j * dt * D_xx






u_ref = gaussian(x, 0.0, x0, sigma_0, k0)

u = u_ref[1:-1]

assert(u.size == Jx-1)


u_complete = np.zeros_like(u_ref)

u_complete[1:-1] = u[:]





n_mod_times_analysis = 25

times_analysis = times[::n_mod_times_analysis]


rel_error_of_times_analysis = np.zeros_like(times_analysis)




fig_1 = Figure1(x, times, screen_size='large')

fig_1.update_u(u_complete, u_ref)

fig_1.redraw()





n_0 = 0
n_1 = 200
n_2 = 400
n_3 = 600

assert(n_0 % n_mod_times_analysis == 0)
assert(n_1 % n_mod_times_analysis == 0)
assert(n_2 % n_mod_times_analysis == 0)
assert(n_3 % n_mod_times_analysis == 0)


t_0 = times[n_0]
t_1 = times[n_1]
t_2 = times[n_2]
t_3 = times[n_3]


print(t_0)
print(t_1)
print(t_2)
print(t_3)





nr_times_analysis = 0

for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        u_ref = gaussian(x, t, x0, sigma_0, k0)
    
        u_complete[1:-1] = u[:]
        
        rel_error_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_complete-u_ref) / np.linalg.norm(u_complete)
        times_analysis[nr_times_analysis] = t 
        
        
        fig_1.update_u(u_complete, u_ref)
        fig_1.update_rel_error(rel_error_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.redraw()
        
        
        nr_times_analysis = nr_times_analysis + 1
        
        
        if n == n_0:
            
            u_ref_0 = u_ref.copy()
            u_0 = u_complete.copy()
            
        if n == n_1:
            
            u_ref_1 = u_ref.copy()
            u_1 = u_complete.copy()
            
        if n == n_2:
            
            u_ref_2 = u_ref.copy()
            u_2 = u_complete.copy()
            
        if n == n_3:
            
            u_ref_3 = u_ref.copy()
            u_3 = u_complete.copy()
            
            
        
    
    u = spsolve(A, B*u)
    



print()
print(t_0)
print(t_1)
print(t_2)
print(t_3)
print()









from style_sheet import mystyle

import matplotlib as mpl
import matplotlib.pyplot as plt

export_pdf = True

if export_pdf == True:

    mpl.use("pgf")
    
    plt.rcParams.update(mystyle.rc_parameters)
    


linestyle_u_ref = mystyle.linestyle_u_ref
linestyle_u     = mystyle.linestyle_u

linewidth_u_ref = mystyle.linewidth_u_ref
linewidth_u     = mystyle.linewidth_u

color_u_ref = mystyle.color_u_ref
color_u     = mystyle.color_u

label_u_ref = r'$u_\mathrm{ref}$'
label_u     = r'$u$'

color_gridlines_major = mystyle.color_gridlines_major
color_gridlines_minor = mystyle.color_gridlines_minor

linestyle_gridlines_major = mystyle.linestyle_gridlines_major
linestyle_gridlines_minor = mystyle.linestyle_gridlines_minor

linewidth_gridlines_major = mystyle.linewidth_gridlines_major
linewidth_gridlines_minor = mystyle.linewidth_gridlines_minor

x_ticks_major = np.array([-8, -4, 0, 4, 8])
x_ticks_minor = np.array([-7, -6, -5, -3, -2, -1, 1, 2, 3, 5, 6, 7])
            
y_ticks_major = np.array([0, 0.5, 1.0])
y_ticks_minor = np.array([0.25, 0.75])

x_min = -8.5
x_max = +8.5

y_min = -0.2
y_max = +1.2


xlabel = r'$x$'

ylabel_00 = r'$|u(x,t_0)|^2$'
ylabel_10 = r'$|u(x,t_1)|^2$'
ylabel_20 = r'$|u(x,t_2)|^2$'
ylabel_30 = r'$|u(x,t_3)|^2$'



width  = 6.5
height = 6.5

fig = plt.figure("figure_wave_packet_dirichlet", figsize=(width, height), facecolor="white", constrained_layout=False)

spacing_x = 0.0
spacing_y = 0.065

gridspec = fig.add_gridspec(ncols=1, nrows=4, left=0.085, right=0.99, bottom=0.06, top=0.99, wspace=spacing_x, hspace=spacing_y, width_ratios=[1], height_ratios=[1, 1, 1, 1])



#==========================================================================================
ax_00 = fig.add_subplot(gridspec[0, 0])
ax_10 = fig.add_subplot(gridspec[1, 0])
ax_20 = fig.add_subplot(gridspec[2, 0])
ax_30 = fig.add_subplot(gridspec[3, 0])
#==========================================================================================

#==========================================================================================
ax_00.plot(x, np.abs(u_ref_0)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_00.plot(x, np.abs(u_0)**2,     linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_00.set_xlim(x_min, x_max)
ax_00.set_ylim(y_min, y_max)

ax_00.set_xticks(x_ticks_major, minor=False)
ax_00.set_xticks(x_ticks_minor, minor=True)

ax_00.set_yticks(y_ticks_major, minor=False)
ax_00.set_yticks(y_ticks_minor, minor=True)

ax_00.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_00.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_00.set_ylabel(ylabel_00)

ax_00.set_xticklabels([])

ax_00.legend(loc='upper left', ncol=2)
#==========================================================================================

#==========================================================================================
ax_10.plot(x, np.abs(u_ref_1)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_10.plot(x, np.abs(u_1)**2,     linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_10.set_xlim(x_min, x_max)
ax_10.set_ylim(y_min, y_max)

ax_10.set_xticks(x_ticks_major, minor=False)
ax_10.set_xticks(x_ticks_minor, minor=True)

ax_10.set_yticks(y_ticks_major, minor=False)
ax_10.set_yticks(y_ticks_minor, minor=True)

ax_10.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_10.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_10.set_ylabel(ylabel_10)
    
ax_10.set_xticklabels([])
#==========================================================================================

#==========================================================================================
ax_20.plot(x, np.abs(u_ref_2)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_20.plot(x, np.abs(u_2)**2,     linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_20.set_xlim(x_min, x_max)
ax_20.set_ylim(y_min, y_max)

ax_20.set_xticks(x_ticks_major, minor=False)
ax_20.set_xticks(x_ticks_minor, minor=True)

ax_20.set_yticks(y_ticks_major, minor=False)
ax_20.set_yticks(y_ticks_minor, minor=True)

ax_20.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_20.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_20.set_ylabel(ylabel_20)
    
ax_20.set_xticklabels([])
#==========================================================================================

#==========================================================================================
ax_30.plot(x, np.abs(u_ref_3)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_30.plot(x, np.abs(u_3)**2,     linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_30.set_xlim(x_min, x_max)
ax_30.set_ylim(y_min, y_max)

ax_30.set_xticks(x_ticks_major, minor=False)
ax_30.set_xticks(x_ticks_minor, minor=True)

ax_30.set_yticks(y_ticks_major, minor=False)
ax_30.set_yticks(y_ticks_minor, minor=True)

ax_30.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_30.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_30.set_xlabel(xlabel)
ax_30.set_ylabel(ylabel_30)
#==========================================================================================



if export_pdf == True:
    
    path = "/home/jfmennemann/git/nls/pdf/"

    filepath = path + "figure_wave_packet_dirichlet.pdf"

    plt.savefig(filepath, backend='pgf')

else:

    plt.show()



































