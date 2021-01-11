from scipy.sparse import eye

from scipy.sparse.linalg import spsolve

import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)

from linear_schroedinger.wave_packet.reference_solution import gaussian

from linear_schroedinger.wave_packet.figure_1 import Figure1


from differentiation import finite_differences_1d



x0 = 0.0

sigma_0 = 0.5

k0 = 8



x_min = -8
x_max = +8

L = x_max - x_min

Jx = 200

x_complete = np.linspace(x_min, x_max, Jx+1, endpoint=True)

x = x_complete[0:-1]

dx = x[1] - x[0]



T = 2


dt = 0.0025

n_mod_times_analysis = 25
# n_mod_times_analysis = 1



n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)





E = eye(Jx)

D2 = finite_differences_1d.get_D2_circulant_2nd_order(Jx, dx)
D4 = finite_differences_1d.get_D2_circulant_4th_order(Jx, dx)
D6 = finite_differences_1d.get_D2_circulant_6th_order(Jx, dx)
D8 = finite_differences_1d.get_D2_circulant_8th_order(Jx, dx)

A_cn_2 = E - 0.25 * 1j * dt * D2
A_cn_4 = E - 0.25 * 1j * dt * D4
A_cn_6 = E - 0.25 * 1j * dt * D6
A_cn_8 = E - 0.25 * 1j * dt * D8

B_cn_2 = E + 0.25 * 1j * dt * D2
B_cn_4 = E + 0.25 * 1j * dt * D4
B_cn_6 = E + 0.25 * 1j * dt * D6
B_cn_8 = E + 0.25 * 1j * dt * D8


A_rk4_2 = 1j * 0.5 * D2
A_rk4_4 = 1j * 0.5 * D4
A_rk4_6 = 1j * 0.5 * D6
A_rk4_8 = 1j * 0.5 * D8







u_ref_0 = gaussian(x, 0.0, x0, sigma_0, k0)

u_cn_2 = u_ref_0.copy()
u_cn_4 = u_ref_0.copy()
u_cn_6 = u_ref_0.copy()
u_cn_8 = u_ref_0.copy()

u_rk4_2 = u_ref_0.copy()
u_rk4_4 = u_ref_0.copy()
u_rk4_6 = u_ref_0.copy()
u_rk4_8 = u_ref_0.copy()


assert(u_cn_2.size == Jx)
assert(u_cn_4.size == Jx)
assert(u_cn_6.size == Jx)
assert(u_cn_8.size == Jx)

assert(u_rk4_2.size == Jx)
assert(u_rk4_4.size == Jx)
assert(u_rk4_6.size == Jx)
assert(u_rk4_8.size == Jx)








times_analysis = times[::n_mod_times_analysis]


rel_error_cn_2_of_times_analysis = np.zeros_like(times_analysis)
rel_error_cn_4_of_times_analysis = np.zeros_like(times_analysis)
rel_error_cn_6_of_times_analysis = np.zeros_like(times_analysis)
rel_error_cn_8_of_times_analysis = np.zeros_like(times_analysis)

rel_error_rk4_2_of_times_analysis = np.zeros_like(times_analysis)
rel_error_rk4_4_of_times_analysis = np.zeros_like(times_analysis)
rel_error_rk4_6_of_times_analysis = np.zeros_like(times_analysis)
rel_error_rk4_8_of_times_analysis = np.zeros_like(times_analysis)


deviation_mass_cn_2_of_times_analysis = np.zeros_like(times_analysis)
deviation_mass_cn_4_of_times_analysis = np.zeros_like(times_analysis)
deviation_mass_cn_6_of_times_analysis = np.zeros_like(times_analysis)
deviation_mass_cn_8_of_times_analysis = np.zeros_like(times_analysis)

deviation_mass_rk4_2_of_times_analysis = np.zeros_like(times_analysis)
deviation_mass_rk4_4_of_times_analysis = np.zeros_like(times_analysis)
deviation_mass_rk4_6_of_times_analysis = np.zeros_like(times_analysis)
deviation_mass_rk4_8_of_times_analysis = np.zeros_like(times_analysis)





fig_1 = Figure1(x, times, screen_size='large')

fig_1.update_u(u_rk4_8, u_ref_0)

fig_1.redraw()



t_0 = 0
t_1 = 0.25
t_2 = 0.5
t_3 = 0.75
t_4 = 1.0
t_5 = 1.25
t_6 = 1.5
t_7 = 1.75
t_8 = 2.0


indices_tmp = (np.abs(times - t_0) < 1e-14) 
n_0 = np.where(indices_tmp)[0][0]

indices_tmp = (np.abs(times - t_1) < 1e-14) 
n_1 = np.where(indices_tmp)[0][0]

indices_tmp = (np.abs(times - t_2) < 1e-14) 
n_2 = np.where(indices_tmp)[0][0]

indices_tmp = (np.abs(times - t_3) < 1e-14) 
n_3 = np.where(indices_tmp)[0][0]

indices_tmp = (np.abs(times - t_4) < 1e-14) 
n_4 = np.where(indices_tmp)[0][0]

indices_tmp = (np.abs(times - t_5) < 1e-14) 
n_5 = np.where(indices_tmp)[0][0]

indices_tmp = (np.abs(times - t_6) < 1e-14) 
n_6 = np.where(indices_tmp)[0][0]

indices_tmp = (np.abs(times - t_7) < 1e-14) 
n_7 = np.where(indices_tmp)[0][0]

indices_tmp = (np.abs(times - t_8) < 1e-14) 
n_8 = np.where(indices_tmp)[0][0]

print(times[n_0])
print(times[n_1])
print(times[n_2])
print(times[n_3])
print(times[n_4])
print(times[n_5])
print(times[n_6])
print(times[n_7])
print(times[n_8])

assert(n_0 % n_mod_times_analysis == 0)
assert(n_1 % n_mod_times_analysis == 0)
assert(n_2 % n_mod_times_analysis == 0)
assert(n_3 % n_mod_times_analysis == 0)
assert(n_4 % n_mod_times_analysis == 0)
assert(n_5 % n_mod_times_analysis == 0)
assert(n_6 % n_mod_times_analysis == 0)
assert(n_7 % n_mod_times_analysis == 0)
assert(n_8 % n_mod_times_analysis == 0)




nr_times_analysis = 0

for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        u_ref = (
                + gaussian(x - 10*L, t, x0, sigma_0, k0) 
                + gaussian(x -  9*L, t, x0, sigma_0, k0) 
                + gaussian(x -  8*L, t, x0, sigma_0, k0) 
                + gaussian(x -  7*L, t, x0, sigma_0, k0) 
                + gaussian(x -  6*L, t, x0, sigma_0, k0) 
                + gaussian(x -  5*L, t, x0, sigma_0, k0) 
                + gaussian(x -  4*L, t, x0, sigma_0, k0) 
                + gaussian(x -  3*L, t, x0, sigma_0, k0) 
                + gaussian(x -  2*L, t, x0, sigma_0, k0) 
                + gaussian(x -  1*L, t, x0, sigma_0, k0)
                + gaussian(x +  0*L, t, x0, sigma_0, k0)
                + gaussian(x +  1*L, t, x0, sigma_0, k0)
                + gaussian(x +  2*L, t, x0, sigma_0, k0)
                + gaussian(x +  3*L, t, x0, sigma_0, k0)
                + gaussian(x +  4*L, t, x0, sigma_0, k0)
                + gaussian(x +  5*L, t, x0, sigma_0, k0)
                + gaussian(x +  6*L, t, x0, sigma_0, k0)
                + gaussian(x +  7*L, t, x0, sigma_0, k0)
                + gaussian(x +  8*L, t, x0, sigma_0, k0)
                + gaussian(x +  9*L, t, x0, sigma_0, k0)
                + gaussian(x + 10*L, t, x0, sigma_0, k0)
                )
        
        rel_error_cn_2_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_cn_2-u_ref) / np.linalg.norm(u_ref)
        rel_error_cn_4_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_cn_4-u_ref) / np.linalg.norm(u_ref)
        rel_error_cn_6_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_cn_6-u_ref) / np.linalg.norm(u_ref)
        rel_error_cn_8_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_cn_8-u_ref) / np.linalg.norm(u_ref)
        
        rel_error_rk4_2_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_rk4_2-u_ref) / np.linalg.norm(u_ref)
        rel_error_rk4_4_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_rk4_4-u_ref) / np.linalg.norm(u_ref)
        rel_error_rk4_6_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_rk4_6-u_ref) / np.linalg.norm(u_ref)
        rel_error_rk4_8_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_rk4_8-u_ref) / np.linalg.norm(u_ref)
        
        
        deviation_mass_cn_2_of_times_analysis[nr_times_analysis] = np.abs( 1 - (np.linalg.norm(u_cn_2)/np.linalg.norm(u_ref_0))**2 )
        deviation_mass_cn_4_of_times_analysis[nr_times_analysis] = np.abs( 1 - (np.linalg.norm(u_cn_4)/np.linalg.norm(u_ref_0))**2 )
        deviation_mass_cn_6_of_times_analysis[nr_times_analysis] = np.abs( 1 - (np.linalg.norm(u_cn_6)/np.linalg.norm(u_ref_0))**2 )
        deviation_mass_cn_8_of_times_analysis[nr_times_analysis] = np.abs( 1 - (np.linalg.norm(u_cn_8)/np.linalg.norm(u_ref_0))**2 )
        
        deviation_mass_rk4_2_of_times_analysis[nr_times_analysis] = np.abs( 1 - (np.linalg.norm(u_rk4_2)/np.linalg.norm(u_ref_0))**2 )
        deviation_mass_rk4_4_of_times_analysis[nr_times_analysis] = np.abs( 1 - (np.linalg.norm(u_rk4_4)/np.linalg.norm(u_ref_0))**2 )
        deviation_mass_rk4_6_of_times_analysis[nr_times_analysis] = np.abs( 1 - (np.linalg.norm(u_rk4_6)/np.linalg.norm(u_ref_0))**2 )
        deviation_mass_rk4_8_of_times_analysis[nr_times_analysis] = np.abs( 1 - (np.linalg.norm(u_rk4_8)/np.linalg.norm(u_ref_0))**2 )
        
        
        
        
        times_analysis[nr_times_analysis] = t 
        
        
        fig_1.update_u(u_rk4_8, u_ref)
        
        
        fig_1.redraw()
        
        
        nr_times_analysis = nr_times_analysis + 1
        
        
        u_snapshot = u_cn_4.copy()
        
        if n == n_0:
            
            u_ref_0 = u_ref.copy()
            u_0 = u_snapshot.copy()
            
        if n == n_1:
            
            u_ref_1 = u_ref.copy()
            u_1 = u_snapshot.copy()
            
        if n == n_2:
            
            u_ref_2 = u_ref.copy()
            u_2 = u_snapshot.copy()
            
        if n == n_3:
            
            u_ref_3 = u_ref.copy()
            u_3 = u_snapshot.copy()
            
        if n == n_4:
            
            u_ref_4 = u_ref.copy()
            u_4 = u_snapshot.copy()
            
        if n == n_5:
            
            u_ref_5 = u_ref.copy()
            u_5 = u_snapshot.copy()
            
        if n == n_6:
            
            u_ref_6 = u_ref.copy()
            u_6 = u_snapshot.copy()
            
        if n == n_7:
            
            u_ref_7 = u_ref.copy()
            u_7 = u_snapshot.copy()
            
        if n == n_8:
            
            u_ref_8 = u_ref.copy()
            u_8 = u_snapshot.copy()
            
            
        
    
    u_cn_2 = spsolve(A_cn_2, B_cn_2 * u_cn_2)
    u_cn_4 = spsolve(A_cn_4, B_cn_4 * u_cn_4)
    u_cn_6 = spsolve(A_cn_6, B_cn_6 * u_cn_6)
    u_cn_8 = spsolve(A_cn_8, B_cn_8 * u_cn_8)
    
    
    k1_rk4_2 = A_rk4_2 * (u_rk4_2                      )
    k2_rk4_2 = A_rk4_2 * (u_rk4_2 + 0.5 * dt * k1_rk4_2)
    k3_rk4_2 = A_rk4_2 * (u_rk4_2 + 0.5 * dt * k2_rk4_2)
    k4_rk4_2 = A_rk4_2 * (u_rk4_2 + 1.0 * dt * k3_rk4_2)
    
    u_rk4_2 = u_rk4_2 + (dt/6.0) * (k1_rk4_2 + 2*k2_rk4_2 + 2*k3_rk4_2 + k4_rk4_2)
    
    
    k1_rk4_4 = A_rk4_4 * (u_rk4_4                      )
    k2_rk4_4 = A_rk4_4 * (u_rk4_4 + 0.5 * dt * k1_rk4_4)
    k3_rk4_4 = A_rk4_4 * (u_rk4_4 + 0.5 * dt * k2_rk4_4)
    k4_rk4_4 = A_rk4_4 * (u_rk4_4 + 1.0 * dt * k3_rk4_4)
    
    u_rk4_4 = u_rk4_4 + (dt/6.0) * (k1_rk4_4 + 2*k2_rk4_4 + 2*k3_rk4_4 + k4_rk4_4)
    
    
    k1_rk4_6 = A_rk4_6 * (u_rk4_6                      )
    k2_rk4_6 = A_rk4_6 * (u_rk4_6 + 0.5 * dt * k1_rk4_6)
    k3_rk4_6 = A_rk4_6 * (u_rk4_6 + 0.5 * dt * k2_rk4_6)
    k4_rk4_6 = A_rk4_6 * (u_rk4_6 + 1.0 * dt * k3_rk4_6)
    
    u_rk4_6 = u_rk4_6 + (dt/6.0) * (k1_rk4_6 + 2*k2_rk4_6 + 2*k3_rk4_6 + k4_rk4_6)
    
    
    k1_rk4_8 = A_rk4_8 * (u_rk4_8                      )
    k2_rk4_8 = A_rk4_8 * (u_rk4_8 + 0.5 * dt * k1_rk4_8)
    k3_rk4_8 = A_rk4_8 * (u_rk4_8 + 0.5 * dt * k2_rk4_8)
    k4_rk4_8 = A_rk4_8 * (u_rk4_8 + 1.0 * dt * k3_rk4_8)
    
    u_rk4_8 = u_rk4_8 + (dt/6.0) * (k1_rk4_8 + 2*k2_rk4_8 + 2*k3_rk4_8 + k4_rk4_8)
    



print()
print(t_0)
print(t_1)
print(t_2)
print(t_3)
print()









from style_sheet import mystyle

import matplotlib as mpl

import matplotlib.pyplot as plt

from matplotlib.ticker import FixedLocator, NullFormatter


export_pdf = True

if export_pdf == True:

    mpl.use("pgf")
    
    plt.rcParams.update(mystyle.rc_parameters)
    


linestyle_u_ref = mystyle.linestyle_u_ref
linestyle_u     = mystyle.linestyle_u

linewidth_u_ref = mystyle.linewidth_u_ref
linewidth_u     = mystyle.linewidth_u

linewidth_rel_error_cn_2 = mystyle.linewidth_rel_error_cn_2
linewidth_rel_error_cn_4 = mystyle.linewidth_rel_error_cn_4
linewidth_rel_error_cn_6 = mystyle.linewidth_rel_error_cn_6
linewidth_rel_error_cn_8 = mystyle.linewidth_rel_error_cn_8

linewidth_rel_error_rk4_2 = mystyle.linewidth_rel_error_rk4_2
linewidth_rel_error_rk4_4 = mystyle.linewidth_rel_error_rk4_4
linewidth_rel_error_rk4_6 = mystyle.linewidth_rel_error_rk4_6
linewidth_rel_error_rk4_8 = mystyle.linewidth_rel_error_rk4_8

linestyle_rel_error_cn_2 = mystyle.linestyle_rel_error_cn_2
linestyle_rel_error_cn_4 = mystyle.linestyle_rel_error_cn_4
linestyle_rel_error_cn_6 = mystyle.linestyle_rel_error_cn_6
linestyle_rel_error_cn_8 = mystyle.linestyle_rel_error_cn_8

linestyle_rel_error_rk4_2 = mystyle.linestyle_rel_error_rk4_2
linestyle_rel_error_rk4_4 = mystyle.linestyle_rel_error_rk4_4
linestyle_rel_error_rk4_6 = mystyle.linestyle_rel_error_rk4_6
linestyle_rel_error_rk4_8 = mystyle.linestyle_rel_error_rk4_8


color_u_ref = mystyle.color_u_ref
color_u     = mystyle.color_u

color_rel_error_cn_2 = mystyle.color_rel_error_cn_2
color_rel_error_cn_4 = mystyle.color_rel_error_cn_4
color_rel_error_cn_6 = mystyle.color_rel_error_cn_6
color_rel_error_cn_8 = mystyle.color_rel_error_cn_8

color_rel_error_rk4_2 = mystyle.color_rel_error_rk4_2
color_rel_error_rk4_4 = mystyle.color_rel_error_rk4_4
color_rel_error_rk4_6 = mystyle.color_rel_error_rk4_6
color_rel_error_rk4_8 = mystyle.color_rel_error_rk4_8

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

t_ticks_major = np.array([0, 0.5, 1.0, 1.5, 2])
t_ticks_minor = np.array([0.25, 0.75, 1.25, 1.75])
          
y_ticks_major_left_column = np.array([0, 0.5, 1.0])
y_ticks_minor_left_column = np.array([0.25, 0.75])

y_ticks_major_right_column = np.array([-1, 0, 1])
y_ticks_minor_right_column = np.array([-0.5, 0.5])

x_min = -8.5
x_max = +8.5

y_min_left_column = -0.1
y_max_left_column = +1.1

y_min_right_column = -1.2
y_max_right_column = +1.2


xlabel = r'$x$'

"""
ylabel_00 = r'$|u_\mathrm{ref}(x,t_0)|^2$'
ylabel_10 = r'$|u_\mathrm{ref}(x,t_1)|^2$'
ylabel_20 = r'$|u_\mathrm{ref}(x,t_2)|^2$'
ylabel_30 = r'$|u_\mathrm{ref}(x,t_3)|^2$'
ylabel_40 = r'$|u_\mathrm{ref}(x,t_4)|^2$'
ylabel_50 = r'$|u_\mathrm{ref}(x,t_5)|^2$'
ylabel_60 = r'$|u_\mathrm{ref}(x,t_6)|^2$'
ylabel_70 = r'$|u_\mathrm{ref}(x,t_7)|^2$'
ylabel_80 = r'$|u_\mathrm{ref}(x,t_8)|^2$'
"""

ylabel_00 = r'$|u(x,t_0)|^2$'
ylabel_10 = r'$|u(x,t_1)|^2$'
ylabel_20 = r'$|u(x,t_2)|^2$'
ylabel_30 = r'$|u(x,t_3)|^2$'
ylabel_40 = r'$|u(x,t_4)|^2$'
ylabel_50 = r'$|u(x,t_5)|^2$'
ylabel_60 = r'$|u(x,t_6)|^2$'
ylabel_70 = r'$|u(x,t_7)|^2$'
ylabel_80 = r'$|u(x,t_8)|^2$'


ylabel_01 = r'$\operatorname{Re}\, u(x,t_0)$'
ylabel_11 = r'$\operatorname{Re}\, u(x,t_1)$'
ylabel_21 = r'$\operatorname{Re}\, u(x,t_2)$'
ylabel_31 = r'$\operatorname{Re}\, u(x,t_3)$'
ylabel_41 = r'$\operatorname{Re}\, u(x,t_4)$'
ylabel_51 = r'$\operatorname{Re}\, u(x,t_5)$'
ylabel_61 = r'$\operatorname{Re}\, u(x,t_6)$'
ylabel_71 = r'$\operatorname{Re}\, u(x,t_7)$'
ylabel_81 = r'$\operatorname{Re}\, u(x,t_8)$'


width  = 8
height = 8

name_fig_1 = "figure_wave_packet_snapshots_cn4"

fig_1 = plt.figure(name_fig_1, figsize=(width, height), facecolor="white", constrained_layout=False)

spacing_x = 0.2
spacing_y = 0.1

gridspec = fig_1.add_gridspec(ncols=2, nrows=9, left=0.065, right=0.99, bottom=0.05, top=0.99, wspace=spacing_x, hspace=spacing_y)



#==========================================================================================
ax_00 = fig_1.add_subplot(gridspec[0, 0])
ax_10 = fig_1.add_subplot(gridspec[1, 0])
ax_20 = fig_1.add_subplot(gridspec[2, 0])
ax_30 = fig_1.add_subplot(gridspec[3, 0])
ax_40 = fig_1.add_subplot(gridspec[4, 0])
ax_50 = fig_1.add_subplot(gridspec[5, 0])
ax_60 = fig_1.add_subplot(gridspec[6, 0])
ax_70 = fig_1.add_subplot(gridspec[7, 0])
ax_80 = fig_1.add_subplot(gridspec[8, 0])

ax_01 = fig_1.add_subplot(gridspec[0, 1])
ax_11 = fig_1.add_subplot(gridspec[1, 1])
ax_21 = fig_1.add_subplot(gridspec[2, 1])
ax_31 = fig_1.add_subplot(gridspec[3, 1])
ax_41 = fig_1.add_subplot(gridspec[4, 1])
ax_51 = fig_1.add_subplot(gridspec[5, 1])
ax_61 = fig_1.add_subplot(gridspec[6, 1])
ax_71 = fig_1.add_subplot(gridspec[7, 1])
ax_81 = fig_1.add_subplot(gridspec[8, 1])
#==========================================================================================

#==========================================================================================
# ax_00.plot(x, np.abs(u_ref_0)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)

print(x.size)
print(u_0.size)

ax_00.plot(x, np.abs(u_0)**2, linewidth=linewidth_u, linestyle=linestyle_u, color=color_u, label=label_u)

ax_00.set_xlim(x_min, x_max)
ax_00.set_ylim(y_min_left_column, y_max_left_column)

ax_00.set_xticks(x_ticks_major, minor=False)
ax_00.set_xticks(x_ticks_minor, minor=True)

ax_00.set_yticks(y_ticks_major_left_column, minor=False)
ax_00.set_yticks(y_ticks_minor_left_column, minor=True)

ax_00.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_00.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_00.set_ylabel(ylabel_00)

ax_00.set_xticklabels([])

# ax_00.legend(loc='upper left', ncol=2)
#==========================================================================================

#==========================================================================================
# ax_10.plot(x, np.abs(u_ref_1)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_10.plot(x, np.abs(u_1)**2, linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_10.set_xlim(x_min, x_max)
ax_10.set_ylim(y_min_left_column, y_max_left_column)


ax_10.set_xticks(x_ticks_major, minor=False)
ax_10.set_xticks(x_ticks_minor, minor=True)

ax_10.set_yticks(y_ticks_major_left_column, minor=False)
ax_10.set_yticks(y_ticks_minor_left_column, minor=True)

ax_10.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_10.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_10.set_ylabel(ylabel_10)
    
ax_10.set_xticklabels([])
#==========================================================================================

#==========================================================================================
# ax_20.plot(x, np.abs(u_ref_2)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_20.plot(x, np.abs(u_2)**2, linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_20.set_xlim(x_min, x_max)
ax_20.set_ylim(y_min_left_column, y_max_left_column)


ax_20.set_xticks(x_ticks_major, minor=False)
ax_20.set_xticks(x_ticks_minor, minor=True)

ax_20.set_yticks(y_ticks_major_left_column, minor=False)
ax_20.set_yticks(y_ticks_minor_left_column, minor=True)

ax_20.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_20.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_20.set_ylabel(ylabel_20)
    
ax_20.set_xticklabels([])
#==========================================================================================

#==========================================================================================
# ax_30.plot(x, np.abs(u_ref_3)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_30.plot(x, np.abs(u_3)**2, linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_30.set_xlim(x_min, x_max)
ax_30.set_ylim(y_min_left_column, y_max_left_column)


ax_30.set_xticks(x_ticks_major, minor=False)
ax_30.set_xticks(x_ticks_minor, minor=True)

ax_30.set_yticks(y_ticks_major_left_column, minor=False)
ax_30.set_yticks(y_ticks_minor_left_column, minor=True)

ax_30.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_30.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_30.set_ylabel(ylabel_30)
    
ax_30.set_xticklabels([])
#==========================================================================================

#==========================================================================================
# ax_40.plot(x, np.abs(u_ref_4)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_40.plot(x, np.abs(u_4)**2, linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_40.set_xlim(x_min, x_max)
ax_40.set_ylim(y_min_left_column, y_max_left_column)


ax_40.set_xticks(x_ticks_major, minor=False)
ax_40.set_xticks(x_ticks_minor, minor=True)

ax_40.set_yticks(y_ticks_major_left_column, minor=False)
ax_40.set_yticks(y_ticks_minor_left_column, minor=True)

ax_40.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_40.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_40.set_ylabel(ylabel_40)
    
ax_40.set_xticklabels([])
#==========================================================================================

#==========================================================================================
# ax_50.plot(x, np.abs(u_ref_5)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_50.plot(x, np.abs(u_5)**2, linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_50.set_xlim(x_min, x_max)
ax_50.set_ylim(y_min_left_column, y_max_left_column)


ax_50.set_xticks(x_ticks_major, minor=False)
ax_50.set_xticks(x_ticks_minor, minor=True)

ax_50.set_yticks(y_ticks_major_left_column, minor=False)
ax_50.set_yticks(y_ticks_minor_left_column, minor=True)

ax_50.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_50.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_50.set_ylabel(ylabel_50)
    
ax_50.set_xticklabels([])
#==========================================================================================

#==========================================================================================
# ax_60.plot(x, np.abs(u_ref_6)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_60.plot(x, np.abs(u_6)**2, linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_60.set_xlim(x_min, x_max)
ax_60.set_ylim(y_min_left_column, y_max_left_column)


ax_60.set_xticks(x_ticks_major, minor=False)
ax_60.set_xticks(x_ticks_minor, minor=True)

ax_60.set_yticks(y_ticks_major_left_column, minor=False)
ax_60.set_yticks(y_ticks_minor_left_column, minor=True)

ax_60.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_60.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_60.set_ylabel(ylabel_60)
    
ax_60.set_xticklabels([])
#==========================================================================================

#==========================================================================================
# ax_70.plot(x, np.abs(u_ref_7)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_70.plot(x, np.abs(u_7)**2, linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_70.set_xlim(x_min, x_max)
ax_70.set_ylim(y_min_left_column, y_max_left_column)


ax_70.set_xticks(x_ticks_major, minor=False)
ax_70.set_xticks(x_ticks_minor, minor=True)

ax_70.set_yticks(y_ticks_major_left_column, minor=False)
ax_70.set_yticks(y_ticks_minor_left_column, minor=True)

ax_70.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_70.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_70.set_ylabel(ylabel_70)
    
ax_70.set_xticklabels([])
#==========================================================================================


#==========================================================================================
# ax_80.plot(x, np.abs(u_ref_8)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_80.plot(x, np.abs(u_8)**2, linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_80.set_xlim(x_min, x_max)
ax_80.set_ylim(y_min_left_column, y_max_left_column)


ax_80.set_xticks(x_ticks_major, minor=False)
ax_80.set_xticks(x_ticks_minor, minor=True)

ax_80.set_yticks(y_ticks_major_left_column, minor=False)
ax_80.set_yticks(y_ticks_minor_left_column, minor=True)

ax_80.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_80.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_80.set_xlabel(xlabel)
ax_80.set_ylabel(ylabel_80)
#==========================================================================================






#==========================================================================================
# ax_01.plot(x, np.abs(u_ref_0)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_01.plot(x, np.real(u_0), linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_01.set_xlim(x_min, x_max)
ax_01.set_ylim(y_min_right_column, y_max_right_column)

ax_01.set_xticks(x_ticks_major, minor=False)
ax_01.set_xticks(x_ticks_minor, minor=True)

ax_01.set_yticks(y_ticks_major_right_column, minor=False)
ax_01.set_yticks(y_ticks_minor_right_column, minor=True)

ax_01.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_01.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_01.set_ylabel(ylabel_01)

ax_01.set_xticklabels([])

# ax_01.legend(loc='upper left', ncol=2)
#==========================================================================================

#==========================================================================================
# ax_11.plot(x, np.abs(u_ref_1)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_11.plot(x, np.real(u_1), linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_11.set_xlim(x_min, x_max)
ax_11.set_ylim(y_min_right_column, y_max_right_column)


ax_11.set_xticks(x_ticks_major, minor=False)
ax_11.set_xticks(x_ticks_minor, minor=True)

ax_11.set_yticks(y_ticks_major_right_column, minor=False)
ax_11.set_yticks(y_ticks_minor_right_column, minor=True)

ax_11.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_11.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_11.set_ylabel(ylabel_11)
    
ax_11.set_xticklabels([])
#==========================================================================================

#==========================================================================================
# ax_21.plot(x, np.abs(u_ref_2)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_21.plot(x, np.real(u_2), linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_21.set_xlim(x_min, x_max)
ax_21.set_ylim(y_min_right_column, y_max_right_column)


ax_21.set_xticks(x_ticks_major, minor=False)
ax_21.set_xticks(x_ticks_minor, minor=True)

ax_21.set_yticks(y_ticks_major_right_column, minor=False)
ax_21.set_yticks(y_ticks_minor_right_column, minor=True)

ax_21.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_21.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_21.set_ylabel(ylabel_21)
    
ax_21.set_xticklabels([])
#==========================================================================================

#==========================================================================================
# ax_31.plot(x, np.abs(u_ref_3)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_31.plot(x, np.real(u_3), linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_31.set_xlim(x_min, x_max)
ax_31.set_ylim(y_min_right_column, y_max_right_column)


ax_31.set_xticks(x_ticks_major, minor=False)
ax_31.set_xticks(x_ticks_minor, minor=True)

ax_31.set_yticks(y_ticks_major_right_column, minor=False)
ax_31.set_yticks(y_ticks_minor_right_column, minor=True)

ax_31.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_31.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_31.set_ylabel(ylabel_31)
    
ax_31.set_xticklabels([])
#==========================================================================================

#==========================================================================================
# ax_41.plot(x, np.abs(u_ref_4)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_41.plot(x, np.real(u_4), linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_41.set_xlim(x_min, x_max)
ax_41.set_ylim(y_min_right_column, y_max_right_column)


ax_41.set_xticks(x_ticks_major, minor=False)
ax_41.set_xticks(x_ticks_minor, minor=True)

ax_41.set_yticks(y_ticks_major_right_column, minor=False)
ax_41.set_yticks(y_ticks_minor_right_column, minor=True)

ax_41.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_41.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_41.set_ylabel(ylabel_41)
    
ax_41.set_xticklabels([])
#==========================================================================================

#==========================================================================================
# ax_51.plot(x, np.abs(u_ref_5)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_51.plot(x, np.real(u_5), linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_51.set_xlim(x_min, x_max)
ax_51.set_ylim(y_min_right_column, y_max_right_column)


ax_51.set_xticks(x_ticks_major, minor=False)
ax_51.set_xticks(x_ticks_minor, minor=True)

ax_51.set_yticks(y_ticks_major_right_column, minor=False)
ax_51.set_yticks(y_ticks_minor_right_column, minor=True)

ax_51.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_51.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_51.set_ylabel(ylabel_51)
    
ax_51.set_xticklabels([])
#==========================================================================================

#==========================================================================================
# ax_61.plot(x, np.abs(u_ref_6)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_61.plot(x, np.real(u_6), linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_61.set_xlim(x_min, x_max)
ax_61.set_ylim(y_min_right_column, y_max_right_column)


ax_61.set_xticks(x_ticks_major, minor=False)
ax_61.set_xticks(x_ticks_minor, minor=True)

ax_61.set_yticks(y_ticks_major_right_column, minor=False)
ax_61.set_yticks(y_ticks_minor_right_column, minor=True)

ax_61.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_61.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_61.set_ylabel(ylabel_61)
    
ax_61.set_xticklabels([])
#==========================================================================================

#==========================================================================================
# ax_71.plot(x, np.abs(u_ref_7)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_71.plot(x, np.real(u_7), linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_71.set_xlim(x_min, x_max)
ax_71.set_ylim(y_min_right_column, y_max_right_column)


ax_71.set_xticks(x_ticks_major, minor=False)
ax_71.set_xticks(x_ticks_minor, minor=True)

ax_71.set_yticks(y_ticks_major_right_column, minor=False)
ax_71.set_yticks(y_ticks_minor_right_column, minor=True)

ax_71.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_71.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_71.set_ylabel(ylabel_71)

ax_71.set_xticklabels([])
#==========================================================================================


#==========================================================================================
# ax_81.plot(x, np.abs(u_ref_8)**2, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=label_u_ref)
ax_81.plot(x, np.real(u_8), linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

ax_81.set_xlim(x_min, x_max)
ax_81.set_ylim(y_min_right_column, y_max_right_column)


ax_81.set_xticks(x_ticks_major, minor=False)
ax_81.set_xticks(x_ticks_minor, minor=True)

ax_81.set_yticks(y_ticks_major_right_column, minor=False)
ax_81.set_yticks(y_ticks_minor_right_column, minor=True)

ax_81.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_81.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_81.set_xlabel(xlabel)
ax_81.set_ylabel(ylabel_81)
#==========================================================================================

















width  = 8
height = 6

name_fig_2 = "figure_wave_packet_time_evolution"

fig_2 = plt.figure(name_fig_2, figsize=(width, height), facecolor="white", constrained_layout=False)

spacing_x = 0.2
spacing_y = 0.1

gridspec = fig_2.add_gridspec(ncols=1, nrows=2, left=0.1, right=0.975, bottom=0.075, top=0.975, wspace=spacing_x, hspace=spacing_y)


ax_00 = fig_2.add_subplot(gridspec[0, 0])
ax_10 = fig_2.add_subplot(gridspec[1, 0])



#==========================================================================================
ax_00.axis([0, T, 1e-5, 1])

ax_00.set_yscale('log')


rel_error_cn_2_of_times_analysis[0] = 1e-6
rel_error_cn_4_of_times_analysis[0] = 1e-6
rel_error_cn_6_of_times_analysis[0] = 1e-6
rel_error_cn_8_of_times_analysis[0] = 1e-6

rel_error_rk4_2_of_times_analysis[0] = 1e-6
rel_error_rk4_4_of_times_analysis[0] = 1e-6
rel_error_rk4_6_of_times_analysis[0] = 1e-6
rel_error_rk4_8_of_times_analysis[0] = 1e-6


ax_00.semilogy(times_analysis, rel_error_cn_2_of_times_analysis,  linewidth=linewidth_rel_error_cn_2,  linestyle=linestyle_rel_error_cn_2,  color=color_rel_error_cn_2,  label=r'$\mathrm{CN}_2$')
ax_00.semilogy(times_analysis, rel_error_cn_4_of_times_analysis,  linewidth=linewidth_rel_error_cn_4,  linestyle=linestyle_rel_error_cn_4,  color=color_rel_error_cn_4,  label=r'$\mathrm{CN}_4$')
ax_00.semilogy(times_analysis, rel_error_cn_6_of_times_analysis,  linewidth=linewidth_rel_error_cn_6,  linestyle=linestyle_rel_error_cn_6,  color=color_rel_error_cn_6,  label=r'$\mathrm{CN}_6$')
ax_00.semilogy(times_analysis, rel_error_cn_8_of_times_analysis,  linewidth=linewidth_rel_error_cn_8,  linestyle=linestyle_rel_error_cn_8,  color=color_rel_error_cn_8,  label=r'$\mathrm{CN}_8$')

ax_00.semilogy(times_analysis, rel_error_rk4_2_of_times_analysis, linewidth=linewidth_rel_error_rk4_2, linestyle=linestyle_rel_error_rk4_2, color=color_rel_error_rk4_2, label=r'$\mathrm{RK4}_2$')
ax_00.semilogy(times_analysis, rel_error_rk4_4_of_times_analysis, linewidth=linewidth_rel_error_rk4_4, linestyle=linestyle_rel_error_rk4_4, color=color_rel_error_rk4_4, label=r'$\mathrm{RK4}_4$')
ax_00.semilogy(times_analysis, rel_error_rk4_6_of_times_analysis, linewidth=linewidth_rel_error_rk4_6, linestyle=linestyle_rel_error_rk4_6, color=color_rel_error_rk4_6, label=r'$\mathrm{RK4}_6$')
ax_00.semilogy(times_analysis, rel_error_rk4_8_of_times_analysis, linewidth=linewidth_rel_error_rk4_8, linestyle=linestyle_rel_error_rk4_8, color=color_rel_error_rk4_8, label=r'$\mathrm{RK4}_8$')

ax_00.set_xticks(t_ticks_major, minor=False)
ax_00.set_xticks(t_ticks_minor, minor=True)


ax_00.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_00.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_00.set_ylabel(r'$\|\bm{u}(t) - \bm{u}_\mathrm{ref}(t) \|_2 / \| \bm{u}_\mathrm{ref}(t) \|_2$')

ax_00.set_xticklabels([])

ax_00.legend(loc='upper right', ncol=2)
#==========================================================================================

#==========================================================================================
ax_10.axis([0, T, 1e-17, 1e-5])

ax_10.set_yscale('log')


deviation_mass_cn_2_of_times_analysis[0] = 1e-18
deviation_mass_cn_4_of_times_analysis[0] = 1e-18
deviation_mass_cn_6_of_times_analysis[0] = 1e-18
deviation_mass_cn_8_of_times_analysis[0] = 1e-18

deviation_mass_rk4_2_of_times_analysis[0] = 1e-18
deviation_mass_rk4_4_of_times_analysis[0] = 1e-18
deviation_mass_rk4_6_of_times_analysis[0] = 1e-18
deviation_mass_rk4_8_of_times_analysis[0] = 1e-18


ax_10.plot(times_analysis, deviation_mass_cn_2_of_times_analysis,  linewidth=linewidth_rel_error_cn_2,  linestyle=linestyle_rel_error_cn_2,  color=color_rel_error_cn_2,  label=label_u)
ax_10.plot(times_analysis, deviation_mass_cn_4_of_times_analysis,  linewidth=linewidth_rel_error_cn_4,  linestyle=linestyle_rel_error_cn_4,  color=color_rel_error_cn_4,  label=label_u)
ax_10.plot(times_analysis, deviation_mass_cn_6_of_times_analysis,  linewidth=linewidth_rel_error_cn_6,  linestyle=linestyle_rel_error_cn_6,  color=color_rel_error_cn_6,  label=label_u)
ax_10.plot(times_analysis, deviation_mass_cn_8_of_times_analysis,  linewidth=linewidth_rel_error_cn_8,  linestyle=linestyle_rel_error_cn_8,  color=color_rel_error_cn_8,  label=label_u)

ax_10.plot(times_analysis, deviation_mass_rk4_2_of_times_analysis, linewidth=linewidth_rel_error_rk4_2, linestyle=linestyle_rel_error_rk4_2, color=color_rel_error_rk4_2, label=label_u)
ax_10.plot(times_analysis, deviation_mass_rk4_4_of_times_analysis, linewidth=linewidth_rel_error_rk4_4, linestyle=linestyle_rel_error_rk4_4, color=color_rel_error_rk4_4, label=label_u)
ax_10.plot(times_analysis, deviation_mass_rk4_6_of_times_analysis, linewidth=linewidth_rel_error_rk4_6, linestyle=linestyle_rel_error_rk4_6, color=color_rel_error_rk4_6, label=label_u)
ax_10.plot(times_analysis, deviation_mass_rk4_8_of_times_analysis, linewidth=linewidth_rel_error_rk4_8, linestyle=linestyle_rel_error_rk4_8, color=color_rel_error_rk4_8, label=label_u)


ax_10.set_xticks(t_ticks_major, minor=False)
ax_10.set_xticks(t_ticks_minor, minor=True)




majorLocator = FixedLocator([1e-17, 1e-15, 1e-13, 1e-11, 1e-9, 1e-7, 1e-5])
minorLocator = mpl.ticker.LogLocator(base=10.0, subs=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks=100)


ax_10.yaxis.set_major_locator(majorLocator)
ax_10.yaxis.set_minor_locator(minorLocator)

ax_10.yaxis.set_minor_formatter(NullFormatter())




ax_10.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_10.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_10.set_xlabel(r'$t$')
ax_10.set_ylabel(r'$\big| 1 - \| \bm{u}(t) \|_2 / \| \bm{u}(0) \|_2 \big|$')
#==========================================================================================




plt.draw()



if export_pdf == True:
    
    path = "/home/jfmennemann/git/nls/pdf/"
    
    
    plt.figure(name_fig_1)
    
    filepath = path + name_fig_1 + ".pdf"

    plt.savefig(filepath, backend='pgf')
    
    
    plt.figure(name_fig_2)
    
    filepath = path + name_fig_2 + ".pdf"

    plt.savefig(filepath, backend='pgf')

else:

    plt.show()









# u_complete = np.zeros_like(x_complete, dtype=np.complex128)

# u_complete[0:-1] = u



