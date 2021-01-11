from scipy.sparse import eye
from scipy.sparse import spdiags

from scipy.sparse.linalg import spsolve

import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)

from nonlinear_schroedinger.solitons.bright_soliton.reference_solution import bright_soliton

from nonlinear_schroedinger.solitons.bright_soliton.figure_1 import Figure1


from differentiation import finite_differences_1d




order_spatial_discretization = 8


a = 4
v = 2
x0 = 0
theta_0 = 0
beta = -1




x_min = -4
x_max = +4

L = x_max - x_min

Jx = 200

x = np.linspace(x_min, x_max, Jx+1, endpoint=True)

dx = x[1] - x[0]





if order_spatial_discretization == 2:
    
    D2 = finite_differences_1d.get_D2_circulant_2nd_order(Jx, dx)


if order_spatial_discretization == 4:
    
    D2 = finite_differences_1d.get_D2_circulant_4th_order(Jx, dx)
    
    
if order_spatial_discretization == 6:
    
    D2 = finite_differences_1d.get_D2_circulant_6th_order(Jx, dx)
    
    
if order_spatial_discretization == 8:
    
    D2 = finite_differences_1d.get_D2_circulant_8th_order(Jx, dx)






dt = 0.001






T = 10

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)










u = bright_soliton(x, 0.0, a, v, -2, theta_0, beta) + bright_soliton(x, 0.0, a, -v, +2, theta_0, beta)

u = u[0:-1]

psi_old = np.abs(u)**2

assert(u.size == Jx)
assert(psi_old.size ==Jx)


u_complete = np.zeros_like(x, dtype=np.complex128)

u_complete[0:-1] = u

u_complete[-1] = u_complete[0]


u_ref = np.nan * np.ones_like(u_complete)



n_mod_times_analysis = 10

times_analysis = times[::n_mod_times_analysis]



norm_u_of_times_analysis = np.zeros_like(times_analysis)
rel_error_of_times_analysis = np.zeros_like(times_analysis)




fig_1 = Figure1(x, times, screen_size='large')

fig_1.update_u(u_complete, u_ref)
   
fig_1.redraw()




nr_times_analysis = 0

for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        u_complete[0:-1] = u
        u_complete[-1] = u_complete[0]
        
        
        norm_u_of_times_analysis[nr_times_analysis] = np.linalg.norm(u)
        
        
        print(norm_u_of_times_analysis[nr_times_analysis] / norm_u_of_times_analysis[0])
        
        defect_of_mass_of_times_analysis = np.abs(1.0 - norm_u_of_times_analysis / norm_u_of_times_analysis[0])
        
        
        
        
        
        times_analysis[nr_times_analysis] = t 
        
        
        fig_1.update_u(u_complete, u_ref)
        
        fig_1.update_rel_error(rel_error_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.update_defect_of_mass(defect_of_mass_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.redraw()
        
        
        nr_times_analysis = nr_times_analysis + 1
    
    
    psi_new = 2 * np.abs(u)**2 - psi_old
    
    b = u + 0.25 * 1j * dt * D2 * u - 0.5 * 1j * dt * beta * psi_new * u
    
    A = eye(Jx) - 0.25 * 1j * dt * D2 + 0.5 * 1j * dt * beta * spdiags(psi_new, 0, Jx, Jx)
    
    u = spsolve(A, b)
    
    psi_old = psi_new
    
    
    

input('press any key ...')






