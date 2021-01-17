"""
Equation: free Schroedinger equation
Initial condition: Gaussian wave packet
Spatial approximation: finite differences
Boundary conditions: homogeneous Dirichlet
Time-integration method: Crank-Nicolson
"""

from scipy.sparse import eye

from scipy.sparse.linalg import spsolve


import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)


from simulations.reference_solutions import gaussian_periodic

from simulations.figure_1 import Figure1


from differentiation import finite_differences_1d



order_spatial_discretization = 2




x_min = -8
x_max = +8

L = x_max - x_min

Jx = 200

x_complete = np.linspace(x_min, x_max, Jx+1, endpoint=True)

x = x_complete[1:-1]

dx = x[1] - x[0]



x0 = 0

sigma_0 = 0.5

k0 = 4


u_ref = gaussian_periodic(x, 0.0, x0, k0, sigma_0, L)

u = u_ref


#------------------------------------------------------------------------------
T = 4

dt = 0.0025

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)

dt_new = times[1] - times[0]

assert(dt_new == dt)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
n_mod_times_analysis = 25
#------------------------------------------------------------------------------



D2 = finite_differences_1d.get_D2_dirichlet_2nd_order(Jx-1, dx)

E = eye(Jx-1)

A = E - 0.25 * 1j * dt * D2
B = E + 0.25 * 1j * dt * D2




u = u_ref



fig_1 = Figure1(x, 1, 0, None, u_ref)

fig_1.update_u(u, u_ref)

fig_1.redraw()




nr_times_analysis = 0

for n in np.arange(times.size+1):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        u_ref = gaussian_periodic(x, t, x0, k0, sigma_0, L)
                
        fig_1.update_u(u, u_ref)
        
        fig_1.redraw()
        
    
    u = spsolve(A, B*u)
    




