"""
Equation: linear Schroedinger equation
Initial condition: coherent state
Spatial approximation: finite differences
Boundary conditions: periodic
Time-integration method: Crank-Nicolson
"""



from scipy.sparse import eye, spdiags
from scipy.sparse.linalg import spsolve


import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)


from simulations.reference_solutions import coherent_state

from simulations.figure_1 import Figure1


from differentiation import finite_differences_1d




order_spatial_discretization = 8



x0 = 1

omega = 5

#------------------------------------------------------------------------------
x_min = -3
x_max = +3

L = x_max - x_min

Jx = 200

x_complete = np.linspace(x_min, x_max, Jx+1, endpoint=True)

x = x_complete[0:-1]

dx = x[1] - x[0]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
T = 2

dt = 0.0025

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)

dt_new = times[1] - times[0]

assert(dt_new == dt)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
n_mod_times_analysis = 25
# n_mod_times_analysis = 1
#------------------------------------------------------------------------------





if order_spatial_discretization == 2:
    
    D2 = finite_differences_1d.get_D2_circulant_2nd_order(Jx, dx)


if order_spatial_discretization == 4:
    
    D2 = finite_differences_1d.get_D2_circulant_4th_order(Jx, dx)
    
    
if order_spatial_discretization == 6:
    
    D2 = finite_differences_1d.get_D2_circulant_6th_order(Jx, dx)
    
    
if order_spatial_discretization == 8:
    
    D2 = finite_differences_1d.get_D2_circulant_8th_order(Jx, dx)






V = 0.5 * omega**2 * x**2

diag_V = spdiags(np.array([V]), np.array([0]), Jx, Jx, 'csr')


E = eye(Jx)

A = E  -  0.25 * dt * 1j * D2  +  0.5 * dt * 1j * diag_V
B = E  +  0.25 * dt * 1j * D2  -  0.5 * dt * 1j * diag_V



u_ref = coherent_state(x, 0.0, x0, omega)

u = u_ref




fig_1 = Figure1(x, 1.5, 200, V, u_ref)

fig_1.update_u(u, u_ref)

fig_1.update_V(V)

fig_1.redraw()


for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        u_ref = coherent_state(x, t, x0, omega)
        
        fig_1.update_u(u, u_ref)
        
        fig_1.redraw()
    
    u = spsolve(A, B*u)
    
input('press any key ...')






