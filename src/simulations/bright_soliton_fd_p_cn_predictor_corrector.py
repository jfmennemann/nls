"""
Equation: cubic nonlinear Schroedinger equation
Initial condition: bright soliton
Spatial approximation: finite differences
Boundary conditions: periodic
Time-integration method: predictor-corrector (Crank-Nicolson)
"""


from scipy.sparse import eye
from scipy.sparse import spdiags

from scipy.sparse.linalg import spsolve


import numpy as np

from numpy import abs

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)


from simulations.reference_solutions import bright_soliton_periodic

from simulations.figure_1 import Figure1


from differentiation import finite_differences_1d



order_spatial_discretization = 8



x_min = -4
x_max = +4

L = x_max - x_min

Jx = 100

x_complete = np.linspace(x_min, x_max, Jx+1, endpoint=True)

x = x_complete[0:-1]

dx = x[1] - x[0]





if order_spatial_discretization == 2:
    
    D2 = finite_differences_1d.get_D2_circulant_2nd_order(Jx, dx)


if order_spatial_discretization == 4:
    
    D2 = finite_differences_1d.get_D2_circulant_4th_order(Jx, dx)
    
    
if order_spatial_discretization == 6:
    
    D2 = finite_differences_1d.get_D2_circulant_6th_order(Jx, dx)
    
    
if order_spatial_discretization == 8:
    
    D2 = finite_differences_1d.get_D2_circulant_8th_order(Jx, dx)




#------------------------------------------------------------------------------
T = 4

dt = 0.001

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)

dt_new = times[1] - times[0]

assert(dt_new == dt)
#------------------------------------------------------------------------------

n_mod_times_analysis = 25



a = 4
v = 2
x0 = 0
theta_0 = 0
beta = -1

u_ref = bright_soliton_periodic(x, 0.0, x0, theta_0, a, v, beta, L)

u = u_ref



fig_1 = Figure1(x, 20, 0, None, u_ref)

fig_1.update_u(u, u_ref)

fig_1.redraw()



nr_times_analysis = 0

for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        u_ref = bright_soliton_periodic(x, t, x0, theta_0, a, v, beta, L)
        
        fig_1.update_u(u, u_ref)
        
        fig_1.redraw()
        
    
    u_abs_squared = abs(u)**2
    
    b =       u + 0.25 * 1j * dt * D2 * u - 0.5 * 1j * dt * beta * u_abs_squared * u
    
    A = eye(Jx) - 0.25 * 1j * dt * D2     + 0.5 * 1j * dt * beta * spdiags(u_abs_squared, np.array([0]), Jx, Jx)
    
    u_star = spsolve(A, b)
    
    
    
    u_star_abs_squared = abs(u_star)**2
    
    A_star = eye(Jx) - 0.25 * 1j * dt * D2     + 0.5 * 1j * dt * beta * spdiags(u_star_abs_squared, np.array([0]), Jx, Jx)
    
    u = spsolve(A_star, b)
    

input('press any key ...')






