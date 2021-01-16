"""
Equation: Gross-Pitaevskii equation
Initial condition: fat soliton
Spatial approximation: finite differences
Boundary conditions: periodic
Time-integration method: RK4
"""


from scipy.sparse import spdiags


import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)


from simulations.reference_solutions import Phi

from simulations.figure_2 import Figure2


from differentiation import finite_differences_1d


order_spatial_discretization = 2





x_min = -20
x_max = +20

L = x_max - x_min

Jx = 200

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






alpha = 4
nu = 0.9

Phi_of_x = Phi(x, 0.0, alpha, nu)


u_ref = Phi_of_x
u     = Phi_of_x


V = alpha * Phi_of_x

diag_V = spdiags(np.array([V[0:-1]]), np.array([0]), Jx, Jx, 'csr')



#------------------------------------------------------------------------------
T = 20

dt = 0.001

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)

dt_new = times[1] - times[0]

assert(dt_new == dt)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
n_mod_times_analysis = 25
# n_mod_times_analysis = 1
#------------------------------------------------------------------------------










fig_1 = Figure2(x, 0.0, 1.5, -4, 1, V, u_ref)

fig_1.update_u(u, u_ref)

fig_1.update_V(V)

fig_1.redraw()



def eval_f(y):
    
    return 1j * D2 * y - 1j * diag_V * y - 1j * np.abs(y)**2 * y

for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        # u_ref = Phi(x, 0.0, alpha, nu)
        
        fig_1.update_u(u, u_ref)
        
        fig_1.redraw()
    
    
    k1 = eval_f(u)
    k2 = eval_f(u + 0.5 * dt * k1)
    k3 = eval_f(u + 0.5 * dt * k2)
    k4 = eval_f(u + 1.0 * dt * k3)
    
    u = u + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
    

input('press any key ...')






