"""
Equation: cubic nonlinear Schroedinger equation
Initial condition: dark soliton
Spatial approximation: finite differences
Boundary conditions: free boundary conditions
Time-integration method: RK4
"""

import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)


from scipy.sparse import diags

from scipy.sparse import csr_matrix



from simulations.reference_solutions import dark_soliton

from simulations.figure_1 import Figure1



x_min = -8
x_max = +8

L = x_max - x_min

Jx = 200

x = np.linspace(x_min, x_max, Jx+1, endpoint=True)

dx = x[1] - x[0]



#------------------------------------------------------------------------------
T = 8

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



# second-order accurate finite difference approximation of the second derivative
D_xx = diags([1, -2, 1], [-1, 0, 1], shape=(Jx+1, Jx+1))

D_xx = D_xx.todense()

# use one-sided finite difference coefficients at the boundaries
D_xx[0, 0] = -1
D_xx[0, 1] = +4
D_xx[0, 2] = -5
D_xx[0, 3] = +2

D_xx[-1, -1] = +2
D_xx[-1, -2] = -5
D_xx[-1, -3] = +4
D_xx[-1, -4] = -1

D_xx = csr_matrix(D_xx)

D_xx = D_xx / dx**2



v = 1

x0 = -5
theta_0 = 0

u0 = 10
beta = 1
phi = 0.0 * np.pi



u_ref = dark_soliton(x, 0, x0, theta_0, u0, v, beta, phi)

u = u_ref
 


def eval_f(y):
    
    return 0.5 * 1j * D_xx * y - 1j * beta * np.abs(y)**2 * y




fig_1 = Figure1(x, 20, 0, None, u_ref)

fig_1.update_u(u, u_ref)
   
fig_1.redraw()


for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        u_ref = dark_soliton(x, t, x0, theta_0, u0, v, beta, phi)
        
        fig_1.update_u(u, u_ref)
        
        fig_1.redraw()
    
    
    k1 = eval_f(u)
    k2 = eval_f(u + 0.5 * dt * k1)
    k3 = eval_f(u + 0.5 * dt * k2)
    k4 = eval_f(u + 1.0 * dt * k3)
    
    u = u + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
    
input('press any key ...')






