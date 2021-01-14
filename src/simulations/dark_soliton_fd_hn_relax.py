"""
Equation: cubic nonlinear Schroedinger equation
Initial condition: dark soliton
Spatial approximation: finite differences
Boundary conditions: homogeneous Neumann
Time-integration method: relaxation method
"""


from scipy.sparse import diags
from scipy.sparse import eye
from scipy.sparse import spdiags
from scipy.sparse import csr_matrix

from scipy.sparse.linalg import spsolve


import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)


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
#------------------------------------------------------------------------------


D_xx = diags([1, -2, 1], [-1, 0, 1], shape=(Jx+1, Jx+1))
    
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

psi_old = np.abs(u_ref)**2




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
    
    
    psi_new = 2 * np.abs(u)**2 - psi_old
    
    b = u + 0.25 * 1j * dt * D_xx * u - 0.5 * 1j * dt * beta * psi_new * u
    
    # create matrix A using psi_new as approximation for the density
    A = eye(Jx+1) - 0.25 * 1j * dt * D_xx + 0.5 * 1j * dt * beta * spdiags(psi_new, 0, Jx+1, Jx+1)
    
    
    # use one-sided finite difference approximation for the 1st derivative at the boundaries
    A = A.todense()
    
    A[0, 0] = -3
    A[0, 1] = +4
    A[0, 2] = -1
    
    A[-1, -1] = +3
    A[-1, -2] = -4
    A[-1, -3] = +1

    # 1st derivative at the boundaries is enforced to be zero
    b[0] = 0
    b[-1] = 0

    A = csr_matrix(A)
    
    
    
    u = spsolve(A, b)
    
    psi_old = psi_new
    
    
    

input('press any key ...')






