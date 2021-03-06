from scipy.sparse import diags, eye
from scipy.sparse import spdiags

from scipy.sparse.linalg import spsolve

import numpy as np

from numpy import array

from numpy import pi
from numpy import exp

from linear_schroedinger.figure_1 import Figure1



def coherent_state(x, t, x0, omega):
    
    return (omega/pi)**0.25 * exp( -(omega/2) * ( x**2 - 2*x*x0*exp(-1j*omega*t) + (x0**2/2.0) * exp(-2*1j*omega*t) + x0**2/2.0 ) - 1j * omega * t / 2.0 )



order_spatial_discretization = 8





x0 = 1

omega = 5









x_min = -5
x_max = +5

Jx = 200

x = np.linspace(x_min, x_max, Jx+1, endpoint=True)

dx = x[1] - x[0]



T = 10

dt = 0.001

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)











if order_spatial_discretization == 2:
    
    D_xx = diags([1, -2, 1], [-1, 0, 1], shape=(Jx-1, Jx-1))
    
    D_xx = D_xx / dx**2
    
if order_spatial_discretization == 4:
    
    D_xx = diags([-1, 16, -30, 16, -1], [-2, -1, 0, 1, 2], shape=(Jx-1, Jx-1))

    D_xx = D_xx / (12 * dx**2)
    
if order_spatial_discretization == 6:
    
    D_xx = diags([2, -27, 270, -490, 270, -27, 2], [-3, -2, -1, 0, 1, 2, 3], shape=(Jx-1, Jx-1))
    
    D_xx = D_xx / (180 * dx**2)
    
if order_spatial_discretization == 8:
    
    D_xx = diags([-9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9], [-4, -3, -2, -1, 0, 1, 2, 3, 4], shape=(Jx-1, Jx-1))
    
    D_xx = D_xx / (5040 * dx**2)






V = 0.5 * omega**2 * x**2

diag_V = spdiags(array([V[1:-1]]), array([0]), Jx-1, Jx-1, 'csr')


E = eye(Jx-1)


A = E  -  0.25 * dt * 1j * D_xx  +  0.5 * dt * 1j * diag_V
B = E  +  0.25 * dt * 1j * D_xx  -  0.5 * dt * 1j * diag_V







u_ref = coherent_state(x, 0.0, x0, omega)


u = u_ref[1:-1]

assert(u.size == Jx-1)


u_complete = np.zeros_like(u_ref)

u_complete[1:-1] = u[:]





n_mod_times_analysis = 25

times_analysis = times[::n_mod_times_analysis]


rel_error_of_times_analysis = np.zeros_like(times_analysis)




fig_1 = Figure1(x, times, screen_size='large')

fig_1.update_u(u_complete, u_ref)

fig_1.update_V(V)
    

fig_1.redraw()



nr_times_analysis = 0

for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        # u_ref = gaussian(x, t, x0, sigma_0, k0)
        u_ref = coherent_state(x, t, x0, omega)

        
        u_complete[1:-1] = u[:]
        
        rel_error_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_complete-u_ref) / np.linalg.norm(u_complete)
        times_analysis[nr_times_analysis] = t 
        
        
        fig_1.update_u(u_complete, u_ref)
        fig_1.update_rel_error(rel_error_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.redraw()
        
        
        nr_times_analysis = nr_times_analysis + 1
    
    u = spsolve(A, B*u)
    

input('press any key ...')






