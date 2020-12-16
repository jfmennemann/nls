from scipy.sparse import diags
from scipy.sparse import eye
from scipy.sparse import spdiags
from scipy.sparse import csr_matrix

from scipy.sparse.linalg import spsolve

import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)



order_spatial_discretization = 2



from nonlinear_schroedinger.solitons.dark_soliton.reference_solutions import dark_soliton

from nonlinear_schroedinger.solitons.dark_soliton.figure_1 import Figure1




v = 2.5

x0 = -4
theta_0 = 0

u0 = 10
beta = 2
phi = 0.0 * np.pi




x_min = -8
x_max = +8

L = x_max - x_min

Jx = 400

x = np.linspace(x_min, x_max, Jx+1, endpoint=True)

dx = x[1] - x[0]



T = 4

dt = 0.001

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)





if order_spatial_discretization == 2:
    
    D_xx = diags([1, -2, 1], [-1, 0, 1], shape=(Jx+1, Jx+1))
    
    D_xx = csr_matrix(D_xx)
    
    D_xx = D_xx / dx**2

"""   
if order_spatial_discretization == 4:
    
    D_xx = diags([16, -1, -1, 16, -30, 16, -1, -1, 16], [-(Jx-1), -(Jx-2), -2, -1, 0, 1, 2, Jx-2, Jx-1], shape=(Jx, Jx))
    
    D_xx = D_xx / (12 * dx**2)
    
if order_spatial_discretization == 6:
    
    D_xx = diags([270, -27, 2, 2, -27, 270, -490, 270, -27, 2, 2, -27, 270], [-(Jx-1), -(Jx-2), -(Jx-3), -3, -2, -1, 0, 1, 2, 3, Jx-3, Jx-2, Jx-1], shape=(Jx, Jx))
    
    D_xx = D_xx / (180 * dx**2)
    
if order_spatial_discretization == 8:
    
    D_xx = diags([8064, -1008, 128, -9, -9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9, -9, 128, -1008, 8064], [-(Jx-1), -(Jx-2), -(Jx-3), -(Jx-4), -4, -3, -2, -1, 0, 1, 2, 3, 4, Jx-4, Jx-3, Jx-2, Jx-1], shape=(Jx, Jx))
        
    D_xx = D_xx / (5040 * dx**2)
"""








u_ref = dark_soliton(x, 0, v, x0, theta_0, u0, beta, phi)
u = u_ref

psi_old = np.abs(u_ref)**2

assert(u.size == Jx+1)
assert(psi_old.size ==Jx+1)







n_mod_times_analysis = 20

times_analysis = times[::n_mod_times_analysis]



norm_u_of_times_analysis = np.zeros_like(times_analysis)
rel_error_of_times_analysis = np.zeros_like(times_analysis)




fig_1 = Figure1(x, times, screen_size='large')

fig_1.update_u(u, u_ref)
   
fig_1.redraw()




nr_times_analysis = 0

for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        
        u_ref = dark_soliton(x, t, v, x0, theta_0, u0, beta, phi)
        
        
        norm_u_of_times_analysis[nr_times_analysis] = np.linalg.norm(u)
        
        
        print(norm_u_of_times_analysis[nr_times_analysis] / norm_u_of_times_analysis[0])
        
        defect_of_mass_of_times_analysis = np.abs(1.0 - norm_u_of_times_analysis / norm_u_of_times_analysis[0])
        
        
        
        rel_error_of_times_analysis[nr_times_analysis] = np.linalg.norm(u-u_ref) / np.linalg.norm(u_ref)
        
        times_analysis[nr_times_analysis] = t 
        
        
        fig_1.update_u(u, u_ref)
        
        fig_1.update_rel_error(rel_error_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.update_defect_of_mass(defect_of_mass_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.redraw()
        
        
        nr_times_analysis = nr_times_analysis + 1
    
    
    psi_new = 2 * np.abs(u)**2 - psi_old
    
    b = u + 0.25 * 1j * dt * D_xx * u - 0.5 * 1j * dt * beta * psi_new * u
    
    A = eye(Jx+1) - 0.25 * 1j * dt * D_xx + 0.5 * 1j * dt * beta * spdiags(psi_new, 0, Jx+1, Jx+1)
    
    
    
    A = A.todense()
    
    A[0, 0] = -3
    A[0, 1] = +4
    A[0, 2] = -1
    
    A[-1, -1] = +3
    A[-1, -2] = -4
    A[-1, -3] = +1

    b[0] = 0
    b[-1] = 0

    A = csr_matrix(A)
    
    
    
    u = spsolve(A, b)
    
    psi_old = psi_new
    
    
    

input('press any key ...')






