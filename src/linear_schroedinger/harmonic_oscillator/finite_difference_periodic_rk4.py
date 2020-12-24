import numpy as np

from numpy import array

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)


from scipy.sparse import diags
from scipy.sparse import spdiags


from linear_schroedinger.harmonic_oscillator.reference_solutions import coherent_state

from linear_schroedinger.harmonic_oscillator.figure_1 import Figure1





order_spatial_discretization = 8



x0 = 1

omega = 5






x_min = -4
x_max = +4

Jx = 200

x = np.linspace(x_min, x_max, Jx+1, endpoint=True)

dx = x[1] - x[0]









if order_spatial_discretization == 2:
     
    D_xx = diags([1, 1, -2, 1, 1], [-(Jx-1), -1, 0, 1, (Jx-1)], shape=(Jx, Jx))
    
    D_xx = D_xx / dx**2
    
if order_spatial_discretization == 4:
    
    D_xx = diags([16, -1, -1, 16, -30, 16, -1, -1, 16], [-(Jx-1), -(Jx-2), -2, -1, 0, 1, 2, Jx-2, Jx-1], shape=(Jx, Jx))

    D_xx = D_xx / (12 * dx**2)
    
if order_spatial_discretization == 6:
    
    D_xx = diags([270, -27, 2, 2, -27, 270, -490, 270, -27, 2, 2, -27, 270], [-(Jx-1), -(Jx-2), -(Jx-3), -3, -2, -1, 0, 1, 2, 3, Jx-3, Jx-2, Jx-1], shape=(Jx, Jx))
    
    D_xx = D_xx / (180 * dx**2)
    
if order_spatial_discretization == 8:
    
    D_xx = diags([8064, -1008, 128, -9, -9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9, -9, 128, -1008, 8064], [-(Jx-1), -(Jx-2), -(Jx-3), -(Jx-4), -4, -3, -2, -1, 0, 1, 2, 3, 4, Jx-4, Jx-3, Jx-2, Jx-1], shape=(Jx, Jx))
    
    D_xx = D_xx / (5040 * dx**2)






V = 0.5 * omega**2 * x**2

diag_V = spdiags(array([V[0:-1]]), array([0]), Jx, Jx, 'csr')


A = 1j * (0.5 * D_xx - diag_V)



eigenvalues_A, eigenvectors_A = np.linalg.eig(A.todense())

max_abs_lambda = np.max(np.abs(eigenvalues_A))




dt = 1.0 * np.sqrt(8) / max_abs_lambda

print(dt)

input('press any key to continue ... ')





T = 100

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)





u_ref = coherent_state(x, 0.0, x0, omega)

u = u_ref[0:-1]

assert(u.size == Jx)



u_complete = np.zeros_like(u_ref)

u_complete[0:-1] = u

u_complete[-1] = u_complete[0]






n_mod_times_analysis = 50

times_analysis = times[::n_mod_times_analysis]



norm_u_of_times_analysis = np.zeros_like(times_analysis)
rel_error_of_times_analysis = np.zeros_like(times_analysis)




fig_1 = Figure1(x, times, screen_size='large')

fig_1.update_u(u_complete, u_ref)

fig_1.update_V(V)
    

fig_1.redraw()



nr_times_analysis = 0

for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        """
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        """
        
        
        u_ref = coherent_state(x, t, x0, omega)
        
        u_complete[0:-1] = u
        u_complete[-1] = u_complete[0]
        
        
        norm_u_of_times_analysis[nr_times_analysis] = np.linalg.norm(u)
        
        defect_of_mass_of_times_analysis = np.abs(1.0 - norm_u_of_times_analysis / norm_u_of_times_analysis[0])
        
        print(defect_of_mass_of_times_analysis[nr_times_analysis])
        
        
        rel_error_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_complete-u_ref) / np.linalg.norm(u_ref)
        
        
        times_analysis[nr_times_analysis] = t 
        
        
        fig_1.update_u(u_complete, u_ref)
        
        fig_1.update_rel_error(rel_error_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.update_defect_of_mass(defect_of_mass_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.redraw()
        
        
        nr_times_analysis = nr_times_analysis + 1
    
    
    k1 = A * u
    k2 = A * (u + 0.5 * dt * k1)
    k3 = A * (u + 0.5 * dt * k2)
    k4 = A * (u + 1.0 * dt * k3)
    
    u = u + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
    

input('press any key ...')






