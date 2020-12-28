from scipy.sparse import eye

from scipy.sparse.linalg import spsolve

import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=8)

from linear_schroedinger.wave_packet.reference_solution import gaussian

from linear_schroedinger.wave_packet.figure_1 import Figure1

from differentiation import finite_differences_1d




x0 = -2.5

sigma_0 = 1.0

k0 = 4.0



x_min = -10
x_max = +10

Jx = 250

x = np.linspace(x_min, x_max, Jx+1, endpoint=True)

dx = x[1] - x[0]



T = 5

dt = 0.001

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)



D2 = finite_differences_1d.get_D2_circulant_2nd_order(Jx, dx)


print(dx**2 * D2.toarray())

input('press any key ...')




E = eye(Jx)


A = E - 0.5 * 1j * dt * D2












u_ref = gaussian(x, 0.0, x0, sigma_0, k0)

u = u_ref[0:-1]

assert(u.size == Jx)


u_complete = np.zeros_like(u_ref)

u_complete[0:-1] = u

u_complete[-1] = u_complete[0]




n_mod_times_analysis = 50

times_analysis = times[::n_mod_times_analysis]


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
        
        u_ref = gaussian(x, t, x0, sigma_0, k0)
    
        u_complete[0:-1] = u
        u_complete[-1] = u_complete[0]
        
        rel_error_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_complete-u_ref) / np.linalg.norm(u_complete)
        times_analysis[nr_times_analysis] = t 
        
        
        fig_1.update_u(u_complete, u_ref)
        # fig_1.update_rel_error(rel_error_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.redraw()
        
        
        nr_times_analysis = nr_times_analysis + 1
    
    u = spsolve(A, u)
    
    

input('press any key ...')






