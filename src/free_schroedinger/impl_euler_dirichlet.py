from scipy.sparse import diags, eye

from scipy.sparse.linalg import spsolve

import numpy as np

from free_schroedinger.figure_1 import Figure1


x_min = -10
x_max = +10

Jx = 250

x = np.linspace(x_min, x_max, Jx+1, endpoint=True)

dx = x[1] - x[0]



T = 5

dt = 0.001

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)


D_xx = diags([1, -2, 1], [-1, 0, 1], shape=(Jx-1, Jx-1))

# print(D_xx.toarray())


D_xx = D_xx / dx**2

E = eye(Jx-1)


A = E - 0.5 * 1j * dt * D_xx



def gaussian(x, t, x0, sigma_0, k0):
    
    tau = 2 * sigma_0**2
    
    alpha = 1.0 + 1j * t / tau
    
    return (1.0/np.sqrt(alpha)) * np.exp( (1.0/alpha) * ( -((x-x0)/(2*sigma_0))**2 + 1j * k0 * (x-x0) - 1j * sigma_0**2 * k0**2 * t/ tau) )





x0 = -4.0

sigma_0 = 1.0

k0 = 4


u_ref = gaussian(x, 0.0, x0, sigma_0, k0)

u = u_ref[1:-1]

assert(u.size == Jx-1)


u_complete = np.zeros_like(u_ref)

u_complete[1:-1] = u[:]





n_mod_times_analysis = 25

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
    
        u_complete[1:-1] = u[:]
        
        rel_error_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_complete-u_ref) / np.linalg.norm(u_complete)
        times_analysis[nr_times_analysis] = t 
        
        
        fig_1.update_u(u_complete, u_ref)
        fig_1.update_rel_error(rel_error_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.redraw()
        
        
        nr_times_analysis = nr_times_analysis + 1
    
    u = spsolve(A, u)
    

input('press any key ...')






