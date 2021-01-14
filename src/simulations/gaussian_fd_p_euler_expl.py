from scipy.sparse import eye


import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)


from simulations.reference_solutions import gaussian

from simulations.figure_1 import Figure1


from differentiation import finite_differences_1d



x_min = -8
x_max = +8

Jx = 200

x_complete = np.linspace(x_min, x_max, Jx+1, endpoint=True)

x = x_complete[0:-1]

dx = x[1] - x[0]




x0 = 0

sigma_0 = 0.5

k0 = 4


u_ref = gaussian(x, 0.0, x0, k0, sigma_0)

u = u_ref




T = 4

dt = 0.001

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)




D2 = finite_differences_1d.get_D2_circulant_2nd_order(Jx, dx)



E = eye(Jx)

A = E - 0.5 * 1j * dt * D2




u_ref = gaussian(x, 0.0, x0, k0, sigma_0)

u = u_ref



n_mod_times_analysis = 50



fig_1 = Figure1(x, 1, 0, None, u_ref)

fig_1.update_u(u, u_ref)

fig_1.redraw()


for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        u_ref = gaussian(x, t, x0, k0, sigma_0)
    
        
        fig_1.update_u(u, u_ref)
        
        fig_1.redraw()
        
        
    u = u + 0.5 * dt * 1j * D2 * u
    
    

input('press any key ...')






