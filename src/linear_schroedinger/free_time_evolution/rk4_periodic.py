from scipy.sparse import diags

import numpy as np

from linear_schroedinger.free_time_evolution.figure_1 import Figure1



def gaussian(x, t, x0, sigma_0, k0):
    
    tau = 2 * sigma_0**2
    
    alpha = 1.0 + 1j * t / tau
    
    return (1.0/np.sqrt(alpha)) * np.exp( (1.0/alpha) * ( -((x-x0)/(2*sigma_0))**2 + 1j * k0 * (x-x0) - 1j * sigma_0**2 * k0**2 * t/ tau) )




order_spatial_discretization = 2


x0 = -2.5

sigma_0 = 0.5

k0 = 4



x_min = -10
x_max = +10

Jx = 200

x = np.linspace(x_min, x_max, Jx+1, endpoint=True)

dx = x[1] - x[0]



T = 5

# dt = 0.001

dt = 0.015
# dt = 0.014

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)






if order_spatial_discretization == 2:
    
    D_xx = diags([1, -2, 1], [-1, 0, 1], shape=(Jx, Jx))
    
    D_xx = D_xx / dx**2
    
    a = 0.5
    
if order_spatial_discretization == 4:
    
    D_xx = diags([-1, 16, -30, 16, -1], [-2, -1, 0, 1, 2], shape=(Jx, Jx))

    D_xx = D_xx / (12 * dx**2)
    
if order_spatial_discretization == 6:
    
    D_xx = diags([2, -27, 270, -490, 270, -27, 2], [-3, -2, -1, 0, 1, 2, 3], shape=(Jx, Jx))
    
    D_xx = D_xx / (180 * dx**2)
    
if order_spatial_discretization == 8:
    
    D_xx = diags([-9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9], [-4, -3, -2, -1, 0, 1, 2, 3, 4], shape=(Jx, Jx))
    
    D_xx = D_xx / (5040 * dx**2)


A = 1j * 0.5 * D_xx


dt_upper_bound = (np.sqrt(8) / 4) * dx**2 / a

print(dt_upper_bound)
input('press any key ...')



u_ref = gaussian(x, 0.0, x0, sigma_0, k0)

u = u_ref[0:-1]

assert(u.size == Jx)


u_complete = np.zeros_like(u_ref)

u_complete[0:-1] = u

u_complete[-1] = u_complete[0]











n_mod_times_analysis = 1

times_analysis = times[::n_mod_times_analysis]



norm_u_of_times_analysis = np.zeros_like(times_analysis)
rel_error_of_times_analysis = np.zeros_like(times_analysis)




fig_1 = Figure1(x, times, screen_size='large')

fig_1.update_u(u_complete, u_ref)

  

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
        
        u_ref = gaussian(x, t, x0, sigma_0, k0)
        
        u_complete[0:-1] = u
        u_complete[-1] = u_complete[0]
        
        
        # norm_u_of_times_analysis[nr_times_analysis] = np.linalg.norm(u)
        
        # defect_of_mass_of_times_analysis = np.abs(1.0 - norm_u_of_times_analysis / norm_u_of_times_analysis[0])
        
        # print(defect_of_mass_of_times_analysis[nr_times_analysis])
        
        rel_error_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_complete-u_ref) / np.linalg.norm(u_ref)
        
        times_analysis[nr_times_analysis] = t 
        
        
        fig_1.update_u(u_complete, u_ref)
        
        fig_1.update_rel_error(rel_error_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.redraw()
        
        
        nr_times_analysis = nr_times_analysis + 1
    
    
    k1 = A * u
    k2 = A * (u + 0.5 * dt * k1)
    k3 = A * (u + 0.5 * dt * k2)
    k4 = A * (u + 1.0 * dt * k3)
    
    u = u + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
    

input('press any key ...')






