import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)

from linear_schroedinger.wave_packet.reference_solution import gaussian

from linear_schroedinger.wave_packet.figure_1 import Figure1


from differentiation import finite_differences_1d



order_spatial_discretization = 8




x0 = 0

sigma_0 = 0.5

k0 = 4



x_min = -8
x_max = +8

L = x_max - x_min

Jx = 200

x = np.linspace(x_min, x_max, Jx+1, endpoint=True)

dx = x[1] - x[0]












if order_spatial_discretization == 2:
    
    D2 = finite_differences_1d.get_D2_circulant_2nd_order(Jx, dx)


if order_spatial_discretization == 4:
    
    D2 = finite_differences_1d.get_D2_circulant_4th_order(Jx, dx)
    
    
if order_spatial_discretization == 6:
    
    D2 = finite_differences_1d.get_D2_circulant_6th_order(Jx, dx)
    
    
if order_spatial_discretization == 8:
    
    D2 = finite_differences_1d.get_D2_circulant_8th_order(Jx, dx)



A = 1j * 0.5 * D2





if order_spatial_discretization == 2:

    # dt_upper_bound_2nd_order = np.sqrt(2) * dx**2
    dt_upper_bound_2nd_order = 1.414 * dx**2
    
    dt = 1.0 * dt_upper_bound_2nd_order
    
    
if order_spatial_discretization == 4:
    
    # dt_upper_bound_4th_order = (3/4) * np.sqrt(2) * dx**2
    dt_upper_bound_4th_order = 1.060 * dx**2

    dt = 1.0 * dt_upper_bound_4th_order
    
    
if order_spatial_discretization == 6:
    
    dt_upper_bound_6th_order = 0.935 * dx**2

    dt = 1.0 * dt_upper_bound_6th_order
    
      
if order_spatial_discretization == 8:
    
    dt_upper_bound_8th_order = 0.870 * dx**2
    
    dt = 1.0 * dt_upper_bound_8th_order
    



T = 2

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)






u_ref = gaussian(x, 0.0, x0, sigma_0, k0)

u = u_ref[0:-1]

assert(u.size == Jx)


u_complete = np.zeros_like(u_ref)

u_complete[0:-1] = u

u_complete[-1] = u_complete[0]











n_mod_times_analysis = 10

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
        
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        
        u_ref = (
                + gaussian(x - 10*L, t, x0, sigma_0, k0) 
                + gaussian(x -  9*L, t, x0, sigma_0, k0) 
                + gaussian(x -  8*L, t, x0, sigma_0, k0) 
                + gaussian(x -  7*L, t, x0, sigma_0, k0) 
                + gaussian(x -  6*L, t, x0, sigma_0, k0) 
                + gaussian(x -  5*L, t, x0, sigma_0, k0) 
                + gaussian(x -  4*L, t, x0, sigma_0, k0) 
                + gaussian(x -  3*L, t, x0, sigma_0, k0) 
                + gaussian(x -  2*L, t, x0, sigma_0, k0) 
                + gaussian(x -  1*L, t, x0, sigma_0, k0)
                + gaussian(x +  0*L, t, x0, sigma_0, k0)
                + gaussian(x +  1*L, t, x0, sigma_0, k0)
                + gaussian(x +  2*L, t, x0, sigma_0, k0)
                + gaussian(x +  3*L, t, x0, sigma_0, k0)
                + gaussian(x +  4*L, t, x0, sigma_0, k0)
                + gaussian(x +  5*L, t, x0, sigma_0, k0)
                + gaussian(x +  6*L, t, x0, sigma_0, k0)
                + gaussian(x +  7*L, t, x0, sigma_0, k0)
                + gaussian(x +  8*L, t, x0, sigma_0, k0)
                + gaussian(x +  9*L, t, x0, sigma_0, k0)
                + gaussian(x + 10*L, t, x0, sigma_0, k0)
                )
        
        u_complete[0:-1] = u
        u_complete[-1] = u_complete[0]
        
        
        # norm_u_of_times_analysis[nr_times_analysis] = np.linalg.norm(u)
        
        # defect_of_mass_of_times_analysis = np.abs(1.0 - norm_u_of_times_analysis / norm_u_of_times_analysis[0])
        
        # print(defect_of_mass_of_times_analysis[nr_times_analysis])
        
        rel_error_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_complete-u_ref) / np.linalg.norm(u_ref)
        
        times_analysis[nr_times_analysis] = t 
        
        
        fig_1.update_u(u_complete, u_ref)
        
        # fig_1.update_rel_error(rel_error_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.redraw()
        
        
        nr_times_analysis = nr_times_analysis + 1
    
    
    k1 = A * u
    k2 = A * (u + 0.5 * dt * k1)
    k3 = A * (u + 0.5 * dt * k2)
    k4 = A * (u + 1.0 * dt * k3)
    
    u = u + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
    

input('press any key ...')






