import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)

from nonlinear_schroedinger.solitons.bright_soliton.reference_solution import bright_soliton

from nonlinear_schroedinger.solitons.bright_soliton.figure_1 import Figure1


from differentiation import finite_differences_1d


order_spatial_discretization = 8


a = 4
v = 2
x0 = 0
theta_0 = 0
beta = -1







x_min = -4
x_max = +4

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
    
    



def eval_f(y):
    
    return 0.5 * 1j * D2 * y - 1j * beta * np.abs(y)**2 * y






A_linear_part = 0.5 * 1j * D2

eigenvalues_A_linear_part, eigenvectors_A_linear_part = np.linalg.eig(A_linear_part.todense())

max_abs_lambda = np.max(np.abs(eigenvalues_A_linear_part))




dt = 1.0 * np.sqrt(8) / max_abs_lambda

# print(dt)
# input('press any key to continue ... ')







T = 25

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)







u_ref = bright_soliton(x, 0.0, a, v, x0, theta_0, beta)

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
   
fig_1.redraw()




nr_times_analysis = 0

for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        u_ref = (
                + bright_soliton(x+0*L, t, a, v, x0, theta_0, beta) 
                + bright_soliton(x+1*L, t, a, v, x0, theta_0, beta)
                + bright_soliton(x+2*L, t, a, v, x0, theta_0, beta)
                + bright_soliton(x+3*L, t, a, v, x0, theta_0, beta)
                + bright_soliton(x+4*L, t, a, v, x0, theta_0, beta)
                + bright_soliton(x+5*L, t, a, v, x0, theta_0, beta)
                + bright_soliton(x+6*L, t, a, v, x0, theta_0, beta)
                + bright_soliton(x+7*L, t, a, v, x0, theta_0, beta)
                + bright_soliton(x+8*L, t, a, v, x0, theta_0, beta)
                )
        
        u_complete[0:-1] = u
        u_complete[-1] = u_complete[0]
        
        
        norm_u_of_times_analysis[nr_times_analysis] = np.linalg.norm(u)
        
        
        print(norm_u_of_times_analysis[nr_times_analysis] / norm_u_of_times_analysis[0])
        
        defect_of_mass_of_times_analysis = np.abs(1.0 - norm_u_of_times_analysis / norm_u_of_times_analysis[0])
        
        
        
        rel_error_of_times_analysis[nr_times_analysis] = np.linalg.norm(u_complete-u_ref) / np.linalg.norm(u_complete)
        
        times_analysis[nr_times_analysis] = t 
        
        
        fig_1.update_u(u_complete, u_ref)
        
        fig_1.update_rel_error(rel_error_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.update_defect_of_mass(defect_of_mass_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.redraw()
        
        
        nr_times_analysis = nr_times_analysis + 1
    
    
    k1 = eval_f(u)
    k2 = eval_f(u + 0.5 * dt * k1)
    k3 = eval_f(u + 0.5 * dt * k2)
    k4 = eval_f(u + 1.0 * dt * k3)
    
    u = u + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
    

input('press any key ...')






