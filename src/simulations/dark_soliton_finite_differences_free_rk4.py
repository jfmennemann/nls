import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)


from scipy.sparse import diags

from scipy.sparse import csr_matrix



from simulations.reference_solutions import dark_soliton

from simulations.figure_1 import Figure1




v = 1

x0 = -5
theta_0 = 0

u0 = 10
beta = 1
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






D_xx = diags([1, -2, 1], [-1, 0, 1], shape=(Jx+1, Jx+1))

D_xx = D_xx.todense()

D_xx[0, 0] = -1
D_xx[0, 1] = +4
D_xx[0, 2] = -5
D_xx[0, 3] = +2

D_xx[-1, -1] = +2
D_xx[-1, -2] = -5
D_xx[-1, -3] = +4
D_xx[-1, -4] = -1

D_xx = csr_matrix(D_xx)

D_xx = D_xx / dx**2
    





def eval_f(y):
    
    return 0.5 * 1j * D_xx * y - 1j * beta * np.abs(y)**2 * y



u_ref = dark_soliton(x, 0, v, x0, theta_0, u0, beta, phi)

u = u_ref

assert(u.size == Jx+1)





n_mod_times_analysis = 20

times_analysis = times[::n_mod_times_analysis]



norm_u_of_times_analysis = np.zeros_like(times_analysis)
rel_error_of_times_analysis = np.zeros_like(times_analysis)





fig_1 = Figure1(x, 20, 0, None, u_ref)

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
        
        defect_of_mass_of_times_analysis = np.abs(1.0 - norm_u_of_times_analysis / norm_u_of_times_analysis[0])
        
        
        
        rel_error_of_times_analysis[nr_times_analysis] = np.linalg.norm(u-u_ref) / np.linalg.norm(u_ref)
        
        times_analysis[nr_times_analysis] = t 
        
        
        fig_1.update_u(u, u_ref)
        
        # fig_1.update_rel_error(rel_error_of_times_analysis, times_analysis, nr_times_analysis)
        
        # fig_1.update_defect_of_mass(defect_of_mass_of_times_analysis, times_analysis, nr_times_analysis)
        
        fig_1.redraw()
        
        
        nr_times_analysis = nr_times_analysis + 1
    
    
    k1 = eval_f(u)
    k2 = eval_f(u + 0.5 * dt * k1)
    k3 = eval_f(u + 0.5 * dt * k2)
    k4 = eval_f(u + 1.0 * dt * k3)
    
    u = u + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
    
input('press any key ...')






