"""
Equation: Gross-Pitaevskii equation
Initial condition: fat soliton
Spatial approximation: Fourier spectral collocation
Boundary conditions: periodic
Time-integration method: Strang-Splitting
"""


from scipy.sparse import spdiags


import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)


from simulations.reference_solutions import Phi


from simulations.figure_2 import Figure2


from differentiation import fourier



x_min = -8
x_max = +8

L = x_max - x_min

Jx = 128

x_complete = np.linspace(x_min, x_max, Jx+1, endpoint=True)

x = x_complete[0:-1]

dx = x[1] - x[0]


mue_x = fourier.get_mue_x(Jx, dx)

mue_squared = mue_x**2 



alpha = 4
nu = 0.9

Phi_of_x = Phi(x, 0.0, alpha, nu)


u_ref = Phi_of_x.copy()
u     = Phi_of_x.copy()


V = -alpha * Phi_of_x

diag_V = spdiags(np.array([V[0:-1]]), np.array([0]), Jx, Jx, 'csr')






density_0 = np.real(u * np.conj(u))
    
N0 = dx * np.sum(density_0)

N0_formula = 4 * np.sqrt(2) * alpha * nu * ( (1.0/(2*nu)) * np.log( (1+nu) / (1-nu) ) - 1)

print(N0)
print(N0_formula)
input()


#------------------------------------------------------------------------------
T = 20

dt = 0.001

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)

dt_new = times[1] - times[0]

assert(dt_new == dt)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
n_mod_times_analysis = 25
# n_mod_times_analysis = 1
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
fig_1 = Figure2(x, 0.0, 4, -8, 2, V, u_ref)

fig_1.update_u(u, u_ref)

fig_1.update_V(V)

fig_1.redraw()
#------------------------------------------------------------------------------



g = 1

hbar = 1
m = 1



compute_ground_state = True

if compute_ground_state == True:

    dt_imag = -1j * dt
    
    exp_mue_squared = np.exp( - 1.0j * dt_imag * mue_squared * hbar / (2.0 * m) )
    
    
    n_iter = 0
    while n_iter < 5000:
        
        if n_iter % n_mod_times_analysis == 0:
                
            rel_error = np.linalg.norm(u-u_ref) / np.linalg.norm(u_ref)
            
            print('n_iter:    {0:d}'.format(n_iter))
            print('rel_error: {0:1.4e}'.format(rel_error))
            print()
            
            # u_ref = Phi(x, 0.0, alpha, nu)
            
            fig_1.update_u(u, u_ref)
            
            fig_1.redraw()
        
        
        density = np.real(u * np.conj(u))
        
        #----------------------------------------------------------------------
        tmp_1 = np.conj(u) * ( -(hbar**2 / (2 * m)) * np.fft.ifft( -mue_squared * np.fft.fft(u) ) )
        E1 = dx * np.sum(tmp_1)
        
        tmp_2 = np.conj(u) * ( V * u )
        E2 = dx * np.sum(tmp_2)
        
        tmp_3 = np.conj(u) * ( (g * density) * u )
        E3 = dx * np.sum(tmp_3)
        #----------------------------------------------------------------------
        
        N = dx * np.sum(density)
        
        mue = np.real((E1 + E2 + E3) / N)
        
        # mue = 0.0
        #======================================================================
        
        #======================================================================
        V_eff = V + g * density - mue
        
        u = np.exp( - 1.0j * V * dt_imag / (2.0 * hbar) ) * u
        #======================================================================
        
        #======================================================================
        u = np.fft.fft(u)
        
        u = exp_mue_squared * u
        
        u = np.fft.ifft(u)
        #======================================================================
        
        #======================================================================
        density = np.real(u * np.conj(u))
        
        V_eff = V + g * density - mue
        
        u = np.exp( - 1.0j * V_eff * dt_imag / (2.0 * hbar) ) * u
        #======================================================================
        
        
        u = np.real(u)
        
        #----------------------------------------------------------------------
        # normalize
        
        N_desired = N0
        
        density = np.real(u * np.conj(u))
        
        N = dx * np.sum(density)
        
        u = np.sqrt(N_desired / N) * u
        #----------------------------------------------------------------------
        
        
        n_iter = n_iter + 1
        
    






exp_mue_squared = np.exp( - 1.0j * dt * mue_squared * hbar / (2.0 * m) )

for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        # u_ref = Phi(x, 0.0, alpha, nu)
        # rel_error = np.linalg.norm(u-u_ref) / np.linalg.norm(u_ref)
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        # print('rel_error: {0:1.4e}'.format(rel_error))
        print()
        
        fig_1.update_u(u, u_ref)
        
        fig_1.redraw()
    
    #==========================================================================
    density = np.real(u * np.conj(u))
    
    """
    #--------------------------------------------------------------------------
    tmp_1 = np.conj(u) * ( -(hbar**2 / (2 * m)) * np.fft.ifft( -mue_squared * np.fft.fft(u) ) )
    E1 = dx * np.sum(tmp_1)
    
    tmp_2 = np.conj(u) * ( V * u )
    E2 = dx * np.sum(tmp_2)
    
    tmp_3 = np.conj(u) * ( (g * density) * u )
    E3 = dx * np.sum(tmp_3)
    #--------------------------------------------------------------------------
    
    N = dx * np.sum(density)
    
    mue = np.real((E1 + E2 + E3) / N)
    """
    
    mue = 0.0
    #==========================================================================
    
    #==========================================================================
    V_eff = V + g * density - mue
    
    u = np.exp( - 1.0j * V * dt / (2.0 * hbar) ) * u
    #==========================================================================
    
    #==========================================================================
    u = np.fft.fft(u)
    
    u = exp_mue_squared * u
    
    u = np.fft.ifft(u)
    #==========================================================================
    
    #==========================================================================
    density = np.real(u * np.conj(u))
    
    V_eff = V + g * density - mue
    
    u = np.exp( - 1.0j * V_eff * dt / (2.0 * hbar) ) * u
    #==========================================================================
    
    
    
    
    
    
