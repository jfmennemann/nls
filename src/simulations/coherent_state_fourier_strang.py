"""
Equation: linear Schroedinger equation
Initial condition: coherent state
Spatial approximation: Fourier spectral collocation
Boundary conditions: periodic
Time-integration method: Strang-Splitting
"""


from scipy.sparse import spdiags


import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)


from simulations.reference_solutions import coherent_state

from simulations.figure_1 import Figure1


from differentiation import fourier



x0 = 1

omega = 5


x_min = -3
x_max = +3

L = x_max - x_min

Jx = 200

x_complete = np.linspace(x_min, x_max, Jx+1, endpoint=True)

x = x_complete[0:-1]

dx = x[1] - x[0]




V = 0.5 * omega**2 * x**2

diag_V = spdiags(np.array([V[0:-1]]), np.array([0]), Jx, Jx, 'csr')





#------------------------------------------------------------------------------
T = 2

dt = 0.0005

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)

dt_new = times[1] - times[0]

assert(dt_new == dt)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
n_mod_times_analysis = 25
# n_mod_times_analysis = 1
#------------------------------------------------------------------------------



u_ref = coherent_state(x, 0.0, x0, omega)

u = u_ref



fig_1 = Figure1(x, 1.5, 200, V, u_ref)

fig_1.update_u(u, u_ref)

fig_1.update_V(V)

fig_1.redraw()






hbar = 1
m = 1

mue_x = fourier.get_mue_x(Jx, dx)

mue_squared = mue_x**2

exp_mue_squared = np.exp( - 1.0j * dt * mue_squared * hbar / (2.0 * m) )


for n in np.arange(times.size):
    
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        print('n: {0:d}'.format(n))
        print('t: {0:1.2f}'.format(t))
        print()
        
        u_ref = coherent_state(x, t, x0, omega)
        
        fig_1.update_u(u, u_ref)
        
        fig_1.redraw()
    
    
    #======================================================================
    density = np.real(u * np.conj(u))
    
    #----------------------------------------------------------------------
    tmp_1 = np.conj(u) * ( -(hbar**2 / (2 * m)) * np.fft.ifft( -mue_squared * np.fft.fft(u) ) )
    E1 = dx * np.sum(tmp_1)
    
    tmp_2 = np.conj(u) * ( V * u )
    E2 = dx * np.sum(tmp_2)
    
    # tmp_3 = np.conj(u) * ( (g * density) * u )
    # E3 = dx * np.sum(tmp_3)
    #----------------------------------------------------------------------
    
    N = dx * np.sum(density)
    
    # mue_batch = np.real((E1 + E2 + E3) / N)
    mue = np.real((E1 + E2) / N)
    #======================================================================
    
    #======================================================================
    # V_eff = V + g * density - mue
    V_eff = V
    
    u = np.exp( - 1.0j * V * dt / (2.0 * hbar) ) * u
    #======================================================================
    
    #======================================================================
    
    u = np.fft.fft(u)
    
    u = exp_mue_squared * u
    
    u = np.fft.ifft(u)
    
    #======================================================================
    
    #======================================================================
    density = np.real(u * np.conj(u))
    
    # V_eff = V + g * density - mue
    V_eff = V
    
    u = np.exp( - 1.0j * V_eff * dt / (2.0 * hbar) ) * u
    #======================================================================
    #========================================================================================== 
    
    
    

input('press any key ...')






