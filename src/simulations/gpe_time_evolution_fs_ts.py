"""
Equation: time-dependent Gross-Pitaevskii equation
Initial condition: stationary solution
Spatial approximation: Fourier spectral collocation
Boundary conditions: periodic
Time-integration method: time-splitting
"""


import numpy as np

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)


from simulations.figure_2 import Figure2



hbar = 1
m = 1


x_min = -8
x_max = +8

L = x_max - x_min

Jx = 128

x_complete = np.linspace(x_min, x_max, Jx+1, endpoint=True)

x = x_complete[0:-1]

dx = x[1] - x[0]

Lx = Jx * dx
   
nue = np.arange(-Jx//2, Jx//2)
# nue[0] = 0

nue = np.fft.fftshift(nue)
    
nue_x = (2 * np.pi / Lx) * nue

nue_x_squared = nue_x**2 




nr_example = 1

if nr_example == 1:
    
    mue_phi = 1
    g = 1
    
    alpha = 1
    
    phi = np.exp(-alpha*x**2/2)
    
    phi_dd = -alpha * phi + alpha**2 * x**2 * phi 
    
    """
    phi_dd_approx = np.fft.ifftn(-nue_x_squared * np.fft.fftn(phi))
    
    rel_error_phi_dd = np.linalg.norm(phi_dd_approx-phi_dd) / np.linalg.norm(phi_dd)
    
    print(rel_error_phi_dd)
    input()
    """
    
    V = mue_phi + 0.5 * (-alpha + alpha**2 * x**2) - g * np.abs(phi)**2
    
    #----------------------------------------------------------------------------------------------
    # settings for the figure
    
    density_min_tmp = 0.0
    density_max_tmp = 1.25
    
    V_min_tmp = 0
    V_max_tmp = 20
    
    density_min = density_min_tmp - 0.1 * (density_max_tmp - density_min_tmp)
    density_max = density_max_tmp + 0.1 * (density_max_tmp - density_min_tmp)
    
    V_min = V_min_tmp - 0.1 * (V_max_tmp - V_min_tmp)
    V_max = V_max_tmp + 0.1 * (V_max_tmp - V_min_tmp)
    #----------------------------------------------------------------------------------------------


elif nr_example == 2:
    
    # in the example of the article g equals 1
    g = 1
    
    from numpy import log, sqrt, tanh
    
    def eval_phi_article(x, alpha, nu):
        
        # the factor sqrt(2) before x is needed to compensate the fact that we consider 1/2 d_xx u
        # instead of d_xx u as in the article
        
        xi = sqrt(2) * x 
        
        theta = 0.25 * log( (1+nu) / (1-nu) )
        
        Delta = 3 * sqrt(2) / (alpha * nu)
        
        phi = (alpha * nu / 3) * ( tanh( (xi/Delta) + theta ) - tanh( (xi/Delta) - theta ) )
    
        return phi
    
    alpha = 4
    nu = 0.95
    
    phi = eval_phi_article(x, alpha, nu)
    
    V = -alpha * phi
    
    #----------------------------------------------------------------------------------------------
    # compute chemical potential of phi
    
    density_phi = np.real(phi * np.conj(phi))
            
    tmp_1 = np.conj(phi) * ( -(hbar**2 / (2 * m)) * np.fft.ifft( -nue_x_squared * np.fft.fft(phi) ) )
    E1 = dx * np.sum(tmp_1)
    
    tmp_2 = np.conj(phi) * ( V * phi )
    E2 = dx * np.sum(tmp_2)
    
    tmp_3 = np.conj(phi) * ( (g * density_phi) * phi )
    E3 = dx * np.sum(tmp_3)
    
    N = dx * np.sum(density_phi)
    
    mue_phi = np.real((E1 + E2 + E3) / N)
    #----------------------------------------------------------------------------------------------
    
    #----------------------------------------------------------------------------------------------
    # settings for the figure
    
    density_min_tmp = 0.0
    density_max_tmp = 4.0
    
    V_min_tmp = -8
    V_max_tmp =  0
    
    density_min = density_min_tmp - 0.1 * (density_max_tmp - density_min_tmp)
    density_max = density_max_tmp + 0.1 * (density_max_tmp - density_min_tmp)
    
    V_min = V_min_tmp - 0.1 * (V_max_tmp - V_min_tmp)
    V_max = V_max_tmp + 0.1 * (V_max_tmp - V_min_tmp)
    #----------------------------------------------------------------------------------------------
    
    
    
    

#------------------------------------------------------------------------------
density_phi = np.real(phi * np.conj(phi))
N_phi = dx * np.sum(density_phi)
#------------------------------------------------------------------------------





u_ref = phi
u     = phi




#------------------------------------------------------------------------------
T = 20

dt = 0.01

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
fig_1 = Figure2(x, density_min, density_max, V_min, V_max, V, u_ref)

fig_1.update_u(u, u_ref)

fig_1.update_V(V)

fig_1.redraw()
#------------------------------------------------------------------------------


exp_nue_squared = np.exp( - 1.0j * dt * nue_x_squared * hbar / (2.0 * m) )

for n in np.arange(times.size):
        
    if n % n_mod_times_analysis == 0:
        
        t = times[n]
        
        u_ref = phi * np.exp(-1j * t * mue_phi / hbar)
        
        rel_error = np.linalg.norm(u-u_ref) / np.linalg.norm(u_ref)
        
        print('n:         {0:d}'.format(n))
        print('t:         {0:1.2f}'.format(t))
        print('rel_error: {0:1.4e}'.format(rel_error))
        print()
        
        fig_1.update_u(u, u_ref)
        
        fig_1.redraw()
    
    #==========================================================================
    density = np.real(u * np.conj(u))
    
    V_eff = V + g * density
    
    u = np.exp( - 1.0j * V_eff * dt / (2.0 * hbar) ) * u
    #==========================================================================
    
    #==========================================================================
    u = np.fft.fft(u)
    
    u = exp_nue_squared * u
    
    u = np.fft.ifft(u)
    #==========================================================================
    
    #==========================================================================
    density = np.real(u * np.conj(u))
    
    V_eff = V + g * density
    
    u = np.exp( - 1.0j * V_eff * dt / (2.0 * hbar) ) * u
    #==========================================================================
    
    
    
    
    
    


    




