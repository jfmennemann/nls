import matplotlib.pyplot as plt

from scipy.sparse import spdiags, eye
# from scipy.sparse import csr_matrix

from scipy.sparse.linalg import spsolve


import numpy as np

np.set_printoptions(edgeitems=8, linewidth=100, precision=8)



x_min = -50
x_max = +50

Jx = 100

x = np.linspace(x_min, x_max, Jx, endpoint=True)

dx = x[1] - x[0]



T = 1000

dt = 0.01

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)



d_minus_1 = +1 * np.ones(Jx)
d_0       = -2 * np.ones(Jx)
d_plus_1  = +1 * np.ones(Jx)

D_xx = spdiags([d_minus_1, d_0, d_plus_1], [-1,0,1], Jx, Jx, 'csr') / dx**2

E = eye(Jx)


A = E + 0.5 * 1j * dt * D_xx






def gaussian(x, t, x0=0, sigma_0=1.0):
    
    tau = 2 * sigma_0**2
    
    alpha = 1 + 1j * t / tau

    return (1.0/np.sqrt(alpha)) * np.exp(-(1.0/alpha) * ((x-x0)/(2*sigma_0))**2 )





x0 = 0.0

sigma_0 = 2.0

k0 = 2.0


u_ref = gaussian(x, 0.0, x0=x0, sigma_0=sigma_0)

u = u_ref




plt.ion()

#======================================================================================
fig_temp = plt.figure("figure_1", figsize=(6, 3), facecolor="white")

#--------------------------------------------------------------------------------------
ax_1 = fig_temp.add_subplot(111)

line_u_abs_squared,     = ax_1.plot(x, np.abs(u)**2,     linewidth=1.0, linestyle='-', color='k')
line_u_ref_abs_squared, = ax_1.plot(x, np.abs(u_ref)**2, linewidth=1.0, linestyle='-', color='tab:green')

ax_1.set_xlabel('x')
ax_1.set_ylabel('|u|^2')

ax_1.set_ylim(-0.1, 1.1)
#--------------------------------------------------------------------------------------

plt.tight_layout()

plt.draw()
fig_temp.canvas.start_event_loop(0.001)
#======================================================================================


for n in np.arange(times.size):
    
    t = times[n]
    
    print(n)
    
    if n % 10 == 0:
    
        u_ref = gaussian(x, t, x0=x0, sigma_0=sigma_0)
    
        line_u_abs_squared.set_ydata(np.abs(u)**2)
        line_u_ref_abs_squared.set_ydata(np.abs(u_ref)**2)
        
        plt.tight_layout()
        
        plt.draw()
        fig_temp.canvas.start_event_loop(0.001)
    
    u = spsolve(A, u)
    




