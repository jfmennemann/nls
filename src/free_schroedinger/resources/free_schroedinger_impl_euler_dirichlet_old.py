import matplotlib.pyplot as plt


from scipy.sparse import spdiags, eye

from scipy.sparse.linalg import spsolve


import numpy as np




x_min = -10
x_max = +10

Jx = 250

x = np.linspace(x_min, x_max, Jx, endpoint=True)

dx = x[1] - x[0]



T = 10

dt = 0.001

n_times = np.int(np.round(T / dt)) + 1
        
times = np.linspace(0, T, n_times, endpoint=True)



d_minus_1 = +1 * np.ones(Jx)
d_0       = -2 * np.ones(Jx)
d_plus_1  = +1 * np.ones(Jx)

D_xx = spdiags([d_minus_1, d_0, d_plus_1], [-1,0,1], Jx, Jx, 'csr') / dx**2

E = eye(Jx)


A = E - 0.5 * 1j * dt * D_xx






def gaussian(x, t, x0, sigma_0, k0):
    
    tau = 2 * sigma_0**2
    
    alpha = 1.0 + 1j * t / tau
    
    return (1.0/np.sqrt(alpha)) * np.exp( (1.0/alpha) * ( -((x-x0)/(2*sigma_0))**2 + 1j * k0 * (x-x0) - 1j * sigma_0**2 * k0**2 * t/ tau) )





x0 = 0.0

sigma_0 = 1.0

k0 = 2


u_ref = gaussian(x, 0.0, x0, sigma_0, k0)

u = u_ref




n_mod_times_analysis = 25

times_analysis = times[::n_mod_times_analysis]




rel_error_of_times_analysis = np.zeros_like(times_analysis)


print(times[-1])
print(times_analysis[-1])



plt.ion()

#======================================================================================
fig_temp = plt.figure("figure_1", figsize=(8, 8), facecolor="white")

#--------------------------------------------------------------------------------------
ax_1 = fig_temp.add_subplot(411)

line_u_abs_squared,     = ax_1.plot(x, np.abs(u)**2,     linewidth=1.0, linestyle='-', color='k',         label='u')
line_u_ref_abs_squared, = ax_1.plot(x, np.abs(u_ref)**2, linewidth=1.0, linestyle='-', color='tab:green', label='u_ref')

ax_1.set_ylim(-0.1, 1.1)

ax_1.set_xlabel('x')
ax_1.set_ylabel('|u|^2')

ax_1.legend(loc='upper right', fancybox=False, ncol=2)
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
ax_2 = fig_temp.add_subplot(412)

line_u_real,     = ax_2.plot(x, np.real(u),     linewidth=1.0, linestyle='-', color='k', label='u')
line_u_ref_real, = ax_2.plot(x, np.real(u_ref), linewidth=1.0, linestyle='-', color='tab:green', label='u_ref')

ax_2.set_ylim(-1.1, 1.1)

ax_2.set_xlabel('x')
ax_2.set_ylabel('Re u')

ax_2.legend(loc='upper right', fancybox=False, ncol=2)
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
ax_3 = fig_temp.add_subplot(413)

line_u_imag,     = ax_3.plot(x, np.imag(u),     linewidth=1.0, linestyle='-', color='k', label='u')
line_u_ref_imag, = ax_3.plot(x, np.imag(u_ref), linewidth=1.0, linestyle='-', color='tab:green', label='u_ref')

ax_3.set_ylim(-1.1, 1.1)

ax_3.set_xlabel('x')
ax_3.set_ylabel('Im u')

ax_3.legend(loc='upper right', fancybox=False, ncol=2)
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
ax_4 = fig_temp.add_subplot(414)

line_rel_error, = ax_4.plot([], [], linewidth=1.0, linestyle='-', color='k')

ax_4.set_xlabel('t')
ax_4.set_ylabel('rel_error')

ax_4.set_xlim(0, T)
ax_4.set_ylim(-0.1, 2.1)
#--------------------------------------------------------------------------------------



plt.tight_layout()

plt.draw()
fig_temp.canvas.start_event_loop(0.001)
#======================================================================================




nr_times_analysis = 0

for n in np.arange(times.size):
    
    t = times[n]
    
    print('n: {0:d}'.format(n))
    print('t: {0:1.2f}'.format(t))
    print()
    
    u_ref = gaussian(x, t, x0, sigma_0, k0)
    
    times_analysis[nr_times_analysis] = t 
    rel_error_of_times_analysis[nr_times_analysis] = np.linalg.norm(u-u_ref) / np.linalg.norm(u)
    
    if n % n_mod_times_analysis == 0:
    
        line_u_abs_squared.set_ydata(np.abs(u)**2)
        line_u_ref_abs_squared.set_ydata(np.abs(u_ref)**2)
        
        line_u_real.set_ydata(np.imag(u))
        line_u_ref_real.set_ydata(np.imag(u_ref))
        
        line_u_imag.set_ydata(np.real(u))
        line_u_ref_imag.set_ydata(np.real(u_ref))
        
        line_rel_error.set_xdata(times_analysis[0:nr_times_analysis])
        line_rel_error.set_ydata(rel_error_of_times_analysis[0:nr_times_analysis])
        
        plt.tight_layout()
        
        plt.draw()
        fig_temp.canvas.start_event_loop(0.001)
        
        nr_times_analysis = nr_times_analysis + 1
    
    u = spsolve(A, u)
    


plt.ioff()

input('press any key ...')





