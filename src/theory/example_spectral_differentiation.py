import numpy as np

from numpy import pi



x_min = -pi
x_max = +pi

J = 8

assert(J % 2 == 0)

x = np.linspace(x_min, x_max, J, endpoint=False)

index_center_x = J//2

assert(np.abs(x[index_center_x]) < 1e-14)

dx = x[1] - x[0]

L = J * dx
             
vec_nue_1st_derivative = np.arange(-J//2, J//2)
vec_nue_1st_derivative[0] = 0

vec_nue_2nd_derivative = np.arange(-J//2, J//2)

vec_nue_1st_derivative = np.fft.fftshift(vec_nue_1st_derivative)
vec_nue_2nd_derivative = np.fft.fftshift(vec_nue_2nd_derivative)

print(vec_nue_1st_derivative)
print(vec_nue_2nd_derivative)
print()
# input('press any key ...')
# print()


vec_lambda_1st_derivative = ( 1j * (2 * pi / L) * vec_nue_1st_derivative )**1
vec_lambda_2nd_derivative = ( 1j * (2 * pi / L) * vec_nue_2nd_derivative )**2



#==============================================================================
u = np.exp(np.sin(x))

u_d_ref  = np.cos(x) * u
u_dd_ref = -np.sin(x) * u + np.cos(x) * u_d_ref
#==============================================================================



u_d  = np.fft.ifftn(vec_lambda_1st_derivative * np.fft.fftn(u))
u_dd = np.fft.ifftn(vec_lambda_2nd_derivative * np.fft.fftn(u))



u_d_imag_max  = np.max(np.imag(u_d))
u_dd_imag_max = np.max(np.imag(u_dd))

error_u_d  = np.linalg.norm(u_d  - u_d_ref,  np.inf)
error_u_dd = np.linalg.norm(u_dd - u_dd_ref, np.inf)


print('u_d_imag_max:  {0:1.4e}'.format(u_d_imag_max))
print('u_dd_imag_max: {0:1.4e}'.format(u_dd_imag_max))
print()
print('error_u_d:     {0:1.4e}'.format(error_u_d))
print('error_u_dd:    {0:1.4e}'.format(error_u_dd))







u_d  = np.real(u_d)
u_dd = np.real(u_dd)


from style_sheet import mystyle

import matplotlib as mpl
import matplotlib.pyplot as plt

export_pdf = False

if export_pdf == True:

    mpl.use("pgf")
    
    plt.rcParams.update(mystyle.rc_parameters)
    


linestyle_u_ref = mystyle.linestyle_u_ref
linestyle_u     = mystyle.linestyle_u

linewidth_u_ref = mystyle.linewidth_u_ref
linewidth_u     = mystyle.linewidth_u

color_u_ref = mystyle.color_u_ref
color_u     = mystyle.color_u

label_u_ref = r'$u_\mathrm{ref}$'
label_u     = r'$u$'

color_gridlines_major = mystyle.color_gridlines_major
color_gridlines_minor = mystyle.color_gridlines_minor

linestyle_gridlines_major = mystyle.linestyle_gridlines_major
linestyle_gridlines_minor = mystyle.linestyle_gridlines_minor

linewidth_gridlines_major = mystyle.linewidth_gridlines_major
linewidth_gridlines_minor = mystyle.linewidth_gridlines_minor

x_ticks_major = np.array([-8, -4, 0, 4, 8])
x_ticks_minor = np.array([-7, -6, -5, -3, -2, -1, 1, 2, 3, 5, 6, 7])
            
y_ticks_major_left_column = np.array([0, 0.5, 1.0])
y_ticks_minor_left_column = np.array([0.25, 0.75])

y_ticks_major_right_column = np.array([-1, 0, 1])
y_ticks_minor_right_column = np.array([-0.5, 0.5])

x_min = -0.5
x_max = +2*pi+0.5

y_min_left_column = -1.1
y_max_left_column = +1.1

y_min_right_column = -1.2
y_max_right_column = +1.2


xlabel = r'$x$'



ylabel_00 = r'$u$'
ylabel_10 = r'$u_prime$'




width  = 10
height = 10

fig = plt.figure("figure_wave_packet_dirichlet", figsize=(width, height), facecolor="white", constrained_layout=False)

spacing_x = 0.2
spacing_y = 0.1

gridspec = fig.add_gridspec(ncols=1, nrows=3, left=0.085, right=0.99, bottom=0.05, top=0.99, wspace=spacing_x, hspace=spacing_y)



#==========================================================================================
ax_00 = fig.add_subplot(gridspec[0, 0])
ax_10 = fig.add_subplot(gridspec[1, 0])
ax_20 = fig.add_subplot(gridspec[2, 0])
#==========================================================================================

#==========================================================================================
ax_00.plot(x, u, linewidth=linewidth_u, linestyle=linestyle_u, color=color_u, label=r'$u$')

# ax_00.set_xlim(x_min, x_max)
# ax_00.set_ylim(y_min_left_column, y_max_left_column)

# ax_00.set_xticks(x_ticks_major, minor=False)
# ax_00.set_xticks(x_ticks_minor, minor=True)

# ax_00.set_yticks(y_ticks_major_left_column, minor=False)
# ax_00.set_yticks(y_ticks_minor_left_column, minor=True)

ax_00.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_00.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_00.set_ylabel(r'$u$')

ax_00.set_xticklabels([])

# ax_00.legend(loc='upper left', ncol=2)
#==========================================================================================

#==========================================================================================
ax_10.plot(x, u_d_ref, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=r'$u_\mathrm{d}$')
ax_10.plot(x, u_d,     linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

# ax_00.set_xlim(x_min, x_max)
# ax_00.set_ylim(y_min_left_column, y_max_left_column)

# ax_00.set_xticks(x_ticks_major, minor=False)
# ax_00.set_xticks(x_ticks_minor, minor=True)

# ax_00.set_yticks(y_ticks_major_left_column, minor=False)
# ax_00.set_yticks(y_ticks_minor_left_column, minor=True)

ax_10.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_10.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_10.set_ylabel(r'$\frac{du}{dx}(x)$')

ax_10.set_xticklabels([])

# ax_00.legend(loc='upper left', ncol=2)
#==========================================================================================

#==========================================================================================
ax_20.plot(x, u_dd_ref, linewidth=linewidth_u_ref, linestyle=linestyle_u_ref, color=color_u_ref, label=r'$u_\mathrm{dd}$')
ax_20.plot(x, u_dd,     linewidth=linewidth_u,     linestyle=linestyle_u,     color=color_u,     label=label_u)

# ax_00.set_xlim(x_min, x_max)
# ax_00.set_ylim(y_min_left_column, y_max_left_column)

# ax_00.set_xticks(x_ticks_major, minor=False)
# ax_00.set_xticks(x_ticks_minor, minor=True)

# ax_00.set_yticks(y_ticks_major_left_column, minor=False)
# ax_00.set_yticks(y_ticks_minor_left_column, minor=True)

ax_20.grid(b=True, which='major', color=color_gridlines_major, linestyle=linestyle_gridlines_major, linewidth=linewidth_gridlines_major)
ax_20.grid(b=True, which='minor', color=color_gridlines_minor, linestyle=linestyle_gridlines_minor, linewidth=linewidth_gridlines_minor)

ax_20.set_xlabel(r'$x$')
ax_20.set_ylabel(r'$\frac{d^2 u}{dx^2}(x)$')

# ax_20.set_xticklabels([])

# ax_00.legend(loc='upper left', ncol=2)
#==========================================================================================


if export_pdf == True:

    path = "/home/jfmennemann/git/nls/pdf/"
    
    filepath = path + "figure_finite_difference_periodic_crank_nicolson.pdf"

    plt.savefig(filepath, backend='pgf')

else:

    plt.show()



  
    