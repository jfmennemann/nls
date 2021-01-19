import numpy as np

from numpy import pi

np.set_printoptions(edgeitems=8, linewidth=200, precision=10)



x_min = -pi
x_max = +pi

J = 4

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

# print(vec_nue_1st_derivative)
# print(vec_nue_2nd_derivative)
# print()
# input('press any key ...')
# print()


vec_lambda_1st_derivative = ( 1j * (2 * pi / L) * vec_nue_1st_derivative )**1
vec_lambda_2nd_derivative = ( 1j * (2 * pi / L) * vec_nue_2nd_derivative )**2







from scipy.linalg import circulant

c = np.zeros((J,))

c[0] = -pi**2/(3*dx**2) - 1.0/6.0

for j in np.arange(1,J):
    
    c[j] = -(-1)**j / ( 2 * np.sin(j * dx / 2)**2 )


D2_ref = circulant(c)





F     = np.zeros((J, J), dtype=np.complex)
F_inv = np.zeros((J, J), dtype=np.complex)

for k in np.arange(J):
    
    for j in np.arange(J):
        
        F[k,j]     = np.exp(-1j * 2 * pi * k * j / J)
        F_inv[j,k] = np.exp(+1j * 2 * pi * k * j / J) / J




Lambda = np.diag(vec_lambda_2nd_derivative)

tmp_1 = np.matmul(Lambda, F)
D2 = np.matmul(F_inv, tmp_1)




print(D2_ref)
print(D2)
print()


rel_error = np.linalg.norm(D2 - D2_ref) / np.linalg.norm(D2_ref)

print(rel_error)



















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

"""
print('u_d_imag_max:  {0:1.4e}'.format(u_d_imag_max))
print('u_dd_imag_max: {0:1.4e}'.format(u_dd_imag_max))
print()
print('error_u_d:     {0:1.4e}'.format(error_u_d))
print('error_u_dd:    {0:1.4e}'.format(error_u_dd))
"""




u_d  = np.real(u_d)
u_dd = np.real(u_dd)




    