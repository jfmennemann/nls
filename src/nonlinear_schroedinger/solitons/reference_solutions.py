import numpy as np

from numpy import exp
from numpy import cos
from numpy import tanh
from numpy import sqrt

def bright_soliton(x, t, a, v, x0, theta_0, beta):
    
    assert(beta < 0)
    
    return (a / sqrt(-beta)) * ( 1.0 / np.cosh( a * (x - v*t - x0) ) ) * np.exp( 1j * ( v*x - 0.5*(v**2-a**2)*t + theta_0 ) )
    
    
def dark_soliton(x, t, v, x0, theta_0, u0, beta, phi):
    
    assert(beta > 0)
    assert(u0 > 0)
    
    a = sqrt(u0 * beta) * cos(phi)
    
    return (1 / sqrt(beta)) * (a * tanh(a*(x-v*t-x0)) + 1j*v) * exp(1j*(-(a**2 + v**2)*t + theta_0 ) )