from numpy import exp
from numpy import cosh
from numpy import tanh
from numpy import cos
from numpy import sqrt

from numpy import pi


def gaussian(x, t, x0, sigma_0, k0):
    
    tau = 2 * sigma_0**2
    
    alpha = 1.0 + 1j * t / tau
    
    return (1.0/sqrt(alpha)) * exp( (1.0/alpha) * ( -((x-x0)/(2*sigma_0))**2 + 1j * k0 * (x-x0) - 1j * sigma_0**2 * k0**2 * t/tau) )



def coherent_state(x, t, x0, omega):
    
    return (omega/pi)**0.25 * exp( -(omega/2) * ( x**2 - 2*x*x0*exp(-1j*omega*t) + (x0**2/2.0) * exp(-2*1j*omega*t) + x0**2/2.0 ) - 1j * omega * t / 2.0 )



def bright_soliton(x, t, a, v, x0, theta_0, beta):
    
    assert(beta < 0)
    
    return (a / sqrt(-beta)) * ( 1.0 / cosh( a * (x - v*t - x0) ) ) * exp( 1j * ( v*x - 0.5*(v**2-a**2)*t + theta_0 ) )



def dark_soliton(x, t, v, x0, theta_0, u0, beta, phi):
    
    assert(beta > 0)
    assert(u0 > 0)
    
    a = sqrt(u0 * beta) * cos(phi)
    
    return (1 / sqrt(beta)) * (a * tanh(a*(x-v*t-x0)) + 1j*v) * exp(1j*(-(a**2 + v**2)*t + theta_0 ) )