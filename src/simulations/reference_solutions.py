from numpy import exp
from numpy import log
from numpy import cosh
from numpy import tanh
from numpy import cos
from numpy import sqrt

from numpy import pi


def Phi(x, t, alpha, nu):
    
    xi = x
    
    theta = 0.25 * log( (1+nu) / (1-nu) )
    
    Delta = 3 * sqrt(2) / (alpha * nu)
    
    Phi = (alpha * nu / 3) * ( tanh( (xi/Delta) + theta ) - tanh( (xi/Delta) - theta ) )

    return Phi



def gaussian(x, t, x0, k0, sigma_0):
    
    tau = 2 * sigma_0**2
    
    alpha = 1.0 + 1j * t / tau
    
    return (1.0/sqrt(alpha)) * exp( (1.0/alpha) * ( -((x-x0)/(2*sigma_0))**2 + 1j * k0 * (x-x0) - 1j * sigma_0**2 * k0**2 * t/tau) )


def gaussian_periodic(x, t, x0, k0, sigma_0, L):
    
    u = (
          gaussian(x - 10*L, t, x0, k0, sigma_0) 
        + gaussian(x -  9*L, t, x0, k0, sigma_0) 
        + gaussian(x -  8*L, t, x0, k0, sigma_0) 
        + gaussian(x -  7*L, t, x0, k0, sigma_0) 
        + gaussian(x -  6*L, t, x0, k0, sigma_0) 
        + gaussian(x -  5*L, t, x0, k0, sigma_0) 
        + gaussian(x -  4*L, t, x0, k0, sigma_0) 
        + gaussian(x -  3*L, t, x0, k0, sigma_0) 
        + gaussian(x -  2*L, t, x0, k0, sigma_0) 
        + gaussian(x -  1*L, t, x0, k0, sigma_0)
        + gaussian(x +  0*L, t, x0, k0, sigma_0)
        + gaussian(x +  1*L, t, x0, k0, sigma_0)
        + gaussian(x +  2*L, t, x0, k0, sigma_0)
        + gaussian(x +  3*L, t, x0, k0, sigma_0)
        + gaussian(x +  4*L, t, x0, k0, sigma_0)
        + gaussian(x +  5*L, t, x0, k0, sigma_0)
        + gaussian(x +  6*L, t, x0, k0, sigma_0)
        + gaussian(x +  7*L, t, x0, k0, sigma_0)
        + gaussian(x +  8*L, t, x0, k0, sigma_0)
        + gaussian(x +  9*L, t, x0, k0, sigma_0)
        + gaussian(x + 10*L, t, x0, k0, sigma_0)
        )
    
    return u


def coherent_state(x, t, x0, omega):
    
    return (omega/pi)**0.25 * exp( -(omega/2) * ( x**2 - 2*x*x0*exp(-1j*omega*t) + (x0**2/2.0) * exp(-2*1j*omega*t) + x0**2/2.0 ) - 1j * omega * t / 2.0 )




def bright_soliton(x, t, x0, theta_0, a, v, beta):
    
    assert(beta < 0)
    
    return (a / sqrt(-beta)) * ( 1.0 / cosh( a * (x - v*t - x0) ) ) * exp( 1j * ( v*x - 0.5*(v**2-a**2)*t + theta_0 ) )


def bright_soliton_periodic(x, t, x0, theta_0, a, v, beta, L):
    
    assert(beta < 0)
        
    u = (
          bright_soliton(x - 10*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x -  9*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x -  8*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x -  7*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x -  6*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x -  5*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x -  4*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x -  3*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x -  2*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x -  1*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x +  0*L, t, x0, theta_0, a, v, beta) 
        + bright_soliton(x +  1*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x +  2*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x +  3*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x +  4*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x +  5*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x +  6*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x +  7*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x +  8*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x +  9*L, t, x0, theta_0, a, v, beta)
        + bright_soliton(x + 10*L, t, x0, theta_0, a, v, beta)
        )
    
    return u



def dark_soliton(x, t, x0, theta_0, u0, v, beta, phi):
    
    assert(beta > 0)
    assert(u0 > 0)
    
    a = sqrt(u0 * beta) * cos(phi)
    
    return (1 / sqrt(beta)) * (a * tanh(a*(x-v*t-x0)) + 1j*v) * exp(1j*(-(a**2 + v**2)*t + theta_0 ) )







