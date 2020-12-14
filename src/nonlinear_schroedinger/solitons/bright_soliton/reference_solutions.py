from numpy import exp
from numpy import cosh
from numpy import sqrt

def bright_soliton(x, t, a, v, x0, theta_0, beta):
    
    assert(beta < 0)
    
    return (a / sqrt(-beta)) * ( 1.0 / cosh( a * (x - v*t - x0) ) ) * exp( 1j * ( v*x - 0.5*(v**2-a**2)*t + theta_0 ) )
    
    