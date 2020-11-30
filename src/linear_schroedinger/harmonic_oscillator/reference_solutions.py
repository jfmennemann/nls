from numpy import exp
from numpy import pi

def coherent_state(x, t, x0, omega):
    
    return (omega/pi)**0.25 * exp( -(omega/2) * ( x**2 - 2*x*x0*exp(-1j*omega*t) + (x0**2/2.0) * exp(-2*1j*omega*t) + x0**2/2.0 ) - 1j * omega * t / 2.0 )

