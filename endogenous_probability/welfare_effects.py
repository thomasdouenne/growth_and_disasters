import numpy as np
from numpy import arange
from pylab import meshgrid

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


# Normalize the color of graphs to 0 :
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


norm = MidpointNormalize(midpoint=0)


# Define functions:
def theta(e,g,l,w):

    return (((1-w**(1-g))*l*u)/(a*(1-g)))**(1/(1-u))


def psi(e,g,l,w):

    return e*r+(1-e)*(
        (1-theta(e,g,l,w))*a - (g*(s**2)/2) - l*(1+d-(theta(e,g,l,w))**u) * (1-w**(1-g))/(1-g)
        )


def lg_growth(e,g,l,w):
    
    return (1-theta(e,g,l,w))*a - psi(e,g,l,w) - l*(1+d-(theta(e,g,l,w)**u))*(1-w)


def trend_growth(e,g,l,w):
    
    return (1-theta(e,g,l,w))*a - psi(e,g,l,w)


def expected_growth(e,g,l,w):
    
    return trend_growth(e,g,l,w) - l*(1+d-(theta(e,g,l,w))**u) * (1-w)


def mrt_sigma_gdp(e,g,l,w):

    return (g*s)/psi(e,g,l,w)

    
def mrt_lambda_gdp(e,g,l,w):
    group_1 = (u**(u/(1-u))-u**(1/(1-u)))/(1-u)
    group_2 = ((1-w**(1-g))/((a**u)*(1-g)))**(1/(1-u))
    group_3 = (1+d)*(1-w**(1-g))/(1-g)
    
    return -(1/psi(e,g,l,w))*(
            l**(u/(1-u)) * group_1 * group_2 - group_3
            )


def mrt_omega_gdp(e,g,l,w):
    group_1 = (u**(u/(1-u))-u**(1/(1-u)))/(1-u)
    group_2 = ((1-w**(1-g))/(a*(1-g)))**(u/(1-u))
    
    return -(w**(-g))/psi(e,g,l,w)*(
            l*(1+d) - l**(1/(1-u)) * group_1 * group_2
            )


def lucas_measure(e,g,l,w):
    psi_optimal = psi(e,g,l,w)
    psi_bau = e*r+(1-e)*(
        (1-0)*a - (g*(s**2)/2) - l*(1+d-(0)**u) * (1-w**(1-g))/(1-g)
        )
    
    return (psi_optimal/psi_bau)**(1/(1-e)) - 1


def rho_to_fit_growth(e,g,l,w):
    
    group_1 = a*(1-theta(e,g,l,w))
    group_2 = l*(1+d-(theta(e,g,l,w))**u) * (1-w)
    group_3 = g*(s**2)/2
    group_4 = l*(1+d-(theta(e,g,l,w))**u) * (1-w**(1-g)) / (1-g)

    return 1/e * (group_1 - group_2 - (1-e)*(group_1 - group_3 - group_4) - 0.02)

# Fix parameters' value :
d = 1
w = 0.6
l = 0.003991
s = 0.02
u = 0.25
a = 0.069
e = 1.5
g = 1+1e-12

r = rho_to_fit_growth(e,g,l,w)
print r

print 'Expected growth with right rho:', expected_growth(e,g,l,w)*100, '%'
print 'Share of production consumed', psi(e,g,l,w)/a*100, '%'

print lucas_measure(e,g,l,w)

#### Calibration rho



print 'Share of production consumed', psi(e,g,l,w)/a*100, '%'
print 'Risk mitigation share:', theta(e,g,l,w)*100, '%'
print 'Disaster ex post probability:', l*(1+d-(theta(e,g,l,w))**u)
print 'Trend growth:', trend_growth(e,g,l,w)*100, '%'
print 'Expected growth:', expected_growth(e,g,l,w)*100, '%'

print mrt_lambda_gdp(e,g,l,w)
print mrt_omega_gdp(e,g,l,w)
print mrt_sigma_gdp(e,g,l,w)
