import numpy as np
import pandas as pd
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


def effect_disasters_expected_growth(e,g,l,w):
    
    return expected_growth(e,g,l,w) - expected_growth(e,g,0,w)


# Fix parameters' value :
d = 1
l = 0.0307
w = 0.948
s = 0.02
u = 0.25
a = 0.069
r = 0.025
e = .5
g = 4

w_1 = 0.948
w_2 = 0.85
w_3 = 0.60
l_1 = 0.0307
l_2 = 0.01064
l_3 = 0.003991

print 'Share of production consumed', psi(e,g,l,w)/a*100, '%'
print 'Risk mitigation share:', theta(e,g,l,w)*100, '%'
print 'Disaster ex post probability:', l*(1+d-(theta(e,g,l,w))**u)
print 'Trend growth:', trend_growth(e,g,l,w)*100, '%'
print 'Expected growth:', expected_growth(e,g,l,w)*100, '%'

print effect_disasters_expected_growth(.2,8,0.001,0.5)

### Heatmaps:

# Create the grid
e_inverse = 1/(arange(4,1,-0.01) - 0.001)
e_normal = arange(1,3,0.01) + 0.001
e = np.concatenate([e_inverse,e_normal])
g = arange(6,1,-0.01) + 0.005
E,G = meshgrid(e, g)

# Build heatmap
values_heatmap = effect_disasters_expected_growth(E, G, 0.001, 0.5) # evaluation of the function on the grid


# Prepare axes
axe_g = ['6', '5', '4', '3', '2', '1']
axe_e = ['1/4', '1/3', '1/2', '1', '2', '3'] # These two graduation are incorrect

fig, ax = plt.subplots()
heatmap_expected_growth = ax.imshow(
    values_heatmap, norm=norm, cmap=plt.cm.seismic, extent = [0,10,10,0], interpolation='none'
    )
ax.set_xticklabels(axe_e)
ax.set_yticklabels(axe_g)
plt.xlabel('IES')
plt.ylabel('RRA')
cbar = plt.colorbar(heatmap_expected_growth)
cbar.set_label('percentage points')
plt.show()

""" Trash:

def mrt_sigma_gdp(e,g,l,w):

    return (g*s)/psi(e,g,l,w)

    
# Alternative calibrations for expected growth
for epsilon in [0.25, 0.5, 1+1e-09, 1.5, 2]:
    print 'with e=', epsilon, ', w=', w_1, 'and l=', l_1, ':', expected_growth(epsilon,g,0.0307,0.948)*100, '%'
    print 'with e=', epsilon, ', w=', w_2, 'and l=', l_2, ':', expected_growth(epsilon,g,0.017,0.71)*100, '%'
    print 'with e=', epsilon, ', w=', w_3, 'and l=', l_3, ':', expected_growth(epsilon,g,0.001,0.5)*100, '%'


print 'Share of production consumed', psi(e,g,l,w)/a*100, '%'
print 'Risk mitigation share:', tau(e,g,l,w)*100, '%'
print 'Disaster ex post probability:', l*(1+d-(tau(e,g,l,w))**u)
print 'Reduction in risk of env. disaster:', tau(e,g,l,w)**u
print 'Trend growth:', trend_growth(e,g,l,w)*100, '%'
print 'Expected growth:', expected_growth(e,g,l,w)*100, '%'

print mrt_lambda_gdp(e,g,l,w)
print mrt_omega_gdp(e,g,l,w)
print mrt_sigma_gdp(e,g,l,w)

"""
