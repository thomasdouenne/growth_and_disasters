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

def lg_growth_by_lambda(e,g):
    first = (1-e)*(1+d)*(1-w**(1-g))/(1-g)
    second = (s**(s/(1-s)) + e*s**(1/(1-s))-e*s**(s/(1-s)))
    second_bis = l**(s/(1-s))/(1-s)
    second_ter = ((1-w**(1-g))/(a**s*(1-g)))**(1/(1-s))
    third = (1+d)*(1-w)
    fourth = (1-w)/(1-s)
    fourth_bis = ((1-w**(1-g))*l*s/(a*(1-g)))**(s/(1-s))
    
    return first - second*second_bis*second_ter - third + fourth*fourth_bis


def lg_growth_by_omega(e,g): # Check the formula, I have weird results
    first = (1-e)*l*(1+d)*w**(-g)
    second = (s**(s/(1-s)) + e*s**(1/(1-s))-e*s**(s/(1-s)))
    second_bis = l**(s/(1-s))/(1-s)
    second_ter = ((1-w**(1-g))/(a*(1-g)))**(s/(1-s))*w**(-g)
    third = l*(1+d)
    fourth = l*((1-w**(1-g))/(a*(1-g)))**(s/(1-s))
    fifth = (1-w)/(1-s)
    fifth_bis = ((1-w**(1-g))*l*s/(a*(1-g)))**(s/(1-s) - 1)
    fifth_ter = (l**2*s**2)/a*w**(-g)
    
    return -first + second*second_bis*second_ter + third - fourth - fifth*fifth_bis*fifth_ter


# Fix parameters' value :

d = 1
w = 0.95
l = 0.02
s = 0.75
a = 0.05

# Create the grid

e = arange(0.05,5,0.025) + 0.01
g = arange(0.05,10,0.05) + 0.01
E,G = meshgrid(e, g)

# Apply function and build graphs
gr_lambda = lg_growth_by_lambda(E, G) # evaluation of the function on the grid
gr_omega = -lg_growth_by_omega(E, G) # evaluation of the function on the grid


# Draw heatmap - lambda
fig, ax = plt.subplots()
heatmap_lambda = ax.imshow(
    gr_lambda, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('ε = IES')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_lambda)
cbar.set_label('dg*/dλ')
plt.show()

# Draw heatmap - omega
fig, ax = plt.subplots()
heatmap_omega = ax.imshow(
    gr_omega, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('ε = IES')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_omega)
cbar.set_label('dg*/dw')
plt.show()

print(lg_growth_by_omega(0.5,2))



