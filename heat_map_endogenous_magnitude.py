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
def omega(e,g):
    return (l*s*de)**(1/g)


def theta(e,g):
    return f/s - (1-dk-omega(e,g))/(de*s*a)


def psi(e,g):
    return r*e + (1-e)*(1-theta(e,g))*a - (1-e)*l*(1-(omega(e,g))**(1-g))/(1-g)


def lr_growth(e,g,l,de):
    return (1-theta(e,g))*a - psi(e,g) - l*(1-omega(e,g))


def dl_lr_growth(e,g,l,de):
    h = 1e-6

    return (lr_growth(e,g,l+h,de) - lr_growth(e,g,l-h,de))/(2*h)


def dde_lr_growth(e,g,l,de):
    h = 1e-6

    return (lr_growth(e,g,l,de+h) - lr_growth(e,g,l,de-h))/(2*h)


# Fix parameters' value :
de = 1
dk = 0.01
l = 0.04
a = 0.05
r = 0.015
s = 20
f = 5 # so spending 25% of GDP in mitigation pull the risk to 0

# Create the grid
e = 1/(arange(0.05,10,0.05) + 0.01)
g = arange(0.05,10,0.05) + 0.01
E,G = meshgrid(e, g)

# Apply function and build graphs
gr_lambda = dl_lr_growth(E, G, l, de) # evaluation of the function on the grid
gr_de = dde_lr_growth(E, G, l, de) # evaluation of the function on the grid


# Draw heatmap - lambda
fig, ax = plt.subplots()
heatmap_lambda = ax.imshow(
    gr_lambda, norm=norm, cmap=plt.cm.seismic, extent = [0,10,10,0], interpolation='none'
    )
plt.xlabel('1/ε = aversion to fluc.')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_lambda)
cbar.set_label('dg*/dλ')
plt.show()

# Draw heatmap - omega
fig, ax = plt.subplots()
heatmap_de = ax.imshow(
    gr_de, norm=norm, cmap=plt.cm.seismic, extent = [0,10,10,0], interpolation='none'
    )
plt.xlabel('1/ε = aversion to fluc.')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_de)
cbar.set_label('dg*/dw')
plt.show()

print(omega(0.5,2)) # Gros problème : omega beaucoup trop faible !
print(theta(0.5,2))
