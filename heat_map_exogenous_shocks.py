# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 15:15:27 2018

@author: t.douenne
"""

import numpy as np
from numpy import exp,arange
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
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
    
    return (1-e)*(1-w**(1-g))/(1-g) - (1-w)


def lg_growth_by_omega(e,g):
    
    return l*(1-(1-e)*w**(-g))



# Fix parameters' value :

d = 1
w = 0.95
l = 0.025

# Create the grid

e = 1/(arange(0.05,5,0.025) + 0.01)
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
plt.xlabel('1/ε = aversion to fluc.')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_lambda)
cbar.set_label('dg*/dλ')
plt.show()

# Draw heatmap - omega
fig, ax = plt.subplots()
heatmap_omega = ax.imshow(
    gr_omega, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('1/ε = aversion to fluc.')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_omega)
cbar.set_label('dg*/dw')
plt.show()
