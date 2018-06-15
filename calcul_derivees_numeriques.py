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
def theta(e,g,l,w):

    return (((1-w**(1-g))*l*s)/(a*(1-g)))**(1/(1-s))


def psi(e,g,l,w):

    return e*r+(1-e)*((1-theta(e,g,l,w))*a - l*(1+d-(theta(e,g,l,w))**s) *(1-w**(1-g))/(1-g))


def lg_growth(e,g,l,w):
    
    return (1-theta(e,g,l,w)*a - psi(e,g,l,w) - l*(1+d-(theta(e,g,l,w)**s))*(1-w))


def dl_lg_growth(e,g,l,w):
    h = 1e-6

    return (lg_growth(e,g,l+h,w) - lg_growth(e,g,l-h,w))/(2*h)


def dw_lg_growth(e,g,l,w):
    h = 1e-6

    return -(lg_growth(e,g,l,w+h) - lg_growth(e,g,l,w-h))/(2*h)


d = 1
w = 0.9
l = 0.025
s = 0.7
a = 0.05
r = 0.015

e = 1/(arange(0.05,5,0.025) + 0.01)
g = arange(0.05,10,0.05) + 0.01
E,G = meshgrid(e, g)

# Apply function and build graphs
effect_lambda_num = dl_lg_growth(E, G, l,w) # evaluation of the function on the grid

effect_omega_num = dw_lg_growth(E, G, l,w) # evaluation of the function on the grid


# Draw heatmap - lambda
fig, ax = plt.subplots()
heatmap_lambda_num = ax.imshow(
    effect_lambda_num, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('1/ε = aversion to fluc.')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_lambda_num)
cbar.set_label('dg*/dλ')
plt.show()

# Draw heatmap - omega
fig, ax = plt.subplots()
heatmap_omega_num = ax.imshow(
    effect_omega_num, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('1/ε = aversion to fluc.')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_omega_num)
cbar.set_label('dg*/dλ')
plt.show()
