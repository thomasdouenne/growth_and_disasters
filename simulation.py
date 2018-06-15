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
    first = (1-e)*(1+d)*(1-w**(1-g))/(1-g)
    second = (s**(s/(1-s)) + e*s**(1/(1-s))-e*s**(s/(1-s)))
    second_bis = l**(s/(1-s))/(1-s)
    second_ter = ((1-w**(1-g))/(a**s*(1-g)))**(1/(1-s))
    third = (1+d)*(1-w)
    fourth = (1-w)/(1-s)
    fourth_bis = ((1-w**(1-g))*l*s/(a*(1-g)))**(s/(1-s))
    
    return first - second*second_bis*second_ter - third + fourth*fourth_bis


# Fix parameters' value :

d = 1
w = 0.8
l = 0.02
s = 1/3
a = 0.05

# Create the grid

e = arange(0.05,5,0.025) + 0.01
g = arange(0.05,10,0.05) + 0.01
E,G = meshgrid(e, g)

# Apply function and build graphs
gr_lambda = lg_growth_by_lambda(E, G) # evaluation of the function on the grid

# Check the function is correct
print(lg_growth_by_lambda(0.5,2))


fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(E, G, gr_lambda, norm=norm, rstride=1, cstride=1, 
                      cmap=cm.RdBu,linewidth=0, antialiased=False)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

fig, ax = plt.subplots()
heatmap = ax.imshow(gr_lambda, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none')
plt.xlabel('ε = IES')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap)
cbar.set_label('dg*/dλ')
plt.show()
