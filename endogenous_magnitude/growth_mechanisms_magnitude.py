import numpy as np
from numpy import arange
from pylab import meshgrid

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


#☼ Comprendre l'effet non-monotone observé sur les calibrations

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


def damage_effect_lambda(e,g):
    
    return -(1-omega(e,g))


def consumption_effect_lambda(e,g):

    return (1-e)*(1-(omega(e,g))**(1-g))/(1-g)


def consumption_effect_deltae(e,g):

    return (1-e)*(1-dk-omega(e,g))/(de**2*s)


def conservation_effect_lambda(e,g):

    return omega(e,g)/g


def conservation_effect_deltae(e,g):

    return omega(e,g)*l/(g*de)


def substitution_from_capital_effect_lambda(e,g):
 
    return omega(e,g)**(1-g)/g


def substitution_from_capital_effect_deltae(e,g):

    return -(1-dk+((1-g)/g)*omega(e,g))/(de**2*s)


# Fix parameters' value :
de = 1
dk = 0.01
l = 0.03
a = 0.05
r = 0.015
s = 20
f = 4 # so spending 20% of GDP in mitigation pull the risk to 0


# Create the grid
e = arange(0.05,5,0.025) + 0.01
g = arange(0.05,10,0.05) + 0.01
E,G = meshgrid(e, g)


# Apply functions and build graphs
damage_lambda = damage_effect_lambda(0.5,g)

print("Damage effect on growth")

plt.title("Damage effect - lambda")
plt.plot(g,damage_lambda)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()


conservation_lambda = conservation_effect_lambda(0.5,g)
conservation_deltae = conservation_effect_deltae(0.5,g)

print("Conservation effect on growth")

plt.title("Conservation effect - lambda")
plt.plot(g,conservation_lambda)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Conservation effect - deltae")
plt.plot(g,conservation_deltae)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()


consumption_lambda = consumption_effect_lambda(0.5,g)
consumption_deltae = consumption_effect_deltae(0.5,g)

print("Consumption effect on growth")

plt.title("Consumption effect - lambda")
plt.plot(g,consumption_lambda)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Consumption effect - deltae")
plt.plot(g,consumption_deltae)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()


substitution_lambda = substitution_from_capital_effect_lambda(0.5,g)
substitution_deltae = substitution_from_capital_effect_deltae(0.5,g)

print("Substitution effect on growth")

plt.title("Substitution effect - lambda")
plt.plot(g,substitution_lambda)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Substitution effect - deltae")
plt.plot(g,substitution_deltae)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()


aggregate_lambda = (
        damage_lambda + conservation_lambda + consumption_lambda + substitution_lambda
        )
dl_growth = dl_lr_growth(0.5,g,l,de)
aggregate_deltae = (
        conservation_deltae + consumption_deltae + substitution_deltae
        )
dw_growth = dde_lr_growth(0.5,g,l,de)

print("Aggregate effect on growth")

plt.title("Aggregate effect - lambda")
plt.plot(g,aggregate_lambda)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Check aggregate effect numerically - lambda")
plt.plot(g,dl_growth)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Aggregate effect - deltae")
plt.plot(g,aggregate_deltae)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Check aggregate effect numerically) - omega")
plt.plot(g,dw_growth)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()


# Build and draw heatmaps
gr_lambda = dl_lr_growth(E, G, l, de) # evaluation of the function on the grid
gr_de = dde_lr_growth(E, G, l, de) # evaluation of the function on the grid

fig, ax = plt.subplots()
heatmap_lambda = ax.imshow(
    gr_lambda, norm=norm, cmap=plt.cm.seismic, extent = [0,10,10,0], interpolation='none'
    )
plt.xlabel('1/ε = aversion to fluc.')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_lambda)
cbar.set_label('dg*/dλ')
plt.show()

fig, ax = plt.subplots()
heatmap_de = ax.imshow(
    gr_de, norm=norm, cmap=plt.cm.seismic, extent = [0,10,10,0], interpolation='none'
    )
plt.xlabel('1/ε = aversion to fluc.')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_de)
cbar.set_label('dg*/dw')
plt.show()
