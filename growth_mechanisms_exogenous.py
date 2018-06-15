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


#############Check sign errors


# Define functions:
def lg_growth_by_lambda(e,g):
    
    return (1-e)*(1-w**(1-g))/(1-g) - (1-w)


def lg_growth_by_omega(e,g):
    
    return -l*(1-(1-e)*w**(-g))


def damage_effect_lambda(e,g): # effect from higher frequency of disasters

    return -(1-w) + g*0


def damage_effect_omega(e,g): # effect from larger damages of disasters

    return -l + g*0


def consumption_effect_lambda(e,g):

    return (1-e)/(1-g)*(1-w**(1-g))


def consumption_effect_omega(e,g):

    return (1-e)*l*w**(-g)



# Fix parameters' value :

d = 1
w = 0.95
l = 0.02

# Create the grid

e = arange(0.05,5,0.025) + 0.01
g = arange(0.05,10,0.05) + 0.01
E,G = meshgrid(e, g)

# Apply function and build graphs
damage_lambda = damage_effect_lambda(0.5,g)
damage_omega = damage_effect_omega(0.5,g)

plt.title("Damage effect - lambda")
plt.plot(g,damage_lambda)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Damage effect - omega")
plt.plot(g,damage_omega)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

consumption_lambda = consumption_effect_lambda(0.5,g)
consumption_omega = consumption_effect_omega(0.5,g)

plt.title("Consumption effect - lambda")
plt.plot(g,consumption_lambda)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Consumption effect - omega")
plt.plot(g,consumption_omega)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()


aggregate_lambda = (
        damage_lambda + consumption_lambda 
        )
aggregate_omega = (
        damage_omega + consumption_omega
        )

plt.title("Aggregate effect - lambda")
plt.plot(g,aggregate_lambda)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Aggregate effect - omega")
plt.plot(g,aggregate_omega)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()



lg_growth_l = lg_growth_by_lambda(0.5,g)
lg_growth_o = lg_growth_by_omega(0.5,g)

plt.title("Check aggregate effect (growth) - lambda")
plt.plot(g,lg_growth_l)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Check aggregate effect (growth) - omega")
plt.plot(g,lg_growth_o)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()
