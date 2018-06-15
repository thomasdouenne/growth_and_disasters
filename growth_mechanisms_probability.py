import numpy as np
from numpy import arange
from pylab import meshgrid

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize



# To do : revoir la gradation des axes verticaux, ainsi que les couleurs (rouge pour négatif)


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
    second = l*(1+d)
    third = (s**(s/(1-s)) + e*s**(1/(1-s))-e*s**(s/(1-s)) - (1-w)*s**(1/(1-s))*(1-g)/(1-w**(1-g)))
    third_bis = w**(-g)/(1-s)
    third_ter = s**(s/(1-s))
    fourth = (l/a**s)**(1/(1-s))
    fifth = ((1-w**(1-g))/((1-g)))**(s/(1-s))
    
    return first - second - (third*third_bis - third_ter)*fourth*fifth


def damage_effect_lambda(e,g): # effect from higher frequency of disasters
    one = (1+d)*(1-w)
    two = (1-w)*(((1-w**(1-g))*l*s)/(a*(1-g)))**(s/(1-s))

    return - one + two


def damage_effect_omega(e,g): # effect from larger damages of disasters
    one = l*(1+d)
    two = l*(((1-w**(1-g))*l*s)/(a*(1-g)))**(s/(1-s))

    return - one + two


def consumption_effect_lambda(e,g):
    one = (1+d)*(1-w**(1-g))/(1-g)
    two = l**(1/(1-s))*(s**(s/(1-s))-s**(1/(1-s)))/(1-s)
    three = ((1-w**(1-g))/(a**s*(1-g)))**(1/(1-s))
    four = 1/l

    return (1-e)*(one - two*three*four)


def consumption_effect_omega(e,g):
    one = (1+d)*l
    two = l**(1/(1-s))*(s**(s/(1-s))-s**(1/(1-s)))/(1-s)
    three = ((1-w**(1-g))/(a**s*(1-g)))**(1/(1-s))
    four = (1-g)/(1-w**(1-g))

    return (1-e)*w**(-g)*(one - two*three*four)


def conservation_effect_lambda(e,g):
    one = s*(1-w)/(1-s)
    two = ((1-w**(1-g))*l*s/(a*(1-g)))**(s/(1-s))

    return one*two


def conservation_effect_omega(e,g):
    one = (1-w)/(1-s)
    two = ((1-w**(1-g))*l*s/(a*(1-g)))**(s/(1-s) - 1)
    three = (l**2*s**2)/a*w**(-g)

    return one*two*three


def substitution_from_capital_effect_lambda(e,g):
    one = a*l**(s/(1-s))/(1-s)
    two = (((1-w**(1-g))*s)/(a*(1-g)))**(1/(1-s))
 
    return - one*two


def substitution_from_capital_effect_omega(e,g):
    one = a*((l*s)/a)**(1/(1-s))
    two = w**(-g)/(1-s)
    three = ((1-w**(1-g))/(1-g))**(s/(1-s))

    return -one*two*three


# Fix parameters' value :

d = 1
w = 0.95
l = 0.02
s = 0.7
a = 0.05
r = 0.15

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

conservation_lambda = conservation_effect_lambda(0.5,g)
conservation_omega = conservation_effect_omega(0.5,g)

plt.title("Conservation effect - lambda")
plt.plot(g,conservation_lambda)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Conservation effect - omega")
plt.plot(g,conservation_omega)
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


substitution_lambda = substitution_from_capital_effect_lambda(0.5,g)
substitution_omega = substitution_from_capital_effect_omega(0.5,g)

plt.title("Substitution effect - lambda")
plt.plot(g,substitution_lambda)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Substitution effect - omega")
plt.plot(g,substitution_omega)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()


aggregate_lambda = (
        damage_lambda + conservation_lambda + consumption_lambda + substitution_lambda
        )
lg_growth_l = lg_growth_by_lambda(0.5,g)
dl_growth = dl_lg_growth(0.5,g,l,w)


plt.title("Aggregate effect - lambda")
plt.plot(g,aggregate_lambda)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Check aggregate effect (growth) - lambda")
plt.plot(g,lg_growth_l)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Check aggregate effect numerically - lambda")
plt.plot(g,dl_growth)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()


aggregate_omega = (
        damage_omega + conservation_omega + consumption_omega + substitution_omega
        )
lg_growth_o = lg_growth_by_omega(0.5,g) # Problème ici !
dw_growth = dw_lg_growth(0.5,g,l,w) # Et ici aussi dans une moindre mesure


plt.title("Aggregate effect - omega")
plt.plot(g,aggregate_omega)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Check aggregate effect (growth) - omega")
plt.plot(g,lg_growth_o)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Check aggregate effect numerically) - omega")
plt.plot(g,dw_growth)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()


# Try to understand why it does not match with aggregate effect

# Vérifier numériquement la calibration de tous les effets...



"""






# Draw heatmap - lambda - conservation
fig, ax = plt.subplots()
heatmap_conservation_lambda = ax.imshow(
    conservation_lambda, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('ε = IES')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_conservation_lambda)
cbar.set_label('dg*/dw')
plt.show()

# Draw heatmap - omega - conservation
fig, ax = plt.subplots()
heatmap_conservation_omega = ax.imshow(
    conservation_omega, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('ε = IES')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_conservation_omega)
cbar.set_label('dg*/dw')
plt.show()


print("Consumption:")

# Draw heatmap - lambda - consumption
fig, ax = plt.subplots()
heatmap_consumption_lambda = ax.imshow(
    consumption_lambda, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('ε = IES')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_consumption_lambda)
cbar.set_label('dg*/dw')
plt.show()

# Draw heatmap - omega - consumption
fig, ax = plt.subplots()
heatmap_consumption_omega = ax.imshow(
    consumption_omega, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('ε = IES')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_consumption_omega)
cbar.set_label('dg*/dw')
plt.show()


print("Damage:")

# Draw heatmap - lambda - damage
fig, ax = plt.subplots()
heatmap_damage_lambda = ax.imshow(
    damage_lambda, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('ε = IES')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_damage_lambda)
cbar.set_label('dg*/dw')
plt.show()

# Draw heatmap - omega - damage
fig, ax = plt.subplots()
heatmap_damage_omega = ax.imshow(
    damage_omega, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('ε = IES')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_damage_omega)
cbar.set_label('dg*/dw')
plt.show()

print("Substitution :")

# Draw heatmap - lambda - substitution
fig, ax = plt.subplots()
heatmap_substitution_lambda = ax.imshow(
    substitution_lambda, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('ε = IES')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_substitution_lambda)
cbar.set_label('dg*/dw')
plt.show()

# Draw heatmap - omega - substitution
fig, ax = plt.subplots()
heatmap_substitution_omega = ax.imshow(
    damage_omega, norm=norm, cmap=plt.cm.seismic, extent = [0,5,10,0], interpolation='none'
    )
plt.xlabel('ε = IES')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_substitution_omega)
cbar.set_label('dg*/dw')
plt.show()

print(consumption_effect_lambda(0.5,20))
print(consumption_effect_omega(0.5,20))



"""

