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
        (1-theta(e,g,l,w))*a - (g*s/2) - l*(1+d-(theta(e,g,l,w))**u) * (1-w**(1-g))/(1-g)
        )


def lg_growth(e,g,l,w):
    
    return (1-theta(e,g,l,w))*a - psi(e,g,l,w) - l*(1+d-(theta(e,g,l,w)**u))*(1-w)


def dl_lg_growth(e,g,l,w):
    h = 1e-6

    return (lg_growth(e,g,l+h,w) - lg_growth(e,g,l-h,w))/(2*h)


def dw_lg_growth(e,g,l,w):
    h = 1e-6

    return -(lg_growth(e,g,l,w+h) - lg_growth(e,g,l,w-h))/(2*h)


def damage_effect_lambda(e,g): # effect from higher frequency of disasters
    one = (1+d)*(1-w)
    two = (1-w)*(((1-w**(1-g))*l*u)/(a*(1-g)))**(u/(1-u))

    return - one + two


def damage_effect_omega(e,g): # effect from larger damages of disasters
    one = l*(1+d)
    two = l*(((1-w**(1-g))*l*u)/(a*(1-g)))**(u/(1-u))

    return - one + two


def consumption_effect_lambda(e,g):
    one = (1+d)*(1-w**(1-g))/(1-g)
    two = l**(1/(1-u))*(u**(u/(1-u))-u**(1/(1-u)))/(1-u)
    three = ((1-w**(1-g))/(a**u*(1-g)))**(1/(1-u))
    four = 1/l

    return (1-e)*(one - two*three*four)


def consumption_effect_omega(e,g):
    one = (1+d)*l
    two = l**(1/(1-u))*(u**(u/(1-u))-u**(1/(1-u)))/(1-u)
    three = ((1-w**(1-g))/(a**u*(1-g)))**(1/(1-u))
    four = (1-g)/(1-w**(1-g))

    return (1-e)*w**(-g)*(one - two*three*four)


def conservation_effect_lambda(e,g):
    one = u*(1-w)/(1-u)
    two = ((1-w**(1-g))*l*u/(a*(1-g)))**(u/(1-u))

    return one*two


def conservation_effect_omega(e,g):
    one = (1-w)/(1-u)
    two = ((1-w**(1-g))*l*u/(a*(1-g)))**(u/(1-u) - 1)
    three = (l**2*u**2)/a*w**(-g)

    return one*two*three


def substitution_from_capital_effect_lambda(e,g):
    one = a*l**(u/(1-u))/(1-u)
    two = (((1-w**(1-g))*u)/(a*(1-g)))**(1/(1-u))
 
    return - one*two


def substitution_from_capital_effect_omega(e,g):
    one = a*((l*u)/a)**(1/(1-u))
    two = w**(-g)/(1-u)
    three = ((1-w**(1-g))/(1-g))**(u/(1-u))

    return -one*two*three


# Fix parameters' value :
d = 3
w = 0.79
l = 0.01
s = 0.02
u = 1/4
a = 0.059
r = 0.029

# Create the grid
e_inverse = 1/(arange(4,1,-0.01) - 0.001)
e_normal = arange(1,3,0.01) + 0.001
e = np.concatenate([e_inverse,e_normal])
g = arange(6,1,-0.01) + 0.005
E,G = meshgrid(e, g)


# Apply functions and build graphs
damage_lambda = damage_effect_lambda(0.5,g)
damage_omega = damage_effect_omega(0.5,g)

print("Damage effect on growth")

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

print("Conservation effect on growth")

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

print("Consumption effect on growth")

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

print("Substitution effect on growth")

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
dl_growth = dl_lg_growth(0.5,g,l,w)

aggregate_omega = (
        damage_omega + conservation_omega + consumption_omega + substitution_omega
        )
dw_growth = dw_lg_growth(0.5,g,l,w) # Et ici aussi dans une moindre mesure

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

plt.title("Aggregate effect - omega")
plt.plot(g,aggregate_omega)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()

plt.title("Check aggregate effect numerically - omega")
plt.plot(g,dw_growth)
plt.xlabel('Risk aversion coefficient')
plt.ylabel('Effect on growth')
plt.show()


# Build heatmap
effect_lambda_num = dl_lg_growth(E, G, l,w) # evaluation of the function on the grid
effect_omega_num = dw_lg_growth(E, G, l,w) # evaluation of the function on the grid

# Prepare axes
axe_g = ['6', '5', '4', '3', '2', '1']
axe_e = ['1/4', '1/3', '1/2', '1', '2', '3'] # These two graduation are incorrect

fig, ax = plt.subplots()
heatmap_lambda_num = ax.imshow(
    effect_lambda_num, norm=norm, cmap=plt.cm.seismic, extent = [0,10,10,0], interpolation='none'
    )
ax.set_xticklabels(axe_e)
ax.set_yticklabels(axe_g)
plt.xlabel('1/ε = aversion to fluc.')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_lambda_num)
cbar.set_label('dg*/dλ')
plt.show()

fig, ax = plt.subplots()
heatmap_omega_num = ax.imshow(
    effect_omega_num, norm=norm, cmap=plt.cm.seismic, extent = [0,10,10,0], interpolation='none'
    )
ax.set_xticklabels(axe_e)
ax.set_yticklabels(axe_g)
plt.xlabel('1/ε = aversion to fluc.')
plt.ylabel('γ = risk aversion coef.')
cbar = plt.colorbar(heatmap_omega_num)
cbar.set_label('dg*/dλ')
plt.show()

print(theta(2,4,l,w)*100)
print(psi(2,3.3,l,w)/a)
print(l*(1+d-(theta(2,3.3,l,w))**u))
print(consumption_effect_lambda(2,4))
