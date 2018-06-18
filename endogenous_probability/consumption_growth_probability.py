import pandas as pd
import numpy as np
from numpy import arange
from pylab import meshgrid
import random

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


def trend_growth(e,g,l,w):
    
    return ((1-theta(e,g,l,w))*a - psi(e,g,l,w))


def stochastic_consumption(e,g,l,w,t):

    return psi(e,g,l,w) * k * (1 + trend_growth(e,g,l,w))**t * (1-w)*(random.randrange(1, 51) == 50)


# Fix parameters' value :
d = 1
w = 0.95
l = 0.02
s = 0.7
a = 0.05
r = 0.015
k = 1 # normalized initial capital stock


#Create vector of consumption
time = arange(0,250,1)
trend_g = trend_growth(0.5,4,l,w)

df = pd.DataFrame(index = time, columns = ['time', 'shock', 'consumption'])
df.time = time
liste = []
for t in time:
    liste.append(1*(random.randrange(1, 51) == 50))
df['shock'] = liste
df['consumption'] = 1

for i in range(0, len(df)-1):
    df.ix[i+1, 'consumption'] = (
            (1+trend_g) * df.ix[i, 'consumption']
            - (1-0.8)*df.ix[i+1, 'shock'] * df.ix[i, 'consumption']
            )


# Create graphs
plt.title("Stochastic consumption growth")
plt.plot(df['time'],df['consumption'])
plt.xlabel('Time')
plt.ylabel('Consumption')
plt.show()

