from math import sqrt
from scipy.stats import norm
import pandas as pd
import numpy as np
from numpy import arange
from pylab import plot, show, grid, xlabel, ylabel
import random

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


"""
brownian() implements one dimensional Brownian motion (i.e. the Wiener process).
"""
def brownian(x0, n, dt, delta, out=None):
    x0 = np.asarray(x0)

    # For each element of x0, generate a sample of n numbers from a
    # normal distribution.
    r = norm.rvs(size=x0.shape + (n,), scale=delta*sqrt(dt))

    # If `out` was not given, create an output array.
    if out is None:
        out = np.empty(r.shape)

    # This computes the Brownian motion by forming the cumulative sum of
    # the random samples. 
    np.cumsum(r, axis=-1, out=out)

    # Add the initial condition.
    out += np.expand_dims(x0, axis=-1)

    return out


# Define functions:
def psi(e,g,l,w):

    return e*r + (1-e)*(a - (g*s/2) - l*(1-w**(1-g))/(1-g))


def trend_growth(e,g,l,w):
    
    return e*(a-r) + (1-e)*((g*s/2) + l*(1-w**(1-g))/(1-g))


# Fix parameters' value :
d = 3
w = 0.79
l = 0.01
s = 0.02
a = 0.059
r = 0.029
T = 200 # Total time.
N = 2400 # Number of steps.
dt = T/N # Time step size
x = np.empty((1,N+1)) # Create an empty array to store the realizations.
x[:, 0] = 1 # Initial values of x.


#Create vector of consumption
time = arange(0,N,1)
trend_g = trend_growth(2,3.3,l,w)

df = pd.DataFrame(index = time, columns = ['time', 'shock', 'consumption'])
df.time = time
liste = []
for t in time:
    liste.append(1* (random.uniform(0, 1) < l*(1+d)*dt))
df['shock'] = liste
df['consumption'] = 1
df['capital_stock'] = 1/psi(2,3.3,l,w) # normalized initial capital stock: normalized consumption over psi
df['brownian'] = brownian(x[:,0], N, dt, s, out=x[:,1:]).transpose()

for i in range(0, len(df)-1):
    df.ix[i+1, 'consumption'] = (
            (1+trend_g)**(dt) * df.ix[i, 'consumption']
            + (df.ix[i+1, 'brownian'] - df.ix[i, 'brownian']) * s * df.ix[i, 'capital_stock']
            - (1-w)*df.ix[i+1, 'shock'] * df.ix[i, 'consumption']
            )
    df.ix[i+1, 'capital_stock'] = (
            df.ix[i+1, 'consumption'] / psi(2,3.3,l,w)
            )


# Create graphs
plt.title("Stochastic consumption growth")
plt.plot(df['time'],df['consumption'])
plt.xlabel('Time')
plt.ylabel('Consumption')
plt.show()

print(psi(2,3.3,l,w)/a)
