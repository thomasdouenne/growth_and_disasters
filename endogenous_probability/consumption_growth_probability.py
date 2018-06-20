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
    """
    Generate an instance of Brownian motion (i.e. the Wiener process):

        X(t) = X(0) + N(0, delta**2 * t; 0, t)

    where N(a,b; t0, t1) is a normally distributed random variable with mean a and
    variance b.  The parameters t0 and t1 make explicit the statistical
    independence of N on different time intervals; that is, if [t0, t1) and
    [t2, t3) are disjoint intervals, then N(a, b; t0, t1) and N(a, b; t2, t3)
    are independent.
    
    Written as an iteration scheme,

        X(t + dt) = X(t) + N(0, delta**2 * dt; t, t+dt)


    If `x0` is an array (or array-like), each value in `x0` is treated as
    an initial condition, and the value returned is a numpy array with one
    more dimension than `x0`.

    Arguments
    ---------
    x0 : float or numpy array (or something that can be converted to a numpy array
         using numpy.asarray(x0)).
        The initial condition(s) (i.e. position(s)) of the Brownian motion.
    n : int
        The number of steps to take.
    dt : float
        The time step.
    delta : float
        delta determines the "speed" of the Brownian motion.  The random variable
        of the position at time t, X(t), has a normal distribution whose mean is
        the position at time t=0 and whose variance is delta**2*t.
    out : numpy array or None
        If `out` is not None, it specifies the array in which to put the
        result.  If `out` is None, a new numpy array is created and returned.

    Returns
    -------
    A numpy array of floats with shape `x0.shape + (n,)`.
    
    Note that the initial value `x0` is not included in the returned array.
    """

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
def theta(e,g,l,w):

    return (((1-w**(1-g))*l*u)/(a*(1-g)))**(1/(1-u))


def psi(e,g,l,w):

    return e*r+(1-e)*(
        (1-theta(e,g,l,w))*a - (g*s/2) - l*(1+d-(theta(e,g,l,w))**u) * (1-w**(1-g))/(1-g)
        )


def trend_growth(e,g,l,w):
    
    return ((1-theta(e,g,l,w))*a - psi(e,g,l,w))


# Fix parameters' value :
d = 3
w = 0.79
l = 0.01
s = 0.02
u = 0.5
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
    liste.append(1* (random.uniform(0, 1) < (l*(1+d-theta(2,3.3,l,w)**u)*dt)))
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
