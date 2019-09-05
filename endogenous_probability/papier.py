import numpy as np
import pandas as pd
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
def tau(e,g,l,w):

    return (((1-w**(1-g))*l*u)/(a*(1-g)))**(1/(1-u))


def psi(e,g,l,w):

    return e*r+(1-e)*(
        (1-tau(e,g,l,w))*a - (g*(s**2)/2) - l*(1+d-(tau(e,g,l,w))**u) * (1-w**(1-g))/(1-g)
        )


def trend_growth(e,g,l,w):
    
    return (1-tau(e,g,l,w))*a - psi(e,g,l,w)


def expected_growth(e,g,l,w):
    
    return trend_growth(e,g,l,w) - l*(1+d-(tau(e,g,l,w))**u) * (1-w)


def effect_disasters_expected_growth(e,g,l,w):
    
    return expected_growth(e,g,l,w) - expected_growth(e,g,0,w)


def lucas_measure(e,g,l,w):
    psi_optimal = psi(e,g,l,w)
    psi_bau = e*r+(1-e)*(
        (1-0)*a - (g*(s**2)/2) - l*(1+d-(0)**u) * (1-w**(1-g))/(1-g)
        )
    
    return (psi_optimal/psi_bau)**(1/(1-e)) - 1


def mrt_lambda_gdp(e,g,l,w):
    group_1 = (u**(u/(1-u))-u**(1/(1-u)))/(1-u)
    group_2 = ((1-w**(1-g))/((a**u)*(1-g)))**(1/(1-u))
    group_3 = (1+d)*(1-w**(1-g))/(1-g)
    
    return -(1/psi(e,g,l,w))*(
            l**(u/(1-u)) * group_1 * group_2 - group_3
            )


def mrt_omega_gdp(e,g,l,w):
    group_1 = (u**(u/(1-u))-u**(1/(1-u)))/(1-u)
    group_2 = ((1-w**(1-g))/(a*(1-g)))**(u/(1-u))
    
    return -(w**(-g))/psi(e,g,l,w)*(
            l*(1+d) - l**(1/(1-u)) * group_1 * group_2
            )
    

def rho_to_fit_growth(e,g,l,w):
    
    group_1 = a*(1-tau(e,g,l,w))
    group_2 = l*(1+d-(tau(e,g,l,w))**u) * (1-w)
    group_3 = g*(s**2)/2
    group_4 = l*(1+d-(tau(e,g,l,w))**u) * (1-w**(1-g)) / (1-g)

    return 1/e * (group_1 - group_2 - (1-e)*(group_1 - group_3 - group_4) - 0.02)


# Fix parameters' value :
d = 1
s = 0.02
u = 0.25
a = 0.069
e = 1.5
g = 4

w_1 = 0.948
w_2 = 0.85
w_3 = 0.60
l_1 = 0.0307
l_2 = 0.01064
l_3 = 0.003991


### Functions to activate :
main_parameters = True
table_tau = False
table_lucas = False
table_mrs_lambda = False
table_mrs_omega = False
heatmap = True
figure_tau = False
figure_lucas = False
figure_mrs = False


### Output values:

# Main output of interest:
if main_parameters == True:
    l = l_1
    w = w_1
    r = rho_to_fit_growth(e,g,l,w)
    print 'Risk mitigation share:', tau(e,g,l,w)*100, '%'
    print 'Reduction in risk of env. disaster:', tau(e,g,l,w)**u*100, '%'
    print 'Share of production consumed', psi(e,g,l,w)/a*100, '%'
    print 'Expected growth:', expected_growth(e,g,l,w)*100, '%'
    print 'Expected damages from env. dis.:', l_1*(1-(tau(e,g,l,w))**u) * (1-w) * 100, '%'
    del l,w

# Alternative calibrations for tau
if table_tau == True:
    for gamma in [1+1e-09, 2, 4, 10]:
        r = rho_to_fit_growth(e,gamma,l_1,w_1)
        print 'with g=', gamma, ', w=', w_1, 'and l=', l_1, ':', tau(e,gamma,l_1,w_1)*100, '%'
        r = rho_to_fit_growth(e,gamma,l_2,w_2)
        print 'with g=', gamma, ', w=', w_2, 'and l=', l_2, ':', tau(e,gamma,l_2,w_2)*100, '%'
        r = rho_to_fit_growth(e,gamma,l_3,w_3)
        print 'with g=', gamma, ', w=', w_3, 'and l=', l_3, ':', tau(e,gamma,l_3,w_3)*100, '%'


### Welfare effects:

# Lucas' measure:
if table_lucas == True:
    for gamma in [1+1e-12, 2, 4, 10]:
        r = rho_to_fit_growth(e,gamma,l_1,w_1)
        print 'with g=', gamma, ', w=', w_1, 'and l=', l_1, ':', lucas_measure(e,gamma,l_1,w_1)*100, '%'
        r = rho_to_fit_growth(e,gamma,l_2,w_2)
        print 'with g=', gamma, ', w=', w_2, 'and l=', l_2, ':', lucas_measure(e,gamma,l_2,w_2)*100, '%'
        r = rho_to_fit_growth(e,gamma,l_3,w_3)
        print 'with g=', gamma, ', w=', w_3, 'and l=', l_3, ':', lucas_measure(e,gamma,l_3,w_3)*100, '%'


# GDP equivalent in welfare of a decrease in lambda:
if table_mrs_lambda == True:
    for gamma in [1+1e-09, 2, 4, 10]:
        r = rho_to_fit_growth(e,gamma,l_1,w_1)
        print 'with g=', gamma, ', w=', w_1, 'and l=', l_1, ':', mrt_lambda_gdp(e,gamma,l_1,w_1)
        r = rho_to_fit_growth(e,gamma,l_2,w_2)
        print 'with g=', gamma, ', w=', w_2, 'and l=', l_2, ':', mrt_lambda_gdp(e,gamma,l_2,w_2)
        r = rho_to_fit_growth(e,gamma,l_3,w_3)
        print 'with g=', gamma, ', w=', w_3, 'and l=', l_3, ':', mrt_lambda_gdp(e,gamma,l_3,w_3)
    
    for epsilon in [0.25, 0.5, 1+1e-09, 1.5, 2.0]:
        r = rho_to_fit_growth(epsilon,g,l_1,w_1)
        print 'with e=', epsilon, ', w=', w_1, 'and l=', l_1, ':', mrt_lambda_gdp(epsilon,g,l_1,w_1)
        r = rho_to_fit_growth(epsilon,g,l_2,w_2)
        print 'with e=', epsilon, ', w=', w_2, 'and l=', l_2, ':', mrt_lambda_gdp(epsilon,g,l_2,w_2)
        r = rho_to_fit_growth(epsilon,g,l_3,w_3)
        print 'with e=', epsilon, ', w=', w_3, 'and l=', l_3, ':', mrt_lambda_gdp(epsilon,g,l_3,w_3)

# GDP equivalent in welfare of an increase in omega:
if table_mrs_omega == True:
    for gamma in [1+1e-09, 2, 4, 10]:
        r = rho_to_fit_growth(e,gamma,l_1,w_1)
        print 'with g=', gamma, ', w=', w_1, 'and l=', l_1, ':', mrt_omega_gdp(e,gamma,l_1,w_1)
        r = rho_to_fit_growth(e,gamma,l_2,w_2)
        print 'with g=', gamma, ', w=', w_2, 'and l=', l_2, ':', mrt_omega_gdp(e,gamma,l_2,w_2)
        r = rho_to_fit_growth(e,gamma,l_3,w_3)
        print 'with g=', gamma, ', w=', w_3, 'and l=', l_3, ':', mrt_omega_gdp(e,gamma,l_3,w_3)
    
    for epsilon in [0.25, 0.5, 1+1e-09, 1.5, 2.0]:
        r = rho_to_fit_growth(epsilon,g,l_1,w_1)
        print 'with e=', epsilon, ', w=', w_1, 'and l=', l_1, ':', mrt_omega_gdp(epsilon,g,l_1,w_1)
        r = rho_to_fit_growth(epsilon,g,l_2,w_2)
        print 'with e=', epsilon, ', w=', w_2, 'and l=', l_2, ':', mrt_omega_gdp(epsilon,g,l_2,w_2)
        r = rho_to_fit_growth(epsilon,g,l_3,w_3)
        print 'with e=', epsilon, ', w=', w_3, 'and l=', l_3, ':', mrt_omega_gdp(epsilon,g,l_3,w_3)


### Heatmaps - Effect of disasters on expected growth:
if heatmap == True:
    # Create the grid
    e_inverse = 1/(arange(4,1,-0.01) - 0.001)
    e_normal = arange(1,3,0.01) + 0.001
    e = np.concatenate([e_inverse,e_normal])
    g = arange(6,1,-0.01) + 0.005
    E,G = meshgrid(e, g)
    
    # Prepare axes
    axe_g = ['6', '5', '4', '3', '2', '1']
    axe_e = ['1/4', '1/3', '1/2', '1', '2', '3']
    
    # Build heatmap
    values_heatmap = effect_disasters_expected_growth(E, G, l_3, w_3) # evaluation of the function on the grid
    
    fig, ax = plt.subplots()
    heatmap_expected_growth = ax.imshow(
        values_heatmap, norm=norm, cmap=plt.cm.seismic, extent = [0,10,10,0], interpolation='none'
        )
    ax.set_xticklabels(axe_e)
    ax.set_yticklabels(axe_g)
    plt.xlabel('IES')
    plt.ylabel('RRA')
    cbar = plt.colorbar(heatmap_expected_growth)
    cbar.set_label('percentage points')
    plt.show()


### Plot Welfare effects:
e = 1.5
dataframe_welfare_effects = pd.DataFrame(index = range(1, 500),
            columns = ['gamma_index', 'gamma',
                       'rho_1', 'rho_2', 'rho_3', 'tau_1', 'tau_2', 'tau_3',
                       'lucas_measure_1', 'mrt_lambda_gdp_1', 'mrt_omega_gdp_1',
                       'lucas_measure_2', 'mrt_lambda_gdp_2', 'mrt_omega_gdp_2',
                       'lucas_measure_3', 'mrt_lambda_gdp_3', 'mrt_omega_gdp_3',
                       ]
            )
for gamma_index in range(1,500):
    gamma_index = 1 + gamma_index
    gamma = 1 + float(gamma_index)/100

    dataframe_welfare_effects['gamma_index'][gamma_index-1] = gamma_index
    dataframe_welfare_effects['gamma'][gamma_index-1] = gamma

    for disaster_type in range(1,4):
        if disaster_type == 1:
            l = l_1
            w = w_1
        elif disaster_type == 2:
            l = l_2
            w = w_2
        else:
            l = l_3
            w = w_3
            
        r = rho_to_fit_growth(e,gamma,l,w)
        dataframe_welfare_effects['tau_{}'.format(disaster_type)][gamma_index-1] = tau(e,gamma,l,w)*100
        dataframe_welfare_effects['rho_{}'.format(disaster_type)][gamma_index-1] = r
        dataframe_welfare_effects['lucas_measure_{}'.format(disaster_type)][gamma_index-1] = lucas_measure(e,gamma,l,w) * 100
        dataframe_welfare_effects['mrt_lambda_gdp_{}'.format(disaster_type)][gamma_index-1] = mrt_lambda_gdp(e,gamma,l,w)
        dataframe_welfare_effects['mrt_omega_gdp_{}'.format(disaster_type)][gamma_index-1] = mrt_omega_gdp(e,gamma,l,w)


if figure_tau == True:
    for disaster_type in range(1,4):
        #plt.title("Optimal share of income spent in policy instrument (in %)")
        #plt.plot(dataframe_welfare_effects['gamma'],
        #         dataframe_welfare_effects[['tau_1'] + ['tau_2'] + ['tau_3']])
        plt.plot(dataframe_welfare_effects['gamma'],
                 dataframe_welfare_effects[['tau_{}'.format(disaster_type)]])
        plt.xlabel('RRA')
        plt.ylabel('Tau')
        plt.show()

if figure_lucas == True:
    for disaster_type in range(1,4):
        #plt.title("Welfare benefits of the policy")
        #plt.plot(dataframe_welfare_effects['gamma'],
        #         dataframe_welfare_effects[['lucas_measure_1'] + ['lucas_measure_2'] + ['lucas_measure_3']])
        plt.plot(dataframe_welfare_effects['gamma'],
                 dataframe_welfare_effects[['lucas_measure_{}'.format(disaster_type)]])
        plt.xlabel('RRA')
        plt.ylabel('Gamma')
        plt.show()

if figure_mrs == True:
    for disaster_type in range(1,4):
        #plt.title("Marginal rate of substitution between proportionate changes in GDP and in disasters probability")
        #plt.plot(dataframe_welfare_effects['gamma'],
        #         dataframe_welfare_effects[['mrt_lambda_gdp_1'] + ['mrt_lambda_gdp_2'] + ['mrt_lambda_gdp_3']])
        plt.plot(dataframe_welfare_effects['gamma'],
                 dataframe_welfare_effects[['mrt_lambda_gdp_{}'.format(disaster_type)]])
        plt.xlabel('RRA')
        plt.ylabel('MRS')
        plt.show()
