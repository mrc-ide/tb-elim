# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 14:19:04 2023

@author: ISA
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import lognorm, beta
from scipy.special import betaln

def get_dist2(prctiles, distribution, visualizing=False):
    # Sort the percentiles
    dat = sorted(prctiles)
    
    # Estimators to help set up initial guesses
    mn = dat[1]
    var = (dat[2] - dat[0])**2 / 4
    
    def objective(x):
        if distribution == 'lognorm':
            cdf_values = [lognorm.cdf(d, x[1], loc=0, scale=np.exp(x[0])) for d in dat]
        elif distribution == 'beta':
            cdf_values = [beta.cdf(d, x[0], x[1], loc=0, scale=1) for d in dat]
        return sum(((np.array(cdf_values) / (np.array([2.5, 50, 97.5]) / 100)) - 1)**2)
    
    # Initial guess for parameters
    mu = np.log(mn)
    sigma = np.sqrt(np.log(var / mn**2 + 1))
    
    # Adjusted initialization for beta
    if mn < 0.5:
        alpha_init = mn * (mn * (1 - mn) / var - 1)
        beta_init = (1 - mn) * (mn * (1 - mn) / var - 1)
    else:
        alpha_init = (2 - mn) * ((1 - mn) / var - 1)
        beta_init = (1 - mn) * ((2 - mn) * (1 - mn) / var - 1)
    
    if distribution == 'lognorm':
        init = [mu, sigma]
    elif distribution == 'beta':
        init = [alpha_init, beta_init]
    
    # Optimize to fit the distribution
    result = minimize(objective, init, method='Nelder-Mead', options={'maxiter': 1e6, 'maxfev': 1e6})
    
    if result.success:
        out = result.x
    else:
        raise Exception('Calibration setup not converged')
    
    if distribution == 'lognorm':
        mu, sigma = out[0], out[1]
        logfn = lambda x: -((np.log(x) - mu)**2 / (2 * sigma**2)) - np.log(x * sigma * np.sqrt(2 * np.pi))
    elif distribution == 'beta':
        a, b = out[0], out[1]
        logfn = lambda x: (a - 1) * np.log(x) + (b - 1) * np.log(1 - x) - betaln(a, b)
    
    aux = {'sim': np.array([lognorm.cdf(d, out[1], loc=0, scale=np.exp(out[0])) if distribution == 'lognorm'
                           else beta.cdf(d, out[0], out[1], loc=0, scale=1) for d in dat]),
           'val': result.fun}
    
    if visualizing:
        x = np.linspace(dat[0] * 0.8, dat[2] * 1.2, 100)
        if distribution == 'lognorm':
            y = lognorm.pdf(x, out[1], loc=0, scale=np.exp(out[0]))
        elif distribution == 'beta':
            y = beta.pdf(x, out[0], out[1], loc=0, scale=1)
        plt.figure()
        plt.plot(x, y)
        for d in dat:
            plt.axvline(x=d, linestyle='--', color='gray')
        plt.show()
    
    return logfn
