# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 21:21:49 2023

@author: ISA
"""

import numpy as np
#from obj import get_objective
#from setup_model import r,p,agg,sel,ref,xi,prm,gps_born,likelihood_function
from tqdm import tqdm
import time

def MCMC_adaptive(F, x0, n, sigma, cov0, displ=True):
    d = len(x0)
    b = 0.05
    sd = sigma * 2.4**2 / d

    # Checks on the initial covariance matrix
    if cov0 is None:
        cov0 = np.eye(d)
        
    n = int(n)  # Convert n to an integer
    #xsto = np.zeros((n, d))
    xsto = np.zeros((n, d))
    outsto = np.zeros(n)
    history = np.zeros((n, d + 1))

    xsto[0] = x0
    xbar = xsto.copy()
    FX = F(x0)
    outsto[0] = FX
    acc = 0

    if displ:
        print("Starting MCMC...")

    start_time = time.time()  # Start the timer

    for t in tqdm(range(1, n), total=n-1, desc="MCMC Progress", ncols=100):
        X = xsto[t - 1]

        # Make a proposal from the distribution
        Y0 = np.random.multivariate_normal(X, 0.1**2 * cov0 * sigma / d)
        if t < 2 * d:
            Y = np.maximum(Y0, 0)
        else:
            covmat = np.cov(xsto[:t, :].T)
            # Incorporate any adjustments to the covariance matrix here, if needed
            covmat = (covmat + covmat.T) / 2
            Y = np.maximum((1-b) * np.random.multivariate_normal(X, sd * covmat) + b * Y0, 0)

        history[t, :d] = Y

        FY = F(Y)
        if np.random.rand() < np.exp(FY - FX) and np.isfinite(FY):
            # Accept
            xsel = Y
            FX = FY
            acc += 1
            history[t, -1] = 1
        else:
            # Reject
            xsel = xsto[t - 1]

        xsto[t] = xsel
        outsto[t] = FX
        xbar[t] = (xbar[t - 1] * t + xsel) / (t + 1)
        
        end_time = time.time()  # End the timer

        elapsed_time = end_time - start_time  # Calculate elapsed time

        # Display options
        if displ:
            print(f"MCMC completed in {elapsed_time:.2f} seconds.")

    accept_rate = acc / n

    return xsto, outsto, history, accept_rate
