# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 16:51:37 2023

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
from make_model import make_model
from allocate import allocate_parameters
from setup_model import likelihood_function
from goveqs_basis import goveqs_basis3
from gov_eqs_scaleup import goveqs_scaleup

def get_objective(x, ref, prm, gps, calfn):
    
    # extract dictionaries
    i = ref['i']
    s = ref['s']
    xi = ref['xi']
    p = prm['p']
    r = prm['r']
    sel = prm['sel']
    agg = prm['agg']
    
    
    # defined final parameters to test code - will be deleted after
    #x = [23.9410, 0.1473, 5.2634, 0.7268, 12.3779]
    # assign parameters to p and r dicts
    p, r = allocate_parameters(x, p, r, xi)
    
    # Initialize lists for first and second values
    first_values = []
    second_values = []
    
    # Iterate through each key in bounds
    for key, bounds in prm['bounds'].items():
        if isinstance(bounds[0], list):
            # Extend lists if bounds are lists (e.g., gamma)
            first_values.extend([b[0] for b in bounds])
            second_values.extend([b[1] for b in bounds])
        else:
            # Append normally if not lists
            first_values.append(bounds[0])
            second_values.append(bounds[1])
    
    # Convert to numpy arrays and ensure dtype is set to float for calculations
    first_values = np.array(first_values, dtype=float)
    second_values = np.array(second_values, dtype=float)
    
    # Reconstruct tmp1 with proper dimensions
    tmp1 = np.array([first_values, second_values])
    
    # Ensure x is reshaped correctly and matches the length of bounds
    x_reshaped = np.array([x])  # x must be the correct length to match tmp1
    
    
    # Attempt the stacking
    tmp1 = np.vstack((tmp1, x_reshaped))
    # Calculate differences between consecutive rows for columns 0, 2, and 1 (assuming 0-based indexing)
    tmp2 = np.diff(tmp1[[0, 2, 1], :], axis=0)
    
    # Check if the minimum value in tmp2 is less than 0
    cond1 = np.min(tmp2) < 0
    
    if cond1:
        out = -np.inf
        aux = np.nan
    else:
        # creates model with params
        M = make_model(p, r, i, s, gps)
        
        # initialise the 'insert' for gov eqs basis
        # (initial state vector with aux and selectors on end)
        init = np.zeros(i['nx'])
        seed = 1e-5
        init[i[('U', 'dom')]] = 1 - 0.168 - seed
        init[i[('U', 'mig_recent')]] = 0.168
        init[i[('I', 'dom')]] = seed
        
        # Adjust parameters as needed
        p0 = p.copy(); r0 = r.copy()
        p0['betadec'] = 0
        r0['gamma'] = r0['gamma_2015']
        M0 = make_model(p0, r0, i, s, gps)
        
        def geq0(y, t, i, s, M0, agg, sel, r0, p0):
            return goveqs_basis3(t, y, i, s, M0, agg, sel, r0, p0)
        
        t0 = np.linspace(0, 5e3, 500)  # Time vector, adjust granularity as needed
        soln0 = odeint(geq0, init, t0, args=(i, s, M0, agg, sel, r0, p0))
       

        # For subsequent periods
        p1 = p0.copy(); r1 = r0.copy()
        r1['TPT'] = [0, r['TPT2020rec'], 0]
        M1 = make_model(p1, r1, i, s, gps)
    
        p2 = p.copy(); r2 = r.copy()
        r2['gamma'] = r1['gamma_2020']
        M2 = make_model(p2, r2, i, s, gps)
        
        def geq1(y, t, M0, M1, M2, timeline, i, s, p2, r2, prm, sel, agg):
            
            return goveqs_scaleup(t, y, M0, M1, M2, timeline, i, s, p2, r2, prm, sel, agg)

        t1 = np.linspace(2010, 2020, 100)  # Adjust as needed for granularity
        soln1 = odeint(geq1, soln0[-1], t1, args=(M0, M1, M2, [2015, 2020, 2010, 2020], i, s, p2, r2, prm, sel, agg))
        #soln1 = odeint(geq1, soln0[-1], t1, args=(M0, M1, M2, [2015, 2020, 2010, 2020], i, s, p2, r2, prm, sel, agg))
        
        # Calculate the derivative of the solution with respect to time
        dsol = np.diff(soln1, axis=0)  # Assuming soln1.y for solve_ivp; adjust if using odeint
        sfin = soln1[-1, :]  # Final state
        
        # Indexing for year 2010 might require finding the closest time point in t1 to 2010
        # This is straightforward if t1 contains 2010 exactly; otherwise, you might need to find the nearest value
        index_2010 = np.where(t1 == 2010)[0][0]  # Adjust accordingly if t1 is not an exact match
        index_2020 = np.where(t1 == 2020)[0][0]  # Assuming 2020 is directly in t1
        
        # Assuming i.aux.inc and i.aux.mort are previously defined, pointing to the right indices
        # For the sake of demonstration, let's say they are defined as:
        # i = {'aux': {'inc': [index_for_inc], 'mort': index_for_mort}}
        
        incd2010 = dsol[index_2010, i['aux']['inc'][0]] * 1e5
        incd2020 = dsol[-1, i['aux']['inc'][0]] * 1e5
        incd = dsol[-1, i['aux']['inc']] * 1e5  # Assuming i['aux']['inc'] is a list of indices
        mort = dsol[-1, i['aux']['mort']] * 1e5
        
        # Calculating proportions, assume 'p' and 's' contain necessary info
        # For p_LTBI, directly using the parameter as no calculation shown
        p_LTBI = p['LTBI_in_migr']
        
        # For p_migrpopn, you need total in migrant states vs. total in all states
        # Assuming s['migr'] contains indices for migrant states and using sfin for final state distribution
        p_migrpopn = np.sum(sfin[s['migr']]) / np.sum(sfin[:i['nstates']])
        p_migrTB = incd[2] / incd[0]

        # Number initiating TPT in 2019
        # Again, this assumes that dsol, i.aux.nTPT, etc. are appropriately defined
        #n_TPT2019 = dsol[-1, i['aux']['nTPT']] * 1e5
    # Iterate over each index in 'nTPT' and calculate the sum
        n_TPT_sum = 0
        for nTPT_index in i['aux']['nTPT']:
            # Adjust index for zero-based Python indexing
            nTPT_index -= 1
            if 0 <= nTPT_index < dsol.shape[1]:
                n_TPT_sum += dsol[-1, nTPT_index] * 1e5
            else:
                raise IndexError(f"Index {nTPT_index} is out of bounds for 'dsol' with shape {dsol.shape}")
    
        
        if np.any(incd > 0.1):
            #out = calfn(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI)
            out = likelihood_function(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI)
            # add to the auxillaries with calculated values
            aux = {
                'soln': soln1,
                'incd': dsol[np.where(t1 == 2010)[0][0]:, i['aux']['inc'][0]] * 1e5,
                'incd2010': incd2010,
                'incd2020': incd2020,
                'mort': mort,
                'p_migrTB': p_migrTB,
                'p_migrpopn': p_migrpopn,
                'p_LTBI': p_LTBI,
                'p_migrect': np.sum(sfin[s['mig_recent']]) / np.sum(sfin[:i['nstates']]),
                'nTPT': n_TPT_sum
            }
        else:
            out = -np.inf
            aux = np.nan
        
    return out, aux


        # # wrapper for gov eqs basis taking only insert and time
        # # def geq(t, in_):
        # #     return goveqs_basis2(t, in_, i, s, M, agg, sel, r, p).flatten()
        
        # # time range for solving equation until
        # t0 = np.arange(2020)
        # soln0 = odeint(geq, init, t0, tfirst=True)
        # #soln0 = solve_ivp(geq, (t0[0], t0[-1]), init, t_eval=t0, vectorized=True)
        # #soln0=soln0.y.T

        
        # # calculates the differences between rows of soln0
        # dsol = np.diff(soln0, axis=0)
        # # filters rows of dsol, the t0 index select rows where the corresponding 
        # # time point in t0 is equal to 2010 
        # # aux index is position of incidence
        # incd2010 = dsol[t0[:-1] == 2010, i['aux']['inc'][0]] * 1e5
        # # final position in row of dsol gives 2020, and aux index for incidence
        # incd2020 = dsol[-1, i['aux']['inc'][0]] * 1e5
        # # incidence in all positions of the aux for 2020
        # incd = dsol[-1, i['aux']['inc']] * 1e5
        # # final mortality (2020) and its index at end of aux
        # mort = dsol[-1, -1] * 1e5
        
        # # extract final row of soln0, which has the solutions at the last time point
        # # (selecting last row 2020, and all 31 columns)
        # sfin = soln0[-1, :]
        # # within this subset, look for latent foreign states, to get proportion of latent
        # p_LTBI = np.sum(sfin[np.intersect1d(s['for'], [s['Lf'], s['Ls']])]) / np.sum(sfin[s['for']])
        # # proportion of migrants is all foreign compartments over total pop
        # p_migrpopn = np.sum(sfin[s['for']]) / np.sum(sfin[:i['nstates']])
        
        # # number of tpt
        # #n_tpt2019 = sfin[s['Pf']] + sfin[s['Ps']] * 1e5
        # n_tpt2019 = (soln0[np.where(t0 == 2019)[0][0], s['Pf']] - soln0[np.where(t0 == 2018)[0][0], s['Pf']]) + \
        #     (soln0[np.where(t0 == 2019)[0][0], s['Ps']] - soln0[np.where(t0 == 2018)[0][0], s['Ps']]) * 1e5
