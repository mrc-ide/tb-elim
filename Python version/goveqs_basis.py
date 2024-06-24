# -*- coding: utf-8 -*-
"""


@author: ia19
"""

import numpy as np


def goveqs_basis3(t, insert, i, s, M, agg, sel, r, p):
    n_states = len(insert)
    out = np.zeros((n_states, 1))
    invec = insert[:i['nstates']]
    invec = invec.reshape(-1,1)

    # Adjust index for zero-based Python indexing
    nTPT_index = i['aux']['nTPT'][0] - 1  # Convert from one-based to zero-based index

    # Lambda with declining beta process over time (scalar)
    lam = M['lam'] @ invec / np.sum(invec) * (1 - p['betadec']) ** np.maximum((t - 2010), 0)
    
    # Full specification of the model
    allmat = M['lin'] + (lam * M['nlin'])
    
    # Model operations
    out[:i['nstates']] = np.dot(allmat, invec)
    morts = M['mort'] * invec
    out[:i['nstates']] -= np.sum(morts, axis=1).reshape(-1, 1)
    dom_morts = np.sum(morts[s['dom'], :])
    out[i[('U', 'dom')]] += dom_morts
    out[s['migr']] -= r['migr'] * invec[s['migr']] / np.sum(invec[s['migr']])
    
    # Migration in
    inmigr = np.sum(morts[s['migr'], :]) + r['migr']
    vec = np.array([
        1 - p['LTBI_in_migr'],
        (1 - p['migrTPT']) * p['LTBI_in_migr'] * 0.02,
        (1 - p['migrTPT']) * p['LTBI_in_migr'] * 0.98,
        p['migrTPT'] * p['LTBI_in_migr'] * 0.02,
        p['migrTPT'] * p['LTBI_in_migr'] * 0.98
    ])
    out[s['migrstates'], 0] += inmigr * vec.reshape(-1)

    # Correctly use the adjusted index for calculating nTPT related output
    if 0 <= nTPT_index < n_states:
        out[nTPT_index] = np.sum((sel['nTPT'] * allmat).dot(invec))
    else:
        raise IndexError(f"Adjusted index {nTPT_index} is out of bounds for 'out' with size {len(out)}.")

    # Auxiliary calculations
    out[i['aux']['inc']] = agg['inc'].dot(sel['inc'] * allmat).dot(invec)
    tmp1 = agg['sources'].dot((sel['Lf2I'] * allmat).dot(invec))
    tmp2 = agg['sources'].dot((sel['Pf2I'] * allmat).dot(invec))
    tmp3 = agg['sources'].dot((sel['Ls2I'] * allmat).dot(invec))
    tmp4 = agg['sources'].dot((sel['Ps2I'] * allmat).dot(invec))
    tmp5 = agg['sources'].dot((sel['R2I'] * allmat).dot(invec))
    out[i['aux']['sources']] = np.concatenate([tmp1, tmp2, tmp3, tmp4, tmp5])
    out[i['aux']['mort']] = np.sum(morts[:, 1])

    return out.flatten()





























# import numpy as np

# def goveqs_basis3(t, insert, i, s, M, agg, sel, r, p):
        
#     # initialise out vector used in odeint
#     out = np.zeros((len(insert), 1))
#     # select just the compartments, no aggregators or selectors
#     invec = insert[:i['nstates']]
#     invec = invec.reshape(-1,1)
    
#     # lambda with declining beta process over time
#     # will be scalar
#     lam = M['lam'] @ (invec) / np.sum(invec) * (1 - p['betadec']) ** np.maximum((t - 2010), 0)
    
#     # full specification of the model
#     # allmat is 33x33, invec.T is 33x1
#     allmat = M['lin'] + (lam * M['nlin'])
#     # again just the compartments of the model
#     # full model times state vec gives 22x1
#     out[:i['nstates']] = np.dot(allmat, invec)
    
#     # mortality
#     morts = M['mort'] * invec
#     # substract mortality from the states
#     out[:i['nstates']] -= np.sum(morts, axis=1).reshape(-1, 1)
    
#     # and births into UK pop
#     dom_morts = np.sum(morts[s['dom'], :])
#     out[i[('U', 'dom')]] += dom_morts
    
    
#     # migration out of the UK
#     out[s['migr']] = out[s['migr']] - r['migr'] * invec[s['migr']] / np.sum(invec[s['migr']])
    
#     # migration in
#     inmigr = np.sum(morts[s['migr'], :]) + r['migr']
#     vec = np.array([
#         1 - p['LTBI_in_migr'],
#         (1 - p['migrTPT']) * p['LTBI_in_migr'] * 0.02,
#         (1 - p['migrTPT']) * p['LTBI_in_migr'] * 0.98,
#         p['migrTPT'] * p['LTBI_in_migr'] * 0.02,
#         p['migrTPT'] * p['LTBI_in_migr'] * 0.98
#     ])
#     out[s['migrstates'], 0] += inmigr * vec.reshape(-1)  # Adjust reshape to match required dimensions

    
    
#     # auxillaries    
#     out[i['aux']['inc']] = agg['inc'].dot(sel['inc'] * allmat).dot(invec)
#     tmp1 = agg['sources'].dot((sel['Lf2I'] * allmat).dot(invec))
#     tmp2 = agg['sources'].dot((sel['Pf2I'] * allmat).dot(invec))
#     tmp3 = agg['sources'].dot((sel['Ls2I'] * allmat).dot(invec))
#     tmp4 = agg['sources'].dot((sel['Ps2I'] * allmat).dot(invec))
#     tmp5 = agg['sources'].dot((sel['R2I'] * allmat).dot(invec))
#     out[i['aux']['sources']] = np.concatenate([tmp1, tmp2, tmp3, tmp4, tmp5])
#     out[i['aux']['mort']] = np.sum(morts[:, 1])
#     print("Indices in i['aux']['nTPT']:", i['aux']['nTPT'])
#     print("Size of array 'out':", out.shape)
#     out[i['aux']['nTPT']] = np.sum((sel['nTPT'] * allmat).dot(invec))
    
#     return out











