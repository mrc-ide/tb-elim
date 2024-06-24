# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 18:54:00 2023

@author: ISA
"""
from get_addresses import get_addresses
import numpy as np
from scipy.sparse import csr_matrix
from get_dist2 import get_dist2


# combining states and born will be used as groups in get_addresses
states = ['U', 'Lf', 'Ls', 'Pf', 'Ps', 'I', 'I2', 'Tx', 'Rlo', 'Rhi', 'R'] # states
gps_born = ['dom', 'mig_recent', 'mig_long'] # where they are born



i, s, d, lim = get_addresses([states, gps_born])
s['everyI'] = np.concatenate((s['I'], s['I2']))
s['migr'] = np.concatenate((s['mig_recent'], s['mig_long']))
s['migrstates'] = [i['U', 'mig_recent'], i['Lf','mig_recent'], i['Ls','mig_recent'], i['Pf','mig_recent'], i['Ps','mig_recent']]

# s['migrL'] = [s['Lf'][1], s['Ls'][1]]
# s['migrP'] = [s['Pf'][1], s['Ps'][1]]

# ~~~~~~~~~~~~ auxillary time
auxillaries = ['inc', 'sources', 'mort', 'nTPT'] # set up to add to end of matrix
lengths = [3, 15, 1, 1]
lim=i['nstates']
i['aux'] = {} # initialise auxes to end of i

for ii in range(len(auxillaries)): # loop over however many auxillaries there are
    # in each iteration of the loop create a list of integers starting 
    # from lim + 1 and ending at lim + lengths[ii] + 1. lim is keeps track of  
    # the last index used, and lengths[ii] is the length associated 
    # with the current auxiliary.
    inds = list(range(lim + 1, lim + lengths[ii] + 1))
    # below assign the indices 'inds' to their corresponding aux 
    # keys in i. e.g. , if auxillaries[ii] at beginning where ii==0
    # is 'inc', the indices associated with 'inc' are stored as a list in i['aux']['inc']
    i['aux'][auxillaries[ii]] = inds
    # lim is updated to last index in inds to make sure the next aux
    # category will start from the next available index.
    lim = inds[-1]
i['nx'] = lim


# store selectors and aggregators
sel = {}
agg = {}


# --- Incidence
# --- Incidence
tmp = np.zeros((3, i['nstates']))
tmp[0, s['everyI']] = 1
tmp[1, np.intersect1d(s['everyI'], s['dom'])] = 1
tmp[2, np.intersect1d(s['everyI'], s['migr'])] = 1
agg['inc'] = tmp


tmp = np.zeros((i['nstates'], i['nstates']))
tmp[s['everyI'], :] = 1
# Remove transitions due to changing migrant status
tmp[s['mig_recent'], s['mig_long']] = 0
tmp[s['mig_long'], s['mig_recent']] = 0

# Subtracting the diagonal
sel['inc'] = tmp - np.diag(np.diag(tmp)) # so that diagonal self to self terms arent counted


#~~~~~~~~~incidence sources

tmp = np.zeros((3, i['nstates']))
tmp[0, np.intersect1d([s['I'], s['I2']], s['dom'])] = 1
tmp[1, np.intersect1d([s['I'], s['I2']], s['mig_recent'])] = 1
tmp[2, np.intersect1d([s['I'], s['I2']], s['mig_long'])] = 1
agg['sources'] = tmp


# From recent infection
tmp = np.zeros((i['nstates'], i['nstates']))
for i_idx in s['everyI']:
    for j_idx in s['Lf']:
        tmp[i_idx, j_idx] = 1
sel['Lf2I'] = tmp - np.diag(np.diag(tmp))

# From recent infection, TPT
tmp = np.zeros((i['nstates'], i['nstates']))
for i_idx in s['everyI']:
    for j_idx in s['Pf']:
        tmp[i_idx, j_idx] = 1
sel['Pf2I'] = tmp - np.diag(np.diag(tmp))

# From remote infection
tmp = np.zeros((i['nstates'], i['nstates']))
for i_idx in s['everyI']:
    for j_idx in s['Ls']:
        tmp[i_idx, j_idx] = 1
sel['Ls2I'] = tmp - np.diag(np.diag(tmp))

# From remote infection, TPT
tmp = np.zeros((i['nstates'], i['nstates']))
for i_idx in s['everyI']:
    for j_idx in s['Ps']:
        tmp[i_idx, j_idx] = 1
sel['Ps2I'] = tmp - np.diag(np.diag(tmp))

# From relapse
tmp = np.zeros((i['nstates'], i['nstates']))
for i_idx in s['everyI']:
    for j_idx in np.concatenate([s['R'], s['Rlo'], s['Rhi']]):
        tmp[i_idx, j_idx] = 1
sel['R2I'] = tmp - np.diag(np.diag(tmp))

#~~~~~~~~~migrant TPT uptake
tmp = np.zeros((i['nstates'], i['nstates']))
for i_idx in np.concatenate([s['Pf'], s['Ps']]):
    for j_idx in np.concatenate([s['Lf'], s['Ls']]):
        tmp[i_idx, j_idx] = 1
sel['nTPT'] = tmp - np.diag(np.diag(tmp))



# ~~~~~~~~~~~ Natural history params
# all rates
r = {
     'gamma': 0.1,
    'progression0': 0.0826,
    'LTBI_stabil': 0.872,
    'reactivation0': 0.0006,
    'Tx': 2,
    'default': 0.01,
    'self_cure': 1/6,
    'relapse': [0.032, 0.14, 0.0015],
    'muTB': 1/6
    }

# all proportions
p = {
    'imm': 0.8,
    'migrTPT': 0,
    'TPTeff': 0.6
    }


r['TPT2020rec']   = 0.004
r['TPT']          = [0,0,0]
r['ACF']          = [0, 0,0]
r['ACF2']         = [0, 0,0]



#~~~~~~~~~~~~~~~~~ parameters

# free_params = ['beta', 'betadec', 'gamma', 'r_TPT2020rec', 'p_relrate', 'r_migr', 'p_LTBI_in_migr']
# param_lengths = [1,         1,         2,       1,          1,         1,               1]

free_params = ['beta', 'betadec', 'gamma', 'p_relrate', 'r_migr', 'p_LTBI_in_migr']
param_lengths = [1,         1,         2,           1,         1,               1]

# give indices to parameters for later allocation
limit = -1
xi = {}

for ii in range(len(free_params)): #iterate over each free param
    # calculate list of indices for current param. in this case, it would 
    # start with zero (limit + 1) and ends with limit + param_lengths[ii] + 1, 
    # where param_lengths[ii] is how many indices it occupies (1 each here). 
    # essentially assigns a range of indices to each parameter, ensuring that they do not overlap.
    indices = list(range(limit + 1, limit + param_lengths[ii] + 1))
    # matches param name with list of indices calculated in previous step. 
    # stores this mapping in the xi dictionary.
    xi[free_params[ii]] = indices
    # updates the limit to the index where it left off. 
    # making sure the next param will receive indices that dont overlap
    limit = indices[-1]

prm = {}
prm['bounds'] = {
    'beta': [0, 40],
    'betadec': [0, 0.15],
    'gamma': [[0, 10], [0, 10]],
    'p_relrate': [1, 20],
    'r_migr': [0, 1],
    'p_LTBI_in_migr': [0, 0.5]}

prm['p'] = p
prm['r'] = r
prm['agg'] = agg
prm['sel'] = sel


ref = {}
ref['i'] = i
ref['s'] = s
ref['xi'] = xi


# ~~~~~~~~~~~~~~~~~~~~~data 

data = {
    'incd2010': [14.1, 14.6, 15.1],
    'incd2020': [6.5, 7, 7.5],
    'mort': [0.28, 0.3, 0.32],
    'p_migrTB': [0.744, 0.764, 0.784],
    'p_migrpopn': [0.138, 0.168, 0.198],
    'p_LTBI': [0.15, 0.2, 0.25],
    'nTPT2019': [1.3*x for x in [0.9, 1, 1.1]]
}

f1a = get_dist2(data['incd2010'], 'lognorm')
f1b = get_dist2(data['incd2020'], 'lognorm')
f2 = get_dist2(data['mort'], 'lognorm')
f3 = get_dist2(data['p_migrTB'], 'beta')
f4 = get_dist2(data['p_migrpopn'], 'beta')
f5 = get_dist2(data['p_LTBI'], 'beta')
f6 = get_dist2(data['nTPT2019'], 'lognorm')

# Define the likelihood function
def likelihood_function(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI):
    likelihood = (
        np.sum(f1a(np.array(incd2010))) +
        np.sum(f1b(np.array(incd2020))) +
        np.sum(f2(np.array(mort))) +
        np.sum(f3(np.array(p_migrTB))) +
        np.sum(f4(np.array(p_migrpopn))) +
        np.sum(f5(np.array(p_LTBI)))
        #np.sum(f6(np.array(nTPT2019)))
    )
    return likelihood

# Example usage with data as lists or arrays
# incd2010 = [14.1, 14.6, 15.1]
# incd2020 = [6.5, 7, 7.5]
# mort = [0.28, 0.3, 0.32]
# p_migrTB = [0.708, 0.728, 0.748]
# p_migrpopn = [0.138, 0.168, 0.198]
# p_LTBI = [0.15, 0.2, 0.25]

likelihood = likelihood_function(data['incd2010'], data['incd2020'], data['mort'], data['p_migrTB'], data['p_migrpopn'], data['p_LTBI'])
#print(likelihood)

