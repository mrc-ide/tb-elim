# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 22:20:56 2023

@author: ISA
"""

from setup_model import states,gps_born
import numpy as np
from scipy.sparse import dia_matrix, csr_matrix


def get_states_for_born(i, born):
    state_values = {} # dict to store state values based on where they're born
    for state in states: # iterate over elements in this dict of states
        key = (state, born) # assign a key to each combo of state and born, to use as lookup in i
        state_values[state] = i.get(key) # for each state, use the key to look it up in i 
    return state_values



def make_model(p, r, i, s, gps):
    M = {}
    m = np.zeros((i['nstates'], i['nstates'])) # construct matrix
    for ib, born in enumerate(gps_born):
        state_values = get_states_for_born(i, born)
        
        # Access the state values for the current 'born'
        U = state_values['U']
        Lf = state_values['Lf']
        Ls = state_values['Ls']
        Pf = state_values['Pf']
        Ps = state_values['Ps']
        I = state_values['I']
        I2 = state_values['I2']
        Tx = state_values['Tx']
        Rlo = state_values['Rlo']
        Rhi = state_values['Rhi']
        R = state_values['R']
    
        # Progression from 'fast' latent
        source = Lf
        destin = I
        rate = r['progression'][ib]  # Replace with a random number
        m[destin, source] += rate
    
        source = Pf
        destin = I2
        rate = r['progression'][ib] * (1 - p['TPTeff'])  # Replace with a random number
        m[destin, source] += rate
    
        # Stabilization of 'fast' to 'slow' latent
        source = Lf
        destin = Ls
        rate = r['LTBI_stabil']  # Replace with a random number
        m[destin, source] += rate
    
        source = Pf
        destin = Ps
        rate = r['LTBI_stabil']  # Replace with a random number
        m[destin, source] += rate
    
        # Reactivation of 'slow' latent
        source = Ls
        destin = I
        rate = r['reactivation'][ib]  # Replace with a random number
        m[destin, source] += rate
    
        source = Ps
        destin = I
        rate = r['reactivation'][ib] * (1 - p['TPTeff'])  # Replace with a random number
        m[destin, source] += rate
    
        # Initiation of treatment
        source = I
        destins = [Tx, Rhi]
        rates = [r['gamma'], r['self_cure']]  # Replace with random numbers
        m[destins, source] += rates
    
        source = I2
        destins = [Tx, Rhi]
        rates = [r['gamma'], r['self_cure']]  # Replace with random numbers
        m[destins, source] += rates
    
        # Treatment completion or interruption
        source = Tx
        destins = [Rlo, Rhi]
        rates = [r['Tx'], r['default']]  # Replace with random numbers
        m[destins, source] += rates
    
        # Relapse
        sources = [Rlo, Rhi, R]
        destin = I2
        rates = r['relapse']  # Replace with a random number
        m[destin, sources] += rates
    
        # Stabilization of relapse risk
        sources = [Rlo, Rhi]
        destin = R
        rates = 0.5  # Replace with a random number
        m[destin, sources] += rates
    
        # Initiation of TPT
        rate = r['TPT'][gps_born.index(born)]
        for src, dst in [(Lf, Pf)]:
           m[dst, src] += rate
        
        rate = r['TPT'][gps_born.index(born)]
        for src, dst in [(Ls, Ps)]:
           m[dst, src] += rate
           
        # Case-finding       
        rate = r['ACF'][gps_born.index(born)]
        for src in [I, I2]:  # Iterate over each source state individually
           m[Tx, src] += rate
    
        rate = r['ACF2'][gps_born.index(born)]
        for src, dst in [(I2, Tx)]:
           m[dst, src] += rate
           
           
    # Transition from recent to long-term migrant status (over 5-year period)
    sources = s['mig_recent']
    destinations = s['mig_long']
    # No direct sub2ind needed, use broadcasting to index directly with arrays
    for srcs, dstn in zip(sources, destinations):
        m[dstn, srcs] += 1/5

        
    # ~~~~~~~~~~~~~~~~ LINEAR COMPONENT
    col_sums = np.sum(m, axis=0)  # sum up each column
    mod_diagonal = m.diagonal() + col_sums  # add column sums to the diagonal
    
    # make a sparse diagonal matrix with the modified diagonal elements
    sparse_matrix = (dia_matrix((mod_diagonal, [0]), shape=(i['nstates'], i['nstates']))).toarray()
    
    #M_lin = m - sparse_matrix
    M['lin'] = m - sparse_matrix
    
    # ~~~~~~~~~~~~~~~~ NON LINEAR COMPONENT

    # Create an empty matrix m with the same shape as i.nstates
    # Define the sizes and parameters

    # make matrix m with zeros
    m = np.zeros((i['nstates'], i['nstates']))

    # iterate over gps_born
    for born in gps_born: # iterate over where they are born e.g. dom and for
        # find indices of specified states in s that intersect with 'born'
        # Find the indices of susceptible states that intersect with 'born'
        susinds = np.intersect1d([s[state] for state in ['U', 'Lf', 'Ls', 'Rlo', 'Rhi', 'R']], s[born])
        # calculate intersection between two sets: the states in s, and where they're born
        # born in s
        # initialise those elements in the matrix with a 1
        # use indices in i as the row index and susinds as the column indices
        m[i[('Lf', born)], susinds] = 1 

    # latent and recovered states 
    L_and_R = [s[state] for state in ['Lf', 'Ls', 'Rlo', 'Rhi', 'R']]
    # multiply these by (1 - imm), take all rows only specified columns
    m[:, L_and_R] *= (1 - p['imm'])
    # column sums
    col_sums = np.sum(m, axis=0)
    # Adjust the diagonal elements to subtract the column sums
    diagonal = np.diag(m).copy()
    diagonal += col_sums

    # same as linear
    sparse_diagonal = (dia_matrix((diagonal, [0]), shape=(i['nstates'], i['nstates']))).toarray()
    m_sparse = csr_matrix(m).toarray()

    #M_nlin = m_sparse - sparse_diagonal
    M['nlin'] = m_sparse - sparse_diagonal
    
    # ~~~~~~~~~~~~~~~~ force of infection// lambda
  
    m = np.zeros(i['nstates'])
    # Set values in m at indices specified for infectious states to beta val
    m[s['everyI']] = r['beta'] # WILL BE BETA ONCE OTHER SCRIPTS ARE COMPLETE

    # Create a sparse diagonal matrix 'M.lam' from 'm'
    #M_lam = csr_matrix(m).toarray()
    M['lam'] = csr_matrix(m).toarray()
    
    # ~~~~~~~~~~~~~~~~ mortality

    m = np.zeros((i['nstates'], 2)) # mortality in each state

    # first column of m to life expectancy
    m[:, 0] = 1/83

    # second column's infectious states have tb mort
    m[s['everyI'], 1] = r['muTB']

    # Create a sparse matrix 'M.mort' from 'm'
    #M_mort = (csr_matrix(m)).toarray()
    M['mort'] = (csr_matrix(m)).toarray()


    return M
    








# file_path = 'matrix.txt'

# # Save the matrix to the text file
# np.savetxt(file_path, m, fmt='%.6f', delimiter='\t')

# print(f"Matrix saved to {file_path}")


# # Get the indices (coordinates) of non-zero elements in 'm'
# non_zero_indices = np.transpose(np.nonzero(result_sparse))

# # Display the non-zero positions and their coordinates
# for coord in non_zero_indices:
#     print(f"Coordinate: {tuple(coord)}, Value: {result_sparse[coord[0], coord[1]]}")







# # Get the indices (coordinates) of non-zero elements in 'm'
# non_zero_indices = np.transpose(np.nonzero(M3['lin']))

# # Display the non-zero positions and their coordinates
# for coord in non_zero_indices:
#     incremented_coord = tuple(coord + 1)
#     print(f"{incremented_coord}, Value: {M3['lin'][coord[0], coord[1]]}")
