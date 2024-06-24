# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 11:05:31 2023

@author: ia19
"""

#~~~~~~~~~~~PYTHON INDEXING STARTS AT 0 (ZERO)

#combining states and born will be used as groups in function argument
# states = ['U', 'Lf', 'Ls', 'Pf', 'Ps', 'I', 'I2', 'Tx', 'Rlo', 'Rhi', 'R'] # states
# gps_born = ['dom', 'for'] # where they are born
# groups = ([states, gps_born])
#initialise lists


def get_addresses(groups, i=None, s=None, d=None, lim=0):
    if i is None:
        i = {}
    if s is None:
        s = {}
    if d is None:
        d = {} #creates empty dictionaries
    
    # here we're iterating over groups, which will be a list of lists
    # initialising empty lists in s for each element in groups
    # setting up data to organise the states by group
    for ig in range(len(groups)): # iterates over the indices of groups. len(groups) gives number of groups in the groups list
        gp = groups[ig] # extracts the n'th element from groups, call this subgroup gp
        for ig2 in range(len(gp)): # inner loop iterates over the indices of gp
            if gp[ig2] not in s: # check if the n-th element of gp is not already a key in s.
                s[gp[ig2]] = [] # if the element from gp isnt in s, add it as a key to s with an empty list


    if len(groups) == 1: #check if there's only one element in groups
        gp1 = groups[0] # if there's only one group, call that group gp1 
        for ig1 in range(len(gp1)): # iterate over elements of gp1
            lim += 1 # for each element in gp1, increment lim
            i[(gp1[ig1])] = lim - 1 # based off lim, assign an index to that element in i
            s[gp1[ig1]].append(lim - 1) # add the current index to the list attached to the state
            d[lim - 1] = [gp1[ig1]] # state name is stored in d with the lim index as the key

    if len(groups) == 2: # check if there's two elements in groups
        gp1 = groups[0] # if there's two groups, call them gp1 and gp2
        gp2 = groups[1]
        for ig1 in range(len(gp1)): # iterate over both groups gp1 and gp2
            for ig2 in range(len(gp2)):
                lim += 1 # for each element in these groups, increment lim
                i[(gp1[ig1], gp2[ig2])] = lim - 1 # assigns a number index to the state key in the dict
                s[gp1[ig1]].append(lim - 1) # adds index to the state in the first group
                s[gp2[ig2]].append(lim - 1) # # adds index to the state in the second group
                d[lim - 1] = [gp1[ig1], ' ', gp2[ig2]]

    i['nstates'] = lim # set total no. of states in i to final val for lim
    return i, s, d, lim
