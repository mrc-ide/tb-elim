# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 13:29:40 2023

@author: ISA
"""

#from setup_model import p,r,xi
import numpy as np


def allocate_parameters(x, p, r, xi):
    r['beta'] = x[xi['beta'][0]]
    p['betadec'] = x[xi['betadec'][0]]
    r['gamma_2015'] = x[xi['gamma'][0]]
    r['gamma_2020'] = x[xi['gamma'][1]]
    r['progression'] = r['progression0'] * np.array([1, x[xi['p_relrate'][0]], 1])
    r['reactivation'] = r['reactivation0'] * np.array([1, x[xi['p_relrate'][0]], 1])
    r['migr'] = x[xi['r_migr'][0]]
    p['LTBI_in_migr'] = x[xi['p_LTBI_in_migr'][0]]
    return p, r
