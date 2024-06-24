# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 13:33:40 2023

@author: ISA
"""

import numpy as np
from goveqs_basis import goveqs_basis3
#from allocate import allocate_parameters
#from setup_model import prm,gps_born
#from make_model import make_model


# times = [2022, 2025]
# t=np.arange(2022, 2037)

# def goveqs_scaleup(t, insert, i, s, M0, M1, p0, p1, times, agg, sel, r):
#     #scale = min(max((t - times[0]) / (times[1] - times[0]), 0), 1)
#     scale = np.minimum(1.0, np.maximum(0.0, (t - times[0]) / (times[1] - times[0])))

    
#     Ms = M1.copy()
#     Ms['lin'] = M0['lin'] + scale * (M1['lin'] - M0['lin'])
    
#     ps = p1.copy()
#     ps['migrTPT'] = p0['migrTPT'] + scale * (p1['migrTPT'] - p0['migrTPT'])
    
#     out = goveqs_basis2(t, insert, i, s, Ms, agg, sel, r, ps)
    
#     return out







# def goveqs_scaleup(t, insert, M0, M1, M2, times, i, s, p, r, prm, sel, agg):
#     if not isinstance(times, list) or len(times) != 2:
#         raise ValueError("Times should be a list containing two tuples [(start1, end1), (start2, end2)]")
    
#     for interval in times:
#         if not isinstance(interval, tuple) or len(interval) != 2:
#             raise ValueError("Each interval in times should be a tuple (start, end)")
    
#     scale = np.maximum((t - np.array([start for start, _ in times])) / np.maximum(np.array([end - start for start, end in times]), 1e-6), 0)
#     scale[0] = min(scale[0], 1)
    
#     # Creating a new dictionary to hold Mt
#     Mt = M1.copy()  # Assuming M1, M0, and M2 are dictionaries with a 'lin' key
#     # Updating 'lin' in Mt based on calculated scales
#     Mt['lin'] = M0['lin'] + scale[0] * (M1['lin'] - M0['lin']) + scale[1] * (M2['lin'] - M0['lin'])
    
#     # Compute output using another function, assuming goveqs_basis3 is also translated to Python
#     out = goveqs_basis3(t, insert, i, s, Mt, agg, sel, r, p)
#     return out






def goveqs_scaleup(t, insert, M0, M1, M2, times, i, s, p, r, prm, sel, agg):
    # Calculate scale factors based on time t and given intervals in times
    # Ensuring no division by zero in the denominator
    times = np.array(times).reshape((2, 2))
    scale = np.maximum((t - times[:, 0]) / np.maximum(times[:, 1] - times[:, 0], 1e-6), 0)
    scale[0] = min(scale[0], 1)
    
    # Creating a new dictionary to hold Mt
    Mt = M1.copy()  # Assuming M1, M0, and M2 are dictionaries with a 'lin' key
    # Updating 'lin' in Mt based on calculated scales
    Mt['lin'] = M0['lin'] + scale[0] * (M1['lin'] - M0['lin']) + scale[1] * (M2['lin'] - M0['lin'])
    
    # Compute output using another function, assuming goveqs_basis3 is also translated to Python
    out = goveqs_basis3(t, insert, i, s, Mt, agg, sel, r, p)
    return out























# def goveqs_scaleup(t, insert, i, s, M0, M1, p0, p1, times, agg, sel, r):
#     scale = np.minimum(np.maximum((t - times[0]) / (times[1] - times[0]), 0), 1)

#     Ms = M1.copy()
#     Ms['lin'] = M0['lin'] + scale * (M1['lin'] - M0['lin'])

#     ps = p1.copy()
#     ps['migrTPT'] = p0['migrTPT'] + scale * (p1['migrTPT'] - p0['migrTPT'])

#     out = goveqs_basis2(t, insert, i, s, Ms, agg, sel, r, ps)

#     return out
