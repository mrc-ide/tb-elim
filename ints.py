# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 14:21:15 2024

@author: ia19
"""

import numpy as np
from scipy.integrate import odeint
import numpy as np
from scipy.integrate import odeint
from make_model import make_model
from allocate import allocate_parameters
from setup_model import likelihood_function
from goveqs_basis import goveqs_basis3
from gov_eqs_scaleup import goveqs_scaleup
import numpy as np
from setup_model import s, ref, prm, gps_born, likelihood
from obj import get_objective
from pyDOE2 import lhs
from mcmc2 import MCMC_adaptive
from tqdm import tqdm
import numpy as np
from scipy.optimize import minimize
from scipy.stats.qmc import LatinHypercube
from scipy.stats import qmc
np.random.seed(42)


# Load necessary data
calibration_res = np.load('calibration_res.npz')
xsto = calibration_res['xsto']
outsto = calibration_res['outsto']

obj = lambda x: get_objective(x, ref, prm, gps_born, likelihood)

ix0 = round(xsto.shape[0] / 2)
dx = round(xsto.shape[0] / (2 * 25))
xs = xsto[ix0::dx]

ix0 = round(xsto.shape[0] / 2)
dx = round(xsto.shape[0] / (2 * 25))
if dx != 0:
    xs = xsto[ix0::dx]
else:
    print("Step size cannot be zero. Unable to perform slicing operation.")
# Sort xs array
sorted_indices = np.argsort(xs[:, 0])  # Assuming time is the first column in xs
xs = xs[sorted_indices]

mk = round(xs.shape[0] / 25)

for ii in range(xs.shape[0]):
    if (ii + 1) % mk == 0:
        print(f"{(ii + 1) / mk} ", end="")

    xx = xs[ii]
    out, aux = obj(xx)

    init = aux['soln'][-1]

    p0, r0 = allocate_parameters(xx, prm['p'], prm['r'], ref['xi'])
    r0['gamma'] = r0['gamma_2020']
    M0 = make_model(p0, r0, ref['i'], ref['s'], gps_born)

    # ---------------------------------------------------------------------
    # --- Model baseline
    
#     geq = lambda t, in_: goveqs_basis2(t, in_, ref['i'], ref['s'], M0, prm['agg'], prm['sel'], r0, p0)
#     t, soln = odeint(geq, [2022, 2031], init)
#     sdiff = np.diff(soln, axis=0)
#     incsto[:, ii, 0] = sdiff[:, ref['i']['aux']['inc'][0]] * 1e5
    
    # ---------------------------------------------------------------------
    # --- Model intervention
    
    p1 = p0.copy()
    r1 = r0.copy()
    r1['TPT'] = 1.4 * np.array([0, 1, 0])
    M1 = make_model(p1, r1, ref['i'], ref['s'], gps_born)

    p2 = p1.copy()
    r2 = r1.copy()
    r2['ACF'] = 0.69 * np.array([1, 1, 1])
    M2 = make_model(p2, r2, ref['i'], ref['s'], gps_born)

    p3 = p2.copy()
    r3 = r2.copy()
    p3['migrTPT'] = 0.25
    M3 = make_model(p3, r3, ref['i'], ref['s'], gps_born)

    p4 = p3.copy()
    r4 = r3.copy()
    r4['TPT'] = 0.2876 * np.array([0, 1, 0])
    M4 = make_model(p4, r4, ref['i'], ref['s'], gps_born)

    models = [M0, M1, M2, M3]
    
    num_models = len(models)
    num_xs = xs.shape[0]
    num_inc = len(ref['i']['aux']['inc'])
    num_mort = len(ref['i']['aux']['mort'])
    num_inc_sources = len(ref['i']['aux']['sources'])
    endsolsto = np.zeros((num_models, len(init)))
    incsto = np.zeros((num_inc, num_xs, num_models))
    mrtsto = np.zeros((num_mort, num_xs, num_models))
    props = np.zeros((num_xs, num_inc_sources, num_models))

    for mi, model in enumerate(models):
        geq = lambda t, in_: goveqs_scaleup(t, in_, ref['i'], ref['s'], M0, models[mi], p0, p1 if mi < 4 else p3, [2022, 2025], prm['agg'], prm['sel'], r0)
        t, soln = odeint(geq, np.arange(2022, 2037), init)
        endsolsto[mi] = soln[-1]
        sdiff = np.diff(soln, axis=0)
        incsto[:, ii, mi] = sdiff[:, ref['i']['aux']['inc'][0]] * 1e5
        mrtsto[:, ii, mi] = sdiff[:, ref['i']['aux']['mort']] * 1e5
        
        # Get proportions from different sources
        vec = sdiff[-1, ref['i']['aux']['incsources']] * 1e5
        props[ii, :, mi] = vec / np.sum(vec)

print()


incmat = np.percentile(incsto, [2.5, 50, 97.5], axis=1).transpose(1, 0, 2)
mrtmat = np.percentile(mrtsto, [2.5, 50, 97.5], axis=1).transpose(1, 0, 2)
