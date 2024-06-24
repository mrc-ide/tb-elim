# -*- coding: utf-8 -*-
"""

@author: ISA
"""

import numpy as np
from obj import get_objective
from allocate import allocate_parameters
from make_model import make_model
from gov_eqs_scaleup import goveqs_scaleup
from setup_model import i,s,r,p,agg,sel,ref,xi,prm,gps_born,likelihood
from scipy.integrate import odeint
from goveqs_basis import goveqs_basis2
import copy
import matplotlib.pyplot as plt


#2017-18 testing - 14,770
# same year shouldve had tpt - 1863
# proportion =  0.12613405551794177

#2018-19 testing - 17,104
# same year shouldve had tpt - 1970
# proportion =  0.115177736202058

#2019-20 testing - 22221
# same year shouldve had tpt - 2049
# proportion =  0.09221007155393547


#2016-17 tpt access should - > 1208

#2015-16 tpt access should - > 317


#tpt needed proportion over 3 years of observed data 0.11117395442464507

#scale up 2015-16 to 2016-17 is 3.8107
#2017-18 to 2018-19 is 1.0574
# 2018-19 to 2019-20 is 1.0401



xsto = np.load('xsto.npy')

xsto = [xsto[:, i][::100] for i in range(5)]
xsto = np.column_stack(xsto)
obj = lambda x: get_objective(x, ref, prm, gps_born,likelihood)
t1 = np.arange(2015, 2021)

for params in xsto:
    # Stage I
    t0 = np.arange(2015)
    soln_stage1 = obj(params)
    final_state = soln_stage1[1]['soln'][-1]

    # Stage II
    t1 = np.arange(2015, 2021)  # Assuming 2023 as the end year, modify as needed
    current_state = final_state
    for year in t1:
        params_with_TPT = update_parameters_with_TPT(year, params)
        soln_stage2 = obj([year, year+1], current_state, params_with_TPT)
        current_state = soln_stage2[-1]
    
    # Store the final state for this parameter set after both stages
    results.append(current_state)








# Assuming you've already written the Python equivalents of your functions
# like get_objective2, allocate_parameters, make_model, goveqs_basis2, and goveqs_scaleup

# Load data
# Assuming calibration_res.mat and Model_setup.mat have been converted to Python-friendly formats

# Function definition

# Selection of samples
ix0 = round(len(xsto) / 2)
dx = round(len(xsto) / 2 / 150)
xs = xsto[ix0::dx, :]
incsto = np.zeros((14, 167, 5))
mrtsto = np.zeros((14, 167, 5))

mk = round(len(xs) / 25)
for ii in range(len(xs)):
    if ii % mk == 0:
        print(f"{ii/mk:.5g}", end=" ")

    xx = xs[ii, :]
    out, aux = obj(xx)
    init = aux['soln'][-1, :]

    p0, r0 = allocate_parameters(xx, p, r, xi)
    M0 = make_model(p0, r0, i, s, gps_born)
    
    p1, r1 = copy.deepcopy(p0), copy.deepcopy(r0)
    p1['migrTPT'] = 1
    M1 = make_model(p1, r1, i, s, gps_born)
    
    p2, r2 = copy.deepcopy(p0), copy.deepcopy(r0)
    p2['migrTPT'] = 1
    r2['TPT'] = [0.69 * x for x in [0, 1]]
    M2 = make_model(p2, r2, i, s, gps_born)
    
    # Models list
    models = [M0, M2]
    

    for mi, model in enumerate(models):
        def geq(t, in_):
            return goveqs_scaleup(t, in_, i, s, M0, models[mi], p0, p1, [2022, 2025] ,agg, sel, r).flatten()
        soln0 = odeint(geq, init, np.arange(2022, 2037), tfirst=True)

        sdiff = np.diff(soln0, axis=0)
        incsto[:, ii, mi] = sdiff[:,i["aux"]["inc"][0]] * 1e5
        mrtsto[:, ii, mi] = sdiff[:,-1] * 1e5

print()





incmat = np.percentile(incsto, [2.5, 50, 97.5], axis=1).transpose(1, 0, 2)
mrtmat = np.percentile(mrtsto, [2.5, 50, 97.5], axis=1).transpose(1, 0, 2)
incmat = np.transpose(incmat, (1, 0, 2))


xx = np.arange(2022, 2036)
cols = [
    (0.8, 0.2, 0.2),  # Light Red
    (0.2, 0.8, 0.2),  # Light Green
    (0.2, 0.2, 0.8),  # Light Blue
    (0.7, 0.3, 0.7),  # Purple
    (0.9, 0.6, 0.0)   # Orange
]
lw = 2  # Line width
fs = 10  # Font size

fig, ax = plt.subplots()
lg = []  # To store the legends

for ii in range(incmat.shape[2]):
    plt_data = incmat[:, :, ii]
    line, = ax.plot(xx, plt_data[1, :], color=cols[ii], linewidth=lw)
    ax.fill_between(xx, plt_data[0, :], plt_data[2, :], color=cols[ii], alpha=0.1)
    lg.append(line)

ax.axhline(y=1.05, color='k', linestyle='--')
ax.set_ylim(bottom=0)
ax.set_xlim(2022, 2035)
ax.set_title('Incidence')
ax.set_ylabel('Rate per 100,000 population')
ax.set_xlabel('Year')
ax.legend(lg, [
    'Baseline', 'ACF', 'ACF + domestic TPT', 'ACF + domestic AND migrant TPT',
    '+ Monthly followup post TPT or Tx', 'Elimination target'
], loc='SouthWest')
ax.tick_params(axis='both', which='major', labelsize=fs)

plt.show()



xx = np.arange(2022, 2036)
lw = 2.0
fs = 12  # Font size
cols = [
    [0.0, 0.0, 0.0],
    [0.7, 0.0, 0.0],
    [0.0, 0.7, 0.0],
    [0.0, 0.0, 0.7],
    [0.7, 0.7, 0.0]
]

fig, ax = plt.subplots()
lines = []

for ii in range(incmat.shape[2]):
    plt_data = incmat[:, :, ii]
    print(f"Plotting scenario {ii + 1}")  # Check iteration
    print(plt_data[1, :])  # Print the data being plotted for verification

    line, = ax.plot(xx, plt_data[1, :], color=cols[ii], linewidth=lw, label=f"Scenario {ii + 1}")
    ax.fill_between(xx, plt_data[0, :], plt_data[2, :], color=cols[ii], alpha=0.1)
    lines.append(line)

ax.axhline(y=1.05, linestyle='--', color='k')
ax.set_ylim(bottom=0)
ax.set_xlim([2022, 2035])
ax.set_title('Incidence')
ax.set_ylabel('Rate per 100,000 population')
ax.set_xlabel('Year')
ax.legend(loc='SouthWest', fontsize=fs)

plt.show()