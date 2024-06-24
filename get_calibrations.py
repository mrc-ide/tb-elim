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

def nobj(x):
    return -obj(x)  # Return the negated objective function value for minimization

obj = lambda x: get_objective(x, ref, prm, gps_born, likelihood)[0]
obj2 = lambda x: get_objective(x, ref, prm, gps_born,likelihood)

# Flattening the bounds for Latin Hypercube Sampling
flat_bounds = []
for key, bounds in prm['bounds'].items():
    if isinstance(bounds[0], list):  # Handling nested lists for bounds
        for bound in bounds:
            flat_bounds.append(bound)
    else:
        flat_bounds.append(bounds)

# Separating lower and upper bounds
lower_bounds, upper_bounds = zip(*flat_bounds)

# Generating samples with Latin Hypercube Sampling
# Generating samples with Latin Hypercube Sampling
nsam = 1000  # Reduced to 10 for debugging
print("Total number of samples:", nsam)  # Debugging
with tqdm(total=nsam, desc="Evaluating Objective") as pbar:
    sampler = LatinHypercube(d=len(flat_bounds))
    sample = sampler.random(n=nsam)
    scaled_sample = qmc.scale(sample, lower_bounds, upper_bounds)

    # Debugging
    print("Length of scaled_sample:", len(scaled_sample))

    # Evaluating the objective function for each sample
    outs = np.zeros(nsam)
    for i, x in enumerate(scaled_sample[:nsam]):  # Limiting the loop to iterate over only 10 samples
        print("Current index:", i)  # Debugging
        try:
            outs[i] = obj(x)
        except Exception as e:
            print(f"Error at index {i}: {e}")
        pbar.update(1)        
# Sorting samples by the objective function values to find the best starting point
ord_indices = np.argsort(-outs)  # Negate for descending order
best_sample = scaled_sample[ord_indices][0]

# Optimization starting from the best sample
with tqdm(total=len(best_sample), desc="Optimizing") as pbar:
    res = minimize(nobj, best_sample, method='Nelder-Mead', callback=lambda x: pbar.update(1))

# Print the optimized parameters and the objective function value
print("Optimized Parameters:", res.x)
print("Objective Function Value:", -res.fun)

# Optimizing starting from the best sample
x0 = res.x  # Use the best sample from previous optimization
res_xord_1 = res.x  # Store the first element of xord (xord[0]) to use as x0

x0_length = len(x0)  # Obtain the length of x0
result = MCMC_adaptive(obj, x0, 1e4, 1, np.eye(x0_length))
xsto, outsto = result[0], result[1]

# Perform MCMC
#xsto, outsto = MCMC_adaptive(obj, x0, 1e2, 1, np.eye(len(x0)))

# Find the index of the maximum value in outsto
inds = np.where(outsto == np.max(outsto))[0]
x0 = xsto[inds[0], :]

# Call obj with x0
out, aux = obj2(x0)
sfin = aux['soln'][-1, :]
migr_indices = np.intersect1d(s['migr'], [s['Lf'], s['Ls']])
migr_sum = np.sum(sfin[migr_indices]) / np.sum(sfin[s['migr']])


cov0 = np.cov(xsto.T)
result = MCMC_adaptive(obj, x0, 1e4, 1, cov0=cov0)
xsto, outsto = result[0], result[1]



niter = [1, 1, 1, 5] * int(1e3)

for ii in range(len(niter)):
    inds = np.where(outsto == np.max(outsto))[0]
    x0 = xsto[inds[0], :]
    cov0 = np.cov(xsto.T)
    result = MCMC_adaptive(obj, x0, niter[ii], 1, cov0=cov0)
    
xsto, outsto = result[0], result[1]


# Save calibration_res
np.savez('calibration_res.npz', xsto=xsto, outsto=outsto)

# Additional part
x2 = xsto[-1, :]
cov0 = np.cov(xsto.T)
results2 = MCMC_adaptive(obj, x2, 1e2, 1, cov0=cov0)
xsto2, outsto2 = results2[0], results2[1]

nx = 200
ix0 = len(xsto) // 2
if len(xsto) < 2 * nx:
    dx = 1  # Set a default value for dx if len(xsto) is too small
else:
    dx = len(xsto) // (2 * nx)

xs = xsto[ix0::dx]

for ii in range(len(xs)):
    if (ii + 1) % (len(xs) // 24) == 0:
        print(f"{(ii + 1) / (len(xs) // 24)} ", end="")
    obj2res = obj2(xs[ii])
    out, aux = obj2res[0], obj2res[1]
    sim = np.array([aux['incd'], aux['mort'], aux['p_migrTB'], aux['p_migrpopn'], aux['p_LTBI']])

# Save calibration_res again
np.savez('calibration_res.npz', xsto=xsto, outsto=outsto, xsto2=xsto2, outsto2=outsto2)

for ii in range(len(xs)):
    if (ii + 1) % (len(xs) // 24) == 0:
        print(f"{(ii + 1) / (len(xs) // 24)} ", end="")
    obj2res = obj2(xs[ii])
    out, aux = obj2res[0], obj2res[1]
    # Ensure each element of sim is a scalar value
    incd = aux['incd'][0] if isinstance(aux['incd'], (list, np.ndarray)) else aux['incd']
    mort = aux['mort'][0] if isinstance(aux['mort'], (list, np.ndarray)) else aux['mort']
    p_migrTB = aux['p_migrTB'][0] if isinstance(aux['p_migrTB'], (list, np.ndarray)) else aux['p_migrTB']
    p_migrpopn = aux['p_migrpopn'][0] if isinstance(aux['p_migrpopn'], (list, np.ndarray)) else aux['p_migrpopn']
    p_LTBI = aux['p_LTBI'][0] if isinstance(aux['p_LTBI'], (list, np.ndarray)) else aux['p_LTBI']
    sim = np.array([incd, mort, p_migrTB, p_migrpopn, p_LTBI])
