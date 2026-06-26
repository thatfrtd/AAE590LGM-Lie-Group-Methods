import numpy as np
from scipy.stats.distributions import chi2

def einsum(indices, func):
    esum = 0
    for k in indices:
        esum = esum + func(k)

    return esum

def sqx2inv(eps, dof):
    val = np.sqrt(chi2.ppf(eps, df = dof))
    
    return val