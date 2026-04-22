import numpy as np
import scipy as sp

def NIS(innovation, S_k_inv):
    '''
    Normalized Innovation Squared

    Nice because it only requires the filter's own quantities so no ground truth is needed
    '''
    nu_k = innovation.T @ S_k_inv @ innovation
    return nu_k

def NEES(xhat_error, P_k):
    '''
    Normalized Estimation Error Squared
    '''
    eps_k = xhat_error.T @ np.linalg.inv(P_k) @ xhat_error
    return eps_k

def averaged_chi_squared_consistency_bands(M, n):
    band_bounds = (sp.stats.chi2.ppf(0.025, M * n) / M, sp.stats.chi2.ppf(0.975, M * n) / M)
    return band_bounds