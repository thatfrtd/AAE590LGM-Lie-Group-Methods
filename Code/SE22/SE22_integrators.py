import numpy as np
from Code.SE22.SE22_maps import se22_compose, se22_exp, se22_wedge

def se22_Lie_group_integrator(X_k, xi_k, dt):
    """
    Integrate using exponential map (assuming twist is constant accross timestep)

    :param X_k: group element at t_k
    :param xi_k: twist at t_k
    :param dt: timestep
    """
    Gamma_dt = np.block([[np.eye(2), np.zeros((2, 2))],
                         [0, 0, 1, dt],
                         [0, 0, 0, 1]]) # Coupling matrix
    
    a = se22_exp(dt * xi_k) @ Gamma_dt
    X_kp1 = se22_compose(X_k, se22_exp(dt * xi_k) @ Gamma_dt) 

    return X_kp1