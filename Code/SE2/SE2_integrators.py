import numpy as np
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_wedge

def se2_euler_integrator(X_k, xi_k, dt):
    """
    Integrate using first order Euler

    :param X_k: group element at t_k
    :param xi_k: twist at t_k
    :param dt: timestep
    """
    X_kp1 = se2_compose(X_k, np.eye(3) + dt * se2_wedge(xi_k)) 

    return X_kp1

def se2_Lie_group_integrator(X_k, xi_k, dt):
    """
    Integrate using exponential map (assuming twist is constant accross timestep)

    :param X_k: group element at t_k
    :param xi_k: twist at t_k
    :param dt: timestep
    """
    X_kp1 = se2_compose(X_k, se2_exp(dt * xi_k)) 

    return X_kp1