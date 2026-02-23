import numpy as np
from Code.SE3.SE3_maps import se3_compose, se3_exp, se3_wedge

def se3_euler_integrator(X_k, xi_k, dt):
    """
    Integrate using first order Euler

    :param X_k: group element at t_k
    :param xi_k: twist at t_k
    :param dt: timestep
    """
    X_kp1 = se3_compose(X_k, np.eye(3) + dt * se3_wedge(xi_k)) 

    return X_kp1

def se3_Lie_group_integrator(X_k, xi_k, dt):
    """
    Integrate using exponential map (assuming twist is constant accross timestep)

    :param X_k: group element at t_k
    :param xi_k: twist at t_k
    :param dt: timestep
    """
    X_kp1 = se3_compose(X_k, se3_exp(dt * xi_k)) 

    return X_kp1