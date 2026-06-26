import numpy as np

def so3_wedge(theta):
    """
    Maps from so3 real number parametrization to so3 Lie algebra
    
    :param theta: angle in R
    """
    Omega = np.block([[0, -theta[2], theta[1]], 
                      [theta[2],  0, -theta[0]],
                      [-theta[1], theta[0], 0]]).reshape((3, 3))

    return Omega

def so3_vee(Omega):
    """
    Maps from so3 Lie algebra to its real number parametrization 
    
    :param Omega: 3x3 skew symmetric matrix (in so3 Lie algebra)
    """
    
    theta = np.array((Omega[2, 1], Omega[0, 2], Omega[1, 0])).reshape((3, 1))

    return theta

def so3_exp(theta):
    """
    Maps from parametrization of so3 Lie algebra to SO3 Lie group
    
    :param theta: angle in R
    """

    theta_mag = np.linalg.norm(theta)
    lambda_vec = theta / theta_mag
    lambda_wedge = so3_wedge(lambda_vec)
    R = np.eye(3) + lambda_wedge * np.sin(theta_mag) + lambda_wedge @ lambda_wedge * (1 - np.cos(theta_mag))

    return R

def so3_log(R):
    """
    Maps from SO3 Lie group to the parametrization of its Lie algebra so3
   
    :param R: 3D rotation matrix
    """

    theta_mag = np.acos((np.linalg.trace(R) - 1) / 2)
    if np.abs(theta_mag) >= 1e-8:
        theta = theta_mag * so3_vee(R - R.T) / (2 * np.sin(theta_mag))
        '''theta = np.array([[R[2, 1] - R[1, 2]],
                          [R[0, 2] - R[2, 0]],
                          [R[1, 0] - R[0, 1]]]) / (2 * np.sin(theta_mag)) * theta_mag'''
    else:
        theta = 0.5 * (1 + 1 / 6 * theta_mag ^ 2 + 7 / 360 * theta_mag ^ 4) * so3_vee(R - R.T)

    return theta

def so3_ljac(theta):
    theta_mag = np.linalg.norm(theta)

    if theta_mag > 1e-8:
        lambda_vec = theta / theta_mag
        V = np.eye(3) + (1 - np.cos(theta_mag)) / theta_mag * so3_wedge(lambda_vec) + (theta_mag - np.sin(theta_mag)) / theta_mag * so3_wedge(lambda_vec) @ so3_wedge(lambda_vec)
    else:
        V = np.eye(3) + 1 / 2 * so3_wedge(theta) + 1 / 6 * so3_wedge(theta) * so3_wedge(theta)

    return V

def so3_inv_ljac(theta):
    theta_mag = np.linalg.norm(theta)

    if theta_mag > 1e-8:
        lambda_vec = theta / theta_mag
        V_inv = np.eye(3) - 1 / 2 * so3_wedge(theta) + (1 - theta_mag * (1 + np.cos(theta_mag)) / (2 * np.sin(theta_mag))) * so3_wedge(lambda_vec) @ so3_wedge(lambda_vec)
    else:
        V_inv = np.eye(3) - 1 / 2 * so3_wedge(theta) + 1 / 12 * so3_wedge(theta) * so3_wedge(theta)

    return V_inv