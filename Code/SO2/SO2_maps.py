import numpy as np

def so2_wedge(theta):
    """
    Maps from so2 real number parametrization to so2 Lie algebra
    
    :param theta: angle in R
    """
    Omega = np.block([[0, -theta], 
                      [theta,  0]]).reshape((2, 2))

    return Omega

def so2_vee(Omega):
    """
    Maps from so2 Lie algebra to its real number parametrization 
    
    :param Omega: 2x2 skew symmetric matrix (in so2 Lie algebra)
    """
    
    theta = Omega[1, 0]

    return theta

def so2_exp(theta):
    """
    Maps from parametrization of so2 Lie algebra to SO2 Lie group
    
    :param theta: angle in R
    """

    R = np.array([[np.cos(theta), -np.sin(theta)], 
                  [np.sin(theta),  np.cos(theta)]]).reshape((2, 2))

    return R

def so2_log(R):
   """
   Maps from SO2 Lie group to the parametrization of its Lie algebra so2
   
   :param R: 2D rotation matrix
   """
   theta = np.arctan2(R[1, 0], R[0, 0])

   return theta

def so2_ljac(theta):
    """

    :param theta: lie algebra element
    """

    if abs(theta) > 1e-8:
        J_l = np.sin(theta) / theta * np.eye(2) + (1 - np.cos(theta)) / theta * np.array([[0, -1], [1, 0]])
    else:
        J_l = ((1 - theta ** 2 / 6 + theta ** 4 / 120) * np.eye(2) 
           + (theta / 2 - theta ** 3 / 24) * np.array([[0, -1], [1, 0]]))
    
    return J_l

def so2_inv_ljac(theta):
    """
    Inverse left Jacobian
    
    :param theta: lie algebra element
    """

    if abs(theta) > 1e-8:
        J_linv = theta / 2 * (np.cos(theta / 2) / np.sin(theta / 2) * np.eye(2) + np.array([[0, 1], [-1, 0]]))
    else:
        J_linv = 1 / 2 * ((2 - theta ** 2 / 6 - theta ** 4 / 360 - 2 * theta ** 6 / 30240 + theta ** 8 / 604800) * np.eye(2) + theta * np.array([[0, -1], [1, 0]]))

    return J_linv