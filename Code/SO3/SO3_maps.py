import numpy as np

def so3_wedge(theta):
    """
    Maps from so3 real number parametrization to so3 Lie algebra
    
    :param theta: angle in R
    """
    Omega = np.block([[0, -theta[3], theta[2]], 
                      [theta[3],  0, -theta[1]],
                      [-theta[2], theta[1], 0]]).reshape((3, 3))

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

    R = np.array([[np.cos(theta), -np.sin(theta)], 
                  [np.sin(theta),  np.cos(theta)]]).reshape((2, 2))

    return R

def so3_log(R):
   """
   Maps from SO3 Lie group to the parametrization of its Lie algebra so3
   
   :param R: 3D rotation matrix
   """
   theta = np.arctan2(R[1, 0], R[0, 0])

   return theta
