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
