import numpy as np
from Code.SO3.SO3_maps import so3_wedge, so3_vee, so3_exp, so3_log

def se3_compose(X1, X2):
    """
    Composition of elements of SE(3)
    
    :param X1: element of SE(3)
    :param X2: element of SE(3)
    """
    X1_t = X1[0:2, 2].reshape((2, 1))
    X1_R = X1[0:2, 0:2]
    
    X2_t = X2[0:2, 2].reshape((2, 1))
    X2_R = X2[0:2, 0:2]

    X12_R = X1_R @ X2_R
    X12_t = X1_R @ X2_t + X1_t
    
    X12 = np.block([[X12_R, X12_t], 
                    [np.zeros((1, 2)), 1]])

    return X12

def se3_inverse(X):
    """
    Inverse of SE(3) element
    
    :param X: element of SE(3)
    """
    X_t = X[0:2, 2].reshape((2, 1))
    X_R = X[0:2, 0:2]

    X_inv_t = -X_R.T @ X_t
    X_inv_R = X_R.T

    X_inv = np.block([[X_inv_R, X_inv_t],
                      [np.zeros((1, 2)), 1]])

    return X_inv

def se3_wedge(xi):
    """
    Maps from se3 real number parametrization to se3 Lie algebra
    
    :param xi: twist vector (v_x, v_y, omega) in R3
    """
    xiwedge = np.block([[so3_wedge(xi[2]), xi[0:2].reshape((2, 1))], 
                        [np.zeros((1, 2)),       0]])

    return xiwedge

def se3_vee(xiwedge):
    """
    Maps from se3 Lie algebra to its real number parametrization 
    
    :param xiwedge: 3x3 se3 Lie algebra element
    """
    
    xi = np.block([xiwedge[[0, 1], [2, 2]], so3_vee(xiwedge[0:2, 0:2])]);

    return xi

def se3_exp(xi):
    """
    Maps from parametrization of se3 Lie algebra to SE3 Lie group
    
    :param xi: twist vector (v_x, v_y, omega) in R3
    """
    # Extract elements
    v = xi[0:2]
    omega = xi[2]

    # SO3 exponential map
    R = so3_exp(omega)

    # Translation-orientation coupling
    if abs(omega) > 1e-8:
        V = np.sin(omega) / omega * np.eye(2) + (1 - np.cos(omega)) / omega * np.array([[0, -1], [1, 0]])
    else:
        V = ((1 - omega ** 2 / 6 + omega ** 4 / 120) * np.eye(2) 
           + (omega / 2 - omega ** 3 / 24) * np.array([[0, -1], [1, 0]]))

    X = np.block([[R,    V @ v.reshape((2, 1))], 
                  [np.zeros((1, 2)),     1]])
    
    return X

def se3_log(xiwedge):
   """
   Maps from SE3 Lie group to the parametrization of its Lie algebra se3
   
   :param xiwedge: 3x3 se3 Lie algebra element
   """
   omega = so3_log(xiwedge[0:2, 0:2])

   if abs(omega) > 1e-8:
       V_inv = omega / 2 * (np.cos(omega / 2) / np.sin(omega / 2) * np.eye(2) + np.array([[0, 1], [-1, 0]]))
   else:
       V_inv = 1 / 2 * ((2 - omega ** 2 / 6 - omega ** 4 / 360 - 2 * omega ** 6 / 30240 + omega ** 8 / 604800) * np.eye(2) + omega * np.array([[0, -1], [1, 0]]))

   v = V_inv @ xiwedge[0:2, 2]

   xi = np.vstack((v.reshape((2, 1)), np.reshape(omega, (1,1))))

   return xi

def se3_Ad(X):
    """
    Maps from se3 to se3 defined by (Ad_X xi)**wedge = X xi**wedge X**-1
    Shifts tangent spaces

    :param X: SE3 Lie group element
    """
    Ad = np.copy(X)
    Ad[0:2, 2] = -np.array([[0, -1], [1, 0]]) @ X[0:2, 2]

    return Ad

def se3_ad(xi):
    """
    Derivative of Adjoint at identity
    
    :param xi: se(2) Lie algebra element parametrization
    """
    ad = se3_wedge(xi)
    ad[0:2, 2] = -np.array([[0, -1], [1, 0]]) @ ad[0:2, 2]
    
    return ad