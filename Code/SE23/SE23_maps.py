import numpy as np
from Code.SO3.SO3_maps import so3_wedge, so3_vee, so3_exp, so3_log

def se23_compose(X1, X2):
    """
    Composition of elements of SE2(3)
    
    :param X1: element of SE2(3)
    :param X2: element of SE2(3)
    """
    X1_t = X1[0:3, 4].reshape((3, 1))
    X1_v = X1[0:3, 3].reshape((3, 1))
    X1_R = X1[0:3, 0:3]
    
    X2_t = X2[0:3, 4].reshape((3, 1))
    X2_v = X2[0:3, 3].reshape((3, 1))
    X2_R = X2[0:3, 0:3]

    X12_R = X1_R @ X2_R
    X12_v = X1_R @ X2_v + X1_v
    X12_t = X1_R @ X2_t + X1_t
    
    X12 = np.block([[X12_R, X12_v, X12_t], 
                    [np.zeros((3, 3)), np.eye(3)]])

    return X12

def se23_inverse(X):
    """
    Inverse of SE2(3) element
    
    :param X: element of SE2(3)
    """
    X_t = X[0:3, 4].reshape((3, 1))
    X_v = X[0:3, 3].reshape((3, 1))
    X_R = X[0:3, 0:3]

    X_inv_t = -X_R.T @ X_t
    X_inv_v = -X_R.T @ X_v
    X_inv_R = X_R.T

    X_inv = np.block([[X_inv_R, X_inv_v, X_inv_t],
                      [np.zeros((3, 3)), np.eye(3)]])

    return X_inv

def se23_wedge(xi):
    """
    Maps from se23 real number parametrization to se23 Lie algebra
    
    :param xi: twist vector (a_x, a_y, b_x, b_y, omega) in R3
    """
    xiwedge = np.block([[so3_wedge(xi[4]), xi[0:2].reshape((2, 1)), xi[2:4].reshape((2, 1))], 
                        [np.zeros((2, 3)), np.zeros((2, 2))]])

    return xiwedge

def se23_vee(xiwedge):
    """
    Maps from se23 Lie algebra to its real number parametrization 
    
    :param xiwedge: 3x3 se23 Lie algebra element
    """
    
    xi = np.block([xiwedge[[0, 1], [2, 2]], xiwedge[[0, 1], [3, 3]], so3_vee(xiwedge[0:3, 0:3])]);

    return xi

def se23_exp(xi):
    """
    Maps from parametrization of se23 Lie algebra to se23 Lie group
    
    :param xi: twist vector (a_x, a_y, b_x, b_y, omega) in R3
    """
    # Extract elements
    a = xi[0:3]
    b = xi[3:6]
    omega = xi[6:9]

    # SO2 exponential map
    R = so3_exp(omega)

    # Translation-orientation coupling
    if abs(omega) > 1e-8:
        V = np.sin(omega) / omega * np.eye(2) + (1 - np.cos(omega)) / omega * np.array([[0, -1], [1, 0]])
    else:
        V = ((1 - omega ** 2 / 6 + omega ** 4 / 120) * np.eye(2) 
           + (omega / 2 - omega ** 3 / 24) * np.array([[0, -1], [1, 0]]))

    X = np.block([[R, V @ a.reshape((2, 1)), V @ b.reshape((2, 1))], 
                  [np.zeros((2, 2)),     np.eye(2)]])
    
    return X

def se23_log(xiwedge):
   """
   Maps from se23 Lie group to the parametrization of its Lie algebra se23
   
   :param xiwedge: 3x3 se23 Lie algebra element
   """
   omega = so3_log(xiwedge[0:3, 0:3])

   if abs(omega) > 1e-8:
       V_inv = omega / 2 * (np.cos(omega / 2) / np.sin(omega / 2) * np.eye(2) + np.array([[0, 1], [-1, 0]]))
   else:
       V_inv = 1 / 2 * ((2 - omega ** 2 / 6 - omega ** 4 / 360 - 2 * omega ** 6 / 30240 + omega ** 8 / 604800) * np.eye(2) + omega * np.array([[0, -1], [1, 0]]))

   a = V_inv @ xiwedge[0:3, 3]
   b = V_inv @ xiwedge[0:3, 4]

   xi = np.vstack((a.reshape((3, 1)), b.reshape((3, 1)), np.reshape(omega, (3,1))))

   return xi

def se23_Ad(X):
    """
    Maps from se23 to se23 defined by (Ad_X xi)**wedge = X xi**wedge X**-1
    Shifts tangent spaces

    :param X: se23 Lie group element
    """
    X_t = X[0:3, 4].reshape((3, 1))
    X_v = X[0:3, 3].reshape((3, 1))
    X_R = X[0:3, 0:3]

    J = np.array([[0, -1], [1, 0]])
    v_circ = -J @ X_v
    t_circ = -J @ X_t

    Ad = np.block([[X_R, np.zeros((2, 2)), v_circ],
                   [np.zeros((2, 2)), X_R, t_circ],
                   [np.zeros((1, 4)),           1]])

    return Ad

def se23_ad(xi):
    """
    Derivative of Adjoint at identity
    
    :param xi: SE2(3) Lie algebra element parametrization
    """
    a = xi[0:3].reshape((3, 1))
    b = xi[3:6].reshape((3, 1))
    omega = xi[6:9]

    J = np.array([[0, -1], [1, 0]])
    a_circ = -J @ a
    b_circ = -J @ b

    omega_wedge = so3_wedge(omega)

    ad = np.block([[omega_wedge, np.zeros((3, 3)), a_circ],
                   [np.zeros((3, 3)), omega_wedge, b_circ],
                   [np.zeros((1, 4)),                   0]])
    
    return ad