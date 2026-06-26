import numpy as np
from Code.SO3.SO3_maps import so3_wedge, so3_vee, so3_exp, so3_log, so3_ljac, so3_inv_ljac
from Code.SO2.SO2_maps import so2_ljac, so2_inv_ljac

def se3_compose(X1, X2):
    """
    Composition of elements of SE(3)
    
    :param X1: element of SE(3)
    :param X2: element of SE(3)
    """
    X1_t = X1[0:3, 3].reshape((3, 1))
    X1_R = X1[0:3, 0:3]
    
    X2_t = X2[0:3, 3].reshape((3, 1))
    X2_R = X2[0:3, 0:3]

    X12_R = X1_R @ X2_R
    X12_t = X1_R @ X2_t + X1_t
    
    X12 = np.block([[X12_R, X12_t], 
                    [np.zeros((1, 3)), 1]])

    return X12

def se3_inverse(X):
    """
    Inverse of SE(3) element
    
    :param X: element of SE(3)
    """
    X_t = X[0:3, 3].reshape((3, 1))
    X_R = X[0:3, 0:3]

    X_inv_t = -X_R.T @ X_t
    X_inv_R = X_R.T

    X_inv = np.block([[X_inv_R, X_inv_t],
                      [np.zeros((1, 3)), 1]])

    return X_inv

def se3_wedge(xi):
    """
    Maps from se3 real number parametrization to se3 Lie algebra
    
    :param xi: twist vector (v_x, v_y, v_z, omega_x, omega_y, omeg_z) in R3
    """
    xiwedge = np.block([[so3_wedge(xi[3]), xi[0:3].reshape((3, 1))], 
                        [np.zeros((1, 3)),       0]])

    return xiwedge

def se3_vee(xiwedge):
    """
    Maps from se3 Lie algebra to its real number parametrization 
    
    :param xiwedge: 3x3 se3 Lie algebra element
    """
    
    xi = np.block([xiwedge[0:3, 3], so3_vee(xiwedge[0:3, 0:3])])

    return xi

def se3_exp(xi):
    """
    Maps from parametrization of se3 Lie algebra to SE3 Lie group
    
    :param xi: twist vector (v_x, v_y, omega) in R3
    """
    # Extract elements
    v = xi[0:3]
    omega = xi[3:6]

    # SO3 exponential map
    R = so3_exp(omega)

    # Translation-orientation coupling
    X = np.block([[R,    so3_ljac(omega) @ v.reshape((3, 1))], 
                  [np.zeros((1, 3)),     1]])
    
    return X

def se3_log(xiwedge):
   """
   Maps from SE3 Lie group to the parametrization of its Lie algebra se3
   
   :param xiwedge: 3x3 se3 Lie algebra element
   """
   omega = so3_log(xiwedge[0:3, 0:3])

   v = so3_inv_ljac(omega) @ xiwedge[0:3, 3]

   xi = np.vstack((v.reshape((3, 1)), np.reshape(omega, (3,1))))

   return xi

def se3_Ad(X):
    """
    Maps from se3 to se3 defined by (Ad_X xi)**wedge = X xi**wedge X**-1
    Shifts tangent spaces

    :param X: SE3 Lie group element
    """
    R = X[0:3, 0:3]
    t = X[0:3, 3]

    Ad = np.copy(X)
    Ad[0:3, 3] = np.block([[R, -so3_wedge(t) @ R],
                           [np.zeros((3, 3)), R]])

    return Ad

def se3_ad(xi):
    """
    Derivative of Adjoint at identity
    
    :param xi: se(3) Lie algebra element parametrization
    """
    rho_wedge = se3_wedge(xi[0:3])
    theta_wedge = se3_wedge(xi[3:6])
    ad = np.block([[theta_wedge, rho_wedge], [np.zeros((3, 3)), theta_wedge]])
    
    return ad

def se3_Ql(xi):
    """
    Jacobian of velocity element w.r.t angular velocity element

    :param xi: se(3) Lie algebra element parametrization
    """
    rho_wedge = se3_wedge(xi[0:3])
    theta_wedge = se3_wedge
    theta = np.linalg.norm(xi[3:6])
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    Q_l = (1 / 2 @ rho_wedge 
         + ((theta - ctheta) / theta ** 3) @ (theta_wedge @ rho_wedge + rho_wedge @ theta_wedge + theta_wedge @ rho_wedge @ theta_wedge)
         + ((theta ** 2 + 2 * ctheta - 2) / (2 * stheta ** 4)) @ (theta_wedge @ theta_wedge @ rho_wedge + rho_wedge @ theta_wedge @ theta_wedge - 3 * theta_wedge @ rho_wedge @ theta_wedge)
         + ((2 * theta - 3 * stheta + theta * ctheta) / (2 * theta ** 5)) @ (theta_wedge @ rho_wedge @ theta_wedge @ theta_wedge + theta_wedge @ theta_wedge @ rho_wedge @ theta_wedge))           

    return Q_l

def se3_ljac(xi):
    """
    Left Jacobian

    :param xi: se(3) Lie algebra element parametrization
    """
    theta = np.linalg.norm(xi[3:6])
    if theta > 1e-10:
        Q_l = se3_Ql(xi)
        J_l_theta = so2_ljac(xi[3:6])
        J_l = np.block([[J_l_theta, Q_l], [np.zeros((3, 3)), J_l_theta]])
    else:
        J_l = se3_ljac_approx(xi)

    return J_l

def se3_ljac_approx(xi):
    """
    Approximate left Jacobian

    :param xi: se(3) Lie algebra element parametrization
    """
    J_l_approx = np.eye(6) + 1 / 2 * se3_wedge(xi)

    return J_l_approx

def se3_inv_ljac(xi):
    """
    Inverse left Jacobian

    :param xi: se(3) Lie algebra element parametrization
    """
    theta = np.linalg.norm(xi[3:6])
    if theta > 1e-10:
        Q_l = se3_Ql(xi)
        J_l_inv_theta = so2_inv_ljac(xi[3:6])
        J_l_inv = np.block([[J_l_inv_theta, -J_l_inv_theta @ Q_l @ J_l_inv_theta], 
                            [np.zeros((3, 3)), J_l_inv_theta]])
    else:
        J_l_inv = se3_inv_ljac_approx(xi)

    return J_l_inv

def se3_inv_ljac_approx(xi):
    """
    Approximate inverse left Jacobian

    :param xi: se(3) Lie algebra element parametrization
    """
    J_l_inv_approx = np.eye(6) - 1 / 2 * se3_wedge(xi)
    
    return J_l_inv_approx