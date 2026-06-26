import numpy as np
from Code.SO3.SO3_maps import so3_wedge, so3_vee, so3_exp, so3_log, so3_ljac, so3_inv_ljac
from Code.SO3.SO3_quat_maps import quat_rot, q_mul, q_conj_cas, q_exp, q_log, q_conj

def se3_quat_compose(X1, X2):
    """
    Composition of elements of SE(3)
    
    :param X1: element of SE(3)
    :param X2: element of SE(3)
    """
    X1_t = X1[0:3].reshape((3, 1))
    X1_q = X1[3:7].reshape((4, 1))
    
    X2_t = X2[0:3].reshape((3, 1))
    X2_q = X2[3:7].reshape((4, 1))

    X12_q = q_mul(X1_q, X2_q)
    X12_t = quat_rot(X1_q, X2_t) + X1_t
    
    X12 = np.vstack((X12_t, X12_q))

    return X12

def se3_quat_inverse(X):
    """
    Inverse of SE(3) element
    
    :param X: element of SE(3)
    """
    X_t = X[0:3].reshape((3, 1))
    X_q = X[3:7].reshape((4, 1))

    X_inv_q = q_conj(X_q)
    X_inv_t = -quat_rot(X_inv_q, X_t)

    X_inv = np.vstack([X_inv_t, X_inv_q])

    return X_inv

def se3_quat_exp(xi):
    """
    Maps from parametrization of se3 Lie algebra to SE3 Lie group
    
    :param xi: twist vector (v_x, v_y, omega) in R3
    """
    # Extract elements
    v = xi[0:3]
    omega = xi[3:6]

    # SO3 exponential map
    q = q_exp(omega.flatten()).reshape((4, 1))

    # Translation-orientation coupling
    V = so3_ljac(omega)
    
    X = np.vstack((V @ v.reshape((3, 1)), q.reshape((4, 1))))
    
    return X

def se3_quat_log(X):
   """
   Maps from SE3 Lie group to the parametrization of its Lie algebra se3
   
   :param X: 7x1 se3 Lie group element
   """
   omega = q_log(X[3:7])

   v = so3_inv_ljac(omega) @ X[0:3]

   xi = np.vstack((v.reshape((3, 1)), np.reshape(omega, (3,1))))

   return xi