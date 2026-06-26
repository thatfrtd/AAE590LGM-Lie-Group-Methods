import numpy as np
import casadi as cas
from Code.SO3.SO3_maps import so3_wedge, so3_vee, so3_exp, so3_log, so3_ljac, so3_inv_ljac
from Code.SO3.SO3_quat_maps import quat_rot, q_mul, q_conj_cas, q_exp, q_log

def se23_quat_compose(X1, X2):
    """
    Composition of elements of quaternion based SE2(3)
    
    :param X1: element of SE2(3)
    :param X2: element of SE2(3)
    """
    X1_t = X1[0:3].reshape((3, 1))
    X1_v = X1[3:6].reshape((3, 1))
    X1_q = X1[6:10].reshape((4, 1))
    
    X2_t = X2[0:3].reshape((3, 1))
    X2_v = X2[3:6].reshape((3, 1))
    X2_q = X2[6:10].reshape((4, 1))

    X12_q = q_mul(X1_q, X2_q)
    X12_v = quat_rot(X1_q, X2_v) + X1_v
    X12_t = quat_rot(X1_q, X2_t) + X1_t
    
    X12 = np.vstack((X12_t, X12_v, X12_q))

    return X12

def se23_quat_compose_cas(X1, X2):
    """
    Composition of elements of quaternion based SE2(3)
    Made using CasADi
    
    :param X1: element of SE2(3)
    :param X2: element of SE2(3)
    """
    X1_t = X1[0:3].reshape((3, 1))
    X1_v = X1[3:6].reshape((3, 1))
    X1_q = X1[6:10].reshape((4, 1))
    
    X2_t = X2[0:3].reshape((3, 1))
    X2_v = X2[3:6].reshape((3, 1))
    X2_q = X2[6:10].reshape((4, 1))

    X12_q = q_mul(X1_q, X2_q)
    X12_v = quat_rot(X1_q, X2_v) + X1_v
    X12_t = quat_rot(X1_q, X2_t) + X1_t
    
    X12 = cas.vertcat(X12_t, X12_v, X12_q)

    return X12

def se23_quat_inverse(X):
    """
    Inverse of quaternion based SE2(3) element
    
    :param X: element of SE2(3)
    """
    X_t = X[0:3].reshape((3, 1))
    X_v = X[3:6].reshape((3, 1))
    X_q = X[6:10].reshape((4, 1))

    X_inv_q = q_conj_cas(X_q)
    X_inv_t = -quat_rot(X_inv_q, X_t)
    X_inv_v = -quat_rot(X_inv_q, X_v)

    X_inv = cas.vcat([X_inv_t, X_inv_v, X_inv_q])

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

def se23_quat_exp(xi):
    """
    Maps from parametrization of se23 Lie algebra to se23 Lie group
    
    :param xi: twist vector (a_x, a_y, b_x, b_y, omega) in R3
    """
    # Extract elements
    a = xi[0:3].reshape((3, 1))
    b = xi[3:6].reshape((3, 1))
    omega = xi[6:9].reshape((3, 1))

    # SO3 exponential map
    q = q_exp(omega.flatten()).reshape((4, 1))

    # Translation-orientation
    V = so3_ljac(omega)

    X = np.vstack((V @ a.reshape((3, 1)), V @ b.reshape((3, 1)), q.reshape((4, 1))))
    
    return X

def se23_quat_log(X):
   """
   Maps from se23 Lie group to the parametrization of its Lie algebra se23
   
   :param X: 10x1 se23 Lie group element
   """
   omega = q_log(X[6:10])

   V_inv = so3_inv_ljac(omega)

   a = V_inv @ X[0:3]
   b = V_inv @ X[3:6]

   xi = np.vstack((a, b, omega))

   return xi