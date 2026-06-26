import numpy as np
import casadi as cas
from Code.SO3.SO3_maps import so3_wedge

# Lie Group functions
def q_exp(theta):
    '''
    Quaternion exponential map

    :param theta: element of parametrization of Lie algebra (3,)
    '''
    theta_mag = np.linalg.norm(theta)

    if theta_mag > 1e-10:
        theta_vec = theta / np.linalg.norm(theta)
        
        q = np.concatenate([theta_vec * np.sin(theta_mag / 2), [np.cos(theta_mag / 2)]])
    else:
        q = np.array([0, 0, 0, 1])

    return q

def q_log(q):
    '''
    Quaternion logarithmic map

    :param q: quaternion, element of S3 group
    '''
    phi_vee = q_vee(q)
    theta_mag = np.linalg.norm(phi_vee) / 2

    if theta_mag > 1e-10:
        w = q[3] * np.sign(q[3])
        v = q[0:3] * np.sign(w)
        theta = 2 * v * np.atan2(theta_mag, w) / theta_mag
    else:
        theta = np.array([0, 0, 0])

    return theta

def q_log_cas(q):
    '''
    Quaternion logarithmic map
    Written using CasADi

    :param q: quaternion, element of S3 group
    '''
    phi_vee = q_vee(q)
    theta_mag = cas.norm_2(phi_vee) / 2

    #if theta_mag > 1e-10: # Can't evaluate symbolic truth statement...
    w = q[3] * cas.sign(q[3])
    v = q[0:3] * cas.sign(w)
    theta = 2 * v * cas.atan2(theta_mag, w) / (theta_mag + 1e-10)
    #else:
    #    theta = cas.GenDM_zeros(3, 1)

    return theta

def q_hat(theta):
    '''
    Quaternion hat: parametrization -> Lie algebra S3

    :param theta: element of parametrization Lie algebra
    '''
    theta_hat = 1 / 2 * q_pure(theta)

    return theta_hat

def q_vee(phi):
    '''
    Quaternion vee: Lie algebra -> parametrization S3

    :param phi: element of Lie algebra
    '''
    phi_vee = 2 * phi[0:3]

    return phi_vee

# General quaternion helpers
def quat_rot(q, v):
    '''
    QUAT_ROT Summary of this function goes here
    '''
    v_rot = cas.blockcat([[cas.DM_eye(3), cas.GenDM_zeros(3, 1)]]) @ q_mul(q, q_mul(q_pure(v), q_conj(q)))

    return v_rot

def quat_rotmatrix(q):
    w = q[3]
    v = q[0:3]

    R = (w ** 2 - v.T @ v) @ np.eye(3) + 2 * v @ v.T + 2 * w @ so3_wedge(v)
    return R

def q_pure(v):
    '''
    Construct pure quaternion (0 scalar part) from vector
    '''
    p = cas.vcat([v, 0])

    return p

def q_conj(q):
    '''
    Q_CONJ Quaternion conjugate
    '''
    q_star = cas.blockcat([[-q[0:3]], [q[3]]])

    return q_star

def q_conj_cas(q):
    '''
    Q_CONJ Quaternion conjugate
    '''
    q_star = cas.blockcat([[-q[0:3]], [q[3]]])

    return q_star

def q_mul_array(q_A, q_B):
    '''
    %Q_MUL Summary of this function goes here
    %   Detailed explanation goes here
    '''

    N = q_A.shape[1]

    q_prod = cas.zeros(4, N)

    for i in range(N):
        q_prod[:, i] = q_mul(q_A[:, i], q_B[:, i])

    return q_prod

def q_conj_mul_matrix(q):
    '''
    Q_CONJ_MUL_MATRIX
    '''
    q_mm_star = cas.blockcat([[q[3] * cas.DM_eye(3) + cas.skew(q[0:3]), q[0:3]], [-q[0:3], q(4)]])

    return q_mm_star

def q_mul(q_A, q_B):
    '''
    Q_MUL Summary of this function goes here
    '''
    q_mm_A =  q_mul_matrix(q_A)
    q_prod = q_mm_A @ q_B

    return q_prod

def q_mul_matrix(q):
    '''
    Q_MUL_MATRIX Summary of this function goes here
    '''
    q_mm = cas.blockcat([[q[3] * cas.DM_eye(3) + cas.skew(q[0:3]), q[0:3]], [-q[0:3].reshape((1, 3)), q[3]]])

    return q_mm