import numpy as np
import casadi as cas
from Code.SO3.SO3_maps import so3_wedge
from Code.SO3.SO3_quat_maps import quat_rot, quat_rotmatrix, q_conj, q_log, q_log_cas
from Code.SE23.SE23_quat_maps import se23_quat_compose, se23_quat_inverse, se23_quat_compose_cas, se23_quat_exp

def full_right_error(x, x_ref):
    '''
    Right error in Lie algebra
    '''
    Psi_R = np.array(se23_quat_compose(cas.DM([x[0:10]]), se23_quat_inverse(cas.DM([x_ref[0:10]]))))
    psi_r_R = Psi_R[0:3].reshape((3, 1))
    psi_v_R = Psi_R[3:6].reshape((3, 1))
    psi_theta_R = q_log(Psi_R[6:10]).reshape((3, 1))

    psi_R = np.vstack((psi_r_R, psi_v_R, psi_theta_R))

    delta_w = (x[10:13] - x_ref[10:13]).reshape((3, 1))
    delta_m = (x[13] - x_ref[13]).reshape((1, 1))

    eps = np.vstack((psi_R, delta_w, delta_m))

    return eps

def full_right_error_cas(x, x_ref):
    '''
    Right error in Lie algebra
    Written using CasADi
    '''
    Psi_R = se23_quat_compose_cas(x[0:10], se23_quat_inverse(x_ref[0:10]))
    psi_r_R = Psi_R[0:3].reshape((3, 1))
    psi_v_R = Psi_R[3:6].reshape((3, 1))
    psi_theta_R = q_log_cas(Psi_R[6:10]).reshape((3, 1))

    psi_R = cas.vertcat(psi_r_R, psi_v_R, psi_theta_R)

    delta_w = (x[10:13] - x_ref[10:13]).reshape((3, 1))
    delta_m = (x[13] - x_ref[13]).reshape((1, 1))

    eps = cas.vertcat(psi_R, delta_w, delta_m)

    return eps

def full_right_error_rate(g_val, lever_val, I_b_val, alpha_val):
    '''
    Derivative of right error in Lie algebra
    '''
    # Constants
    g_vec = cas.SX([[0], [0], [-g_val]])
    lever = cas.SX(lever_val)
    I_b = cas.SX(I_b_val)
    alpha = cas.SX([alpha_val])

    eps = cas.SX.sym('eps', 13, 1)
    psi_r_R = eps[0:3]
    psi_v_R = eps[3:6]
    psi_theta_R = eps[6:9]
    delta_w = eps[9:12]
    delta_m = eps[12]

    x = cas.SX.sym('x', 14, 1)
    r = x[0:3]
    q = x[6:10]

    x_ref = cas.SX.sym('x_ref', 14, 1)
    r_ref = x_ref[0:3]
    v_ref = x_ref[3:6]
    q_ref = x_ref[6:10]
    w_ref = x_ref[10:13]
    m_ref = x_ref[13]

    I_b_m = I_b * m_ref # Recover current moment of inertia

    u_ref = cas.SX.sym('u_ref', 4, 1)
    T_ref = u_ref[0:3]
    tau_ref = u_ref[3]

    delta_u = cas.SX.sym('delta_u', 4, 1)
    delta_T = delta_u[0:3]
    delta_tau = delta_u[3]

    w_a = delta_T / m_ref - T_ref / (m_ref ** 2) * delta_m

    psi_r_R_dot = psi_v_R + quat_rot(q, cas.skew(quat_rot(q_conj(q_ref), r_ref)) @ delta_w) 
    psi_v_R_dot = cas.skew(g_vec) @ psi_theta_R + quat_rot(q, w_a) + quat_rot(q, cas.skew(quat_rot(q_conj(q_ref), v_ref)) @ delta_w) 
    psi_theta_R_dot = (cas.DM_eye(3) + cas.skew(psi_theta_R)) @ quat_rot(q_ref, delta_w)

    psi_R_dot = cas.vertcat(psi_r_R_dot, psi_v_R_dot, psi_theta_R_dot)

    M = cas.cross(lever, T_ref) + cas.blockcat([[tau_ref*1e-3], [0], [0]]) # Moment
    w_dot = 1 / I_b_m * (cas.cross(w_ref, I_b_m * w_ref) + M) # Assuming principle MoI frame
    w_dot_jac_w = cas.jacobian(w_dot, w_ref)
    w_dot_jac_m = cas.jacobian(w_dot, m_ref)
    w_dot_jac_u = cas.jacobian(w_dot, u_ref)
    delta_w_dot = w_dot_jac_w @ delta_w + w_dot_jac_m @ delta_m + w_dot_jac_u @ delta_u
    
    T_norm = cas.norm_2(T_ref)
    m_dot = -(T_norm + cas.norm_2(tau_ref)) * alpha
    m_dot_jac_u = cas.jacobian(m_dot, u_ref)
    delta_m_dot = m_dot_jac_u @ delta_u

    # Package error state derivative
    eps_dot = cas.vertcat(psi_R_dot, delta_w_dot, delta_m_dot)

    # Generate functions
    eps_dot_func = cas.Function("eps_dot_func", [eps, delta_u, x, x_ref, u_ref], [eps_dot], ["eps", "delta_u", "x", "x_ref", "u_ref"], ["eps_dot"])

    return eps_dot_func

def full_right_error_rate_linearized_ref(g_val, lever_val, I_b_val, alpha_val):
    '''
    Derivative of right error in Lie algebra linearized around reference
    '''
    # Constants
    g_vec = cas.SX([[0], [0], [-g_val]])
    lever = cas.SX(lever_val)
    I_b = cas.SX(I_b_val)
    alpha = cas.SX([alpha_val])

    eps = cas.SX.sym('eps', 13, 1)
    psi_r_R = eps[0:3]
    psi_v_R = eps[3:6]
    psi_theta_R = eps[6:9]
    delta_w = eps[9:12]
    delta_m = eps[12]

    x_ref = cas.SX.sym('x_ref', 14, 1)
    r_ref = x_ref[0:3]
    v_ref = x_ref[3:6]
    q_ref = x_ref[6:10]
    w_ref = x_ref[10:13]
    m_ref = x_ref[13]

    I_b_m = I_b * m_ref # Recover current moment of inertia

    u_ref = cas.SX.sym('u_ref', 4, 1)
    T_ref = u_ref[0:3]
    tau_ref = u_ref[3]

    delta_u = cas.SX.sym('delta_u', 4, 1)
    delta_T = delta_u[0:3]
    delta_tau = delta_u[3]

    w_a = delta_T / m_ref - T_ref / (m_ref ** 2) * delta_m

    psi_r_R_dot = psi_v_R + cas.skew(r_ref) @ quat_rot(q_ref, delta_w)
    psi_v_R_dot = cas.skew(g_vec) @ psi_theta_R + quat_rot(q_ref, w_a) + cas.skew(v_ref) @ quat_rot(q_ref, delta_w)
    psi_theta_R_dot = quat_rot(q_ref, delta_w)

    psi_R_dot = cas.vertcat(psi_r_R_dot, psi_v_R_dot, psi_theta_R_dot)

    M = cas.cross(lever, T_ref) + cas.blockcat([[tau_ref*1e-3], [0], [0]]) # Moment
    w_dot = 1 / I_b_m * (cas.cross(w_ref, I_b_m * w_ref) + M) # Assuming principle MoI frame
    w_dot_jac_w = cas.jacobian(w_dot, w_ref)
    w_dot_jac_m = cas.jacobian(w_dot, m_ref)
    w_dot_jac_u = cas.jacobian(w_dot, u_ref)
    delta_w_dot = w_dot_jac_w @ delta_w + w_dot_jac_m @ delta_m + w_dot_jac_u @ delta_u
    
    T_norm = cas.norm_2(T_ref)
    m_dot = -(T_norm + cas.norm_2(tau_ref)) * alpha
    m_dot_jac_u = cas.jacobian(m_dot, u_ref)
    delta_m_dot = m_dot_jac_u @ delta_u

    # Package error state derivative
    eps_dot = cas.vertcat(psi_R_dot, delta_w_dot, delta_m_dot)
    A = cas.jacobian(eps_dot, eps)
    B = cas.jacobian(eps_dot, delta_u)

    # Generate functions
    eps_dot_func = cas.Function("eps_dot_func", [eps, delta_u, x_ref, u_ref], [eps_dot], ["eps", "delta_u", "x_ref", "u_ref"], ["eps_dot"])
    A_func = cas.Function("A_func", [x_ref, u_ref], [A], ["x_ref", "u_ref"], ["A"])
    B_func = cas.Function("B_func", [x_ref, u_ref], [B], ["x_ref", "u_ref"], ["B"])

    # Package
    cont_dyn = {"f":eps_dot, "A":A, "B":B, "f_func":eps_dot_func, "A_func":A_func, "B_func":B_func}

    return cont_dyn

def full_right_error_rate_linearized_est(eps, x, x_ref, u, u_ref, g):
    '''
    Derivative of right error in Lie algebra linearized around estimate
    '''
    psi_r_R = eps[0:3]
    psi_v_R = eps[3:6]
    psi_theta_R = eps[6:10]
    delta_w = eps[10:13]
    delta_m = eps[13]

    r = x[0:3]
    q = x[6:10]

    r_ref = x_ref[0:3]
    v_ref = x_ref[3:6]
    q_ref = x_ref[6:10]
    m_ref = x[13]

    T_ref = u_ref[0:3]
    delta_T = u[0:3] - T_ref
    delta_tau = u[3] - u_ref[3]
    w_a = delta_T / m_ref - T_ref / m_ref ** 2 * delta_m

    psi_r_R_dot = psi_v_R + quat_rot(q, cas.skew(quat_rot(q_conj(q_ref), r_ref)), delta_w) 
    psi_v_R_dot = cas.skew(g) @ psi_theta_R + quat_rot(q, w_a) + quat_rot(q, cas.skew(quat_rot(q_conj(q_ref), v_ref)), delta_w) 
    psi_theta_R_dot = (np.eye(3) + cas.skew(psi_theta_R)) @ quat_rot(q_ref, delta_w)

    psi_R_dot = np.vstack((psi_r_R_dot, psi_v_R_dot, psi_theta_R_dot))

    #delta_w_dot = 
    #delta_m_dot = 

    eps_dot = np.vstack((psi_R_dot, delta_w_dot, delta_m_dot))

    return eps_dot

def full_compose(x1, x2):
    '''
    Composition of full state
    '''
    x12 = np.zeros((14,))
    x12[0:10] = se23_quat_compose(x1, x2).flatten()
    x12[10:14] = x1[10:14] + x2[10:14]

    return x12

def full_exp(eps):
    '''
    Exponential of full state
    '''
    x = np.zeros((14,))
    x[0:10] = se23_quat_exp(eps[0:9]).reshape((10,))
    x[10:14] = eps[9:13]

    return x

def full_add(x, eps):
    '''
    Combined exponentiation, composition operation
    '''
    x_new = full_compose(full_exp(eps), x)

    return x_new