import numpy as np
import sympy as sym
import casadi as cas
from Code.SO3.SO3_maps import so3_wedge, so3_exp
from Code.SO3.SO3_quat_maps import quat_rot, q_mul, q_pure, q_hat
from Code.SE3.SE3_maps import se3_wedge, se3_ad, se3_inverse
from Code.Dynamics.EulerPoincare import euler_poincare

def rocket_landing_6DoF_dynamics(g_val, lever_val, I_b_val, alpha_val, jit):
    # Constants
    g_vec = cas.SX([[0], [0], [-g_val]])
    lever = cas.SX(lever_val)
    I_b = cas.SX(I_b_val)
    alpha = cas.SX([alpha_val])

    # States
    x = cas.SX.sym('x', 14, 1)
    r = x[0:3]
    v = x[3:6]
    q = x[6:10]
    w = x[10:13]
    m = x[13]

    I_b_m = I_b * m # Recover current moment of inertia

    # Controls
    u = cas.SX.sym('u', 4, 1)
    T = u[0:3] # Thrust force
    tau = u[3] # Torque

    # Kinematics
    r_dot = v
    q_dot = q_mul(q, q_hat(w))
    #q_dot = 1 / 2 * q_mul(q, q_pure(w))

    # Dynamics
    v_dot = g_vec + quat_rot(q, T) / m
    M = cas.cross(lever, T) + cas.blockcat([[tau*1e-3], [0], [0]]) # Moment
    w_dot = 1 / I_b_m * (cas.cross(w, I_b_m * w) + M) # Assuming principle MoI frame
    T_norm = cas.norm_2(T)
    m_dot = -(T_norm + cas.norm_2(tau*1e-3)) * alpha

    # Package derivatives
    x_dot = cas.vertcat(r_dot, v_dot, q_dot, w_dot, m_dot)

    f = x_dot
    A = cas.jacobian(x_dot, x)
    B = cas.jacobian(x_dot, u)

    # Generate functions
    opts_jit = {'jit':jit}
    f_func = cas.Function("f_func", [x, u], [x_dot], ["x", "u"], ["f"], opts_jit)
    A_func = cas.Function("A_func", [x, u], [A], ["x", "u"], ["A"], opts_jit)
    B_func = cas.Function("B_func", [x, u], [B], ["x", "u"], ["B"], opts_jit)

    # Hackily get sparsity structure of discretized A matrix (B is dense)
    A2 = cas.SX_eye(14) + A @ (cas.SX_eye(14) + A)
    A3 = A2 + A @ A2
    A_sp = np.nonzero(np.array(cas.GenDM_ones(A3.sparsity())))

    # Package
    cont_dyn = {"f":f, "A":A, "B":B, "f_func":f_func, "A_func":A_func, "B_func":B_func, "A_sp":A_sp}

    return cont_dyn

def rocket_landing_6DoF_dynamics_LieGroup(x, u, J_b, g, lever): # ADD MASS
    # Extract states
    X = x[0:9].reshape((3, 3))
    xi = x[9:12]

    # Kinematics
    X_dot = X @ se3_wedge(xi)

    # Dynamics
    u_applied = np.vstack([[u], [-u[1] * lever]]) # [force; torque]
    u_external = se3_inverse(X) @ g # gravity
    u_net = u_applied + u_external
    xi_dot = euler_poincare(xi, J_b, se3_ad(xi), u_net)

    # Package derivatives
    x_dot = np.vstack((X_dot.reshape((9,1)), xi_dot))

    return x_dot

def so3_exp_cas(theta):
    '''
    Exponential map for SO(3) written using CasADi
    '''
    theta_mag = cas.norm_2(theta)
    lambda_vec = theta / theta_mag
    lambda_wedge = so3_wedge(lambda_vec)
    R = cas.SX_eye(3) + lambda_wedge * cas.sin(theta_mag) + lambda_wedge @ lambda_wedge * (1 - cas.cos(theta_mag))

    return R