import numpy as np
import sympy as sym
import casadi as cas
from Code.SO2.SO2_maps import so2_wedge, so2_exp
from Code.SE2.SE2_maps import se2_wedge, se2_ad, se2_inverse
from Code.Dynamics.EulerPoincare import euler_poincare

def rocket_landing_3DoF_dynamics(g_val, lever_val, I_b_val, alpha_val, jit):
    # Constants
    g_vec = cas.SX([[0], [-g_val]])
    lever = cas.SX([lever_val])
    I_b = cas.SX([I_b_val])
    alpha = cas.SX([alpha_val])

    # States
    x = cas.SX.sym('x', 7, 1)
    r = x[0:2]
    v = x[2:4]
    theta = x[4:5]
    w = x[5:6]
    m = x[6:7]

    # Controls
    u = cas.SX.sym('u', 2, 1)
    T = u[0:2]

    # Kinematics
    r_dot = v
    theta_dot = w

    # Dynamics
    R = so2_exp_cas(theta)
    v_dot = g_vec + R @ T / m
    w_dot = -T[1, 0] * lever / (I_b * m)
    T_norm = cas.norm_2(T)
    m_dot = -T_norm * alpha

    # Package derivatives
    x_dot = cas.vertcat(r_dot, v_dot, theta_dot, w_dot, m_dot)

    f = x_dot
    A = cas.jacobian(x_dot, x)
    B = cas.jacobian(x_dot, u)

    # Generate functions
    opts_jit = {'jit':jit}
    f_func = cas.Function("f_func", [x, u], [x_dot], ["x", "u"], ["f"], opts_jit)
    A_func = cas.Function("A_func", [x, u], [A], ["x", "u"], ["A"], opts_jit)
    B_func = cas.Function("B_func", [x, u], [B], ["x", "u"], ["B"], opts_jit)

    # Hackily get sparsity structure of discretized A matrix (B is dense)
    A2 = cas.SX_eye(7) + A @ (cas.SX_eye(7) + A)
    A3 = A2 + A @ A2
    A_sp = np.nonzero(np.array(cas.GenDM_ones(A3.sparsity())))

    # Package
    cont_dyn = {"f":f, "A":A, "B":B, "f_func":f_func, "A_func":A_func, "B_func":B_func, "A_sp":A_sp}

    return cont_dyn

def rocket_landing_3DoF_dynamics_LieGroup(x, u, J_b, g, lever): # ADD MASS
    # Extract states
    X = x[0:9].reshape((3, 3))
    xi = x[9:12]

    # Kinematics
    X_dot = X @ se2_wedge(xi)

    # Dynamics
    u_applied = np.vstack([[u], [-u[1] * lever]]) # [force; torque]
    u_external = se2_inverse(X) @ g # gravity
    u_net = u_applied + u_external
    xi_dot = euler_poincare(xi, J_b, se2_ad(xi), u_net)

    # Package derivatives
    x_dot = np.vstack((X_dot.reshape((9,1)), xi_dot))

    return x_dot

def so2_exp_cas(theta):
    R = cas.blockcat([[cas.cos(theta), -cas.sin(theta)], 
                      [cas.sin(theta),  cas.cos(theta)]])
    return R