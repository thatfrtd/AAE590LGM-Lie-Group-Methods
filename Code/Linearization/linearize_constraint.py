import numpy as np
import casadi as cas
import cvxpy as cp

def linearize_constraint(constraint, nx, nu, var_type, var_ind):
    #LINEARIZE_CONSTRAINT Summary of this function goes here
    #   Detailed explanation goes here

    x_sym = cas.SX.sym("x", nx, 1)
    u_sym = cas.SX.sym("u", nu, 1)

    if var_type == "x":
        x_selected = cas.blockcat([[x_sym[i, 0] for i in var_ind]])
        constraint_jac = cas.jacobian(constraint(x_sym, u_sym), x_sym)
        constraint_jacobian = cas.Function("constraint_jac_x", [x_sym, u_sym], [constraint_jac], ["x", "u"], ["jac"])
        taylor_ceof_func = lambda x_ref, u_ref, k : np.concatenate((np.array(constraint(x_ref[:, k], u_ref[:, k])).flatten(), np.array(constraint_jacobian(x_ref[:, k], u_ref[:, k])).flatten(), np.array(constraint_jacobian(x_ref[:, k], u_ref[:, k]) @ x_ref[:, k]).flatten()))
        taylor_func = lambda x, u, x_ref, u_ref, coef, k : coef[0] + cp.reshape(coef[1:(np.prod(constraint_jac.shape) + 1)], constraint_jac.shape, 'F') @ x[:, k] - coef[(np.prod(constraint_jac.shape)  + 1) : (np.prod(constraint_jac.shape) + 1 + 1)]
        linearized_constraint_func = lambda x, u, x_ref, u_ref, k : np.array(constraint(x_ref[:, k], u_ref[:, k])) + np.array(constraint_jacobian(x_ref[:, k], u_ref[:, k])) @ (x[var_ind, k] - x_ref[var_ind, k])
        linearized_constraint = {"val":linearized_constraint_func, "tay":taylor_func, "tay_coef":taylor_ceof_func}
    elif var_type == "u":
        u_selected = cas.blockcat([[u_sym[i, 0] for i in var_ind]])
        constraint_jac = cas.jacobian(constraint(x_sym, u_sym), u_sym)
        constraint_jacobian = cas.Function("constraint_jac_u", [x_sym, u_sym], [constraint_jac], ["x", "u"], ["jac"])
        taylor_ceof_func = lambda x_ref, u_ref, k : np.concatenate((np.array(constraint(x_ref[:, k], u_ref[:, k])).flatten(), np.array(constraint_jacobian(x_ref[:, k], u_ref[:, k])).flatten(), np.array(constraint_jacobian(x_ref[:, k], u_ref[:, k]) @ u_ref[:, k]).flatten()))
        taylor_func = lambda x, u, x_ref, u_ref, coef, k : coef[0] + cp.reshape(coef[1:(np.prod(constraint_jac.shape) + 1)], constraint_jac.shape, 'F') @ u[:, k] - coef[(np.prod(constraint_jac.shape)  + 1) : (np.prod(constraint_jac.shape) + 1 + 1)]
        linearized_constraint_func = lambda x, u, x_ref, u_ref, k : np.array(constraint(x_ref[:, k], u_ref[:, k])) + np.array(constraint_jacobian(x_ref[:, k], u_ref[:, k])) @ (u[var_ind, k] - u_ref[var_ind, k])
        linearized_constraint = {"val":linearized_constraint_func, "tay":taylor_func, "tay_coef":taylor_ceof_func}

    return linearized_constraint