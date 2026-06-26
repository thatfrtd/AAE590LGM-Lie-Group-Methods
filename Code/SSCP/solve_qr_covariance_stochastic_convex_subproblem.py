import numpy as np
import cvxpy as cp
import cvxpygen as cpg
import time
from Code.SSCP.qr_covariance import qr_pos_diag, X_Psqrt_noG, qr_derivative, tril, triu, tril_cp, triu_cp, tril_to_mat, tril_to_mat_cp, triu_to_mat

def solve_qr_covariance_stochastic_convex_subproblem(opt_prob, ptr_ops, x_ref, u_ref, S_ref, L_ref, X_ref):
    '''
    Solve qr covariance based stochastic convex subproblem for Stochastic Sequential Convex Programming

    No parameter variables
    Too many cvxpy parameters so doesn't use disciplined convex programming :(

    :param opt_prob: DeterministicProblem object
    :param ptr_ops: ptr algorithm options {w_vc, w_tr, alpha_x, alpha_u, solver}
    :param x_ref: reference state trajectory
    :param u_ref: reference controls
    '''

    # Extract problem info
    nx = opt_prob.n["x"]
    neps = opt_prob.n["eps"]
    nu = opt_prob.n["u"]
    N = opt_prob.n["N"]

    # Adjust for ZOH-ness of input
    u_ref_adj = np.hstack((u_ref, np.zeros((nu, 1))))

    ## Create Convex Problem
    # Define variables
    x = cp.Variable((nx, N), name= 'x') # Nominal state
    S = cp.Variable((neps * int((neps + 1) / 2), N), name = 'S') # State dispersion - telling it is PSD makes it defined over lower dimensional space of PSD matrices
    u = cp.Variable((nu, N - 1), name='u') # Nominal control
    L = cp.Variable((nu, neps, N - 1), name = 'L') # K @ P_hat
    #L = np.zeros((nu, neps, N - 1))
    eta = cp.Variable((N,), name = 'eta')
    V = cp.Variable((nx, N - 1), name = 'V')
    V_S = cp.Variable((neps * int((neps + 1) / 2), N - 1), name = 'V_S')
    v_0 = cp.Variable((nx,), name = 'v_0')
    v_N = cp.Variable((nx,), name = 'v_N')
    v_S_0 = cp.Variable((neps * int((neps + 1) / 2),), name = 'v_S_0')
    v_S_N = cp.Variable((1,), name = "v_S_N")
    if opt_prob.n["ncvx_k"] > 0:
        v_prime = cp.Variable((opt_prob.n["ncvx_k"]), name = "v_prime")
    else:
        v_prime = 0

    # Define objective
    objective = cp.Minimize(opt_prob.objective(x, u, S, L, x_ref, u_ref))
    virtual_control_cost = cp.Minimize(virtual_control_cost_func(ptr_ops["w_vc"], V, V_S, v_0, v_N, v_S_0, v_S_N, v_prime))
    trust_region_cost = cp.Minimize(trust_region_cost_func(ptr_ops["w_tr"], eta))
    augmented_objective = objective + virtual_control_cost + trust_region_cost

    ## Define constraints
    dynamics_constraints = []
    convex_constraints = []
    nonconvex_constraints = []
    trust_region_constraints = []
    nc_count = 0
    for k in range(N):
        if k < N - 1:
            # Nominal dynamics
            dynamics_constraints.append(x[:, k + 1] == opt_prob.disc["A_k"][:, :, k] @ x[:, k] + opt_prob.disc["B_k"][:, :, k] @ u[:, k] + opt_prob.disc["c_k"][:, k] + V[:, k])
            # Sqrt covariance dynamics
            S_k = tril_to_mat_cp(S[:, k])
            S_ref_k = tril_to_mat(S_ref[:, k])
            L_k = L[:, :, k]
            L_ref_k = L_ref[:, :, k]
            X_kp1 = X_Psqrt_noG(opt_prob.disc["A_err_k"][:, :, k], S_k, opt_prob.disc["B_err_k"][:, :, k], L_k)
            X_ref_kp1 = X_Psqrt_noG(opt_prob.disc["A_err_k"][:, :, k], S_ref_k, opt_prob.disc["B_err_k"][:, :, k], L_ref_k)
            Q_kp1, R_kp1 = qr_pos_diag(X_ref_kp1)
            dX_kp1 = X_kp1 - X_ref_kp1
            Q_kp1_ref, S_kp1_ref = qr_pos_diag(X_ref[:, :, k].T)
            S_kp1_ref = S_kp1_ref.T

            dynamics_constraints.append(S[:, k + 1] == triu(S_kp1_ref).flatten() + triu_cp(qr_derivative(dX_kp1.T, Q_kp1, R_kp1)).flatten() + V_S[:, k])

            U_k = u[:, k]
            u_adj = u
        else:
            U_k = np.zeros((nu, 1))
            u_adj = np.zeros((nu, N))
            dX_kp1 = np.zeros_like(dX_kp1)
            S_k = tril_to_mat_cp(S[:, k])
            S_ref_k = tril_to_mat(S_ref[:, k])
            L_k = np.zeros_like(L_k)
            L_ref_k = np.zeros_like(L_k)

        # Convex constraints
        for cc in range(opt_prob.n["cvx"]):
            cc_k = opt_prob.convex_constraints[cc]["k"]
            if np.any(np.isin(cc_k, k)):
                cvx_constraint_func = opt_prob.convex_constraints[cc]["func"]   
                if "trigger" in opt_prob.convex_constraints[cc]:
                    trigger = np.array([opt_prob.convex_constraints[cc]["trigger"](x_ref[:, k], u_ref_adj[:, k])])
                else:
                    trigger = 1
                if opt_prob.convex_constraints[cc]["type"] == "<=":
                    convex_constraints.append(cvx_constraint_func(x[:, k], U_k, S_k, L_k) * trigger <= 0)  
                elif opt_prob.convex_constraints[cc]["type"] == "==":
                    convex_constraints.append(cvx_constraint_func(x[:, k], U_k, S_k, L_k) * trigger == 0)                         

        # Nonconvex constraints
        for nc in range(opt_prob.n["ncvx"]):
            nc_k = opt_prob.nonconvex_constraints[nc]["k"]
            if np.any(np.isin(nc_k, k)):
                ncvx_constraint_coef_func = opt_prob.nonconvex_constraints[nc]["func"]["tay_coef"]
                ncvx_constraint_tay_func = opt_prob.nonconvex_constraints[nc]["func"]["tay"]
                ncvx_coef = ncvx_constraint_coef_func(x_ref, u_ref_adj, k)
                if opt_prob.nonconvex_constraints[nc]["type"] == "<=":
                    nonconvex_constraints.append(ncvx_constraint_tay_func(x, u_adj, x_ref, u_ref_adj, ncvx_coef, k) <= v_prime[nc_count])
                    nonconvex_constraints.append(v_prime[nc_count] <= 0)
                elif opt_prob.nonconvex_constraints[nc]["type"] == "==":
                    nonconvex_constraints.append(ncvx_constraint_tay_func(x, u_adj, x_ref, u_ref_adj, ncvx_coef, k) == v_prime[nc_count])
                nc_count = nc_count + 1

        trust_region_constraints.append(ptr_ops["alpha_x"] * cp.sum_squares(x[:, k] - x_ref[:, k]) 
                                      + ptr_ops["alpha_u"] * cp.sum_squares(U_k - u_ref_adj[:, k])
                                      + ptr_ops["alpha_X"] * cp.norm2(cp.vec(dX_kp1)) <= eta[k])
    
    constraints = dynamics_constraints + \
                    convex_constraints + \
                    nonconvex_constraints + \
                    trust_region_constraints + [
                    opt_prob.initial_bc["x"](x[:, 0], x_ref[:, 0]) == v_0, # Initial condition constraint
                    opt_prob.terminal_bc["x"](x[:, -1], x_ref[:, -1]) == v_N, # Terminal condition constraint
                    opt_prob.initial_bc["S"](S[:, 0]) == v_S_0, # Initial covariance constraint
                    opt_prob.terminal_bc["S"](S[:, -1]) <= v_S_N] # Terminal covariance constraint

    # Define problem
    cvx_prob = cp.Problem(augmented_objective, constraints)

    # Solve
    t0 = time.time()
    val = cvx_prob.solve(solver = ptr_ops["solver"], verbose = ptr_ops["solver_verbose"])
    t1 = time.time()
    if ptr_ops["print_solve_time"]:
        solver_name = ptr_ops["solver"]
        print(f"CVXPY {solver_name} Solve time: {1000 * (t1 - t0):.3f} ms with {val:.3f}")

    # Extract Solution
    x_sol = cvx_prob.var_dict['x'].value
    u_sol = cvx_prob.var_dict['u'].value
    S_sol = cvx_prob.var_dict['S'].value
    L_sol = cvx_prob.var_dict['L'].value

    eta = cvx_prob.var_dict['eta'].value
    V = cvx_prob.var_dict['V'].value
    V_S = cvx_prob.var_dict['V_S'].value
    v_0 = cvx_prob.var_dict['v_0'].value
    v_N = cvx_prob.var_dict['v_N'].value
    v_S_0 = cvx_prob.var_dict['v_S_0'].value
    v_S_N = cvx_prob.var_dict['v_S_N'].value
    if opt_prob.n["ncvx_k"] > 0:
        v_prime = cvx_prob.var_dict['v_prime'].value
    else:
        v_prime = 0

    solve_status = cvx_prob.status

    J = opt_prob.objective(x_sol, u_sol, S_sol, L_sol, x_ref, u_ref).value
    J_ref = opt_prob.objective(x_ref, u_ref, S_ref, L_ref, x_ref, u_ref).value
    J_tr = trust_region_cost_func(ptr_ops["w_tr"], eta)[0]
    J_vc = virtual_control_cost_func(ptr_ops["w_vc"], V, V_S, v_0, v_N, v_S_0, v_S_N, v_prime).value

    delta_X_inf_norm_sum = 0
    dX = 0
    for k in range(N - 1):
        X_kp1 = X_Psqrt_noG(opt_prob.disc["A_err_k"][:, :, k], tril_to_mat(S_sol[:, k]), opt_prob.disc["B_err_k"][:, :, k], L_sol[:, :, k])
        dX_kp1 = X_kp1 - X_ref[:, :, k]
        delta_X_inf_norm_sum = delta_X_inf_norm_sum + np.linalg.norm(dX_kp1.flatten(), np.inf)
        dX = dX + np.linalg.norm(dX_kp1, 'fro')

    dx = np.sum(np.linalg.norm(x_sol - x_ref, 2, 0))
    du = np.sum(np.linalg.norm(u_sol - u_ref, 2, 0))

    # Calculate solution information
    sol_info = {'status':solve_status,
                'vd':V,
                'vd_S':V_S,
                'vs':v_prime,
                'vbc_0':v_0,
                'vbc_N':v_N,
                'J':J,
                'J_tr':J_tr,
                'J_vc':J_vc,
                'dJ':100 * (J - J_ref) / J_ref,
                'dx':dx,
                'du':du,
                'dX':dX,
                'delta':dx + du + dX,
                'eta':np.sum(eta),
                'eta_x':ptr_ops["alpha_x"] * np.sum(np.linalg.norm(x_sol - x_ref, 2, 0) ** 2),
                'eta_u':ptr_ops["alpha_u"] * np.sum(np.linalg.norm(u_sol - u_ref, 2, 0) ** 2),
                'eta_X':ptr_ops["alpha_X"] * delta_X_inf_norm_sum}
    
    return x_sol, u_sol, S_sol, L_sol, sol_info

def virtual_control_cost_func(w_vc, V, V_S, v_0, v_N, v_S_0, v_S_N, v_prime):
        J_vc = w_vc * (cp.sum(cp.norm(V, 1, 0)) + cp.sum(cp.norm(V_S, 1, 0)) + cp.norm(v_0, 1) + cp.norm(v_N, 1) + cp.norm(v_S_0, 1) + cp.norm(v_S_N, 1) + cp.norm(v_prime, 1))
        return J_vc

def trust_region_cost_func(w_tr, eta):
    J_tr = w_tr @ eta.T
    return J_tr