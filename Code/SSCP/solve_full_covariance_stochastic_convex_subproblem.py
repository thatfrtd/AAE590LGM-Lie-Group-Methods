import numpy as np
import cvxpy as cp
import cvxpygen as cpg
import time

def solve_full_covariance_stochastic_convex_subproblem(opt_prob, cvx_prob, ptr_ops, x_ref, u_ref):
    '''
    Solve full covariance based stochastic convex subproblem for Stochastic Sequential Convex Programming

    No parameter variables

    :param opt_prob: DeterministicProblem object
    :param cvx_prob: CVXPy problem object (can be list)
    :param ptr_ops: ptr algorithm options {w_vc, w_tr, alpha_x, alpha_u, solver}
    :param x_ref: reference state trajectory
    :param u_ref: reference controls
    '''

    # Extract problem info
    nx = opt_prob.n["x"]
    nu = opt_prob.n["u"]
    N = opt_prob.n["N"]

    # Adjust for ZOH-ness of input
    u_ref_adj = np.hstack((u_ref, np.zeros((nu, 1))))

    ## Create Convex Problem
    # Define variables
    if str(type(cvx_prob)) == "<class 'list'>":
        #print("Creating problem")
        x = cp.Variable((nx, N), name='x') # Nominal state
        P_hat = cp.Variable((nx * nx, N), PSD = True) # State dispersion - telling it is PSD makes it defined over lower dimensional space of PSD matrices
        U = cp.Variable((nu * nx, N - 1), name = 'U') # K @ P_hat
        u = cp.Variable((nu, N - 1), name='u') # Nominal control
        Y = cp.Variable((nu * nu, N - 1), name = 'Y', PSD = True) # Control covariance - telling it is PSD makes it defined over lower dimensional space of PSD matrices
        tau = cp.Variable((N,), name = 'tau', nonneg = True) # Upper bound on sqrt of largest eigenvalue of Y
        eta = cp.Variable((N,), name = 'eta')
        V = cp.Variable((nx, N - 1), name = 'V')
        v_tau = cp.Variable((N,), name = 'v_tau')
        v_0 = cp.Variable((nx, 1), name = 'v_0')
        v_N = cp.Variable((nx, 1), name = 'v_N')
        if opt_prob.n["ncvx_k"] > 0:
            v_prime = cp.Variable((opt_prob.n["ncvx_k"]), name = "v_prime")

        # Define parameters
        x_ref_param = cp.Parameter((nx, N), name = 'x_ref')
        u_ref_param = cp.Parameter((nu, N), name = 'u_ref')
        tau_ref_param = cp.Parameter((N,), name = "tau_ref")

        # Define objective
        objective = cp.Minimize(opt_prob.objective(x, u, x_ref, u_ref))
        virtual_control_cost = cp.Minimize(virtual_control_cost_func(ptr_ops["w_vc"], V, v_tau, v_0, v_N, v_prime))
        trust_region_cost = cp.Minimize(trust_region_cost_func(ptr_ops["w_tr"], eta))
        control_cov_cost = cp.Minimize(control_cov_cost_func(ptr_ops["eps_Y"], Y, tau, opt_prob.quant))
        augmented_objective = objective + virtual_control_cost + trust_region_cost

        ## Define constraints
        Ak_params = []
        Bk_params = []
        ck_param = cp.Parameter((nx, N - 1), name = 'ck')
        cvx_trigger_param = []
        dynamics_constraints = []
        convexification_constraints = []
        convex_constraints = []
        nonconvex_constraints = []
        trust_region_constraints = []
        nc_count = 0
        for k in range(N):
            if k < N - 1:
                # Dynamics constraints
                Ak_params.append(cp.Parameter((nx, nx), name = "Ak_" + str(k)))
                Bk_params.append(cp.Parameter((nx, nu), name = "Bk_" + str(k)))

                dynamics_constraints.append(x[:, k + 1] == Ak_params[k] @ x[:, k] + Bk_params[k] @ u[:, k]  + ck_param[:, k] + V[:, k])
                dynamics_constraints.append(P_hat[:, k + 1] == Ak_params[k] @ P_hat[:, k] @ Ak_params[k].T + Bk_params[k] @ U[:, k] @ Ak_params[k].T + Ak_params[k] @ U[:, k].T @ Bk_params[k].T + Bk_params[k] @ Y[:, k] @ Bk_params[k].T + opt_prob.Q_hat[:, :, k + 1])

                U_k = u[:, k]
                u_adj = u
            else:
                U_k = np.zeros((nu, 1))
                u_adj = np.zeros((nu, N))

            # Convexification constraints
            convexification_constraints.append([[P_hat[:, k], U[:, k]], 
                                                [U[:, k],     Y[:, k]]] >> 0) # Schur's lemma constraint making Y the upper bound on control covariance
            convexification_constraints.append(cp.lambda_max(Y[:, k]) - tau_ref_param[k] ^ 2 - 2 * tau_ref_param[k] * (tau[k] - tau_ref_param[k]) <= v_tau[k]) # linearized constraint on tau being upper bound on sqrt of max Y eigenvalue

            # Convex constraints
            for cc in range(opt_prob.n["cvx"]):
                cc_k = opt_prob.convex_constraints[cc]["k"]
                if np.any(np.isin(cc_k, k)):
                    cvx_constraint_func = opt_prob.convex_constraints[cc]["func"]   
                    if "trigger" in opt_prob.convex_constraints[cc]:
                        trigger = cp.Parameter((1,), name = f"cvx_trigger_{cc}_{k}")
                        cvx_trigger_param.append(trigger)
                    else:
                        trigger = 1
                    if opt_prob.convex_constraints[cc]["type"] == "<=":
                        convex_constraints.append(cvx_constraint_func(x[:, k], U_k, tau[k]) * trigger <= 0)  
                    elif opt_prob.convex_constraints[cc]["type"] == "==":
                        convex_constraints.append(cvx_constraint_func(x[:, k], U_k, tau[k]) * trigger == 0)                         

            # Nonconvex constraints
            for nc in range(opt_prob.n["ncvx"]):
                nc_k = opt_prob.nonconvex_constraints[nc]["k"]
                if np.any(np.isin(nc_k, k)):
                    ncvx_constraint_coef_func = opt_prob.nonconvex_constraints[nc]["func"]["tay_coef"]
                    ncvx_constraint_tay_func = opt_prob.nonconvex_constraints[nc]["func"]["tay"]
                    coef_ex = ncvx_constraint_coef_func(x_ref, u_ref_adj, k)
                    ncvx_coef = cp.Parameter(coef_ex.shape, "ncvx_constraint_coef_" + str(nc_count))
                    if opt_prob.nonconvex_constraints[nc]["type"] == "<=":
                        nonconvex_constraints.append(ncvx_constraint_tay_func(x, u_adj, tau, x_ref_param, u_ref_param, tau_ref_param, ncvx_coef, k) <= v_prime[nc_count])
                        nonconvex_constraints.append(v_prime[nc_count] <= 0)
                    elif opt_prob.nonconvex_constraints[nc]["type"] == "==":
                        nonconvex_constraints.append(ncvx_constraint_tay_func(x, u_adj, tau, x_ref_param, u_ref_param, tau_ref_param, ncvx_coef, k) == v_prime[nc_count])
                    nc_count = nc_count + 1

            trust_region_constraints.append(ptr_ops["alpha_x"] * cp.sum_squares(x[:, k] - x_ref_param[:, k]) 
                                          + ptr_ops["alpha_u"] * cp.sum_squares(U_k - u_ref_param[:, k])
                                          + ptr_ops["alpha_tau"] * cp.norm_inf(tau[k] - tau_ref_param[:, k]) <= eta[k])
        
        constraints = dynamics_constraints + \
                        convex_constraints + \
                        nonconvex_constraints + \
                        trust_region_constraints + [
                        opt_prob.initial_bc(x[:, 0], x_ref_param[:, 0]) == v_0, # Initial condition constraint
                        opt_prob.terminal_bc(x[:, -1], x_ref_param[:, -1]) == v_N, # Terminal condition constraint
                        P_hat[:, 0] == opt_prob.P_hat_minus_0 + opt_prob.Q_hat_0, # Initial covariance constraint
                        opt_prob.P_f - P_hat[:, -1] - opt_prob.P_tilde[:, :, -1] >> 0]
    
        # Define problem
        cvx_prob = cp.Problem(augmented_objective, constraints)

    # Set parameters
    nc_count = 0
    for k in range(N):
        if k < N - 1:
            cvx_prob.param_dict["Ak_" + str(k)].value = opt_prob.disc["A_k"][:, :, k]
            cvx_prob.param_dict["Bk_" + str(k)].value = opt_prob.disc["B_k"][:, :, k]
        
        for cc in range(opt_prob.n["cvx"]):
            cc_k = opt_prob.convex_constraints[cc]["k"]
            if np.any(np.isin(cc_k, k)):
                cvx_constraint_func = opt_prob.convex_constraints[cc]["func"]   
                if "trigger" in opt_prob.convex_constraints[cc]:
                    cvx_prob.param_dict[f"cvx_trigger_{cc}_{k}"].value = opt_prob.convex_constraints[cc]["trigger"](x_ref[:, k], u_ref_adj[:, k])

        for nc in range(opt_prob.n["ncvx"]):
            nc_k = opt_prob.nonconvex_constraints[nc]["k"]
            if np.any(np.isin(nc_k, k)):
                ncvx_constraint_coef_func = opt_prob.nonconvex_constraints[nc]["func"]["tay_coef"]
                if "trigger" in opt_prob.nonconvex_constraints[nc]:
                    trigger = opt_prob.nonconvex_constraints[nc]["trigger"](x_ref[:, k], u_ref_adj[:, k])
                else:
                    trigger = 1
                cvx_prob.param_dict["ncvx_constraint_coef_" + str(nc_count)].value = ncvx_constraint_coef_func(x_ref, u_ref_adj, k) * trigger
                nc_count = nc_count + 1

    cvx_prob.param_dict["ck"].value = opt_prob.disc["c_k"]

    cvx_prob.param_dict["x_ref"].value = x_ref
    cvx_prob.param_dict["u_ref"].value = u_ref_adj

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

    eta = cvx_prob.var_dict['eta'].value
    V = cvx_prob.var_dict['V'].value
    v_0 = cvx_prob.var_dict['v_0'].value
    v_N = cvx_prob.var_dict['v_N'].value
    if opt_prob.n["ncvx_k"] > 0:
        v_prime = cvx_prob.var_dict['v_prime'].value
    else:
        v_prime = 0

    solve_status = cvx_prob.status

    J = opt_prob.objective(x_sol, u_sol, x_ref, u_ref).value
    J_ref = opt_prob.objective(x_ref, u_ref, x_ref, u_ref).value
    J_tr = trust_region_cost_func(ptr_ops["w_tr"], eta)[0]
    J_vc = virtual_control_cost_func(ptr_ops["w_vc"], V, v_0, v_N, v_prime).value

    dx = np.sum(np.linalg.norm(x_sol - x_ref, 2))
    du = np.sum(np.linalg.norm(u_sol - u_ref, 2))

    # Calculate solution information
    sol_info = {'status':solve_status,
                'vd':V,
                'vs':v_prime,
                'vbc_0':v_0,
                'vbc_N':v_N,
                'J':J,
                'J_tr':J_tr,
                'J_vc':J_vc,
                'dJ':100 * (J - J_ref) / J_ref,
                'dx':dx,
                'du':du,
                'delta':dx + du,
                'eta':np.linalg.norm(eta),
                'eta_x':ptr_ops["alpha_x"] * np.sum(np.linalg.norm(x_sol - x_ref, 2)),
                'eta_u':ptr_ops["alpha_u"] * np.sum(np.linalg.norm(u_sol - u_ref, 2))}
    
    return x_sol, u_sol, cvx_prob, sol_info

def virtual_control_cost_func(w_vc, V, v_tau, v_0, v_N, v_prime):
        #J_vc = w_vc * (cp.sum(cp.sum_squares(V)) + cp.sum(cp.sum_squares(v_0)) + cp.sum(cp.sum_squares(v_N)))
        J_vc = w_vc * (cp.sum(cp.norm(V, 1, 0)) + cp.sum(v_tau) + cp.norm(v_0, 1) + cp.norm(v_N, 1) + cp.norm(v_prime, 1))
        return J_vc

def trust_region_cost_func(w_tr, eta):
    J_tr = w_tr @ eta.T
    return J_tr

def trust_region_cost_func_explicit(w_tr, x, x_ref, u, u_ref):
    J_tr = (cp.sum(cp.square(x - x_ref), 0) + cp.sum(cp.square(u - u_ref), 0)) @ w_tr.T
    return J_tr

def virtual_control_cost_func_np(w_vc, V, v_0, v_N, v_prime):
    J_vc = w_vc * (np.sum(np.linalg.norm(V, 1, 0)) + np.linalg.norm(v_0, 1) + np.linalg.norm(v_N, 1) + np.linalg.norm(v_prime, 1))
    return J_vc

def control_cov_cost_func(eps_Y, Y, tau, quant):
    J_sigma = 0
    for k in range(Y.shape[2]):
        J_sigma = J_sigma + quant * tau[k] + eps_Y * cp.trace(Y[:, :, k])
    return J_sigma