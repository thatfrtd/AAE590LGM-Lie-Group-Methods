import numpy as np
import cvxpy as cp
import cvxpygen as cpg
import time

def solve_geodesic_convex_subproblem(opt_prob, cvx_prob, ptr_ops, x_ref, u_ref):
    '''
    Solve geodesic convex subproblem for Sequential Convex Programming

    No parameter variables

    :param opt_prob: DeterministicProblem object
    :param cvx_prob: CVXPy problem object (can be list)
    :param ptr_ops: ptr algorithm options {w_vc, w_tr, alpha_x, alpha_u, solver}
    :param x_ref: reference state trajectory
    :param eps_ref: reference error state trajectory
    :param u_ref: reference controls
    '''

    # Extract problem info
    nx = opt_prob.n["x"]
    neps = opt_prob.n["eps"]
    nu = opt_prob.n["u"]
    N = opt_prob.n["N"]
    Nu = opt_prob.n["Nu"]

    # Adjust for FOH-ness of input
    u_ref_adj = np.hstack((u_ref, np.zeros((nu, 1))))

    ## Create Convex Problem
    # Define variables
    if str(type(cvx_prob)) == "<class 'list'>":
        #print("Creating problem")
        eps = cp.Variable((neps, N), name='eps')
        delta_U = cp.Variable((nu, N - 1), name='delta_U')
        eta = cp.Variable((N,), name = 'eta')
        V = cp.Variable((neps, N - 1), name = 'V')
        v_0 = cp.Variable((neps, 1), name = 'v_0')
        v_N = cp.Variable((neps, 1), name = 'v_N')
        if opt_prob.n["ncvx_k"] > 0:
            v_prime = cp.Variable((opt_prob.n["ncvx_k"]), name = "v_prime")
        else:
            v_prime = 0

        # Define parameters
        x_ref_param = cp.Parameter((nx, N), name = 'x_ref')
        u_ref_param = cp.Parameter((nu, N), name = 'u_ref')
        eps_0_param = cp.Parameter((neps,), name = 'eps_0')
        eps_N_param = cp.Parameter((neps,), name = 'eps_N')

        # Define objective
        objective = cp.Minimize(opt_prob.objective(eps, delta_U + u_ref_param[:, :Nu], x_ref_param, u_ref_param))
        if opt_prob.n["ncvx_k"] > 0:
            virtual_control_cost = cp.Minimize(virtual_control_cost_func(ptr_ops["w_vc"], V, v_0, v_N, v_prime))
        else:
            virtual_control_cost = cp.Minimize(virtual_control_cost_func(ptr_ops["w_vc"], V, v_0, v_N, 0))
        trust_region_cost = cp.Minimize(trust_region_cost_func(ptr_ops["w_tr"], eta))
        augmented_objective = objective + virtual_control_cost + trust_region_cost

        ## Define constraints
        Ak_params = []
        Bk_params = []
        ck_param = cp.Parameter((neps, N - 1), name = 'ck')
        cvx_trigger_param = []
        dynamics_constraints = []
        convex_constraints = []
        nonconvex_constraints = []
        trust_region_constraints = []
        nc_count = 0
        for k in range(N):
            if k < N - 1:
                # Dynamics constraints
                Ak_params.append(cp.Parameter((neps, neps), name = "Ak_" + str(k)))
                Bk_params.append(cp.Parameter((neps, nu), name = "Bk_" + str(k)))

                dynamics_constraints.append(eps[:, k + 1] == Ak_params[k] @ eps[:, k] + Bk_params[k] @ delta_U[:, k]  + ck_param[:, k] + V[:, k])
            
                delta_U_k = delta_U[:, k]
                U_k = delta_U_k + u_ref_param[:, k]
                u_adj = delta_U_k + u_ref_param[:, k]
            else:
                U_k = np.zeros((nu, 1))
                u_adj = np.zeros((nu, N))

            # Convex constraints
            for cc in range(opt_prob.n["cvx"]):
                cc_k = opt_prob.convex_constraints[cc]["k"]
                if np.any(np.isin(cc_k, k)):
                    cvx_constraint_func = opt_prob.convex_constraints[cc]["func"]   
                    if "trigger" in opt_prob.convex_constraints[cc]:
                        trigger = opt_prob.convex_constraints[cc]["trigger"](x_ref[:, k], u_ref_adj[:, k])
                    else:
                        trigger = 1
                    if trigger > 0:
                        if opt_prob.convex_constraints[cc]["type"] == "<=":
                            convex_constraints.append(cvx_constraint_func(eps[:, k], U_k, x_ref_param[:, k]) <= 0)  
                        elif opt_prob.convex_constraints[cc]["type"] == "==":
                            convex_constraints.append(cvx_constraint_func(eps[:, k], U_k, x_ref_param[:, k]) == 0)                         

            # Nonconvex constraints
            for nc in range(opt_prob.n["ncvx"]):
                nc_k = opt_prob.nonconvex_constraints[nc]["k"]
                if np.any(np.isin(nc_k, k)):
                    ncvx_constraint_coef_func = opt_prob.nonconvex_constraints[nc]["func"]["tay_coef"]
                    ncvx_constraint_tay_func = opt_prob.nonconvex_constraints[nc]["func"]["tay"]
                    coef_ex = ncvx_constraint_coef_func(x_ref, u_ref_adj, k)
                    ncvx_coef = cp.Parameter(coef_ex.shape, "ncvx_constraint_coef_" + str(nc_count))
                    if opt_prob.nonconvex_constraints[nc]["type"] == "<=":
                        nonconvex_constraints.append(ncvx_constraint_tay_func(eps, u_adj, x_ref_param, u_ref_param, ncvx_coef, k) <= v_prime[nc_count])
                        nonconvex_constraints.append(v_prime[nc_count] <= 0)
                    elif opt_prob.nonconvex_constraints[nc]["type"] == "==":
                        nonconvex_constraints.append(ncvx_constraint_tay_func(eps, u_adj, x_ref_param, u_ref_param, ncvx_coef, k) == v_prime[nc_count])
                    nc_count = nc_count + 1

            trust_region_constraints.append(ptr_ops["alpha_x"] * cp.norm_inf(eps[:, k]) + ptr_ops["alpha_u"] * cp.norm_inf(delta_U_k) <= eta[k])
            #trust_region_constraints.append(ptr_ops["alpha_x"] * cp.norm_inf(eps[:, k]) + ptr_ops["alpha_u"] * cp.norm_inf(delta_U_k) <= 1e-6)
        
        constraints = dynamics_constraints + \
                        convex_constraints + \
                        nonconvex_constraints + \
                        trust_region_constraints + [
                        opt_prob.initial_bc(eps[:, 0], eps_0_param) == v_0, # Initial condition constraint
                        opt_prob.terminal_bc(eps[:, -1], eps_N_param) == v_N] # Terminal condition constraint
    
        # Define problem
        cvx_prob = cp.Problem(augmented_objective, constraints)

    # Set parameters
    nc_count = 0
    for k in range(N):
        if k < N - 1:
            cvx_prob.param_dict["Ak_" + str(k)].value = opt_prob.disc["A_k"][:, :, k]
            cvx_prob.param_dict["Bk_" + str(k)].value = opt_prob.disc["B_k"][:, :, k]
        
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
    cvx_prob.param_dict["eps_0"].value = np.array(opt_prob.err_func(opt_prob.x_0.flatten(), x_ref[:, 0])).flatten()
    cvx_prob.param_dict["eps_N"].value = np.array(opt_prob.err_func(opt_prob.x_f.flatten(), x_ref[:, -1])).flatten()

    # Solve
    t0 = time.time()
    val = cvx_prob.solve(solver = ptr_ops["solver"], verbose = ptr_ops["solver_verbose"])
    t1 = time.time()
    if ptr_ops["print_solve_time"]:
        solver_name = ptr_ops["solver"]
        print(f"CVXPY {solver_name} Solve time: {1000 * (t1 - t0):.3f} ms with {val:.3f}")

    # Extract Solution
    eps_sol = cvx_prob.var_dict['eps'].value
    delta_u_sol = cvx_prob.var_dict['delta_U'].value
    x_sol = np.zeros_like(x_ref)
    for k in range(N):
        x_sol[:, k] = opt_prob.eps_add(x_ref[:, k], eps_sol[:, k]) # Add change in Lie algebra to reference by composing reference with the exponential of eps

    eta = cvx_prob.var_dict['eta'].value
    V = cvx_prob.var_dict['V'].value
    v_0 = cvx_prob.var_dict['v_0'].value
    v_N = cvx_prob.var_dict['v_N'].value
    if opt_prob.n["ncvx_k"] > 0:
        v_prime = cvx_prob.var_dict['v_prime'].value
    else:
        v_prime = 0

    solve_status = cvx_prob.status

    J = opt_prob.objective(eps_sol, delta_u_sol + u_ref, x_ref, u_ref).value
    J_ref = opt_prob.objective(eps_sol * 0, u_ref, x_ref, u_ref).value
    J_tr = trust_region_cost_func(ptr_ops["w_tr"], eta)[0]
    J_vc = virtual_control_cost_func(ptr_ops["w_vc"], V, v_0, v_N, v_prime).value

    dx = np.sum(np.linalg.norm(eps_sol, 2))
    du = np.sum(np.linalg.norm(delta_u_sol, 2))

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
                'eta_x':ptr_ops["alpha_x"] * np.sum(np.linalg.norm(eps_sol, 2)),
                'eta_u':ptr_ops["alpha_u"] * np.sum(np.linalg.norm(delta_u_sol, 2))}
    
    return x_sol, eps_sol, delta_u_sol, cvx_prob, sol_info

def virtual_control_cost_func(w_vc, V, v_0, v_N, v_prime):
        #J_vc = w_vc * (cp.sum(cp.sum_squares(V)) + cp.sum(cp.sum_squares(v_0)) + cp.sum(cp.sum_squares(v_N)))
        J_vc = w_vc * (cp.sum(cp.norm(V, 1, 0)) + cp.norm(v_0, 1) + cp.norm(v_N, 1) + cp.norm(v_prime, 1))
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