import numpy as np
import cvxpy as cp
import cvxpygen
from scipy.sparse import coo_array
import time

def solve_convex_subproblem(opt_prob, cvx_prob, ptr_ops, x_ref, u_ref):
    '''
    Solve convex subproblem for Sequential Convex Programming

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
        X = cp.Variable((nx, N), name='X')
        U = cp.Variable((nu, N - 1), name='U')
        eta = cp.Variable((N,), name = 'eta')
        V = cp.Variable((nx, N - 1), name = 'V')
        v_0 = cp.Variable((nx, 1), name = 'v_0')
        v_N = cp.Variable((nx, 1), name = 'v_N')
        if opt_prob.n["ncvx_k"] > 0:
            v_prime = cp.Variable((opt_prob.n["ncvx_k"]), name = "v_prime")
        else:
            v_prime = 0

        # Define parameters
        x_ref_param = cp.Parameter((nx, N), name = 'x_ref')
        u_ref_param = cp.Parameter((nu, N), name = 'u_ref')

        # Define objective
        objective = cp.Minimize(opt_prob.objective(X, U, x_ref_param, u_ref_param))
        virtual_control_cost = cp.Minimize(virtual_control_cost_func(ptr_ops["w_vc"], V, v_0, v_N, v_prime))
        trust_region_cost = cp.Minimize(trust_region_cost_func(ptr_ops["w_tr"], eta))
        augmented_objective = objective + virtual_control_cost + trust_region_cost

        ## Define constraints
        A_param = cp.Parameter((opt_prob.cont_dyn["A_sp"][0].size, N - 1), name = 'Ak')
        B_param = cp.Parameter((nx * nu, N - 1), name = 'Bk')
        ck_param = cp.Parameter((nx, N - 1), name = 'ck')
        cvx_trigger_param = []
        dynamics_constraints = []
        convex_constraints = []
        nonconvex_constraints = []
        trust_region_constraints = []
        nc_count = 0
        for k in range(N):
            if k < N - 1:
                # Dynamics constraints
                Bk_param = B_param[:, k].reshape((nx, nu), 'C')
                dynamics_constraints.append(X[:, k + 1] == manual_sparse_multiply(A_param[:, k], X[:, k], opt_prob.cont_dyn["A_sp"]) + Bk_param @ U[:, k]  + ck_param[:, k] + V[:, k])
            
                U_k = U[:, k]
                u_adj = U
            else:
                U_k = np.zeros((nu, 1))
                u_adj = np.zeros((nu, N))

            # Convex constraints
            for cc in range(opt_prob.n["cvx"]):
                cc_k = opt_prob.convex_constraints[cc]["k"]
                if np.any(np.isin(cc_k, k)):
                    cvx_constraint_func = opt_prob.convex_constraints[cc]["func"]   
                    if "trigger" in opt_prob.convex_constraints[cc]:
                        trigger = cp.Parameter((1,), name = f"cvx_trigger_{cc}_{k}", nonneg = True)
                        cvx_trigger_param.append(trigger)
                    else:
                        trigger = 1
                    if opt_prob.convex_constraints[cc]["type"] == "<=":
                        convex_constraints.append(cvx_constraint_func(X[:, k], U_k) * trigger <= 0)  
                    elif opt_prob.convex_constraints[cc]["type"] == "==":
                        convex_constraints.append(cvx_constraint_func(X[:, k], U_k) * trigger == 0)                         

            # Nonconvex constraints
            for nc in range(opt_prob.n["ncvx"]):
                nc_k = opt_prob.nonconvex_constraints[nc]["k"]
                if np.any(np.isin(nc_k, k)):
                    ncvx_constraint_coef_func = opt_prob.nonconvex_constraints[nc]["func"]["tay_coef"]
                    ncvx_constraint_tay_func = opt_prob.nonconvex_constraints[nc]["func"]["tay"]
                    coef_ex = ncvx_constraint_coef_func(x_ref, u_ref_adj, k)
                    ncvx_coef = cp.Parameter(coef_ex.shape, "ncvx_constraint_coef_" + str(nc_count))
                    if opt_prob.nonconvex_constraints[nc]["type"] == "<=":
                        nonconvex_constraints.append(ncvx_constraint_tay_func(X, u_adj, x_ref_param, u_ref_param, ncvx_coef, k) <= v_prime[nc_count])
                        nonconvex_constraints.append(v_prime[nc_count] <= 0)
                    elif opt_prob.nonconvex_constraints[nc]["type"] == "==":
                        nonconvex_constraints.append(ncvx_constraint_tay_func(X, u_adj, x_ref_param, u_ref_param, ncvx_coef, k) == v_prime[nc_count])
                    nc_count = nc_count + 1
         
            trust_region_constraints.append(ptr_ops["alpha_x"] * cp.sum_squares(X[:, k] - x_ref_param[:, k], 0) + ptr_ops["alpha_u"] * cp.sum_squares(U_k - u_ref_param[:, k], 0) <= eta[k])
        
        constraints = dynamics_constraints + \
                        convex_constraints + \
                        nonconvex_constraints + \
                        trust_region_constraints + [
                        opt_prob.initial_bc(X[:, 0], x_ref_param[:, 0]) == v_0, # Initial condition constraint
                        opt_prob.terminal_bc(X[:, -1], x_ref_param[:, -1]) == v_N] # Terminal condition constraint
    
        # Define problem
        cvx_prob = cp.Problem(augmented_objective, constraints)

        if ptr_ops["cvxpygen"]: # VS code is exploding and tries to start Julia for cpg....... ahhhhhh
            cvxpygen.cpg.generate_code(cvx_prob, code_dir='Deterministic_6DoF_fixed_FOH_QOCO', solver = ptr_ops["solver"])
            from Deterministic_6DoF_fixed_FOH_QOCO.cpg_solver import cpg_solve
            cvx_prob.register_solve(ptr_ops["solver"] + "_cpg", cpg_solve)
            ptr_ops["solver"] = ptr_ops["solver"] + "_cpg"

    # Set parameters
    A_k_sparse = np.zeros((opt_prob.cont_dyn["A_sp"][0].size, N - 1))
    nc_count = 0
    for k in range(N):
        if k < N - 1:
            A_k_sparse[:, k] = opt_prob.disc["A_k"][:, :, k][opt_prob.cont_dyn["A_sp"]]
        
        for cc in range(opt_prob.n["cvx"]):
            cc_k = opt_prob.convex_constraints[cc]["k"]
            if np.any(np.isin(cc_k, k)):
                cvx_constraint_func = opt_prob.convex_constraints[cc]["func"]   
                if "trigger" in opt_prob.convex_constraints[cc]:
                    cvx_prob.param_dict[f"cvx_trigger_{cc}_{k}"].value = np.array([opt_prob.convex_constraints[cc]["trigger"](x_ref[:, k], u_ref_adj[:, k])])
                    
        for nc in range(opt_prob.n["ncvx"]):
            nc_k = opt_prob.nonconvex_constraints[nc]["k"]
            if np.any(np.isin(nc_k, k)):
                ncvx_constraint_coef_func = opt_prob.nonconvex_constraints[nc]["func"]["tay_coef"]
                ncoef = ncvx_constraint_coef_func(x_ref, u_ref_adj, k)
                if "trigger" in opt_prob.nonconvex_constraints[nc]:
                    trigger = opt_prob.nonconvex_constraints[nc]["trigger"](x_ref[:, k], u_ref_adj[:, k])
                    if trigger == 0:
                        ncoef = np.zeros_like(ncoef)
                else:
                    trigger = 1
                cvx_prob.param_dict["ncvx_constraint_coef_" + str(nc_count)].value = ncoef * trigger
                nc_count = nc_count + 1

    cvx_prob.param_dict["Ak"].value = A_k_sparse
    cvx_prob.param_dict["Bk"].value = opt_prob.disc["B_k"].reshape((nx * nu, N - 1))
    cvx_prob.param_dict["ck"].value = opt_prob.disc["c_k"]

    cvx_prob.param_dict["x_ref"].value = x_ref
    cvx_prob.param_dict["u_ref"].value = u_ref_adj

    # Solve
    t0 = time.time()
    val = cvx_prob.solve(solver = ptr_ops["solver"], verbose = ptr_ops["solver_verbose"], ignore_dpp = ptr_ops["ignore_dpp"])
    t1 = time.time()
    if ptr_ops["print_solve_time"]:
        solver_name = ptr_ops["solver"]
        print(f"CVXPY {solver_name} Solve time: {1000 * (t1 - t0):.3f} ms with {val:.3f}")

    # Extract Solution
    x_sol = cvx_prob.var_dict['X'].value
    u_sol = cvx_prob.var_dict['U'].value

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

def manual_sparse_multiply(sp_data, vec, sp_shape):
    
    result = 0
    ind_0 = np.flatnonzero(sp_shape[0] == 0)
    for i in ind_0:
        result = result + sp_data[i] * vec[sp_shape[1][i]]
    result = result.reshape(1,'C')

    for j in range(1, np.max(sp_shape[0]) + 1):
        result_j = 0
        ind_j = np.flatnonzero(sp_shape[0] == j)
        for i in ind_j:
            result_j = result_j + sp_data[i] * vec[sp_shape[1][i]]
        result = cp.concatenate([result, result_j.reshape(1,'C')])
    
    return result