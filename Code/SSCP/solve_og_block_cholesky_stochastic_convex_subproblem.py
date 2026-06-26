import numpy as np
import cvxpy as cp
import cvxpygen as cpg
import time
from Code.SSCP.qr_covariance import qr_pos_diag, X_Psqrt_noG, qr_derivative, tril, triu, tril_cp, triu_cp, tril_to_mat, tril_to_mat_cp, triu_to_mat
from Code.Helpers import einsum

def solve_block_cholesky_stochastic_convex_subproblem(opt_prob, cvx_prob, ptr_ops, x_ref, u_ref, S_x_ref, S_u_ref):
    '''
    Solve block Cholesky based stochastic convex subproblem for Stochastic Sequential Convex Programming

    No parameter variables

    :param opt_prob: DeterministicProblem object
    :param cvx_prob: CVXPy problem object (can be list)
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
    #print("Creating problem")
    X = cp.Variable((nx, N), name='X')
    U = cp.Variable((nu, N - 1), name='U')
    S_x = cp.Variable((neps, neps, N), name = 'S_x')
    S_u = cp.Variable((nu, neps, N - 1), name = 'S_u')
    eta = cp.Variable((N,), name = 'eta')
    V = cp.Variable((nx, N - 1), name = 'V')
    v_0 = cp.Variable((nx, 1), name = 'v_0')
    v_N = cp.Variable((nx, 1), name = 'v_N')
    if opt_prob.n["ncvx_k"] > 0:
        v_prime = cp.Variable((opt_prob.n["ncvx_k"]), name = "v_prime")

    # Define objective
    obj_sum = 0
    virtual_control_cost = cp.Minimize(virtual_control_cost_func(ptr_ops["w_vc"], V, v_0, v_N, v_prime))
    trust_region_cost = cp.Minimize(trust_region_cost_func(ptr_ops["w_tr"], eta))

    ## Define constraints
    cvx_trigger_param = []
    dynamics_constraints = []
    stoch_constraints = []
    convex_constraints = []
    nonconvex_constraints = []
    trust_region_constraints = []
    S_unc = np.zeros((neps, neps, N))
    S_unc[:, :, 0] = np.linalg.cholesky(opt_prob.P_0)
    nc_count = 0
    for k in range(N):    
        if k < N - 1:
            # Dynamics constraints
            dynamics_constraints.append(X[:, k + 1] == opt_prob.disc["A_k"][:, :, k] @ X[:, k] + opt_prob.disc["B_k"][:, :, k] @ U[:, k] + opt_prob.disc["c_k"][:, k] + V[:, k])
        
            # Covariance dynamics
            dynamics_constraints.append(S_x[:, :, k + 1] == opt_prob.disc["A_err_k"][:, :, k] @ S_x[:, :, k] + opt_prob.disc["B_err_k"][:, :, k] @ S_u[:, :, k])

            S_x_k = S_x[:, :, k]
            S_u_k = S_u[:, :, k]
            U_k = U[:, k]
            u_adj = U

            obj_sum = obj_sum + opt_prob.objective(X[:, k], U_k, S_u_k, x_ref[:, k], u_ref[:, k])
        else:
            U_k = np.zeros((nu, 1))
            u_adj = np.zeros((nu, N))
            S_x_N = S_x_k
            S_u_k = np.zeros((nu, neps))

        # Convex constraints
        for cc in range(opt_prob.n["cvx"]):
            cc_k = opt_prob.convex_constraints[cc]["k"]
            if np.any(np.isin(cc_k, k)):
                cvx_constraint_func = opt_prob.convex_constraints[cc]["func"]   
                if "trigger" in opt_prob.convex_constraints[cc]:
                    trigger = opt_prob.convex_constraints[cc]["trigger"](x_ref[:, k], u_ref_adj[:, k])
                    cvx_trigger_param.append(trigger)
                else:
                    trigger = 1
                if opt_prob.convex_constraints[cc]["type"] == "<=":
                    convex_constraints.append(cvx_constraint_func(X[:, k], U_k, S_x_k, S_u_k) * trigger <= 0)  
                elif opt_prob.convex_constraints[cc]["type"] == "==":
                    convex_constraints.append(cvx_constraint_func(X[:, k], U_k, S_x_k, S_u_k) * trigger == 0)                         

        # Nonconvex constraints
        for nc in range(opt_prob.n["ncvx"]):
            nc_k = opt_prob.nonconvex_constraints[nc]["k"]
            if np.any(np.isin(nc_k, k)):
                ncvx_constraint_coef_func = opt_prob.nonconvex_constraints[nc]["func"]["tay_coef"]
                ncvx_constraint_tay_func = opt_prob.nonconvex_constraints[nc]["func"]["tay"]
                if "trigger" in opt_prob.nonconvex_constraints[nc]:
                    trigger = opt_prob.nonconvex_constraints[nc]["trigger"](x_ref[:, k], u_ref_adj[:, k])
                else:
                    trigger = 1
                ncvx_coef = ncvx_constraint_coef_func(x_ref, u_ref_adj, k) * trigger
                if opt_prob.nonconvex_constraints[nc]["type"] == "<=":
                    nonconvex_constraints.append(ncvx_constraint_tay_func(X, u_adj, x_ref, u_ref_adj, ncvx_coef, k) <= v_prime[nc_count])
                    nonconvex_constraints.append(v_prime[nc_count] <= 0)
                elif opt_prob.nonconvex_constraints[nc]["type"] == "==":
                    nonconvex_constraints.append(ncvx_constraint_tay_func(X, u_adj, x_ref, u_ref_adj, ncvx_coef, k) == v_prime[nc_count])
                nc_count = nc_count + 1

        trust_region_constraints.append(ptr_ops["alpha_x"] * cp.sum_squares(X[:, k] - x_ref[:, k]) + ptr_ops["alpha_u"] * cp.sum_squares(U_k - u_ref_adj[:, k]) <= eta[k])
    
    constraints = dynamics_constraints + \
                    stoch_constraints + \
                    convex_constraints + \
                    nonconvex_constraints + \
                    trust_region_constraints + [
                    opt_prob.initial_bc["x"](X[:, 0], x_ref[:, 0]) == v_0, # Initial condition constraint
                    opt_prob.terminal_bc["x"](X[:, -1], x_ref[:, -1]) == v_N, # Terminal condition constraint
                    opt_prob.initial_bc["S"](S_x[:, :, 0]) == 0, # Initial covariance constraint
                    opt_prob.terminal_bc["S"](S_x_N) <= 0] # Terminal covariance constraint 

    objective = cp.Minimize(obj_sum)# + cp.sum((np.eye(neps - 1, neps) @ cp.norm2(S_x_N, 1)) - np.diag(np.sqrt(opt_prob.P_f))))
    augmented_objective = objective + virtual_control_cost + trust_region_cost

    # Define problem
    cvx_prob = cp.Problem(augmented_objective, constraints)

    # Solve
    t0 = time.time()
    val = cvx_prob.solve(solver = ptr_ops["solver"], verbose = ptr_ops["solver_verbose"], ignore_dpp = True)#, mosek_params = {'MSK_DPAR_INTPNT_CO_TOL_PFEAS': 1e-12,
                                                                                                            #                 'MSK_DPAR_INTPNT_CO_TOL_DFEAS': 1e-12,
                                                                                                            #                 'MSK_DPAR_INTPNT_CO_TOL_INFEAS':1e-12,
                                                                                                            #                 'MSK_DPAR_INTPNT_CO_TOL_PFEAS':1e-12,
                                                                                                            #                 'MSK_DPAR_INTPNT_CO_TOL_REL_GAP':1e-12,
                                                                                                            #                 'MSK_DPAR_INTPNT_QO_TOL_DFEAS':1e-12,
                                                                                                            #                 'MSK_DPAR_INTPNT_QO_TOL_DFEAS':1e-12,
                                                                                                            #                 'MSK_DPAR_INTPNT_TOL_PFEAS':1e-12})
    t1 = time.time()
    if ptr_ops["print_solve_time"]:
        solver_name = ptr_ops["solver"]
        print(f"CVXPY {solver_name} Solve time: {1000 * (t1 - t0):.3f} ms with {val:.3f}")

    # Extract Solution
    x_sol = cvx_prob.var_dict['X'].value
    u_sol = cvx_prob.var_dict['U'].value
    S_x_sol = cvx_prob.var_dict["S_x"].value
    S_u_sol = cvx_prob.var_dict["S_u"].value

    eta = cvx_prob.var_dict['eta'].value
    V = cvx_prob.var_dict['V'].value
    v_0 = cvx_prob.var_dict['v_0'].value
    v_N = cvx_prob.var_dict['v_N'].value
    if opt_prob.n["ncvx_k"] > 0:
        v_prime = cvx_prob.var_dict['v_prime'].value
    else:
        v_prime = 0

    solve_status = cvx_prob.status

    J = einsum(range(opt_prob.n["Nu"]), lambda k : opt_prob.objective(x_sol[:, k], u_sol[:, k], S_u_sol[:, :, k], x_ref[:, k], u_ref[:, k]).value)
    J_ref = einsum(range(opt_prob.n["Nu"]), lambda k :opt_prob.objective(x_ref[:, k], u_ref[:, k], S_u_ref[:, :, k], x_ref[:, k], u_ref[:, k]).value)
    J_tr = trust_region_cost_func(ptr_ops["w_tr"], eta)[0]
    J_vc = virtual_control_cost_func(ptr_ops["w_vc"], V, v_0, v_N, v_prime).value

    dx = np.sum(np.linalg.norm(x_sol - x_ref, 2))
    du = np.sum(np.linalg.norm(u_sol - u_ref, 2))
    dS_x = 0
    dS_u = 0
    for k in range(N - 1):
        dS_x = dS_x + np.linalg.norm(S_x_sol[:, :, k] - S_x_ref[:, :, k], "fro")
        dS_u = dS_u + np.linalg.norm(S_u_sol[:, :, k] - S_u_ref[:, :, k], "fro")

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
                'dS_x':dS_x,
                'dS_u':dS_u,
                'delta':dx + du,
                'eta':np.linalg.norm(eta),
                'eta_x':ptr_ops["alpha_x"] * np.sum(np.linalg.norm(x_sol - x_ref, 2)),
                'eta_u':ptr_ops["alpha_u"] * np.sum(np.linalg.norm(u_sol - u_ref, 2))}
    
    return x_sol, u_sol, S_x_sol, S_u_sol, [], sol_info

def virtual_control_cost_func(w_vc, V, v_0, v_N, v_prime):
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