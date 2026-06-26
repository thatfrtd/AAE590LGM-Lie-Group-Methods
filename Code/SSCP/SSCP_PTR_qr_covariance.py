import numpy as np
from numpy.linalg import norm as norm
import time
from Code.SSCP.solve_qr_covariance_stochastic_convex_subproblem import solve_qr_covariance_stochastic_convex_subproblem
from Code.SSCP.qr_covariance import calculate_X_S_defect_noG, tril_to_mat

def SSCP_PTR_qr_covariance(opt_prob, ptr_ops):
    '''
    Penalized Trust Region Stochastic Sequential Convex Programming algorithm using qr covariance formulation
    
    Only set up for ZOH discretization right now

    :param opt_prob: StochasticProblem object that describes optimization problem
    :param ptr_ops: {iter_min, iter_max, solver, defect_tol, delta_tol, w_tr, w_vc, alpha_x, alpha_u}
    '''
    
    ## Initial Discretization
    defect_k = np.zeros((opt_prob.n["x"], opt_prob.n["N"] + 1, ptr_ops["iter_max"] + 1))
    defect_k[:, :, 0] = opt_prob.discretize(opt_prob.guess["x"], opt_prob.guess["u"])
    defect = np.zeros((ptr_ops["iter_max"] + 1,))
    defect[0] = np.sum(np.linalg.norm(defect_k[:, :, 0], 2))
    defect_S_k = np.zeros((opt_prob.n["N"] + 1, ptr_ops["iter_max"] + 1))
    X, defect_S_k[:, 0] = calculate_X_S_defect_noG_with_bc(opt_prob.disc["A_err_k"], opt_prob.disc["B_err_k"], opt_prob.guess["S"], opt_prob.guess["L"], opt_prob.initial_bc["S"], opt_prob.terminal_bc["S"])
    defect_S = np.zeros((ptr_ops["iter_max"] + 1,))
    defect_S[0] = np.sum(defect_S_k[:, 0])

    ## Prepare printout
    print(" k |       status      |   vd  |  vd_S |   vs  | vbc_0 | vbc_N |    J    |   J_tr  |   J_vc   |   dJ %  |   dx   |   du   |   dX   | delta |  dyn  |  dyn_X |  eta  | eta_x | eta_u | eta_X ")
    
    # Set up arrays and initialize to initial guess
    x = np.zeros((opt_prob.n["x"], opt_prob.n["N"], ptr_ops["iter_max"] + 1))
    u = np.zeros((opt_prob.n["u"], opt_prob.n["Nu"], ptr_ops["iter_max"] + 1))
    x[:, :, 0] = opt_prob.guess["x"]
    u[:, :, 0] = opt_prob.guess["u"]
    delta = np.zeros((ptr_ops["iter_max"],))

    S = np.zeros((opt_prob.n["eps"] * int((opt_prob.n["eps"] + 1) / 2), opt_prob.n["N"], ptr_ops["iter_max"] + 1))
    L = np.zeros((opt_prob.n["u"], opt_prob.n["eps"], opt_prob.n["N"] - 1, ptr_ops["iter_max"] + 1))
    S[:, :, 0] = opt_prob.guess["S"]
    L[:, :, :, 0] = opt_prob.guess["L"]

    ## PTR
    info = [{}]
    info[0]["J"] = opt_prob.objective(x[:, :, 0], u[:, :, 0], S[:, :, 0], L[:, :, :, 0], x[:, :, 0], u[:, :, 0]).value
    i = 0
    converged = False
    cvx_prob = []
    while i < ptr_ops["iter_min"] or (i < ptr_ops["iter_max"] and not converged):
        t0 = time.time()
        ## Solve Convex Subproblem
        x[:, :, i + 1], u[:, :, i + 1], S[:, :, i + 1], L[:, :, :, i + 1], info_i = solve_qr_covariance_stochastic_convex_subproblem(opt_prob, ptr_ops, x[:, :, i], u[:, :, i], S[:, :, i], L[:, :, :, i], X)
        info.append(info_i)
        
        ## Discretize
        defect_k[:, :, i + 1] = opt_prob.discretize(x[:, :, i + 1], u[:, :, i + 1])
        # Calculate covariance defect
        X, defect_S_k[:, i + 1] = calculate_X_S_defect_noG_with_bc(opt_prob.disc["A_err_k"], opt_prob.disc["B_err_k"], S[:, :, i + 1], L[:, :, :, i + 1], opt_prob.initial_bc["S"], opt_prob.terminal_bc["S"])

        ## Calculate Stopping Conditions
        defect[i + 1] = np.sum(np.linalg.norm(defect_k[:, :, i + 1], 2, 0))
        defect_S[i + 1] = np.sum(defect_S_k[:-1, i + 1])
        delta[i] = info_i["delta"]
        
        # Create mean dynamics satisfaction string
        if defect[i + 1] <= ptr_ops["defect_tol"]:
            dyn_str = "true"
        else:
            dyn_str = f"{defect[i + 1]:>5.1g}"

        # Create covariance dynamics satisfaction string
        if defect_S[i + 1] <= ptr_ops["defect_S_tol"]:
            dyn_S_str = "true"
        else:
            dyn_S_str = f"{defect_S[i + 1]:>5.1g}"
        
        # Check PTR convergence 
        converged = delta[i] <= ptr_ops["delta_tol"] and defect_S[i + 1] <= ptr_ops["defect_S_tol"] and defect[i + 1] <= ptr_ops["defect_tol"]

        ## Display Iteration Results 
        print(f"{i:>2} | {info_i['status']:>17} | {norm(info_i['vd']):>5.1g} | {norm(info_i['vd_S']):>5.1g} | {norm(info_i['vs']):>5.1g} | {norm(info_i['vbc_0']):>5.1g} | {norm(info_i['vbc_N']):>5.1g} | {info_i['J']:>5.3g} | {info_i['J_tr']:>5.3g} | {info_i['J_vc']:>5.3g} | {info_i['dJ']:>5.3g} | {np.sum(info_i['dx']):>6.1g} | {np.sum(info_i['du']):>6.1g} | {np.sum(info_i['dX']):>6.1g} | {info_i['delta']:>5.1g} | {dyn_str:>5} | {dyn_S_str:>5} | {norm(info_i['eta']):>5.1g} | {norm(info_i['eta_x']):>5.1g} | {norm(info_i['eta_u']):>5.1g} | {norm(info_i['eta_X']):>5.1g}")
    
        ## Iterate
        i = i + 1
        t1 = time.time()
        if ptr_ops["print_solve_time"]:
            print(f"PTR Solve time: {1000 * (t1 - t0):.3f} ms")

    defect_end = defect[i]

    K = recover_gain_matrices(S, L)
    P = recover_covariance_matrices(S)
    
    ptr_sol = {"x":x, "u":u, "K":K, "P":P, "S":S, "L":L, "X":X, "info":info, "converged":converged, "converged_i":i, "delta":delta, "defect":defect, "defect_k":defect_k, "defect_S_k":defect_S_k}
    
    return ptr_sol

def calculate_X_S_defect_noG_with_bc(A_err_k, B_err_k, S_k, L_k, ic, tc):
    defect_S_k = np.zeros(S_k.shape[1] + 1)
    defect_S_k[0] = np.linalg.norm(ic(S_k[:, 0]))
    X, defect_S_k[1:-1] = calculate_X_S_defect_noG(A_err_k, B_err_k, S_k, L_k)
    a = tc(S_k[:, -1]).value
    defect_S_k[-1] = np.max([tc(S_k[:, -1]).value, 0]) # tc should return scalar

    return X, defect_S_k

def recover_gain_matrices(S, L):
    K = np.zeros_like(L.shape[1], L.shape[0], L.shape[2])
    for k in range(L.shape[2]):
        K[:, :, k] = L[:, :, k] @ np.linalg.inv(tril_to_mat(S[:, :, k]))

    return K

def recover_covariance_matrices(S):
    P = np.zeros_like(S)
    for k in range(S.shape[2]):
        P[:, :, k] = S[:, :, k] @ S[:, :, k].T

    return P