import numpy as np
from numpy.linalg import norm as norm
import time
from Code.SSCP.solve_og_block_cholesky_stochastic_convex_subproblem import solve_block_cholesky_stochastic_convex_subproblem
from Code.Helpers import einsum

def SSCP_PTR_block_cholesky(opt_prob, ptr_ops):
    '''
    Penalized Trust Region Stochastic Sequential Convex Programming algorithm using block Cholesky formulation
    
    Only set up for ZOH discretization right now

    :param opt_prob: StochasticProblem object that describes optimization problem
    :param ptr_ops: {iter_min, iter_max, solver, defect_tol, delta_tol, w_tr, w_vc, alpha_x, alpha_u}
    '''
    
    ## Initial Discretization
    defect_k = np.zeros((opt_prob.n["x"], opt_prob.n["N"] + 1, ptr_ops["iter_max"] + 1))
    defect_k[:, :, 0] = opt_prob.discretize(opt_prob.guess["x"], opt_prob.guess["u"])
    defect = np.zeros((ptr_ops["iter_max"] + 1,))
    defect[0] = np.sum(np.linalg.norm(defect_k[:, :, 0], 2))
    
    ## Prepare printout
    print(" k |       status      |   vd  |   vs  | vbc_0 | vbc_N |    J    |   J_tr  |   J_vc   |   dJ %  |   dx   |   du   |  dS_x  |  dS_u  | delta |  dyn  |  eta  | eta_x | eta_u")
    
    x = np.zeros((opt_prob.n["x"], opt_prob.n["N"], ptr_ops["iter_max"] + 1))
    u = np.zeros((opt_prob.n["u"], opt_prob.n["Nu"], ptr_ops["iter_max"] + 1))
    x[:, :, 0] = opt_prob.guess["x"]
    u[:, :, 0] = opt_prob.guess["u"]
    delta = np.zeros((ptr_ops["iter_max"],))
    
    S_x = np.zeros((opt_prob.n["eps"], opt_prob.n["eps"], opt_prob.n["N"], ptr_ops["iter_max"] + 1))
    S_u = np.zeros((opt_prob.n["u"], opt_prob.n["eps"], opt_prob.n["N"] - 1, ptr_ops["iter_max"] + 1))

    info = [{}]
    info[0]["J"] = einsum(range(opt_prob.n["Nu"]), lambda k : opt_prob.objective(x[:, k, 0], u[:, k, 0], np.zeros((opt_prob.n["u"], opt_prob.n["eps"])), x[:, k, 0], u[:, k, 0]).value)
    i = 0
    converged = False
    cvx_prob = []
    while i < ptr_ops["iter_min"] or (i < ptr_ops["iter_max"] and not converged):
        t0 = time.time()
        ## Solve Convex Subproblem
        x[:, :, i + 1], u[:, :, i + 1], S_x[:, :, :, i + 1], S_u[:, :, :, i + 1], cvx_prob, info_i = solve_block_cholesky_stochastic_convex_subproblem(opt_prob, cvx_prob, ptr_ops, x[:, :, i], u[:, :, i], S_x[:, :, :, i], S_u[:, :, :, i])
        info.append(info_i)
        
        ## Discretize
        defect_k[:, :, i + 1] = opt_prob.discretize(x[:, :, i + 1], u[:, :, i + 1])
   
        ## Calculate Stopping Conditions
        defect[i + 1] = np.sum(np.linalg.norm(defect_k[:, :, i + 1], 2))
        delta[i] = info_i["delta"]
        
        # Create dynamics satisfaction string
        if defect[i + 1] <= ptr_ops["defect_tol"]:
            dyn_str = "true"
        else:
            dyn_str = f"{defect[i + 1]:>5.1g}"
        
        # Check PTR convergence 
        converged = delta[i] <= ptr_ops["delta_tol"] and defect[i + 1] <= ptr_ops["defect_tol"]

        ## Display Iteration Results 
        print(f"{i:>2} | {info_i['status']:>17} | {norm(info_i['vd']):>5.1g} | {norm(info_i['vs']):>5.1g} | {norm(info_i['vbc_0']):>5.1g} | {norm(info_i['vbc_N']):>5.1g} | {info_i['J']:>5.3g} | {info_i['J_tr']:>5.3g} | {info_i['J_vc']:>5.3g} | {info_i['dJ']:>5.3g} | {np.sum(info_i['dx']):>6.1g} | {np.sum(info_i['du']):>6.1g} | {np.sum(info_i['dS_x']):>6.1g} | {np.sum(info_i['dS_u']):>6.1g} | {info_i['delta']:>5.1g} | {dyn_str:>5} | {norm(info_i['eta']):>5.1g} | {norm(info_i['eta_x']):>5.1g} | {norm(info_i['eta_u']):>5.1g}")
    
        ## Iterate
        i = i + 1
        t1 = time.time()
        if ptr_ops["print_solve_time"]:
            print(f"PTR Solve time: {1000 * (t1 - t0):.3f} ms")
        print(f"Terminal Cov constraint: {opt_prob.terminal_bc["S"](S_x[:, :, -1, i]).value}")

    defect_end = defect[i]

    # Recover K
    K = recover_sqrt_covariances_from_S(opt_prob, S_x[:, :, :, i], S_u[:, :, :, i])
    
    ptr_sol = {"x":x, "u":u, "K":K, "S_x":S_x[:, :, :, i], "S_u":S_u[:, :, :, i], "info":info, "converged":converged, "converged_i":i, "delta":delta, "defect":defect, "defect_k":defect_k}
    
    return ptr_sol

def recover_sqrt_covariances_from_S(opt_prob, S_x, S_u):
    '''
    Recover gains from sqrt covariances
    '''
    neps = opt_prob.n["eps"]
    nu = opt_prob.n["u"]
    N = opt_prob.n["N"]

    K = np.zeros((nu, neps, N - 1))
    for k in range(N - 1):
        K[:, :, k] = S_u[:, :, k] @ np.linalg.pinv(S_x[:, :, k]) # Only if S_x is square
    return K