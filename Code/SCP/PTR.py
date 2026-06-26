import numpy as np
from numpy.linalg import norm as norm
import time
from Code.SCP.solve_convex_subproblem_faster import solve_convex_subproblem

def PTR(opt_prob, ptr_ops):
    '''
    Penalized Trust Region Sequential Convex Programming algorithm
    
    Only set up for ZOH discretization right now

    :param opt_prob: DeterministicProblem object that describes optimization problem
    :param ptr_ops: {iter_min, iter_max, solver, defect_tol, delta_tol, w_tr, w_vc, alpha_x, alpha_u}
    '''
    
    ## Initial Discretization
    defect_k = np.zeros((opt_prob.n["x"], opt_prob.n["N"] + 1, ptr_ops["iter_max"] + 1))
    defect_k[:, :, 0] = opt_prob.discretize(opt_prob.guess["x"], opt_prob.guess["u"])
    defect = np.zeros((ptr_ops["iter_max"] + 1,))
    defect[0] = np.sum(np.linalg.norm(defect_k[:, :, 0], 2))
    
    ## Prepare printout
    print(" k |       status      |   vd  |   vs  | vbc_0 | vbc_N |    J    |   J_tr  |   J_vc   |   dJ %  |   dx   |   du   | delta |  dyn  |  eta  | eta_x | eta_u")
    
    x = np.zeros((opt_prob.n["x"], opt_prob.n["N"], ptr_ops["iter_max"] + 1))
    u = np.zeros((opt_prob.n["u"], opt_prob.n["Nu"], ptr_ops["iter_max"] + 1))
    x[:, :, 0] = opt_prob.guess["x"]
    u[:, :, 0] = opt_prob.guess["u"]
    delta = np.zeros((ptr_ops["iter_max"],))
    
    info = [{}]
    info[0]["J"] = opt_prob.objective(x[:, :, 0], u[:, :, 0], x[:, :, 0], u[:, :, 0]).value
    i = 0
    converged = False
    cvx_prob = []
    while i < ptr_ops["iter_min"] or (i < ptr_ops["iter_max"] and not converged):
        t0 = time.time()
        ## Solve Convex Subproblem
        x[:, :, i + 1], u[:, :, i + 1], cvx_prob, info_i = solve_convex_subproblem(opt_prob, cvx_prob, ptr_ops, x[:, :, i], u[:, :, i])
        info.append(info_i)
        
        ## Discretize
        td0 = time.time()
        defect_k[:, :, i + 1] = opt_prob.discretize(x[:, :, i + 1], u[:, :, i + 1])
        td1 = time.time()
        if ptr_ops["print_solve_time"]:
            print(f"Discretization time: {1000 * (td1 - td0):.3f} ms")

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
        print(f"{i:>2} | {info_i['status']:>17} | {norm(info_i['vd']):>5.1g} | {norm(info_i['vs']):>5.1g} | {norm(info_i['vbc_0']):>5.1g} | {norm(info_i['vbc_N']):>5.1g} | {info_i['J']:>5.3g} | {info_i['J_tr']:>5.3g} | {info_i['J_vc']:>5.3g} | {info_i['dJ']:>5.3g} | {np.sum(info_i['dx']):>6.1g} | {np.sum(info_i['du']):>6.1g} | {info_i['delta']:>5.1g} | {dyn_str:>5} | {norm(info_i['eta']):>5.1g} | {norm(info_i['eta_x']):>5.1g} | {norm(info_i['eta_u']):>5.1g}")
    
        ## Iterate
        i = i + 1
        t1 = time.time()
        if ptr_ops["print_solve_time"]:
            print(f"PTR Solve time: {1000 * (t1 - t0):.3f} ms")

    defect_end = defect[i]
    
    ptr_sol = {"x":x, "u":u, "info":info, "converged":converged, "converged_i":i, "delta":delta, "defect":defect, "defect_k":defect_k}
    
    return ptr_sol