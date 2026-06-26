import numpy as np
import casadi as cas
import scipy as sp
from Code.Discretization.discretize_ZOH_manifold import discretize_error_dynamics_ZOH_manifold

class DeterministicLieProblem:
    n = {} # x, u, eps, N, Nu, cvx, ncvx, ncvx_k
    x_0 = []
    x_f = []
    t_k = []
    tf = []
    err_func = [] # Function to compute error state
    eps_add = [] # Function for the combined operation of composing group elememnts with the exponential of Lie algebra elements
    cont_dyn = {} # f_func, A_func, B_func
    cont_manifold_dyn = {} # f_func, A_func, B_func
    disc_dyn = {} # disc_dyn(x, u, t) -> (A_ks, B_ks, c_ks, defect_ks)
    disc = {}
    guess = {} # x, u
    objective = []
    convex_constraints = {}
    nonconvex_constraints = {}
    initial_bc = []
    terminal_bc = []
    integration_tolerance = []

    def __init__(self, x_0, x_f, t_k, neps, cont_manifold_dyn, err_func, eps_add, cont_dyn, guess, objective, convex_constraints, nonconvex_constraints, initial_bc, terminal_bc, integration_tolerance, disc_threads):
        self.x_0 = x_0
        self.x_f = x_f
        self.t_k = t_k
        self.tf = t_k[-1]
        self.err_func = err_func
        self.eps_add = eps_add
        self.cont_manifold_dyn = cont_manifold_dyn
        self.cont_dyn = cont_dyn
        self.guess = guess
        nc_count = 0
        for nc in range(len(nonconvex_constraints)):
            nc_count = nc_count + len(nonconvex_constraints[nc]["k"])
        self.n = {"x":guess["x"].shape[0], "eps":neps, "u":guess["u"].shape[0], "N":guess["x"].shape[1], "Nu":guess["u"].shape[1], "cvx":len(convex_constraints), "ncvx":len(nonconvex_constraints), "ncvx_k":nc_count}
        self.disc_dyn = discretize_error_dynamics_ZOH_manifold(cont_dyn["A_func"], cont_dyn["B_func"], cont_manifold_dyn["f_func"], self.n["x"], self.n["eps"], self.n["u"], self.n["N"], err_func, integration_tolerance, disc_threads)
        self.disc = {"A_k":np.zeros([self.n["eps"], self.n["eps"], self.n["N"] - 1]), 
                     "B_k":np.zeros([self.n["eps"], self.n["u"], self.n["N"] - 1]), 
                     "c_k":np.zeros([self.n["eps"], self.n["N"] - 1]), 
                     "defects":np.zeros([self.n["x"], self.n["N"] - 1])}
        self.objective = objective
        self.convex_constraints = convex_constraints
        self.nonconvex_constraints = nonconvex_constraints
        self.initial_bc = initial_bc
        self.terminal_bc = terminal_bc
        self.integration_tolerance = integration_tolerance
    
    def discretize(self, x_k, u_k):
        A_ks, B_ks, Delta_ks = self.disc_dyn(x_k, u_k, self.t_k)
        
        for k in range(self.n["N"] - 1):
            self.disc["A_k"][:, :, k] = np.array(A_ks[:, k]).reshape((self.n["eps"], self.n["eps"])).T
            self.disc["B_k"][:, :, k] = np.array(B_ks[:, k]).reshape((self.n["u"], self.n["eps"])).T
        self.disc["c_k"] = np.array(Delta_ks)
    
        Delta_k = np.hstack((self.initial_bc(np.zeros((self.n["eps"],)), np.array(self.err_func(x_k[:, 0], self.x_0.flatten())).flatten()).reshape((self.n["eps"], 1)), 
                             np.array(Delta_ks), 
                             self.terminal_bc(np.zeros((self.n["eps"],)), np.array(self.err_func(x_k[:, -1], self.x_f.flatten())).flatten()).value.reshape((self.n["eps"], 1))))

        return Delta_k
    
    def cont_prop(self, u, t_eval):
        # Replace with casadi integrator for faster performance (or generate code from CasADi)
        
        x_sol = np.zeros((self.n["x"], t_eval.size))
        x_sol[:, 0] = self.x_0.flatten()
        u_sol = np.zeros((self.n["u"], t_eval.size))
        for k in range(self.n["N"] - 1):
            t_span_k = self.t_k[k:(k+2)]
            k_eval = np.where(np.logical_and(t_eval >= t_span_k[0], t_eval <= t_span_k[1]))[0]
            sol = sp.integrate.solve_ivp(lambda t, x : np.array(self.cont_manifold_dyn["f_func"](x, u[:, k])).flatten(), t_span_k, x_sol[:, k_eval[0]], t_eval = t_eval[k_eval], rtol = self.integration_tolerance, atol = self.integration_tolerance)
            x_sol[:, k_eval] = sol.y
            u_sol[:, k_eval] = np.repeat(u[:, k].reshape((-1, 1)), len(k_eval), 1)
        
        return x_sol, u_sol

    