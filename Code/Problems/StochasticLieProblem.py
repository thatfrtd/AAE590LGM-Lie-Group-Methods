import numpy as np
import casadi as cas
import scipy as sp
from Code.Discretization.discretize_ZOH import discretize_error_dynamics_ZOH
from Code.Discretization.discretize_ZOH_manifold import discretize_error_dynamics_ZOH_manifold

class StochasticLieProblem:
    n = {} # x, eps, u, N, Nu, cvx, ncvx, ncvx_k
    x_0 = []
    x_f = []
    P_0 = []
    P_f = []
    t_k = []
    tf = []
    err_func = [] # Function to compute error state
    eps_add = [] # Function for the combined operation of composing group elememnts with the exponential of Lie algebra elements
    cont_dyn = {} # f_func, A_func, B_func
    cont_err_dyn = {} # f_func, A_func, B_func
    disc_dyn = {} # disc_dyn(x, u, t) -> (A_ks, B_ks, c_ks, defect_ks)
    disc_err_dyn = {} # disc_dyn(x, u, t) -> (A_ks, B_ks, defect_ks)
    disc = {}
    guess = {} # x, u
    objective = []
    convex_constraints = {}
    nonconvex_constraints = {}
    initial_bc = [] # x for state, S for covariance
    terminal_bc = [] # x for state, S for covariance
    integration_tolerance = []

    def __init__(self, x_0, x_f, P_0, P_f, t_k, neps, err_func, eps_add, cont_dyn, cont_err_dyn, guess, objective, convex_constraints, nonconvex_constraints, initial_bc, terminal_bc, integration_tolerance, disc_threads):
        self.x_0 = x_0
        self.x_f = x_f
        self.P_0 = P_0
        self.P_f = P_f
        self.t_k = t_k
        self.tf = t_k[-1]
        self.err_func = err_func
        self.eps_add = eps_add
        self.cont_dyn = cont_dyn
        self.cont_err_dyn = cont_err_dyn
        self.guess = guess
        nc_count = 0
        for nc in range(len(nonconvex_constraints)):
            nc_count = nc_count + len(nonconvex_constraints[nc]["k"])
        self.n = {"x":guess["x"].shape[0], "eps":neps, "u":guess["u"].shape[0], "N":guess["x"].shape[1], "Nu":guess["u"].shape[1], "cvx":len(convex_constraints), "ncvx":len(nonconvex_constraints), "ncvx_k":nc_count}
        self.disc_dyn = discretize_error_dynamics_ZOH(cont_dyn["A_func"], cont_dyn["B_func"], cont_dyn["f_func"], self.n["x"], self.n["u"], self.n["N"], integration_tolerance, disc_threads)
        self.disc = {"A_k":np.zeros([self.n["x"], self.n["x"], self.n["N"] - 1]), 
                     "B_k":np.zeros([self.n["x"], self.n["u"], self.n["N"] - 1]), 
                     "c_k":np.zeros([self.n["x"], self.n["N"] - 1]), 
                     "defects":np.zeros([self.n["x"], self.n["N"] - 1]),
                     "A_err_k":np.zeros([self.n["eps"], self.n["eps"], self.n["N"] - 1]), 
                     "B_err_k":np.zeros([self.n["eps"], self.n["u"], self.n["N"] - 1])}
        self.disc_err_dyn = discretize_error_dynamics_ZOH_manifold(cont_err_dyn["A_func"], cont_err_dyn["B_func"], cont_dyn["f_func"], self.n["x"], self.n["eps"], self.n["u"], self.n["N"], err_func, integration_tolerance, disc_threads)
        self.objective = objective
        self.convex_constraints = convex_constraints
        self.nonconvex_constraints = nonconvex_constraints
        self.initial_bc = initial_bc
        self.terminal_bc = terminal_bc
        self.integration_tolerance = integration_tolerance
    
    def discretize(self, x_k, u_k):
        A_ks, B_ks, c_ks, Delta_ks = self.disc_dyn(x_k, u_k, self.t_k)
        A_err_ks, B_err_ks, Delta_err_ks = self.disc_err_dyn(x_k, u_k, self.t_k)

        for k in range(self.n["N"] - 1):
            self.disc["A_k"][:, :, k] = np.array(A_ks[:, k]).reshape((self.n["x"], self.n["x"])).T
            self.disc["B_k"][:, :, k] = np.array(B_ks[:, k]).reshape((self.n["u"], self.n["x"])).T
            self.disc["c_k"][:, k] = np.array(c_ks[:, k]).reshape((self.n["x"],))
            self.disc["A_err_k"][:, :, k] = np.array(A_err_ks[:, k]).reshape((self.n["eps"], self.n["eps"])).T
            self.disc["B_err_k"][:, :, k] = np.array(B_err_ks[:, k]).reshape((self.n["u"], self.n["eps"])).T
     
        Delta_k = np.hstack((self.initial_bc["x"](x_k[:, 0], x_k[:, 0]).reshape((self.n["x"], 1)), 
                             np.array(Delta_ks), 
                             self.terminal_bc["x"](x_k[:, -1], x_k[:, -1]).value.reshape((self.n["x"], 1))))

        return Delta_k
    
    def cont_prop(self, x_ref, u_ref, K, N_mc, mc_threads):
        # Monte Carlo propagation of solution to stochastic trajectory optimization problem
        # ADD STOCHASTIC PERTURBATIONS
        # ADD ESTIMATOR
        # Generate code from CasADi?
        
        integ_traj_fb_func = self.integrate_traj_fb(x_ref, u_ref, K, self.t_k, self.integration_tolerance)
        integrate_fb_func_vec = integ_traj_fb_func.map(N_mc, "thread", mc_threads)
       
        rng = np.random.default_rng(2)
        x_0s = np.zeros((self.n["x"], N_mc))
        for i in range(N_mc):
            x_0s[:, i] = self.eps_add(self.x_0.flatten(), rng.multivariate_normal(np.zeros((self.n["eps"],)), self.P_0))

        mc_sol = integrate_fb_func_vec(x_0 = x_0s)
        x_mc = np.array(mc_sol["x_k"]).reshape((self.n["x"], N_mc, self.n["N"])).transpose((0, 2, 1))
        u_mc = np.array(mc_sol["u_k"]).reshape((self.n["u"], N_mc, self.n["N"] - 1)).transpose((0, 2, 1))

        return x_mc, u_mc
    
    def integrate_traj_fb(self, x_ref_k, u_ref_k, K_k, t_k, tolerance):
        integrate_func = cas_integrate(self.cont_dyn["f_func"], self.n["x"], self.n["u"], tolerance)
        x_0 = cas.SX.sym('x_0', self.n["x"], 1)
        x_k = x_0

        u_k = cas.GenDM_zeros((self.n["u"], self.n["Nu"]))
        for k in range(len(t_k) - 1): # Should be smaller steps for stochastic perturbations
            eps = self.err_func(x_k[:, k], x_ref_k[:, k])
            u_k_fb = K_k[:, :, k] @ eps
            u_k[:, k] = u_ref_k[:, k] + u_k_fb
            sol_k = integrate_func(x0 = x_k[:, k], p = cas.vertcat(cas.DM([t_k[k + 1] - t_k[k]]).reshape((1, 1)), u_k[:, k]))
            x_k = cas.horzcat(x_k, sol_k["xf"])
        
        integ_traj_fb_func = cas.Function("integrate_traj_fb",
                                    [x_0], [x_k, u_k],
                                    ["x_0"], ["x_k", "u_k"])
        
        return integ_traj_fb_func
    
def cas_integrate(f_func, nx, nu, tolerance):

    x = cas.SX.sym('x', nx, 1)
    u = cas.SX.sym('u', nu, 1)
    tf = cas.SX.sym('tf')
    p = cas.vertcat(tf, u)
    opts = {'reltol':tolerance, 'abstol':tolerance, 'jit':False}
    ivp = {'x':x, 'p':p, 'ode':f_func(x, u) * tf}
    F = cas.integrator('F', 'cvodes', ivp, opts)

    return F
