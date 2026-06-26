import numpy as np
import casadi as cas
import pickle
import time
import matplotlib.pyplot as plt
from Code.Dynamics.rocket_landing_6DoF_dynamics import rocket_landing_6DoF_dynamics
from Code.Problems.DeterministicProblem import DeterministicProblem
from Code.Linearization.linearize_constraint import linearize_constraint
from Code.Discretization.discretize_ZOH_manifold import discretize_error_dynamics_ZOH_manifold
from Code.SCP.PTR import PTR
from Code.Dynamics.FullRightError.full_right_error import full_right_error, full_add, full_right_error_rate, full_right_error_rate_linearized_ref, full_right_error_cas
from Code.Plotting.plot_convergence import plot_convergence
from Code.Plotting.plot_rocket_landing_6DoF import plot_rocket_landing_6DoF_trajectory, plot_rocket_landing_6DoF_history, comparison_plot_rocket_landing_6DoF_trajectory
from Code.SO3.SO3_quat_maps import q_exp, quat_rot
from Code.SE3.SE3_quat_maps import se3_quat_compose, se3_quat_exp
from Code.Plotting.covariance_plotting import plot_cov_ellipse, calc_and_plot_cov_ellipse, cov_ellipse_3D

## Initialize Problem
# Constants
m_0 = 2100 # [kg]
I_b = np.array([[1915], [14000], [14000]])*(1e-3)**2 / m_0 # [kg m2]
g = 9.81*1e-3 # [m / s2] Gravity
lever = np.array([[-0.25], [0], [0]])*1e-3 # [m] distance from CoM to rocket engine
Isp = 225 # [s]
alpha = 1 / (g * Isp)
T_max = 25000 / m_0 * 1e-3 # [N]
T_min = 6000 / m_0 * 1e-3 # [N]
tau_max = 50 / m_0 * 1e-3 # [Nm]
gimbal_max_angle = np.deg2rad(8) # [rad]

# Initialize problem parameters
nx = 14
neps = 13
nu = 4
N = 25
Nu = N - 1 # ZOH
tf = 45
t_k = np.linspace(0, tf, N)
delta_t = t_k[1] - t_k[0]

# Initialize solver parameters
ptr_ops = {"iter_min":3,
           "iter_max":50,
           "defect_tol":5e-4,
           "delta_tol":2e-2,
           "w_tr":1e3 * np.ones((1,N)),
           "w_vc":1e4,
           "alpha_x":1,
           "alpha_u":1,
           "solver":"QOCO", 
           "cvxpygen":False, # Want to but VS Code being stupid with opening Julia kernal ????
           "sparse":False,
           "ignore_dpp":False,
           "solver_verbose":False,
           "print_solve_time":False}
integration_tolerance = 1e-6
disc_threads = 1 # Number of threads to have CasADi use during discretization - makes difference?? - should also compile this code...

# Initial conditions
r_0 = np.array([250, -100, 433]).reshape((3, 1)) * 1e-3 # [km]
#r_0 = np.array([0, 0, 433]).reshape((3, 1)) * 1e-3 # [km]
v_0 = np.array([-5, 40, -25]).reshape((3, 1)) * 1e-3 * 0.8 # [km / s]
theta_0 = np.deg2rad(np.array([0, -70, 0]))
q_0 = q_exp(theta_0).reshape((4, 1))
w_0 = np.array([0.0, 0.0, 0.0]).reshape((3, 1))

x_0 = np.vstack((r_0, v_0, q_0, w_0, 1))
initial_bc = lambda x, x_ref : x - x_0.flatten()

# Terminal conditions
r_f = np.array([0, 0, 30]).reshape((3, 1)) * 1e-3 # [km]
v_f = np.array([0, 0, 0]).reshape((3, 1)) * 1e-3 # [km]
theta_f = np.deg2rad(np.array([0, -90, 0]))
q_f = q_exp(theta_f).reshape((4, 1))
w_f = np.array([0, 0, 0]).reshape((3, 1))

x_f = np.vstack((r_f, v_f, q_f, w_f, 1)) # Mass is guess!

## Create dynamics
cont_dyn = rocket_landing_6DoF_dynamics(g, lever, I_b, alpha, False)
Rerr_cont_dyn = full_right_error_rate_linearized_ref(g, lever, I_b, alpha)

with open('rocket_6DoF_sol.pickle', 'rb') as handle:
    sol = pickle.load(handle)

def integrate(f_func, nx, nu, tolerance):

    x = cas.SX.sym('x', nx, 1)
    u = cas.SX.sym('u', nu, 1)
    tf = cas.SX.sym('tf')
    p = cas.vertcat(tf, u)
    opts = {'reltol':tolerance, 'abstol':tolerance, 'jit':False}
    ivp = {'x':x, 'p':p, 'ode':f_func(x, u) * tf}
    F = cas.integrator('F', 'cvodes', ivp, opts)

    return F

def integrate_traj(f_func, nx, nu, u_k, t_k, tolerance):
    integrate_func = integrate(f_func, nx, nu, tolerance)
    x_0 = cas.SX.sym('x_0', nx, 1)
    x_k = x_0

    for k in range(len(t_k) - 1):
        sol_k = integrate_func(x0 = x_k[:, k], p = np.vstack((np.array((t_k[k + 1] - t_k[k])).reshape((1, 1)), u_k[:, k].reshape((nu, 1)))))
        x_k = cas.horzcat(x_k, sol_k["xf"])
    
    integ_traj_func = cas.Function("integrate_traj",
                                   [x_0], [x_k],
                                   ["x_0"], ["x_k"])
    
    return integ_traj_func

f_func = cont_dyn["f_func"]
integ_traj_func = integrate_traj(f_func, nx, nu, sol["u"], t_k, integration_tolerance)

N_mc = 100
integrate_func_vec = integ_traj_func.map(N_mc, "thread", 4)

disc_err_dyn = discretize_error_dynamics_ZOH_manifold(Rerr_cont_dyn["A_func"], Rerr_cont_dyn["B_func"], cont_dyn["f_func"], nx, neps, nu, N, full_right_error_cas, integration_tolerance, disc_threads)
A_ks, B_ks, Delta_ks = disc_err_dyn(sol["x"], sol["u"], t_k)
disc = {"A_k":np.zeros([neps, neps, N - 1]), 
        "B_k":np.zeros([neps, nu, N - 1]), 
        "defects":np.zeros([neps, N - 1])}
for k in range(N - 1):
    disc["A_k"][:, :, k] = np.array(A_ks[:, k]).reshape((neps, neps)).T
    disc["B_k"][:, :, k] = np.array(B_ks[:, k]).reshape((nu, neps)).T
    
rng = np.random.default_rng(2)
r_std = np.array([[0.001], [0.001], [0.001]])
v_std = np.array([[0.001], [0.001], [0.001]])
theta_std = np.deg2rad(np.array([[1], [1], [1]]))
w_std = np.deg2rad(np.array([[0.01], [0.01], [0.01]]))
m_std = np.array([[0.001]])
P_0 = np.diag(np.concatenate((r_std, v_std, theta_std, w_std, m_std)).flatten() / 3) ** 2
P_k = np.zeros((neps, neps, N))
P_k[:, :, 0] = P_0

x_0s = np.zeros((nx, N_mc))
for i in range(N_mc):
    x_0s[:, i] = full_add(x_0.flatten(), rng.multivariate_normal(np.zeros((neps,)), P_0))

t0 = time.time()
x_mc = np.array(integrate_func_vec(x_0 = x_0s)["x_k"]).reshape((nx, N_mc, N)).transpose((0, 2, 1))
t1 = time.time()
print(f"MC time: {1000 * (t1 - t0):.3f} ms")

# Plot
def cov_ellipse_6D(P):
    P_eigval, P_eigvec = np.linalg.eig(P)
    a1 = np.reshape(np.linspace(0, np.pi, 5), (5,))
    a2 = np.reshape(np.linspace(0, np.pi, 5), (5,))
    a3 = np.reshape(np.linspace(0, np.pi, 5), (5,))
    a4 = np.reshape(np.linspace(0, np.pi, 5), (5,))
    a5 = np.reshape(np.linspace(0, 2 * np.pi, 10), (10,))

    ellipse_3sigma_data = np.zeros((6, a1.size, a2.size, a3.size, a4.size, a5.size))
    for i1 in range(a1.size):
        for i2 in range(a2.size):
            for i3 in range(a3.size):
                for i4 in range(a4.size):
                    for i5 in range(a5.size):
                        ellipse_3sigma_data[:, i1, i2, i3, i4, i5] = np.reshape(P_eigvec @ np.reshape(np.vstack((3 * np.sqrt(P_eigval[0]) * np.sin(a1[i1]) * np.sin(a2[i2]) * np.sin(a3[i3]) * np.sin(a4[i4]) * np.sin(a5[i5]), 
                                                                                                                 3 * np.sqrt(P_eigval[1]) * np.sin(a1[i1]) * np.sin(a2[i2]) * np.sin(a3[i3]) * np.sin(a4[i4]) * np.cos(a5[i5]),
                                                                                                                 3 * np.sqrt(P_eigval[2]) * np.sin(a1[i1]) * np.sin(a2[i2]) * np.sin(a3[i3]) * np.cos(a4[i4]),
                                                                                                                 3 * np.sqrt(P_eigval[3]) * np.sin(a1[i1]) * np.sin(a2[i2]) * np.cos(a3[i3]),
                                                                                                                 3 * np.sqrt(P_eigval[4]) * np.sin(a1[i1]) * np.cos(a2[i2]), 
                                                                                                                 3 * np.sqrt(P_eigval[5]) * np.cos(a1[i1]))), (6, 1)), (6,))
    return ellipse_3sigma_data.reshape((6, -1))

def plot_Lie_algebra_cov_ellipse(P, Xbar, color, linestyle, ax): # Should change to be more general and go in 'covariance_plotting.py'
    ellipse_3sgima_data = cov_ellipse_6D(P)
    r_boundary = np.zeros((3, np.size(ellipse_3sgima_data, 1)))
    for i in range(np.size(ellipse_3sgima_data, 1)):
        b = se3_quat_compose(se3_quat_exp(ellipse_3sgima_data[:, i]), Xbar)
        r_boundary[:, i] = b[0:3, 0]
    ax.scatter(r_boundary[0, :], r_boundary[1, :], r_boundary[2, :], color = color, linestyle = linestyle, alpha = 0.2) 
    
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

for i in np.arange(10):
    error = np.zeros((neps, N_mc))
    for r in range(N_mc):
        error[:, r] = full_right_error(sol["x"][:, i], x_mc[:, i, r]).flatten()
    error_cov = np.cov(error)
    error_cov = error_cov[[0, 1, 2, 6, 7, 8], :]
    error_cov = error_cov[:, [0, 1, 2, 6, 7, 8]]
    P_k[:, :, i + 1] = disc["A_k"][:, :, i] @ P_k[:, :, i] @ disc["A_k"][:, :, i].T
    #plot_Lie_algebra_cov_ellipse(error_cov, sol["x"][[0, 1, 2, 6, 7, 8, 9], i], "green", "-", ax)
    P_k_se3 = P_k[[0, 1, 2, 6, 7, 8], :, i]
    P_k_se3 = P_k_se3[:, [0, 1, 2, 6, 7, 8]]
    plot_Lie_algebra_cov_ellipse(P_k_se3, sol["x"][[0, 1, 2, 6, 7, 8, 9], i], "blue", "-", ax)
    ax.scatter(x_mc[0, i, :].flatten(), x_mc[1, i, :].flatten(), x_mc[2, i, :].flatten(), alpha = 0.5, marker='o')

#ax.scatter(x_mc[0, :, :].flatten(), x_mc[1, :, :].flatten(), x_mc[2, :, :].flatten(), alpha = 0.5, marker='.')
ax.plot(sol["x"][0, :10], sol["x"][1, :10], sol["x"][2, :10])

ax.set_xlabel("X [km]")
ax.set_ylabel("Y [km]")
ax.set_zlabel("Z [km]")
plt.gca().set_aspect('equal', adjustable = 'box')
plt.show()
