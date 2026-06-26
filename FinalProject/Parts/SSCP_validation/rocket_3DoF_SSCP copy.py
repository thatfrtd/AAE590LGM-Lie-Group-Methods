import numpy as np
import cvxpy as cp
import casadi as cas
import pickle
import matplotlib.pyplot as plt
from Code.Helpers import einsum, sqx2inv
from Code.Dynamics.rocket_landing_3DoF_dynamics import rocket_landing_3DoF_dynamics
from Code.Problems.StochasticProblem import StochasticProblem
from Code.Linearization.linearize_constraint import linearize_constraint
from Code.Plotting.plot_convergence import plot_convergence
from Code.Plotting.plot_rocket_landing_3DoF import plot_rocket_landing_3DoF_trajectory, plot_rocket_landing_3DoF_history, comparison_plot_rocket_landing_3DoF_trajectory
from Code.SO2.SO2_maps import so2_exp
from Code.SSCP.SSCP_PTR_og_block_cholesky import SSCP_PTR_block_cholesky
from Code.SPD.LP_maps import interp_LP_vec, geodesic_dist_LP
from Code.SSCP.qr_covariance import triu, tril, tril_to_mat, tril_to_mat_cp
from Code.SSCP.SSCP_PTR_qr_covariance import SSCP_PTR_qr_covariance

## Initialize Problem
# Constants
m_0 = 2100 # [kg]
I_b = 150000*(1e-3)**2 / m_0 # [kg m2]
g = 9.81*1e-3 # [m / s2] Gravity
lever = -3*1e-3 # [m] distance from CoM to rocket engine
Isp = 200 # [s]
alpha = 1 / (g * Isp)
T_max = 3 * g # [N]
T_min = 0.55 * T_max # [N]
gimbal_max_angle = np.deg2rad(6) # [rad]
c = np.array([[g], [lever], [I_b], [alpha]])

# Initialize problem parameters
nx = 7
neps = nx
nu = 2
N = 25
Nu = N - 1 # ZOH
tf = 35
t_k = np.linspace(0, tf, N)
delta_t = t_k[1] - t_k[0]

# Initialize solver parameters
ptr_ops = {"iter_min":1,
           "iter_max":30,
           "defect_tol":1e-4,
           "defect_S_tol":5e-4,
           "delta_tol":1e-2,
           "w_tr":1e-1 * np.ones((1,N)),
           "w_vc":1e2,
           "alpha_x":1,
           "alpha_u":1,
           "alpha_X":1,
           "solver":"Mosek", 
           "ignore_dpp":False,
           "solver_verbose":True,
           "print_solve_time":False}
integration_tolerance = 1e-10
disc_threads = 3 # Number of threads to have CasADi use during discretization - makes difference?? - should also compile this code...

# Initial conditions
r_0 = np.array([0, 4600]).reshape((2, 1)) * 1e-3
v_0 = so2_exp(np.deg2rad(-60)) @ np.array([306, 0]).reshape((2, 1)) * 1e-3
theta_0 = np.array(np.deg2rad(120)).reshape((1, 1))
w_0 = np.array(0).reshape((1, 1))

x_0 = np.vstack((r_0, v_0, theta_0, w_0, 1))

r_std = np.array([[0.0173], [0.0265]])
v_std = np.array([[0.004], [0.054]])
theta_std = np.deg2rad(np.array([[1]]))
w_std = np.deg2rad(np.array([[0.5]]))
m_std = np.array([[0.001]])
P_0 = np.diag(np.concatenate((r_std, v_std, theta_std, w_std, m_std)).flatten()) ** 2
P_0[0, 1] = -(0.0186 ** 2)
P_0[1, 0] = -(0.0186 ** 2)
P_0[2, 3] = -(0.0034 ** 2)
P_0[3, 2] = -(0.0034 ** 2)
S_0 = tril(np.linalg.cholesky(P_0)) / 3

initial_bc = {"x":lambda x, x_ref : x - x_0.flatten(), 
              "S":lambda S : S - tril(np.linalg.cholesky(P_0))}

# Terminal conditions
r_f = np.array([0, 0.030]).reshape((2, 1))
v_f = np.array([0, 0]).reshape((2, 1))
theta_f = np.array(np.pi/2).reshape((1, 1))
w_f = np.array(0).reshape((1, 1))

x_f = np.vstack((r_f, v_f, theta_f, w_f, 1)) # Mass is guess!

r_std = np.array([[0.001], [0.0003]])
v_std = np.array([[0.001], [0.001]])
theta_std = np.deg2rad(np.array([[1]]))
w_std = np.deg2rad(np.array([[0.5]]))
m_std = np.array([[0.3]])
P_N_guess = np.diag(np.concatenate((r_std, v_std, theta_std, w_std, m_std)).flatten()) ** 2 # Mass is guess
S_N_guess = tril(np.linalg.cholesky(P_N_guess))
P_N = np.diag(np.concatenate((r_std, v_std, theta_std, w_std)).flatten()) ** 2
S_N = tril(np.linalg.cholesky(P_N_guess))

terminal_bc = {"x":lambda x, x_ref : cp.concatenate([x[0:-1] - x_f[0:-1].flatten(), [0]]),
               "S":lambda S : (cp.norm2(np.linalg.inv(tril_to_mat(S_N)) @ tril_to_mat_cp(S)) - 1)}

## Create dynamics
cont_dyn = rocket_landing_3DoF_dynamics(g, lever, I_b, alpha, False)

## Specify Constraints
# State convex constraints
glideslope_max_angle = np.deg2rad(65)
glideslope_constraint = {"k":np.arange(N), "type":"<=", "func":lambda x, u, S_x, S_u : cp.norm2(x[0:2]) - x[1] / np.cos(glideslope_max_angle)}
angvel_max = np.deg2rad(20)
angvel_constraint = {"k":np.arange(N), "type":"<=", "func":lambda x, u, S_x, S_u : cp.norm2(x[5]) - angvel_max}

state_convex_constraints = [glideslope_constraint, angvel_constraint]

# Control convex constraints
thrust_max_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u, S_x, S_u : cp.norm2(u) - T_max}
thrust_max_constraint_chance = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u, S_x, S_u : cp.norm2(u) + sqx2inv(0.95, 3) * cp.norm2(S_u) - T_max}
gimbal_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u, S_x, S_u : cp.norm2(u[1]) - u[0] * np.tan(gimbal_max_angle)}
gimbal_constraint_chance = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u, S_x, S_u : (cp.norm2(u[1]) + sqx2inv(0.95, 1) * cp.norm2(S_u[1, :]) - (u[0] - sqx2inv(0.95, 1) * cp.norm2(S_u[0, :])) * np.tan(gimbal_max_angle))}
bad_thrust_min_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u, S_x, S_u : T_min - u[0]}
bad_thrust_min_constraint_chance = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u, S_x, S_u : T_min - u[0] + sqx2inv(0.95, 1) * cp.norm2(S_u[0, :])} # Replace with nonconvex one? need to be able to nonconvex anyway...

control_convex_constraints = [thrust_max_constraint_chance, gimbal_constraint_chance, bad_thrust_min_constraint_chance]

# Combine convex constraints
convex_constraints = state_convex_constraints + control_convex_constraints

# State nonconvex constraints
state_nonconvex_constraints = []

# Control nonconvex constraints
thrust_min_constraint = {"k":np.arange(Nu), "type":"<=", "func":linearize_constraint(lambda x, u : cas.DM([T_min]) - cas.norm_2(u), nx, nu, "u", [0, 1])}

control_nonconvex_constraints = [thrust_min_constraint]

# Combine nonconvex constraints
nonconvex_constraints = state_nonconvex_constraints + control_nonconvex_constraints

## Specify Objective 
#objective = lambda x, u, x_ref, u_ref : cp.sum(-x[6, -1])
objective = lambda x, u, S_x, S_u, x_ref, u_ref : cp.sum(cp.norm2(u[0:2]) + sqx2inv(0.99, 2) * cp.norm2(S_u)) * delta_t

## Generate Initial Guess
# Load saved deterministic solution
with open('rocket_3DoF_sol.pickle', 'rb') as handle:
    sol_deterministic = pickle.load(handle)
# Interpolate covariance
S_guess = interp_LP_vec(S_0, S_N_guess, N)
# Guess zero feedback
L_guess = np.zeros((nu, neps, Nu))
# Combine guesses
guess = {"x":sol_deterministic["x"], "u":sol_deterministic["u"], "S":S_guess, "L":L_guess}

## Create Problem Object
prob_3DoF = StochasticProblem(x_0, x_f, P_0, P_N, t_k, cont_dyn, guess, objective, convex_constraints, nonconvex_constraints, initial_bc, terminal_bc, integration_tolerance, disc_threads)

## Solve Problem
ptr_sol = SSCP_PTR_qr_covariance(prob_3DoF, ptr_ops)
#ptr_sol = SSCP_PTR_block_cholesky(prob_3DoF, ptr_ops)

## Plot Results
# Convergence plot
plot_convergence(ptr_sol)

# Trajectory plot
converged_i = ptr_sol["converged_i"]
x_disc = ptr_sol["x"][:, :, converged_i]
u_disc = ptr_sol["u"][:, :, converged_i]
plot_rocket_landing_3DoF_trajectory(t_k, x_disc, u_disc, glideslope_max_angle, gimbal_max_angle, T_min, T_max, "Discrete")

t_cont = t_k #np.linspace(0, tf, 100)
x_cont, u_cont = prob_3DoF.cont_prop(u_disc, t_cont)
plot_rocket_landing_3DoF_trajectory(t_cont, x_cont, u_cont, glideslope_max_angle, gimbal_max_angle, T_min, T_max, "Continuous")

# State, control histories plots
#plot_rocket_landing_3DoF_history(x_disc, u_disc, gimbal_max_angle)

plt.stairs(np.linalg.norm(u_cont[0:2, 0:-1], 2, 0) / T_max, t_cont)
plt.hlines(T_min / T_max, 0, tf, color = "gray", linestyles = "--")
plt.hlines(T_max / T_max, 0, tf, color = "gray", linestyles = "--")
plt.grid()
plt.title("Throttle Level vs Time")
plt.show()

plt.stairs(np.rad2deg(np.acos(u_cont[0, 0:Nu] / np.linalg.norm(u_cont[0:2, 0:Nu], 2, 0))), t_cont)
plt.hlines(np.rad2deg(gimbal_max_angle), 0, tf, color = "gray", linestyles = "--")
plt.title("Gimbal Satisfaction")
plt.grid()
plt.show()

# Plot iteration trajectories
#comparison_plot_rocket_landing_3DoF_trajectory(x_iters, u_iters, glideslope_max_angle, gimbal_max_angle, T_min, T_max)
