import numpy as np
import cvxpy as cp
import casadi as cas
import pickle
import matplotlib.pyplot as plt
from Code.Dynamics.rocket_landing_3DoF_dynamics import rocket_landing_3DoF_dynamics
from Code.Problems.DeterministicProblem import DeterministicProblem
from Code.Linearization.linearize_constraint import linearize_constraint
from Code.SCP.PTR import PTR
from Code.Plotting.plot_convergence import plot_convergence
from Code.Plotting.plot_rocket_landing_3DoF import plot_rocket_landing_3DoF_trajectory, plot_rocket_landing_3DoF_history, comparison_plot_rocket_landing_3DoF_trajectory
from Code.SO2.SO2_maps import so2_exp

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
nu = 2
N = 25
Nu = N - 1 # ZOH
tf = 35
t_k = np.linspace(0, tf, N)
delta_t = t_k[1] - t_k[0]

# Initialize solver parameters
ptr_ops = {"iter_min":1,
           "iter_max":100,
           "defect_tol":1e-4,
           "delta_tol":1e-2,
           "w_tr":1e-1 * np.ones((1,N)),
           "w_vc":1e2,
           "alpha_x":1,
           "alpha_u":1,
           "solver":"Clarabel", 
           "cvxpygen":False,
           "solver_verbose":False,
           "ignore_dpp":False,
           "print_solve_time":True}
integration_tolerance = 1e-10
disc_threads = 3 # Number of threads to have CasADi use during discretization - makes difference?? - should also compile this code...

# Initial conditions
r_0 = np.array([0, 4600]).reshape((2, 1)) * 1e-3
v_0 = so2_exp(np.deg2rad(-60)) @ np.array([306, 0]).reshape((2, 1)) * 1e-3
theta_0 = np.array(np.deg2rad(120)).reshape((1, 1))
w_0 = np.array(0).reshape((1, 1))

x_0 = np.vstack((r_0, v_0, theta_0, w_0, 1))
initial_bc = lambda x, x_ref : cp.concatenate([x - x_0.flatten()])
#initial_bc = lambda x, x_ref : cp.concatenate([x[0:4] - x_0[0:4].flatten(), [0], x[5:7] - x_0[5:7].flatten()])

# Terminal conditions
r_f = np.array([0, 0.03]).reshape((2, 1))
v_f = np.array([0, 0]).reshape((2, 1))
theta_f = np.array(np.pi/2).reshape((1, 1))
w_f = np.array(0).reshape((1, 1))

x_f = np.vstack((r_f, v_f, theta_f, w_f, 1)) # Mass is guess!
terminal_bc = lambda x, x_ref : cp.concatenate([x[0:-1] - x_f[0:-1].flatten(), [0]])

## Create dynamics
cont_dyn = rocket_landing_3DoF_dynamics(g, lever, I_b, alpha, False)

## Specify Constraints
# State convex constraints
glideslope_max_angle = np.deg2rad(65)
glideslope_constraint = {"k":np.arange(N), "type":"<=", "func":lambda x, u : cp.norm2(x[0:2]) - x[1] / np.cos(glideslope_max_angle)}
angvel_max = np.deg2rad(20)
angvel_constraint = {"k":np.arange(N), "type":"<=", "func":lambda x, u : cp.norm2(x[5]) - angvel_max}

state_convex_constraints = [glideslope_constraint, angvel_constraint]

# Control convex constraints
thrust_max_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u : cp.norm2(u) - T_max}
gimbal_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u : cp.norm2(u) - u[0] / np.cos(gimbal_max_angle)}
bad_thrust_min_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u : T_min - u[0]}

control_convex_constraints = [thrust_max_constraint, gimbal_constraint]

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
objective = lambda x, u, x_ref, u_ref : cp.sum(cp.norm2(u[0:2, :], 0)) * delta_t

## Generate Initial Guess
x_guess = x_0 + (x_f - x_0) * (t_k) / tf
u_guess = np.array([[T_min], [0]]) * np.ones((nu, Nu))
guess = {"x":x_guess, "u":u_guess}

#plot_rocket_landing_3DoF_trajectory(t_k, x_guess, u_guess, glideslope_max_angle, gimbal_max_angle, T_min, T_max, "Guess")

## Create Problem Object
prob_3DoF = DeterministicProblem(x_0, x_f, t_k, cont_dyn, guess, objective, convex_constraints, nonconvex_constraints, initial_bc, terminal_bc, integration_tolerance, disc_threads)

## Solve Problem
ptr_sol = PTR(prob_3DoF, ptr_ops)

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

sol = {"x":x_disc, "u":u_disc}

## Save Trajectory
with open('rocket_3DoF_sol.pickle', 'wb') as handle:
    pickle.dump(sol, handle)