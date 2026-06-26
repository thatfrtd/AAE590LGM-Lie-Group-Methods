import numpy as np
import cvxpy as cp
import casadi as cas
import pickle
from Code.Dynamics.rocket_landing_6DoF_dynamics import rocket_landing_6DoF_dynamics
from Code.Problems.DeterministicLieProblem import DeterministicLieProblem
from Code.Linearization.linearize_constraint import linearize_constraint
from Code.iSCP.iPTR import iPTR
from Code.Dynamics.FullRightError.full_right_error import full_right_error, full_add, full_right_error_rate, full_right_error_rate_linearized_ref, full_right_error_cas
from Code.Plotting.plot_convergence import plot_convergence
from Code.Plotting.plot_rocket_landing_6DoF import plot_rocket_landing_6DoF_trajectory, plot_rocket_landing_6DoF_history, comparison_plot_rocket_landing_6DoF_trajectory
from Code.SO3.SO3_quat_maps import q_exp, quat_rot
import matplotlib.pyplot as plt

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
gimbal_max_angle = np.deg2rad(5) # [rad]

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
ptr_ops = {"iter_min":10,
           "iter_max":50,
           "defect_tol":5e-4,
           "delta_tol":2e-2,
           "w_tr":1e3 * np.ones((1,N)),
           "w_vc":1e4,
           "alpha_x":10,
           "alpha_u":10,
           "solver":"QOCO", 
           "solver_verbose":False,
           "print_solve_time":False}
integration_tolerance = 1e-6
disc_threads = 2 # Number of threads to have CasADi use during discretization - makes difference?? - should also compile this code...

# Initial conditions
r_0 = np.array([250, -100, 433]).reshape((3, 1)) * 1e-3 # [km]
#r_0 = np.array([0, 0, 433]).reshape((3, 1)) * 1e-3 # [km]
v_0 = np.array([-5, 40, -25]).reshape((3, 1)) * 1e-3 * 0.8 # [km / s]
theta_0 = np.deg2rad(np.array([0, -70, 0]))
q_0 = q_exp(theta_0).reshape((4, 1))
w_0 = np.array([0, 0, 0]).reshape((3, 1))

x_0 = np.vstack((r_0, v_0, q_0, w_0, 1))
initial_bc = lambda x, x_ref : x - x_0.flatten()
initial_bc_Lie = lambda eps, eps_ref : eps - eps_ref

# Terminal conditions
r_f = np.array([0, 0, 30]).reshape((3, 1)) * 1e-3 # [km]
v_f = np.array([0, 0, 0]).reshape((3, 1)) * 1e-3 # [km]
theta_f = np.deg2rad(np.array([0, -90, 0]))
q_f = q_exp(theta_f).reshape((4, 1))
w_f = np.array([0, 0, 0]).reshape((3, 1))

x_f = np.vstack((r_f, v_f, q_f, w_f, 1)) # Mass is guess!
terminal_bc = lambda x, x_ref : cp.concatenate([x[0:-1] - x_f[0:-1].flatten(), [0]])
terminal_bc_Lie = lambda eps, eps_ref : cp.concatenate([eps[0:-1] - eps_ref[0:-1].flatten(), [0]])

## Create dynamics
cont_dyn = rocket_landing_6DoF_dynamics(g, lever, I_b, alpha, False)
Rerr_cont_dyn = full_right_error_rate_linearized_ref(g, lever, I_b, alpha)

## Specify Constraints
# Create LoS trigger function
r_LoS_min = 0.2 # [km]
r_LoS_max = 0.3 # [km]
LoS_distance_trigger = lambda x, u : 1e3*-np.minimum(np.linalg.norm(x[0:3]) - r_LoS_max, 0) * -np.minimum(r_LoS_min - np.linalg.norm(x[0:3]), 0)

# State convex constraints
glideslope_max_angle = np.deg2rad(65)
glideslope_constraint = {"k":np.arange(N), "type":"<=", "func":lambda eta, u, x_ref : cp.norm2(x_ref[0:3] + eta[0:3]) - (x_ref[2] + eta[2]) / np.cos(glideslope_max_angle)}
mass_constraint = {"k":np.arange(N), "type":"<=", "func":lambda eta, u, x_ref : eta[12] + x_ref[13] - 1}
angvel_max = np.deg2rad(5)
angvel_constraint = {"k":np.arange(N), "type":"<=", "func":lambda eta, u, x_ref : cp.norm2(x_ref[10:13] + eta[9:12]) - angvel_max}
angvel_max_STC = np.deg2rad(2)
angvel_constraint_STC = {"k":np.arange(N), "type":"<=", "func":lambda eta, u, x_ref : cp.norm2(x_ref[10:13] + eta[9:12]) - angvel_max_STC, "trigger":LoS_distance_trigger}
state_convex_constraints = [glideslope_constraint, mass_constraint]

# Control convex constraints
thrust_max_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda eta, u, x_ref : cp.norm2(u[0:3]) - T_max}
gimbal_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda eta, u, x_ref : 1e6*(cp.norm2(u[0:3]) - u[0] / np.cos(gimbal_max_angle))}
max_torque_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda eta, u, x_ref : cp.abs(u[3]) - tau_max}
bad_thrust_min_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda eta, u, x_ref : T_min - u[0]} # Replace with nonconvex one? need to be able to nonconvex anyway...
control_convex_constraints = [max_torque_constraint, thrust_max_constraint, gimbal_constraint]

# Combine convex constraints
convex_constraints = state_convex_constraints + control_convex_constraints

# State nonconvex constraints
pitch_max = np.deg2rad(65)
pitch_max_constraint = {"k":np.arange(N), "type":"<=", "func":linearize_constraint(lambda x, u : np.cos(pitch_max) - cas.DM([[0, 0, 1]]) @ quat_rot(x[6:10], cas.DM([[1], [0], [0]])), nx, nu, "x", range(nx))}
LoS_max = np.deg2rad(5)
p_b = np.array([0.5, 0, -np.sqrt(3) / 2]) # Direction of sensor in body frame
d_b = np.array([0, 0, 0]) # Position of sensor in body frame
LoS_constraint = {"k":np.arange(N), "type":"<=", "func":linearize_constraint(lambda x, u : (quat_rot(x[6:10], x[0:3]).T @ p_b + cas.norm_2(x[0:3]) * np.cos(LoS_max) + d_b.T @ p_b), nx, nu, "x", range(nx)), "trigger":LoS_distance_trigger}
state_nonconvex_constraints = []#[pitch_max_constraint, LoS_constraint]

# Control nonconvex constraints
thrust_min_constraint = {"k":np.arange(Nu), "type":"<=", "func":linearize_constraint(lambda x, u : cas.DM([T_min]) - cas.norm_2(u), nx, nu, "u", [0, 1, 2])}

control_nonconvex_constraints = []#[thrust_min_constraint]

# Combine nonconvex constraints
nonconvex_constraints = state_nonconvex_constraints + control_nonconvex_constraints

## Specify Objective 
#objective = lambda eta, u, x_ref, u_ref : cp.sum(-(x_ref[13, -1] + eta[12, -1])) + cp.sum(cp.abs(u[3, :]))
objective = lambda x, u, x_ref, u_ref : cp.sum(cp.norm2(u[0:3, :], 0) + cp.abs(u[3, :])) * delta_t

## Generate Initial Guess
eps_ic_tc = full_right_error(x_f.flatten(), x_0.flatten())
x_guess = np.zeros((nx, N))
for k in range(N):
    x_guess[:, k] = full_add(x_0.flatten(), eps_ic_tc.flatten() * t_k[k] / tf).flatten()

x_guess_truth = x_0 + (x_f - x_0) * (t_k) / tf

u_guess = np.array([[g], [0], [0], [0]]) * np.ones((nu, Nu))
guess = {"x":x_guess, "u":u_guess}

with open('rocket_6DoF_sol.pickle', 'rb') as handle:
    sol = pickle.load(handle)
guess = sol

## Create Problem Object
prob_6DoF = DeterministicLieProblem(x_0, x_f, t_k, neps, cont_dyn, full_right_error_cas, full_add, Rerr_cont_dyn, guess, objective, convex_constraints, nonconvex_constraints, initial_bc_Lie, terminal_bc_Lie, integration_tolerance, disc_threads)

# Propagate and plot guess
x_guess_cont, u_guess_cont = prob_6DoF.cont_prop(u_guess, t_k)

#plot_rocket_landing_6DoF_trajectory(t_k, x_guess, u_guess, glideslope_max_angle, gimbal_max_angle, T_min, T_max, p_b, r_LoS_min, r_LoS_max, "Guess_Lie")
#plot_rocket_landing_6DoF_trajectory(t_k, x_guess_truth, u_guess, glideslope_max_angle, gimbal_max_angle, T_min, T_max, p_b, r_LoS_min, r_LoS_max, "Guess_Truth")
#plot_rocket_landing_6DoF_trajectory(t_k, x_guess_cont, u_guess_cont, glideslope_max_angle, gimbal_max_angle, T_min, T_max, "Guess_cont")

## Solve Problem
ptr_sol = iPTR(prob_6DoF, ptr_ops)

## Plot Results
# Convergence plot
plot_convergence(ptr_sol)

# Trajectory plot
converged_i = ptr_sol["converged_i"]
x_disc = ptr_sol["x"][:, :, converged_i]
u_disc = ptr_sol["u"][:, :, converged_i]
plot_rocket_landing_6DoF_trajectory(t_k, x_disc, u_disc, glideslope_max_angle, gimbal_max_angle, T_min, T_max, p_b, r_LoS_min, r_LoS_max, "Discrete")

t_cont = t_k #np.linspace(0, tf, 100)
x_cont, u_cont = prob_6DoF.cont_prop(u_disc, t_cont)
plot_rocket_landing_6DoF_trajectory(t_cont, x_cont, u_cont, glideslope_max_angle, gimbal_max_angle, T_min, T_max, p_b, r_LoS_min, r_LoS_max, "Continuous")

plt.plot(t_cont, np.rad2deg(np.acos(x_cont[2, :] / np.linalg.norm(x_cont[0:3, :], 2, 0))))
plt.hlines(np.rad2deg(glideslope_max_angle), 0, tf, color = "gray", linestyles = "--")
plt.title("Glideslope Satisfaction")
plt.show()

plt.stairs(np.linalg.norm(u_cont[0:3, 0:-1], 2, 0) / T_max, t_cont)
plt.hlines(T_min / T_max, 0, tf, color = "gray", linestyles = "--")
plt.hlines(T_max / T_max, 0, tf, color = "gray", linestyles = "--")
plt.grid()
plt.title("Throttle Level vs Time")
plt.show()

plt.stairs(np.rad2deg(np.acos(u_cont[0, 0:-1] / np.linalg.norm(u_cont[0:3, 0:-1], 2, 0))), t_cont)
plt.hlines(np.rad2deg(gimbal_max_angle), 0, tf, color = "gray", linestyles = "--")
plt.title("Gimbal Satisfaction")
plt.grid()
plt.show()

plt.stairs(u_cont[3, 0:-1] * m_0 * 1e3, t_cont)
plt.hlines(tau_max * m_0 * 1e3, 0, tf, color = "gray", linestyles = "--")
plt.title("Roll Torque")
plt.grid()
plt.show()

plt.stairs(cp.abs(u_disc[3, :] * m_0 *1e3).value, t_k)
plt.hlines(tau_max * m_0 * 1e3, 0, tf, color = "gray", linestyles = "--")
plt.title("Roll Torque")
plt.grid()
plt.show()

pitch = np.squeeze(np.rad2deg(np.acos(np.array([cas.DM([[0, 0, 1]]) @ quat_rot(x_disc[6:10, k], cas.DM([[1], [0], [0]])) for k in range(Nu)]))))
pitches = np.squeeze(np.array([pitch_max_constraint["func"]["val"](x_disc, u_disc, x_disc, u_disc, k) for k in range(Nu)]))
real_pitches = np.squeeze(np.rad2deg(np.acos([np.array([[0, 0, 1]]) @ quat_rot(x_disc[6:10, k], np.array([[1], [0], [0]])) for k in range(Nu)])))

plt.plot(t_k[:Nu], np.rad2deg(pitches))
plt.plot(t_k[:Nu], real_pitches)
plt.plot(t_k[:Nu], pitch)
plt.plot(t_k, np.linalg.norm(x_disc[0:3], 2, 0))
plt.grid()
plt.title("Pitches")
plt.show()

LoS_vals = np.squeeze(np.rad2deg(np.acos([-(np.array(quat_rot(cas.DM(x_disc[6:10, k]), cas.DM(x_disc[0:3, k]))).T @ p_b + d_b.T @ p_b) / np.linalg.norm(x_disc[0:3, k]) for k in range(N)])))
plt.plot(LoS_vals, np.linalg.norm(x_disc[0:3, :], 2, 0)*1e3)
plt.hlines(r_LoS_min*1e3, 0, 80, color = "gray", linestyles = "--")
plt.hlines(r_LoS_max*1e3, 0, 80, color = "gray", linestyles = "--")
plt.vlines(np.rad2deg(LoS_max), 0, np.linalg.norm(x_disc[0:3, 0], 2, 0)*1e3, color = "gray", linestyles = "--")
plt.title("State Triggered Line of Sight Constraint")
plt.grid()
plt.show()

plt.plot(np.rad2deg(np.linalg.norm(x_disc[10:13, :], 2, 0)), np.linalg.norm(x_disc[0:3, :], 2, 0)*1e3)
plt.hlines(r_LoS_min*1e3, 0, np.rad2deg(angvel_max)*1.5, color = "gray", linestyles = "--")
plt.hlines(r_LoS_max*1e3, 0, np.rad2deg(angvel_max)*1.5, color = "gray", linestyles = "--")
plt.vlines(np.rad2deg(angvel_max_STC), r_LoS_min*1e3, r_LoS_max*1e3, color = "gray", linestyles = "--")
plt.vlines(np.rad2deg(angvel_max), 0, np.linalg.norm(x_disc[0:3, 0], 2, 0)*1e3, color = "gray", linestyles = "--")
plt.title("State Triggered Angular Velocity Constraint")
plt.grid()
plt.show()

# State, control histories plots
#plot_rocket_landing_6DoF_history(x_disc, u_disc, gimbal_max_angle)

# Plot iteration trajectories
#comparison_plot_rocket_landing_6DoF_trajectory(x_iters, u_iters, glideslope_max_angle, gimbal_max_angle, T_min, T_max)
