import numpy as np
import cvxpy as cp
import casadi as cas
import mosek
import pickle
import matplotlib.pyplot as plt
from Code.Helpers import einsum, sqx2inv
from Code.Dynamics.rocket_landing_6DoF_dynamics import rocket_landing_6DoF_dynamics
from Code.Problems.StochasticLieProblem import StochasticLieProblem
from Code.Linearization.linearize_constraint import linearize_constraint
from Code.Dynamics.FullRightError.full_right_error import full_right_error, full_add, full_right_error_rate, full_right_error_rate_linearized_ref, full_right_error_cas
from Code.Plotting.plot_convergence import plot_convergence
from Code.Plotting.plot_rocket_landing_6DoF import plot_rocket_landing_6DoF_trajectory, plot_rocket_landing_6DoF_history, comparison_plot_rocket_landing_6DoF_trajectory
from Code.SO3.SO3_quat_maps import q_exp, quat_rot
from Code.SPD.LP_maps import interp_LP_vec, geodesic_dist_LP
from Code.SSCP.qr_covariance import triu, tril, tril_to_mat, tril_to_mat_cp
from Code.SSCP.SSCP_PTR_block_cholesky import SSCP_PTR_block_cholesky

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
N = 50
Nu = N - 1 # ZOH
tf = 35
t_k = np.linspace(0, tf, N)
delta_t = t_k[1] - t_k[0]

# Initialize solver parameters
ptr_ops = {"iter_min":2,
           "iter_max":50,
           "defect_tol":5e-4,
           "delta_tol":2e-2,
           "w_tr":1e2 * np.ones((1,N)),
           "w_vc":1e4,
           "alpha_x":1,
           "alpha_u":0,
           "solver":"Mosek", 
           "solver_verbose":True,
           "print_solve_time":True}
integration_tolerance = 1e-6
disc_threads = 1 # Number of threads to have CasADi use during discretization - makes difference?? - should also compile this code...

# Initial conditions
r_0 = np.array([250, -100, 433]).reshape((3, 1)) * 1e-3 # [km]
#r_0 = np.array([0, 0, 433]).reshape((3, 1)) * 1e-3 # [km]
v_0 = np.array([-5, 40, -25]).reshape((3, 1)) * 1e-3 * 0.8 # [km / s]
#v_0 = np.array([-5, 4, -5]).reshape((3, 1)) * 1e-3 * 0.8 # [km / s]
theta_0 = np.deg2rad(np.array([0, -70, 0]))
#theta_0 = np.deg2rad(np.array([0, -90, 0]))
q_0 = q_exp(theta_0).reshape((4, 1))
w_0 = np.array([0.0, 0.0, 0.0]).reshape((3, 1))

x_0 = np.vstack((r_0, v_0, q_0, w_0, 1))

r_std = np.array([[0.01], [0.01], [0.02]])
v_std = np.array([[0.005], [0.005], [0.015]])
theta_std = np.deg2rad(np.array([[6], [6], [6]]))
w_std = np.deg2rad(np.array([[3], [3], [3]]))
m_std = np.array([[0.001]])
P_0 = np.diag(np.concatenate((r_std, v_std, theta_std, w_std, m_std)).flatten() / 3) ** 2
S_0 = tril(np.linalg.cholesky(P_0))

initial_bc = {"x":lambda x, x_ref : x - x_0.flatten(), 
              "S":lambda S : S - S_0.flatten()}

# Terminal conditions
r_f = np.array([0, 0, 30]).reshape((3, 1)) * 1e-3 # [km]
v_f = np.array([0, 0, 0]).reshape((3, 1)) * 1e-3 # [km]
theta_f = np.deg2rad(np.array([0, -90, 0]))
q_f = q_exp(theta_f).reshape((4, 1))
w_f = np.array([0, 0, 0]).reshape((3, 1))

x_f = np.vstack((r_f, v_f, q_f, w_f, 1)) # Mass is guess!

r_std = np.array([[0.009], [0.009], [0.003]])
v_std = np.array([[0.003], [0.003], [0.003]])
theta_std = np.deg2rad(np.array([[5], [5], [5]]))
w_std = np.deg2rad(np.array([[1.5], [1.5], [1.5]]))
m_std = np.array([[0.10]])
P_N_guess = np.diag(np.concatenate((r_std, v_std, theta_std, w_std, m_std)).flatten() / 3) ** 2 # Mass is guess
S_N_guess = tril(np.linalg.cholesky(P_N_guess))
P_N = np.diag(np.concatenate((r_std, v_std, theta_std, w_std)).flatten() / 3) ** 2
S_N = tril(np.linalg.cholesky(P_N_guess))

# terminal_bc = {"x":lambda x, x_ref : cp.concatenate([x[0:-1] - x_f[0:-1].flatten(), [0]]),
#                "S":lambda S : cp.norm2(np.linalg.inv(tril_to_mat(S_N)) @ np.eye(neps - 1, neps) @ S @ np.eye(neps, neps - 1)) - 1}

terminal_bc = {"x":lambda x, x_ref : cp.concatenate([x[0:-1] - x_f[0:-1].flatten(), [0]]),
               "S":lambda S : (cp.norm2(np.linalg.inv(tril_to_mat(S_N)) @ S) - 1)}

## Create dynamics
cont_dyn = rocket_landing_6DoF_dynamics(g, lever, I_b, alpha, False)
Rerr_cont_dyn = full_right_error_rate_linearized_ref(g, lever, I_b, alpha)

## Specify Constraints
# GRADIENT OF CONSTRAINT w.r.t. LIE ALGEBRA PERTURBATION GIVES (approx) WORST CASE PERTURBATION
# - Should do full derivative by using reference S/L

# Create LoS trigger function
r_LoS_min = 0.2 # [km]
r_LoS_max = 0.3 # [km]
LoS_distance_trigger = lambda x, u : 1e3*-np.minimum(np.linalg.norm(x[0:3]) - r_LoS_max, 0) * -np.minimum(r_LoS_min - np.linalg.norm(x[0:3]), 0)

# State convex constraints
glideslope_max_angle = np.deg2rad(65)
glideslope_constraint = {"k":np.arange(N), "type":"<=", "func":lambda x, u, S, L : cp.norm2(x[0:2]) - x[2] * np.tan(glideslope_max_angle)}
glideslope_constraint_chance = {"k":np.arange(1, N - 1), "type":"<=", "func":lambda x, u, S, L : cp.norm2(x[0:2]) + sqx2inv(0.95, 2) * cp.norm2(S[0:2, :]) - (x[2] - sqx2inv(0.95, 1) * cp.norm2(S[2, :])) * np.tan(glideslope_max_angle)}
angvel_max = np.deg2rad(6)
angvel_constraint = {"k":np.arange(N), "type":"<=", "func":lambda x, u, S, L : cp.norm2(x[10:13]) - angvel_max}
angvel_constraint_chance = {"k":np.arange(1, N - 1), "type":"<=", "func":lambda x, u, S, L : cp.norm2(x[10:13]) + sqx2inv(0.95, 3) * cp.norm2(S[9:12, :]) - angvel_max}
angvel_max_STC = np.deg2rad(2)
angvel_constraint_STC = {"k":np.arange(N), "type":"<=", "func":lambda x, u, S, L : cp.norm2(x[10:13]) - angvel_max_STC, "trigger":LoS_distance_trigger}
angvel_constraint_STC_chance = {"k":np.arange(1, N - 1), "type":"<=", "func":lambda x, u, S, L : cp.norm2(x[10:13]) + sqx2inv(0.95, 3) * cp.norm2(S[9:12, :]) - angvel_max_STC, "trigger":LoS_distance_trigger}
max_r_std = np.array([[0.1], [0.1], [0.1]])
max_v_std = np.array([[0.03], [0.03], [0.03]])
max_theta_std = np.array(np.deg2rad([[30], [30], [30]]))
max_w_std = np.array(np.deg2rad([[10], [10], [10]]))
max_x_std = np.concatenate((max_r_std, max_v_std, max_theta_std, max_w_std)).flatten() / 3
S_x_max_constraint_row = {"k":np.arange(1, N - 1), "type":"<=", "func":lambda x, u, S, L : cp.norm2(S[:-1, :], 1) - max_x_std}
S_x_max_constraint = {"k":np.arange(1, N - 1), "type":"<=", "func":lambda x, u, S, L : cp.norm2(cp.diag(cp.vstack([1 / max_x_std.reshape((-1, 1)), [[0]]])) @ S) - 1}
state_convex_constraints = [glideslope_constraint, angvel_constraint_chance, angvel_constraint_STC]

# Control convex constraints
thrust_max_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u, S, L : cp.norm2(u[0:3]) - T_max}
gimbal_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u, S, L : 1e3*(cp.norm2(u[0:3]) - u[0] / np.cos(gimbal_max_angle))}
max_torque_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u, S, L : cp.abs(u[3]) - tau_max}
bad_thrust_min_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u, S, L : T_min - u[0]} # Replace with nonconvex one? need to be able to nonconvex anyway...
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
state_nonconvex_constraints = [pitch_max_constraint, LoS_constraint]

# Control nonconvex constraints
thrust_min_constraint = {"k":np.arange(Nu), "type":"<=", "func":linearize_constraint(lambda x, u : cas.DM([T_min]) - cas.norm_2(u), nx, nu, "u", [0, 1, 2])}

control_nonconvex_constraints = [thrust_min_constraint]

# Combine nonconvex constraints
nonconvex_constraints = state_nonconvex_constraints + control_nonconvex_constraints

## Specify Objective 
#objective = lambda x, u, S_u, x_ref, u_ref : (cp.norm2(u[0:3]) + sqx2inv(0.99, 3) * cp.norm2(S_u[0:3, :]) + cp.abs(u[3]) + sqx2inv(0.99, 1) * cp.norm2(S_u[3, :])) * delta_t
objective = lambda x, u, S_u, x_ref, u_ref : (cp.norm2(u[0:3]) + cp.abs(u[3])) * delta_t

## Generate Initial Guess
# Load saved deterministic solution
with open('rocket_6DoF_sol.pickle', 'rb') as handle:
    sol_deterministic = pickle.load(handle)
# Interpolate covariance
S_guess = interp_LP_vec(S_0, S_N_guess, N)
# Guess zero feedback
L_guess = np.zeros((nu, neps, Nu))
# Combine guesses
guess = {"x":sol_deterministic["x"], "u":sol_deterministic["u"], "S":S_guess, "L":L_guess}

## Create Problem Object
prob_6DoF = StochasticLieProblem(x_0, x_f, P_0, P_N, t_k, neps, full_right_error_cas, full_add, cont_dyn, Rerr_cont_dyn, guess, objective, convex_constraints, nonconvex_constraints, initial_bc, terminal_bc, integration_tolerance, disc_threads)

## Solve Problem
ptr_sol = SSCP_PTR_block_cholesky(prob_6DoF, ptr_ops)

# with np.printoptions(precision=4, suppress=True, linewidth = 200):
#     for k in range(Nu):
#         print(ptr_sol["K"][:, :, k, ptr_sol["converged_i"]])
#     m = np.max(ptr_sol["K"][:, :, :, ptr_sol["converged_i"]])
#     print(f"Max K val{m}")"

P_N_sol = ptr_sol["S_x"][:, :, N - 1] @ ptr_sol["S_x"][:, :, N - 1].T
final_eig_ck = np.linalg.eigvals(P_N_sol - P_N_guess)
print(f"Final Eig ck: {final_eig_ck}")
'''
x_cov_mag = np.zeros((N,))
x_cov_trace = np.zeros((N,))
u_cov_mag = np.zeros((Nu,))
u_cov_trace = np.zeros((Nu,))
for k in range(N):
    x_cov_mag[k] = np.linalg.norm(ptr_sol["S_x"][:, :, k], 'nuc')
    x_cov_trace[k] = np.linalg.trace(np.sqrt(np.abs(ptr_sol["S_x"][:, :, k] @ ptr_sol["S_x"][:, :, k].T)))
    if k < Nu:
        u_cov_mag[k] = np.linalg.norm(ptr_sol["BK"][:, :, k], 'nuc')
        u_cov_trace[k] = np.linalg.trace(np.sqrt(np.abs(ptr_sol["BK"][:, :, k] @ ptr_sol["BK"][:, :, k].T)))

plt.subplot(1, 2, 1)
plt.plot(x_cov_mag, label = "Mag")
plt.plot(x_cov_trace, label = "Trace")
plt.legend()
plt.grid()
plt.ylabel("State Covariance Mag")

plt.subplot(1, 2, 2)
plt.plot(u_cov_mag, label = "Mag")
plt.plot(u_cov_trace, label = "Trace")
plt.legend()
plt.grid()
plt.ylabel("Control Covariance Mag")

plt.show()
'''
state_stds_normalized = np.zeros((neps - 1, N))
for k in range(N):
    state_stds_normalized[:, k] = np.sqrt(np.diag(ptr_sol["S_x"][:-1, :, k] @ ptr_sol["S_x"][:-1, :, k].T)) / max_x_std#np.diag(np.sqrt(P_N_guess))

sl = cp.norm2(ptr_sol["S_x"][:, :, -1], 1).value - np.diag(np.sqrt(P_N_guess))
sl_t = state_stds_normalized[:, -1]
print(f"Std N ck:{sl}")
print(f"Std N normalized ck:{np.sqrt(np.diag(ptr_sol["S_x"][:-1, :, -1] @ ptr_sol["S_x"][:-1, :, -1].T)) / np.diag(np.sqrt(P_N))}")
print(f"Terminal Cov constraint: {(cp.norm2(np.linalg.inv(tril_to_mat(S_N)) @ ptr_sol["S_x"][:, :, N - 1]).value - 1)}")

plt.plot(state_stds_normalized.T, label = ["r_x", "r_y", "r_z", "v_x", "v_y", "v_z", "theta_x", "theta_y", "theta_z", "w_x", "w_y", "w_z"])#, "m"])
plt.legend()
plt.grid()
plt.ylabel("State Normalized Standard Deviations")
plt.show()
''''''
## Plot Results
# Convergence plot
#plot_convergence(ptr_sol)

# Trajectory plot
converged_i = ptr_sol["converged_i"]
x_disc = ptr_sol["x"][:, :, converged_i]
u_disc = ptr_sol["u"][:, :, converged_i]
plot_rocket_landing_6DoF_trajectory(t_k, x_disc, u_disc, glideslope_max_angle, gimbal_max_angle, T_min, T_max, p_b, r_LoS_min, r_LoS_max, "Discrete")

t_cont = t_k #np.linspace(0, tf, 100)
x_cont = x_disc
u_cont = u_disc
#x_cont, u_cont = prob_6DoF.cont_prop(u_disc, t_cont)
#plot_rocket_landing_6DoF_trajectory(t_cont, x_cont, u_cont, glideslope_max_angle, gimbal_max_angle, T_min, T_max, p_b, r_LoS_min, r_LoS_max, "Continuous")

plt.plot(t_cont, np.rad2deg(np.acos(x_cont[2, :] / np.linalg.norm(x_cont[0:3, :], 2, 0))))
plt.hlines(np.rad2deg(glideslope_max_angle), 0, tf, color = "gray", linestyles = "--")
plt.title("Glideslope Satisfaction")
plt.show()

plt.stairs(np.linalg.norm(u_cont[0:3, 0:Nu], 2, 0) / T_max, t_cont)
plt.hlines(T_min / T_max, 0, tf, color = "gray", linestyles = "--")
plt.hlines(T_max / T_max, 0, tf, color = "gray", linestyles = "--")
plt.grid()
plt.title("Throttle Level vs Time")
plt.show()

plt.stairs(np.rad2deg(np.acos(u_cont[0, 0:Nu] / np.linalg.norm(u_cont[0:3, 0:Nu], 2, 0))), t_cont)
plt.hlines(np.rad2deg(gimbal_max_angle), 0, tf, color = "gray", linestyles = "--")
plt.title("Gimbal Satisfaction")
plt.grid()
plt.show()

plt.stairs(u_cont[3, 0:Nu] * m_0 * 1e3, t_cont)
plt.hlines(tau_max * m_0 * 1e3, 0, tf, color = "gray", linestyles = "--")
plt.hlines(-tau_max * m_0 * 1e3, 0, tf, color = "gray", linestyles = "--")
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
plt.xlabel("Sensor LoS Angle [deg]")
plt.ylabel("Distance from Landing Site [m]")
plt.grid()
plt.show()

plt.plot(np.rad2deg(np.linalg.norm(x_disc[10:13, :], 2, 0)), np.linalg.norm(x_disc[0:3, :], 2, 0)*1e3)
plt.hlines(r_LoS_min*1e3, 0, np.rad2deg(angvel_max)*1.5, color = "gray", linestyles = "--")
plt.hlines(r_LoS_max*1e3, 0, np.rad2deg(angvel_max)*1.5, color = "gray", linestyles = "--")
plt.vlines(np.rad2deg(angvel_max_STC), r_LoS_min*1e3, r_LoS_max*1e3, color = "gray", linestyles = "--")
plt.vlines(np.rad2deg(angvel_max), 0, np.linalg.norm(x_disc[0:3, 0], 2, 0)*1e3, color = "gray", linestyles = "--")
plt.title("State Triggered Angular Velocity Constraint")
plt.xlabel("Angular Velocity [deg / s]")
plt.ylabel("Distance from Landing Site [m]")
plt.grid()
plt.show()

# State, control histories plots
#plot_rocket_landing_6DoF_history(x_disc, u_disc, gimbal_max_angle)

# Plot iteration trajectories
#comparison_plot_rocket_landing_6DoF_trajectory(x_iters, u_iters, glideslope_max_angle, gimbal_max_angle, T_min, T_max)


#sol = {"x":x_disc, "u":u_disc}

## Save Trajectory
#with open('rocket_6DoF_sol.pickle', 'wb') as handle:
#    pickle.dump(sol, handle)