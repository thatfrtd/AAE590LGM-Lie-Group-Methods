import numpy as np
import cvxpy as cp
import casadi as cas
from Code.Dynamics.rocket_landing_6DoF_dynamics import rocket_landing_6DoF_dynamics
from Code.Problems.DeterministicProblem import DeterministicProblem
from Code.Linearization.linearize_constraint import linearize_constraint
from Code.Discretization.discretize_ZOH_manifold import discretize_error_dynamics_ZOH_manifold
from Code.SCP.PTR import PTR
from Code.Plotting.plot_convergence import plot_convergence
from Code.Plotting.plot_rocket_landing_6DoF import plot_rocket_landing_6DoF_trajectory, plot_rocket_landing_6DoF_history, comparison_plot_rocket_landing_6DoF_trajectory
from Code.SO3.SO3_quat_maps import q_exp, quat_rot, q_conj, q_mul, q_log
from Code.Dynamics.FullRightError.full_right_error import full_right_error, full_add, full_right_error_rate, full_right_error_rate_linearized_ref, full_right_error_cas
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
nx = 14 # dimension of parametrization
neps = 13 # degrees of freedom
nu = 4
N = 20
Nu = N - 1 # ZOH
tf = 45
t_k = np.linspace(0, tf, N)
delta_t = t_k[1] - t_k[0]

# Initialize solver parameters
ptr_ops = {"iter_min":4,
           "iter_max":200,
           "defect_tol":5e-3,
           "delta_tol":2e2,
           "w_tr":4e3 * np.ones((1,N)),
           "w_vc":1e4,
           "alpha_x":1,
           "alpha_u":1,
           "solver":"Clarabel", 
           "solver_verbose":False,
           "print_solve_time":True}
integration_tolerance = 1e-12
disc_threads = 2 # Number of threads to have CasADi use during discretization - makes difference?? - should also compile this code...

# Initial conditions
r_0 = np.array([250, -100, 433]).reshape((3, 1)) * 1e-3 # [km]
#r_0 = np.array([0, 0, 433]).reshape((3, 1)) * 1e-3 # [km]
v_0 = np.array([0, 40, -25]).reshape((3, 1)) * 1e-3 * 0.8 # [km / s]
theta_0 = np.deg2rad(np.array([0, -90, 0]))
q_0 = q_exp(theta_0).reshape((4, 1))
w_0 = np.array([0, 0, 0]).reshape((3, 1))


q_0_log = q_log(q_0.flatten())
q_0_ck = q_exp(q_0_log)

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

## Specify Constraints
# Create LoS trigger function
r_LoS_min = 0.2 # [km]
r_LoS_max = 0.3 # [km]
LoS_distance_trigger = lambda x, u : 1e3*-np.minimum(np.linalg.norm(x[0:3]) - r_LoS_max, 0) * -np.minimum(r_LoS_min - np.linalg.norm(x[0:3]), 0)

# State convex constraints
glideslope_max_angle = np.deg2rad(65)
glideslope_constraint = {"k":np.arange(N), "type":"<=", "func":lambda x, u : cp.norm2(x[0:3]) - x[2] / np.cos(glideslope_max_angle)}
angvel_max = np.deg2rad(5)
angvel_constraint = {"k":np.arange(N), "type":"<=", "func":lambda x, u : cp.norm2(x[10:13]) - angvel_max}
angvel_max_STC = np.deg2rad(2)
angvel_constraint_STC = {"k":np.arange(N), "type":"<=", "func":lambda x, u : cp.norm2(x[10:13]) - angvel_max_STC, "trigger":LoS_distance_trigger}
state_convex_constraints = [glideslope_constraint, angvel_constraint, angvel_constraint_STC]

# Control convex constraints
thrust_max_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u : cp.norm2(u[0:3]) - T_max}
gimbal_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u : cp.norm2(u[0:3]) - u[0] / np.cos(gimbal_max_angle)}
max_torque_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u : cp.abs(u[3]) - tau_max}
bad_thrust_min_constraint = {"k":np.arange(Nu), "type":"<=", "func":lambda x, u : T_min - u[0]} # Replace with nonconvex one? need to be able to nonconvex anyway...
control_convex_constraints = [max_torque_constraint, thrust_max_constraint, gimbal_constraint]

# Combine convex constraints
convex_constraints = state_convex_constraints + control_convex_constraints

# State nonconvex constraints
pitch_max = np.deg2rad(65)
pitch_max_constraint = {"k":np.arange(N), "type":"<=", "func":linearize_constraint(lambda x, u : np.cos(pitch_max) - cas.DM([[0, 0, 1]]) @ quat_rot(x[6:10], cas.DM([[1], [0], [0]])), nx, nu, "x", range(nx))}
'''r_pitch_min = 0.0 # [km]
r_pitch_max = 0.2 # [km]
pitch_max_STC = np.deg2rad(3)
pitch_distance_trigger = lambda x, u : 1000*-np.minimum((np.linalg.norm(x[0:3]) - r_pitch_max), 0) * -np.minimum((r_pitch_min - np.linalg.norm(x[0:3])), 0)
pitch_max_constraint_STC = {"k":np.arange(N), "type":"<=", "func":linearize_constraint(lambda x, u : np.cos(pitch_max_STC) - cas.DM([[0, 0, 1]]) @ quat_rot(x[6:10], cas.DM([[1], [0], [0]])), nx, nu, "x", range(nx)), "trigger":pitch_distance_trigger}'''
LoS_max = np.deg2rad(5)
p_b = np.array([0.5, 0, -np.sqrt(3) / 2]) # Direction of sensor in body frame
d_b = np.array([0, 0, 0]) # Position of sensor in body frame
LoS_constraint = {"k":np.arange(N), "type":"<=", "func":linearize_constraint(lambda x, u : (quat_rot(x[6:10], x[0:3]).T @ p_b + cas.norm_2(x[0:3]) * np.cos(LoS_max) + d_b.T @ p_b), nx, nu, "x", range(nx)), "trigger":LoS_distance_trigger}
state_nonconvex_constraints = [pitch_max_constraint,LoS_constraint]

# Control nonconvex constraints
thrust_min_constraint = {"k":np.arange(Nu), "type":"<=", "func":linearize_constraint(lambda x, u : cas.DM([T_min]) - cas.norm_2(u), nx, nu, "u", [0, 1, 2])}

control_nonconvex_constraints = [thrust_min_constraint]

# Combine nonconvex constraints
nonconvex_constraints = state_nonconvex_constraints + control_nonconvex_constraints

## Specify Objective 
objective = lambda x, u, x_ref, u_ref : cp.sum(-x[13, -1]) + cp.sum(cp.abs(u[3, :]))
#objective = lambda x, u, x_ref, u_ref : cp.sum(cp.norm2(u[0:3, :], 0) + cp.abs(u[3, :])) * delta_t

## Generate Initial Guess
x_guess = x_0 + (x_f - x_0) * (t_k) / tf
u_guess = np.array([[g], [-0.0001], [0.0], [tau_max / 100]]) * np.ones((nu, Nu))
guess = {"x":x_guess, "u":u_guess}

## Create Problem Object
prob_6DoF = DeterministicProblem(x_0, x_f, t_k, cont_dyn, guess, objective, convex_constraints, nonconvex_constraints, initial_bc, terminal_bc, integration_tolerance, disc_threads)

# Propagate and plot guess
x_guess_cont, u_guess_cont = prob_6DoF.cont_prop(u_guess, t_k)
#prob_6DoF.x_0 = full_add(prob_6DoF.x_0.reshape((nx,)), np.array([0, 0, 0, 0, 0, 0, 0.02, -0.02, 0.01, 0, 0, 0.002, 0]).reshape((13,)))
#x_0_diff = x_0 - prob_6DoF.x_0.reshape((nx,1))
x_guess_cont2, u_guess_cont2 = prob_6DoF.cont_prop(u_guess/2, t_k)
#x_guess_cont = x_guess
#u_guess_cont = u_guess

eps_guess = np.zeros((neps, N))
x_reconstruct = np.zeros((nx, N))
for k in range(N):
    eps_guess[:, k] = full_right_error(x_guess_cont2[:, k], x_guess_cont[:, k]).flatten()
    x_reconstruct[:, k] = full_add(x_guess_cont[:, k], eps_guess[:, k])

x_err = x_reconstruct - x_guess_cont2

## TEST FULL ERROR RATE CALCULATION
fRerr_func = full_right_error_rate(g, lever, I_b, alpha)

fRrate = np.array(fRerr_func(eps_guess[:, 0], u_guess_cont2[:, 0] - u_guess_cont[:, 0], x_guess_cont2[:, 0], x_guess_cont[:, 0], u_guess_cont[:, 0]))

def integrate_full_error(fd_func, ferr_func, nx, neps, nu, tolerance):

    X = cas.SX.sym('x', nx + nx + neps, 1)
    p = cas.SX.sym('p', 1 + nu + nu)
    tf = p[0]
    u = p[1:(nu + 1)]
    u_ref = p[(nu + 1):]

    x = X[0:nx]
    x_ref = X[nx:(nx * 2)]
    eps = X[(nx * 2):]
    x_dot = fd_func(x, u)
    x_ref_dot = fd_func(x_ref, u_ref)
    eps_dot = ferr_func(eps, u - u_ref, x, x_ref, u_ref)
    X_dot = cas.vertcat(x_dot, x_ref_dot, eps_dot)
    f_func = cas.Function("f_func", [X, p], [X_dot], ["X", "p"], ["X_dot"])
    opts = {'reltol':tolerance, 'abstol':tolerance, 'jit':False}
    ivp = {'x':X, 'p':p, 'ode':f_func(X, p) * tf}
    F = cas.integrator('F', 'cvodes', ivp, opts)

    return F

eps_err = np.zeros((neps, N - 1))
integrate_func = integrate_full_error(cont_dyn["f_func"], fRerr_func, nx, neps, nu, 1e-8)
for k in range(N - 1):
    r = integrate_func(x0 = np.concat((x_guess_cont2[:, k], x_guess_cont[:, k], eps_guess[:, k])), p = np.concat(([delta_t], u_guess_cont2[:, k], u_guess_cont[:, k])))
    eps_guess_int = np.array(r["xf"][(nx * 2):]).flatten()

    eps_guess_kp1 = eps_guess[:, k + 1]
    eps_err[:, k] = eps_guess_int - eps_guess_kp1
    theta_eps = np.rad2deg(eps_guess_kp1[6:9])

    a = 1

plt.plot(t_k[:-1], eps_err[0:3, :].T)
plt.grid()
plt.title("Position error propagation error")
plt.show()

plt.plot(t_k[:-1], eps_err[3:6, :].T)
plt.grid()
plt.title("Velocity error propagation error")
plt.show()

plt.plot(t_k[:-1], (eps_err[6:9, :].T + np.pi) % (2 * np.pi) - np.pi)
plt.grid()
plt.title("Orientation error propagation error")
plt.show()

plt.plot(t_k[:-1], eps_err[9:12, :].T)
plt.grid()
plt.title("Angular velocity error propagation error")
plt.show()

## TEST FULL ERROR RATE CALCULATION LINEARIZED AT REFERENCE
Rerr_cont_dyn = full_right_error_rate_linearized_ref(g, lever, I_b, alpha)
fRerr_ref_func = Rerr_cont_dyn["f_func"]

fRrate = np.array(fRerr_ref_func(eps_guess[:, 0], u_guess_cont2[:, 0] - u_guess_cont[:, 0], x_guess_cont[:, 0], u_guess_cont[:, 0]))

disc_err_dyn = discretize_error_dynamics_ZOH_manifold(Rerr_cont_dyn["A_func"], Rerr_cont_dyn["B_func"], cont_dyn["f_func"], nx, neps, nu, N, full_right_error_cas, integration_tolerance, disc_threads)
A_ks, B_ks, Delta_ks = disc_err_dyn(x_guess_cont, u_guess_cont[:, :Nu], t_k)
        
disc = {"A_k":np.zeros([neps, neps, N - 1]), 
        "B_k":np.zeros([neps, nu, N - 1]), 
        "defects":np.zeros([neps, N - 1])}
for k in range(N - 1):
    disc["A_k"][:, :, k] = np.array(A_ks[:, k]).reshape((neps, neps)).T
    disc["B_k"][:, :, k] = np.array(B_ks[:, k]).reshape((nu, neps)).T
'''
B_k_0 = disc["B_k"][:, :, 0]
B_k_Nm1 = disc["B_k"][:, :, -1]

with np.printoptions(precision=3, suppress=True, linewidth = 200):
    print(B_k_0)
    print(B_k_Nm1)
'''
Delta_k = np.hstack((initial_bc_Lie(np.zeros((neps,)), full_right_error(x_guess_cont[:, 0], x_0.flatten()).flatten()).reshape((neps, 1)), 
                    np.array(Delta_ks), 
                    terminal_bc_Lie(np.zeros((neps,)), full_right_error(x_guess_cont[:, -1], x_f.flatten()).flatten()).value.reshape((neps, 1))))

def integrate_full_error_ref(fd_func, ferr_func, nx, neps, nu, tolerance):

    X = cas.SX.sym('x', nx + neps, 1)
    p = cas.SX.sym('p', 1 + nu + nu)
    tf = p[0]
    delta_u = p[1:(nu + 1)]
    u_ref = p[(nu + 1):]

    x_ref = X[0:nx]
    eps = X[nx:]
    x_ref_dot = fd_func(x_ref, u_ref)
    eps_dot = ferr_func(eps, delta_u, x_ref, u_ref)
    X_dot = cas.vertcat(x_ref_dot, eps_dot)
    f_func = cas.Function("f_func", [X, p], [X_dot], ["X", "p"], ["X_dot"])
    opts = {'reltol':tolerance, 'abstol':tolerance, 'jit':False}
    ivp = {'x':X, 'p':p, 'ode':f_func(X, p) * tf}
    F = cas.integrator('F', 'cvodes', ivp, opts)

    return F

prob_6DoF.discretize(x_guess_cont, u_guess_cont[:, :Nu])
x_guess_disc = prob_6DoF.disc_prop(u_guess_cont2)

eps_guess_int = np.zeros((neps, N))
eps_guess_disc = np.zeros((neps, N))
eps_guess_xdisc = np.zeros((neps, N))
eps_guess_int[:, 0] = eps_guess[:, 0]
eps_guess_disc[:, 0] = eps_guess[:, 0]
eps_guess_xdisc[:, 0] = eps_guess[:, 0]
eps_err = np.zeros((neps, N - 1))
eps_disc_err = np.zeros((neps, N - 1))
eps_xdisc_err = np.zeros((neps, N - 1))
integrate_func = integrate_full_error_ref(cont_dyn["f_func"], fRerr_ref_func, nx, neps, nu, 1e-8)
for k in range(N - 1):
    delta_u_k = u_guess_cont2[:, k] - u_guess_cont[:, k]
    r = integrate_func(x0 = np.concat((x_guess_cont[:, k], eps_guess_int[:, k])), p = np.concat(([delta_t], delta_u_k, u_guess_cont[:, k])))
    eps_guess_int[:, k + 1] = np.array(r["xf"][(nx):]).flatten()
    eps_guess_disc[:, k + 1] = disc["A_k"][:, :, k] @ eps_guess_disc[:, k] + disc["B_k"][:, :, k] @ delta_u_k + Delta_k[:, k + 1].flatten()
    eps_guess_xdisc[:, k + 1] = full_right_error(x_guess_disc[:, k + 1], x_guess_cont[:, k + 1]).flatten()

    eps_guess_kp1 = eps_guess[:, k + 1]
    eps_err[:, k] = eps_guess_int[:, k + 1] - eps_guess_kp1
    eps_disc_err[:, k] = eps_guess_disc[:, k + 1] - eps_guess_kp1
    eps_xdisc_err[:, k] = eps_guess_xdisc[:, k + 1] - eps_guess_kp1
    theta_err = np.rad2deg(eps_err[6:9, k])
    theta_eps = np.rad2deg(eps_guess_kp1[6:9])

    a = 1



plt.plot(t_k[:-1], eps_err[0:3, :].T)
plt.plot(t_k[:-1], eps_disc_err[0:3, :].T, linestyle = "--")
plt.plot(t_k[:-1], eps_xdisc_err[0:3, :].T, linestyle = ":")
plt.grid()
plt.title("Position error propagation error")
plt.show()

plt.plot(t_k[:-1], eps_err[3:6, :].T)
plt.plot(t_k[:-1], eps_disc_err[3:6, :].T, linestyle = "--")
plt.plot(t_k[:-1], eps_xdisc_err[3:6, :].T, linestyle = ":")
plt.grid()
plt.title("Velocity error propagation error")
plt.show()

plt.plot(t_k[:-1], (eps_err[6:9, :].T + np.pi) % (2 * np.pi) - np.pi)
plt.plot(t_k[:-1], (eps_disc_err[6:9, :].T + np.pi) % (2 * np.pi) - np.pi, linestyle = "--")
plt.plot(t_k[:-1], (eps_xdisc_err[6:9, :].T + np.pi) % (2 * np.pi) - np.pi, linestyle = ":")
plt.grid()
plt.title("Orientation error propagation error")
plt.show()

plt.plot(t_k[:-1], eps_err[9:12, :].T)
plt.plot(t_k[:-1], eps_disc_err[9:12, :].T, linestyle = "--")
plt.plot(t_k[:-1], eps_xdisc_err[9:12, :].T, linestyle = ":")
plt.grid()
plt.title("Angular velocity error propagation error")
plt.show()

plt.plot(t_k[:-1], eps_err[12, :].T)
plt.plot(t_k[:-1], eps_disc_err[12, :].T, linestyle = "--")
plt.plot(t_k[:-1], eps_xdisc_err[12, :].T, linestyle = ":")
plt.grid()
plt.title("Mass error propagation error")
plt.show()

#plot_rocket_landing_6DoF_trajectory(t_k, x_guess, u_guess, glideslope_max_angle, gimbal_max_angle, T_min, T_max, "Guess")
#plot_rocket_landing_6DoF_trajectory(t_k, x_guess_cont, u_guess_cont, glideslope_max_angle, gimbal_max_angle, T_min, T_max, "Guess_cont")

#plot_rocket_landing_6DoF_trajectory(t_k, x_reconstruct, u_guess, glideslope_max_angle, gimbal_max_angle, T_min, T_max, "Guess_reconstruct")

## Solve Problem
ptr_sol = PTR(prob_6DoF, ptr_ops)

## Plot Results
# Convergence plot
#plot_convergence(ptr_sol)

# Trajectory plot

converged_i = ptr_sol["converged_i"]
x_disc = ptr_sol["x"][:, :, converged_i]
u_disc = ptr_sol["u"][:, :, converged_i]
#plot_rocket_landing_6DoF_trajectory(t_k, x_disc, u_disc, glideslope_max_angle, gimbal_max_angle, T_min, T_max, "Discrete")

t_cont = t_k #np.linspace(0, tf, 100)
x_cont, u_cont = prob_6DoF.cont_prop(u_disc, t_cont)
plot_rocket_landing_6DoF_trajectory(t_cont, x_cont, u_cont, glideslope_max_angle, gimbal_max_angle, T_min, T_max, p_b, r_LoS_min, r_LoS_max, "Continuous")

u_disc_old = ptr_sol["u"][:, :, converged_i - 1] # Get a trajectory that isn't quite converged
x_cont_old, u_cont_old = prob_6DoF.cont_prop(u_disc_old, t_cont)
plot_rocket_landing_6DoF_trajectory(t_cont, x_cont_old, u_cont_old, glideslope_max_angle, gimbal_max_angle, T_min, T_max, p_b, r_LoS_min, r_LoS_max, "Continuous_old")

eps = np.zeros((neps, N))
x_reconstruct = np.zeros((nx, N))
for k in range(N):
    eps[:, k] = full_right_error(x_cont_old[:, k], x_cont[:, k]).flatten()
    x_reconstruct[:, k] = full_add(x_cont[:, k], eps[:, k])

prob_6DoF.discretize(x_cont, u_cont[:, :-1])
x_disc_old = prob_6DoF.disc_prop(u_cont_old)

plot_rocket_landing_6DoF_trajectory(t_cont, x_reconstruct, u_cont_old, glideslope_max_angle, gimbal_max_angle, T_min, T_max, p_b, r_LoS_min, r_LoS_max, "Reconstruction")

A_ks, B_ks, Delta_ks = disc_err_dyn(x_cont, u_cont[:, :-1], t_k)
        
disc = {"A_k":np.zeros([neps, neps, N - 1]), 
        "B_k":np.zeros([neps, nu, N - 1]), 
        "defects":np.zeros([neps, N - 1])}
for k in range(N - 1):
    disc["A_k"][:, :, k] = np.array(A_ks[:, k]).reshape((neps, neps)).T

Delta_k = np.hstack((initial_bc_Lie(np.zeros((neps,)), full_right_error(x_cont[:, 0], x_0.flatten()).flatten()).reshape((neps, 1)), 
                    np.array(Delta_ks), 
                    terminal_bc_Lie(np.zeros((neps,)), full_right_error(x_cont[:, -1], x_f.flatten()).flatten()).value.reshape((neps, 1))))


with np.printoptions(precision=3, suppress=True, linewidth = 200):
    print(Delta_k)

eps_int = np.zeros((neps, N))
eps_disc = np.zeros((neps, N))
eps_xdisc = np.zeros((neps, N))
eps_int[:, 0] = eps[:, 0]
eps_disc[:, 0] = eps[:, 0]
eps_xdisc[:, 0] = eps[:, 0]
eps_err = np.zeros((neps, N - 1))
eps_disc_err = np.zeros((neps, N - 1))
eps_xdisc_err = np.zeros((neps, N - 1))
for k in range(N - 1):
    delta_u_k = u_cont_old[:, k] - u_cont[:, k]
    r = integrate_func(x0 = np.concat((x_cont[:, k], eps_int[:, k])), p = np.concat(([delta_t], delta_u_k, u_cont[:, k])))
    eps_int[:, k + 1] = np.array(r["xf"][(nx):]).flatten()
    eps_disc[:, k + 1] = disc["A_k"][:, :, k] @ eps_disc[:, k] + disc["B_k"][:, :, k] @ delta_u_k + Delta_k[:, k + 1].flatten()
    eps_xdisc[:, k + 1] = full_right_error(x_disc_old[:, k + 1], x_cont[:, k + 1]).flatten()

    eps_kp1 = eps[:, k + 1]
    eps_err[:, k] = eps_int[:, k + 1] - eps_kp1
    eps_disc_err[:, k] = eps_disc[:, k + 1] - eps_kp1
    eps_xdisc_err[:, k] = eps_xdisc[:, k + 1] - eps_kp1    
    theta_err = np.rad2deg(eps_err[6:9, k])
    theta_eps = np.rad2deg(eps_kp1[6:9])

    a = 1

plt.plot(t_k[:-1], eps_err[0:3, :].T)
plt.plot(t_k[:-1], eps_disc_err[0:3, :].T, linestyle = "--")
plt.plot(t_k[:-1], eps_xdisc_err[0:3, :].T, linestyle = ":")
plt.grid()
plt.title("Position error propagation error")
plt.show()

plt.plot(t_k[:-1], eps_err[3:6, :].T)
plt.plot(t_k[:-1], eps_disc_err[3:6, :].T, linestyle = "--")
plt.plot(t_k[:-1], eps_xdisc_err[3:6, :].T, linestyle = ":")
plt.grid()
plt.title("Velocity error propagation error")
plt.show()

plt.plot(t_k[:-1], eps_err[6:9, :].T)
plt.plot(t_k[:-1], eps_disc_err[6:9, :].T, linestyle = "--")
plt.plot(t_k[:-1], eps_xdisc_err[6:9, :].T, linestyle = ":")
plt.grid()
plt.title("Orientation error propagation error")
plt.show()

plt.plot(t_k[:-1], eps_err[9:12, :].T)
plt.plot(t_k[:-1], eps_disc_err[9:12, :].T, linestyle = "--")
plt.plot(t_k[:-1], eps_xdisc_err[9:12, :].T, linestyle = ":")
plt.grid()
plt.title("Angular velocity error propagation error")
plt.show()

plt.plot(t_k[:-1], eps_err[12, :].T)
plt.plot(t_k[:-1], eps_disc_err[12, :].T, linestyle = "--")
plt.plot(t_k[:-1], eps_xdisc_err[12, :].T, linestyle = ":")
plt.grid()
plt.title("Mass error propagation error")
plt.show()