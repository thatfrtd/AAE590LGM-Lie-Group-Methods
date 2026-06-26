import numpy as np
import matplotlib.pyplot as plt
import time
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_log, se2_wedge, se2_vee, se2_inverse, se2_Ad, se2_ad
from Code.SE2.SE2_integrators import se2_euler_integrator, se2_Lie_group_integrator
from Code.Plotting.covariance_plotting import plot_cov_ellipse, calc_and_plot_cov_ellipse, cov_ellipse_3D
from Code.Plotting.SE2_plotting import plot_scatter_unicycles, plot_scatter_unicycles_with_cov, plot_Lie_algebra_cov_ellipse
from LieGroupMPC.PD_test_second_problem import PD_test_problem
from Code.SE2.SE2_integrators import se2_kinodynamic_integrator

# Second order system on SE2

# Problem 1: The Banana Distribution

## Unicycle SE2 Dynamics
# Xdot = X xi^hat
# xi = [v; 0; omega]
# xidot = euler poincare

## Initialize Problem Parameters
prob = PD_test_problem()
prob.B_d = np.zeros((3, 3, prob.N_nodes))
prob.B_d_inv = np.zeros((3, 3, prob.N_nodes))
for k in range(prob.N_nodes):
    prob.B_d[:, :, k] = np.eye(3)
    prob.B_d_inv[:, :, k] = np.linalg.pinv(prob.B_d[:, :, k])

## Create nominal reference trajectory (figure 8)    
## 1a) Simulate the Unicycle
# Sim
t_ref = prob.t_k
X_ref = prob.simulate_no_fb(np.zeros([3, prob.N_steps]), prob.X_0)

X_bar_nou = prob.simulate(np.zeros([3, prob.N_steps]), lambda xi_e, eta, xi, xi_bar, k : np.zeros([3, 1]), X_ref, prob.X_0, prob.xi_0)

K_p = 3e-2 * np.eye(3)
K_d = 3e-2 * np.eye(3)
PD_func = lambda xi_e, eta, xi, xi_bar, k : prob.PD(xi_e, eta, xi, xi_bar, K_p, K_d)
dyn_inv_func = lambda xi_e, eta, xi, xi_bar, k : prob.log_linear_dynamic_inversion(xi_e, xi_bar, K_p, k)
X_bar_P = prob.simulate(np.zeros([3, prob.N_steps]), PD_func, X_ref, prob.X_0, prob.xi_0)
X_bar_DI = prob.simulate(np.zeros([3, prob.N_steps]), dyn_inv_func, X_ref, prob.X_0, prob.xi_0)

# Monte Carlo
N_runs = 50
start = time.perf_counter()
X_runs, xi_runs = prob.monte_carlo(N_runs, PD_func, X_ref)
#X_runs = prob.monte_carlo_multithreaded(N_runs)
end = time.perf_counter()
print(f"Elapsed: {(end - start) * 1e3:.3f} ms")

## 1b) Plot the Distribution
'''
plot_scatter_unicycles(X_ref, X_runs, False)
plt.title("Reference Trajectory in (x, y) Plane")
plt.gca().set_aspect('equal', adjustable = 'box')
plt.grid()
#plt.savefig("./HW4/Parts/Q1/AAE590LGM_HW4_Q1b_ref.png")
plt.show()

## 1c) Fit a Gaussian

# Compute sample mean and covariance of final positions
# Plot 3 sigma ellipse on top of scatter

plot_scatter_unicycles_with_cov(X_ref, X_runs, False)
plt.title("Trajectory in (x, y) Plane with 3 Sigma Ellipse at T = 35 s")
plt.gca().set_aspect('equal', adjustable = 'box')
plt.grid()
#plt.savefig("./HW4/Parts/Q1/AAE590LGM_HW4_Q1c_gfit.png")
plt.show()


## 1e) Intermediate Snapshots

plot_scatter_unicycles_with_cov(X_ref, X_runs[:, :, :50, :], False)
plt.title("Trajectory in (x, y) Plane with 3 Sigma Ellipse at T = 5 s")
plt.gca().set_aspect('equal', adjustable = 'box')
plt.grid()
#plt.savefig("./HW4/Parts/Q1/AAE590LGM_HW4_Q1e_t5.png")
plt.show()

plot_scatter_unicycles_with_cov(X_ref, X_runs[:, :, :30, :], False)
plt.title("Trajectory in (x, y) Plane with 3 Sigma Ellipse at T = 3 s")
plt.gca().set_aspect('equal', adjustable = 'box')
plt.grid()
#plt.savefig("./HW4/Parts/Q1/AAE590LGM_HW4_Q1e_t15.png")
plt.show()

plot_scatter_unicycles_with_cov(X_ref, X_runs[:, :, :1, :], False)
plt.title("Trajectory in (x, y) Plane with 3 Sigma Ellipse at T = 0 s")
plt.gca().set_aspect('equal', adjustable = 'box')
plt.grid()
#plt.savefig("./HW4/Parts/Q1/AAE590LGM_HW4_Q1e_t30.png")
plt.show()


## 1f) Lie Algebra Gaussian Fit
# Calculate error w.r.t. reference trajectory
error = np.zeros((3, N_runs))
for i in range(N_runs):
    error[:, i] = np.reshape(se2_log(se2_compose(se2_inverse(X_ref[:, :, -1]), X_runs[:, :, -1, i])), (3,))


# Plot e_1-e_2 and e_1-e_3 projections
plt.subplot(1, 2, 1)
e12_mean, e12_std, e12_cov = calc_and_plot_cov_ellipse(error[0:2, :])
plt.scatter(error[0, :], error[1, :])
plt.xlabel("epsilon_1")
plt.ylabel("epsilon_2")
plt.title("e_1-e_2 Error Projection")
plt.grid()

plt.subplot(1, 2, 2)
e13_mean, e13_std, e13_cov = calc_and_plot_cov_ellipse(error[[0, 2], :])
plt.scatter(error[0, :], error[2, :])
plt.xlabel("epsilon_1")
plt.ylabel("epsilon_3")
plt.title("e_1-e_3 Error Projection")
plt.grid()

plt.suptitle("Lie Algebra Error Projections at T = 35 s")
#plt.savefig("./HW4/Parts/Q1/AAE590LGM_HW4_Q1f_liealg.png")
plt.show()
'''
## 1g) Exponetiate the Ellipsoid Boundary
# Propagate covariance under control
def F_func(X_k, xi_k):
    A = se2_Ad(se2_exp(-xi_k * prob.delta_t))
    return A
Phat_k = np.zeros((3, 3, prob.N_nodes))
Phat_k[:, :, 0] = prob.Phat_0[0:3, 0:3]
Q_k = lambda X : prob.delta_t ** 2 * prob.Q_body

for k in range(prob.N_steps):
    F_k = F_func(X_ref[:, :, k], prob.xi_k[:, k])
    Phat_k[:, :, k + 1] = (F_k - K_p * prob.delta_t) @ Phat_k[:, :, k] @ (F_k - K_p * prob.delta_t).T + Q_k(X_ref[:, :, k])

# Plot
for i in np.arange(1, prob.N_nodes, 30):
    error = np.zeros((3, N_runs))
    for r in range(N_runs):
        error[:, r] = np.reshape(se2_log(se2_compose(se2_inverse(X_ref[:, :, i]), X_runs[:, :, i, r])), (3,))
    error_cov = np.cov(error)
    plot_Lie_algebra_cov_ellipse(error_cov, X_ref[:, :, i], "green", "-")
    plot_Lie_algebra_cov_ellipse(Phat_k[:, :, i], X_ref[:, :, i], "blue", "-")
    plot_scatter_unicycles_with_cov(X_ref[:, :, :(i + 1)], X_runs[:, :, :(i + 1), :], False)
plt.title("Trajectory in (x, y) Plane with 3 Sigma Ellipse at T = 35 s")
plt.grid()
plt.gca().set_aspect('equal', adjustable = 'box')
#plt.savefig("./HW4/Parts/Q1/AAE590LGM_HW4_Q1g_t35.png")
plt.show()