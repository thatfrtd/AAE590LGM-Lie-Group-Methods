import numpy as np
import matplotlib.pyplot as plt
import time
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_log, se2_wedge, se2_vee, se2_inverse
from Code.SE2.SE2_integrators import se2_euler_integrator, se2_Lie_group_integrator
from Code.Plotting.covariance_plotting import plot_cov_ellipse, calc_and_plot_cov_ellipse, cov_ellipse_3D
from Code.Plotting.SE2_plotting import plot_scatter_unicycles, plot_scatter_unicycles_with_cov, plot_Lie_algebra_cov_ellipse
from HW4.Parts.AAE590LGM_HW4_Problem import AAE590LGM_HW4_Problem

# Problem 1: The Banana Distribution

## Unicycle SE2 Dynamics
# Xdot = X xi^hat
# xi = [v; 0; omega]

## Initialize Problem Parameters
prob = AAE590LGM_HW4_Problem()

## 1a) Simulate the Unicycle
# Sim
t_ref = prob.t_k
X_ref = prob.simulate(np.zeros((3, prob.N_nodes)))

# Monte Carlo
N_runs = 100
start = time.perf_counter()
X_runs = prob.monte_carlo(N_runs)
#X_runs = prob.monte_carlo_multithreaded(N_runs)
end = time.perf_counter()
print(f"Elapsed: {(end - start) * 1e3:.3f} ms")

## 1b) Plot the Distribution

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

plot_scatter_unicycles_with_cov(X_ref, X_runs[:, :, :150, :], False)
plt.title("Trajectory in (x, y) Plane with 3 Sigma Ellipse at T = 15 s")
plt.gca().set_aspect('equal', adjustable = 'box')
plt.grid()
#plt.savefig("./HW4/Parts/Q1/AAE590LGM_HW4_Q1e_t15.png")
plt.show()

plot_scatter_unicycles_with_cov(X_ref, X_runs[:, :, :300, :], False)
plt.title("Trajectory in (x, y) Plane with 3 Sigma Ellipse at T = 30 s")
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
plt.xlabel("\epsilon_1")
plt.ylabel("\epsilon_2")
plt.title("e_1-e_2 Error Projection")
plt.grid()

plt.subplot(1, 2, 2)
e13_mean, e13_std, e13_cov = calc_and_plot_cov_ellipse(error[[0, 2], :])
plt.scatter(error[0, :], error[2, :])
plt.xlabel("\epsilon_1")
plt.ylabel("\epsilon_3")
plt.title("e_1-e_3 Error Projection")
plt.grid()

plt.suptitle("Lie Algebra Error Projections at T = 35 s")
#plt.savefig("./HW4/Parts/Q1/AAE590LGM_HW4_Q1f_liealg.png")
plt.show()

## 1g) Exponetiate the Ellipsoid Boundary
error_cov = np.cov(error)
plot_Lie_algebra_cov_ellipse(error_cov, X_ref[:, :, -1], "green", "-")
plot_scatter_unicycles_with_cov(X_ref, X_runs, False)
plt.title("Trajectory in (x, y) Plane with 3 Sigma Ellipse at T = 35 s")
plt.grid()
plt.gca().set_aspect('equal', adjustable = 'box')
#plt.savefig("./HW4/Parts/Q1/AAE590LGM_HW4_Q1g_t35.png")
plt.show()
