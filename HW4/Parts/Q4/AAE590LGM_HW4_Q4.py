import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from Code.SO2.SO2_maps import so2_exp
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_log, se2_wedge, se2_vee, se2_inverse, se2_mat_to_tuplevec, se2_tuplevec_to_mat, se2_Ad
from Code.SE2.SE2_integrators import se2_euler_integrator, se2_Lie_group_integrator
from Code.Plotting.covariance_plotting import plot_cov_ellipse, calc_and_plot_cov_ellipse, cov_ellipse_3D, create_cov_ellipse, plot_cov_ellipse_ax
from Code.Plotting.SE2_plotting import plot_scatter_unicycles, plot_scatter_unicycles_with_cov, plot_Lie_algebra_cov_ellipse, plot_scatter_unicycles_ax
from HW4.Parts.AAE590LGM_HW4_Problem import AAE590LGM_HW4_Problem
from Code.Filtering.KalmanFilter import ESKF
from Code.Filtering.filter_performance_metrics import averaged_chi_squared_consistency_bands, NEES

# Problem 4: Invariant EKF

## Unicycle SE2 Dynamics
# Xdot = X xi^hat
# xi = [v; 0; omega]

## Initialize Problem Parameters
# Simulation
prob = AAE590LGM_HW4_Problem()

# Monte Carlo
N_runs = 100

def IEKF_SE2(prob):
    t_ref = prob.t_k
    
    def f_disc_func(X_k, xi_k):
        X_kp1 = se2_compose(X_k, se2_exp(xi_k * prob.delta_t))
        return X_kp1

    def F_func(X_k, xi_k):
        A = se2_Ad(se2_exp(-xi_k * prob.delta_t))
        return A
    
    def x_update_func(X_k_km1, delta_x_k):
        X_k_k = se2_compose(X_k_km1, se2_exp(delta_x_k)) 
        return X_k_k

    # 2a) Prediction Step    
    Q_k = lambda X : prob.delta_t ** 2 * prob.Q_body
    
    # 2b) Measurement Update
    h_GPS = lambda X : X[0:2, 2]
    h_GPS_noise = lambda X : prob.rng.normal(h_GPS(X), prob.sigma_GPS)
    R_z_GPS = prob.sigma_GPS ** 2 * np.eye(2)
    H_GPS_func = lambda X : np.hstack((X[0:2, 0:2], np.zeros((2, 1))))
    
    # 2c) Implementation
    xhat_0 = np.zeros((3, 1))
    P_0 = np.diag((0.01, 0.01, 0.01))
    
    X_true = prob.monte_carlo(1)
    x_true = np.zeros((3, prob.N_nodes))
    for k in range(prob.N_nodes):
        x_true[:, k] = se2_mat_to_tuplevec(X_true[:, :, k]).reshape((3,))
    
    xhat_k = np.zeros((3, prob.N_nodes))
    Xhat_k = np.zeros((3, 3, prob.N_nodes))
    Phat_k = np.zeros((3, 3, prob.N_nodes))
    xhat_k[:, 0] = xhat_0.reshape((3,))
    Xhat_k[:, :, 0] = se2_tuplevec_to_mat(xhat_0).reshape((3, 3))
    Phat_k[:, :, 0] = P_0
    nu_k = []
    K_k = []
    innovation = []
    t_k_GPS = []
    z_k_GPS = []
    k_GPS = []
    for k in range(prob.N_steps):
        if np.mod(t_ref[k], 1 / prob.GPS_rate) == 0 and k > 1:
            z_GPS = h_GPS_noise(X_true[:, :, k])
            z_k_GPS.append(z_GPS)
            t_k_GPS.append(t_ref[k])
            k_GPS.append(k)
    
            Xhat_k[:, :, k + 1], Phat_k[:, :, k + 1], nu_k_i, K_k_i, innovation_i = ESKF(Xhat_k[:, :, k], Phat_k[:, :, k], prob.xi_k[:, k], h_GPS, z_GPS, f_disc_func, F_func, Q_k, H_GPS_func, R_z_GPS, x_update_func)
    
            nu_k.append(nu_k_i)
            K_k.append(K_k_i)
            innovation.append(innovation_i)
        else:
            z_GPS = np.nan    
    
            Xhat_k[:, :, k + 1], Phat_k[:, :, k + 1] = ESKF(Xhat_k[:, :, k], Phat_k[:, :, k], prob.xi_k[:, k], h_GPS, z_GPS, f_disc_func, F_func, Q_k, H_GPS_func, R_z_GPS, x_update_func)
        xhat_k[:, k + 1] = se2_mat_to_tuplevec(Xhat_k[:, :, k + 1]).reshape((3,))
 
    nu_k = np.array(nu_k)
    K_k = np.array(K_k)
    innovation = np.array(innovation)
    t_k_GPS = np.array(t_k_GPS)
    z_k_GPS = np.array(z_k_GPS)
    k_GPS = np.array(k_GPS)

    return x_true, X_true, xhat_k, Xhat_k, Phat_k, nu_k, K_k, innovation, t_k_GPS, z_k_GPS, k_GPS
x_true, X_true, xhat_k, Xhat_k, Phat_k, nu_k, K_k, innovation, t_k_GPS, z_k_GPS, k_GPS = IEKF_SE2(prob)

# Rotate Phat_k into world frame
def rotate_covariance_to_world(Xhat_k, Phat_k):
    Rhat_k_mat = np.block([[Xhat_k[0:2, 0:2], np.zeros((2, 1))], [np.zeros((1, 2)), 1]])
    Phat_k_world = Rhat_k_mat @ Phat_k @ Rhat_k_mat.T
    return Phat_k_world

Phat_k_world = np.zeros_like(Phat_k)
for k in range(prob.N_steps):
    Phat_k_world[:, :, k] = rotate_covariance_to_world(Xhat_k[:, :, k], Phat_k[:, :, k])

# 2d) Animation

fig, ax = plt.subplots()
ax.set_title("World Frame IEKF")
ax.set_xlabel("X Position [m]")
ax.set_ylabel("Y Position [m]")
ax.grid()
plt.gca().set_aspect('equal', adjustable = 'box')
run_plot, run_quiver, ref_plot, ref_quiver = plot_scatter_unicycles_ax(X_true[:, :, 0], Xhat_k[:, :, 0].reshape((3, 3, 1, 1)), True, ax)
cov_plot = plot_cov_ellipse_ax(Phat_k_world[0:2, 0:2, 0], xhat_k[0:2, 0], "red", "-", ax)
gps_plot = ax.scatter(z_k_GPS[k_GPS <= prob.N_nodes, 0], z_k_GPS[k_GPS <= prob.N_nodes, 1], color = "green")

def update(k):
    ref_plot[0].set_xdata(x_true[0, :k])
    ref_plot[0].set_ydata(x_true[1, :k])
    run_plot[0].set_xdata(xhat_k[0, :k])
    run_plot[0].set_ydata(xhat_k[1, :k])
    e_data = create_cov_ellipse(Phat_k_world[0:2, 0:2, k], xhat_k[0:2, k])
    cov_plot[0].set_xdata(e_data[0, :] + xhat_k[0, k])
    cov_plot[0].set_ydata(e_data[1, :] + xhat_k[1, k])
    gps_plot.set_offsets(z_k_GPS[k_GPS <= k, :])
    return (gps_plot, cov_plot[0], ref_plot, run_plot)

ani = animation.FuncAnimation(fig=fig, func=update, frames = prob.N_nodes, interval=10)

plt.show()
#plt.savefig("./HW4/Parts/Q1/AAE590LGM_HW4_Q2_ref.png")


# 2e) Monte Carlo Evaluation
nu_k_runs = np.zeros((N_runs, t_k_GPS.size))
eps_k_runs = np.zeros((N_runs, prob.N_nodes))
for i in range(N_runs):
    x_true, X_true, xhat_k, Xhat_k, Phat_k, nu_k_runs[i, :], K_k, innovation, t_k_GPS, z_k_GPS, k_GPS = IEKF_SE2(prob)
    
    Phat_k_world = np.zeros_like(Phat_k)
    for k in range(prob.N_steps):
        Phat_k_world[:, :, k] = rotate_covariance_to_world(Xhat_k[:, :, k], Phat_k[:, :, k])

    '''xhat_error = x_true - xhat_k
    xhat_error[2] = np.mod(xhat_error[2] + np.pi, 2 * np.pi) - np.pi
    eps_k = np.zeros((1, prob.N_nodes))
    for k in range(2, prob.N_nodes - 4):
        eps_k_runs[i, k] = NEES(xhat_error[:, k], Phat_k_world[:, :, k])'''

    for k in range(2, prob.N_nodes - 4):
        Xhat_error = se2_log(se2_compose(se2_inverse(Xhat_k[:, :, k]), X_true[:, :, k, 0]))
        eps_k_runs[i, k] = NEES(Xhat_error, Phat_k[:, :, k])

    print(i)
    
nu_k_avg = np.average(nu_k_runs, 0)
eps_k_avg = np.average(eps_k_runs, 0); 

eps_bounds = averaged_chi_squared_consistency_bands(N_runs, 3)
nu_bounds = averaged_chi_squared_consistency_bands(N_runs, 2)

plt.subplot(1, 2, 1)
plt.plot(prob.t_k, eps_k_avg, color = "red")
plt.title("Time Averaged NEES")
plt.hlines(eps_bounds, 0, prob.T, linestyles = "--", label="95% bounds")
plt.grid()
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(t_k_GPS, nu_k_avg, color = "red")
plt.title("Time Averaged NIS")
plt.hlines(nu_bounds, 0, prob.T, linestyles = "--", label="95% bounds")
plt.grid()
plt.legend()

plt.savefig("./HW4/Parts/Q4/AAE590LGM_HW4_Q4d_mc.png")

plt.show()