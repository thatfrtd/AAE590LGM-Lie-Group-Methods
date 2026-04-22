import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_log, se2_wedge, se2_vee, se2_inverse, se2_mat_to_tuplevec, se2_tuplevec_to_mat
from Code.SE2.SE2_integrators import se2_euler_integrator, se2_Lie_group_integrator
from Code.Plotting.covariance_plotting import plot_cov_ellipse, calc_and_plot_cov_ellipse, cov_ellipse_3D, create_cov_ellipse, plot_cov_ellipse_ax
from Code.Plotting.SE2_plotting import plot_scatter_unicycles, plot_scatter_unicycles_with_cov, plot_Lie_algebra_cov_ellipse, plot_scatter_unicycles_ax
from HW4.Parts.AAE590LGM_HW4_Problem import AAE590LGM_HW4_Problem
from Code.Filtering.KalmanFilter import EKF
from Code.Filtering.filter_performance_metrics import averaged_chi_squared_consistency_bands, NEES

# Problem 2: World-Frame EKF

## Unicycle SE2 Dynamics
# Xdot = X xi^hat
# xi = [v; 0; omega]

## Initialize Problem Parameters
# Simulation
prob = AAE590LGM_HW4_Problem()

# Monte Carlo
N_runs = 100

def EKF_SE2(prob):
    t_ref = prob.t_k
    
    def f_disc(x_k, xi_k, dt):
        x_kp1 = [x_k[0] + xi_k[0] * sp.cos(x_k[2]) * dt, 
                 x_k[1] + xi_k[0] * sp.sin(x_k[2]) * dt, 
                 x_k[2] + xi_k[2] * dt]
        return x_kp1

    # 2a) Prediction Step
    x_sym = sp.symbols("x:3")
    xi_sym = sp.symbols("xi:3")
    f_disc_sym = sp.Matrix(f_disc(x_sym, xi_sym, prob.delta_t))
    F = f_disc_sym.jacobian(x_sym)
    F_func = sp.lambdify((x_sym, xi_sym), F, "numpy")
    
    f_disc_func = sp.lambdify((x_sym, xi_sym), f_disc_sym, "numpy")
    
    L_k = lambda theta : prob.delta_t * np.array([[np.cos(theta), 0, 0], [np.sin(theta), 0, 0], [0, 0, 1]])
    Q_k_world = lambda x : L_k(x[2]) @ prob.Q_body @ L_k(x[2]).T
    
    # 2b) Measurement Update
    h_GPS = lambda x : x[0:2]
    h_GPS_noise = lambda x : prob.rng.normal(h_GPS(x), prob.sigma_GPS)    
    R_z_GPS = prob.sigma_GPS ** 2 * np.eye(2)
    H_GPS = np.hstack((np.eye(2), np.zeros((2, 1))))
    H_GPS_func = lambda x : H_GPS
    
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
            z_GPS = h_GPS_noise(x_true[:, k])
            z_k_GPS.append(z_GPS)
            t_k_GPS.append(t_ref[k])
            k_GPS.append(k)
    
            xhat_k[:, k + 1], Phat_k[:, :, k + 1], nu_k_i, K_k_i, innovation_i = EKF(xhat_k[:, k], Phat_k[:, :, k], prob.xi_k[:, k], h_GPS, z_GPS, f_disc_func, F_func, Q_k_world, H_GPS_func, R_z_GPS)
    
            nu_k.append(nu_k_i)
            K_k.append(K_k_i)
            innovation.append(innovation_i)
        else:
            z_GPS = np.nan    
    
            xhat_k[:, k + 1], Phat_k[:, :, k + 1] = EKF(xhat_k[:, k], Phat_k[:, :, k], prob.xi_k[:, k], h_GPS, z_GPS, f_disc_func, F_func, Q_k_world, H_GPS_func, R_z_GPS)
        Xhat_k[:, :, k + 1] = se2_tuplevec_to_mat(xhat_k[:, k + 1]).reshape((3, 3))
 
    nu_k = np.array(nu_k)
    K_k = np.array(K_k)
    innovation = np.array(innovation)
    t_k_GPS = np.array(t_k_GPS)
    z_k_GPS = np.array(z_k_GPS)
    k_GPS = np.array(k_GPS)

    return x_true, X_true, xhat_k, Xhat_k, Phat_k, nu_k, K_k, innovation, t_k_GPS, z_k_GPS, k_GPS
x_true, X_true, xhat_k, Xhat_k, Phat_k, nu_k, K_k, innovation, t_k_GPS, z_k_GPS, k_GPS = EKF_SE2(prob)

# 2d) Animation

fig, ax = plt.subplots()
ax.set_title("World Frame EKF")
ax.set_xlabel("X Position [m]")
ax.set_ylabel("Y Position [m]")
ax.grid()
plt.gca().set_aspect('equal', adjustable = 'box')
run_plot, run_quiver, ref_plot, ref_quiver = plot_scatter_unicycles_ax(X_true[:, :, 0], Xhat_k[:, :, 0].reshape((3, 3, 1, 1)), True, ax)
cov_plot = plot_cov_ellipse_ax(Phat_k[0:2, 0:2, 0], xhat_k[0:2, 0], "red", "-", ax)
gps_plot = ax.scatter(z_k_GPS[k_GPS <= prob.N_nodes, 0], z_k_GPS[k_GPS <= prob.N_nodes, 1], color = "green")

def update(k):
    ref_plot[0].set_xdata(x_true[0, :k])
    ref_plot[0].set_ydata(x_true[1, :k])
    run_plot[0].set_xdata(xhat_k[0, :k])
    run_plot[0].set_ydata(xhat_k[1, :k])
    e_data = create_cov_ellipse(Phat_k[0:2, 0:2, k], xhat_k[0:2, k])
    cov_plot[0].set_xdata(e_data[0, :] + xhat_k[0, k])
    cov_plot[0].set_ydata(e_data[1, :] + xhat_k[1, k])
    gps_plot.set_offsets(z_k_GPS[k_GPS <= k, :])
    return (gps_plot, cov_plot[0], ref_plot, run_plot)

ani = animation.FuncAnimation(fig=fig, func=update, frames = prob.N_nodes, interval=10)

#plt.savefig("./HW4/Parts/Q1/AAE590LGM_HW4_Q2_ref.png")
plt.show()

# 2e) Monte Carlo Evaluation
nu_k_runs = np.zeros((N_runs, t_k_GPS.size))
eps_k_runs = np.zeros((N_runs, prob.N_nodes))
for i in range(N_runs):
    x_true, X_true, xhat_k, Xhat_k, Phat_k, nu_k_runs[i, :], K_k, innovation, t_k_GPS, z_k_GPS, k_GPS = EKF_SE2(prob)
    
    xhat_error = x_true[:, 0:] - xhat_k[:, 0:]
    xhat_error[2] = np.mod(xhat_error[2] + np.pi, 2 * np.pi) - np.pi

    eps_k = np.zeros((1, prob.N_nodes))
    for k in range(2, prob.N_nodes - 1):
        eps_k_runs[i, k] = NEES(xhat_error[:, k], Phat_k[:, :, k])

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

plt.savefig("./HW4/Parts/Q2/AAE590LGM_HW4_Q2e_mc.png")

plt.show()