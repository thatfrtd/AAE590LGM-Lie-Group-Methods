import numpy as np
import matplotlib.pyplot as plt
from Code.Plotting.covariance_plotting import cov_ellipse_3D, calc_and_plot_cov_ellipse
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_log, se2_wedge, se2_vee, se2_inverse

def plot_scatter_unicycles(X_ref, X_runs, plot_lines):
    N_runs = np.size(X_runs, 3)

    arrow_scale = 0.1
    plt.xlabel("X Position [m]")
    plt.ylabel("Y Position [m]")
    plt.grid()
    plt.gca().set_aspect('equal', adjustable = 'box')
    if plot_lines:
        for i in range(N_runs):
            plt.plot(X_runs[0, 2, :, i], X_runs[1, 2, :, i], linestyle = "--")
    else:
        for i in range(N_runs):
            plt.scatter(X_runs[0, 2, -1, i], X_runs[1, 2, -1, i], alpha = 0.5)
    for i in range(N_runs):
        plt.quiver(X_runs[0, 2, -1, i], X_runs[1, 2, -1, i], X_runs[0, 0, -1, i], X_runs[1, 0, -1, i], scale = 1/arrow_scale, alpha = 0.25)
    
    plt.plot(X_ref[0, 2, :], X_ref[1, 2, :],'k',linewidth = 0.7)
    plt.quiver(X_ref[0, 2, -1], X_ref[1, 2, -1], X_ref[0, 0, -1], X_ref[1, 0, -1], scale = 1/arrow_scale)

def plot_scatter_unicycles_ax(X_ref, X_runs, plot_lines, ax):
    N_runs = np.size(X_runs, 3)

    arrow_scale = 0.1
    if plot_lines:
        for i in range(N_runs):
            run_plot = ax.plot(X_runs[0, 2, :, i], X_runs[1, 2, :, i], linestyle = "--")
    else:
        for i in range(N_runs):
            ax.scatter(X_runs[0, 2, -1, i], X_runs[1, 2, -1, i], alpha = 0.5)
    for i in range(N_runs):
        run_quiver = ax.quiver(X_runs[0, 2, -1, i], X_runs[1, 2, -1, i], X_runs[0, 0, -1, i], X_runs[1, 0, -1, i], scale = 1/arrow_scale, alpha = 0.25)
    
    ref_plot = ax.plot(X_ref[0, 2, :], X_ref[1, 2, :],'k',linewidth = 0.7)
    ref_quiver = ax.quiver(X_ref[0, 2, -1], X_ref[1, 2, -1], X_ref[0, 0, -1], X_ref[1, 0, -1], scale = 1/arrow_scale)
    return (run_plot, run_quiver, ref_plot, ref_quiver)

def plot_scatter_unicycles_with_cov(X_ref, X_runs, plot_lines):
    N_runs = np.size(X_runs, 3)
    
    arrow_scale = 0.1
    plt.xlabel("X Position [m]")
    plt.ylabel("Y Position [m]")
    
    if plot_lines:
        for i in range(N_runs):
            plt.plot(X_runs[0, 2, :, i], X_runs[1, 2, :, i])
    else:
        for i in range(N_runs):
            plt.scatter(X_runs[0, 2, -1, i], X_runs[1, 2, -1, i], alpha = 0.5)
    for i in range(N_runs):
        plt.quiver(X_runs[0, 2, -1, i], X_runs[1, 2, -1, i], X_runs[0, 0, -1, i], X_runs[1, 0, -1, i], scale = 1/arrow_scale, alpha = 0.25)
    
    plt.plot(X_ref[0, 2, :], X_ref[1, 2, :],'k')
    plt.quiver(X_ref[0, 2, -1], X_ref[1, 2, -1], X_ref[0, 0, -1], X_ref[1, 0, -1], scale = 1/arrow_scale)

    r_f = np.squeeze(X_runs[0:2, 2, -1, :])
    calc_and_plot_cov_ellipse(r_f)

def plot_Lie_algebra_cov_ellipse(P, Xbar, color, linestyle): # Should change to be more general and go in 'covariance_plotting.py'
    if Xbar.ndim > 2:
        N_runs = np.size(Xbar, 3)
    else:
        N_runs = 1

    ellipse_3sgima_data = cov_ellipse_3D(P)
    X_boundary = np.zeros((3, 3, np.size(ellipse_3sgima_data, 1)))
    for i in range(np.size(ellipse_3sgima_data, 1)):
        X_boundary[:, :, i] = se2_compose(Xbar, se2_exp(ellipse_3sgima_data[:, i]))
    r_boundary = X_boundary[0:2, 2, :]
    plt.scatter(r_boundary[0, :], r_boundary[1, :], 2, color = color, linestyle = linestyle, alpha = 0.5) 
    