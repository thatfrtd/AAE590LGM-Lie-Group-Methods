import numpy as np
import matplotlib.pyplot as plt

def create_cov_ellipse(P, xbar):
    P_eigval, P_eigvec = np.linalg.eig(P)
    thetas = np.reshape(np.linspace(0, 2 * np.pi, 100), (1, 100))
    ellipse_3sigma_data = P_eigvec @ np.reshape(np.vstack((3 * np.sqrt(P_eigval[0]) * np.cos(thetas), 3 * np.sqrt(P_eigval[1]) * np.sin(thetas))), (2, 100))
    return ellipse_3sigma_data

def plot_cov_ellipse_ax(P, xbar, color, linestyle, ax):
    ellipse_3sigma_data = create_cov_ellipse(P, xbar)
    eplot = ax.plot(xbar[0] + ellipse_3sigma_data[0, :], xbar[1] + ellipse_3sigma_data[1, :], color = color, linestyle = linestyle)
    return eplot

def plot_cov_ellipse(P, xbar, color, linestyle):
    ellipse_3sigma_data = create_cov_ellipse(P, xbar)
    eplot = plt.plot(xbar[0] + ellipse_3sigma_data[0, :], xbar[1] + ellipse_3sigma_data[1, :], color = color, linestyle = linestyle)
    return eplot

def calc_and_plot_cov_ellipse(X):
    X_cov = np.cov(X)
    X_std = np.sqrt(np.diag(X_cov))
    X_mean = np.mean(X, 1)
    plot_cov_ellipse(X_cov, X_mean, "k", "-")
    return X_mean, X_std, X_cov

def cov_ellipse_3D(P):
    P_eigval, P_eigvec = np.linalg.eig(P)
    phi = np.reshape(np.linspace(0, np.pi, 50), (50,))
    psi = np.reshape(np.linspace(0, 2 * np.pi, 100), (100,))

    ellipse_3sigma_data = np.zeros((3, phi.size, psi.size))
    for i in range(phi.size):
        for j in range(psi.size):
            ellipse_3sigma_data[:, i, j] = np.reshape(P_eigvec @ np.reshape(np.vstack((3 * np.sqrt(P_eigval[0]) * np.sin(phi[i]) * np.cos(psi[j]), 
                                                                                       3 * np.sqrt(P_eigval[1]) * np.sin(phi[i]) * np.sin(psi[j]), 
                                                                                       3 * np.sqrt(P_eigval[2]) * np.cos(phi[i]))), (3, 1)), (3,))
    return ellipse_3sigma_data.reshape((3, -1))