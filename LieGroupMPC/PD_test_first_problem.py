import numpy as np
import matplotlib.pyplot as plt
from parfor import pmap
from os import environ
environ['OMP_NUM_THREADS'] = '4'
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_log, se2_wedge, se2_vee, se2_inverse, se2_ljac, se2_inv_ljac, se2_Ad, se2_ad
from Code.SE2.SE2_integrators import se2_euler_integrator, se2_Lie_group_integrator

class PD_test_problem:
    # Simulation
    delta_t = 0.1 # [s] Time step
    T = 15 # [s] Simulation duration
    
    # Reference trajectory
    v_nom = 2 # [m / s] Commanded forward speed
    r_fig8 = 5 # [m] Radius of each lobe of the figure-eight
    omega_turn = 0.4 # [rad / s] Magnitude of the commanded heading rate on each lobe
    T_fig8 = 10 * np.pi # [s] Period of one figure-eight cycle
    T_lobe = T_fig8 / 2 # [s] Period of one lobe
    X_0 = se2_exp(np.array([0, 0, 0]).T)

    # Initial dispersion
    Phat_0 = np.diag((0.01, 0.01, 0.03))

    # Initial filter covariance
    Ptilde_0 = np.diag((0.01, 0.01, 0.01))

    # Process noise
    sigma_v = 0.001 # [m / s] Forward speed standard deviation
    sigma_lat = 0 # [m / s] Lateral speed standard deviation (0 so no sideslip)
    sigma_omega = 0.003 # [rad / s] Heading rate standard deviation
    
    # GPS measurement model parameters
    sigma_GPS = 0.3 # [m] GPS measurement noise standard deviation
    GPS_rate = 1 # [Hz] GPS measurements arrive once per second

    rng = np.random.default_rng(2)
    N_steps = np.round(T / delta_t).astype(int)
    N_nodes = N_steps + 1
    t_k = np.linspace(0, T, N_nodes)

    v_x = v_nom * np.ones((1, N_nodes))
    v_y = np.zeros((1, N_nodes)) # No sideslip
    tau = np.mod(t_k, T_fig8)
    omega = np.where(tau < T_lobe, omega_turn, -omega_turn)
    xi_k = np.vstack((v_x, v_y, omega))
    Q_body = np.diag(np.square([sigma_v, 0, sigma_omega]))
    B_cont = np.array([[1, 0], 
                       [0, 0], 
                       [0, 1.0]])
    A = np.zeros((3, 3, N_nodes))
    A_d = np.zeros((3, 3, N_nodes))
    B_d = np.zeros((3, 2, N_nodes))
    B_d_inv = np.zeros((2, 3, N_nodes))
    
    def __init__(self):
        # Discretize (ONLY WORKS FOR CONSTANT INPUTS (so really only first order systems))
        A = lambda xi_k : -se2_ad(xi_k)
        A_d = lambda xi_k : se2_Ad(se2_exp(-xi_k * self.delta_t))
        B_d = lambda xi_k : np.linalg.pinv(A(xi_k)) @ (A_d(xi_k) - np.eye(3)) @ self.B_cont
        B_d_inv = lambda xi_k : np.linalg.pinv(B_d(xi_k))
        for k in range(self.N_nodes):
            self.A[:, :, k] = A(self.xi_k[:, k])
            self.A_d[:, :, k] = A_d(self.xi_k[:, k])
            self.B_d[:, :, k] = B_d(self.xi_k[:, k])
            self.B_d_inv[:, :, k] = B_d_inv(self.xi_k[:, k])

    def simulate_no_fb(self, xi_noise, X_0): 
        X_k = np.zeros([3, 3, self.N_nodes])
        X_k[:, :, 0] = X_0

        for k in range(self.N_steps):
            X_k[:, :, k + 1] = se2_Lie_group_integrator(X_k[:, :, k], self.xi_k[:, k] + xi_noise[:, k], self.delta_t)

        return X_k

    def simulate(self, xi_noise, u_func, X_bar, X_0): 
        X_k = np.zeros([3, 3, self.N_nodes])
        X_k[:, :, 0] = X_0
        u_fb_k = np.zeros([self.B_d.shape[1], self.N_steps])
        eta_k = np.zeros([3, 3, self.N_nodes])
        xi_e_k = np.zeros([3, self.N_nodes])

        for k in range(self.N_steps):
            eta_k[:, :, k] = se2_inverse(X_bar[:, :, k]) @ X_k[:, :, k]
            xi_e_k[:, k] = se2_log(eta_k[:, :, k]).reshape((3,))
            u_ff = self.xi_k[:, k]
            u_fb_k[:, k] = u_func(xi_e_k[:, k], eta_k[:, :, k], np.zeros([3, 1]), u_ff, k).reshape((-1,))
            X_k[:, :, k + 1] = se2_Lie_group_integrator(X_k[:, :, k], u_ff + self.B_d[:, :, k] @ u_fb_k[:, k] + xi_noise[:, k], self.delta_t)

        return X_k

    def monte_carlo(self, N_runs, u_func, X_bar):
        X_k_runs = np.zeros((3, 3, self.N_nodes, N_runs))

        xi_noise = np.vstack((self.rng.normal(0, self.sigma_v, size = (1, self.N_steps, N_runs)), 
                              np.zeros((1, self.N_steps, N_runs)), 
                              self.rng.normal(0, self.sigma_omega, size = (1, self.N_steps, N_runs))))
        
        a = self.rng.multivariate_normal(np.zeros((3,)), self.Phat_0)
        X_0 = np.zeros((3, 3, N_runs))
        for i in range(N_runs):
            X_0[:, :, i] = se2_compose(self.X_0, se2_exp(self.rng.multivariate_normal(np.zeros((3,)), self.Phat_0)))

        for i in range(N_runs):
            X_k_runs[:, :, :, i] = self.simulate(np.squeeze(xi_noise[:, :, i]), u_func, X_bar, X_0[:, :, i])

        return X_k_runs
    
    def monte_carlo_multithreaded(self, N_runs):
        X_k_runs = np.zeros((3, 3, self.N_nodes, N_runs))

        xi_noise = np.vstack((self.rng.normal(0, self.sigma_v, size = (1, self.N_steps, N_runs)), 
                              np.zeros((1, self.N_steps, N_runs)), 
                              self.rng.normal(0, self.sigma_omega, size = (1, self.N_steps, N_runs))))
        
        # issue tasks and gather results
        sim_func = lambda i: self.simulate(np.squeeze(xi_noise[:, :, i]))
        results = pmap(sim_func, range(N_runs))

        for i in range(N_runs):
            X_k_runs[:, :, :, i] = results[i]

        return X_k_runs
    
    @staticmethod
    def error_rate(xi_e, eta, xi, xi_bar):
        xi_e_dot = se2_inv_ljac(xi_e) @ xi_e @ (xi - se2_Ad(se2_inverse(eta)) @ xi_bar)
        return xi_e_dot
    
    @staticmethod
    def PD(xi_e, eta, xi, xi_bar, K_p, K_d):
        #xi_e_dot = error_rate(xi_e, eta, xi, xi_bar)
        xi_e_dot = se2_inv_ljac(xi_e) @ xi_e @ (xi - se2_Ad(se2_inverse(eta)) @ xi_bar)
        u_fb = -K_p @ xi_e - K_d @ xi_e_dot
        return u_fb
    
    @staticmethod
    def PD_approx(xi_e, eta, xi, xi_bar, K_p, K_d):
        xi_e_dot = xi - se2_Ad(se2_inverse(eta)) @ xi_bar
        u_fb = -K_p @ xi_e - K_d @ xi_e_dot
        return u_fb 
    
    def log_linear_dynamic_inversion(self, xi_e, xi_bar, K_p, k):
        u_fb = self.B_d_inv[:, :, k] @ (se2_ljac(xi_e) @ (-self.B_d[:, :, k] @ K_p @ xi_e) - se2_ad(xi_e) @ xi_bar)
        return u_fb