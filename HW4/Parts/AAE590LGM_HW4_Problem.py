import numpy as np
import matplotlib.pyplot as plt
from parfor import pmap
from os import environ
environ['OMP_NUM_THREADS'] = '4'
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_log, se2_wedge, se2_vee, se2_inverse
from Code.SE2.SE2_integrators import se2_euler_integrator, se2_Lie_group_integrator

class AAE590LGM_HW4_Problem:
    # Simulation
    delta_t = 0.1 # [s] Time step
    T = 35 # [s] Simulation duration
    
    # Reference trajectory
    v_nom = 2 # [m / s] Commanded forward speed
    r_fig8 = 5 # [m] Radius of each lobe of the figure-eight
    omega_turn = 0.4 # [rad / s] Magnitude of the commanded heading rate on each lobe
    T_fig8 = 10 * np.pi # [s] Period of one figure-eight cycle
    T_lobe = T_fig8 / 2 # [s] Period of one lobe

    # Process noise
    sigma_v = 1 # [m / s] Forward speed standard deviation
    sigma_lat = 0 # [m / s] Lateral speed standard deviation (0 so no sideslip)
    sigma_omega = 0.05 # [rad / s] Heading rate standard deviation
    
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

    ## Create nominal reference trajectory (figure 8)
    def simulate(self, xi_noise):
        X_0 = np.eye(3)
 
        X_k = np.zeros([3, 3, self.N_nodes])
        X_k[:, :, 0] = X_0

        for k in range(self.N_steps):
            X_k[:, :, k + 1] = se2_Lie_group_integrator(X_k[:, :, k], self.xi_k[:, k] + xi_noise[:, k], self.delta_t)

        return X_k
    
    def monte_carlo(self, N_runs):
        X_k_runs = np.zeros((3, 3, self.N_nodes, N_runs))

        xi_noise = np.vstack((self.rng.normal(0, self.sigma_v, size = (1, self.N_steps, N_runs)), 
                              np.zeros((1, self.N_steps, N_runs)), 
                              self.rng.normal(0, self.sigma_omega, size = (1, self.N_steps, N_runs))))
        
        for i in range(N_runs):
            X_k_runs[:, :, :, i] = self.simulate(np.squeeze(xi_noise[:, :, i]))

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