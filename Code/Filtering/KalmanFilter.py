import numpy as np
from Code.Filtering.filter_performance_metrics import NIS

def discrete_covariance_prediction(P_k_k, F_k, Q_k):
    P_kp1_k = F_k @ P_k_k @ F_k.T + Q_k
    return P_kp1_k

def Kalman_gain(P_k_km1, H, R_z):
    S_k = H @ P_k_km1 @ H.T + R_z
    S_k_inv = np.linalg.inv(S_k)
    K_k = P_k_km1 @ (H.T @ S_k_inv)
    return K_k, S_k_inv

def covariance_correction(P_k_km1, K_k, H, R_z):
    P_k_k = (np.eye(3) - K_k @ H) @ P_k_km1 @ (np.eye(3) - K_k @ H).T + K_k @ R_z @ K_k.T
    return P_k_k

def EKF(x_km1_km1, P_km1_km1, xi_k, h, z_k, f_disc, F, Q_k, H, R_z):
    n_x = x_km1_km1.size

    # Predict
    x_k_km1 = f_disc(x_km1_km1, xi_k)
    P_k_km1 = discrete_covariance_prediction(P_km1_km1, F(x_km1_km1, xi_k), Q_k(x_km1_km1))

    H_k = H(x_k_km1)

    # Correct if Measurement Exists
    if not np.all(np.isnan(z_k)):
        # Correct
        K_k, S_k_inv = Kalman_gain(P_k_km1, H_k, R_z)
        zhat_k = h(x_k_km1).reshape(z_k.shape) # Expected measurement
        innovation = z_k - zhat_k
    
        x_k_k = x_k_km1.reshape(n_x,) + K_k @ innovation
        P_k_k = covariance_correction(P_k_km1, K_k, H_k, R_z)

        # Filter performance metric (Normalized Innovation Squared - NIS)
        nu_k = NIS(innovation, S_k_inv)
        return x_k_k.reshape((n_x,)), P_k_k.reshape((n_x, n_x)), nu_k, K_k, innovation
    else:
        return x_k_km1.reshape((n_x,)), P_k_km1.reshape((n_x, n_x))

# def ESKF_update():
def ESKF(x_km1_km1, P_km1_km1, xi_k, h, z_k, f_disc, F, Q_k, H, R_z, x_update_func):
    n_x = np.size(x_km1_km1, axis = 0)

    # Predict
    x_k_km1 = f_disc(x_km1_km1, xi_k)
    P_k_km1 = discrete_covariance_prediction(P_km1_km1, F(x_km1_km1, xi_k), Q_k(x_km1_km1))
    
    H_k = H(x_k_km1)

    # Correct if Measurement Exists
    if not np.all(np.isnan(z_k)):
        # Correct
        K_k, S_k_inv = Kalman_gain(P_k_km1, H_k, R_z)
        zhat_k = h(x_k_km1).reshape(z_k.shape) # Expected measurement
        innovation = z_k - zhat_k

        delta_x = K_k @ innovation
        x_k_k = x_update_func(x_k_km1.reshape(x_km1_km1.shape), delta_x)
        P_k_k = covariance_correction(P_k_km1, K_k, H_k, R_z)

        # Filter performance metric (Normalized Innovation Squared - NIS)
        nu_k = NIS(innovation, S_k_inv)
        return x_k_k.reshape(x_km1_km1.shape), P_k_k.reshape((n_x, n_x)), np.squeeze(nu_k), K_k, innovation
    else:
        return np.squeeze(x_k_km1.reshape(x_km1_km1.shape)), P_k_km1.reshape((n_x, n_x))

# def LIEKF_update():