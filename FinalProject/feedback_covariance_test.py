import numpy as np

# Feedback Covariance Test
# Test how feedback in the Lie algebra affects the state covariance
# Need to use Dynamic Inversion for the feedback?

def covariance_update(P_k, A_k, B_k, K_k):
    P_kp1 = (A_k - B_k @ K_k) @ P_k @ (A_k - B_k @ K_k).T
    return P_kp1

def linear_feedback(ubar_k, error_k, K_k):
    u_k = ubar_k + K_k @ error_k
    return u_k

def dyn_inv_feedback(ubar_k, error_k, K_k, B_k, U_l_B_k): # good or not? how does U_l_B_k get discretized?
    u_k = ubar_k + np.linalg.inv(U_l_B_k) @ B_k @ K_k @ error_k
    return u_k