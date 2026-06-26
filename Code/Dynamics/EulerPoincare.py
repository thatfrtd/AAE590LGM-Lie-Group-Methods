import numpy as np

def euler_poincare(xi, J_b, ad_xi, u):
    xi_dot = np.linalg.inv(J_b) @ (ad_xi.T @ J_b @ xi + u)
    return xi_dot