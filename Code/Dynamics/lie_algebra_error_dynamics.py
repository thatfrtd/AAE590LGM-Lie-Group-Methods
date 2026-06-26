import numpy as np
from Code.Dynamics.EulerPoincare import euler_poincare

def lie_algebra_2nd_error_dynamics(eps, J_b, ad, ad_dot, u, xi_bar):
    # Extract error states
    psi = eps[0:6] # Group error
    zeta = eps[6:12] # Twist error (group error rate)

    # Group error dynamics
    xi = xi_bar + zeta

    # Twist error dynamics
    xi_dot = euler_poincare(xi, J_b, ad(xi), u)
    xi_bar_dot = euler_poincare(xi_bar, J_b, ad(xi_bar), u)
    zeta_dot = -ad_dot(xi_bar) @ psi - ad(xi_bar) @ psi_dot + (xi_dot - xi_bar_dot)

    # Package
    eps_dot = [psi_dot, zeta_dot]

    return eps_dot
