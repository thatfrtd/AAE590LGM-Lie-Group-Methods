import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
#import numba as nb
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_log, se2_wedge, se2_vee, se2_inverse, se2_ad
from Code.SE2.SE2_integrators import se2_Lie_group_integrator

# Simulate PID controlled aircraft on SE(2) Lie group
# No estimator, just controller (maybe noisy state measurements)
# Start with 2D then move to 3D (SE(3))
# Tuning through convex optimization?
# Try numba for speeding up

## Initialize
# Define vehicle (including input constraints)

# Define simulation parameters


# Get track midline reference trajectory

# Define PID controller
K_p = np.eye(3)
K_i = np.eye(3)
K_d = np.eye(3)

# Define dynamics


## Simulation loop
# Compute PID control input
psi_k = left_invariant_error(X_k, X_d_k)
psidot_k = left_invariant_error_derivative(psi, xi_k, xi_d_k)
u_k = PID(psi_k, psidot_k, K_p, K_i, K_d)

# Simulate forward
X_kp1 = se2_Lie_group_integrator(X_k, u_k, dt)

## Plot
# Trajectory plot

# State and control histories

## Helper Functions
def left_invariant_error(X, X_d):

    Psi_l = se2_compose(se2_inverse(X_d), X)

    return Psi_l

def right_invariant_error(X, X_d):

    Psi_r = se2_compose(X, se2_inverse(X_d))

    return Psi_r

def left_invariant_error_derivative(psi, xi, xi_d):
    '''analytically calculate left invariant error derivative'''
    psidot = -se2_ad(xi_d) * psi + (xi - xi_d)

    return psidot

def PID(psi, psidot, K_p, K_i, K_d):

    u = K_p @ psi + K_d @ psidot

    return u