import numpy as np
import matplotlib.pyplot as plt
# import numba as nb
from Code.SE22.SE22_maps import SE22_compose, SE22_exp, SE22_log, SE22_wedge, SE22_vee, SE22_inverse
from Code.SE22.SE22_integrators import SE22_Lie_group_integrator

# Test Invariant Extended Kalman filter on a flying object in 2D
# IMU for propagation, GPS
# Left and Right IEKFs
# Test against EKF?
# Try numba for speeding up

## Initialize
X_0 = np.eye(4)

a = SE22_compose(X_0, X_0)

# Define vehicle

# Define simulation parameters

# Get reference trajectory

# Define dynamics

## Define Measurements
# Accelerometer


# Gyrometer


# GPS measurements


# Initialize filter

## Simulation loop
# Filter

# Simulate forward

## Plot
# Trajectory plot

# State error history
