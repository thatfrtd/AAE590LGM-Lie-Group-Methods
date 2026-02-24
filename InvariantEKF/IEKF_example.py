import numpy as np
import matplotlib.pyplot as plt
# import numba as nb
from Code.SE22.SE22_maps import SE22_compose, SE22_exp, SE22_log, SE22_wedge, SE22_vee, SE22_inverse
from Code.SE22.SE22_integrators import SE22_Lie_group_integrator

# Test equivariant Kalman filter on a flying object in 2D
# IMU for propagation, GPS
# Test against EKF?
# Try numba for speeding up

## Initialize
# Define vehicle

# Define simulation parameters

# Get reference trajectory

# Define dynamics

# Define measurements

# Initialize filter

## Simulation loop
# Filter

# Simulate forward

## Plot
# Trajectory plot

# State error history

# Bias convergence
