import numpy as np
import matplotlib.pyplot as plt
# import numba as nb
from Code.SE23.SE23_maps import se23_compose, se23_exp, se23_log, se23_wedge, se23_vee, se23_inverse
from Code.SE23.SE23_integrators import se23_Lie_group_integrator

# Test equivariant Kalman filter on a flying object with accelerometer and gyroscopic biases
# IMU for propagation, GPS, barometer, magnetometer
# Test against MEKF? - Pablo's code...
# Test against (imperfect) InEKF?
# Test against TFG-InEKF?
# Try numba for speeding up
# Add angular velocity as a state to track

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
