import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
#import numba as nb
from Code.SE22.SE22_maps import se22_compose, se22_exp, se22_log, se22_wedge, se22_vee, se22_inverse
from Code.SE22.SE22_integrators import se22_Lie_group_integrator

# Simulate MPC controlled aircraft on SE2(2) Lie group
# No estimator, just controller (maybe noisy state measurements)
# Should compare to PID controller
# Start with 2D then move to 3D (SE2(3))
# Tuning with differentiating through convex program?
# Try numba for speeding up discretization

## Initialize
# Define vehicle

# Define simulation parameters

# Get track midline reference trajectory

# Define baseline PID controller?

# Define dynamics

## Define convex problem - should be QP

# Define constraints

# Define objective

## Simulation loop
# Discretize dynamics

# Solve convex problem

# Simulate forward

## Plot
# Trajectory plot with PID comparison

# State and control histories