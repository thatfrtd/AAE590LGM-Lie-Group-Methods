import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
#import numba as nb
from Code.SE3.SE3_maps import se3_compose, se3_exp, se3_log, se3_wedge, se3_vee, se3_inverse
from Code.SE3.SE3_integrators import se3_Lie_group_integrator

# Test intrinsic SCP on Lie groups
# Compare with standard SCP? - would be good for stochastic comparison...
# Try numba for speeding up
# Think about how to make stochastic

## Initialize
# Define vehicle

# Define simulation parameters

# Get track midline reference trajectory

# Define dynamics

## Define convex problem

# Define constraints

# Define objective

# Solve with iSCP

## Plot
# Trajectory plot

# State and control histories