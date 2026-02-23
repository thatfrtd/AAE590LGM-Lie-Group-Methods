import numpy as np
import matplotlib.pyplot as plt
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_log, se2_wedge, se2_vee, se2_inverse
from Code.SE2.SE2_integrators import se2_euler_integrator, se2_Lie_group_integrator

## Unicycle sim

# Initialize
X_0 = np.eye(3)
v = 1 # [m / s]
omega = 0.5 # [rad / s]
xi_k = np.array([[v], [0], [omega]])
dt = 0.1 # [s]
tf = 20
time = np.arange(0, tf + dt, dt)

# Simulation loop
X_k_euler = X_0
X_k_Lie = X_0

positions_euler = np.zeros([2, time.size])
positions_Lie = np.zeros([2, time.size])

error_euler = np.zeros_like(time)
error_Lie = np.zeros_like(time)

for k in range(time.size - 1):
   X_k_euler = se2_euler_integrator(X_k_euler, xi_k, dt)
   X_k_Lie = se2_Lie_group_integrator(X_k_Lie, xi_k, dt)

   positions_euler[:, k + 1] = X_k_euler[0:2, 2]
   positions_Lie[:, k + 1] = X_k_Lie[0:2, 2]

   R_euler = X_k_euler[0:2, 0:2]
   R_Lie = X_k_Lie[0:2, 0:2]

   error_euler[k + 1] = np.linalg.norm(R_euler.T @ R_euler - np.eye(2), "fro")
   error_Lie[k + 1] = np.linalg.norm(R_Lie.T @ R_Lie - np.eye(2), "fro")

# Plot
# Trajectories
plt.plot(positions_euler[0, :], positions_euler[1, :], label = "Euler")
plt.plot(positions_Lie[0, :], positions_Lie[1, :], label = "Lie Group")
plt.xlabel("X Position [m]")
plt.ylabel("Y Position [m]")
plt.title("Unicycle Trajectories vs Integrator")
plt.legend()
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("./HW2/Parts/Q4/AAE590LGM_HW2_Q4_traj.png")
plt.show()

# Error
plt.semilogy(time, error_euler, "o", label = "Euler")
plt.semilogy(time, error_Lie, "o", label =  "Lie Group")
plt.xlabel("Time [s]")
plt.ylabel("Rotation Error")
plt.title("Rotation Error vs Time for Both Integrators")
plt.legend()
plt.grid(True)
plt.savefig("./HW2/Parts/Q4/AAE590LGM_HW2_Q4_error.png")
plt.show()

# Compare with Exp map
X_f = se2_compose(X_0, se2_exp(tf * xi_k))
print(np.linalg.norm(X_f - X_k_Lie, "fro"))
