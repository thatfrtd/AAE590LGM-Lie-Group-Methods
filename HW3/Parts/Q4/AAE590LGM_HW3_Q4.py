import numpy as np
import matplotlib.pyplot as plt
from Code.SE22.SE22_maps import se22_compose, se22_exp, se22_log, se22_wedge, se22_vee, se22_inverse
from Code.SE22.SE22_integrators import se22_Lie_group_integrator

## 1b.2) Propagate waypoints
X_0 = np.eye(4)

# Segment info
N_seg = 4
a = np.array([[2, 0, 0,-1],
              [0, 0, 0, 0]]) # [m / s]
omega = np.array([0, 0, 0.3, 0]) # [rad / s]
T = np.array([5, 5, np.pi / 0.3, 5]) # [s]

# xi_k = np.array([[v_x], [v_y], [omega]])
# dt = 0.04 # [s]

# # Simulation loop
# X_k = np.zeros([3, 3, N_seg + 1]) # Waypoint group elements
# X_k[:, :, 0] = X_0

# T_seg = np.cumsum(T)
# time = np.arange(0, T_seg[-1], dt)
# k_seg = np.round(T / dt).astype(int)
# k_seg_total = sum(k_seg) + 1
# X_seg = np.zeros([3, 3, k_seg_total])
# X_seg[:, :, 0] = X_0
# xi_seg = np.hstack((np.matlib.repmat(xi_k[:, :, 0], 1, k_seg[0]),
#                     np.matlib.repmat(xi_k[:, :, 1], 1, k_seg[1]),
#                     np.matlib.repmat(xi_k[:, :, 2], 1, k_seg[2]),
#                     np.matlib.repmat(xi_k[:, :, 3], 1, k_seg[3]),
#                     np.matlib.repmat(xi_k[:, :, 4], 1, k_seg[4]),
#                     np.matlib.repmat(xi_k[:, :, 5], 1, k_seg[5]),
#                     np.matlib.repmat(xi_k[:, :, 6], 1, k_seg[6]),
#                     np.matlib.repmat(xi_k[:, :, 7], 1, k_seg[7])))

# for k in range(N_seg):
#    X_k[:, :, k + 1] = se22_Lie_group_integrator(X_k[:, :, k], xi_k[:, :, k], T[k])
   
# for s in range(k_seg_total - 1):
#     X_seg[:, :, s + 1] = se22_Lie_group_integrator(X_seg[:, :, s], xi_seg[:, s], dt)

# # Plot
# # Trajectories
# plt.plot(X_seg[0, 2, :], X_seg[1, 2, :])
# plt.quiver(X_k[0, 2, :], X_k[1, 2, :], X_k[0, 0, :], X_k[1, 0, :])
# for k in range(N_seg):
#     plt.text(X_k[0, 2, k], X_k[1, 2, k], f"Segment {k}", fontsize=10, 
#               bbox=dict(facecolor='red', alpha=0.5))
# plt.xlabel("X Position [m]")
# plt.ylabel("Y Position [m]")
# plt.title("Reference Trajectory in (x, y) Plane")
# plt.legend()
# plt.grid()
# plt.gca().set_aspect('equal', adjustable='box')
# #plt.savefig("./HW2/Parts/Q4/AAE590LGM_HW2_Q4_traj.png")
# plt.show()

# Simulate spiral
dt_spiral = 0.01 # [s]
tf = 10 # [s]
N_t = np.round(tf / dt_spiral).astype(int)
X_k_spiral = np.zeros([4, 4, N_t + 1])
omega_spiral = 0.5 # [rad / s]
a_spiral = np.array([[1], [0]])
xi_spiral = np.block([[a_spiral], [np.zeros((2, 1))], [omega_spiral]])
X_k_spiral[:, :, 0] = X_0

for k in range(N_t):
    X_k_spiral[:, :, k + 1] = se22_Lie_group_integrator(X_k_spiral[:, :, k], xi_spiral, dt_spiral)

plt.plot(X_k_spiral[0, 3, :], X_k_spiral[1, 3, :], label = "Line")
v_ind = np.round(np.linspace(0, N_t, 10)).astype(int)
plt.quiver(X_k_spiral[0, 3, v_ind], X_k_spiral[1, 3, v_ind], X_k_spiral[0, 2, v_ind], X_k_spiral[1, 2, v_ind])
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.title("Spiral")
plt.legend()
plt.gca().set_aspect('equal', adjustable='box')
plt.grid()
plt.show()

a = X_k_spiral[0:2, 2, :]
vel_mag = np.linalg.vector_norm(X_k_spiral[[0, 1], 2, :], axis = 0)
plt.plot(np.linspace(0, tf, N_t + 1), vel_mag)
plt.grid()
plt.show()