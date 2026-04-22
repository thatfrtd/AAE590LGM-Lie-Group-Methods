import numpy as np
import matplotlib.pyplot as plt
from Code.SE22.SE22_maps import se22_compose, se22_exp, se22_log, se22_wedge, se22_vee, se22_inverse
from Code.SE22.SE22_integrators import se22_Lie_group_integrator

def turning_radius(velocity, desired_turn_rate):
    g = 9.81 # [m / s2]

    turn_radius = velocity ** 2 / (g * np.tan(desired_turn_rate))

    return turn_radius

## 1b.2) Propagate waypoints
X_0 = np.eye(4)

# Segment info
dt = 0.01 # [s]
N_seg = 4
a = np.array([[2, 0, 0,-1],
              [0, 0, 0, 0]]) # [m / s]
omega = np.array([0, 0, 0.3, 0]).reshape([1, 4]) # [rad / s]
T = np.array([5, 5, np.round(np.pi / 0.3 / dt) * dt, 5]) # [s]

xi_k = np.vstack((a, np.zeros_like(a), omega)).reshape([5, 1, -1])

# Simulation loop
X_k = np.zeros([4, 4, N_seg + 1]) # Waypoint group elements
X_k[:, :, 0] = X_0

T_seg = np.cumsum(T)
time = np.arange(0, T_seg[-1], dt)
k_seg = np.round(T / dt).astype(int)
k_seg_total = sum(k_seg) + 1
X_seg = np.zeros([4, 4, k_seg_total])
X_seg[:, :, 0] = X_0
xi_seg = np.hstack((np.matlib.repmat(xi_k[:, :, 0], 1, k_seg[0]),
                    np.matlib.repmat(xi_k[:, :, 1], 1, k_seg[1]),
                    np.matlib.repmat(xi_k[:, :, 2], 1, k_seg[2]),
                    np.matlib.repmat(xi_k[:, :, 3], 1, k_seg[3])))

for k in range(N_seg):
   xi = xi_k[:, :, k]
   if xi[4] != 0: # Add centripital acceleration
       rvec = X_k[0:2, 3, k].reshape([2, 1])
       r = np.linalg.norm(rvec)
       vvec = X_k[0:2, 2, k].reshape([2, 1])
       v = np.linalg.norm(vvec)
       vhat = vvec / v
       omega = xi[4]
       r_turn = turning_radius(v, omega)
       R = X_k[0:2, 0:2, k].T
       a_centrip = np.sign(xi[4]) * v ** 2 / r_turn * np.array([[0, -1], [1, 0]]) @ R @ vhat
       xi[0:2] = xi[0:2] + a_centrip
   X_k[:, :, k + 1] = se22_Lie_group_integrator(X_k[:, :, k], xi, T[k])
   xi_kc = xi_k[:, :, k]
   X_kc = X_k[:, :, k]
   X_kp1 = X_k[:, :, k + 1]

   hiii = 1
   
for s in range(k_seg_total - 1):
    xi = xi_seg[:, s].reshape([5, 1])
    if xi[4] != 0: # Add centripital acceleration
       rvec = X_seg[0:2, 3, s].reshape([2, 1])
       r = np.linalg.norm(rvec)
       vvec = X_seg[0:2, 2, s].reshape([2, 1])
       v = np.linalg.norm(vvec)
       vhat = vvec / v
       omega = xi[4]
       r_turn = turning_radius(v, omega)
       R = X_seg[0:2, 0:2, s].T
       a_centrip = np.sign(xi[4]) * v ** 2 / r_turn * np.array([[0, -1], [1, 0]]) @ R @ vhat
       xi[0:2] = xi[0:2] + a_centrip
    X_seg[:, :, s + 1] = se22_Lie_group_integrator(X_seg[:, :, s], xi, dt)

# Plot
# Trajectories
plt.plot(X_seg[0, 3, :], X_seg[1, 3, :])
# plt.quiver(X_k[0, 3, :], X_k[1, 3, :], X_k[0, 2, :], X_k[1, 2, :])
# for k in range(N_seg):
#     plt.text(X_k[0, 3, k], X_k[1, 3, k], f"Segment {k}", fontsize=10, 
#               bbox=dict(facecolor='red', alpha=0.5))
plt.xlabel("X Position [m]")
plt.ylabel("Y Position [m]")
plt.title("Flight Profile in (x, y) Plane")
plt.legend()
plt.grid()
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("./HW3/Parts/Q4/AAE590LGM_HW3_Q4_traj.png")
plt.show()

vel_mag = np.linalg.vector_norm(X_seg[0:2, 2, :], axis = 0)
plt.plot(time, vel_mag[0:-1])
plt.grid()
plt.title("Velocity vs Time")
plt.xlabel("Time [s]")
plt.ylabel("Velocity [m / s]")
plt.savefig("./HW3/Parts/Q4/AAE590LGM_HW3_Q4_vel.png")
plt.show()

heading = np.zeros([k_seg_total - 1, 1])
for k in range(k_seg_total - 1):
    xi_recover = se22_log(X_seg[:, :, k])
    heading[k] = xi_recover[4]
plt.plot(time, np.rad2deg(heading))
plt.grid()
plt.title("Heading vs Time")
plt.xlabel("Time [s]")
plt.ylabel("Heading [deg]")
plt.savefig("./HW3/Parts/Q4/AAE590LGM_HW3_Q4_head.png")
plt.show()

plt.plot(time, X_seg[0, 3, 0:-1])
plt.plot(time, X_seg[1, 3, 0:-1])
plt.grid()
plt.title("Position vs Time")
plt.xlabel("Time")
plt.xlabel("Position")
plt.savefig("./HW3/Parts/Q4/AAE590LGM_HW3_Q4_pos.png")
plt.show()

# Simulate spiral
dt_spiral = 0.01 # [s]
tf = 10 # [s]
N_t = np.round(tf / dt_spiral).astype(int)
X_k_spiral = np.zeros([4, 4, N_t + 1])
omega_spiral = 0.5 # [rad / s]
a_spiral = np.array([[1], [0]])
xi_spiral = np.block([[a_spiral], [np.zeros((2, 1))], [omega_spiral]]).reshape([5, 1])
X_k_spiral[:, :, 0] = X_0

for k in range(N_t):
    xi = xi_spiral
    if xi[4] != 0: # Add centripital acceleration
       rvec = X_k_spiral[0:2, 3, k].reshape([2, 1])
       r = np.linalg.norm(rvec)
       vvec = X_k_spiral[0:2, 2, k].reshape([2, 1])
       v = np.linalg.norm(vvec)
       if v < 1e-5:
           vhat = np.zeros([2, 1])
           r_turn = 1
       else:
           vhat = vvec / v
           omega = xi[4]
           r_turn = turning_radius(v, omega)
       R = X_k_spiral[0:2, 0:2, k].T
       a_centrip = np.sign(xi[4]) * v ** 2 / r_turn * np.array([[0, -1], [1, 0]]) @ R @ vhat
       xi[0:2] = xi[0:2] + a_centrip
    X_k_spiral[:, :, k + 1] = se22_Lie_group_integrator(X_k_spiral[:, :, k], xi, dt_spiral)

plt.plot(X_k_spiral[0, 3, :], X_k_spiral[1, 3, :], label = "Line")
v_ind = np.round(np.linspace(0, N_t, 10)).astype(int)
plt.quiver(X_k_spiral[0, 3, v_ind], X_k_spiral[1, 3, v_ind], X_k_spiral[0, 2, v_ind], X_k_spiral[1, 2, v_ind])
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.title("Spiral")
plt.legend()
plt.gca().set_aspect('equal', adjustable='box')
plt.grid()
plt.savefig("./HW3/Parts/Q4/AAE590LGM_HW3_Q4_spiral.png")
plt.show()

vel_mag = np.linalg.vector_norm(X_k_spiral[0:2, 2, :], axis = 0)
plt.plot(np.linspace(0, tf, N_t + 1), vel_mag)
plt.grid()
plt.title("Velocity vs Time for Spiral")
plt.xlabel("Time [s]")
plt.xlabel("Velocity [m / s]")
plt.savefig("./HW3/Parts/Q4/AAE590LGM_HW3_Q4_spiralvel.png")
plt.show()

hi = 1