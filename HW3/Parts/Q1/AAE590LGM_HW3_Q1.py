import numpy as np
import matplotlib.pyplot as plt
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_log, se2_wedge, se2_vee, se2_inverse
from Code.SE2.SE2_integrators import se2_euler_integrator, se2_Lie_group_integrator

## 1b.1) Plot turning radius vs velocity
def turning_radius(velocity, desired_turn_rate):
    g = 9.81 # [m / s2]

    turn_radius = velocity ** 2 / (g * np.tan(desired_turn_rate))

    return turn_radius

V = np.linspace(10, 40, 100) # [m / s]
phi_0 = np.deg2rad(30) # Desired turn rate [rad]
turn_radius = turning_radius(V, phi_0)

#plt.plot(V, turn_radius)
#plt.xlabel("Velocity [m / s]")
#plt.ylabel("Turning radius [m]")
#plt.title("Turning Radius vs Velocity for 30 degree Turn Rate")
#plt.grid()
#plt.show()

## 1b.2) Propagate waypoints
X_0 = np.eye(3)

# Segment info
N_seg = 8
v_x = np.array([20, 15, 20, 15, 20, 15, 20, 15]) # [m / s]
v_y = np.zeros_like(v_x)
omega = np.array([0, 0.3, 0, 0.3, 0, 0.3, 0, 0.3]) # [rad / s]
T = np.array([10, 5.24, 5, 5.24, 10, 5.24, 5, 5.24]) # [s]

xi_k = np.array([[v_x], [v_y], [omega]])
dt = 0.04 # [s]

# Simulation loop
X_k = np.zeros([3, 3, N_seg + 1]) # Waypoint group elements
X_k[:, :, 0] = X_0

T_seg = np.cumsum(T)
time = np.arange(0, T_seg[-1], dt)
k_seg = np.round(T / dt).astype(int)
k_seg_total = sum(k_seg) + 1
X_seg = np.zeros([3, 3, k_seg_total])
X_seg[:, :, 0] = X_0
xi_seg = np.hstack((np.matlib.repmat(xi_k[:, :, 0], 1, k_seg[0]),
                    np.matlib.repmat(xi_k[:, :, 1], 1, k_seg[1]),
                    np.matlib.repmat(xi_k[:, :, 2], 1, k_seg[2]),
                    np.matlib.repmat(xi_k[:, :, 3], 1, k_seg[3]),
                    np.matlib.repmat(xi_k[:, :, 4], 1, k_seg[4]),
                    np.matlib.repmat(xi_k[:, :, 5], 1, k_seg[5]),
                    np.matlib.repmat(xi_k[:, :, 6], 1, k_seg[6]),
                    np.matlib.repmat(xi_k[:, :, 7], 1, k_seg[7])))

for k in range(N_seg):
   X_k[:, :, k + 1] = se2_Lie_group_integrator(X_k[:, :, k], xi_k[:, :, k], T[k])
   
for s in range(k_seg_total - 1):
    X_seg[:, :, s + 1] = se2_Lie_group_integrator(X_seg[:, :, s], xi_seg[:, s], dt)

# Plot
# Trajectories
plt.plot(X_seg[0, 2, :], X_seg[1, 2, :])
plt.quiver(X_k[0, 2, :], X_k[1, 2, :], X_k[0, 0, :], X_k[1, 0, :])
for k in range(N_seg):
    plt.text(X_k[0, 2, k], X_k[1, 2, k], f"Segment {k}", fontsize=10, 
              bbox=dict(facecolor='red', alpha=0.5))
plt.xlabel("X Position [m]")
plt.ylabel("Y Position [m]")
plt.title("Reference Trajectory in (x, y) Plane")
plt.legend()
plt.grid()
plt.gca().set_aspect('equal', adjustable='box')
#plt.savefig("./HW2/Parts/Q4/AAE590LGM_HW2_Q4_traj.png")
plt.show()

def required_bank_angle(velocity, omega):
    g = 9.81 # [m / s2]

    bank_angle = np.atan(omega * velocity / g)

    return bank_angle

bank_angles = required_bank_angle(xi_seg[0, :], xi_seg[2, :])
plt.plot(time[0:-1], np.rad2deg(bank_angles))
plt.xlabel("Time [s]")
plt.ylabel("Required Bank Angle [deg]")
plt.title("Required Bank Angle vs Time")
plt.grid()
plt.show()

# Simulate loop and line
dt_loopline = 0.01 # [s]
tf = 10 # [s]
N_t = np.round(tf / dt_loopline).astype(int)
X_k_loop = np.zeros([3, 3, N_t + 1])
X_k_line = np.zeros_like(X_k_loop)
xi_line = np.array([[20], [0], [0]])
xi_loop = np.array([[15], [0], [0.3]])
X_k_loop[:, :, 0] = X_0
X_k_line[:, :, 0] = X_0

for k in range(N_t):
    X_k_loop[:, :, k + 1] = se2_Lie_group_integrator(X_k_loop[:, :, k], xi_loop, dt_loopline)
    a = se2_Lie_group_integrator(X_k_loop[:, :, k], xi_loop, dt)
    X_k_line[:, :, k + 1] = se2_Lie_group_integrator(X_k_line[:, :, k], xi_line, dt_loopline)

plt.plot(X_k_line[0, 2, :], X_k_line[1, 2, :], label = "Line")
plt.plot(X_k_loop[0, 2, :], X_k_loop[1, 2, :], label = "Loop")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.title("Linus and Loopy")
plt.legend()
plt.gca().set_aspect('equal', adjustable='box')
plt.grid()
plt.show()