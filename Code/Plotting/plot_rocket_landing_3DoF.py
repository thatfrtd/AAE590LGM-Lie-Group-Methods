import numpy as np
import matplotlib.pyplot as plt
from Code.SO2.SO2_maps import so2_exp

def plot_rocket_landing_3DoF_trajectory(t, x, u, glideslope_max_angle, max_gimbal_angle, T_min, T_max, pretag):
    
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212, projection = 'polar')
    
    Nu = u.shape[1]

    T = np.zeros_like(u)
    for k in range(len(t) - 1):
        T[:, k] = so2_exp(x[4, k]) @ u[:, k]

    ax1.plot(x[0, :], x[1, :])
    ax1.quiver(x[0, :Nu], x[1, :Nu], T[0, :Nu], T[1, :Nu], color = "red", scale = 1e-1)
    ax1.set_xlabel("X [m]")
    ax1.set_ylabel("Y [m]")
    ax1.set_title("State Trajectory")
    ax1.set_aspect('equal')
    ax1.grid()

    thrust = np.linalg.norm(u, 2)
    gimbal = np.atan2(u[1, :], u[0, :])

    #ax2.polar(gimbal, thrust)
    #ax2.set_title("Control Trajectory")

    plt.suptitle(f"{pretag} 3DoF Rocket Landing Trajectory")
    plt.show()

    return 

def plot_rocket_landing_3DoF_history(t, x, u, max_gimbal_angle, T_min, T_max):
    
    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(3, 2)

    ax1.plot(t, x[0:2, :])
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Position [m]")
    ax1.set_title("Position vs Time")
    ax1.grid(True)

    ax2.plot(t, x[2:4, :])
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Velocity [m / s]")
    ax2.set_title("Velocity vs Time")
    ax2.grid(True)

    ax3.plot(t, np.rad2deg(x[4, :]))
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("Orientation [deg]")
    ax3.set_title("Orientation vs Time")
    ax3.grid(True)

    ax4.plot(t, np.rad2deg(x[5, :]))
    ax4.set_xlabel("Time [s]")
    ax4.set_ylabel("Angular Velocity [deg]")
    ax4.set_title("Angvel vs Time")
    ax4.grid(True)

    ax5.plot(t, np.linalg.norm(u, 2))
    ax5.hline(T_min)
    ax5.hline(T_max)
    ax5.set_xlabel("Time [s]")
    ax5.set_ylabel("Thrust [N]")
    ax5.set_title("Thrust vs Time")
    ax5.grid(True)

    ax6.plot(t, np.atan2(u[1, :], u[0, :]))
    ax6.hline(max_gimbal_angle)
    ax6.set_xlabel("Time [s]")
    ax6.set_ylabel("Gimbal [deg]")
    ax6.set_title("Gimbal Angle vs Time")
    ax6.grid(True)

    return

def comparison_plot_rocket_landing_3DoF_trajectory():

    return