import numpy as np
import matplotlib.pyplot as plt
from Code.SO3.SO3_quat_maps import quat_rot

def plot_rocket_landing_6DoF_trajectory(t, x, u, glideslope_max_angle, max_gimbal_angle, T_min, T_max, p_b, r_LoS_min, r_LoS_max, pretag):
    
    ax1 = plt.subplot(111, projection = "3d")
    #ax1 = plt.subplot(211, projection = "3d")
    #ax2 = plt.subplot(212, projection = 'polar')
    
    Nu = u.shape[1]

    T = np.zeros_like(u[0:3, :])
    sensor = np.zeros_like(u[0:3, :])
    for k in range(len(t) - 1):
        T[:, k] = np.array(quat_rot(x[6:10, k], u[0:3, k])).reshape((3,))
        r_mag = np.linalg.norm(x[0:3, k])
        if r_mag > r_LoS_min and r_mag < r_LoS_max:
            sensor[:, k] = np.array(quat_rot(x[6:10, k], p_b)).reshape((3,))

    u_cone, v_cone = np.mgrid[0:2*np.pi:100j, 0:np.pi:80j]
    x_cone = np.cos(u_cone)*np.sin(v_cone)
    y_cone = np.sin(u_cone)*np.sin(v_cone)
    z_cone = np.sqrt(x_cone ** 2 +  y_cone ** 2) * 1 / np.tan(glideslope_max_angle)
    a = np.sin(glideslope_max_angle)

    '''cone_scale = np.max(x[2, :]) / np.max(z_cone)
    x_cone = x_cone * cone_scale
    y_cone = y_cone * cone_scale
    z_cone = z_cone * cone_scale'''

    ax1.plot_surface(x_cone, y_cone, z_cone, alpha = 0.2)

    ax1.plot(x[0, :], x[1, :], x[2, :])
    ax1.quiver(x[0, :Nu], x[1, :Nu], x[2, :Nu], -T[0, :Nu], -T[1, :Nu], -T[2, :Nu], color = "red", length=10, normalize=False)
    ax1.quiver(x[0, :Nu], x[1, :Nu], x[2, :Nu], -sensor[0, :Nu], -sensor[1, :Nu], -sensor[2, :Nu], color = "blue", length=0.2, normalize=False)
    ax1.set_xlabel("X [km]")
    ax1.set_ylabel("Y [km]")
    ax1.set_zlabel("Z [km]")
    #ax1.set_title("State Trajectory")
    ax1.set_aspect('equal')
    ax1.legend(["Glideslope", "Trajectory", "Thrust", "Sensor LoS"])
    ax1.grid()

    thrust = np.linalg.norm(u, 2)
    gimbal = np.atan2(u[1, :], u[0, :])

    #ax2.polar(gimbal, thrust)
    #ax2.set_title("Control Trajectory")

    plt.title(f"{pretag} 6DoF Rocket Landing Trajectory")
    plt.show()

    return 

def plot_rocket_landing_6DoF_history(t, x, u, max_gimbal_angle, T_min, T_max):
    
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

def comparison_plot_rocket_landing_6DoF_trajectory():

    return