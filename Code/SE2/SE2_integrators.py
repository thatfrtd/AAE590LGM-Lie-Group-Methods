import numpy as np
from scipy.integrate import solve_ivp
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_wedge, se2_ad
from Code.Dynamics.EulerPoincare import euler_poincare

def se2_euler_integrator(X_k, xi_k, dt):
    """
    Integrate using first order Euler

    :param X_k: group element at t_k
    :param xi_k: twist at t_k
    :param dt: timestep
    """
    X_kp1 = se2_compose(X_k, np.eye(3) + dt * se2_wedge(xi_k)) 

    return X_kp1

def se2_Lie_group_integrator(X_k, xi_k, dt):
    """
    Integrate using exponential map (assuming twist is constant accross timestep)

    :param X_k: group element at t_k
    :param xi_k: twist at t_k
    :param dt: timestep
    """
    X_kp1 = se2_compose(X_k, se2_exp(dt * xi_k)) 

    return X_kp1

def se2_kinodynamic_integrator(X_k, xi_k, J_b, u, dt):
    """
    Integrate initial value problem using Euler-Poincare dynamics and reconstruction equation
    """
    def kinodynamic_equation(t, y, J_b, u):
        X = y[0:9].reshape([3, 3])
        xi = y[9:12]
        
        X_dot = X @ se2_wedge(xi)
        xi_dot = euler_poincare(xi, J_b, se2_ad(xi), u)

        ydot = np.vstack((X_dot.reshape((-1, 1)), xi_dot.reshape((-1, 1)))).reshape((-1,))
        return ydot
    
    y0 = np.vstack((X_k.reshape([-1, 1]), xi_k.reshape([-1, 1]))).reshape((-1,))
    sol = solve_ivp(lambda t, y : kinodynamic_equation(t, y, J_b, u), [0, dt], y0, rtol = 1e-6, atol = 1e-6)
    X_kp1 = sol.y[0:9, -1].reshape((3, 3))
    xi_kp1 = sol.y[9:12, -1]

    return X_kp1, xi_kp1