import numpy as np
import numba as nb
import casadi as cas
from scipy.integrate import solve_ivp

## State based discretization
def discretize_dynamics_ZOH(f, A, B, c, N, tspan, x_ref, u_ref, tolerances):
    '''Discretization of a dynamical system assuming ZOH control'''
    # Make c optional? c = f - Ax - Bu

    nx = x_ref(0).size
    nu = u_ref(0).size

    t_k = np.linspace(tspan[0], tspan[1], N + 1)
    A_k = np.zeros([nx, nx, N])
    B_k = np.zeros([nx, nu, N])
    c_k = np.zeros([nx, 1, N])
    
    for k in range(N):
        A_k[:, :, k], B_k[:, :, k], c_k[:, :, k] = integrate_discrete_ZOH(x_ref(t_k[k]), A, B, c, f, u_ref, (t_k[k], t_k[k + 1]), tolerances)

    return A_k, B_k, c_k

def integrate_discrete_ZOH(x0, A, B, c, f, u, tspan, tolerances):
    '''Integrates STM and state with Bk and ck'''
    #   Uses scipy.solve_ivp to integrate the State Transition Matrix and the state using
    #   the given A matrix and dynamics f over the time period in tspan using
    #   the specified error tolerance.

    # Create initial condition
    nx = x0.size

    STM0 = np.eye(nx).reshape(-1, 1)
    B0 = np.zeros_like(B(0, x0, u(0))).reshape(-1, 1)
    c0 = np.zeros([nx, 1])

    nu = u(0).size

    y0 = np.vstack([x0, STM0, B0, c0])
     
    # Simulate    
    sol = solve_ivp(STM_diff_eq_ZOH, tspan, y0.flatten(), t_eval = [tspan[1]], method="RK45", rtol = tolerances["rtol"], atol = tolerances["atol"], args = (A, B, c, f, u, nx))

    # Unpack solution
    y = sol.y
    x = y[0:nx]
    A_k = np.reshape(y[nx : (nx * (nx + 1))], (nx, nx))
    B_k = A_k @ np.reshape(y[(nx * (nx + 1)) : (nx * (nx + 1) + nx * nu)], (nx, nu))
    c_k = A_k @ y[(nx * (nx + 1) + nx * nu) : (nx * (nx + 1) + nx * nu + nx)]

    return A_k, B_k, c_k

def STM_diff_eq_ZOH(t, y, A, B, c, f, u, n):
    x = y[0:n][:, np.newaxis]
    STM = np.reshape(y[n : (n * (n + 1))], (n, n))

    STM_inverse = np.linalg.inv(STM)

    xdot = f(t, x, u(t))
    A_kdot = A(t, x, u(t)) @ STM
    B_kdot = STM_inverse @ B(t, x, u(t))
    c_kdot = STM_inverse @ c(t, x, u(t))

    A_kdot_flat = A_kdot.flatten()[:, np.newaxis]
    B_kdot_flat = B_kdot.flatten()[:, np.newaxis]

    ydot = np.vstack((xdot, A_kdot_flat, B_kdot_flat, c_kdot))

    return ydot.flatten()

## Error state based discretization

def STM_diff_eq_ZOH(y, A, B, f, u, nx, nu):
    # Extract variables
    x = y[0 : nx]
    Phi_A = cas.reshape(y[nx : (nx * (nx + 1))], (nx, nx)) # STM
    Phi_B = cas.reshape(y[(nx * (nx + 1)) : (nx * (nx + 1) + nx * nu)], (nx, nu))

    # Calculate Jacobians
    A_t = A(x, u)
    B_t = B(x, u)

    # Calculate derivatives
    x_dot = f(x, u)
    A_k_dot = A_t @ Phi_A
    B_k_dot = A_t @ Phi_B + B_t

    # Package derivatives
    y_dot = cas.vertcat(x_dot, 
                        cas.reshape(A_k_dot, (nx * nx, 1)), 
                        cas.reshape(B_k_dot, (nx * nu, 1)))

    return y_dot

def integrate_error_discrete_ZOH(x0, u, delta_t, A, B, f, nx, nu, tolerance):
    # Integrates STM and state with B, and ck
    #   Uses CVODES to integrate the State Transition Matrix and the state using
    #   the given A matrix and dynamics f over the time period in tspan using
    #   the specified error tolerance.

    # Create initial condition
    STM0 = cas.SX.eye(nx)
    B0 = cas.SX.zeros(nx, nu)    

    y0 = cas.vertcat(x0, cas.reshape(STM0, (nx * nx, 1)), cas.reshape(B0, (nx * nu, 1)))

    # Create integrator
    y = cas.SX.sym('y', (nx + nx * nx + nx * nu))
    p = cas.vertcat(delta_t, u)
    opts = {'reltol':tolerance, 'abstol':tolerance}
    ivp = {'x':y, 'p':p, 'ode':STM_diff_eq_ZOH(y, A, B, f, u, nx, nu) * delta_t}
    Fd = cas.integrator('Fd', 'cvodes', ivp, opts)

    # Solve
    sol = Fd(x0 = y0, p = p)

    y_f = sol["xf"]

    # Unpack solution
    x_kp1 = y_f[0:nx, :]

    A_k = y_f[nx : (nx * (nx + 1)), :]
    B_k = y_f[nx * (nx + 1) : (nx * (nx + 1) + nx * nu)]

    A_k_matrix = cas.reshape(A_k, (nx, nx))
    B_k_matrix = cas.reshape(B_k, (nx, nu))
    
    c_k = x_kp1 - (A_k_matrix @ x0 + B_k_matrix @ u)

    return x_kp1, A_k, B_k, c_k

def discretize_error_dynamics_ZOH(A, B, f, nx, nu, N, tolerance, num_threads):

    x_ref_k = cas.SX.sym("x_ref_k", nx, 1)
    x_ref_kp1 = cas.SX.sym("x_ref_kp1", nx, 1)
    u_ref_k = cas.SX.sym('u_ref_k', nu, 1)
    delta_t = cas.SX.sym('delta_t')

    x_kp1, A_k, B_k, c_k = integrate_error_discrete_ZOH(x_ref_k, u_ref_k, delta_t, A, B, f, nx, nu, tolerance)    

    Delta_k = x_kp1 - x_ref_kp1

    integrate_error_ZOH = cas.Function("integrate_error_ZOH", 
                                       [x_ref_k, x_ref_kp1, u_ref_k, delta_t], [A_k, B_k, c_k, Delta_k], 
                                       ["x_ref_k", "x_ref_kp1", "u_ref_k", "delta_t"], ["A_k", "B_k", "c_k", "Delta_k"])
    
    integrate_error_func = integrate_error_ZOH.map(N - 1, "thread", num_threads)
    
    x_ref = cas.SX.sym("x_ref", nx, N)
    u_ref = cas.SX.sym("u_ref", nu, N - 1)
    t_k = cas.SX.sym("t_k", 1, N)

    integ_map = integrate_error_func(x_ref_k = x_ref[:, 0:-1], x_ref_kp1 = x_ref[:, 1:], u_ref_k = u_ref, delta_t = t_k[1:] - t_k[0:-1])
    
    opts = {'jit':False}
    discretize_error_func = cas.Function("discretize_error", 
                                         [x_ref, u_ref, t_k], [integ_map["A_k"], integ_map["B_k"], integ_map["c_k"], integ_map["Delta_k"]], 
                                         ["x_ref", "u_ref", "t_k"], ["A_ks", "B_ks", "c_ks", "Delta_ks"], opts)#.expand()

    return discretize_error_func


def discrete_ZOH(f_func, x0, u, delta_t, tolerance, jit):

    p = cas.vertcat(delta_t, u)
    opts = {'reltol':tolerance, 'abstol':tolerance, 'jit':False}
    ivp = {'x':x0, 'p':p, 'ode':f_func(x0, u) * delta_t}
    F = cas.integrator('F', 'cvodes', ivp, opts)

    sol = F(x0 = x0, p = p)

    x_kp1 = sol["xf"]
    A_k = cas.jacobian(x_kp1, x0)
    B_k = cas.jacobian(x_kp1, u)
    c_k = x_kp1 - (A_k @ x0 - B_k @ u) 

    return x_kp1, A_k, B_k, c_k

def discretize_dynamics_ZOH(f, nx, nu, N, tolerance, num_threads):

    x_ref_k = cas.SX.sym("x_ref_k", nx, 1)
    x_ref_kp1 = cas.SX.sym("x_ref_kp1", nx, 1)
    u_ref_k = cas.SX.sym('u_ref_k', nu, 1)
    delta_t = cas.SX.sym('delta_t')

    jit = False
    x_kp1, A_k, B_k, c_k = discrete_ZOH(f, x_ref_k, u_ref_k, delta_t, tolerance, jit)
    
    Delta_k = x_kp1 - x_ref_kp1

    integrate_error_ZOH = cas.Function("integrate_error_ZOH", 
                                       [x_ref_k, x_ref_kp1, u_ref_k, delta_t], [A_k, B_k, c_k, Delta_k], 
                                       ["x_ref_k", "x_ref_kp1", "u_ref_k", "delta_t"], ["A_k", "B_k", "c_k", "Delta_k"])
    
    integrate_error_func = integrate_error_ZOH.map(N - 1, "thread", num_threads)
    
    x_ref = cas.SX.sym("x_ref", nx, N)
    u_ref = cas.SX.sym("u_ref", nu, N - 1)
    t_k = cas.SX.sym("t_k", 1, N)

    integ_map = integrate_error_func(x_ref_k = x_ref[:, 0:-1], x_ref_kp1 = x_ref[:, 1:], u_ref_k = u_ref, delta_t = t_k[1:] - t_k[0:-1])
    
    opts = {'jit':jit}
    discretize_error_func = cas.Function("discretize_error", 
                                         [x_ref, u_ref, t_k], [integ_map["A_k"], integ_map["B_k"], integ_map["c_k"], integ_map["Delta_k"]], 
                                         ["x_ref", "u_ref", "t_k"], ["A_ks", "B_ks", "c_ks", "Delta_ks"], opts)#.expand()

    return discretize_error_func
