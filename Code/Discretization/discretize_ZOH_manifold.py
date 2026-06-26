import casadi as cas

def STM_diff_eq_ZOH_manifold(y, A_err, B_err, f, u_ref, nx, neps, nu):
    # Extract variables
    x_ref = y[0 : nx]
    Phi_A = cas.reshape(y[(nx) : (nx + neps * neps)], (neps, neps)) # STM
    Phi_B = cas.reshape(y[(nx + neps * neps) : (nx + neps * neps + neps * nu)], (neps, nu))

    # Calculate Jacobians
    A_t = A_err(x_ref, u_ref)
    B_t = B_err(x_ref, u_ref)

    # Calculate derivatives
    x_ref_dot = f(x_ref, u_ref)
    A_k_dot = A_t @ Phi_A
    B_k_dot = A_t @ Phi_B + B_t

    # Package derivatives
    y_dot = cas.vertcat(x_ref_dot, 
                        cas.reshape(A_k_dot, (neps * neps, 1)), 
                        cas.reshape(B_k_dot, (neps * nu, 1)))

    return y_dot

def integrate_error_discrete_ZOH_manifold(x0_ref, u_ref, delta_t, A_err, B_err, f, nx, neps, nu, tolerance):
    # Integrates STM and state with B_err
    #   Uses CVODES to integrate the State Transition Matrix and the state using
    #   the given A_err matrix and dynamics f over the time period in tspan using
    #   the specified error tolerance.

    # Create initial condition
    STM0 = cas.SX.eye(neps)
    B0 = cas.SX.zeros(neps, nu)    

    y0 = cas.vertcat(x0_ref, cas.reshape(STM0, (neps * neps, 1)), cas.reshape(B0, (neps * nu, 1)))

    # Create integrator
    y = cas.SX.sym('y', (nx + neps * neps + neps * nu))
    p = cas.vertcat(delta_t, u_ref)
    opts = {'reltol':tolerance, 'abstol':tolerance}
    ivp = {'x':y, 'p':p, 'ode':STM_diff_eq_ZOH_manifold(y, A_err, B_err, f, u_ref, nx, neps, nu) * delta_t}
    Fd = cas.integrator('Fd', 'cvodes', ivp, opts)

    # Solve
    sol = Fd(x0 = y0, p = p)

    y_f = sol["xf"]

    # Unpack solution
    x_kp1 = y_f[0:nx, :]

    A_k = y_f[nx : (nx + neps * neps), :]
    B_k = y_f[(nx + neps * neps) : (nx + neps * neps + neps * nu)]

    return x_kp1, A_k, B_k

def discretize_error_dynamics_ZOH_manifold(A_err, B_err, f, nx, neps, nu, N, err_func, tolerance, num_threads):

    x_ref_k = cas.SX.sym("x_ref_k", nx, 1)
    x_ref_kp1 = cas.SX.sym("x_ref_kp1", nx, 1)
    u_ref_k = cas.SX.sym('u_ref_k', nu, 1)
    delta_t = cas.SX.sym('delta_t')

    x_ref_kp1_prop, A_k, B_k = integrate_error_discrete_ZOH_manifold(x_ref_k, u_ref_k, delta_t, A_err, B_err, f, nx, neps, nu, tolerance)    
    
    Delta_k = -err_func(x_ref_kp1_prop, x_ref_kp1) # Same as c_k for error state

    integrate_error_ZOH = cas.Function("integrate_error_ZOH_manifold", 
                                       [x_ref_k, x_ref_kp1, u_ref_k, delta_t], [A_k, B_k, Delta_k], 
                                       ["x_ref_k", "x_ref_kp1", "u_ref_k", "delta_t"], ["A_k", "B_k", "Delta_k"])
    
    integrate_error_func = integrate_error_ZOH.map(N - 1, "thread", num_threads)
    
    x_ref = cas.SX.sym("x_ref", nx, N)
    u_ref = cas.SX.sym("u_ref", nu, N - 1)
    t_k = cas.SX.sym("t_k", 1, N)

    integ_map = integrate_error_func(x_ref_k = x_ref[:, 0:-1], x_ref_kp1 = x_ref[:, 1:], u_ref_k = u_ref, delta_t = t_k[1:] - t_k[0:-1])
    
    opts = {'jit':False}
    discretize_error_func = cas.Function("discretize_error_manifold", 
                                         [x_ref, u_ref, t_k], [integ_map["A_k"], integ_map["B_k"], integ_map["Delta_k"]], 
                                         ["x_ref", "u_ref", "t_k"], ["A_ks", "B_ks", "Delta_ks"], opts)#.expand()

    return discretize_error_func
