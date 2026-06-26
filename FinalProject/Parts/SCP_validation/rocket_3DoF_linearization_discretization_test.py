import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sympy as sym
import casadi as cas
import CyRK as cy
import time
from Code.SE2.SE2_maps import se2_compose, se2_exp, se2_log, se2_wedge, se2_vee, se2_inverse
from Code.Dynamics.rocket_landing_3DoF_dynamics import rocket_landing_3DoF_dynamics
from Code.Discretization.discretize_ZOH import discretize_error_dynamics_ZOH, integrate_error_discrete_ZOH, discretize_dynamics_ZOH

def integrate(f_func, nx, nu, tolerance):

    x = cas.SX.sym('x', nx, 1)
    u = cas.SX.sym('u', nu, 1)
    tf = cas.SX.sym('tf')
    p = cas.vertcat(tf, u)
    opts = {'reltol':tolerance, 'abstol':tolerance, 'jit':False}
    ivp = {'x':x, 'p':p, 'ode':f_func(x, u) * tf}
    F = cas.integrator('F', 'cvodes', ivp, opts)

    return F

print(str(type([])))

# Constants
m_0 = 1 # [kg]
I_b = 1 # [kg m2]
g = 9.81 # [m / s2] Gravity
lever = -1 # [m] distance from CoM to rocket engine
Isp = 800 # [s]
alpha = 1 / (g * Isp)
c = np.array([[g], [lever], [I_b], [alpha]])

# Initial conditions
r_0 = np.array([20, 100]).reshape((2, 1))
v_0 = np.array([-10, -20]).reshape((2, 1))
theta_0 = np.array(np.pi/2).reshape((1, 1))
w_0 = np.array(0).reshape((1, 1))

x_0 = np.vstack((r_0, v_0, theta_0, w_0, m_0))

# Terminal conditions
r_f = np.array([0, 0]).reshape((2, 1))
v_f = np.array([0, 0]).reshape((2, 1))
theta_f = np.array(np.pi/2).reshape((1, 1))
w_f = np.array(0).reshape((1, 1))

x_f = np.vstack((r_f, v_f, theta_f, w_f, m_0 - 0.1))

# Initialize problem parameters
nx = 7
nu = 2
N = 40
Nu = N - 1 # ZOH
tf = 35
t_k = np.linspace(0, tf, N)
delta_t = t_k[1] - t_k[0]

u = np.array([20, 0]).reshape((2, 1))

jit = False

# Test dynamics
cont_dyn = rocket_landing_3DoF_dynamics(g, lever, I_b, alpha, jit)
f_func = cont_dyn["f_func"]
A_func = cont_dyn["A_func"]
B_func = cont_dyn["B_func"]

integrate_func = integrate(f_func, nx, nu, 1e-8)
r = integrate_func(x0 = x_0, p = np.vstack((10, u)))

print(r["xf"])

integrate_func_vec = integrate_func.map(N, "thread", 2)

start = time.perf_counter()

a = integrate_func_vec(x0 = np.ones((7, N)), p = np.ones((3, N)))

end = time.perf_counter()
print(f"Integrate Elapsed: {(end - start) * 1e3:.3f} ms")

a = integrate_error_discrete_ZOH(cas.SX.sym("x0", nx, 1), cas.SX.sym("u", nu, 1), cas.SX.sym("dt"), A_func, B_func, f_func, nx, nu, 1e-8)

discretize_error_func = discretize_error_dynamics_ZOH(A_func, B_func, f_func, nx, nu, N, 1e-8, 2)

start = time.perf_counter()

A_ks, B_ks, c_ks, Delta_ks = discretize_error_func(np.ones([nx, N]) * x_0, np.ones([nu, Nu]), t_k)

end = time.perf_counter()
print(f"Error Elapsed: {(end - start) * 1e3:.3f} ms")

discretize_func = discretize_dynamics_ZOH(f_func, nx, nu, N, 1e-8, 2)

start = time.perf_counter()

A_ks_ad, B_ks_ad, c_ks_ad, Delta_ks_ad = discretize_func(np.ones([nx, N]) * x_0, np.ones([nu, Nu]), t_k)


print(A_ks_ad.shape)
print(B_ks_ad.shape)
print(c_ks_ad.shape)
print(Delta_ks_ad.shape)


end = time.perf_counter()
print(f"CVODES Elapsed: {(end - start) * 1e3:.3f} ms")

print(Delta_ks.shape)

a = np.array(f_func(x_0, u))
print(a)

start = time.perf_counter()
sol = sp.integrate.solve_ivp(lambda t, x : np.array(f_func(x, u)).flatten(), [t_k[0], t_k[-1]], x_0.flatten(), t_eval = t_k, rtol = 1e-8, atol = 1e-8)
y_scipy = sol.y
end = time.perf_counter()
print(f"Elapsed: {(end - start) * 1e3:.3f} ms")


x_k_guess = x_0 + (x_f - x_0) * (t_k) / tf

A_ks, B_ks, c_ks, Delta_ks = discretize_error_func(x_k_guess, u * np.ones([nu, Nu]), t_k)

#print(np.array(Delta_ks))
start = time.perf_counter()
a = integrate_func_vec(x0 = np.ones((7, N)) * x_0, p = np.vstack((t_k, np.ones((nu, N)) * u)))
x_k_cas = np.array(a["xf"])
end = time.perf_counter()
print(f"Elapsed: {(end - start) * 1e3:.3f} ms")

x_k = np.zeros([nx, N])
x_k[:, 0] = x_0.flatten()
for k in range(N - 1):
    A_k = np.array(A_ks[:, k]).reshape((nx, nx)).T
    B_k = np.array(B_ks[:, k]).reshape((nu, nx)).T
    c_k = np.array(c_ks[:, k]).reshape((nx,))

    x_k[:, k + 1] = A_k @ x_k[:, k] + B_k @ u.flatten() + c_k

x_k_ad = np.zeros([nx, N])
x_k_ad[:, 0] = x_0.flatten()
for k in range(N - 1):
    A_k = np.array(A_ks_ad[:, k]).reshape((nx, nx)).T
    B_k = np.array(B_ks_ad[:, k]).reshape((nu, nx)).T
    c_k = np.array(c_ks_ad[:, k]).reshape((nx,))

    x_k_ad[:, k + 1] = A_k @ x_k_ad[:, k] + B_k @ u.flatten() + c_k


print(B_k)

plt.plot(y_scipy[0,:], y_scipy[1,:])
plt.plot(x_k[0, :], x_k[1, :])
plt.plot(x_k_cas[0, :], x_k_cas[1, :])
plt.legend(["scipy", "STM", "CasADi"])
plt.show()

plt.plot(x_k[0, :] - y_scipy[0,:], x_k[1, :] - y_scipy[1,:])
plt.plot(x_k_ad[0, :] - y_scipy[0,:], x_k_ad[1, :] - y_scipy[1,:])
plt.plot(x_k_cas[0, :] - y_scipy[0,:], x_k_cas[1, :] - y_scipy[1,:])
plt.legend(["STM", "AD", "Cas"])
plt.show()

#print(np.array(A_ks[:, 0]).reshape((nx, nx)))

#print(np.array(B_ks[:, 0]).reshape((nx, nu)))
#print(np.array(c_ks[:, 0]).reshape((nx, 1)))
#discretize_error_func.generate('gen.c')
#print(open('gen.c','r').read())
