import pytest as pt
import numpy as np
from Code.SO2.SO2_maps import so2_wedge
from Code.SE22.SE22_maps import se22_compose, se22_inverse, se22_exp, se22_log, se22_vee, se22_wedge, se22_Ad, se22_ad

def test_se22_group_affine():
    """
    Test that group affine property holds for random SE2 elements with mixed invariant dynamics f_u(X) = xi_R_wedge X + X xi_L_wedge
    f_u(XY) = f_u(X) Y + X f_u(Y) - X f_u(I) Y
    """
    n_samples = 100

    rng = np.random.default_rng()
    xi_X_random = np.reshape(rng.uniform(-10, 10, 5 * n_samples), (5, n_samples))
    xi_Y_random = np.reshape(rng.uniform(-10, 10, 5 * n_samples), (5, n_samples))
    a = np.reshape(rng.uniform(-10, 10, 2 * n_samples), (2, n_samples))
    omega = rng.uniform(-10, 10, n_samples)

    checks = np.zeros([n_samples, 1])
    for i in range(n_samples):
        X_random = se22_exp(xi_X_random[:, i])
        Y_random = se22_exp(xi_Y_random[:, i])

        M = np.zeros((4, 4))
        N = np.block([[so2_wedge(omega[i]), a[:, i].reshape((2, 1)), np.zeros((2, 1))], [np.zeros((2, 4))]])
        C = np.block([np.zeros((4, 3)), np.array((0, 0, 1, 0)).reshape((4, 1))])
        xi_R_wedge = M - C
        xi_L_wedge = N + C

        lhs = xi_R_wedge @ X_random @ Y_random + X_random @ Y_random @ xi_L_wedge # f_u(XY)
        rhs = (xi_R_wedge @ X_random + X_random @ xi_L_wedge) @ Y_random + X_random @ (xi_R_wedge @ Y_random + Y_random @ xi_L_wedge) - X_random @ (xi_R_wedge + xi_L_wedge) @ Y_random # f_u(X) Y + X f_u(Y) - X f_u(I) Y

        checks[i] = np.all(lhs == pt.approx(rhs))

    assert np.all(checks) 

def test_se22_group_affine_wind():
    """
    Test that group affine property holds for random SE2 elements with mixed invariant dynamics with wind f_u(X) = xi_R_wedge X + X xi_L_wedge
    f_u(XY) = f_u(X) Y + X f_u(Y) - X f_u(I) Y
    """
    n_samples = 100

    rng = np.random.default_rng()
    xi_X_random = np.reshape(rng.uniform(-10, 10, 5 * n_samples), (5, n_samples))
    xi_Y_random = np.reshape(rng.uniform(-10, 10, 5 * n_samples), (5, n_samples))
    a = np.reshape(rng.uniform(-10, 10, 2 * n_samples), (2, n_samples))
    omega = rng.uniform(-10, 10, n_samples)
    
    w_a = np.array((3, -2,)).reshape((2, 1)) # [m / s2]

    checks = np.zeros([n_samples, 1])
    for i in range(n_samples):
        X_random = se22_exp(xi_X_random[:, i])
        Y_random = se22_exp(xi_Y_random[:, i])

        M = np.block([np.zeros((4, 2)), np.vstack((w_a, np.zeros((2, 1)))), np.zeros((4, 1))])
        N = np.block([[so2_wedge(omega[i]), a[:, i].reshape((2, 1)), np.zeros((2, 1))], [np.zeros((2, 4))]])
        C = np.block([np.zeros((4, 3)), np.array((0, 0, 1, 0)).reshape((4, 1))])
        xi_R_wedge = M - C
        xi_L_wedge = N + C

        lhs = xi_R_wedge @ X_random @ Y_random + X_random @ Y_random @ xi_L_wedge # f_u(XY)
        rhs = (xi_R_wedge @ X_random + X_random @ xi_L_wedge) @ Y_random + X_random @ (xi_R_wedge @ Y_random + Y_random @ xi_L_wedge) - X_random @ (xi_R_wedge + xi_L_wedge) @ Y_random # f_u(X) Y + X f_u(Y) - X f_u(I) Y

        checks[i] = np.all(lhs == pt.approx(rhs))

    assert np.all(checks) 