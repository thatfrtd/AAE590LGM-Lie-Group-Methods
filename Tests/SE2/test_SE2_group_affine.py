import pytest as pt
import numpy as np
from Code.SE2.SE2_maps import se2_compose, se2_inverse, se2_exp, se2_log, se2_vee, se2_wedge, se2_Ad, se2_ad

def test_se2_group_affine():
    """
    Test that group affine property holds for random SE2 elements with unicycle dynamics f_u(X) = X xi^
    f_u(XY) = f_u(X) Y + X f_u(Y) - X f_u(I) Y
    """
    n_samples = 100

    rng = np.random.default_rng()
    xi_X_random = np.reshape(rng.uniform(-10, 10, 3 * n_samples), (3, n_samples))
    xi_Y_random = np.reshape(rng.uniform(-10, 10, 3 * n_samples), (3, n_samples))

    xi = np.array([[15], [0], [0.3]])

    checks = np.zeros([n_samples, 1])
    for i in range(n_samples):
        X_random = se2_exp(xi_X_random[:, i])
        Y_random = se2_exp(xi_Y_random[:, i])

        lhs = X_random @ Y_random @ se2_wedge(xi)  # f_u(XY)
        rhs = X_random @ se2_wedge(xi) @ Y_random + X_random @ Y_random @ se2_wedge(xi) - X_random @ se2_wedge(xi) @ Y_random # f_u(X) Y + X f_u(Y) - X f_u(I) Y

        checks[i] = np.all(lhs == pt.approx(rhs))

    assert np.all(checks) 
