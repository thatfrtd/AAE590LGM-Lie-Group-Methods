import pytest as pt
import numpy as np
from Code.SE3.SE3_maps import se3_compose, se3_inverse, se3_exp, se3_log, se3_vee, se3_wedge, se3_Ad, se3_ad

def test_se3_group_affine():
    """
    Test that group affine property holds for random SE3 elements with unicycle dynamics f_u(X) = X xi^
    f_u(XY) = f_u(X) Y + X f_u(Y) - X f_u(I) Y
    """
    n_samples = 100

    rng = np.random.default_rng()
    xi_X_random = np.reshape(rng.uniform(-10, 10, 6 * n_samples), (6, n_samples))
    xi_Y_random = np.reshape(rng.uniform(-10, 10, 6 * n_samples), (6, n_samples))

    xi = np.array([[15], [0], [-2], [0.1], [-0.2], [0.3]])

    checks = np.zeros([n_samples, 1])
    for i in range(n_samples):
        X_random = se3_exp(xi_X_random[:, i])
        Y_random = se3_exp(xi_Y_random[:, i])

        lhs = X_random @ Y_random @ se3_wedge(xi)  # f_u(XY)
        rhs = X_random @ se3_wedge(xi) @ Y_random + X_random @ Y_random @ se3_wedge(xi) - X_random @ se3_wedge(xi) @ Y_random # f_u(X) Y + X f_u(Y) - X f_u(I) Y

        checks[i] = np.all(lhs == pt.approx(rhs))

    assert np.all(checks) 
