import pytest as pt
import numpy as np
from Code.SO3.SO3_maps import so3_wedge, so3_vee, so3_exp, so3_log

def test_so3_exp_log():
    """
    Test that exp(log(R)) = R for random rotation matrices
    """
    rng = np.random.default_rng()
    theta_random = rng.uniform(-10, 10, 3 * 100).reshape((3, 100))

    checks = np.zeros([theta_random.size, 1])
    for i in range(theta_random.size):
        R_random = so3_exp(theta_random[:, i])
        checks[i] = np.all(so3_exp(so3_log(R_random)) == pt.approx(R_random))

    assert np.all(checks) 

def test_so3_log_exp():
    """
    Test that log(exp(theta)) = theta for all theta in (-pi, pi]
    """
    theta_array = -np.linspace(-np.pi, np.pi, 3 * 100).reshape((3, 100))

    checks = np.zeros([theta_array.size, 1])
    for i in range(theta_array.size):
        checks[i] = np.all(so3_log(so3_exp(theta_array[:, i])) == pt.approx(theta_array[:, i]))

    assert np.all(checks) 

def test_so3_commutativity():
    """
    Test that R(theta_1)R(theta_2) = R(theta_1 + theta_2)
    """
    iter = 100

    rng = np.random.default_rng()
    theta_1 = rng.uniform(-10, 10, 3 * iter).reshape((3, 100))
    theta_2 = rng.uniform(-10, 10, 3 * iter).reshape((3, 100))

    checks = np.zeros([iter, 1])
    for i in range(iter):
        checks[i] = so3_exp(theta_1[:, i]) @ so3_exp(theta_2[:, i]) == pt.approx(so3_exp(theta_1[:, i] + theta_2[:, i]))

    assert np.all(checks) 