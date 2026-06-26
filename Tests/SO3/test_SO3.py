import pytest as pt
import numpy as np
from Code.SO3.SO3_maps import so3_wedge, so3_vee, so3_exp, so3_log, so3_ljac, so3_inv_ljac

def test_so3_exp_log():
    """
    Test that exp(log(R)) = R for random rotation matrices
    """
    rng = np.random.default_rng()
    theta_random = rng.uniform(-10, 10, 3 * 100).reshape((3, 100))

    checks = np.zeros([theta_random.shape[1], 1])
    for i in range(theta_random.shape[1]):
        R_random = so3_exp(theta_random[:, i])
        checks[i] = np.all(so3_exp(so3_log(R_random)) == pt.approx(R_random))

    assert np.all(checks) 

def test_so3_log_exp():
    """
    Test that log(exp(theta)) = theta for all theta in (-pi, pi]
    """
    theta_array = np.linspace(-np.pi, np.pi, 3 * 100).reshape((3, 100))

    checks = np.zeros([theta_array.shape[1], 1])
    for i in range(theta_array.shape[1]):
        a = so3_log(so3_exp(theta_array[:, i]))
        b = theta_array[:, i]
        checks[i] = np.all(so3_log(so3_exp(theta_array[:, i])) == pt.approx(theta_array[:, i]))

    assert np.all(checks) 

def test_so3_constraint():
    '''
    Test that R * R.T = I
    '''
    theta_array = np.linspace(-np.pi, np.pi, 3 * 100).reshape((3, 100))

    checks = np.zeros([theta_array.shape[1], 1])
    for i in range(theta_array.shape[1]):
        checks[i] = np.all(so3_exp(theta_array[:, i]).T @ so3_exp(theta_array[:, i]) == pt.approx(np.eye(3)))

    assert np.all(checks) 
#test_so3_log_exp()

def test_so3_jacs():
    """
    Test that ljac(xi) @ ljac_inv(xi) = I for all theta in (-pi, pi]
    """                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    n_samples = 100

    rng = np.random.default_rng()
    xi_random = np.reshape(rng.uniform(-10, 10, 3 * n_samples), (3, n_samples))

    checks = np.zeros([n_samples, 1])
    for i in range(n_samples):
        checks[i] = np.all(so3_ljac(xi_random[:, i]) @ so3_inv_ljac(xi_random[:, i]) == pt.approx(np.eye(3)))

    assert np.all(checks)

def test_so2_inv_jac():
    """
    Test that ljac_inv(xi) = inv(ljac(xi)) for all theta in (-pi, pi]
    """
    n_samples = 100

    rng = np.random.default_rng()
    xi_random = np.reshape(rng.uniform(-10, 10, 3 * n_samples), (3, n_samples))

    checks = np.zeros([n_samples, 1])
    for i in range(n_samples):
        checks[i] = np.all(so3_inv_ljac(xi_random[:, i]) == pt.approx(np.linalg.inv(so3_ljac(xi_random[:, i]))))

    assert np.all(checks)