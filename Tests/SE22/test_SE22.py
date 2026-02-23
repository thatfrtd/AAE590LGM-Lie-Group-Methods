import pytest as pt
import numpy as np
from Code.SE22.SE22_maps import se22_compose, se22_inverse, se22_exp, se22_log, se22_vee, se22_wedge, se22_Ad, se22_ad

def test_se22_exp_log():
    """
    Test that exp(log(X)) = X for random SE22 elements
    """
    n_samples = 100

    rng = np.random.default_rng()
    xi_random = np.reshape(rng.uniform(-10, 10, 5 * n_samples), (5, n_samples))

    checks = np.zeros([n_samples, 1])
    for i in range(n_samples):
        X_random = se22_exp(xi_random[:, i])

        X_log = se22_log(X_random)
        X_tilde = se22_exp(X_log)

        checks[i] = np.all(X_tilde == pt.approx(X_random))

    assert np.all(checks) 

def test_se22_Ad():
    """
    Test that Ad_X(xi) = (X xi^wedge X^-1)^vee for random SE22, se22 elements
    """
    n_samples = 100

    rng = np.random.default_rng()
    xi_X_random = np.reshape(rng.uniform(-10, 10, 5 * n_samples), (5, n_samples))
    xi_random = np.reshape(rng.uniform(-10, 10, 5 * n_samples), (5, n_samples))

    checks = np.zeros([n_samples, 1])
    for i in range(n_samples):
        X_random = se22_exp(xi_X_random[:, i])

        xi_prime = se22_vee(X_random @ se22_wedge(xi_random[:, i]) @ se22_inverse(X_random))
        xi_prime_ad = se22_Ad(X_random) @ xi_random[:, i]

        checks[i] = np.all(xi_prime_ad == pt.approx(xi_prime))

    assert np.all(checks) 

def test_se22_Ad_composition():
    """
    Test that Ad_{X_1 X_2} = Ad_{X_1}Ad_{X_2} for random SE22 elements
    """
    n_samples = 100

    rng = np.random.default_rng()
    xi1_random = np.reshape(rng.uniform(-10, 10, 5 * n_samples), (5, n_samples))
    xi2_random = np.reshape(rng.uniform(-10, 10, 5 * n_samples), (5, n_samples))

    checks = np.zeros([n_samples, 1])
    for i in range(n_samples):
        X1_random = se22_exp(xi1_random[:, i])
        X2_random = se22_exp(xi2_random[:, i])

        Ad_X1X2 = se22_Ad(se22_compose(X1_random, X2_random))
        Ad_X1_Ad_X2 = se22_Ad(X1_random) @ se22_Ad(X2_random)

        checks[i] = np.all(Ad_X1X2 == pt.approx(Ad_X1_Ad_X2))

    assert np.all(checks) 

def test_se22_Ad_inv():
    """
    Test that (Ad_X)^-1 = Ad_{X^-1} for random SE22 elements
    """
    n_samples = 100

    rng = np.random.default_rng()
    xi_random = np.reshape(rng.uniform(-10, 10, 5 * n_samples), (5, n_samples))

    checks = np.zeros([n_samples, 1])
    for i in range(n_samples):
        X_random = se22_exp(xi_random[:, i])

        Ad_X_inv = np.linalg.inv(se22_Ad(X_random))
        Ad_Xinv = se22_Ad(se22_inverse(X_random))

        a = se22_Ad(X_random) @ Ad_X_inv
        b = se22_Ad(X_random) @ Ad_Xinv

        checks[i] = np.all(Ad_X_inv == pt.approx(Ad_Xinv))

    assert np.all(checks) 

def test_se22_ad():
    """
    Test that ad_{xi_1}(xi_2) = [xi_1, xi_2] for random se22 elements
    """
    n_samples = 100

    rng = np.random.default_rng()
    xi1_random = np.reshape(rng.uniform(-10, 10, 5 * n_samples), (5, n_samples))
    xi2_random = np.reshape(rng.uniform(-10, 10, 5 * n_samples), (5, n_samples))

    checks = np.zeros([n_samples, 1])
    for i in range(n_samples):
        Lie_bracket = se22_vee(se22_wedge(xi1_random[:, i]) @ se22_wedge(xi2_random[:, i])
                            - se22_wedge(xi2_random[:, i]) @ se22_wedge(xi1_random[:, i]))
        Lie_bracket_ad = se22_ad(xi1_random[:, i]) @ xi2_random[:, i]

        checks[i] = np.all(Lie_bracket_ad == pt.approx(Lie_bracket))

    assert np.all(checks) 