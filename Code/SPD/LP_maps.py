import numpy as np
from Code.SSCP.qr_covariance import triu, tril, tril_to_mat

# Group of lower triangular matrices with positive diagonals
# They form Cholesky factors of SPD matrices

def Exp_LP(L, X):
    Exp_X = np.tril(L, -1) + np.tril(X, -1) + np.diag(np.diag(L) * np.exp(np.diag(X) / np.diag(L)))

    return Exp_X

def Log_LP(L, K):
    Log_K = np.tril(K, -1) - np.tril(L, -1) + np.diag(np.diag(L) * np.log(np.diag(K) / np.diag(L)))

    return Log_K 

def compose_LP(L1, L2):
    L12 = np.tril(L1, -1) + np.tril(L2, -1) + np.diag(np.diag(L1) * np.diag(L2))

    return L12

def inverse_LP(L):
    L_inv = np.diag(1 / np.diag(L)) - np.tril(L, -1)

    return L_inv

def err_LP(P, Q):
    eye = np.eye(P.shape[0])
    err = Log_LP(eye, compose_LP(inverse_LP(P), Q))

    return err

def interp_LP(L0, LN, N):
    err = err_LP(L0, LN)
    eye = np.eye(L0.shape[0])
    
    t = np.linspace(0, 1, N)
    L_interp = np.zeros((L0.shape[0], L0.shape[1], N))
    for k in range(N):
        L_interp[:, :, k] = compose_LP(L0, Exp_LP(eye, err * t[k]))

    return L_interp

def interp_LP_vec(L0_vec, LN_vec, N):
    L0 = tril_to_mat(L0_vec)
    LN = tril_to_mat(LN_vec)
    err = err_LP(L0, LN)
    eye = np.eye(L0.shape[0])
    
    t = np.linspace(0, 1, N)
    L_interp = np.zeros((L0_vec.shape[0], N))
    for k in range(N):
        L_interp[:, k] = tril(compose_LP(L0, Exp_LP(eye, err * t[k]))).flatten()

    return L_interp

def geodesic_dist_LP(L, K):
    dist = np.sqrt(np.linalg.norm(np.tril(L, -1) - np.tril(K, -1), 'fro') ** 2 + np.linalg.norm(np.diag(np.log(np.diag(L)) - np.log(np.diag(K))), "fro") ** 2)
    
    return dist