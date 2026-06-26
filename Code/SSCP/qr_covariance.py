import numpy as np
import casadi as cas
import cvxpy as cp

def qr_pos_diag(M):
    '''
    QR decomposition of matrix with enforcement of positive diagonal condition
    '''
    Q, R = np.linalg.qr(M)

    # Fix negative diagonal elemnts of R while maintaining decomposition
    for k in range(min(M.shape[0], M.shape[1])):
        if R[k, k] < 0:
            R[k, :] = -R[k, :]
            Q[:, k] = -Q[:, k]

    return Q, R

def qr_derivative(dM, Q, R):
    '''
    Derivative of R matrix from QR decomposition w.r.t. input M matrix

    :param dM: M - M_ref
    :param Q: Q_ref matrix from QR decomposition
    :param R: R_ref matrix from QR decomposition
    '''
    A = Q.T @ dM @ np.linalg.inv(R)
    B = tril_to_mat_cp(tril_cp(A))
    Q_TdQ = B - B.T # Anti-symmetric
    dR = Q.T @ dM - Q_TdQ @ R

    return dR

def qr_derivative_wrong(dM, Q, R):
    '''
    Derivative of R matrix from QR decomposition w.r.t. input M matrix
    From Naoya's paper but wrong... (correct version is in his code)
    
    :param dM: M - M_ref
    :param Q: Q_ref matrix from QR decomposition
    :param R: R_ref matrix from QR decomposition
    '''
    dR = Q.T @ dM - tril_to_mat_cp(tril_cp(Q.T @ dM @ np.linalg.inv(R))) @ R

    return dR

def X_Psqrt_noG(A, S, B, L):
    X = A @ S + B @ L

    return X

def X_Psqrt_cp(A, S, B, L, G):
    X = cp.hstack([A @ S + B @ L, G])

    return X

def calculate_X_S_defect_noG(A_k, B_k, S_k, L_k):
    Nm1 = A_k.shape[2]

    X_kp1 = np.zeros_like(A_k)
    defect_S_k = np.zeros(Nm1)
    for k in range(Nm1):
        X_kp1[:, :, k] = X_Psqrt_noG(A_k[:, :, k], tril_to_mat(S_k[:, k]), B_k[:, :, k], L_k[:, :, k])
        defect_S_k[k] = np.linalg.norm(S_k[:, k + 1] - triu(np.linalg.qr(X_kp1[:, :, k].T, 'r')), 'fro')

    return X_kp1, defect_S_k

def triu(X):
    '''
    Vectorize upper triangular (including diagonal) part of matrix
    '''
    triu_vec = X[0, :].reshape((X.shape[0], 1))
    for k in range(1, X.shape[0]):
        triu_vec = np.vstack([triu_vec, X[k, k:].reshape((-1, 1))])

    return triu_vec

def tril(X):
    '''
    Vectorize lower triangular (including diagonal) part of matrix
    '''
    tril_vec = X[:, 0].reshape((X.shape[0], 1))
    for k in range(1, X.shape[0]):
        tril_vec = np.vstack([tril_vec, X[k:, k].reshape((-1, 1))])

    return tril_vec

def triu_to_mat(X_vec):
    '''
    Matricize vector of upper triangular (including diagonal) part of matrix
    '''
    m = int((-1 + np.sqrt(1 + 8 * X_vec.size)) / 2) # inverse of triangular numbers
    X = np.zeros((m, m))
    i = 0
    for k in range(m):
        X[k, k:] = X_vec[i : i + (m - k)].flatten()
        i = i + (m - k)

    return X

def tril_to_mat(X_vec):
    '''
    Matricize vector of lower triangular (including diagonal) part of matrix
    '''
    m = int((-1 + np.sqrt(1 + 8 * X_vec.size)) / 2) # inverse of triangular numbers
    X = np.zeros((m, m))
    i = 0
    for k in range(m):
        a = X_vec[i : i + (m - k)]
        X[k:, k] = X_vec[i : i + (m - k)].flatten()
        i = i + (m - k)

    return X

def tril_to_mat_cp(X_vec):
    '''
    Matricize vector of lower triangular (including diagonal) part of matrix
    using cvxpy functions
    '''
    X = cp.vec_to_upper_tri(X_vec).T

    return X

def triu_cp(X):
    '''
    Vectorize upper triangular (including diagonal) part of matrix
    using cvxpy functions
    '''
    triu_vec = X[0, :].reshape((X.shape[0], 1), 'C')
    for k in range(1, X.shape[0]):
        triu_vec = cp.vstack([triu_vec, X[k, k:].reshape((-1, 1), 'C')])

    return triu_vec

def tril_cp(X):
    '''
    Vectorize lower triangular (including diagonal) part of matrix
    using cvxpy functions
    '''
    tril_vec = X[:, 0].reshape((X.shape[0], 1), 'C')
    for k in range(1, X.shape[0]):
        tril_vec = cp.vstack([tril_vec, X[k:, k].reshape((-1, 1), 'C')])

    return tril_vec