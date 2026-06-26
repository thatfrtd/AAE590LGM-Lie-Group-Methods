import numpy as np
from Code.SSCP.qr_covariance import qr_pos_diag, X_Psqrt_noG, qr_derivative, qr_derivative_wrong, tril, triu, tril_cp, triu_cp, tril_to_mat, tril_to_mat_cp, triu_to_mat
from Code.SPD.LP_maps import geodesic_dist_LP

# Intiailize
x = np.random.standard_normal((10, 3))
dx = np.random.standard_normal((10, 3)) / 100
x2 = x + dx

# Get truth
Q, R = qr_pos_diag(x)
Q2, R2 = qr_pos_diag(x2)
dR = R2 - R
dQ = Q2 - Q

# Approximate
dR_approx1 = qr_derivative_wrong(dx, Q, R).value
LT_ck_1 = tril(dR_approx1 - np.diag(np.diag(dR_approx1)))

dR_approx2 = qr_derivative(dx, Q, R).value
LT_ck_2 = tril(dR_approx2 - np.diag(np.diag(dR_approx2)))

# Calculate Errors
err_0 = np.sum(np.abs(dR).flatten())
err_0_gd = geodesic_dist_LP(R2.T, R.T)
approx1_diff = dR - dR_approx1
err_1_1 = np.sum(np.abs(approx1_diff).flatten())
err_1_1_gd = geodesic_dist_LP(R2.T, R.T + dR_approx1.T)
approx2_diff = dR - dR_approx2
err_1_2 = np.sum(np.abs(approx2_diff).flatten())
err_1_2_gd = geodesic_dist_LP(R2.T, R.T + dR_approx2.T)

print(f"L+ Geodesic Error Zeroth Order QR R Deriv Approx: {err_0_gd}")
print(f"L+ Geodesic Error First Order QR R Deriv Approx (wrong): {err_1_1_gd}")
print(f"L+ Geodesic Error First Order QR R Deriv Approx: {err_1_2_gd}")
