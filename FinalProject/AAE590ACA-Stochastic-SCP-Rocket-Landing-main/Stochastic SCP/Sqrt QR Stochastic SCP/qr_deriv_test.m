
% Initialize
x = randn([10, 3]);
dx = randn([10, 3]) / 10;
x2 = x + dx;

% Get truth
[Q, R] = qr_pos_diag(x);
[Q2, R2] = qr_pos_diag(x2);
dR = R2 - R;
dQ = Q2 - Q;
dx_gd = geodesic_dist_Lp(R, R2);

% Approximate
dR_approx1 = qr_r_deriv(dx, Q, R);
LT_ck_1 = tril(dR_approx1 - diag(diag(dR_approx1)));

% Calculate errors
err_0 = sum(abs(dR));
err_0_gd = geodesic_dist_Lp(R2.', R.');
approx1_diff = dR - dR_approx1;

err_1 = sum(abs(approx1_diff));
err_1_gd = geodesic_dist_Lp(R2.', R.' + dR_approx1.');

fprintf("L+ dx geo size %.5f\n", dx_gd)
fprintf("L+ Geodesic Error Zeroth Order QR R Deriv Approx: %.5f\n", err_0_gd)
fprintf("L+ Geodesic Error First Order QR R Deriv Approx %.5f\n", err_1_gd)
