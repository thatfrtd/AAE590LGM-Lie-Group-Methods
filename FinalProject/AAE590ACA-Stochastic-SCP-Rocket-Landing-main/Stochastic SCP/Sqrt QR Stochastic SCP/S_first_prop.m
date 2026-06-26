function [S_kp1, diff_mag, delta_X_sum_square] = S_first_prop(S, L, S_ref, L_ref, A, B, G)
%S_FIRST_PROP Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    S
    L
    S_ref
    L_ref
    A
    B
    G
end

X = X_Psqrt(A, S, B, L, G);
X_ref = X_Psqrt(A, S_ref, B, L_ref, G);
X_diff = X - X_ref;

[Q_ref, R_ref] = qr_pos_diag(X_ref.');

S_kp1_ref = R_ref.';
S_kp1 = S_kp1_ref + qr_r_deriv(X_diff.', Q_ref, R_ref).';

diff_mag = norm(X_diff, "fro");
delta_X_sum_square = sum_square(X_diff(:));

end