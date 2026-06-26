function [dR] = qr_r_deriv(dM, Q, R)
%QR_R_DERIV Summary of this function goes here
%   Derivative of R matrix from QR decomposition w.r.t. input M matrix
arguments (Input)
    dM
    Q
    R
end

A = Q.' * dM / R;
B = tril(A);
Q_TdQ = B - B.';
dR = Q.' * dM - Q_TdQ * R;

end