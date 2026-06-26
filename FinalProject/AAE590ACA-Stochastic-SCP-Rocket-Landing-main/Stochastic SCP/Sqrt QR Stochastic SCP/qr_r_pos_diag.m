function [R] = qr_r_pos_diag(M)
%QR_R_POS_DIAG Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    M
end

[~, R] = qr_pos_diag(M);

end