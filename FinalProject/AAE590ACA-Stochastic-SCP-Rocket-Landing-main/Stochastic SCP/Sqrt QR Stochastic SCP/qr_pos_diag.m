function [Q, R] = qr_pos_diag(M)
%QR_POS_DIAG Summary of this function goes here
%   QR decomposition of matrix with enforcement of positive diagonal condition
arguments (Input)
    M
end

[Q, R] = qr(M, "econ");

% Fix negative diagonal elemnts of R while maintaining decomposition
for k = 1 : min(size(M))
    if R(k, k) < 0
        R(k, :) = -R(k, :);
        Q(:, k) = -Q(:, k);
    end
end

end