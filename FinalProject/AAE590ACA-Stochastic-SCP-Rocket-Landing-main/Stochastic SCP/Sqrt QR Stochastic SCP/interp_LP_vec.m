function [L_interp] = interp_LP_vec(L0_vec, LN_vec, N)
%INTERP_LP_VEC Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    L0_vec
    LN_vec
    N
end

L0 = trilmat(L0_vec);
LN = trilmat(LN_vec);
err = err_LP(L0, LN);
Leye = eye(size(L0, 1));

t = linspace(0, 1, N);
L_interp = zeros(size(L0_vec, 1), N);
for k = 1 : N
    L_interp(:, k) = trilvec(compose_LP(L0, Exp_LP(Leye, err * t(k))));
end

end

function [Exp_X] = Exp_LP(L, X)
Exp_X = tril(L, -1) + tril(X, -1) + diag(diag(L) .* exp(diag(X) ./ diag(L)));
end

function [Log_K] = Log_LP(L, K)
Log_K = tril(K, -1) - tril(L, -1) + diag(diag(L) .* log(diag(K) ./ diag(L)));
end

function [L12] = compose_LP(L1, L2)
L12 = tril(L1, -1) + tril(L2, -1) + diag(diag(L1) .* diag(L2));
end

function [L_inv] = inverse_LP(L)
L_inv = diag(1 ./ diag(L)) - tril(L, -1);
end

function [err] = err_LP(P, Q)
err = Log_LP(eye(size(P, 1)), compose_LP(inverse_LP(P), Q));
end