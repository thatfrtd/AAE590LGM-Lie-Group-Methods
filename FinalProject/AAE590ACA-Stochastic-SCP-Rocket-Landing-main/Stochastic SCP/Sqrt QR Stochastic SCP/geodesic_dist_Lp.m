function [dist] = geodesic_dist_Lp(L, K)
%GEODESIC_DIST_LP Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    L
    K
end

dist = sqrt(norm(tril(L) - tril(K), 'fro') ^ 2 + norm(log(diag(L)) - log(diag(K)), "fro") ^ 2);

end