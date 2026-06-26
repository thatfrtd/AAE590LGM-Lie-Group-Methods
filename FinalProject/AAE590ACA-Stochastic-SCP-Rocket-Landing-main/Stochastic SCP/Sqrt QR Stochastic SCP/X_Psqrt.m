function [X] = X_Psqrt(A, S, B, L, G)
%X_PSQRT Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    A
    S
    B
    L
    G
end

X = [A * S + B * L, G];

end