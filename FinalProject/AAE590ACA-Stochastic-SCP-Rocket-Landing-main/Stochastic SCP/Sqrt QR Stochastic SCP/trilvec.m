function [S_x_vec] = trilvec(S_x)
%TRILVEC Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    S_x
end

nx = size(S_x, 1);

S_x_vec = S_x(find(tril(ones(nx))));

end