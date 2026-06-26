% function [S_x] = trilmat(S_x_vec)
% %TRILMAT Summary of this function goes here
% %   Detailed explanation goes here
% arguments (Input)
%     S_x_vec
% end
% 
% nx = (-1 + sqrt(1 + 8 * numel(S_x_vec))) / 2; % inverse of triangular numbers
% 
% S_x = zeros(nx);
% S_x(find(tril(ones(nx)))) = S_x_vec;
% 
% end

function [S_x] = trilmat(S_x_vec)
%TRILMAT Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    S_x_vec
end

nx = (-1 + sqrt(1 + 8 * numel(S_x_vec))) / 2; % inverse of triangular numbers

S_x = [];
i = 1;
for k = 1 : nx
    S_x = [S_x, [zeros([k - 1, 1]); S_x_vec(i : (i + (nx - k)))]];
    i = i + (nx - k) + 1;
end

end