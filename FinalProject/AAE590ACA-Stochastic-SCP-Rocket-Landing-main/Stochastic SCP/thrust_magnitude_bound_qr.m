function [Gamma_k] = thrust_magnitude_bound_qr(S_u_ref, u_ref, k, t_k, T_max, m_0, alpha, nu, nx)
%GAMMA_K Summary of this function goes here
%   Detailed explanation goes here
tri = @(k) k * (k + 1) / 2 * nx;

delta_t = t_k(2) - t_k(1);

TWR_0 = T_max / m_0;

TWR_multiplier_k = exp(einsum(@(i) alpha * delta_t ...
    * max(0, norm(u_ref(1:2, i)) ...
           - sigma_mag_confidence(1e-3 / (2 * k), nu) * norm(S_u_ref(:, :, k))), 1:(k - 1)));

Gamma_k = TWR_0 * TWR_multiplier_k;
end

