function [x_sol, u_sol, S_x_sol, S_u_sol, sol_info] = solve_qr_stochastic_ptr_convex_subproblem_no_p_2(prob, ptr_ops, x_ref, u_ref, S_x_ref, S_u_ref)
%SOLVE_PTR_CONVEX_SUBPROBLEM Summary of this function goes here
%   Detailed explanation goes here
P_yk_sqrt = sqrtm_array(prob.disc.Ptilde_minus_k);

tri = @(k) k * (k + 1) / 2;

cvx_begin %quiet
    variable X(prob.n.x, prob.N)
    variable U(prob.n.u, prob.Nu)
    variable eta(1, prob.Nu)
    variable V(prob.n.x, prob.N - 1)
    variable V_S(tri(prob.n.x), prob.N - 1)
    variable v_prime(prob.n.ncvx)
    variable v_0(prob.n.x, 1)
    variable v_N(prob.n.x, 1)
    variable v_NP(1,1)
    variable S_x(tri(prob.n.x), prob.N)
    variable S_u(prob.n.u, prob.n.x, prob.N - 1)
    minimize( prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0, S_x, S_u) ...
        + virtual_control_cost(V, V_S, v_prime, v_0, v_N, v_NP, ptr_ops.w_vc) ...
        + trust_region_cost(eta, 0, ptr_ops.w_tr, 0))
    subject to
        % Dynamics
        for k = 1:(prob.N - 1)
            X(:, k + 1) == prob.disc.A_k(:, :, k) * prob.unscale_x(X(:, k)) ...
                         + prob.disc.B_k(:, :, k) * prob.unscale_u(U(:, k)) ...
                         + prob.disc.c_k(:, :, k) ...
                         + V(:, k);

            [S_x_kp1_approx, ~, delta_X_sum_square(k)] = S_first_prop(trilmat(S_x(:, k)), S_u(:, :, k), trilmat(S_x_ref(:, k)), S_u_ref(:, :, k), prob.disc.A_k(:, :, k), prob.disc.B_k(:, :, k), prob.disc.L_k(:, :, k + 1) * P_yk_sqrt(:, :, k));
            S_x(:, k + 1) == trilvec(S_x_kp1_approx) ... 
                           + V_S(:, k);
        end

        % Constraints
        for k = 1:prob.Nu
            % Convex Constraints
            for cc = 1:prob.n.cvx
                prob.convex_constraints{cc}(prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), 0, trilmat(S_x(:, k)), S_u(:, :, k)) <= 0;
            end
            % Nonconvex Constraints
            for nc = 1:prob.n.ncvx
                prob.nonconvex_constraints{nc}(prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), 0, trilmat(S_x(:, k)), S_u(:, :, k), x_ref, u_ref, 0, trilmat(S_x_ref(:, k)), S_u_ref, k) ...
                    - v_prime(nc) <= 0;
            end
        end
        v_prime >= 0;

        % Boundary Conditions
        prob.initial_bc(prob.unscale_x(X(:, 1)), 0) + v_0 == 0; % Initial mean state
        prob.S_0 == S_x(:, 1); % Initial estimated state covariance
        prob.terminal_bc(prob.unscale_x(X(:, prob.N)), 0) + v_N == 0; % Final mean state
        norm(trilmat(prob.S_f) \ trilmat(S_x(:, end))) - 1 - v_NP <= 0; % Final estimated state covariance

        v_NP >= 0;

        % Trust Region Constraints
        if ptr_ops.q == 22
            ptr_ops.alpha_x * sum_square(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu), 1) + ptr_ops.alpha_u * sum_square(U - u_ref, 1) + ptr_ops.alpha_u * sum_square(U - u_ref, 1) + ptr_ops.alpha_X * sum(delta_X_sum_square) <= eta;
        else
            ptr_ops.alpha_x * norms(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu), ptr_ops.q, 1) + ptr_ops.alpha_u * norms(U - u_ref, ptr_ops.q, 1) <= eta;
        end
cvx_end

x_sol = X;
u_sol = U;
S_x_sol = S_x;
S_u_sol = S_u;

if ptr_ops.q == 22
    eta_x = sum_square(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu), 1);
    eta_u = sum_square(U - u_ref, 1);
else
    eta_x = vecnorm(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu), ptr_ops.q, 1);
    eta_u = vecnorm(U - u_ref, ptr_ops.q, 1);
end

dx = vecnorm(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu), 2, 1);
du = vecnorm(U - u_ref, 2, 1);
dX = zeros([1, prob.N - 1]);
eta_X = zeros([1, prob.N - 1]);
Delta_S = zeros([1, prob.N + 1]);
for k = 1 : prob.N - 1
    [S_x_kp1_approx, dX(k), delta_X_sum_square_k] = S_first_prop(trilmat(S_x_sol(:, k)), S_u_sol(:, :, k), trilmat(S_x_ref(:, k)), S_u_ref(:, :, k), prob.disc.A_k(:, :, k), prob.disc.B_k(:, :, k), prob.disc.L_k(:, :, k + 1) * P_yk_sqrt(:, :, k));
    eta_X(k) = delta_X_sum_square_k;
    Delta_S(k + 1) = norm(S_x_sol(:, k + 1) - trilvec(S_x_kp1_approx));
end
Delta_S(1) = norm(prob.S_0 - S_x_sol(:, 1)); % Should always be 0...
Delta_S(end) = max(norm(trilmat(prob.S_f) \ trilmat(S_x_sol(:, end))) - 1, 0);

sol_info.status = cvx_status;
sol_info.vd = V;
sol_info.vd_S = V_S;
sol_info.vs = v_prime;
sol_info.vbc_0 = v_0;
sol_info.vbc_N = v_N;
sol_info.vbc_NP = norm(trilmat(prob.S_f) \ trilmat(S_x_sol(:, end))) - 1;
sol_info.J = prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0, S_x_sol, S_u_sol);
sol_info.J_tr = trust_region_cost(eta, 0, ptr_ops.w_tr, 0);
sol_info.J_vc = virtual_control_cost(V, V_S, v_prime, v_0, v_N, v_NP, ptr_ops.w_vc);
sol_info.dJ = 100 * (prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0, S_x_sol, S_u_sol) - prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), 0, S_x_ref, S_u_ref)) / prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), 0, S_x_ref, S_u_ref);
sol_info.dx = dx;
sol_info.du = du;
sol_info.dp = 0;
sol_info.dX = dX;
sol_info.delta = dx + du + dX;
sol_info.Delta_S = Delta_S;
sol_info.eta = eta;
sol_info.eta_x = eta_x;
sol_info.eta_u = eta_u;
sol_info.eta_p = 0;
sol_info.eta_X = eta_X;
end

function [J_tr] = trust_region_cost(eta, eta_p, w_tr, w_tr_p)
    J_tr = w_tr * eta' + w_tr_p * eta_p;
end

function [J_vc] = virtual_control_cost(V, V_S, v_prime, v_0, v_N, v_NP, w_vc)
    J_vc = w_vc * (norm(v_prime, 1) + v_NP + sum(norms(V, 1, 1)) + sum(norms(V_S, 1, 1)) + norm(v_0, 1) + norm(v_N, 1));
end