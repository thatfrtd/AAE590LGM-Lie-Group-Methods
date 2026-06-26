function [ptr_sol] = Stochastic_ptr(prob, ptr_ops, options)
arguments
    prob
    ptr_ops
    options.slack_control = false
end
%PTR Sequential Convex Programming algorithm
%   If converged, solution satisfies the nonlinear continuous-time equations of motion
% to within a tolerance on the order of eps_feasible feasible, satisfies all algebraic constraints at each
% temporal node, and approximates (local) optimality of the original optimal control problem.

x_ref = zeros([prob.n.x, prob.N, ptr_ops.iter_max + 1]);
u_ref = zeros([prob.n.u, prob.Nu, ptr_ops.iter_max + 1]);
p_ref = zeros([prob.n.p, ptr_ops.iter_max + 1]);
X_ref = zeros([prob.n.x, prob.n.x * prob.N * (prob.N + 1) / 2, ptr_ops.iter_max + 1]);
S_ref = zeros([prob.n.u, prob.n.x * prob.N * (prob.N + 1) / 2, ptr_ops.iter_max + 1]);

x_ref(:, :, 1) = prob.scale_x(prob.guess.x);
u_ref(:, :, 1) = prob.scale_u(prob.guess.u);
p_ref(:, 1) = prob.scale_p(prob.guess.p);

ptr_sol.converged = false;
ptr_sol.objective = zeros([1, ptr_ops.iter_max + 1]);
ptr_sol.Delta = zeros([prob.n.x, prob.N + 1, ptr_ops.iter_max + 1]);
ptr_sol.delta_xp = zeros([1, ptr_ops.iter_max]);

% Convexify along initial guess
[prob, ptr_sol.Delta(:, :, 1)] = convexify_along_reference(prob, prob.guess.x, prob.guess.u, prob.guess.p);

disp(" k |       status      |   vd  |   vs  | vbc_NP |  vbc_N |   J   |   J_tr  |   J_vc   |   dJ %  |   dx   |   du   |   dp   | delta |  dyn  |  eta  | eta_x | eta_u | eta_p")

for i = 1:(ptr_ops.iter_max)
    % Solve convex subproblem and update reference
    if prob.n.p == 0
        if options.slack_control
            [x_ref(:, :, i + 1), u_ref(:, :, i + 1), X_ref(:, :, i + 1), S_ref(:, :, i + 1), sol_info] = solve_stochastic_ptr_convex_subproblem_no_p_slackcontrol(prob, ptr_ops, x_ref(:, :, i), u_ref(:, :, i), X_ref(:, :, i), S_ref(:, :, i));
        else
            [x_ref(:, :, i + 1), u_ref(:, :, i + 1), X_ref(:, :, i + 1), S_ref(:, :, i + 1), sol_info] = solve_stochastic_ptr_convex_subproblem_no_p_2(prob, ptr_ops, x_ref(:, :, i), u_ref(:, :, i), X_ref(:, :, i), S_ref(:, :, i));
        end
    else
        [x_ref(:, :, i + 1), u_ref(:, :, i + 1), p_ref(:, i + 1), sol_info] = solve_stochastic_ptr_convex_subproblem(prob, ptr_ops, x_ref(:, :, i), u_ref(:, :, i), p_ref(:, i));
    end

    if i > 1
        K = recover_gain_matrix(X_ref(:, :, i + 1), S_ref(:, :, i + 1));
        K_ref = recover_gain_matrix(X_ref(:, :, i), S_ref(:, :, i));

        % Look at qr covariance propagation accuracy
        tri = @(k) k * (k + 1) / 2 * prob.n.x;
        P_yk_sqrt = sqrtm_array(prob.disc.Ptilde_minus_k);
        S_x_rel_err = zeros([prob.N - 1, 2]);
        diff_mag = zeros([prob.N - 1, 1]);
        for k = 1 : (prob.N - 1)
            k_ind = (tri(k - 1) + 1):tri(k);
            kp1_ind = (tri(k) + 1):tri(k + 1);
    
            S_x_k = qr_r_pos_diag(X_ref(:, k_ind, i + 1).').';
    
            S_x_ck = norm(S_x_k * S_x_k.' - X_ref(:, k_ind, i + 1) * X_ref(:, k_ind, i + 1).', "fro"); % Passes
            S_u_k = K(:, :, k) * S_x_k;

            S_x_k_ref = qr_r_pos_diag(X_ref(:, k_ind, i).').';
            S_x_kp1 = qr_r_pos_diag(X_Psqrt(prob.disc.A_k(:, :, k), S_x_k, prob.disc.B_k(:, :, k), S_u_k, prob.disc.L_k(:, :, k + 1) * P_yk_sqrt(:, :, k)).').';%qr_r_pos_diag(X_ref(:, kp1_ind, i + 1).').';
            S_x_kp1_ref = qr_r_pos_diag(X_ref(:, kp1_ind, i).').';
            
            S_x_ck2 = norm(S_x_kp1 * S_x_kp1.' - X_ref(:, kp1_ind, i + 1) * X_ref(:, kp1_ind, i + 1).', "fro"); % Passes
            
            S_u_k_ref = K_ref(:, :, k) * S_x_k_ref;
    
            S_u_ck = norm(S_u_k * S_u_k.' - S_ref(:, k_ind, i + 1) * S_ref(:, k_ind, i + 1).', "fro"); % Passes

            [S_x_kp1_approx, diff_mag(k)] = S_first_prop(S_x_k, S_u_k, S_x_k_ref, S_u_k_ref, prob.disc.A_k(:, :, k), prob.disc.B_k(:, :, k), prob.disc.L_k(:, :, k + 1) * P_yk_sqrt(:, :, k));
            S_x_rel_err(k, 1) = geodesic_dist_Lp(S_x_kp1, S_x_kp1_approx) / geodesic_dist_Lp(S_x_kp1, S_x_kp1_ref);
            S_x_rel_err(k, 2) = norm(S_x_kp1 - S_x_kp1_approx, "fro") / norm(S_x_kp1 - S_x_kp1_ref, "fro");
            a = (S_x_kp1 - S_x_kp1_approx) / S_x_kp1;
        end
    end

    % Convexify along reference trajectory
    [prob, ptr_sol.Delta(:, :, i + 1)] = convexify_along_reference(prob, prob.unscale_x(x_ref(:, :, i + 1)), prob.unscale_u(u_ref(:, :, i + 1)), prob.unscale_p(p_ref(:, i + 1)));

    % Update algorithm weights (4.24)
    %ptr_ops.w_tr = update_trust_region_weights(ptr_sol.Delta(:, i + 1)', ptr_ops.update_w_tr, ptr_ops.w_tr, ptr_ops.Delta_min);

    % Check stopping criteria (4.30)
    ptr_sol.delta_xp(i) = ptr_stopping(x_ref(:, :, i + 1), p_ref(:, i + 1), x_ref(:, :, i), p_ref(:, i), 2);

    % Display results of iteration
    if i ~= 1
        ptr_sol.info(i).dJ = ptr_sol.info(i - 1);
    end
    sol_info.delta = ptr_sol.delta_xp(i);
    if sum(norms(ptr_sol.Delta(:, :, i + 1), 2, 1)) <= ptr_ops.Delta_min
        sol_info.dyn = "True";
    else
        sol_info.dyn = sprintf("%.1g", norm(ptr_sol.Delta(:, i + 1)));
    end

    ptr_sol.info(i) = sol_info;

    if i == 1
        fprintf("%2.f | %17s | %5.1g | %5.1g | %5.1g | %5.1g | %5.3g | %5.3g | %5.3g |         | %6.1g | %6.1g | %5.1g | %5.1g | %5s | %5.1g | %5.1g | %5.1g | %5.1g\n", i, ptr_sol.info(i).status, norm(ptr_sol.info(i).vd), norm(ptr_sol.info(i).vs), ptr_sol.info(i).vbc_NP, norm(ptr_sol.info(i).vbc_N), ptr_sol.info(i).J, ptr_sol.info(i).J_tr, ptr_sol.info(i).J_vc, sum(ptr_sol.info(i).dx), sum(ptr_sol.info(i).du), sum(ptr_sol.info(i).dp), ptr_sol.info(i).delta, string(ptr_sol.info(i).dyn), norm(ptr_sol.info(i).eta), norm(ptr_sol.info(i).eta_x), norm(ptr_sol.info(i).eta_u), norm(ptr_sol.info(i).eta_p))
    else
        fprintf("%2.f | %17s | %5.1g | %5.1g | %5.1g | %5.1g | %5.3g | %5.3g | %5.3g | %6.3f | %6.1g | %6.1g | %5.1g | %5.1g | %5s | %5.1g | %5.1g | %5.1g | %5.1g\n", i, ptr_sol.info(i).status, norm(ptr_sol.info(i).vd), norm(ptr_sol.info(i).vs), ptr_sol.info(i).vbc_NP, norm(ptr_sol.info(i).vbc_N), ptr_sol.info(i).J, ptr_sol.info(i).J_tr, ptr_sol.info(i).J_vc, ptr_sol.info(i).dJ, sum(ptr_sol.info(i).dx), sum(ptr_sol.info(i).du), sum(ptr_sol.info(i).dp), ptr_sol.info(i).delta, string(ptr_sol.info(i).dyn), norm(ptr_sol.info(i).eta), norm(ptr_sol.info(i).eta_x), norm(ptr_sol.info(i).eta_u), norm(ptr_sol.info(i).eta_p))
    end

    if i >= ptr_ops.iter_min
        if ptr_sol.delta_xp(i) < ptr_ops.delta_tol && ~ptr_sol.converged && sol_info.dyn == "True"
            ptr_sol.converged = true;
            ptr_sol.converged_i = i + 1;
    
            break
        end
    end
end

if ptr_sol.converged == false
    warning("PTR did not converge after %g iterations. delta_xp = %.3f. norm(Delta) = %.3f\n", i, ptr_sol.delta_xp(i), norm(ptr_sol.Delta(:, end)))
end

ptr_sol.x = prob.unscale_x(x_ref);
ptr_sol.u = prob.unscale_u(u_ref);
ptr_sol.p = prob.unscale_p(p_ref);
ptr_sol.X = X_ref;
ptr_sol.S = S_ref;
end

function [w_tr] = update_trust_region_weights(Delta, update, w_tr, Delta_min)
    w_tr_update = 1 ./ abs(Delta) .* (abs(Delta) >= Delta_min) ...
             + 1 / Delta_min * (abs(Delta) < Delta_min);

    w_tr = update .* w_tr_update + ~update .* w_tr;
end

function [delta_xp] = ptr_stopping(X, p, x_ref, p_ref, q)
    delta_xp = zero_if_empty(vecnorm(p - p_ref, q, 1)) + max(vecnorm(X - x_ref, q, 1));
end