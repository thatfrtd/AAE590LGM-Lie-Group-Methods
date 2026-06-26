function [ptr_sol] = Stochastic_qr_ptr(prob, ptr_ops)
arguments
    prob
    ptr_ops
end
%PTR Sequential Convex Programming algorithm
%   If converged, solution satisfies the nonlinear continuous-time equations of motion
% to within a tolerance on the order of eps_feasible feasible, satisfies all algebraic constraints at each
% temporal node, and approximates (local) optimality of the original optimal control problem.

prob.S_0 = trilvec(qr_r_pos_diag(chol(prob.Phat0))');
prob.S_f = trilvec(qr_r_pos_diag(chol(prob.Pf))');

x_ref = zeros([prob.n.x, prob.N, ptr_ops.iter_max + 1]);
u_ref = zeros([prob.n.u, prob.Nu, ptr_ops.iter_max + 1]);
p_ref = zeros([prob.n.p, ptr_ops.iter_max + 1]);
S_x_ref = zeros([prob.n.x * (prob.n.x + 1) / 2, prob.N, ptr_ops.iter_max + 1]);
S_u_ref = zeros([prob.n.u, prob.n.x, prob.N - 1, ptr_ops.iter_max + 1]);

x_ref(:, :, 1) = prob.scale_x(prob.guess.x);
u_ref(:, :, 1) = prob.scale_u(prob.guess.u);
p_ref(:, 1) = prob.scale_p(prob.guess.p);
S_x_ref(:, :, 1) = interp_LP_vec(prob.S_0, trilvec(qr_r_pos_diag(chol(prob.Pf))'), prob.N);
% S_u_ref guess is zeros (no feedback) - could solve deterministic cov
% steering first which doesn't need trust region...

ptr_sol.converged = false;
ptr_sol.objective = zeros([1, ptr_ops.iter_max + 1]);
ptr_sol.Delta = zeros([prob.n.x, prob.N + 1, ptr_ops.iter_max + 1]);
ptr_sol.delta = zeros([1, ptr_ops.iter_max]);

% Convexify along initial guess
[prob, ptr_sol.Delta(:, :, 1)] = convexify_along_reference(prob, prob.guess.x, prob.guess.u, prob.guess.p);

disp(" k |       status      |   vd  |   vs  | vbc_NP |  vbc_N |   J   |  J_tr  |   J_vc   |   dJ %  |   dx   |   du   |   dX   | delta |  dyn  |  dyn_S  |  eta  | eta_x | eta_u | eta_X")

for i = 1:(ptr_ops.iter_max)
    % Solve convex subproblem and update reference
    [x_ref(:, :, i + 1), u_ref(:, :, i + 1), S_x_ref(:, :, i + 1), S_u_ref(:, :, :, i + 1), sol_info] = solve_qr_stochastic_ptr_convex_subproblem_no_p_2(prob, ptr_ops, x_ref(:, :, i), u_ref(:, :, i), S_x_ref(:, :, i), S_u_ref(:, :, :, i));

    % Convexify along reference trajectory
    [prob, ptr_sol.Delta(:, :, i + 1)] = convexify_along_reference(prob, prob.unscale_x(x_ref(:, :, i + 1)), prob.unscale_u(u_ref(:, :, i + 1)), prob.unscale_p(p_ref(:, i + 1)));

    % Display results of iteration
    total_defect = sum(norms(ptr_sol.Delta(:, :, i + 1), 2, 1));
    if total_defect <= ptr_ops.Delta_min
        sol_info.dyn = "True";
    else
        sol_info.dyn = sprintf("%.1g", total_defect);
    end
    total_defect_S = sum(sol_info.Delta_S);
    if total_defect_S <= ptr_ops.Delta_S_min
        sol_info.dyn_S = "True";
    else
        sol_info.dyn_S = sprintf("%.1g", total_defect_S);
    end

    ptr_sol.info(i) = sol_info;
    ptr_sol.delta(i) = sum(sol_info.delta);

    if i == 1
        fprintf("%2.f | %17s | %5.1g | %5.1g | %5.1g  | %5.1g  | %5.3g | %5.3g | %5.3g |         | %6.1g | %6.1g | %5.1g | %5.1g | %5s | %5s | %5.1g | %5.1g | %5.1g | %5.1g\n", i, ptr_sol.info(i).status, norm(ptr_sol.info(i).vd), norm(ptr_sol.info(i).vs), ptr_sol.info(i).vbc_NP, norm(ptr_sol.info(i).vbc_N), ptr_sol.info(i).J, ptr_sol.info(i).J_tr, ptr_sol.info(i).J_vc, sum(ptr_sol.info(i).dx), sum(ptr_sol.info(i).du), sum(ptr_sol.info(i).dX), sum(ptr_sol.info(i).delta), string(ptr_sol.info(i).dyn), string(ptr_sol.info(i).dyn_S), norm(ptr_sol.info(i).eta, 1), norm(ptr_sol.info(i).eta_x, 1), norm(ptr_sol.info(i).eta_u, 1), norm(ptr_sol.info(i).eta_X, 1))
    else
        fprintf("%2.f | %17s | %5.1g | %5.1g | %5.1g  | %5.1g  | %5.3g | %5.3g | %5.3g | %6.3f | %6.1g | %6.1g | %5.1g | %5.1g | %5s | %5s | %5.1g | %5.1g | %5.1g | %5.1g\n", i, ptr_sol.info(i).status, norm(ptr_sol.info(i).vd), norm(ptr_sol.info(i).vs), ptr_sol.info(i).vbc_NP, norm(ptr_sol.info(i).vbc_N), ptr_sol.info(i).J, ptr_sol.info(i).J_tr, ptr_sol.info(i).J_vc, ptr_sol.info(i).dJ, sum(ptr_sol.info(i).dx), sum(ptr_sol.info(i).du), sum(ptr_sol.info(i).dX), sum(ptr_sol.info(i).delta), string(ptr_sol.info(i).dyn), string(ptr_sol.info(i).dyn_S), norm(ptr_sol.info(i).eta, 1), norm(ptr_sol.info(i).eta_x, 1), norm(ptr_sol.info(i).eta_u, 1), norm(ptr_sol.info(i).eta_X, 1))
    end

    if i >= ptr_ops.iter_min
        if ptr_sol.delta(i) < ptr_ops.delta_tol && ~ptr_sol.converged && sol_info.dyn == "True"
            ptr_sol.converged = true;
            ptr_sol.converged_i = i + 1;
    
            break
        end
    end
end

if ptr_sol.converged == false
    warning("PTR did not converge after %g iterations. delta_xp = %.3f. norm(Delta) = %.3f\n", i, ptr_sol.delta(i), norm(ptr_sol.Delta(:, end)))
end

ptr_sol.x = prob.unscale_x(x_ref);
ptr_sol.u = prob.unscale_u(u_ref);
ptr_sol.p = prob.unscale_p(p_ref);
ptr_sol.S_x = S_x_ref;
ptr_sol.S_u = S_u_ref;
end
