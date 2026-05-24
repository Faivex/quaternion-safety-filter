%% conservatism_analysis.m
%
% Quantifies the conservatism of the proposed sector-bound LMI synthesis
% (Theorem 3) against a linearization-based LMI baseline (Remark 8) over
% a 2D grid in design-space (theta_max, omega_max).
%
% For each (theta_max, omega_max) pair the script:
%   - Solves the full sector-bound LMI (Theorem 3)
%   - Solves the linearization-based LMI with two-multiplier S-procedure
%     (faithful Remark 8: phi treated as ||phi|| <= gamma*omega_max,
%      independent of the disturbance d_bar)
%   - Records log det Q for both (or marks infeasibility)
%
% Outputs:
%   - Console block with nominal-point comparison
%   - CSV table with the full 2D sweep
%   - PDF figure: feasibility map + log det Q surface for the proposed
%
% Reframing (Option A from the discussion):
%   The point is no longer "X times bigger ellipsoid" but
%     "linearized formulation is infeasible across most of design-space;
%      proposed is feasible everywhere."

clc; clear; close all;

%% Common setup
params = loadSystemParameters();
J = diag([params.Ixx, params.Iyy, params.Izz]);
J_inv = J \ eye(3);

A = [zeros(3,3), 0.5*eye(3); zeros(3,3), zeros(3,3)];
B = [zeros(3,3); J_inv];
B_phi = [eye(3); zeros(3,3)];
E_sec = [zeros(3,3), eye(3)];

n = 6; m = 3;

d_max     = 1e-5;
E_bar     = [zeros(3,3); d_max * J_inv];
u_max     = [1e-4; 1e-4; 1e-4];

alpha_vals = logspace(-1, 2, 5);
tau_vals   = logspace(-1, 2, 5);
mu_vals    = logspace(-1, 2, 5);

%% Nominal-point comparison
fprintf('\n========================================================\n');
fprintf('  Nominal point: theta_max = 40 deg, omega_max = 5 rad/s\n');
fprintf('========================================================\n');

[r_full_nom, r_lin_nom] = compare_at(...
    40, 5.0, u_max, A, B, B_phi, E_sec, E_bar, ...
    alpha_vals, tau_vals, mu_vals, n, m);

print_nominal(40, 5.0, r_full_nom, r_lin_nom);

%% 2D sweep
fprintf('\n========================================================\n');
fprintf('  2D sweep: theta_max in [10, 60] deg, omega_max in [1, 10] rad/s\n');
fprintf('========================================================\n');

theta_sweep = 5:5:60;
omega_sweep = 1:1:10;           

N_th = length(theta_sweep);
N_om = length(omega_sweep);

logdet_full = nan(N_om, N_th);   % rows = omega, cols = theta
logdet_lin  = nan(N_om, N_th);
feas_full   = false(N_om, N_th);
feas_lin    = false(N_om, N_th);

t_start = tic;
for i = 1:N_om
    for j = 1:N_th
        th = theta_sweep(j); om = omega_sweep(i);
        fprintf('  (theta=%2d deg, omega=%4.1f rad/s) ... ', th, om);
        try
            [rF, rL] = compare_at(th, om, u_max, A, B, B_phi, E_sec, E_bar, ...
                alpha_vals, tau_vals, mu_vals, n, m);
            if ~isempty(rF.Q)
                logdet_full(i,j) = rF.logdetQ; feas_full(i,j) = true;
            end
            if ~isempty(rL.Q)
                logdet_lin(i,j)  = rL.logdetQ; feas_lin(i,j)  = true;
            end
            fprintf('proposed: %s   linearized: %s\n', ...
                fmt(logdet_full(i,j)), fmt(logdet_lin(i,j)));
        catch ME
            fprintf('FAILED (%s)\n', ME.message);
        end
    end
end
fprintf('Total sweep time: %.1f minutes\n', toc(t_start)/60);

%% Save raw sweep data for standalone re-plotting
results_dir = 'results/conservatism';
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
sweep_data = struct( ...
    'theta_sweep', theta_sweep, ...
    'omega_sweep', omega_sweep, ...
    'logdet_full', logdet_full, ...
    'logdet_lin',  logdet_lin,  ...
    'feas_full',   feas_full,   ...
    'feas_lin',    feas_lin     );
save(fullfile(results_dir, 'conservatism_sweep_data.mat'), '-struct', 'sweep_data');
fprintf('Sweep data saved to %s\n', fullfile(results_dir, 'conservatism_sweep_data.mat'));

%% Save results
results_dir = 'results/conservatism';
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

%% Summary statistics
n_grid = N_th * N_om;
fprintf('\n--- Summary ---\n');
fprintf('Grid points: %d\n', n_grid);
fprintf('Proposed   feasible at: %d / %d (%.0f%%)\n', ...
    sum(feas_full(:)), n_grid, 100*sum(feas_full(:))/n_grid);
fprintf('Linearized feasible at: %d / %d (%.0f%%)\n', ...
    sum(feas_lin(:)),  n_grid, 100*sum(feas_lin(:))/n_grid);

%% Plot via standalone function (also runnable independently)
plot_conservatism_sweep(fullfile(results_dir, 'conservatism_sweep_data.mat'), ...
                       results_dir);


%% ========================================================================
%%                              FUNCTIONS
%% ========================================================================

function [res_full, res_lin] = compare_at(theta_max_deg, omega_max, u_max, ...
        A, B, B_phi, E_sec, E_bar, alpha_vals, tau_vals, mu_vals, n, m)
    q_v_max = sin(theta_max_deg * pi / 180 / 2);
    gamma   = sqrt(0.5 * (1 - sqrt(1 - q_v_max^2)));

    res_full = solve_full_lmi(A, B, B_phi, E_sec, E_bar, gamma, ...
        q_v_max, omega_max, u_max, alpha_vals, tau_vals, n, m);
    res_full.q_v_max = q_v_max;
    res_full.theta_max_deg = theta_max_deg;

    Delta_max = gamma * omega_max;
    res_lin = solve_linear_lmi(A, B, B_phi, E_bar, Delta_max, ...
        q_v_max, omega_max, u_max, alpha_vals, mu_vals, n, m);
    res_lin.q_v_max = q_v_max;
end

function res = solve_full_lmi(A, B, B_phi, E_sec, E_bar, gamma, ...
        q_v_max, omega_max, u_max, alpha_vals, tau_vals, n, m)
    best_vol = -inf; best_Q = []; best_Y = [];
    best_alpha = nan; best_tau = nan;
    for alpha = alpha_vals
        for tau = tau_vals
            cvx_begin sdp quiet
                cvx_solver mosek
                variable Q(n,n) symmetric
                variable Y(m,n)
                Q >= 1e-4 * eye(n);
                Psi_11 = A*Q + Q*A' + B*Y + Y'*B' + alpha*Q;
                M = [Psi_11,                  B_phi,         E_bar,        gamma*sqrt(tau)*Q*E_sec';
                     B_phi',                  -tau*eye(3),   zeros(3),     zeros(3);
                     E_bar',                  zeros(3),      -alpha*eye(3),zeros(3);
                     gamma*sqrt(tau)*E_sec*Q, zeros(3),      zeros(3),     -eye(3)];
                M <= 0;
                for i = 1:m
                    [Q, Y(i,:)'; Y(i,:), u_max(i)^2] >= 0;
                end
                C1 = [eye(3), zeros(3)];
                [Q, Q*C1'; C1*Q, q_v_max^2*eye(3)] >= 0;
                C2 = [zeros(3), eye(3)];
                [Q, Q*C2'; C2*Q, omega_max^2*eye(3)] >= 0;
                maximize(log_det(Q))
            cvx_end
            if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
                v = log(det(Q));
                if v > best_vol
                    best_vol = v; best_Q = Q; best_Y = Y;
                    best_alpha = alpha; best_tau = tau;
                end
            end
        end
    end
    res.Q = best_Q; res.Y = best_Y;
    res.logdetQ = best_vol;
    res.alpha = best_alpha; res.tau = best_tau;
end

function res = solve_linear_lmi(A, B, B_phi, E_bar, Delta_max, ...
        q_v_max, omega_max, u_max, alpha_vals, mu_vals, n, m)
    % Two-multiplier Remark 8: alpha for d_bar, mu for phi.
    best_vol = -inf; best_Q = []; best_Y = [];
    best_alpha = nan; best_mu = nan;
    for alpha = alpha_vals
        for mu = mu_vals
            cvx_begin sdp quiet
                cvx_solver mosek
                variable Q(n,n) symmetric
                variable Y(m,n)
                Q >= 1e-4 * eye(n);
                Psi_11 = A*Q + Q*A' + B*Y + Y'*B' + (alpha + mu) * Q;
                M = [Psi_11,   B_phi,                         E_bar;
                     B_phi',   -(mu / Delta_max^2) * eye(3),  zeros(3);
                     E_bar',   zeros(3),                      -alpha * eye(3)];
                M <= 0;
                for i = 1:m
                    [Q, Y(i,:)'; Y(i,:), u_max(i)^2] >= 0;
                end
                C1 = [eye(3), zeros(3)];
                [Q, Q*C1'; C1*Q, q_v_max^2*eye(3)] >= 0;
                C2 = [zeros(3), eye(3)];
                [Q, Q*C2'; C2*Q, omega_max^2*eye(3)] >= 0;
                maximize(log_det(Q))
            cvx_end
            if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
                v = log(det(Q));
                if v > best_vol
                    best_vol = v; best_Q = Q; best_Y = Y;
                    best_alpha = alpha; best_mu = mu;
                end
            end
        end
    end
    res.Q = best_Q; res.Y = best_Y;
    res.logdetQ = best_vol;
    res.alpha = best_alpha; res.tau = best_mu;
end

function print_nominal(th, om, rF, rL)
    fprintf('\n  --------------------------------------------------------\n');
    if ~isempty(rF.Q)
        fprintf('  Proposed (Thm 3)   log det Q = %+8.3f   alpha* = %.3f   tau* = %.3f\n', ...
            rF.logdetQ, rF.alpha, rF.tau);
    else
        fprintf('  Proposed (Thm 3)   INFEASIBLE\n');
    end
    if ~isempty(rL.Q)
        fprintf('  Linearized (Rem 8) log det Q = %+8.3f   alpha* = %.3f   mu* = %.3f\n', ...
            rL.logdetQ, rL.alpha, rL.tau);
    else
        fprintf('  Linearized (Rem 8) INFEASIBLE\n');
    end
end

function s = fmt(x)
    if isnan(x); s = '  infeas.'; else; s = sprintf('%+8.3f', x); end
end
