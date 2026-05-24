function [u_cbf, qp_infeasible, h_val, slack_used] = safety_filter_cbf(z, u_nom, params, cbf_opts)
%SAFETY_FILTER_CBF  Higher-Order CBF (HOCBF) safety filter with slack relaxation
%
%   [u_cbf, qp_infeasible, h_val, slack_used] = safety_filter_cbf(z, u_nom, params, cbf_opts)
%
%   The actual safety constraint
%       h_0(z) = q_v_max^2 - ||q_v||^2
%   has relative degree 2 with respect to u, so we use the HOCBF
%   construction (Xiao & Belta 2019, Ames et al. 2019):
%       h_1(z) = h_0_dot(z) + kappa_0 * h_0(z)
%       d/dt(h_1) + kappa_1 * h_1 >= 0    <- enforced by the QP
%
%   After algebraic simplification:
%       h_0_dot       = -q_0 * (q_v' * omega)
%       d/dt(h_0_dot) = 0.5*(q_v'omega)^2 - 0.5*q_0^2*||omega||^2
%                       - q_0 * q_v' * J^{-1} * u
%       d/dt(h_1)     = d/dt(h_0_dot) + kappa_0 * h_0_dot
%
%   SLACK-VARIABLE RELAXATION
%   -------------------------
%   Under input saturation, the hard HOCBF condition can be infeasible.
%   We use the standard slack-relaxed QP (penalty p on slack):
%       min_{u, xi}  0.5*||u - u_nom||^2 + 0.5 * p * xi^2
%       s.t.   L_g h_1 * u + xi >= -drift_h_1 - kappa_1 * h_1
%              -u_max <= u <= u_max
%              xi >= 0
%   This always returns a feasible u, while reporting slack_used > 0 as a
%   diagnostic for when the CBF could not enforce safety.
%
%   The disturbance is OMITTED by design (standard non-robust baseline).
%   Gyroscopic compensation u_applied = u + omega x J*omega is added by
%   the caller, matching the LMI branch.
%
%   Inputs
%   ------
%   z         6x1  reduced state [q_v; omega]
%   u_nom     3x1  nominal control torque (already saturated)
%   params    struct from loadSystemParameters()
%   cbf_opts  struct with fields:
%       .q_v_max         scalar, e.g. sin(theta_max/2)
%       .u_max           3x1, per-axis torque limit
%       .kappa_0         class-K gain for inner level (default 20)
%       .kappa_1         class-K gain for outer level (default 5)
%       .slack_penalty   penalty p on slack^2 (default 1e6)
%
%   Outputs
%   -------
%   u_cbf          3x1, filtered control torque (BEFORE gyroscopic comp)
%   qp_infeasible  logical, true if QP solver failed entirely
%   h_val          scalar, h_0(z) at current state (actual safety margin)
%   slack_used     scalar, slack variable value (>0 means CBF condition
%                  could not be enforced under input limits)

    % ---- Unpack ----
    qv    = z(1:3);
    omega = z(4:6);

    q_v_max   = cbf_opts.q_v_max;
    u_max_vec = cbf_opts.u_max;
    kappa_0   = cbf_opts.kappa_0;
    kappa_1   = cbf_opts.kappa_1;
    if isfield(cbf_opts, 'slack_penalty')
        p_slack = cbf_opts.slack_penalty;
    else
        p_slack = 1e6;
    end

    J     = diag([params.Ixx, params.Iyy, params.Izz]);
    J_inv = J \ eye(3);

    % ---- Useful scalars ----
    qv_sq   = qv' * qv;
    q0      = sqrt(max(1 - qv_sq, 1e-9));      % guarded sqrt
    qv_om   = qv' * omega;
    om_sq   = omega' * omega;

    % ---- Safety barrier and its derivatives ----
    h_0      = q_v_max^2 - qv_sq;
    h_0_dot  = -q0 * qv_om;
    h_1      = h_0_dot + kappa_0 * h_0;

    drift_h0_dot = 0.5 * qv_om^2 - 0.5 * q0^2 * om_sq;
    drift_h1     = drift_h0_dot + kappa_0 * h_0_dot;
    L_g_h1       = -q0 * (qv' * J_inv);    % 1x3

    % ---- Outside-safe-set guard ----
    % Once h_0 < 0 the system is already in the unsafe set and the HOCBF
    % construction no longer provides any guarantee.  Continuing the QP
    % can produce wild control as the formulation is mathematically out
    % of its domain (h_1 dynamics become meaningless).  Switch to a
    % bounded emergency-braking controller that simply decelerates omega.
    % This makes the CBF's failure mode "give up and brake" rather than
    % "thrash unboundedly", which is more representative of how CBFs are
    % deployed in practice and cleaner for visualization.
    if h_0 < 0
        u_brake = -1e-3 * (J * omega);    % proportional to angular momentum
        u_cbf = max(min(u_brake, u_max_vec), -u_max_vec);
        slack_used = abs(h_0);            % report severity of violation
        qp_infeasible = true;
        h_val = h_0;
        return;
    end

    % ---- Slack-relaxed QP setup ----
    % Decision variable: [u; xi], where xi >= 0 is the slack
    % Objective: 0.5*u'u - u_nom'u + 0.5*p*xi^2
    H = blkdiag(eye(3), p_slack);
    f = [-u_nom; 0];
    % Inequality: -L_g_h1*u - xi <= drift_h1 + kappa_1*h_1
    A_ineq = [-L_g_h1, -1];
    b_ineq = drift_h1 + kappa_1 * h_1;
    % Box: -u_max <= u <= u_max, 0 <= xi <= inf
    lb = [-u_max_vec; 0];
    ub = [ u_max_vec; inf];

    persistent qp_opts
    if isempty(qp_opts)
        qp_opts = optimoptions('quadprog', 'Display', 'off', ...
                               'Algorithm', 'interior-point-convex');
    end

    [sol, ~, exitflag] = quadprog(H, f, A_ineq, b_ineq, [], [], lb, ub, [], qp_opts);

    if exitflag > 0 && ~isempty(sol)
        u_cbf      = sol(1:3);
        slack_used = sol(4);
        % Threshold of 1e-3 ignores floating-point noise from quadprog;
        % only flag steps where the CBF condition was non-trivially
        % relaxed (slack >= 1e-3 corresponds to a meaningful violation
        % of the ḣ_1 + κ_1 h_1 >= 0 condition).
        qp_infeasible = (slack_used > 1e-3);
    else
        % True solver failure: emergency braking (decelerate omega)
        u_brake = -1e-3 * (J * omega);     % proportional braking
        u_cbf = max(min(u_brake, u_max_vec), -u_max_vec);
        slack_used = inf;
        qp_infeasible = true;
    end

    h_val = h_0;
end
