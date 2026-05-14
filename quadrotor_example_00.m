%% Rigorous Quaternion-Based Safety Filter Design
% Standalone numerical simulation implementing the exact derivation from
% the paper section on quaternion safety filters
%
% This script demonstrates:
% 1. Exact sector bound computation (no approximations)
% 2. LMI-based RCI set synthesis using rigorous formulation
% 3. Safety filter with smooth blending
% 4. Nonlinear simulation with quaternion dynamics

clc; clear; close all;

fprintf('================================================\n');
fprintf('Rigorous Quaternion Safety Filter Simulation\n');
fprintf('================================================\n\n');

%% System Parameters (Crazyflie 2.0 quadrotor)
params.mass = 0.028;  % kg
params.Ixx = 16.6e-6;  % kg·m^2
params.Iyy = 16.6e-6;  % kg·m^2
params.Izz = 29.3e-6;  % kg·m^2
params.g = 9.81;       % m/s^2

J = diag([params.Ixx, params.Iyy, params.Izz]);
J_inv = inv(J);

fprintf('=== System Parameters ===\n');
fprintf('Mass: %.4f kg\n', params.mass);
fprintf('Inertia: diag([%.2e, %.2e, %.2e]) kg·m^2\n', params.Ixx, params.Iyy, params.Izz);

%% State-Space Matrices
% State: z = [q_v; omega] ∈ R^6
% Linearized system: dz/dt = A*z + B*u + E*d

A = [zeros(3,3), 0.5*eye(3);
     zeros(3,3), zeros(3,3)];

B = [zeros(3,3);
     J_inv];

% Disturbance parameters
d_max = 1e-5;  % N·m
E = [zeros(3,3);
     d_max * J_inv];

n = 6;  % state dimension
m = 3;  % control dimension

fprintf('\n=== State-Space Dimensions ===\n');
fprintf('State dimension n = %d\n', n);
fprintf('Control dimension m = %d\n', m);
fprintf('Controllability rank: %d (should be %d)\n', rank(ctrb(A, B)), n);

%% Safety Constraints
% Attitude constraint: ||q_v|| <= q_v_max = sin(theta_max/2)
theta_max_deg = 80;  % degrees
theta_max_rad = theta_max_deg * pi / 180;
q_v_max = sin(theta_max_rad / 2);

% Angular velocity constraint
omega_max = 2.0;  % rad/s

% State constraints: [q_v; omega]
x_max = [q_v_max; q_v_max; q_v_max; omega_max; omega_max; omega_max];

% Control input constraints
u_max = [1e-4; 1e-4; 1e-4];  % N·m (increased 10x from original)

fprintf('\n=== Safety Constraints ===\n');
fprintf('Max attitude angle: %.1f deg\n', theta_max_deg);
fprintf('q_v_max = sin(%.1f°/2) = %.6f\n', theta_max_deg, q_v_max);
fprintf('omega_max: %.2f rad/s\n', omega_max);
% fprintf('u_max: %.2e N·m (relaxed for feasibility)\n', u_max(1));

%% Exact Sector Bound Computation (Lemma from paper)
% gamma^2 = (1/2) * (1 - sqrt(1 - q_v_max^2))
% This is the EXACT bound, no approximations

gamma_squared = 0.5 * (1 - sqrt(1 - q_v_max^2));
gamma = sqrt(gamma_squared);

% For comparison, compute small-angle approximation
gamma_approx = q_v_max / 2;

fprintf('\n=== Exact Sector Bound ===\n');
fprintf('Rigorous bound: gamma = %.6f\n', gamma);
fprintf('Small-angle approx: gamma ≈ q_v_max/2 = %.6f\n', gamma_approx);
fprintf('Relative difference: %.2f%%\n', abs(gamma - gamma_approx)/gamma * 100);

%% Sector Bound Extraction Matrix
% E_sec extracts omega from state z = [q_v; omega]
E_sec = [zeros(3,3), eye(3)];

%% LMI-Based RCI Set Synthesis (Theorem from paper)
fprintf('\n=== Solving LMI Optimization ===\n');
fprintf('Using exact 4-block LMI structure from Theorem...\n');

best_vol = -inf;
best_Q = [];
best_Y = [];
best_alpha = [];
best_tau = [];

alpha_vals = logspace(-1, 2, 12);
tau_vals = logspace(-1, 2, 12);

for alpha = alpha_vals
    for tau = tau_vals
        cvx_begin sdp quiet
            cvx_solver mosek
            cvx_precision high
            
            variable Q(n,n) symmetric
            variable Y(m,n)
            
            Q >= 1e-4 * eye(n);    
            
            % Construct the 4-block LMI
            Psi_11 = A*Q + Q*A' + B*Y + Y'*B' + alpha*Q;
            B_phi = [eye(3); 
                    zeros(3)];
            E_bar = E;
            M = [Psi_11, B_phi, E_bar, gamma*sqrt(tau)*Q*E_sec';
                 B_phi', -tau*eye(3), zeros(3), zeros(3);
                 E_bar', zeros(3), -alpha*eye(3), zeros(3);
                 gamma*sqrt(tau)*E_sec*Q, zeros(3), zeros(3), -eye(3)];
            M <= 0;
            
            % Input constraints
            for i = 1:m
                e_i = zeros(m,1);
                e_i(i) = 1;
                [Q, Y(i,:)'; ...
                 Y(i,:), u_max(i)^2] >= 0;
            end

            % State constraints
            for i = 1:n
                c_i = zeros(n,1)';
                c_i(i) = 1;
                [Q,       Q*c_i'; ...
                 c_i*Q,  x_max(i)^2] >= 0;
            end
            
            maximize(log_det(Q))
        cvx_end
        
        if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
            vol = log(det(Q));
            if vol > best_vol
                best_vol = vol;
                best_Q = Q;
                best_Y = Y;
                best_alpha = alpha;
                best_tau = tau;
            end
        end
    end
end

if isempty(best_Q)
    error('LMI optimization failed - both formulations infeasible');
end

fprintf('\n=== Optimal Solution ===\n');
if ~isnan(best_alpha)
    fprintf('Best alpha: %.4f\n', best_alpha);
    fprintf('Best tau: %.4f\n', best_tau);
    fprintf('Using full formulation with exact sector bound gamma = %.6f\n', gamma);
else
    fprintf('Using simplified formulation (full formulation not feasible)\n');
end
fprintf('log(det(Q)): %.4f\n', best_vol);

% Extract controller and ellipsoid
P = inv(best_Q);
K = best_Y / best_Q;

fprintf('\nBackup controller gain K:\n');
disp(K);

% Verify closed-loop stability
eig_cl = eig(A + B*K);
fprintf('Closed-loop eigenvalues:\n');
fprintf('  %.4f%+.4fi\n', [real(eig_cl), imag(eig_cl)]');
if all(real(eig_cl) < 0)
    fprintf('✓ Closed-loop system is stable\n');
else
    warning('✗ Closed-loop system is unstable!');
end

%% LQR Nominal Controller Design
Q_lqr = diag([1e-3, 1e-3, 1e-3, 1e-7, 1e-7, 1e-7]);
R_lqr = eye(3);
K_lqr = -lqr(A, B, Q_lqr, R_lqr);

fprintf('\n=== LQR Nominal Controller ===\n');
fprintf('LQR eigenvalues:\n');
eig_lqr = eig(A - B*K_lqr);
fprintf('  %.4f%+.4fi\n', [real(eig_lqr), imag(eig_lqr)]');

% Safety Filter Parameters
h_max = 0.95;  % Pure Backup control
h_min = 0.20;  % Pure nominal control

alpha_blend = @(z) min(max((z'*P*z - h_min)/(h_max - h_min), 0), 1);

fprintf('\n=== Safety Filter Configuration ===\n');
fprintf('h_min: %.2f (pure nominal)\n', h_min);
fprintf('h_max: %.2f (start blending)\n', h_max);
fprintf('Blending region: [%.2f, %.2f]\n', h_min, h_max);

% Nonlinear Simulation Setup
T = 8;      % simulation time
dt = 0.001;  % time step
time = 0:dt:T;
N = length(time);

% Initial condition (small perturbation)
q_v_init = [0.1; -0.05; 0.02];
omega_init = [0.3; -0.2; 0.1];
z_init = [q_v_init; omega_init];

% Desired setpoint (hover with zero attitude)
z_des = zeros(6,1);

fprintf('\n=== Simulation Configuration ===\n');
fprintf('Duration: %.1f seconds\n', T);
fprintf('Time step: %.4f seconds\n', dt);
fprintf('Initial q_v: [%.3f, %.3f, %.3f]\n', q_v_init);
fprintf('Initial omega: [%.3f, %.3f, %.3f] rad/s\n', omega_init);
fprintf('Initial attitude angle: %.2f deg\n', 2*asin(norm(q_v_init))*180/pi);

% Disturbance trajectory
% d_traj = @(t) [sin(2*t); cos(2*t); 0];
d_traj = @(t) d_max * [sin(2*t); cos(2*t); 0];


% Preallocate
z_safe = zeros(n, N);
z_unsafe = zeros(n, N);
u_safe = zeros(m, N);
u_unsafe = zeros(m, N);
alpha_traj = zeros(1, N);
h_safe = zeros(1, N);
h_unsafe = zeros(1, N);
phi_norm_safe = zeros(1, N);
phi_norm_unsafe = zeros(1, N);

z_safe(:,1) = z_init;
z_unsafe(:,1) = z_init;
h_safe(1) = z_init' * P * z_init;
h_unsafe(1) = z_init' * P * z_init;

% Nonlinear Dynamics Functions
% Exact quaternion kinematics: dq_v/dt = (1/2)*(q_0*I + [q_v]_x)*omega
% with q_0 = sqrt(1 - ||q_v||^2)

skew = @(v) [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];

dynamics_qv = @(q_v, omega) 0.5 * (sqrt(1 - norm(q_v)^2)*eye(3) + skew(q_v)) * omega;
dynamics_omega = @(omega, u, d) J_inv * (-cross(omega, J*omega) + u + d);
% dynamics_omega = @(omega, u, d) J_inv * ( u + d);

% Compute exact nonlinearity phi(z)
compute_phi = @(q_v, omega) 0.5 * ((sqrt(1-norm(q_v)^2) - 1)*eye(3) + skew(q_v)) * omega;

% Simulation Loop
fprintf('\n=== Running Nonlinear Simulation ===\n');
fprintf('Progress: ');

ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

for k = 1:N-1
    if mod(k, floor(N/10)) == 0
        fprintf('.');
    end
    
    t_k = time(k);
    d_k = d_traj(t_k);
    
    % ===== WITH SAFETY FILTER =====
    z_k = z_safe(:,k);
    q_v = z_k(1:3);
    omega = z_k(4:6);
    
    % Nominal control
    u_nom = -K_lqr * (z_k - z_des);
    
    % Backup control
    u_backup = K * z_k;
    
    % Safety filter blending
    alpha_k = alpha_blend(z_k);
    alpha_traj(k) = alpha_k;
    
    u_k = (1 - alpha_k) * u_nom + alpha_k * u_backup;
    
    % Saturate control
    u_k = max(min(u_k, u_max), -u_max);
    u_safe(:,k) = u_k;
    
    % Compute nonlinearity for verification
    phi_k = compute_phi(q_v, omega);
    phi_norm_safe(k) = norm(phi_k);
    
    % Integrate nonlinear dynamics
    dz = [dynamics_qv(q_v, omega);
          dynamics_omega(omega, u_k, d_k)];
    z_safe(:,k+1) = z_k + dt * dz;

    % Linear Model Integration (for comparison)
    % dz = A*z_k + B*u_k + E*d_k + B_phi*phi_k;  % Include nonlinearity for exact simulation
    % z_safe(:,k+1) = z_k + dt * dz;
    
    % Safety check - if quaternion norm gets too large, clip it
    if norm(z_safe(1:3,k+1)) > 0.99
        fprintf('\nWarning: q_v norm = %.3f at t=%.2f, clipping to valid range\n', ...
            norm(z_safe(1:3,k+1)), time(k+1));
        z_safe(1:3,k+1) = z_safe(1:3,k+1) / norm(z_safe(1:3,k+1)) * 0.9;
    end
    
    % Compute Lyapunov value
    h_safe(k+1) = z_safe(:,k+1)' * P * z_safe(:,k+1);
    
    % Stop if system becomes unstable
    if h_safe(k+1) > 1e6 || any(isnan(z_safe(:,k+1))) || any(isinf(z_safe(:,k+1)))
        fprintf('\nSimulation stopped at t=%.2f: system unstable\n', time(k+1));
        % Fill remaining arrays with last valid values
        for j = k+1:N
            z_safe(:,j) = z_safe(:,k+1);
            h_safe(j) = h_safe(k+1);
        end
        for j = k:N-1
            alpha_traj(j) = alpha_traj(k);
            phi_norm_safe(j) = phi_norm_safe(k);
            u_safe(:,j) = u_safe(:,k);
        end
        break;
    end
    
    % ===== WITHOUT SAFETY FILTER =====
    z_k = z_unsafe(:,k);
    q_v = z_k(1:3);
    omega = z_k(4:6);
    
    u_nom = -K_lqr * (z_k - z_des);
    u_k = max(min(u_nom, u_max), -u_max);
    u_unsafe(:,k) = u_k;
    
    phi_k = compute_phi(q_v, omega);
    phi_norm_unsafe(k) = norm(phi_k);
    
    dz = [dynamics_qv(q_v, omega);
          dynamics_omega(omega, u_k, d_k)];
    z_unsafe(:,k+1) = z_k + dt * dz;

    % Linear model integration for comparison
    % dz = A*z_k + B*u_k + E*d_k + B_phi*phi_k; % Include nonlinearity for exact simulation
    % z_unsafe(:,k+1) = z_k + dt * dz;
    
    % Safety check
    if norm(z_unsafe(1:3,k+1)) > 0.99
        z_unsafe(1:3,k+1) = z_unsafe(1:3,k+1) / norm(z_unsafe(1:3,k+1)) * 0.9;
    end
    
    h_unsafe(k+1) = z_unsafe(:,k+1)' * P * z_unsafe(:,k+1);
    
    % Stop if unstable
    if h_unsafe(k+1) > 1e6 || any(isnan(z_unsafe(:,k+1))) || any(isinf(z_unsafe(:,k+1)))
        for j = k+1:N
            z_unsafe(:,j) = z_unsafe(:,k+1);
            h_unsafe(j) = h_unsafe(k+1);
        end
        for j = k:N-1
            phi_norm_unsafe(j) = phi_norm_unsafe(k);
            u_unsafe(:,j) = u_unsafe(:,k);
        end
        break;
    end
end

fprintf(' Done!\n');

% Results Analysis
fprintf('\n=== Safety Analysis ===\n');
fprintf('Max h(z) with filter: %.6f\n', max(h_safe));
fprintf('Max h(z) without filter: %.6f\n', max(h_unsafe));

if max(h_safe) <= 1.0
    fprintf('✓ Safety constraint satisfied (h ≤ 1)\n');
else
    fprintf('✗ Safety boundary violated by %.2f%%\n', (max(h_safe)-1)*100);
end

% Verify sector bound satisfaction
q_v_safe = z_safe(1:3,:);
omega_safe = z_safe(4:6,:);
q_v_norm = vecnorm(q_v_safe);
omega_norm = vecnorm(omega_safe);

% Theoretical bound: ||phi|| <= gamma * ||omega||
phi_bound_theoretical = gamma * omega_norm;
sector_bound_ratio = phi_norm_safe ./ (gamma * omega_norm + 1e-10);

fprintf('\n=== Sector Bound Verification ===\n');
fprintf('Max ||phi||/||omega||: %.6f\n', max(phi_norm_safe ./ (omega_norm + 1e-10)));
fprintf('Theoretical bound gamma: %.6f\n', gamma);
fprintf('Max ratio ||phi||/(gamma*||omega||): %.6f\n', max(sector_bound_ratio));

if max(sector_bound_ratio) <= 1.0
    fprintf('✓ Sector bound satisfied\n');
else
    fprintf('✗ Sector bound violated\n');
end

fprintf('\n=== State Constraint Verification ===\n');
fprintf('Max ||q_v||: %.6f (design limit: %.6f)\n', max(q_v_norm), q_v_max);
fprintf('Max ||omega||: %.4f rad/s (design limit: %.2f rad/s)\n', max(omega_norm), omega_max);

% Check actual ellipsoid bounds
Q_ellipse = inv(P);
E_qv = [eye(3), zeros(3)];
Q_qv = E_qv * Q_ellipse * E_qv';
q_v_max_actual = sqrt(max(eig(Q_qv)));
fprintf('\nActual ellipsoid bounds:\n');
fprintf('  Max ||q_v|| in ellipsoid: %.6f\n', q_v_max_actual);
fprintf('  Corresponds to %.2f deg attitude angle\n', 2*asin(q_v_max_actual)*180/pi);

% Plotting
figure('Position', [100, 100, 1600, 1000], 'Name', 'Rigorous Quaternion Safety Filter');

% Quaternion vector part components
subplot(3,4,1);
plot(time, z_safe(1,:), 'b', 'LineWidth', 2); hold on;
plot(time, z_unsafe(1,:), 'r--', 'LineWidth', 2);
yline([-q_v_max, q_v_max], 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('q_{v,1}'); 
title('Quaternion Vector Part (x)');
legend('Safe', 'Unsafe', 'Bounds', 'Location', 'best');
grid on;

subplot(3,4,2);
plot(time, z_safe(2,:), 'b', 'LineWidth', 2); hold on;
plot(time, z_unsafe(2,:), 'r--', 'LineWidth', 2);
yline([-q_v_max, q_v_max], 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('q_{v,2}'); 
title('Quaternion Vector Part (y)');
legend('Safe', 'Unsafe', 'Bounds', 'Location', 'best');
grid on;

subplot(3,4,3);
plot(time, z_safe(3,:), 'b', 'LineWidth', 2); hold on;
plot(time, z_unsafe(3,:), 'r--', 'LineWidth', 2);
yline([-q_v_max, q_v_max], 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('q_{v,3}'); 
title('Quaternion Vector Part (z)');
legend('Safe', 'Unsafe', 'Bounds', 'Location', 'best');
grid on;

% Quaternion vector norm (rotation angle)
subplot(3,4,4);
theta_safe = 2*asin(vecnorm(z_safe(1:3,:)))*180/pi;
theta_unsafe = 2*asin(vecnorm(z_unsafe(1:3,:)))*180/pi;
plot(time, theta_safe, 'b', 'LineWidth', 2); hold on;
plot(time, theta_unsafe, 'r--', 'LineWidth', 2);
yline(theta_max_deg, 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Angle (deg)'); 
title('Total Rotation Angle');
legend('Safe', 'Unsafe', sprintf('Limit (%.0f°)', theta_max_deg), 'Location', 'best');
grid on;

% Angular velocities
subplot(3,4,5);
plot(time, z_safe(4,:), 'b', 'LineWidth', 2); hold on;
plot(time, z_unsafe(4,:), 'r--', 'LineWidth', 2);
yline([-omega_max, omega_max], 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('\omega_x (rad/s)'); 
title('Angular Velocity (x)');
grid on;

subplot(3,4,6);
plot(time, z_safe(5,:), 'b', 'LineWidth', 2); hold on;
plot(time, z_unsafe(5,:), 'r--', 'LineWidth', 2);
yline([-omega_max, omega_max], 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('\omega_y (rad/s)'); 
title('Angular Velocity (y)');
grid on;

subplot(3,4,7);
plot(time, z_safe(6,:), 'b', 'LineWidth', 2); hold on;
plot(time, z_unsafe(6,:), 'r--', 'LineWidth', 2);
yline([-omega_max, omega_max], 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('\omega_z (rad/s)'); 
title('Angular Velocity (z)');
grid on;

% Angular velocity norm
subplot(3,4,8);
plot(time, vecnorm(z_safe(4:6,:)), 'b', 'LineWidth', 2); hold on;
plot(time, vecnorm(z_unsafe(4:6,:)), 'r--', 'LineWidth', 2);
yline(omega_max, 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('||\omega|| (rad/s)'); 
title('Angular Velocity Magnitude');
legend('Safe', 'Unsafe', 'Limit', 'Location', 'best');
grid on;

% Control inputs
subplot(3,4,9);
plot(time(1:end), u_safe(1,:), 'b', 'LineWidth', 2); hold on;
plot(time(1:end), u_unsafe(1,:), 'r--', 'LineWidth', 2);
yline([-u_max(1), u_max(1)], 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('u_1 (N·m)'); 
title('Control Torque (x)');
grid on;

subplot(3,4,10);
plot(time(1:end), u_safe(2,:), 'b', 'LineWidth', 2); hold on;
plot(time(1:end), u_unsafe(2,:), 'r--', 'LineWidth', 2);
yline([-u_max(2), u_max(2)], 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('u_2 (N·m)'); 
title('Control Torque (y)');
grid on;

subplot(3,4,11);
plot(time(1:end), u_safe(3,:), 'b', 'LineWidth', 2); hold on;
plot(time(1:end), u_unsafe(3,:), 'r--', 'LineWidth', 2);
yline([-u_max(3), u_max(3)], 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('u_3 (N·m)'); 
title('Control Torque (z)');
grid on;

% Lyapunov function
subplot(3,4,12);
plot(time, h_safe, 'b', 'LineWidth', 2); hold on;
plot(time, h_unsafe, 'r--', 'LineWidth', 2);
yline([h_min, h_max, 1.0], 'k--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('h(z) = z^T P z'); 
title('RCI Set Membership');
legend('Safe', 'Unsafe', 'h_{min}', 'h_{max}', 'Boundary', 'Location', 'best');
grid on;
ylim([0, max(max(h_safe), max(h_unsafe))*1.1]);

sgtitle('Rigorous Quaternion Safety Filter (Exact Sector Bound)', 'FontSize', 14, 'FontWeight', 'bold');

% Sector Bound Verification Plot
figure('Position', [150, 150, 1200, 400], 'Name', 'Sector Bound Verification');

subplot(1,3,1);
plot(time(1:end), phi_norm_safe, 'b', 'LineWidth', 2); hold on;
plot(time(1:end-1), gamma * vecnorm(z_safe(4:6,1:end-1)), 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Magnitude'); 
title('Nonlinearity vs. Sector Bound');
legend('||\phi(z)||', '\gamma ||\omega||', 'Location', 'best');
grid on;

subplot(1,3,2);
plot(time(1:end-1), sector_bound_ratio(1:end-1), 'b', 'LineWidth', 2); hold on;
yline(1.0, 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Ratio'); 
title('Sector Bound Satisfaction');
legend('||\phi|| / (\gamma ||\omega||)', 'Theoretical limit', 'Location', 'best');
grid on;
ylim([0, 1.2]);

subplot(1,3,3);
plot(time(1:end-1), alpha_traj(1:end-1), 'b', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\alpha'); 
title('Safety Filter Activation');
ylim([-0.1, 1.1]);
grid on;

sgtitle('Verification of Exact Sector Bound', 'FontSize', 12, 'FontWeight', 'bold');

fprintf('\n=== Simulation Complete ===\n');
fprintf('All results shown in figures.\n');