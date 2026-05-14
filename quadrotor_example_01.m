clc;
clear;
close all;

%% Linearize the system
% Load system parameters
params = loadSystemParameters();

J = diag([params.Ixx, params.Iyy, params.Izz]);
J_inv = inv(J);

fprintf('=== System Parameters ===\n');
fprintf('Mass: %.4f kg\n', params.m);
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
theta_max_deg = 40;  % degrees
theta_max_rad = theta_max_deg * pi / 180;
q_v_max = sin(theta_max_rad / 2);

% Angular velocity constraint
omega_max = 5.0;  % rad/s

% State constraints: [q_v; omega]
x_max = [q_v_max; q_v_max; q_v_max; omega_max; omega_max; omega_max];

% Control input constraints
u_max = [1e-4; 1e-4; 1e-4];  % N·m

fprintf('\n=== Safety Constraints ===\n');
fprintf('Max attitude angle: %.1f deg\n', theta_max_deg);
fprintf('q_v_max = sin(%.1f°/2) = %.6f\n', theta_max_deg, q_v_max);
fprintf('omega_max: %.2f rad/s\n', omega_max);

%% Exact Sector Bound Computation
% gamma^2 = (1/2) * (1 - sqrt(1 - q_v_max^2))
% This is the EXACT bound, no approximations

gamma_squared = 0.5 * (1 - sqrt(1 - q_v_max^2));
gamma = sqrt(gamma_squared);

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

alpha_vals = logspace(-1, 2, 5);
tau_vals = logspace(-1, 2, 5);

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
                [Q, Y(i,:)'; ...
                 Y(i,:), u_max(i)^2] >= 0;
            end

            % State constraints
            C_1 = [eye(3), zeros(3)];  % Extract q_v part of state
            [Q,       Q*C_1'; ...
            C_1*Q,  q_v_max^2*eye(3)] >= 0;

            C_2 = [zeros(3), eye(3)];  % Extract omega part of state
            [Q,       Q*C_2'; ...
            C_2*Q,  omega_max^2*eye(3)] >= 0;

                 
            
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

%% Simulation
% Initial condition
x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % initial state
% Simulation time
T = 5; % total time
dt = 0.01; % time step
time = 0:dt:T; % time vector
% Preallocate state and control input arrays for both scenarios
x_safe = zeros(length(x0), length(time));
u_safe = zeros(m, length(time));
x_unsafe = zeros(length(x0), length(time));
u_unsafe = zeros(m, length(time));
x_safe(:,1) = x0; % set initial condition
x_unsafe(:,1) = x0; % set initial condition

% Define alpha as a function handle before simulation loop
h_max = 0.9;
h_min = 0.1;
alpha_fun = @(x) min(max((x'*P*x - h_min)/(h_max - h_min), 0), 1);

% Define Goal (Explicitly for plotting)
x_des = 2;
y_des = 1;
x_goal = [x_des; y_des; 0]; % Target Position

% Simple PD controller for position to generate desired roll and pitch
Kp_y = -0.2; % position proportional gain
Kd_y = -0.2; % position derivative gain
Kp_x = 0.2;  % position proportional gain
Kd_x = 0.2;  % position derivative gain

% Simulation loop
for k = 1:length(time)-1
    % Run both scenarios in parallel
    for scenario = 1:2
        if scenario == 1
            x_current = x_safe(:,k);
        else
            x_current = x_unsafe(:,k);
        end
        
        vx = x_current(4);
        vy = x_current(5);
        ex = x_des - x_current(1);
        ey = y_des - x_current(2);
        theta_des = Kp_x*ex - Kd_x*vx;
        phi_des = Kp_y*ey - Kd_y*vy;

        % Limit desired angles
        phi_des = max(min(phi_des, 60*pi/180), -60*pi/180);
        theta_des = max(min(theta_des, 60*pi/180), -60*pi/180);

        % Attitude control
        eta = x_current(7:9); % current angles [phi; theta; psi]
        eta_des = [phi_des; theta_des; 0]; % desired angles
        e_eta = eta - eta_des; % angle error
        omega = x_current(10:12); % current angular rates [p; q; r]
        k_p_phi = -1e-3; k_d_phi = 2e-4;
        k_p_theta = -1e-3; k_d_theta = 2e-4;
        k_p_psi = -3e-4; k_d_psi = 1e-4;
        u_nominal = [k_p_phi*e_eta(1) + k_d_phi*(0 - omega(1));
                    k_p_theta*e_eta(2) + k_d_theta*(0 - omega(2));
                    k_p_psi*e_eta(3) + k_d_psi*(0 - omega(3))];

        % Saturate control inputs
        for i = 1:length(u_nominal)
            if abs(u_nominal(i)) > u_max(i)
                u_nominal(i) = sign(u_nominal(i)) * u_max(i);
            end
        end
        
        if scenario == 1
            % With safety filter
            % Compute quaternion and angular velocity part of state for alpha computation
            % quaternion = eul2quat(x_current(7:9)');
            quat = eul2quat(x_current(7:9)', 'XYZ');
            x_filtered = [quat(2:4)'; x_current(10:12)]; % extract q_v and omega
            alpha = alpha_fun(x_filtered);
            u_current = (alpha)*K*x_filtered + (1 - alpha)*u_nominal;
        else
            % Without safety filter
            u_current = u_nominal;
        end

        % Add coriolis and gyroscopic effects compensation to the control input
        p = x_current(10);
        q = x_current(11);
        r = x_current(12);
        coriolis_comp = [params.Iyy - params.Izz; params.Izz - params.Ixx; params.Ixx - params.Iyy] .* [q*r; p*r; p*q];
        u_current = u_current + coriolis_comp;

        % % Saturate control inputs
        % for i = 1:m
        %     if abs(u_current(i)) > u_max(i)
        %         u_current(i) = sign(u_current(i)) * u_max(i);
        %     end
        % end

        % Store saturated control inputs
        if scenario == 1
            u_safe(:,k) = u_current;
        else
            u_unsafe(:,k) = u_current;
        end

        % Define disturbance values
        w = @(t) [sin(2*t); cos(2*t); 0]; % define w as a function of time
        w_t = w(time(k)); % evaluate w at current time

        % Simulate system
        ode_fun = @(t, state) F(state, params) + G(state, params)*u_current + [zeros(6,3);E]*w_t;
        [~, x_next] = ode45(ode_fun, [0 dt], x_current);
        
        if scenario == 1
            x_safe(:,k+1) = x_next(end,:)';
        else
            x_unsafe(:,k+1) = x_next(end,:)';
        end
    end
end

plotSafetyResults(time, x_safe, x_unsafe, u_safe, u_unsafe, x_goal, q_v_max, u_max, P, h_min, h_max, 'results/quadrotor_example_01');

% Functions
function F_x = F(state, params)
    % Dynamics function F
    % state = [x; y; z; vx; vy; vz; phi; theta; psi; p; q; r]

    % Extract state variables
    vx = state(4);
    vy = state(5);
    vz = state(6);
    phi = state(7);
    theta = state(8);
    psi = state(9);
    p = state(10);
    q = state(11);
    r = state(12);

    % Gravitational acceleration
    g = params.g;
    m = params.m;

    % State derivatives
    % dot{p} = v 
    p_dot = [vx; vy; vz];

    % Compute trigonometric functions
    c_phi = cos(phi);
    s_phi = sin(phi);
    c_theta = cos(theta);
    s_theta = sin(theta);
    c_psi = cos(psi);
    s_psi = sin(psi);

    % Rotation matrix (body to inertial)
    C_IB = [c_theta*c_psi, c_theta*s_psi, -s_theta;
            -c_phi*s_psi + s_phi*s_theta*c_psi, c_phi*c_psi + s_phi*s_theta*s_psi, s_phi*c_theta;
            s_phi*s_psi + c_phi*s_theta*c_psi, -s_phi*c_psi + c_phi*s_theta*s_psi, c_phi*c_theta];
    % dot{v} = g*e_3 - \frac{1}{m}R(\eta)(Fe_3)
    e3 = [0; 0; 1];
    v_dot = g*e3 - (1/m)*C_IB*(m*g*e3); % Assuming Fe_3 = mg*e_3 at hover

    % dot{\eta} = W(\eta)\omega
    W = [1, s_phi*tan(theta), c_phi*tan(theta);
         0, c_phi, -s_phi;
         0, s_phi/cos(theta), c_phi/cos(theta)];
    eta_dot = W * [p; q; r];

    % dot{\omega} = J^{-1}(-\omega \times J\omega)
    J = diag([params.Ixx, params.Iyy, params.Izz]);
    omega = [p; q; r];
    omega_dot = J \ (-cross(omega, J*omega)); % Assuming no control

    % Combine all derivatives
    F_x = [p_dot; v_dot; eta_dot; omega_dot];
end

function G_x = G(~, params)
    % Control input matrix G
    Ixx = params.Ixx;
    Iyy = params.Iyy;
    Izz = params.Izz;

    % Control inputs affect acceleration and angular acceleration
    G_x = zeros(12, 3);

    % Torques affect angular accelerations
    G_x(10, 1) = 1/Ixx; % Roll torque
    G_x(11, 2) = 1/Iyy; % Pitch torque
    G_x(12, 3) = 1/Izz; % Yaw torque
end