function plotSafetyResults(time, x_safe, x_unsafe, u_safe, u_unsafe, x_ref, q_v_max, u_max, P, h_min, h_max, savePath)
    % plotSafetyResults - Visualizes simulation results for quaternion safety filter
    % x_ref:    3x1 static goal OR 2xN/3xN trajectory history
    % savePath: (optional) directory for PDF export of the paper figure

    if nargin < 12, savePath = ''; end

    %% Color palette (shared by all figures)
    c1 = [0,    0.447, 0.741];   % blue
    c2 = [0.85, 0.325, 0.098];   % red
    c3 = [0.929, 0.694, 0.125];  % gold
    c4 = [0.494, 0.184, 0.556];  % purple
    c5 = [0.466, 0.674, 0.188];  % green
    lc = [0.3, 0.3, 0.3];        % constraint/limit lines

    % Scope default overrides to avoid leaking into unrelated figures
    old_lw = get(0, 'DefaultLineLineWidth');
    old_fs = get(0, 'DefaultAxesFontSize');
    set(0, 'DefaultLineLineWidth', 1.5);
    set(0, 'DefaultAxesFontSize', 12);
    cleanup = onCleanup(@() restore_defaults(old_lw, old_fs));

    %% Pre-compute all derived signals (single vectorized pass)

    % Quaternion conversion: input is Nx3 (each row = [phi theta psi])
    Q_safe   = eul2quat(x_safe(7:9,:)',   'XYZ');  % Nx4, columns: [w x y z]
    Q_unsafe = eul2quat(x_unsafe(7:9,:)', 'XYZ');

    qv_safe        = Q_safe(:,2:4)';    % 3xN
    qv_unsafe      = Q_unsafe(:,2:4)';  % 3xN
    qv_norm_safe   = vecnorm(qv_safe,   2, 1);  % 1xN
    qv_norm_unsafe = vecnorm(qv_unsafe, 2, 1);  % 1xN

    omega_safe = x_safe(10:12, :);      % 3xN
    Z_safe     = [qv_safe; omega_safe]; % 6xN — reduced state

    % Lyapunov function V(z) = z'Pz and mixing coefficient alpha
    Lyap_safe     = sum((P * Z_safe) .* Z_safe, 1);
    alpha_monitor = min(max((Lyap_safe - h_min) / (h_max - h_min), 0), 1);

    % Sector-bounded nonlinearity phi(z) = 0.5*((q0-1)*omega + qv x omega)
    gamma   = sqrt(0.5 * (1 - sqrt(1 - q_v_max^2)));
    q0_safe = Q_safe(:,1)';  % 1xN
    phi_safe        = 0.5 * ((q0_safe - 1) .* omega_safe + cross(qv_safe, omega_safe));
    phi_norm_safe   = vecnorm(phi_safe,   2, 1);  % 1xN
    omega_norm_safe = vecnorm(omega_safe, 2, 1);  % 1xN

    % Sector bound ratio (epsilon avoids 0/0 at rest)
    sector_bound_ratio = phi_norm_safe ./ (gamma * omega_norm_safe + 1e-9);

    % Euler angle constraint limit derived from q_v_max (total rotation bound)
    theta_lim_deg = 2 * asin(q_v_max) * 180/pi;

    %% COMBINED PAPER FIGURE (2x3)
    % Row 1: Position  |  Euler angles  |  Control inputs (tau1, tau2)
    % Row 2: ||q_v||   |  V(z) + alpha  |  Sector bound ratio

    fig_paper = figure('Color', 'w', 'Name', 'Paper Figure');
    fig_paper.Units = 'inches';
    fig_paper.Position(3:4) = [7, 3.2];
    set(fig_paper, 'DefaultAxesFontSize', 9);
    tiledlayout(2, 3, 'TileSpacing', 'tight', 'Padding', 'tight');

    % (1,1) 2D X-Y trajectory
    nexttile;
    pos_safe   = x_safe(1:2, :);
    pos_unsafe = x_unsafe(1:2, :);
    cx = mean(pos_safe(1,:));
    cy = mean(pos_safe(2,:));
    half_x = 1.5 * max(max(abs(pos_safe(1,:) - cx)), 0.1);
    half_y = 1.5 * max(max(abs(pos_safe(2,:) - cy)), 0.1);
    plot(pos_unsafe(1,:), pos_unsafe(2,:), ':', 'Color', [c2, 0.5], 'LineWidth', 0.8, 'DisplayName', 'Unsafe'); hold on;
    if size(x_ref, 2) > 1
        plot(x_ref(1,:), x_ref(2,:), '-', 'Color', lc, 'LineWidth', 1, 'DisplayName', 'Ref.');
    else
        plot(x_ref(1), x_ref(2), 'p', 'MarkerSize', 10, 'MarkerFaceColor', 'y', ...
            'MarkerEdgeColor', 'k', 'DisplayName', 'Goal');
    end
    plot(pos_safe(1,:), pos_safe(2,:), '-', 'Color', c1, 'LineWidth', 1.5, 'DisplayName', 'Safe');
    plot(pos_safe(1,1), pos_safe(2,1), 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.2 0.7 0.2], ...
        'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    xlabel('$x$ (m)', 'Interpreter', 'latex');
    ylabel('$y$ (m)', 'Interpreter', 'latex');
    xlim([cx - half_x, cx + half_x]);
    ylim([cy - half_y, cy + half_y]);
    lgd = legend('Location', 'southoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'Box', 'off');
    lgd.ItemTokenSize = [10, 12];
    grid on;

    % (1,2) Euler angles phi, theta, psi
    nexttile;
    angles_safe   = x_safe(7:9, :);
    angles_unsafe = x_unsafe(7:9, :);
    plot(time, angles_unsafe(1,:)*180/pi, ':', 'Color', [c1, 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off'); hold on;
    plot(time, angles_unsafe(2,:)*180/pi, ':', 'Color', [c2, 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(time, angles_unsafe(3,:)*180/pi, ':', 'Color', [c3, 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    yline( theta_lim_deg, '--', 'Color', lc, 'LineWidth', 1, 'DisplayName', 'Constraint');
    yline(-theta_lim_deg, '--', 'Color', lc, 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(time, angles_safe(1,:)*180/pi, '-', 'Color', c1, 'LineWidth', 1.5, 'DisplayName', '$\phi$');
    plot(time, angles_safe(2,:)*180/pi, '-', 'Color', c2, 'LineWidth', 1.5, 'DisplayName', '$\theta$');
    plot(time, angles_safe(3,:)*180/pi, '-', 'Color', c3, 'LineWidth', 1.5, 'DisplayName', '$\psi$');
    safe_ang_max = 1.5 * max(max(abs(angles_safe * 180/pi)));
    ang_lim = max(safe_ang_max, theta_lim_deg * 1.1);
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    xlim([time(1), time(end)]);
    ylim([-ang_lim, ang_lim]);
    lgd = legend('Location', 'southoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'Box', 'off');
    lgd.ItemTokenSize = [10, 12];
    grid on;

    % (1,3) Control inputs tau1, tau2, tau3
    nexttile;
    plot(time, u_unsafe(1,:), ':', 'Color', [c1, 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off'); hold on;
    plot(time, u_unsafe(2,:), ':', 'Color', [c2, 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(time, u_unsafe(3,:), ':', 'Color', [c3, 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    yline( u_max(1), '--', 'Color', lc, 'LineWidth', 1, 'DisplayName', 'Constraint');
    yline(-u_max(1), '--', 'Color', lc, 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(time, u_safe(1,:), '-', 'Color', c1, 'LineWidth', 1.5, 'DisplayName', '$\tau_1$');
    plot(time, u_safe(2,:), '-', 'Color', c2, 'LineWidth', 1.5, 'DisplayName', '$\tau_2$');
    plot(time, u_safe(3,:), '-', 'Color', c3, 'LineWidth', 1.5, 'DisplayName', '$\tau_3$');
    xlabel('Time (s)');
    ylabel('N$\cdot$m', 'Interpreter', 'latex');
    ylim([-u_max(1)*1.3, u_max(1)*1.3]);
    xlim([time(1), time(end)]);
    lgd = legend('Location', 'southoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'Box', 'off');
    lgd.ItemTokenSize = [10, 12];
    grid on;

    % (2,1) Attitude constraint ||q_v||
    nexttile;
    plot(time, qv_norm_unsafe, ':', 'Color', [c2, 0.6], 'LineWidth', 0.8, 'DisplayName', 'Unsafe'); hold on;
    plot(time, qv_norm_safe,   '-', 'Color', c1,        'LineWidth', 1.5, 'DisplayName', 'Safe');
    yline(q_v_max, '--', 'Color', lc, 'LineWidth', 1.5, 'DisplayName', '$\bar{q}_v$');
    xlabel('Time (s)');
    ylabel('$\|\mathbf{q}_v\|$', 'Interpreter', 'latex');
    ylim([0, q_v_max * 1.5]);
    xlim([time(1), time(end)]);
    lgd = legend('Location', 'southoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'Box', 'off');
    lgd.ItemTokenSize = [10, 12];
    grid on;

    % (2,2) Lyapunov V(z) and mixing alpha
    nexttile;
    plot(time, Lyap_safe,     '-', 'Color', c4, 'LineWidth', 1.5, 'DisplayName', '$V(\mathbf{z})$'); hold on;
    plot(time, alpha_monitor, '-', 'Color', c5, 'LineWidth', 1.5, 'DisplayName', '$\alpha$');
    yline(1.0, '--', 'Color', lc, 'LineWidth', 1.5, 'DisplayName', 'RCI boundary');
    xlabel('Time (s)');
    ylabel('Value');
    ylim([0, max(max(Lyap_safe), 1.0) * 1.3]);
    xlim([time(1), time(end)]);
    lgd = legend('Location', 'southoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'Box', 'off');
    lgd.ItemTokenSize = [10, 12];
    grid on;

    % (2,3) Sector bound ratio
    nexttile;
    plot(time, sector_bound_ratio, '-', 'Color', c1, 'LineWidth', 1.5, ...
        'DisplayName', '$\|\phi\|/(\gamma\|\omega\|)$'); hold on;
    yline(1.0, '--', 'Color', lc, 'LineWidth', 1.5, 'DisplayName', 'Bound');
    xlabel('Time (s)');
    ylabel('Ratio');
    ylim([0, 1.25]);
    xlim([time(1), time(end)]);
    lgd = legend('Location', 'southoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'Box', 'off');
    lgd.ItemTokenSize = [10, 12];
    grid on;

    formatFigureIEEE();

    if ~isempty(savePath)
        if ~exist(savePath, 'dir'), mkdir(savePath); end
        exportgraphics(fig_paper, fullfile(savePath, 'quaternion_safety_results.pdf'), ...
            'ContentType', 'vector', 'Resolution', 1000);
        fprintf('Paper figure saved to: %s\n', fullfile(savePath, 'quaternion_safety_results.pdf'));
    end
end

function restore_defaults(lw, fs)
    set(0, 'DefaultLineLineWidth', lw);
    set(0, 'DefaultAxesFontSize', fs);
end
