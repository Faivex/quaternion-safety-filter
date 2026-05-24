function plotSafetyResults(time, x_safe, x_unsafe, x_cbf, ...
                                u_safe, u_unsafe, u_cbf, ...
                                x_ref, q_v_max, omega_max, u_max, ...
                                qp_infeas_log, slack_log, savePath)
%PLOTSAFETYRESULTSTHREE  Three-way comparison figure: LMI vs Unsafe vs HOCBF.
%
%   Six-panel layout (2x3):
%   (1,1) 2D trajectory with time-gradient colouring
%   (1,2) ||q_v||      vs constraint q_v_max
%   (1,3) ||omega||    vs constraint omega_max
%   (2,1) ||tau||_inf  vs constraint u_max
%   (2,2) h_0(z) = q_v_max^2 - ||q_v||^2   (HOCBF barrier; <0 = unsafe)
%   (2,3) slack variable xi(t) with infeasibility markers

    if nargin < 14, savePath = ''; end

    %% Colour palette (Paul Tol "bright", colourblind-safe, print-friendly)
    col_lmi    = [0.220, 0.424, 0.690];   % blue
    col_cbf    = [0.945, 0.490, 0.196];   % vivid orange (was green)
    col_unsafe = [0.847, 0.110, 0.337];   % crimson
    col_ref    = [0.345, 0.345, 0.345];   % medium grey for refs/goals
    col_lim    = [0.000, 0.000, 0.000];   % black for constraint lines
    col_thresh = [0.620, 0.000, 0.620];   % purple (slack threshold)
    col_start  = [0.450, 0.450, 0.450];   % medium grey (paired with col_ref)

    old_lw = get(0, 'DefaultLineLineWidth');
    old_fs = get(0, 'DefaultAxesFontSize');
    set(0, 'DefaultLineLineWidth', 1.5);
    set(0, 'DefaultAxesFontSize', 12);
    cleanup = onCleanup(@() restore_defaults(old_lw, old_fs));

    %% Derived signals
    Q_safe   = eul2quat(x_safe(7:9,:)',   'XYZ');
    Q_unsafe = eul2quat(x_unsafe(7:9,:)', 'XYZ');
    Q_cbf    = eul2quat(x_cbf(7:9,:)',    'XYZ');

    qv_safe   = Q_safe(:,2:4)';
    qv_unsafe = Q_unsafe(:,2:4)';
    qv_cbf    = Q_cbf(:,2:4)';
    qv_norm_safe   = vecnorm(qv_safe,   2, 1);
    qv_norm_unsafe = vecnorm(qv_unsafe, 2, 1);
    qv_norm_cbf    = vecnorm(qv_cbf,    2, 1);

    om_norm_safe   = vecnorm(x_safe(10:12,:),   2, 1);
    om_norm_unsafe = vecnorm(x_unsafe(10:12,:), 2, 1);
    om_norm_cbf    = vecnorm(x_cbf(10:12,:),    2, 1);

    u_inf_safe   = max(abs(u_safe),   [], 1);
    u_inf_unsafe = max(abs(u_unsafe), [], 1);
    u_inf_cbf    = max(abs(u_cbf),    [], 1);

    h0_safe = q_v_max^2 - qv_norm_safe.^2;
    h0_cbf  = q_v_max^2 - qv_norm_cbf.^2;

    %% Slack: use the same threshold for marker logic and dashed line
    SLACK_THRESH = 1e-3;
    if isempty(slack_log)
        slack_log = zeros(size(time));
    end
    slack_log = slack_log(:)';
    slack_plot = slack_log(1:length(time));
    finite_max = max([slack_plot(isfinite(slack_plot)), SLACK_THRESH * 10]);
    slack_plot(~isfinite(slack_plot)) = finite_max;

    %% Line widths and styles
    lw_main   = 1.7;   % LMI
    lw_cbf    = 1.5;   % CBF
    lw_unsafe = 1.4;   % Unsafe (bumped up so it's visible)
    style_unsafe = '-.';  % dash-dot (more visible than dotted)

    %% Figure setup
    fig = figure('Color', 'w', 'Name', 'Three-Way Comparison', 'Visible', 'on');
    fig.Units = 'inches';
    fig.Position(3:4) = [7.2, 3.9];
    set(fig, 'DefaultAxesFontSize', 9);
    tl = tiledlayout(2, 3, 'TileSpacing', 'tight', 'Padding', 'tight');

    % --- (1,1) 2D trajectory with time-gradient colouring ---
    nexttile;
    pos_safe   = x_safe(1:2,:);
    pos_unsafe = x_unsafe(1:2,:);
    pos_cbf    = x_cbf(1:2,:);
    cx = mean(pos_safe(1,:));
    cy = mean(pos_safe(2,:));
    half_x = 1.5 * max(max(abs(pos_safe(1,:) - cx)), 0.1);
    half_y = 1.5 * max(max(abs(pos_safe(2,:) - cy)), 0.1);

    % Unsafe reference
    h_nominal_traj = plot(pos_unsafe(1,:), pos_unsafe(2,:), style_unsafe, ...
        'Color', [col_unsafe, 0.7], 'LineWidth', lw_unsafe, ...
        'HandleVisibility', 'off'); hold on;

    % Reference (goal or trajectory)
    if size(x_ref, 2) > 1
        h_ref = plot(x_ref(1,:), x_ref(2,:), '-', 'Color', col_ref, ...
            'LineWidth', 1.0, 'DisplayName', 'Ref.');
    else
        h_ref = plot(x_ref(1), x_ref(2), 'p', 'MarkerSize', 12, ...
            'MarkerFaceColor', col_ref, 'MarkerEdgeColor', 'k', ...
            'LineWidth', 1, 'DisplayName', 'Goal');
    end

    % Time-gradient trajectories (~200 segments)
    N = length(time);
    stride = max(1, floor(N/200));
    idx = 1:stride:N;
    if idx(end) < N, idx(end+1) = N; end

    % CBF gradient: light orange -> deep orange-brown
    cbf_cmap = make_gradient([0.99, 0.75, 0.45], col_cbf * 0.85, length(idx));
    for k = 1:length(idx)-1
        seg = idx(k):idx(k+1);
        plot(pos_cbf(1, seg), pos_cbf(2, seg), '-', ...
            'Color', cbf_cmap(k,:), 'LineWidth', lw_cbf, 'HandleVisibility','off');
    end
    % LMI gradient: light blue -> deep navy
    lmi_cmap = make_gradient([0.55, 0.72, 0.92], col_lmi * 0.75, length(idx));
    for k = 1:length(idx)-1
        seg = idx(k):idx(k+1);
        plot(pos_safe(1, seg), pos_safe(2, seg), '-', ...
            'Color', lmi_cmap(k,:), 'LineWidth', lw_main, 'HandleVisibility','off');
    end

    % Legend proxies (mid-tone) for the figure-wide legend
    h_proxy_proposed = plot(nan, nan, '-', 'Color', col_lmi, 'LineWidth', lw_main, ...
        'DisplayName', 'Proposed');
    h_proxy_hocbf    = plot(nan, nan, '-', 'Color', col_cbf, 'LineWidth', lw_cbf, ...
        'DisplayName', 'HOCBF');
    h_proxy_nominal  = plot(nan, nan, style_unsafe, 'Color', col_unsafe, ...
        'LineWidth', lw_unsafe, 'DisplayName', 'Nominal');

    % Start marker
    plot(pos_safe(1,1), pos_safe(2,1), 'o', 'MarkerSize', 7, ...
        'MarkerFaceColor', col_start, 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.8, 'HandleVisibility','off');

    xlabel('$x$ (m)', 'Interpreter', 'latex');
    ylabel('$y$ (m)', 'Interpreter', 'latex');
    xlim([cx - half_x, cx + half_x]);
    ylim([cy - half_y, cy + half_y]);
    grid on;

    % --- (1,2) ||q_v|| ---
    nexttile;
    plot(time, qv_norm_unsafe, style_unsafe, 'Color', [col_unsafe, 0.7], ...
        'LineWidth', lw_unsafe, 'HandleVisibility', 'off'); hold on;
    plot(time, qv_norm_cbf, '-', 'Color', col_cbf, 'LineWidth', lw_cbf, ...
        'HandleVisibility', 'off');
    plot(time, qv_norm_safe, '-', 'Color', col_lmi, 'LineWidth', lw_main, ...
        'HandleVisibility', 'off');
    yl = yline(q_v_max, '--', '$\bar{q}_v$', 'Color', col_lim, 'LineWidth', 2.0);
    yl.Interpreter = 'latex'; yl.LabelHorizontalAlignment = 'left';
    yl.LabelVerticalAlignment = 'bottom'; yl.HandleVisibility = 'off';
    xlabel('Time (s)');
    ylabel('$\|q_v\|$', 'Interpreter', 'latex');
    ylim_top = min(max([qv_norm_safe, qv_norm_cbf]) * 1.15, 1.05);
    ylim_top = max(ylim_top, q_v_max * 1.5);
    ylim([0, ylim_top]);
    xlim([time(1), time(end)]);
    grid on;

    % --- (1,3) ||omega|| ---
    nexttile;
    plot(time, om_norm_unsafe, style_unsafe, 'Color', [col_unsafe, 0.7], ...
        'LineWidth', lw_unsafe, 'HandleVisibility', 'off'); hold on;
    plot(time, om_norm_cbf, '-', 'Color', col_cbf, 'LineWidth', lw_cbf, ...
        'HandleVisibility', 'off');
    plot(time, om_norm_safe, '-', 'Color', col_lmi, 'LineWidth', lw_main, ...
        'HandleVisibility', 'off');
    yl = yline(omega_max, '--', '$\omega_{\max}$', 'Color', col_lim, 'LineWidth', 2.0);
    yl.Interpreter = 'latex'; yl.LabelHorizontalAlignment = 'left';
    yl.LabelVerticalAlignment = 'bottom'; yl.HandleVisibility = 'off';
    xlabel('Time (s)');
    ylabel('$\|\omega\|$ (rad/s)', 'Interpreter', 'latex');
    om_top = min(max([om_norm_safe, om_norm_cbf]) * 1.2, omega_max * 2);
    om_top = max(om_top, omega_max * 1.3);
    ylim([0, om_top]);
    xlim([time(1), time(end)]);
    grid on;

    % --- (2,1) ||tau||_inf ---
    nexttile;
    t_u = time(1:end-1);
    plot(t_u, u_inf_unsafe(1:end-1), style_unsafe, 'Color', [col_unsafe, 0.7], ...
        'LineWidth', lw_unsafe, 'HandleVisibility', 'off'); hold on;
    plot(t_u, u_inf_cbf(1:end-1), '-', 'Color', col_cbf, 'LineWidth', lw_cbf, ...
        'HandleVisibility', 'off');
    plot(t_u, u_inf_safe(1:end-1), '-', 'Color', col_lmi, 'LineWidth', lw_main, ...
        'HandleVisibility', 'off');
    yl = yline(u_max(1), '--', '$\tau_{\max}$', 'Color', col_lim, 'LineWidth', 2.0);
    yl.Interpreter = 'latex'; yl.LabelHorizontalAlignment = 'left';
    yl.LabelVerticalAlignment = 'bottom'; yl.HandleVisibility = 'off';
    xlabel('Time (s)');
    ylabel('$\|\tau\|_\infty$ (N$\cdot$m)', 'Interpreter', 'latex');
    ylim([0, u_max(1)*1.5]);
    xlim([t_u(1), t_u(end)]);
    grid on;

    % --- (2,2) HOCBF barrier h_0 ---
    nexttile;
    plot(time, h0_cbf, '-', 'Color', col_cbf, 'LineWidth', lw_cbf, ...
        'HandleVisibility', 'off'); hold on;
    plot(time, h0_safe, '-', 'Color', col_lmi, 'LineWidth', lw_main, ...
        'HandleVisibility', 'off');
    yl = yline(0, '--', '$h_0 = 0$', 'Color', col_lim, 'LineWidth', 2.0);
    yl.Interpreter = 'latex'; yl.LabelHorizontalAlignment = 'left';
    yl.LabelVerticalAlignment = 'bottom'; yl.HandleVisibility = 'off';
    xlabel('Time (s)');
    ylabel('$h_0(z)$', 'Interpreter', 'latex');
    all_h0 = [h0_safe, h0_cbf];
    h0_lo = min(min(all_h0), -0.05);
    h0_hi = max(max(all_h0), q_v_max^2) * 1.1;
    ylim([h0_lo, h0_hi]);
    xlim([time(1), time(end)]);
    grid on;

    % --- (2,3) Slack variable xi(t) ---
    ax_slack = nexttile;
    slack_for_plot = max(slack_plot, 1e-10);
    t_s   = time(1:end-1);
    sfp_s = slack_for_plot(1:end-1);
    h_xi = semilogy(t_s, sfp_s, '-', 'Color', col_cbf, 'LineWidth', lw_cbf, ...
        'DisplayName', '$\xi(t)$'); hold on;
    h_thr = yline(SLACK_THRESH, '--', 'Color', col_thresh, 'LineWidth', 1.4, ...
        'DisplayName', 'Infeas. Threshold');
    panel_handles = [h_xi, h_thr];
    xlabel('Time (s)');
    ylabel('$\xi$', 'Interpreter', 'latex');
    active = sfp_s > SLACK_THRESH * 0.1;
    if any(active)
        t_lo = t_s(find(active, 1, 'first'));
        t_hi = t_s(find(active, 1, 'last'));
        margin = 0.05 * (t_s(end) - t_s(1));
        xlim([max(t_s(1), t_lo - margin), min(t_s(end), t_hi + margin)]);
    else
        xlim([t_s(1), t_s(end)]);
    end
    ylim_lo = 1e-10;
    ylim_hi = max([finite_max * 5, SLACK_THRESH * 100]);
    ylim([ylim_lo, ylim_hi]);
    lgd_slack = legend(ax_slack, panel_handles, ...
        'Location', 'southeast', 'Interpreter', 'latex', ...
        'Box', 'off', 'NumColumns', 1);
    lgd_slack.ItemTokenSize = [8, 10];
    grid on;

    %% Figure-wide shared legend (method curves + reference/goal)
    method_handles = [h_proxy_proposed, h_proxy_hocbf, h_proxy_nominal];
    if exist('h_ref', 'var') && ishghandle(h_ref)
        method_handles(end+1) = h_ref;
    end
    % Build legend on the trajectory axes (where the proxy handles live),
    % then reposition it across the bottom of the figure in normalized
    % figure coordinates.  This is more reliable across MATLAB versions
    % than tl.Layout.Tile = 'south'.
    ax_first = ancestor(h_proxy_proposed, 'axes');
    lgd_main = legend(ax_first, method_handles, ...
        'Orientation', 'horizontal', 'Interpreter', 'latex', 'Box', 'off');
    lgd_main.ItemTokenSize = [14, 14];
    lgd_main.Units = 'normalized';
    drawnow;   % force layout pass so actual Position is finalized
    lgd_pos = lgd_main.Position;
    lgd_main.Position = [(1 - lgd_pos(3))/2, 0.005, lgd_pos(3), lgd_pos(4)];
    % Make room at the bottom of the tiled layout for the legend
    tl.Padding = 'compact';
    set(tl, 'OuterPosition', [0, 0.07, 1, 0.93]);

    %% Save BEFORE optional formatting hook so PDF always gets written
    if ~isempty(savePath)
        if ~exist(savePath, 'dir')
            mkdir(savePath);
        end
        outfile = fullfile(savePath, 'quaternion_safety_results_three.pdf');
        try
            exportgraphics(fig, outfile, 'ContentType', 'vector', 'Resolution', 1000);
            fprintf('Three-way comparison figure saved to: %s\n', outfile);
        catch ME
            warning('Failed to save figure: %s');
            % Fallback: try the older print mechanism
            try
                set(fig, 'PaperPositionMode', 'auto');
                print(fig, outfile, '-dpdf', '-vector');
                fprintf('Saved (fallback): %s\n', outfile);
            catch
                fprintf('Could not save figure to %s\n', outfile);
            end
        end
    end

    % if exist('formatFigureIEEE', 'file')
    %     try, formatFigureIEEE(); catch, end
    % end
end

%% ---------- helpers ----------
function g = make_gradient(c_start, c_end, n)
    % Linear RGB interpolation
    g = [linspace(c_start(1), c_end(1), n)', ...
         linspace(c_start(2), c_end(2), n)', ...
         linspace(c_start(3), c_end(3), n)'];
end

function restore_defaults(lw, fs)
    set(0, 'DefaultLineLineWidth', lw);
    set(0, 'DefaultAxesFontSize', fs);
end
