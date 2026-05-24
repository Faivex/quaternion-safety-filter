function plot_conservatism_sweep(data_path, save_dir)
%PLOT_CONSERVATISM_SWEEP  Two-panel figure summarizing the 2D LMI sweep.
%
%   plot_conservatism_sweep('results/conservatism/conservatism_sweep_data.mat')
%   plot_conservatism_sweep(data_path, save_dir)
%
%   Panels:
%     (a) RCI volume coverage (Proposed): vol(E(P)) / vol(box) as a
%         percentage of the constraint-box volume.  Infeasible cells are
%         shown in grey.
%     (b) Volume ratio  vol_Proposed / vol_Linearized  over the overlap
%         region only.  Non-overlap cells are grey (= linearized
%         infeasible).  Colormap is logarithmic but tick labels show real
%         multipliers (e.g. 1x, 5x).
%
%   Reads precomputed data; no SDP solve required.

    if nargin < 1 || isempty(data_path)
        data_path = fullfile('results','conservatism','conservatism_sweep_data.mat');
    end
    if nargin < 2 || isempty(save_dir)
        save_dir = fileparts(data_path);
    end
    if ~exist(data_path, 'file')
        error('Sweep data not found: %s\nRun conservatism_analysis.m first.', data_path);
    end

    S = load(data_path);
    theta_sweep = S.theta_sweep;
    omega_sweep = S.omega_sweep;
    logdet_full = S.logdet_full;
    logdet_lin  = S.logdet_lin;
    feas_full   = S.feas_full;
    feas_lin    = S.feas_lin;

    %% Compute real-unit quantities
    n_state = 6;
    ellipsoid_const = pi^(n_state/2) / gamma(n_state/2 + 1);

    [TG, OG] = meshgrid(theta_sweep, omega_sweep);
    qv_max_grid  = sin(TG * pi/180 / 2);
    vol_box_grid = (2*qv_max_grid).^3 .* (2*OG).^3;

    vol_full = ellipsoid_const * exp(logdet_full/2);
    vol_lin  = ellipsoid_const * exp(logdet_lin /2);

    coverage_full = 100 * vol_full ./ vol_box_grid;
    ratio_overlap = vol_full ./ vol_lin;

    %% Palette
    col_lmi  = [0.220, 0.424, 0.690];
    col_grey = [0.85,  0.85,  0.85];

    n_cmap   = 256;
    cmap_lmi = make_sequential([0.95 0.97 1.00], col_lmi*0.85, n_cmap);

    %% Figure
    fig = figure('Color', 'w'); fig.Units = 'inches';
    fig.Position(3:4) = [7.0, 3.0];
    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');

    % ----- Panel (a): Proposed RCI coverage -----
    ax1 = nexttile;
    h1 = imagesc(theta_sweep, omega_sweep, coverage_full);
    set(h1, 'AlphaData', double(feas_full));
    set(ax1, 'Color', col_grey, 'YDir', 'normal');
    colormap(ax1, cmap_lmi);
    cb1 = colorbar(ax1);
    cb1.Label.Interpreter = 'latex';
    cb1.Label.String = 'RCI coverage of constraint box (\%)';
    xlabel('$\theta_{\max}$ (deg)', 'Interpreter', 'latex');
    ylabel('$\omega_{\max}$ (rad/s)', 'Interpreter', 'latex');
    title('(a) Proposed RCI coverage', ...
        'Interpreter', 'latex', 'FontWeight', 'normal');
    set_axes_ticks(ax1, theta_sweep, omega_sweep);
    axis tight;

    % ----- Panel (b): volume ratio (multiplicative) -----
    ax2 = nexttile;
    overlap = feas_full & feas_lin;
    ratio_plot = ratio_overlap;
    ratio_plot(~overlap) = NaN;

    log_ratio = log10(ratio_plot);
    rmax = max(log_ratio(overlap), [], 'all');
    if isempty(rmax) || isnan(rmax) || rmax <= 0
        rmax = 1;
    end

    h2 = imagesc(theta_sweep, omega_sweep, log_ratio);
    set(h2, 'AlphaData', double(overlap));
    set(ax2, 'Color', col_grey, 'YDir', 'normal');
    colormap(ax2, cmap_lmi);
    caxis(ax2, [0, rmax]);
    cb2 = colorbar(ax2);
    cb2.TickLabelInterpreter = 'latex';
    tick_log    = linspace(0, rmax, 5);
    tick_real   = 10.^tick_log;
    cb2.Ticks   = tick_log;
    cb2.TickLabels = arrayfun(@(r) format_ratio(r), tick_real, 'UniformOutput', false);
    cb2.Label.Interpreter = 'latex';
    cb2.Label.String = 'Volume ratio  (Proposed / Linearized)';

    xlabel('$\theta_{\max}$ (deg)', 'Interpreter', 'latex');
    ylabel('$\omega_{\max}$ (rad/s)', 'Interpreter', 'latex');
    title('(b) Volume gain (overlap region)', ...
        'Interpreter', 'latex', 'FontWeight', 'normal');
    set_axes_ticks(ax2, theta_sweep, omega_sweep);
    axis tight;

    % Legend in panel (b): grey cells = linearized infeasible
    hold(ax2, 'on');
    h_lin_inf = patch(NaN, NaN, col_grey, 'EdgeColor', 'none', ...
        'DisplayName', 'Linearized infeasible');
    legend(ax2, h_lin_inf, ...
        'Location', 'southoutside', 'Orientation', 'horizontal', ...
        'Box', 'off', 'Interpreter', 'latex');

    set(findall(fig, '-property', 'FontSize'), 'FontSize', 9);

    %% Save
    if ~isempty(save_dir)
        if ~exist(save_dir, 'dir'); mkdir(save_dir); end
        outfile = fullfile(save_dir, 'conservatism_2D.pdf');
        try
            exportgraphics(fig, outfile, 'ContentType', 'vector', 'Resolution', 1000);
            fprintf('Figure saved to %s\n', outfile);
        catch ME
            warning('Could not save figure: %s');
        end
    end

    %% Console summary
    n_grid = numel(feas_full);
    fprintf('\n--- Conservatism summary ---\n');
    fprintf('Grid size: %d points (theta %d-%d deg, omega %g-%g rad/s)\n', ...
        n_grid, min(theta_sweep), max(theta_sweep), ...
        min(omega_sweep), max(omega_sweep));
    fprintf('Proposed   feasible: %d / %d (%.0f%%)\n', ...
        sum(feas_full(:)), n_grid, 100*sum(feas_full(:))/n_grid);
    fprintf('Linearized feasible: %d / %d (%.0f%%)\n', ...
        sum(feas_lin(:)),  n_grid, 100*sum(feas_lin(:))/n_grid);
    if any(overlap(:))
        fprintf('Volume ratio (Proposed/Lin) over overlap:  median %sx, max %sx\n', ...
            format_ratio(median(ratio_overlap(overlap))), ...
            format_ratio(max(ratio_overlap(overlap), [], 'all')));
    end
end

%% ---- helpers ----
function s = format_ratio(r)
    if r >= 100
        s = sprintf('%.0f', r);
    elseif r >= 10
        s = sprintf('%.1f', r);
    else
        s = sprintf('%.2f', r);
    end
end

function cmap = make_sequential(c_lo, c_hi, n)
    cmap = [linspace(c_lo(1), c_hi(1), n)', ...
            linspace(c_lo(2), c_hi(2), n)', ...
            linspace(c_lo(3), c_hi(3), n)'];
end

function set_axes_ticks(ax, xs, ys)
    set(ax, 'XTick', xs, 'YTick', ys);
end