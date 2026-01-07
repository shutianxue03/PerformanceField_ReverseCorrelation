function basicFxn_drawCorr_permutation( ...
    X_allIter_allSubj, Y_allIter_allSubj, colors, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, type_corr, type_tail, str_title, markers_allSubj, nIter)
% =========================================================================
% basicFxn_drawCorr_permutation
% Inputs
%   X_allIter_allSubj : 3D array 
%   Y_allIter_allSubj : 3D array 
%       * May be in any dimension order, but must contain (nIter, nSubj, nCond)
%       * Internally reshaped to canonical: [nIter x nSubj x nCond]
%   colors            : [nCond x 3] RGB color per condition
%   x_ticks, y_ticks  : tick locations (or NaN to skip setting)
%   x_ticklabels, y_ticklabels : tick labels (or NaN to skip setting)
%   flag_zeroMean     : 0 = no centering
%                       1 = subtract mean across observers within each condition
%                       2 = subtract mean across conditions within each observer
%   flag_plotIdvdCI   : 1 = draw per-point CI bars (X horizontal, Y vertical)
%   flag_plotUnikSymbol : 1 = each observer uses a unique marker from markers_allSubj
%                         0 = all observers use 'o'
%   type_corr         : correlation type for corr / partialcorr (e.g., 'Pearson','Spearman')
%   type_tail         : tail for significance tests ('both','right','left')
%   str_title         : base title string (will be appended with results)
%   markers_allSubj   : cell array (nSubj x 1) of marker symbols per observer (e.g., 'o','s','+','x')
%   nIter             : number of iterations 
%
% Output
%   str_sig           : string label based on median partial-corr p-value:
%                       '_sig' (<=.05), '_mg' (<=.1), '_ns' otherwise
% =========================================================================

%% 0) Plot style defaults (edit once, reuse forever)
wd_border = 3;     % axes + marker linewidth
sz_ticks  = 45;    % tick label size
sz_marker = 20;    % marker size

CI_level = .95; % CI range for stats only! Default for plotting is .68

%% 1) Validate / canonicalize shapes: enforce [nIter x nSubj x nCond]
assert(ndims(X_allIter_allSubj) == 3 && ndims(Y_allIter_allSubj) == 3, ...
    'ALERT: X_allIter_allSubj and Y_allIter_allSubj must be 3D arrays.');

nCond = size(colors, 1);
nSubj = length(markers_allSubj);

X_allIter_allSubj = fxn_reshape(X_allIter_allSubj, nSubj, nCond, 'X_allIter_allSubj');
Y_allIter_allSubj = fxn_reshape(Y_allIter_allSubj, nSubj, nCond, 'Y_allIter_allSubj');

[nIter2, nSubj2, nCond2] = size(X_allIter_allSubj);
assert(all(size(Y_allIter_allSubj) == [nIter2, nSubj2, nCond2]), ...
    'ALERT: X and Y must match after reshape/permute.');
assert(nSubj2 == nSubj && nCond2 == nCond, ...
    'ALERT: Reshaped X has wrong nSubj/nCond.');
assert(nIter2 == nIter, ...
    'ALERT: Input nIter=%d but inferred nIter=%d from X.', nIter, nIter2);

%% 2) Bootstrap summary for plotting: median + CI per (observer × condition)
% Output arrays are [nSubj x nCond]
[X_med_allSubj, X_lb_allSubj, X_ub_allSubj] = getCI(X_allIter_allSubj, 1, 1);
[Y_med_allSubj, Y_lb_allSubj, Y_ub_allSubj] = getCI(Y_allIter_allSubj, 1, 1);

%% 3) Optional centering (flag_zeroMean): define plotting-space variables
% We apply the same shift to medians and CI bounds so errorbars match points.
switch flag_zeroMean
    case 1 % subtract mean across observers within each condition
        shiftX = mean(X_med_allSubj, 1);   % 1 x nCond
        shiftY = mean(Y_med_allSubj, 1);   % 1 x nCond

        X_med_plot = X_med_allSubj - shiftX;
        Y_med_plot = Y_med_allSubj - shiftY;

        X_lb_plot = X_lb_allSubj - shiftX;
        X_ub_plot = X_ub_allSubj - shiftX;
        Y_lb_plot = Y_lb_allSubj - shiftY;
        Y_ub_plot = Y_ub_allSubj - shiftY;

    case 2 % subtract mean across conditions within each observer
        shiftX = mean(X_med_allSubj, 2);   % nSubj x 1
        shiftY = mean(Y_med_allSubj, 2);   % nSubj x 1

        X_med_plot = X_med_allSubj - shiftX;
        Y_med_plot = Y_med_allSubj - shiftY;

        X_lb_plot = X_lb_allSubj - shiftX;
        X_ub_plot = X_ub_allSubj - shiftX;
        Y_lb_plot = Y_lb_allSubj - shiftY;
        Y_ub_plot = Y_ub_allSubj - shiftY;

    otherwise % 0 (no centering)
        X_med_plot = X_med_allSubj;
        Y_med_plot = Y_med_allSubj;

        X_lb_plot = X_lb_allSubj;  X_ub_plot = X_ub_allSubj;
        Y_lb_plot = Y_lb_allSubj;  Y_ub_plot = Y_ub_allSubj;
end

%% 4) Scatter plot of medians (optionally with per-point CI bars)
figure('Position', [0 200 1e3 1e3]);
hold on; box on

for iCond = 1:nCond
    for iSubj = 1:nSubj

        % Choose marker style (unique per subject vs. uniform circles)
        if flag_plotUnikSymbol
            style_marker = markers_allSubj{iSubj};
            % Some markers do not have a fillable face (e.g., '+', 'x')
            if ismember(style_marker, {'+','x','*','.'})
                color_face = 'none';
                color_edge = colors(iCond, :);
            else
                color_face = colors(iCond, :);
                color_edge = 'w';
            end
        else
            style_marker = 'o';
            color_face = colors(iCond, :);
            color_edge = 'w';
        end

        x0 = X_med_plot(iSubj, iCond);
        y0 = Y_med_plot(iSubj, iCond);

        % Optional: individual CIs (horizontal for X, vertical for Y)
        if flag_plotIdvdCI == 1
            plot([X_lb_plot(iSubj,iCond), X_ub_plot(iSubj,iCond)], [y0, y0], '-', ...
                'Color', ones(1,3)/2, 'LineWidth', 2, 'HandleVisibility','off');
            plot([x0, x0], [Y_lb_plot(iSubj,iCond), Y_ub_plot(iSubj,iCond)], '-', ...
                'Color', ones(1,3)/2, 'LineWidth', 2, 'HandleVisibility','off');
        end

        % The point itself
        plot(x0, y0, style_marker, ...
            'MarkerFaceColor', color_face, ...
            'MarkerEdgeColor', color_edge, ...
            'MarkerSize', sz_marker, ...
            'LineWidth', wd_border);
    end
end

%% 5) Prepare indices for "controlling for condition" analyses
ANOVA_indCond = repmat(1:nCond, nSubj, 1);  % nSubj x nCond

%% 6) Regression on MEDIANS
nSamples = 200;

% Global regression on medians (all conditions pooled)
x_global = linspace(min(X_med_plot(:)), max(X_med_plot(:)), nSamples);
polyfit_global_OnMed = polyfit(X_med_plot(:), Y_med_plot(:), 1);
yfit_global_OnMed = polyval(polyfit_global_OnMed, x_global);

% Per-condition regression on medians (within each condition)
x_allCond = nan(nCond, nSamples);
yfit_allCond_OnMed = nan(nCond, nSamples);
for iCond = 1:nCond
    xmin = min(X_med_plot(:, iCond));
    xmax = max(X_med_plot(:, iCond));
    if xmin == xmax
        x_allCond(iCond, :) = linspace(xmin - 1e-6, xmax + 1e-6, nSamples);
    else
        x_allCond(iCond, :) = linspace(xmin, xmax, nSamples);
    end
    polyfit_perCond_OnMed = polyfit(X_med_plot(:, iCond), Y_med_plot(:, iCond), 1);
    yfit_allCond_OnMed(iCond, :) = polyval(polyfit_perCond_OnMed, x_allCond(iCond, :));
end

%% 7) Regression PER BOOT
yfit_global_allIter   = nan(nIter, nSamples);          % global regression line per iteration
yfit_perCond_allIter  = nan(nIter, nSamples, nCond);   % per-condition regression line per iteration
eta2_global_allIter   = nan(nIter, 1);                 % eta^2 per iteration (global)

r_partial_allIter = nan(nIter, 1); % partial corr (global association controlling for condition)
p_partial_allIter = nan(nIter, 1);

r_perCond_allIter = nan(nIter, nCond);
p_perCond_allIter = nan(nIter, nCond);

for iIter = 1:nIter
    x_perBoot = squeeze(X_allIter_allSubj(iIter, :, :));   % nSubj x nCond
    y_perBoot = squeeze(Y_allIter_allSubj(iIter, :, :));   % nSubj x nCond

    % Regression uses the same coordinate system as the scatter (centering applied).
    switch flag_zeroMean
        case 1
            x_reg = x_perBoot - shiftX;   % shiftX: 1 x nCond
            y_reg = y_perBoot - shiftY;   % shiftY: 1 x nCond
        case 2
            x_reg = x_perBoot - shiftX;   % shiftX: nSubj x 1
            y_reg = y_perBoot - shiftY;   % shiftY: nSubj x 1
        otherwise
            x_reg = x_perBoot;
            y_reg = y_perBoot;
    end

    % --- Correlations ---
    % Partial correlation controlling for condition.
    [r_partial_allIter(iIter), p_partial_allIter(iIter)] = partialcorr( ...
        x_perBoot(:), y_perBoot(:), ANOVA_indCond(:), 'type', type_corr, 'tail', type_tail);

    % Correlation within each condition (across observers)
    for iCond = 1:nCond
        [r_perCond_allIter(iIter, iCond), p_perCond_allIter(iIter, iCond)] = corr( ...
            x_perBoot(:, iCond), y_perBoot(:, iCond), 'type', type_corr, 'tail', type_tail);
    end

    % --- Regressions ---
    % Global
    if numel(unique(x_reg(:))) >= 2
        polyfit_global_perBoot = polyfit(x_reg(:), y_reg(:), 1);
        yfit_global_allIter(iIter, :) = polyval(polyfit_global_perBoot, x_global);
        eta2_global_allIter(iIter) = var(polyval(polyfit_global_perBoot, x_reg(:))) / var(y_reg(:));
    end

    % Per condition
    for iCond = 1:nCond
        polyfit_perCond_perBoot = polyfit(x_reg(:, iCond), y_reg(:, iCond), 1);
        yfit_perCond_allIter(iIter, :, iCond) = polyval(polyfit_perCond_perBoot, x_allCond(iCond, :));
    end
end

%% 8) Summarize distributions (median + CI)
[r_med_partial, r_lb_partial, r_ub_partial] = getCI(r_partial_allIter, 1, 1, CI_level);
[p_med_partial, p_lb_partial, p_ub_partial] = getCI(p_partial_allIter, 1, 1, CI_level);
[eta2_med, eta2_lb, eta2_ub] = getCI(eta2_global_allIter, 1, 1, CI_level);

[r_med_perCond, r_lb_perCond, r_ub_perCond] = getCI(r_perCond_allIter, 1, 1, CI_level);
[p_med_perCond, p_lb_perCond, p_ub_perCond] = getCI(p_perCond_allIter, 1, 1, CI_level);

%% 9) Plot regression lines (from medians) + CI bands (from per iteration)
% Global
% if p_med_partial < .1
    [~, yfit_global_lb, yfit_global_ub] = getCI(yfit_global_allIter, 1, 1);
    patch([x_global fliplr(x_global)], [yfit_global_lb fliplr(yfit_global_ub)], ...
        ones(1,3)*.6, 'EdgeColor','none', 'FaceAlpha', 0.20, 'HandleVisibility','off');
    plot(x_global, yfit_global_OnMed, '-', 'color', ones(1,3)*.4, ...
        'HandleVisibility','off', 'LineWidth', wd_border*1.5);
% end

% Per condition
for iCond = 1:nCond
    if p_med_perCond(iCond) < .1
        lineStyle = '-';
    else
        lineStyle = '--';
    end
    [~, yfit_perCond_lb, yfit_perCond_ub] = getCI(yfit_perCond_allIter(:, :, iCond), 1, 1);
    patch([x_allCond(iCond, :) fliplr(x_allCond(iCond, :))], [yfit_perCond_lb fliplr(yfit_perCond_ub)], ...
        colors(iCond,:), 'EdgeColor','none', 'FaceAlpha', 0.12, 'HandleVisibility','off');
    plot(x_allCond(iCond, :), yfit_allCond_OnMed(iCond, :), lineStyle, 'Color', colors(iCond,:), ...
        'LineWidth', wd_border, 'HandleVisibility','off');
end % iCond

%% 10) Axes formatting
if ~isnan(x_ticks), xticks(x_ticks), xlim(x_ticks([1, end])), end
if ~isnan(y_ticks), yticks(y_ticks), ylim(y_ticks([1, end])), end
if ~isnan(x_ticklabels), xticks(x_ticklabels), end
if ~isnan(y_ticklabels), yticks(y_ticklabels), end

axis square
ax = gca;
ax.XAxis.FontSize = sz_ticks;
ax.YAxis.FontSize = sz_ticks;
ax.LineWidth = wd_border;

%% 11) Title text
% str_sig = '_ns';
% if p_med_partial <= .05
%     str_sig = '_sig';
% elseif p_med_partial <= .1
%     str_sig = '_mg';
% end

str_partial = sprintf('Partial r=%.2f [%.2f, %.2f], eta2=%.2f [%.2f, %.2f]', ...
    r_med_partial, r_lb_partial, r_ub_partial, eta2_med, eta2_lb, eta2_ub);

str_perCond = "";
for iCond = 1:nCond
    if mod(iCond-1, 3) == 0 && iCond > 1
        str_perCond = str_perCond + sprintf('\n');
    end
    str_perCond = str_perCond + sprintf('Cond#%d: r=%.2f [%.2f, %.2f] | ', ...
        iCond, r_med_perCond(iCond), r_lb_perCond(iCond), r_ub_perCond(iCond));
end

title(sprintf('%s\n%s\n%s', str_title, str_partial, str_perCond));

end
