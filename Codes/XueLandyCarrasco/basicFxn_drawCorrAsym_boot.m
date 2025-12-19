function basicFxn_drawCorrAsym_boot(X_allBoot_allSubj, Y_allBoot_allSubj, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)

% =========================================================================
% basicFxn_drawCorrAsym_boot.m
%
% Purpose
%   Plot and compare *asymmetry* metrics between two measurements (X vs Y)
%   across observers, while propagating bootstrap uncertainty.
%
% Inputs (canonical / recommended)
%   asymX_allBoot_allSubj : [nBoot x nSubj] bootstrap samples for X asymmetry
%   asymY_allBoot_allSubj : [nBoot x nSubj] bootstrap samples for Y asymmetry
%
%   x_ticks, y_ticks       : vectors (e.g., 5 values) for tick locations
%   x_ticklabels, y_ticklabels : tick labels (numeric or cellstr)
%
%   flag_plotIdvdCI      : 1 plot per-point CI bars; 0 plot points only
%   flag_plotUnikSymbol  : 1 use subject-specific markers; 0 all circles
%   type_corr            : 'pearson' | 'spearman' | 'kendall'
%   type_tail            : 'both' | 'left' | 'right'
%   text_title           : title prefix string
%   markers_allSubj      : {nSubj x 1} marker strings (e.g., {'o','s',...})
%
% Outputs
%   str_sig   : '_sig' / '_mg' / '_ns' based on median p-value across boots
%   str_title : full title string including r/p medians and CIs
%   stats     : struct with medians, CIs, and bootstrap distributions
%
% Notes
%   - Uses getCI.m to summarize bootstrap distributions (median + CI).
%   - If you pass [nSubj x nBoot], the function will auto-permute.
% =========================================================================

%% ---- Figure settings ----
wd_ref    = 3;   % reference line width
wd_border = 4;   % axis/marker line width
fsz_ticks = 45;
sz_marker = 40;
CI_level = .90; % CI range for stats only! Default for plotting is .68; use 90% because of oone-tailed correlation

%% ---- Validate / canonicalize shapes to [nBoot x nSubj] ----
assert(ndims(X_allBoot_allSubj) <= 2 && ndims(Y_allBoot_allSubj) <= 2, ...
    'ALERT: X/Y asym inputs must be 2D: [nBoot x nSubj] (or [nSubj x nBoot]).');

% Determine nSubj from markers
nSubj = numel(markers_allSubj);

% Reshape to ensure both X and Y are [nBoot x nSubj]
if size(X_allBoot_allSubj, 2) ~= nSubj
    X_allBoot_allSubj = X_allBoot_allSubj';
end
if size(Y_allBoot_allSubj, 2) ~= nSubj
    Y_allBoot_allSubj = Y_allBoot_allSubj';
end

[nBoot, nSubj2] = size(X_allBoot_allSubj);
assert(nSubj2 == nSubj, 'ALERT: Cannot reconcile nSubj from markers with X/Y size.');

%% ---- Bootstrap summary for plotting: per-subject median and CI ----
% Outputs are 1 x nSubj by default; squeeze/transpose to nSubj x 1
[X_med, X_lb, X_ub] = getCI(X_allBoot_allSubj, 1, 1);
[Y_med, Y_lb, Y_ub] = getCI(Y_allBoot_allSubj, 1, 1);

% X_med = X_med(:); X_lb = X_lb(:); X_ub = X_ub(:);
% Y_med = Y_med(:); Y_lb = Y_lb(:); Y_ub = Y_ub(:);

%% ---- Plot ----
figure('Position', [0 200 1e3 1e3]); hold on; box on

for iSubj = 1:nSubj

    if flag_plotUnikSymbol == 1
        mk = markers_allSubj{iSubj};
    else
        mk = 'o';
    end

    % Marker face handling: symbols without face (e.g., '+','x')
    if ismember(mk, {'+','x','*','.'})
        faceColor = 'none';
        edgeColor = 'k';
    else
        faceColor = 'w';
        edgeColor = 'k';
    end

    % Optional per-point CI: draw horizontal (X) and vertical (Y) bars
    if flag_plotIdvdCI == 1
        plot([X_lb(iSubj), X_ub(iSubj)], [Y_med(iSubj), Y_med(iSubj)], '-', ...
            'Color', edgeColor, 'LineWidth', 2, 'HandleVisibility', 'off');
        plot([X_med(iSubj), X_med(iSubj)], [Y_lb(iSubj), Y_ub(iSubj)], '-', ...
            'Color', edgeColor, 'LineWidth', 2, 'HandleVisibility', 'off');
    end
    plot(X_med(iSubj), Y_med(iSubj), mk, ...
        'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor, ...
        'MarkerSize', sz_marker, 'LineWidth', wd_border);
end % iSubj

%% Linear regression on median points (reference line)
lm_med = polyfit(X_med, Y_med, 1);
x_lm = linspace(min(X_med), max(X_med), 200);
yfit_OnMed = polyval(lm_med, x_lm);

%% Linear regression & correlation analysis PER BOOT
yfit_allBoot = nan(nBoot, numel(x_lm));
eta2_allBoot = nan(nBoot, 1);

r_allBoot = nan(nBoot, 1);
p_r_allBoot = r_allBoot;
rho_allBoot = r_allBoot;
p_rho_allBoot = r_allBoot;
tau_allBoot = r_allBoot;
p_tau_allBoot = r_allBoot;

for iBoot = 1:nBoot
    x_perBoot = X_allBoot_allSubj(iBoot, :).';
    y_perBoot = Y_allBoot_allSubj(iBoot, :).';

    % Obtain linear regression fits
    if numel(unique(x_perBoot)) >= 2
        polyfit_perBoot = polyfit(x_perBoot, y_perBoot, 1);
        yfit_allBoot(iBoot, :) = polyval(polyfit_perBoot, x_lm);
        eta2_allBoot(iBoot) = var(polyval(polyfit_perBoot, x_perBoot)) / var(y_perBoot);
    end

    % Conduct correlation analysis of different types
    [r_allBoot(iBoot), p_r_allBoot(iBoot)] = corr(x_perBoot, y_perBoot, 'Type', 'Pearson', 'tail', 'both'); % left: assume slope <0
    [rho_allBoot(iBoot), p_rho_allBoot(iBoot)] = corr(x_perBoot, y_perBoot, 'Type', 'spearman', 'tail', 'both'); % left: assume slope <0
    [tau_allBoot(iBoot), p_tau_allBoot(iBoot)] = corr(x_perBoot, y_perBoot, 'Type', 'Kendall', 'tail', 'both'); % left: assume slope <0
end % iBoot

%% Obtain medians and CI of stats
% Linear regression
[~, yfit_lb, yfit_ub] = getCI(yfit_allBoot, 1, 1);
[eta2_med, eta2_lb, eta2_ub] = getCI(eta2_allBoot, 1, 1, CI_level);

% Correlation analysis
[r_med, r_lb, r_ub] = getCI(r_allBoot, 1, 1, CI_level);
% [p_r_med, p_r_lb, p_r_ub] = getCI(p_r_allBoot, 1, 1);
[rho_med, rho_lb, rho_ub] = getCI(rho_allBoot, 1, 1, CI_level);
% [p_rho_med, p_rho_lb, p_rho_ub] = getCI(p_rho_allBoot, 1, 1);
[tau_med, tau_lb, tau_ub] = getCI(tau_allBoot, 1, 1, CI_level);
% [p_tau_med, p_tau_lb, p_tau_ub] = getCI(p_tau_allBoot, 1, 1);

% Significance annotation based on CI excluding zero
str_sig_pearson = '';
str_sig_spearman = '';
str_sig_kendall = '';

if r_lb*r_ub > 0
    str_sig_pearson = '*';
end

if rho_lb*rho_ub > 0
    str_sig_spearman = '*';
end

if tau_lb*tau_ub > 0
    str_sig_kendall = '*';  
end

%% Plot regression fits
patch([x_lm fliplr(x_lm)], [yfit_lb fliplr(yfit_ub)], ones(1,3)*.6, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.20, 'HandleVisibility', 'off');
plot(x_lm, yfit_OnMed, '-', 'Color', ones(1,3)*.4, ...
    'HandleVisibility', 'off', 'LineWidth', wd_border*1.5);

%% ---- Reference lines (zero asymmetry) ----
xline(0, 'Color', ones(1,3)/2, 'LineWidth', wd_ref);
yline(0, 'Color', ones(1,3)/2, 'LineWidth', wd_ref);

%% ---- Axes ticks/limits ----
if ~isnan(x_ticks), xticks(x_ticks); xlim(x_ticks([1 end])); end
if ~isnan(y_ticks), yticks(y_ticks); ylim(y_ticks([1 end])); end
if ~isnan(x_ticklabels), xticklabels(x_ticklabels); end
if ~isnan(y_ticklabels), yticklabels(y_ticklabels); end

axis square
ax = gca;
ax.XAxis.FontSize = fsz_ticks;
ax.YAxis.FontSize = fsz_ticks;
ax.LineWidth = wd_border;

%% ---- Title / annotation ----
str_corr = sprintf(' Pearson%s: r=%.2f [%.2f, %.2f]\n Spearman%s: rho=%.2f [%.2f, %.2f]\n Kendall%s: tau=%.2f [%.2f, %.2f]', ...
    str_sig_pearson, r_med, r_lb, r_ub, ...
    str_sig_spearman, rho_med, rho_lb, rho_ub, ...
    str_sig_kendall, tau_med, tau_lb, tau_ub);

str_eta2 = sprintf('eta2=%.2f [%.2f, %.2f]', eta2_med, eta2_lb, eta2_ub);

str_title = sprintf('%s | %s\n%s', str_title, str_eta2, str_corr);
title(str_title);

end
