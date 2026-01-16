function basicFxn_drawCorrAsym_permutation(X_allIter_allSubj, Y_allIter_allSubj, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)

% =========================================================================
% basicFxn_drawCorrAsym_permutation.m
%
% Purpose
%   Plot and compare *asymmetry* metrics between two measurements (X vs Y) across observers
%
% Inputs (canonical / recommended)
%   asymX_allIter_allSubj : [nIter x nSubj] 
%   asymY_allIter_allSubj : [nIter x nSubj] 
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
%   str_sig   : '_sig' / '_mg' / '_ns' based on median p-value across iterations
%   str_title : full title string including r/p medians and CIs
%   stats     : struct with medians, CIs, and distributions
%
% =========================================================================

%% ---- Figure settings ----
wd_ref    = 3;   % reference line width
wd_border = 4;   % axis/marker line width
fsz_ticks = 45;
sz_marker = 40;
CI_level = .90; % CI range for stats only! Default for plotting is .68; use 90% because of oone-tailed correlation

%% ---- Validate / canonicalize shapes to [nIter x nSubj] ----
assert(ndims(X_allIter_allSubj) <= 2 && ndims(Y_allIter_allSubj) <= 2, ...
    'ALERT: X/Y asym inputs must be 2D: [nIter x nSubj] (or [nSubj x nIter]).');

% Determine nSubj from markers
nSubj = numel(markers_allSubj);

% Reshape to ensure both X and Y are [nIter x nSubj]
if size(X_allIter_allSubj, 2) ~= nSubj
    X_allIter_allSubj = X_allIter_allSubj';
end
if size(Y_allIter_allSubj, 2) ~= nSubj
    Y_allIter_allSubj = Y_allIter_allSubj';
end

[nIter, nSubj2] = size(X_allIter_allSubj);
assert(nSubj2 == nSubj, 'ALERT: Cannot reconcile nSubj from markers with X/Y size.');

%% Obtain median and CI
% Outputs are 1 x nSubj by default; squeeze/transpose to nSubj x 1
[X_med, X_lb, X_ub] = getCI(X_allIter_allSubj, 1, 1);
[Y_med, Y_lb, Y_ub] = getCI(Y_allIter_allSubj, 1, 1);

%% Permutation test on MEDIAN asymmetry points (p-values)
nPerm = 1e4;

x_obs = X_med(:);
y_obs = Y_med(:);

% observed correlations
r_obs_pearson  = corr(x_obs, y_obs, 'Type', 'Pearson',  'Rows', 'complete');
r_obs_spearman = corr(x_obs, y_obs, 'Type', 'Spearman', 'Rows', 'complete');
r_obs_kendall  = corr(x_obs, y_obs, 'Type', 'Kendall',  'Rows', 'complete');

% null distributions (permute observer labels of y; keeps x fixed)
r_perm_pearson  = nan(nPerm,1);
r_perm_spearman = nan(nPerm,1);
r_perm_kendall  = nan(nPerm,1);

for iPerm = 1:nPerm
    idx = randperm(nSubj);
    yp  = y_obs(idx);

    r_perm_pearson(iPerm)  = corr(x_obs, yp, 'Type', 'Pearson',  'Rows', 'complete');
    r_perm_spearman(iPerm) = corr(x_obs, yp, 'Type', 'Spearman', 'Rows', 'complete');
    r_perm_kendall(iPerm)  = corr(x_obs, yp, 'Type', 'Kendall',  'Rows', 'complete');
end

% ---- permutation p-values for ALL tails ----
tail_all = {'left','both','right'};     % L / two-sided / R
tail_lab = {'L','2','R'};

p_perm_pearson_all  = nan(1,3);
p_perm_spearman_all = nan(1,3);
p_perm_kendall_all  = nan(1,3);

for iTail = 1:3
    p_perm_pearson_all(iTail)  = fxn_perm_pval(r_perm_pearson,  r_obs_pearson,  tail_all{iTail});
    p_perm_spearman_all(iTail) = fxn_perm_pval(r_perm_spearman, r_obs_spearman, tail_all{iTail});
    p_perm_kendall_all(iTail)  = fxn_perm_pval(r_perm_kendall,  r_obs_kendall,  tail_all{iTail});
end


%% Plot individual data
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

%% Plot group averages
[X_ave, ~, ~, X_sem] = getCI(X_med(:), 2, 1);
[Y_ave, ~, ~, Y_sem] = getCI(Y_med(:), 2, 1);
errorbar(X_ave, Y_ave, X_sem, 'horizontal', 'k.', 'LineWidth', wd_border*2);
errorbar(X_ave, Y_ave, Y_sem, 'vertical', 'k.', 'LineWidth', wd_border*2);

%% Linear regression on median points (reference line)
lm_med = polyfit(X_med, Y_med, 1);
x_lm = linspace(min(X_med), max(X_med), 200);
yfit_OnMed = polyval(lm_med, x_lm);

%% Linear regression & correlation analysis PER Iteration
yfit_allIter = nan(nIter, numel(x_lm));

r_allIter = nan(nIter, 1);
p_r_allIter = r_allIter;
rho_allIter = r_allIter;
p_rho_allIter = r_allIter;
tau_allIter = r_allIter;
p_tau_allIter = r_allIter;

for iIter = 1:nIter
    x_perIter = X_allIter_allSubj(iIter, :).';
    y_perIter = Y_allIter_allSubj(iIter, :).';

    % Obtain linear regression fits
    if numel(unique(x_perIter)) >= 2
        polyfit_perIter = polyfit(x_perIter, y_perIter, 1);
        yfit_allIter(iIter, :) = polyval(polyfit_perIter, x_lm);
    end

    % Conduct correlation analysis of different types
    [r_allIter(iIter), p_r_allIter(iIter)] = corr(x_perIter, y_perIter, 'Type', 'Pearson', 'tail', 'both'); % left: assume slope <0
    [rho_allIter(iIter), p_rho_allIter(iIter)] = corr(x_perIter, y_perIter, 'Type', 'spearman', 'tail', 'both'); % left: assume slope <0
    [tau_allIter(iIter), p_tau_allIter(iIter)] = corr(x_perIter, y_perIter, 'Type', 'Kendall', 'tail', 'both'); % left: assume slope <0
end % iIter

%% Obtain medians and CI of stats
% Linear regression
[~, yfit_lb, yfit_ub] = getCI(yfit_allIter, 1, 1);

% Correlation analysis
[r_med, r_lb, r_ub] = getCI(r_allIter, 1, 1, CI_level);
% [p_r_med, p_r_lb, p_r_ub] = getCI(p_r_allIter, 1, 1);
[rho_med, rho_lb, rho_ub] = getCI(rho_allIter, 1, 1, CI_level);
% [p_rho_med, p_rho_lb, p_rho_ub] = getCI(p_rho_allIter, 1, 1);
[tau_med, tau_lb, tau_ub] = getCI(tau_allIter, 1, 1, CI_level);
% [p_tau_med, p_tau_lb, p_tau_ub] = getCI(p_tau_allIter, 1, 1);

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
% optional: make the p-value triplet compact in the title
pstr_pearson  = sprintf('p(%s/%s/%s)=[%.3f, %.3f , %.3f]', tail_lab{:}, p_perm_pearson_all);
pstr_spearman = sprintf('p(%s/%s/%s)=[%.3f, %.3f , %.3f]', tail_lab{:}, p_perm_spearman_all);
pstr_kendall  = sprintf('p(%s/%s/%s)=[%.3f, %.3f , %.3f]', tail_lab{:}, p_perm_kendall_all);

str_corr = sprintf([' Pearson%s: r=%.2f [%.2f, %.2f], %s\n' ...
    ' Spearman%s: rho=%.2f [%.2f, %.2f], %s\n' ...
    ' Kendall%s: tau=%.2f [%.2f, %.2f], %s'], ...
    str_sig_pearson,  r_med,   r_lb,   r_ub,   pstr_pearson, ...
    str_sig_spearman, rho_med, rho_lb, rho_ub, pstr_spearman, ...
    str_sig_kendall,  tau_med, tau_lb, tau_ub, pstr_kendall);

str_title = sprintf('%s\n%s\n', str_title, str_corr);
title(str_title);

end

function p = fxn_perm_pval(nullStats, obsStat, type_tail)
% Permutation p-value with +1 correction.
% type_tail: 'both' (two-sided), 'left', or 'right'

nullStats = nullStats(~isnan(nullStats));
n = numel(nullStats);

switch lower(type_tail)
    case {'both','two','two-sided','twosided'}
        p = (1 + sum(abs(nullStats) >= abs(obsStat))) / (n + 1);
    case {'right','greater'}
        p = (1 + sum(nullStats >= obsStat)) / (n + 1);
    case {'left','less'}
        p = (1 + sum(nullStats <= obsStat)) / (n + 1);
    otherwise
        error('Unknown type_tail: %s', type_tail);
end
end
