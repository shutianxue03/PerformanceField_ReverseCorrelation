function basicFxn_drawCorrAsym_permutation( ...
    X_allIter_allSubj, Y_allIter_allSubj, colors, flag_UseRUseRho, flag_CIrange, x_ticks, y_ticks, ...
    x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, ...
    str_title, markers_allSubj)

% =========================================================================
% basicFxn_drawCorrAsym_permutation (REVISED)
%
% Changes vs. your original:
%   - REMOVE all per-iteration correlation/regression calculations
%   - KEEP permutation test on MEDIANS (NHST-friendly p-values)
%   - ADD bootstrap over SUBJECTS for CIs of r/rho/tau + regression band
%   - REVISE string generation to report: point estimate (on medians) + bootstrap CI + permutation p
%
% Inputs
%   X_allIter_allSubj : [nIter x nSubj] or [nSubj x nIter]
%   Y_allIter_allSubj : [nIter x nSubj] or [nSubj x nIter]
% =========================================================================

%% ---- Figure settings ----
wd_ref  = 5;   % reference line width
wd_border = 4;   % axis/marker line width
fsz_ticks = 45;
sz_marker = 35;

CI68 = 0.68;     % per-point CI bars in plot (from iterations)
nPerm = 1e4;
nBoot = 1e4;
seedPerm = 1;
seedBoot = 2;

% Determine conducting one- or two-tailed corr analysis
if flag_CIrange==1, CI_level=.90; % one-tailed corr
else, CI_level=.95;% two-tailed corr 
end

%% ---- Validate / canonicalize shapes to [nIter x nSubj] ----
assert(ndims(X_allIter_allSubj) <= 2 && ndims(Y_allIter_allSubj) <= 2, ...
    'ALERT: X/Y asym inputs must be 2D: [nIter x nSubj] (or [nSubj x nIter]).');

nSubj = numel(markers_allSubj);

% Ensure [nIter x nSubj]
if size(X_allIter_allSubj, 2) ~= nSubj
    X_allIter_allSubj = X_allIter_allSubj';
end
if size(Y_allIter_allSubj, 2) ~= nSubj
    Y_allIter_allSubj = Y_allIter_allSubj';
end

[nIter, nSubj2] = size(X_allIter_allSubj);
assert(nSubj2 == nSubj, 'ALERT: Cannot reconcile nSubj from markers with X/Y size.');
assert(all(size(Y_allIter_allSubj) == [nIter, nSubj]), 'ALERT: X and Y must have the same size.');

%% ---- 1) Summaries across iterations for plotting points ----
% Median across iterations for each subject; CI half-lengths for per-point bars (68% by default)
[X_med, ~, ~, X_sem_neg68, X_sem_pos68] = getCI(X_allIter_allSubj, 1, 1, CI68);
[Y_med, ~, ~, Y_sem_neg68, Y_sem_pos68] = getCI(Y_allIter_allSubj, 1, 1, CI68);

X_med = X_med(:);
Y_med = Y_med(:);

X_sem_neg68 = X_sem_neg68(:);
X_sem_pos68 = X_sem_pos68(:);
Y_sem_neg68 = Y_sem_neg68(:);
Y_sem_pos68 = Y_sem_pos68(:);

x_obs = X_med;
y_obs = Y_med;

%% ---- 2) Observed correlations on MEDIANS ----
r_obs_pearson  = corr(x_obs, y_obs, 'Type', 'Pearson',  'Rows', 'complete');
r_obs_spearman = corr(x_obs, y_obs, 'Type', 'Spearman', 'Rows', 'complete');
r_obs_kendall  = corr(x_obs, y_obs, 'Type', 'Kendall',  'Rows', 'complete');

%% ---- 3) Permutation test on MEDIANS (permute observer labels of Y) ----
r_perm_pearson  = nan(nPerm,1);
r_perm_spearman = nan(nPerm,1);
r_perm_kendall  = nan(nPerm,1);

% Pregenerate subj indices
rng(seedPerm, 'twister');
indRand_allPerm = zeros(nPerm, nSubj, 'uint16');
for iPerm = 1:nPerm
    indRand_allPerm(iPerm, :) = uint16(randperm(nSubj));
end

parfor iPerm = 1:nPerm
    indRandPerm = double(indRand_allPerm(iPerm, :));
    yRand  = y_obs(indRandPerm); % only shuffle y, no need for x

    r_perm_pearson(iPerm)  = corr(x_obs, yRand, 'Type', 'Pearson',  'Rows', 'complete');
    r_perm_spearman(iPerm) = corr(x_obs, yRand, 'Type', 'Spearman', 'Rows', 'complete');
    r_perm_kendall(iPerm)  = corr(x_obs, yRand, 'Type', 'Kendall',  'Rows', 'complete');
end % iPerm

% permutation p-values for all tails (L / 2 / R)
tail_all = {'left','both','right'};
tail_lab = {'L','2','R'};
nTails = numel(tail_all);

p_perm_pearson_all  = nan(1,nTails);
p_perm_spearman_all = nan(1,nTails);
p_perm_kendall_all  = nan(1,nTails);

for iTail = 1:nTails
    p_perm_pearson_all(iTail)  = fxn_perm_pval(r_perm_pearson,  r_obs_pearson,  tail_all{iTail});
    p_perm_spearman_all(iTail) = fxn_perm_pval(r_perm_spearman, r_obs_spearman, tail_all{iTail});
    p_perm_kendall_all(iTail)  = fxn_perm_pval(r_perm_kendall,  r_obs_kendall,  tail_all{iTail});
end % iTail

%% ---- 4) Bootstrap over SUBJECTS for CIs of correlations + regression band ----
groupAveX_allBoot = nan(nBoot,1);
groupAveY_allBoot = nan(nBoot,1);
r_pearson_allBoot  = nan(nBoot,1);
r_spearman_allBoot = nan(nBoot,1);
r_kendall_allBoot  = nan(nBoot,1);

% Regression band: predicted y at x_lm for each bootstrap
xmin = min(x_obs); xmax = max(x_obs);
if xmin == xmax
    xmin = xmin - 1e-6; xmax = xmax + 1e-6;
end
x_lm = linspace(xmin, xmax, 200);
yfit_allBoot = nan(nBoot, numel(x_lm));

% Pregenerate subj indices
rng(seedBoot, 'twister');
indRand_allBoot = randi(nSubj, [nBoot, nSubj], 'uint16');

parfor iBoot = 1:nBoot
    indRandBoot = double(indRand_allBoot(iBoot, :));
    xRand = x_obs(indRandBoot);
    yRand = y_obs(indRandBoot);

    groupAveX_allBoot(iBoot) = mean(xRand, 'omitnan');
    groupAveY_allBoot(iBoot) = mean(yRand, 'omitnan');

    r_pearson_allBoot(iBoot)  = corr(xRand, yRand, 'Type', 'Pearson',  'Rows', 'complete');
    r_spearman_allBoot(iBoot) = corr(xRand, yRand, 'Type', 'Spearman', 'Rows', 'complete');
    r_kendall_allBoot(iBoot)  = corr(xRand, yRand, 'Type', 'Kendall',  'Rows', 'complete');

    indOK = ~isnan(xRand) & ~isnan(yRand);
    if sum(indOK) >= 2 && numel(unique(xRand(indOK))) >= 2
        p = polyfit(xRand(indOK), yRand(indOK), 1);
        yfit_allBoot(iBoot,:) = polyval(p, x_lm);
    end
end % iBoot

% Summarize point and interval estimates (bootstrap)
[groupAveX_med, ~, ~, groupAveX_sem_neg68, groupAveX_sem_pos68] = getCI(groupAveX_allBoot, 1, 1, CI68);
[groupAveY_med, ~, ~, groupAveY_sem_neg68, groupAveY_sem_pos68] = getCI(groupAveY_allBoot, 1, 1, CI68);

[r_med,   r_lb,   r_ub]   = getCI(r_pearson_allBoot,  1, 1, CI_level);
[rho_med, rho_lb, rho_ub] = getCI(r_spearman_allBoot, 1, 1, CI_level);
[tau_med, tau_lb, tau_ub] = getCI(r_kendall_allBoot,  1, 1, CI_level);

[~, yfit_lb68, yfit_ub68] = getCI(yfit_allBoot, 1, 1, CI68);

% Regression line from MEDIANS (point estimate)
indOK_med = ~isnan(x_obs) & ~isnan(y_obs);
lm_med = polyfit(x_obs(indOK_med), y_obs(indOK_med), 1);
yfit_OnMed = polyval(lm_med, x_lm);

%% Generate color gradient
% --- 1) Decide reference axis direction u ---
switch flag_UseRUseRho
    case 'useR'
        sgnCorr = sign(r_med);
    case 'useRho'
        sgnCorr = sign(rho_med);
end
u = [1, sgnCorr] / sqrt(2);

% --- 2) Project each point onto the reference axis ---
t = X_med * u(1) + Y_med * u(2);   % Nx1 signed coordinate along axis

% --- 3) FIX hue by side ---
% bottom-left (typically t<0) -> colors(1,:)
color0 = repmat(colors(1,:), numel(t), 1);   % default
color0(t < 0, :) = repmat(colors(2,:), sum(t < 0), 1);

% --- 4) Brightness only, normalized SEPARATELY for t<0 and t>0 ---
s = zeros(size(t));

neg = (t < 0);
pos = (t > 0);

tmaxNeg = max(abs(t(neg)));   % farthest on negative side
tmaxPos = max(abs(t(pos)));   % farthest on positive side

if any(neg) && tmaxNeg > 0
    s(neg) = abs(t(neg)) ./ tmaxNeg;   % 0 near origin -> 1 at neg extreme
end
if any(pos) && tmaxPos > 0
    s(pos) = abs(t(pos)) ./ tmaxPos;   % 0 near origin -> 1 at pos extreme
end

s = max(0, min(1, s));  % clamp
sMin=.25; % 
s=sMin+(1-sMin)*s;
% Mix with white: s=0 -> white-ish, s=1 -> pure hue
edgeColor_allSubj = (1 - s).*ones(size(color0)) + s.*color0;  % Nx3

%% ---- 5) Plot individual data (medians + optional per-point CI bars) ----
figure('Position', [0 200 1e3 1e3]); hold on; box on

for iSubj = 1:nSubj

    if flag_plotUnikSymbol == 1
        mk = markers_allSubj{iSubj};
    else
        mk = 'o';
    end

    if ismember(mk, {'+','x','*','.'})
        faceColor = 'none';
        edgeColor = 'k';
    else
        faceColor = 'w';
        edgeColor = 'k';
    end

    edgeColor = edgeColor_allSubj(iSubj, :);

    % Plot 68% CI
    if flag_plotIdvdCI == 1
        % horizontal (X) and vertical (Y) per-point CI bars (from iterations)
        errorbar(X_med(iSubj), Y_med(iSubj), Y_sem_neg68(iSubj), Y_sem_pos68(iSubj), '.', 'Color', edgeColor, 'LineWidth', wd_border/1.5, 'HandleVisibility', 'off', 'CapSize', 0);
        errorbar(X_med(iSubj), Y_med(iSubj), X_sem_neg68(iSubj), X_sem_pos68(iSubj), '.', 'horizontal', 'Color', edgeColor, 'LineWidth', wd_border/1.5, 'HandleVisibility', 'off', 'CapSize', 0);
    end

    % Plot data point
    % plot(X_med(iSubj), Y_med(iSubj), mk, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', sz_marker, 'LineWidth', wd_border);
    plot(X_med(iSubj), Y_med(iSubj), mk, 'MarkerFaceColor', edgeColor, 'MarkerEdgeColor', 'w', 'MarkerSize', sz_marker, 'LineWidth', wd_border);
end % iSubj

%% ---- 6) Plot group average (bootstrap over subjects on medians) ----
errorbar(groupAveX_med, groupAveY_med, groupAveX_sem_neg68, groupAveX_sem_pos68, 'horizontal', 'k.', 'LineWidth', wd_border*2, 'HandleVisibility','off', 'CapSize', 0);
errorbar(groupAveX_med, groupAveY_med, groupAveY_sem_neg68, groupAveY_sem_pos68, 'vertical',   'k.', 'LineWidth', wd_border*2, 'HandleVisibility','off', 'CapSize', 0);
plot(groupAveX_med, groupAveY_med, 's', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'linewidth', wd_border*1.5, 'markersize', sz_marker, 'HandleVisibility','off')

%% ---- 7) Plot regression fits: bootstrap CI band + median-fit line ----
patch([x_lm fliplr(x_lm)], [yfit_lb68 fliplr(yfit_ub68)], ones(1,3)*.6, 'EdgeColor', 'none', 'FaceAlpha', 0.20, 'HandleVisibility', 'off');
plot(x_lm, yfit_OnMed, '-', 'Color', ones(1,3)*.4, 'HandleVisibility', 'off', 'LineWidth', wd_border*1.5);

%% ---- 8) Reference lines (zero asymmetry) ----
xline(0, '--', 'Color', ones(1,3)/2, 'LineWidth', wd_ref);
yline(0, '--', 'Color', ones(1,3)/2, 'LineWidth', wd_ref);

%% ---- 9) Axes ticks/limits ----
if ~isnan(x_ticks), xticks(x_ticks); xlim(x_ticks([1 end])); end
if ~isnan(y_ticks), yticks(y_ticks); ylim(y_ticks([1 end])); end
if ~isnan(x_ticklabels), xticklabels(x_ticklabels); end
if ~isnan(y_ticklabels), yticklabels(y_ticklabels); end

axis square
xabs = gca;
xabs.XAxis.FontSize = fsz_ticks;
xabs.YAxis.FontSize = fsz_ticks;
xabs.LineWidth = wd_border;

%% ---- 10) Title / annotation ----
pstr_pearson  = sprintf('p(%s/%s/%s)=[%.3f, %.3f, %.3f]', tail_lab{:}, p_perm_pearson_all);
pstr_spearman = sprintf('p(%s/%s/%s)=[%.3f, %.3f, %.3f]', tail_lab{:}, p_perm_spearman_all);
pstr_kendall  = sprintf('p(%s/%s/%s)=[%.3f, %.3f, %.3f]', tail_lab{:}, p_perm_kendall_all);

sigMark = @(lb,ub) ternary((lb*ub>0), '*', '');

str_corr = sprintf([ ...
    'Pearson%s: r=%.2f [%.2f, %.2f], %s\n' ...
    'Spearman%s: \\rho=%.2f [%.2f, %.2f], %s\n' ...
    'Kendall%s: \\tau=%.2f [%.2f, %.2f], %s'], ...
    sigMark(r_lb,   r_ub),   r_med,  r_lb,   r_ub,   pstr_pearson, ...
    sigMark(rho_lb, rho_ub), rho_med, rho_lb, rho_ub, pstr_spearman, ...
    sigMark(tau_lb, tau_ub), tau_med,  tau_lb, tau_ub, pstr_kendall);

title(sprintf('%s %.0f%% CI \n%s\n', str_title, CI_level*100, str_corr));

%% Print string at the bottom-center
switch flag_UseRUseRho
    case 'useR'
        str_print = sprintf('r = %.2f , CI_{%.0f}=[%.2f, %.2f]', r_med, CI_level*100, r_lb, r_ub);
    case 'useRho'
        str_print = sprintf('\\rho = %.2f, CI_{%.0f}=[%.2f, %.2f]', rho_med, CI_level*100, rho_lb, rho_ub);
end

xabs = gca;
text(xabs, 0.5, 0.02, str_print, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 50, ...
    'Color', 'k', ...
    'Interpreter', 'tex', ...
    'Clipping', 'off');

end

%% ===== helper: permutation p-value with +1 correction =====
function p = fxn_perm_pval(nullStats, obsStat, type_tail)
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

%% ===== helper: simple ternary =====
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end