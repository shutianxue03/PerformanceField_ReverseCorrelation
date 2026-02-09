function basicFxn_compAsym_permutation(asymX_allIter_allSubj, asymY_allIter_allSubj, nBins, y_ticks, sz_fig, str_title, str_ylabel)
% basicFxn_compAsym_permutation
%
% asymX_allIter_allSubj: delta contrast [nSubj x nIter]
% asymY_allIter_allSubj: delta parameter [nSubj x nIter]
% nBins: number of bins (based on X) used to group subjects on each iteration
%
% What it does:
% 1) For each iteration:
%    - bin subjects by X
%    - compute mean(Y) within each bin (for plotting)
%    - run independent-groups one-way ANOVA on RAW Y with group=bin label
%    - if nBins==2, run 2-sample t-test and compute Cohen's d via fxn_getES()
% 2) Summarize stats across iterations (median + CI_level_stats CI)
% 3) Permutation test for p-values: shuffle bin labels within each iteration, recompute stats,
%    compare to observed (median across iterations)
% 4) Plot: bar chart of median binned mean(Y) across iterations + CI_level_plot CI

%% Settings
[~, nIter] = size(asymX_allIter_allSubj);

colors = repmat(linspace(0, .5, nBins)', 1, 3);
fsz_title = 10;
wd = 2;
wd_bar = .5;
wd_ref = wd;
ref = 0;

CI_level_stats = 0.95;   % for stats strings
CI_level_plot  = 0.68;   % for plotting

nPerm = 1e4;

%% Preallocate binned means (for plotting)
Y_meanBin_allIter = nan(nBins, nIter);
X_meanBin_allIter = nan(nBins, nIter);   %#ok<NASGU> % optional
N_bin_allIter     = zeros(nBins, nIter); %#ok<NASGU> % optional

% Statistical analysis
tscore_allIter = nan(nIter, 1);
CohenD_allIter = nan(nIter, 1);
pT_allIter     = nan(nIter, 1);

%% (1) Per-iteration stats: independent-groups ANOVA across bins; and (if nBins==2) ttest2 + Cohen's d
% Assumes inside your main function you already have:
%   [nSubj, nIter] = size(asymX_allIter_allSubj);
%   nBins, etc.
%
% Outputs (preallocate before this block):
%   Fvalue_allIter, eta2_allIter, pA_allIter
%   tscore_allIter, CohenD_allIter, pT_allIter  (only meaningful if nBins==2)

for iIter = 1:nIter
    asymX_med_allSubj = asymX_allIter_allSubj(:, iIter);
    asymY_med_allSubj = asymY_allIter_allSubj(:, iIter);

    [~, ~, binIdx] = histcounts(asymX_med_allSubj, nBins);

    % validXY = ~isnan(x) & ~isnan(y);
    validXY = (binIdx > 0) & ~isnan(asymY_med_allSubj);
    if sum(validXY) < nBins
        continue
    end

    % Quantile edges -> roughly equal N per bin
    edges = quantile(asymX_med_allSubj(validXY), linspace(0, 1, nBins+1));
    edges(1)   = edges(1) - eps;      % include min
    edges(end) = edges(end) + eps;    % include max

    % If ties make repeated edges, discretize can fail; fall back to histcounts
    if numel(unique(edges)) < numel(edges)
        [~, ~, binIdx] = histcounts(asymX_med_allSubj, nBins);
    else
        binIdx = discretize(asymX_med_allSubj, edges);
    end

    for iBin = 1:nBins
        idx = (binIdx == iBin) & validXY;
        N_bin_allIter(iBin, iIter) = sum(idx);

        if any(idx)
            Y_meanBin_allIter(iBin, iIter) = mean(asymY_med_allSubj(idx), 'omitnan');
            X_meanBin_allIter(iBin, iIter) = mean(asymX_med_allSubj(idx), 'omitnan');
        end
    end

    % --- one-way ANOVA (independent groups) ---
    % if numel(unique(binIdx(validXY))) >= 2
    %     [pA, tbl] = anova1(y(validXY), binIdx(validXY), 'off');
    %     pA_allIter(iIter)     = pA;
    %     Fvalue_allIter(iIter) = tbl{2,5};
    %
    %     % eta^2 for independent-groups ANOVA: SS_between / SS_total
    %     SSb = tbl{2,2};
    %     SSt = tbl{4,2};
    %     eta2_allIter(iIter) = SSb / SSt;
    % end

    % --- if 2 bins: 2-sample t-test (independent) + Cohen's d ---
    if nBins == 2
        y1 = asymY_med_allSubj((binIdx == 1) & ~isnan(asymY_med_allSubj));
        y2 = asymY_med_allSubj((binIdx == 2) & ~isnan(asymY_med_allSubj));

        % Conduct two-sample t-test
        if ~isempty(y1) && ~isempty(y2)
            [~, pT, ~, stats] = ttest2(y1, y2);  % or: 'Vartype','unequal'
            pT_allIter(iIter)     = pT;
            tscore_allIter(iIter) = stats.tstat;

            % use your helper
            CohenD_allIter(iIter) = fxn_getES(y1, y2, 'independent');
        end
    end
end % iIter

%% Calculate bin difference
if nBins==2
    yDiff_allIter = Y_meanBin_allIter(1, :) - Y_meanBin_allIter(2, :);

    % Obtain 95% CI for bin difference across iterations
    [yDiff_med, yDiff_lb95, yDiff_ub95] = getCI(yDiff_allIter(:), 1, 1, CI_level_stats);
    [yDiff_med, yDiff_lb68, yDiff_ub68] = getCI(yDiff_allIter(:), 1, 1, .68);
    str_diff = sprintf('Bin1 - Bin2 = %.2f [%.2f, %.2f]\n', yDiff_med, yDiff_lb95, yDiff_ub95);
end

%% (2) Permutation p-values on MEDIAN x/y across iterations (no per-iter permutation)

% Median across iterations per subject (nSubj x 1)
x_med = getCI(asymX_allIter_allSubj, 1, 2);  % median over dim-2
y_med = getCI(asymY_allIter_allSubj, 1, 2);

x_med = x_med(:);
y_med = y_med(:);

% --- bin subjects by x_med (quantile bins) ---
validXY = ~isnan(x_med) & ~isnan(y_med);
if sum(validXY) < nBins
    error('Not enough valid subjects to form %d bins.', nBins);
end

edges = quantile(x_med(validXY), linspace(0, 1, nBins+1));
edges(1)   = edges(1) - eps;
edges(end) = edges(end) + eps;

if numel(unique(edges)) < numel(edges)
    [~, ~, binIdx] = histcounts(x_med, nBins);
else
    binIdx = discretize(x_med, edges);
end

valid = validXY & (binIdx > 0);
if numel(unique(binIdx(valid))) < 2
    error('Binning failed: fewer than 2 non-empty bins.');
end

% --- observed statistic on median y ---
if nBins ~= 2
    error('This block currently implemented for nBins==2 only (ttest2).');
end

y1 = y_med(valid & binIdx == 1);
y2 = y_med(valid & binIdx == 2);
if isempty(y1) || isempty(y2)
    error('Empty bin after binning (nBins=2).');
end

% Signed observed t-stat (important for left vs right)
[~, ~, ~, stats_obs] = ttest2(y1, y2);
t_obs_signed = stats_obs.tstat;

% --- permutation: shuffle y across subjects (bins fixed) ---
t_perm_signed = nan(nPerm, 1);

parfor iPerm = 1:nPerm
    y_perm = y_med;  %#ok<PFBNS> copy

    % permute y among VALID subjects only (FIXED)
    tmp = y_perm(valid);
    tmp = tmp(randperm(numel(tmp)));
    y_perm(valid) = tmp;

    yp1 = y_perm(valid & binIdx == 1);
    yp2 = y_perm(valid & binIdx == 2);

    if ~isempty(yp1) && ~isempty(yp2)
        [~, ~, ~, stp] = ttest2(yp1, yp2);
        t_perm_signed(iPerm) = stp.tstat;
    end
end

% Drop failed perms (if any)
t_perm_signed = t_perm_signed(~isnan(t_perm_signed));
nEff = numel(t_perm_signed);

% --- permutation p-values (L / 2 / R) with +1 correction ---
p_left  = (1 + sum(t_perm_signed <= t_obs_signed)) / (nEff + 1);
p_right = (1 + sum(t_perm_signed >= t_obs_signed)) / (nEff + 1);
p_two   = (1 + sum(abs(t_perm_signed) >= abs(t_obs_signed))) / (nEff + 1);

% for your title formatting, mimic the correlation script
p_ttest_all = [p_left, p_two, p_right];     % order = L / 2 / R
tail_lab = {'L','2','R'};
pstr_ttest = sprintf('p(%s/%s/%s)=[%.3f, %.3f, %.3f]', ...
    tail_lab{1}, tail_lab{2}, tail_lab{3}, p_ttest_all(1), p_ttest_all(2), p_ttest_all(3));

%% 3) Median + CI for stats across iterations
% [F_med, F_lb, F_ub]       = getCI(Fvalue_allIter, 1, 1, CI_level_stats);
% [eta_med, eta_lb, eta_ub] = getCI(eta2_allIter,   1, 1, CI_level_stats);

% str_ANOVA = sprintf('ANOVA: F=%.2f [%.2f, %.2f], p=%.3f, eta2=%.2f [%.2f, %.2f]', ...
%     F_med, F_lb, F_ub, p_ANOVA_perm, eta_med, eta_lb, eta_ub);

if nBins == 2
    [t_med, t_lb, t_ub] = getCI(tscore_allIter, 1, 1, CI_level_stats);
    [d_med, d_lb, d_ub] = getCI(CohenD_allIter, 1, 1, CI_level_stats);

    str_ttest = sprintf('2-sample t: t=%.2f [%.2f, %.2f], %s, CohenD=%.2f [%.2f, %.2f]', ...
        t_med, t_lb, t_ub, pstr_ttest, d_med, d_lb, d_ub);
else
    str_ttest = '';
end

%% 4) Plot: bars show median binned mean(Y) across iterations + CI for plotting
[Y_med_plot, Y_lb_plot, Y_ub_plot] = getCI(Y_meanBin_allIter', 1, 1, CI_level_plot);

figure('Position', [0, 200, sz_fig])
% PLOT 1
subplot(1,2,1), hold on

for iBin = 1:nBins
    bar(iBin, Y_med_plot(iBin), 'FaceColor', colors(iBin,:), 'EdgeColor', colors(iBin,:), 'BarWidth', wd_bar, 'HandleVisibility', 'off');
    errorbar(iBin, Y_med_plot(iBin), Y_med_plot(iBin)-Y_lb_plot(iBin), Y_ub_plot(iBin)-Y_med_plot(iBin), '.', 'color', 'k', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off');

    % % Fancy white overlay
    % if Y_med_plot(iBin) > 0
    %     errorbar(iBin, Y_med_plot(iBin), Y_med_plot(iBin)-Y_lb_plot(iBin), 0, ...
    %         '.w', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off');
    % else
    %     errorbar(iBin, Y_med_plot(iBin), 0, Y_ub_plot(iBin)-Y_med_plot(iBin), ...
    %         '.w', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off');
    % end
end

% reference line
if ~isnan(ref)
    yline(ref, 'color', ones(1,3)/2, 'linewidth', wd_ref, 'HandleVisibility', 'off');
end

% ticks / limits
if ~isnan(y_ticks)
    yticks(y_ticks)
    ylim(y_ticks([1 end]))
end

xticks(1:nBins)

buffer = 0.6;
xlim([1-buffer, nBins+buffer])

ylabel(str_ylabel)

% PLOT 2: distribution of difference
subplot(1,2,2), hold on
histogram(yDiff_allIter, 'Normalization', 'probability', 'FaceColor', ones(1,3)*0.7, 'EdgeColor', 'k');
xline(yDiff_med, 'r-', 'LineWidth', wd, 'Label', 'Median', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle', 'HandleVisibility', 'off');
xline(yDiff_lb95, 'r--', 'LineWidth', wd, 'Label', '95% CI', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle', 'HandleVisibility', 'off');
xline(yDiff_ub95, 'r--', 'LineWidth', wd, 'HandleVisibility', 'off');
xline(yDiff_lb68, 'b--', 'LineWidth', wd, 'Label', '95% CI', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle', 'HandleVisibility', 'off');
xline(yDiff_ub68, 'b--', 'LineWidth', wd, 'HandleVisibility', 'off');
xline(0, 'k-', 'LineWidth', wd*2, 'Label', 'No difference', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle', 'HandleVisibility', 'off');
xlabel('Bin1 - Bin2');
ylabel('Proportion of iterations');

% ax = gca;
% ax.XAxis.FontSize = fsz_ticks;
% ax.YAxis.FontSize = fsz_ticks;
% ax.LineWidth = wd;

% title
% if nBins == 2
%     title(sprintf('%s\n%s\n%s', str_title, str_ANOVA, str_ttest), 'fontsize', fsz_title)
% else
%     title(sprintf('%s\n%s', str_title, str_ANOVA), 'fontsize', fsz_title)
% end

if nBins == 2
    sgtitle(sprintf('%s\n%s\n%s', str_title, str_ttest, str_diff), 'fontsize', fsz_title)
else
    sgtitle(str_title, 'fontsize', fsz_title)
end

end