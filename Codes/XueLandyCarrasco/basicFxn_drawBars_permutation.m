function basicFxn_compAsym_permutation(asymX_allIter_allSubj, asymY_allIter_allSubj, nBins, y_ticks, sz_fig, str_title)
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
[nSubj, nIter] = size(asymX_allIter_allSubj);

colors = repmat(linspace(0, .5, nBins)', 1, 3);

fsz_ticks = 15;
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

%% Preallocate per-iteration stats on RAW Y grouped by bin label
Fvalue_allIter = nan(nIter, 1);
eta2_allIter   = nan(nIter, 1);

tscore_allIter = nan(nIter, 1);
CohenD_allIter = nan(nIter, 1);

%% 1) Per-iteration binning + stats
for iIter = 1:nIter
    x = asymX_allIter_allSubj(:, iIter);
    y = asymY_allIter_allSubj(:, iIter);

    % Bin subjects by X
    [~, ~, binIdx] = histcounts(x, nBins);

    % --- binned means for plotting ---
    for iBin = 1:nBins
        idx = (binIdx == iBin);
        N_bin_allIter(iBin, iIter) = sum(idx);

        if any(idx)
            Y_meanBin_allIter(iBin, iIter) = mean(y(idx), 'omitnan');
            X_meanBin_allIter(iBin, iIter) = mean(x(idx), 'omitnan');
        end
    end

    % --- stats on RAW y by bin groups (independent groups!) ---
    valid = (binIdx > 0) & ~isnan(y);
    if numel(unique(binIdx(valid))) >= 2
        [Fvalue_allIter(iIter), eta2_allIter(iIter)] = oneway_anova_indep(y(valid), binIdx(valid));
    end

    % --- if 2 bins: 2-sample t-test + Cohen's d (via fxn_getES) ---
    if nBins == 2
        idx1 = (binIdx == 1) & ~isnan(y);
        idx2 = (binIdx == 2) & ~isnan(y);
        if any(idx1) && any(idx2)
            [~, ~, ~, stats] = ttest2(y(idx1), y(idx2)); % independent samples
            tscore_allIter(iIter) = stats.tstat;

            % Use your existing effect size helper
            CohenD_allIter(iIter) = fxn_getES(y(idx1), y(idx2));
        end
    end
end

%% 2) Permutation p-value (shuffle bin labels within each iteration)
F_obs = median(Fvalue_allIter, 'omitnan');
F_med_allPerm = nan(nPerm, 1);

if nBins == 2
    t_obs = median(tscore_allIter, 'omitnan');
    t_med_allPerm = nan(nPerm, 1);
end

for iPerm = 1:nPerm
    F_perm_iter = nan(nIter, 1);

    if nBins == 2
        t_perm_iter = nan(nIter, 1);
    end

    for iIter = 1:nIter
        x = asymX_allIter_allSubj(:, iIter);
        y = asymY_allIter_allSubj(:, iIter);

        [~, ~, binIdx] = histcounts(x, nBins);

        valid = (binIdx > 0) & ~isnan(y);
        if numel(unique(binIdx(valid))) < 2
            continue
        end

        % Shuffle bin labels among valid subjects (break X->bin association)
        bin_shuf = binIdx(valid);
        bin_shuf = bin_shuf(randperm(numel(bin_shuf)));

        % ANOVA
        [F_perm_iter(iIter), ~] = oneway_anova_indep(y(valid), bin_shuf);

        % 2-sample t-test if 2 bins
        if nBins == 2
            yv = y(valid);
            idx1 = (bin_shuf == 1);
            idx2 = (bin_shuf == 2);
            if any(idx1) && any(idx2)
                [~, ~, ~, stats] = ttest2(yv(idx1), yv(idx2));
                t_perm_iter(iIter) = stats.tstat;
            end
        end
    end

    F_med_allPerm(iPerm) = median(F_perm_iter, 'omitnan');

    if nBins == 2
        t_med_allPerm(iPerm) = median(t_perm_iter, 'omitnan');
    end
end

p_ANOVA_perm = (1 + sum(F_med_allPerm >= F_obs)) / (numel(F_med_allPerm) + 1);

if nBins == 2
    p_ttest_perm = (1 + sum(abs(t_med_allPerm) >= abs(t_obs))) / (numel(t_med_allPerm) + 1);
end

%% 3) Median + CI for stats across iterations
[F_med, F_lb, F_ub]       = getCI(Fvalue_allIter, 1, 1, CI_level_stats);
[eta_med, eta_lb, eta_ub] = getCI(eta2_allIter,   1, 1, CI_level_stats);

str_ANOVA = sprintf('ANOVA (across bins): F=%.2f [%.2f, %.2f], p=%.3f, eta2=%.2f [%.2f, %.2f]', ...
    F_med, F_lb, F_ub, p_ANOVA_perm, eta_med, eta_lb, eta_ub);

if nBins == 2
    [t_med, t_lb, t_ub] = getCI(tscore_allIter, 1, 1, CI_level_stats);
    [d_med, d_lb, d_ub] = getCI(CohenD_allIter, 1, 1, CI_level_stats);

    str_ttest = sprintf('2-sample t: t=%.2f [%.2f, %.2f], p=%.3f, CohenD=%.2f [%.2f, %.2f]', ...
        t_med, t_lb, t_ub, p_ttest_perm, d_med, d_lb, d_ub);
else
    str_ttest = '';
end

%% 4) Plot: bars show median binned mean(Y) across iterations + CI for plotting
[Y_med_plot, Y_lb_plot, Y_ub_plot] = getCI(Y_meanBin_allIter', 1, 1, CI_level_plot);

figure('Position', [0, 200, sz_fig])
hold on

for iBin = 1:nBins
    bar(iBin, Y_med_plot(iBin), ...
        'FaceColor', colors(iBin,:), 'EdgeColor', colors(iBin,:), ...
        'BarWidth', wd_bar, 'HandleVisibility', 'off');

    errorbar(iBin, Y_med_plot(iBin), ...
        Y_med_plot(iBin)-Y_lb_plot(iBin), Y_ub_plot(iBin)-Y_med_plot(iBin), ...
        '.', 'color', colors(iBin,:), 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off');

    % Fancy white overlay
    if Y_med_plot(iBin) > 0
        errorbar(iBin, Y_med_plot(iBin), Y_med_plot(iBin)-Y_lb_plot(iBin), 0, ...
            '.w', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off');
    else
        errorbar(iBin, Y_med_plot(iBin), 0, Y_ub_plot(iBin)-Y_med_plot(iBin), ...
            '.w', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off');
    end
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

ax = gca;
ax.XAxis.FontSize = fsz_ticks;
ax.YAxis.FontSize = fsz_ticks;
ax.LineWidth = wd;

% title
if nBins == 2
    title(sprintf('%s\n%s\n%s', str_title, str_ANOVA, str_ttest), 'fontsize', fsz_title)
else
    title(sprintf('%s\n%s', str_title, str_ANOVA), 'fontsize', fsz_title)
end

end

%% ===== helper: one-way ANOVA for independent groups (bin label) =====
function [Fvalue, eta2] = oneway_anova_indep(y, g)
% Independent-groups one-way ANOVA
% y: Nx1 values, g: Nx1 integer group labels

y = y(:);
g = g(:);

valid = ~isnan(y) & ~isnan(g);
y = y(valid);
g = g(valid);

groups = unique(g);
k = numel(groups);
N = numel(y);

if k < 2 || N <= k
    Fvalue = NaN;
    eta2 = NaN;
    return
end

grand = mean(y, 'omitnan');

SS_between = 0;
SS_within  = 0;

for i = 1:k
    idx = (g == groups(i));
    yi  = y(idx);
    ni  = numel(yi);
    if ni == 0, continue; end

    mi = mean(yi, 'omitnan');
    SS_between = SS_between + ni * (mi - grand).^2;
    SS_within  = SS_within  + sum((yi - mi).^2, 'omitnan');
end

df_between = k - 1;
df_within  = N - k;

MS_between = SS_between / df_between;
MS_within  = SS_within  / df_within;

Fvalue = MS_between / MS_within;

SS_total = SS_between + SS_within;
eta2 = SS_between / SS_total;  % eta-squared
end
