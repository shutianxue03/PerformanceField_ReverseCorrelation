function basicFxn_drawBars_permutation( ...
    data_allIter_allSubj, ref, colors, x_ticklabels, y_ticks, y_ticklabels, ...
    flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIter, markers_allSubj)

% basicFxn_drawBars_permutation
% - Omnibus RM one-way ANOVA (per iter -> CI; perm p via within-subj label shuffle on median)
% - Pairwise posthoc paired t-tests for ALL pairs (per iter -> CI; perm p via sign-flip on median diffs)
% - Bonferroni correction for pairwise perm p-values
% - Only print significant pairs (p_perm < pThresh_print) in the title

%% ---- settings ----
sz_marker_idvd = 10;
sz_ticks = 25; % behav and sep: 35; tunC: 25
sz_title = 10;
wd = 3; % behav: 3; sep: 2; tunC: 3
wd_bar = .5;

CI95 = .95;
nPerm = 1e4;
% pThresh_print = 0.1;   % print only pairs with perm p < this (after Bonferroni)
nBoot=1e3;
fprintf('\n\n *** %s ***\n', str_title);

%% ---- shape checks: force to [nIter x nSubj x nCond] ----
assert(ndims(data_allIter_allSubj) == 3, 'ALERT: data_allIter_allSubj must be 3D arrays.');

nCond = size(colors, 1);
nSubj = numel(markers_allSubj);

data_allIter_allSubj = fxn_reshape(data_allIter_allSubj, nSubj, nCond, 'data_allIter_allSubj');

[nIter2, nSubj2, nCond2] = size(data_allIter_allSubj);
assert(nIter2 == nIter, 'ALERT: Input nIter=%d but inferred nIter=%d from X.', nIter, nIter2);
assert(nSubj2 == nSubj && nCond2 == nCond, 'ALERT: Reshaped matrix has wrong nSubj/nCond.');

pairs  = nchoosek(1:nCond, 2);
nPairs = size(pairs, 1);

%% Obtain median and CI of the data
[data_med_allSubj, ~, ~] = getCI(data_allIter_allSubj, 1, 1, CI95); % [nSubj x nCond] median over iters

%% [NHST] Permutation p-values for ANOVA and paired contrasts (subject-level)

% (A) Omnibus RM one-way ANOVA: within-subject label shuffle
[Fvalue_obs, eta2p_obs]     = rm_oneway(data_med_allSubj);   % observed F on subject-level table
Fvalue_allPerm = nan(nPerm, 1);

parfor iPerm = 1:nPerm
    data_perPerm = data_med_allSubj;
    for iSubj = 1:nSubj
        data_perPerm(iSubj, :) = data_perPerm(iSubj, randperm(nCond));  % shuffle condition labels within subject
    end
    Fvalue_allPerm(iPerm) = rm_oneway(data_perPerm);
end

pperm_ANOVA = (1 + sum(Fvalue_allPerm >= Fvalue_obs)) / (nPerm + 1);

% (B) Planned pairwise contrasts: sign-flip permutation on within-subject differences
t_obs_allPairs    = nan(1, nPairs);   % observed |t| for each pair
d_obs_allPairs = t_obs_allPairs;
pperm_allPairs    = t_obs_allPairs;   % uncorrected permutation p-values

% sgnMat will be used to generate null distributions consistently across pairs.
sgnMat = (rand(nSubj, nPerm) > 0.5) * 2 - 1; % [nSubj x nPerm]

for iPair = 1:nPairs
    iCondA = pairs(iPair, 1);
    iCondB = pairs(iPair, 2);

    diffk = data_med_allSubj(:, iCondA) - data_med_allSubj(:, iCondB);  % [nSubj x 1]
    indOK = ~isnan(diffk);
    diffk = diffk(indOK);
    nk    = numel(diffk);

    if nk < 2 || std(diffk, 0) == 0
        t_obs_allPairs(iPair) = NaN;
        pperm_allPairs(iPair) = NaN;
        continue;
    end

    % Observed paired-t statistic magnitude (two-tailed)
    denom = std(diffk, 0) / sqrt(nk);
    t_obs = abs(mean(diffk) / denom);
    t_obs_allPairs(iPair) = t_obs;
    d_obs_allPairs(iPair) = mean(diffk)/std(diffk);

    % Sign-flip null: flip within-subject diffs by +/-1
    % Use the first nk rows of sgnMat restricted to valid subjects.
    sgn = sgnMat(indOK, :);                     % [nk x nPerm]
    num_perm = mean(diffk .* sgn, 1);           % [1 x nPerm]
    t_perm   = abs(num_perm / denom);           % [1 x nPerm]

    pperm_allPairs(iPair) = (1 + sum(t_perm >= t_obs)) / (nPerm + 1);
end % iPair

% Strings for NHST
str_ANOVA = sprintf('ANOVA: F=%.2f, p=%.3f, eta2p=%.2f', Fvalue_obs, pperm_ANOVA, eta2p_obs);
% string for t-tests will be shown with bootstraps

%% [Bootstrap] CIs for planned contrasts (row bootstrap over subjects)

CI_level_diff = 1 - 0.05;   % planned -> 95%
CI_level_ref  = 1 - 0.05;   % planned -> 95%

diffCond_allBoot = nan(nBoot, nPairs); % difference among conditions
diffRef_allBoot = nan(nBoot, nCond); % difference between each cond and the ref

for iBoot = 1:nBoot
    indResampled = randi(nSubj, [1, nSubj]);   % resample rows (subjects) with replacement

    dataRand = data_med_allSubj(indResampled, :);   % [nSubj x nCond] bootstrap sample

    % Pairwise planned differences
    for iPair = 1:nPairs
        iCondA = pairs(iPair, 1);
        iCondB = pairs(iPair, 2);
        diffCond_allBoot(iBoot, iPair) = mean(dataRand(:, iCondA) - dataRand(:, iCondB));
    end

    % Comparisons to reference
    if ~isnan(ref)
        for iCond = 1:nCond
            diffRef_allBoot(iBoot, iCond) = mean(dataRand(:, iCond) - ref);
        end
    end
end % iBoot

[diffCond_med, diff_lb, diff_ub] = getCI(diffCond_allBoot, 1, 1, CI_level_diff);
[diffRef_med, diffRef_lb, diffRef_ub] = getCI(diffRef_allBoot, 1, 1, CI_level_ref);

%% Bootstrap group averages
% We bootstrap subjects (rows) on data_med_allSubj (subject-level summaries).
CI_plot = 0.68;          % 68% for plotting (visual readability)
nBootPlot = 5000;

ave_allBoot = nan(nBootPlot, nCond);
for iBoot = 1:nBootPlot
    idx = randi(nSubj, [1 nSubj]);                 % resample subjects with replacement
    ave_allBoot(iBoot, :) = mean(data_med_allSubj(idx, :), 1);  % bootstrap mean per condition
end

% Obtain average and sem of group averages
[data_ave, ~, ~, data_sem_neg, data_sem_pos] = getCI(ave_allBoot, 1, 1, CI_plot);  % across bootstraps

%% [Plot] Group averages
figure('Position', [0, 200, sz_fig]); hold on

for iCond = 1:nCond
    % Bootstrapped group average
    bar(iCond, data_ave(iCond), 'FaceColor', colors(iCond,:), 'EdgeColor', colors(iCond,:), 'BarWidth', wd_bar, 'HandleVisibility', 'off');

    % Error bars for bootstrap CI
    errorbar(iCond, data_ave(iCond), data_sem_neg(iCond), data_sem_pos(iCond), '.', 'color', colors(iCond,:), 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off');

    % Fancy white overlay (pure aesthetics)
    if data_ave(iCond) > 0
        errorbar(iCond, data_ave(iCond), data_sem_neg(iCond), 0, '.w', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off');
    else
        errorbar(iCond, data_ave(iCond), 0, data_sem_pos(iCond), '.w', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off');
    end
end % iCond

%% [Plot] Idvd medians
if flag_plotIDVD
    if nCond == 2
        xIDVD = [1.3, 1.7];
    else
        xIDVD = 1:nCond;
    end

    for iSubj = 1:nSubj
        plot(xIDVD, data_med_allSubj(iSubj,:), '-', ...
            'color', ones(1,3)*.7, 'markerfacecolor','w', 'markeredgecolor', ones(1,3)*.7, 'markersize', sz_marker_idvd, 'linewidth', wd);
    end
end

%% [Plot] Reference
if ~isnan(ref)
    yline(ref, '--', 'color', ones(1,3)/2, 'handlevisibility', 'off', 'linewidth', wd);
end

%% [Plot] ticks, limits
xticks(1:nCond); xticklabels([])
% if ~isempty(x_ticklabels), xticklabels(x_ticklabels), end
if ~isnan(y_ticks)
    yticks(y_ticks);
    ylim(y_ticks([1, end]));
end
if ~isnan(y_ticklabels)
    yticklabels(y_ticklabels);
end

if flag_plotIDVD
    buffer = .6;
else
    buffer = .5;
end
xlim([1-buffer, nCond+buffer]);

ax = gca;
ax.XAxis.FontSize = sz_ticks;
ax.YAxis.FontSize = sz_ticks;
ax.LineWidth = wd;

%% [Plot] "diff bar" (nCond==2 only) (must be placed after axis are set)
% Use bootstrap CI of the mean paired difference (not SEM).
if nCond == 2 && flag_plotDiff

    % --- bootstrap CI for mean paired difference (subject bootstrap) ---
    diff_allBoot = nan(nBootPlot, 1);
    for iBoot = 1:nBootPlot
        idx = randi(nSubj, [1 nSubj]);
        diff_allBoot(iBoot) = mean(data_med_allSubj(idx, 1) - data_med_allSubj(idx, 2), 'omitnan');
    end

    % 95% CI for reporting on plot
    CI_plot_diff = 0.95;
    [diff_med_plot, diff_lb_plot, diff_ub_plot] = getCI(diff_allBoot, 1, 1, CI_plot_diff);

    % For the plotted vertical errorbar, use half-lengths around median
    diff_sem_neg = diff_med_plot - diff_lb_plot;   % lower half-length
    diff_sem_pos = diff_ub_plot - diff_med_plot;   % upper half-length

    % --- place the comparison line using YLIM (not data max) ---
    yl = ylim;                      % [ymin ymax]
    yMin = yl(1);
    yMax = yl(2);

    % Base height at 80% of y-axis span
    yBar = yMin + 0.80 * (yMax - yMin);

    % If CI would exceed yMax, nudge downward
    yTop = yBar + diff_sem_pos;
    if yTop > yMax
        yBar = yMax - diff_sem_pos - 0.02 * (yMax - yMin);
    end
    % --- draw horizontal line between the two conditions ---
    plot([1, 2], [yBar, yBar], 'k-', 'LineWidth', wd, 'HandleVisibility', 'off');

    % --- draw CI errorbar at the center ---
    errorbar(1.5, yBar, diff_sem_neg, diff_sem_pos, 'k-', 'LineWidth', wd, 'HandleVisibility', 'off', 'CapSize', 0);

    % --- print point + interval estimate above the line (aligned left) ---
    % left anchor slightly to the right of x=1 to avoid overlap with the bar
    % xText = 1.5;
    % yText = yBar + 0.03 * (yMax - yMin);

    xText = 1.5;                                  % left aligned near bar 1
    yPad  = 0.02 * (yMax - yMin);                  % padding in axis units
    yText = (yBar + diff_sem_pos) + yPad;          % above upper CI

    % If that would exceed yMax, clamp a bit
    if yText > yMax
        yText = yMax - 0.01 * (yMax - yMin);
    end

    str_delta = sprintf('\\Delta=%.2f [%.2f, %.2f]', diff_med_plot, diff_lb_plot, diff_ub_plot);

    text(xText, yText, str_delta, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 20, ...
        'Color', 'k', ...
        'Interpreter', 'tex', ...
        'Clipping', 'off');
end

%% Strings
% (1) Pairwise comparisons
str_diffCond = sprintf('Planned pairwise contrasts (%.1f%%CI)\n', CI_level_diff*100);
for iPair = 1:nPairs
    str_cross0 = '';
    if diff_lb(iPair)*diff_ub(iPair)>0, str_cross0 = '*'; end
    str_diffCond = [str_diffCond, sprintf('%s-%s%s=%.2f [%.2f, %.2f] | t=%.2f, p=%.3f, d=%.2f\n', ...
        x_ticklabels{pairs(iPair,1)}, x_ticklabels{pairs(iPair,2)}, str_cross0, diffCond_med(iPair), diff_lb(iPair), diff_ub(iPair), t_obs_allPairs(iPair), pperm_allPairs(iPair), d_obs_allPairs(iPair))];
end

% (2) Compare to the reference
str_ref = 'No ref';
if ~isnan(ref)
    str_ref = sprintf('Ref=%.1f (%.2f%% CI)\n', ref, CI_level_ref);

    % NHST for each condition vs ref (sign-flip on (X-ref))
    t_ref  = nan(1, nCond);    % signed t
    d_ref = nan(1, nCond);    % dz = mean(diff)/sd(diff)
    pperm_ref = nan(1, nCond);

    % pre-generate sign flips once (consistent null draws)
    sgnMat_ref = (rand(nSubj, nPerm) > 0.5) * 2 - 1;% [nk x nPerm]

    for iCond = 1:nCond
        diffk = data_med_allSubj(:, iCond) - ref;   % [nSubj x 1]
        nk = numel(diffk);
        sd_diff = std(diffk, 0);
        denom = sd_diff / sqrt(nk);

        % observed stats
        t_ref(iCond)  = mean(diffk) / denom;      % signed
        d_ref(iCond) = mean(diffk) / sd_diff;    % signed

        % permutation p (two-tailed)
        num_perm = mean(diffk .* sgnMat_ref, 1);
        t_perm   = abs(num_perm / denom);

        pperm_ref(iCond) = (1 + sum(t_perm >= abs(t_ref(iCond)))) / (nPerm + 1);

        % Build string
        str_cross0 = '';
        if diffRef_lb(iCond)*diffRef_ub(iCond)>0, str_cross0 = '*'; end
        str_ref = [str_ref, sprintf('%s%s: %.2f [%.2f, %.2f] | t=%.2f, p=%.3f, d=%.2f\n', ...
            x_ticklabels{iCond}, str_cross0, data_ave(iCond), data_ave(iCond)-data_sem_neg(iCond), data_ave(iCond)+data_sem_pos(iCond), ...
            t_ref(iCond), pperm_ref(iCond), d_ref(iCond))];

    end % iCond


end

%% [Plot] Title
title(sprintf('%s\n%s\n%s%s\n', str_title, str_ANOVA, str_diffCond, str_ref), 'fontsize', sz_title);


end
