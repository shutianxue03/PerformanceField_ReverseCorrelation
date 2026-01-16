
function basicFxn_drawBars_permutation(data_allIter_allSubj, ref, colors, x_ticks, y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIter, markers_allSubj)

% Inputs:
%    med_allSubj: nsubj x nBars, median of boostrapping
%    ref: one reference line, ylines
%    colors: nBars x 3, each rows indicates a bar
%    x_ticks: cell, containing strings indicating location names
%    y_ticks: vector, containing 3 values
%    y_ticklabels: vector, containing 3 value
%    flag_plotIdvd: 1=plot idvd data; 0=NOT
%    flag_plotDiff: ~isnan=plot SEM of difference
%    text_title: string
%    flag_pairwiseComp: 1=print pairwise comparison in the title
%    sz_fig

%% define sizes
sz_marker_idvd = 10;
fsz_ticks = 15; % tunC: 20; BEHAV: xx; NOM: 15
fsz_title = 10; % font size of titlte
wd = 2; % line width of axis
wd_bar = .5; % the width of the bar (not the bar edge!!)
wd_ref = wd; % line width of the reference line
interval_diffBar = 8; % higher, closer the comparison bar is to the top of the figure
fprintf('\n\n *** %s ***\n', str_title)

CI_level = .95; % CI range for stats only! Default for plotting is .68

%% Enforce shape of matrix to [nIter x nSubj x nCond]
assert(ndims(data_allIter_allSubj) == 3, 'ALERT: data_allIter_allSubj must be 3D arrays.');

nCond = size(colors, 1);
nSubj = length(markers_allSubj);

data_allIter_allSubj = fxn_reshape(data_allIter_allSubj, nSubj, nCond, 'data_allIter_allSubj');

[nIter2, nSubj2, nCond2] = size(data_allIter_allSubj);
assert(nSubj2 == nSubj && nCond2 == nCond, ...
    'ALERT: Reshaped matrix has wrong nSubj/nCond.');
assert(nIter2 == nIter, ...
    'ALERT: Input nIter=%d but inferred nIter=%d from X.', nIter, nIter2);

%% Conduct ANOVA and paired t-test per iteration
% To obttain test statistics and effect size

% ANOVA
Fvalue_allIter = nan(nIter, 1);
eta2p_allIter = nan(nIter, 1);
% Paired t-test
tscore_allIter = nan(nIter, 1);
CohenD_allIter = nan(nIter, 1);

for iIter = 1:nIter
    data_perIter = squeeze(data_allIter_allSubj(iIter, :, :));

    % ANOVA
    [Fvalue_allIter(iIter), eta2p_allIter(iIter)] = rm_oneway(data_perIter);

    % Paired t-test
    if nCond == 2
        [~, ~, ~, stats] = ttest(data_perIter(:, 1), data_perIter(:, 2));
        CohenD = fxn_getES(data_perIter(:, 1), data_perIter(:, 2));
        tscore_allIter(iIter) = stats.tstat;
        CohenD_allIter(iIter) = CohenD;
    else % >2 levels: ONLY specified pairs, no exhaustive pairwise
        % think of this later
    end

end % iIter

%% Permutation to obtain p-value
% Obtain median and CI for stats
[data_med_allSubj] = getCI(data_allIter_allSubj, 1, 1);

nPerm = 1e4; % number of permutation runs
Fvalue_allPerm = nan(nPerm,1);
tscore_allPerm = Fvalue_allPerm;

% ANOVA: observed F-value
Fvalue_obs = rm_oneway(data_med_allSubj);

% t-test: observed t-score
if nCond==2
    diff_med_allSubj  = data_med_allSubj(:, 1)-data_med_allSubj(:, 2);
    tscore_obs = mean(diff_med_allSubj,'omitnan') / (std(diff_med_allSubj,0,'omitnan') / sqrt( sum(~isnan(diff_med_allSubj))));
end

for iPerm = 1:nPerm
    data_shuffled = data_med_allSubj;

    % shuffle condition labels within each subject (keeps each subject's distribution)
    for iSubj = 1:nSubj
        data_shuffled(iSubj, :) = data_shuffled(iSubj, randperm(nCond));
    end

    % ANOVA
    Fvalue_allPerm(iPerm) = rm_oneway(data_shuffled);

    % t-test
    if nCond==2
        sgn = (rand(nSubj,1) > 0.5)*2 - 1;     % random +/-1
        dp  = diff_med_allSubj .* sgn;
        tscore_allPerm(iPerm) = mean(dp,'omitnan') / (std(dp,0,'omitnan')/sqrt(sum(~isnan(dp))));
    end
end % iPerm

p_ANOVA_perm = (1 + sum(Fvalue_allPerm >= Fvalue_obs)) / (numel(Fvalue_allPerm) + 1);
if nCond==2
    p_ttest_perm = (1 + sum(abs(tscore_allPerm) >= abs(tscore_obs))) / (nPerm + 1);
end

%% Obtain median and CI for stats
% F value and partial eta2
[Fvalue_med, Fvalue_lb, Fvalue_ub] = getCI(Fvalue_allIter, 1, 1, CI_level);
[eta2p_med, eta2p_lb, eta2p_ub] = getCI(eta2p_allIter, 1, 1, CI_level);

% t-score and cohen's D
[tscore_med, tscore_lb, tscore_ub] = getCI(tscore_allIter, 1, 1, CI_level);
[cohenD_med, cohenD_lb, cohenD_ub] = getCI(CohenD_allIter, 1, 1, CI_level);

%% Put up strings
str_ANOVA = sprintf('ANOVA: F=%.2f [%.2f, %.2f], p=%.3f, eta2p=%.2f [%.2f, %.2f]', ...
    Fvalue_med, Fvalue_lb, Fvalue_ub, ...
    p_ANOVA_perm, ...
    eta2p_med, eta2p_lb, eta2p_ub);
if nCond==2
    str_ttest = sprintf('T-test: t=%.2f [%.2f, %.2f], p=%.3f, CohenD=%.2f [%.2f, %.2f]', ...
        tscore_med, tscore_lb, tscore_ub, ...
        p_ttest_perm, ...
        cohenD_med, cohenD_lb, cohenD_ub);
else
    str_ttest = '';
end

% Difference
if nCond==2
    diff_allIter_allSubj = squeeze(data_allIter_allSubj(:, :, 1) - data_allIter_allSubj(:, :, 2));
    [diff_ave_allIter] = getCI(diff_allIter_allSubj, 2, 2);
    [diff_ave_med, diff_ave_lb, diff_ave_ub] = getCI(diff_ave_allIter, 1, 1, CI_level);
    str_sig = ''; if diff_ave_lb*diff_ave_ub>0, str_sig = '*'; end
    str_diff = sprintf('Diff%s=%.2f [%.2f, %.2f]', str_sig, diff_ave_med, diff_ave_lb, diff_ave_ub);
else
    str_diff = '';
end

%% Obtain group ave and sem for plotting
[data_ave, ~, ~, data_sem] = getCI(data_med_allSubj, 2, 1);

%% PLOT
figure('Position', [0, 200, sz_fig])
hold on

%% Plot group means (bars) and SEM (errorbars)
for iCond = 1:nCond

    bar(iCond, data_ave(iCond), 'FaceColor', colors(iCond, :),  'EdgeColor', colors(iCond, :), 'barwidth', wd_bar, 'HandleVisibility', 'off')
    errorbar(iCond, data_ave(iCond), data_sem(iCond), '.', 'color', colors(iCond, :), 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off')
    % Fancy errorbars
    if data_ave(iCond)>0
        errorbar(iCond, data_ave(iCond), data_sem(iCond), 0, '.w', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off')
    else
        errorbar(iCond, data_ave(iCond), 0, data_sem(iCond), '.w', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off')
    end

    % if plot HVA vs. VMA
    %     bar(iCond, ave(iCond), 'FaceColor', 'w',  'EdgeColor', 'k', 'barwidth', wd_bar, 'linewidth', wd)
    %     errorbar(iCond, ave(iCond), sem(iCond), '.', 'color', 'k', 'CapSize', 0, 'linewidth', wd)
end % iCond

%% Plot idvd data (grey lines)
buffer = .2;
if flag_plotIDVD

    if nCond == 2 % for nBars=2, need space to plot idvd data between two bars
        x = [1+buffer, 2-buffer];
    else % for nBars>2, plot idvd data at the center of each bar
        x = 1:nCond;
    end

    for iSubj = 1:nSubj
        plot(x, data_med_allSubj(iSubj, :), '-',  ...
            'color',ones(1,3)*.7, 'markerfacecolor', 'w', 'markeredgecolor', ones(1,3)*.7, ...
            'markersize', sz_marker_idvd, 'linewidth', wd)
    end
end

%% Draw ref and compare bars with ref
% Get group means of all iterations
[data_allIter_ave] = getCI(data_allIter_allSubj, 2, 2);
[~, data_ave_lb, data_ave_ub] = getCI(data_allIter_ave, 1, 1, CI_level);

str_Comp2Ref = '';
if ~isnan(ref)
    if length(ref)>1; error('ALERT: there are more than one REF!!!'), end
    yline(ref, 'color', ones(1,3)/2, 'linewidth', wd_ref, 'HandleVisibility', 'off');

    str_Comp2Ref = sprintf('Ref.=%.1f: ', ref);

    for iCond = 1:nCond
        str_sig = ''; if (data_ave_lb(iCond)-ref) * (data_ave_ub(iCond)-ref)>0, str_sig = '*'; end
        if iCond==4, str_Comp2Ref = sprintf('%s\n', str_Comp2Ref); end
        str_Comp2Ref = [str_Comp2Ref, sprintf('%s%s: [%.2f, %.2f] | ', x_ticks{iCond}, str_sig, data_ave_lb(iCond), data_ave_ub(iCond))];
    end
end

%% Plot SEM of difference (one errorbar)
if nCond==2 && flag_plotDiff
    yDiffSEM = y_ticks(end) - (y_ticks(end) - y_ticks(1))/interval_diffBar;
    [~, ~, ~, diff_sem] = getCI(getCI(diff_allIter_allSubj, 1, 1), 2, 2);
    errorbar(1.5, yDiffSEM, diff_sem, 'k', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off')
    plot([1,2], [yDiffSEM, yDiffSEM], 'k-', 'linewidth', wd, 'HandleVisibility', 'off')
    %     string_s = getString_starts(p);
end

%% ticks, limits and labels
if ~isnan(y_ticks), yticks(y_ticks), ylim(y_ticks([1, end])), end
if ~isnan(y_ticklabels), yticklabels(y_ticklabels), end
xticks(1:nCond), xticklabels([]) % manually add tick labels /symbols on the slide

% if plot idvd data, leave more space for the middle
if flag_plotIDVD, buffer = .6; else, buffer = .5; end
xlim([1-buffer, nCond+buffer])

%% size
ax = gca;
ax.XAxis.FontSize = fsz_ticks;
ax.YAxis.FontSize = fsz_ticks;
ax.LineWidth = wd;

%% title
title(sprintf('%s\n%s\n%s\n%s\n%s', str_title, str_ANOVA, str_ttest, str_diff, str_Comp2Ref), 'fontsize', fsz_title)

end

function [Fvalue, eta2_partial] = rm_oneway(data)
% Repeated-measures one-way ANOVA F for Condition (within-subject)
% Standard partition: SS_cond, SS_error (subject x condition interaction)

[nSubj, nCond] = size(data);

grand = mean(data(:), 'omitnan');
subj_mean = mean(data, 2, 'omitnan');     % nSubj x 1
cond_mean = mean(data, 1, 'omitnan');     % 1 x nCond

SS_cond = nSubj * sum((cond_mean - grand).^2, 'omitnan');
SS_subj = nCond * sum((subj_mean - grand).^2, 'omitnan');

% Total SS
SS_tot  = sum((data - grand).^2, 'all', 'omitnan');

% Error term in RM one-way: subject x condition interaction
SS_err = SS_tot - SS_cond - SS_subj;

df_cond = nCond - 1;
df_err  = (nSubj - 1) * (nCond - 1);

MS_cond = SS_cond / df_cond;
MS_err  = SS_err  / df_err;

Fvalue = MS_cond / MS_err;

% Compute the effect size: partial eta-squared
eta2_partial = SS_cond/(SS_cond+SS_err);
end