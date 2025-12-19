
function basicFxn_drawBars_boot(data_allBoot_allSubj, data_obs_allSubj, ref, colors, x_ticks, y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nBoot, markers_allSubj)

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
sz_marker_ave = 30;
fsz_ticks = 15; % tunC: 20; BEHAV: xx; NOM: 15
fsz_title = 10; % font size of titlte
wd = 2; % line width of axis
wd_bar = .5; % the width of the bar (not the bar edge!!)
wd_ref = wd; % line width of the reference line
% markers_allSubj = {'o', 's', 'd', '^','v',  '<', '+','p', 'h', 'x', '>',      'o', 's', 'd', '^'}; % for each subj
interval_diffBar = 8; % higher, closer the comparison bar is to the top of the figure
% markers_allSubj = {'o', 'o', 'o', 'o','o',  'o', 'o','o', 'o', 'o', 'o',      'o', 'o', 'o'};
nsubj_max = 11;
fprintf('\n\n *** %s ***\n', str_title)

CI_level = .95; % CI range for stats only! Default for plotting is .68

%% 1) Validate / canonicalize shapes: enforce [nBoot x nSubj x nCond]
assert(ndims(data_allBoot_allSubj) == 3, 'ALERT: data_allBoot_allSubj must be 3D arrays.');
assert(ndims(data_obs_allSubj) == 2, 'ALERT: data_obs_allSubj must be 2D arrays.');

nCond = size(colors, 1);
nSubj = length(markers_allSubj);

data_allBoot_allSubj = fxn_reshape(data_allBoot_allSubj, nSubj, nCond, 'data_allBoot_allSubj');

[nBoot2, nSubj2, nCond2] = size(data_allBoot_allSubj);
assert(nSubj2 == nSubj && nCond2 == nCond, ...
    'ALERT: Reshaped matrix has wrong nSubj/nCond.');
assert(nBoot2 == nBoot, ...
    'ALERT: Input nBoot=%d but inferred nBoot=%d from X.', nBoot, nBoot2);

if size(data_obs_allSubj, 1) ~= nSubj, data_obs_allSubj = data_obs_allSubj'; end
[nSubj2, nCond2] = size(data_obs_allSubj);
assert(nSubj2 == nSubj && nCond2 == nCond, ...
    'ALERT: Reshaped matrix has wrong nSubj/nCond.');

%% Conduct ANOVA and pairwise comparison PER BOOT
% ANOVA
% p_ANOVA_allBoot = nan(nBoot, 1);
F_allBoot = nan(nBoot, 1);
eta2_allBoot = nan(nBoot, 1);

% Paired t-test per boot
% p_ttest_allBoot = nan(nBoot, 1);
t_allBoot = nan(nBoot, 1);
CohenD_allBoot = nan(nBoot, 1);

ANOVA_indCond = repmat(1:nCond, nSubj, 1);  % nSubj x nCond

for iBoot = 1:nBoot
    data_perBoot = squeeze(data_allBoot_allSubj(iBoot, :, :));

    % ANOVA
    [~, tbl] = anovan(data_perBoot(:), ANOVA_indCond(:), 'varnames', 'Condition', 'display', 'off');
    F_allBoot(iBoot) = cell2mat(tbl(2, 6)); % F-statistic

    % eta-squared
    SS_between = cell2mat(tbl(2, 2));
    SS_total = cell2mat(tbl(end, 2));
    eta2_allBoot(iBoot) = SS_between / SS_total;

    % Pairwise t-test
    if nCond == 2
        [~, ~, ~, stats] = ttest(data_perBoot(:, 1), data_perBoot(:, 2));
        CohenD = fxn_getES(data_perBoot(:, 1), data_perBoot(:, 2));
        % p_ttest_allBoot(iBoot) = p;
        t_allBoot(iBoot) = stats.tstat;
        CohenD_allBoot(iBoot) = CohenD;
    else % >2 levels: ONLY specified pairs, no exhaustive pairwise
        % think of this later
    end

end % iBoot

%% Get group means of all boots (for comparing to ref)
[data_allBoot_ave] = getCI(data_allBoot_allSubj, 2, 2);
[~, data_ave_lb, data_ave_ub] = getCI(data_allBoot_ave, 1, 1, CI_level);

%% Obtain baseline F and eta2 by shuffling labels
% ANOVA
F_obs = rm_oneway_F(data_obs_allSubj);

% t-test
diff  = data_obs_allSubj(:, 1)-data_obs_allSubj(:, 2);

% observed t
t_obs = mean(diff, 'omitnan') / (std(diff, 'omitnan') / sqrt(nSubj));

nPerm = 1e4;
F_allPerm = nan(nPerm,1);
t_allPerm = F_allPerm;

for iPerm = 1:nPerm
    data_shuffled = data_obs_allSubj;
    % shuffle condition labels within each subject (keeps each subject's distribution)
    for iSubj = 1:nSubj
        data_shuffled(iSubj, :) = data_shuffled(iSubj, randperm(nCond));
    end
    
    % ANOVA
    F_allPerm(iPerm) = rm_oneway_F(data_shuffled);

    % t-test
    sgn = (rand(nSubj,1) > 0.5)*2 - 1;     % random +/-1
    dp  = diff .* sgn;
    t_allPerm(iPerm) = mean(dp,'omitnan') / (std(dp,'omitnan') / sqrt(nSubj));

end % iPerm

p_ANOVA_perm = (1 + sum(F_allPerm >= F_obs)) / (numel(F_allPerm) + 1);
p_ttest_perm = (1 + sum(abs(t_allPerm) >= abs(t_obs))) / (nPerm + 1);

%% Obtain median and CI
[data_med, data_lb, data_ub] = getCI(data_allBoot_allSubj, 1, 1);

% [ttest_p_med, ttest_p_lb, ttest_p_ub] = getCI(p_ttest_allBoot, 1, 1, CI_level);
[ttest_t_med, ttest_t_lb, ttest_t_ub] = getCI(t_allBoot, 1, 1, CI_level);
[ttest_d_med, ttest_d_lb, ttest_d_ub] = getCI(CohenD_allBoot, 1, 1, CI_level);

% [p_ANOVA_med, p_ANOVA_lb, p_ANOVA_ub] = getCI(p_ANOVA_allBoot, 1, 1, CI_level);
[F_med, F_lb, F_ub] = getCI(F_allBoot, 1, 1, CI_level);
[eta2_med, eta2_lb, eta2_ub] = getCI(eta2_allBoot, 1, 1, CI_level);

str_ANOVA = sprintf('ANOVA: F=%.2f [%.2f, %.2f], p=%.3f, eta2=%.2f [%.2f, %.2f]', ...
    F_med, F_lb, F_ub, ...
    p_ANOVA_perm, ...
    eta2_med, eta2_lb, eta2_ub);
if nCond==2
    str_ttest = sprintf('T-test: t=%.2f [%.2f, %.2f], p=%.3f, CohenD=%.2f [%.2f, %.2f]', ...
        ttest_t_med, ttest_t_lb, ttest_t_ub, ...
        p_ttest_perm, ...
        ttest_d_med, ttest_d_lb, ttest_d_ub);
else
    str_ttest = '';
end

%% Difference
diff_allBoot_allSubj = squeeze(data_allBoot_allSubj(:, :, 1) - data_allBoot_allSubj(:, :, 2));
[diff_ave_allBoot] = getCI(diff_allBoot_allSubj, 2, 2);
[diff_ave_med, diff_ave_lb, diff_ave_ub] = getCI(diff_ave_allBoot, 1, 1, CI_level);
str_sig = ''; if diff_ave_lb*diff_ave_ub>0, str_sig = '*'; end
if nCond==2
    str_diff = sprintf('Diff%s=%.2f [%.2f, %.2f]', str_sig, diff_ave_med, diff_ave_lb, diff_ave_ub);
else
    str_diff = '';
end

%% Obtain group ave and sem for plotting
[data_ave, ~, ~, data_sem] = getCI(data_med, 2, 1);

%%
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

%% Plot idvd data
buffer = .2;
if flag_plotIDVD

    if nCond == 2 % for nBars=2, need space to plot idvd data between two bars
        x = [1+buffer, 2-buffer];
    else % for nBars>2, plot idvd data at the center of each bar
        x = 1:nCond;
    end

    for iSubj = 1:nSubj
        plot(x, data_med(iSubj, :), '-',  ...
            'color',ones(1,3)*.7, 'markerfacecolor', 'w', 'markeredgecolor', ones(1,3)*.7, ...
            'markersize', sz_marker_idvd, 'linewidth', wd)
    end
end

%% Draw ref and compare bars with ref (print in command window)
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

%% Compare two loc (ttest, draw sem of diff)
if nCond==2 && flag_plotDiff
    [~, ~, ~, diff_sem] = getCI(getCI(diff_allBoot_allSubj, 1, 1), 2, 2);

    yDiffSEM = y_ticks(end) - (y_ticks(end) - y_ticks(1))/interval_diffBar;
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

function F = rm_oneway_F(Y)
% Repeated-measures one-way ANOVA F for Condition (within-subject)
% Standard partition: SS_cond, SS_error (subject x condition interaction)

[nSubj, nCond] = size(Y);

grand = mean(Y(:), 'omitnan');
subj_mean = mean(Y, 2, 'omitnan');     % nSubj x 1
cond_mean = mean(Y, 1, 'omitnan');     % 1 x nCond

SS_cond = nSubj * sum((cond_mean - grand).^2, 'omitnan');
SS_subj = nCond * sum((subj_mean - grand).^2, 'omitnan');

% Total SS
SS_tot  = sum((Y - grand).^2, 'all', 'omitnan');

% Error term in RM one-way: subject x condition interaction
SS_err = SS_tot - SS_cond - SS_subj;

df_cond = nCond - 1;
df_err  = (nSubj - 1) * (nCond - 1);

MS_cond = SS_cond / df_cond;
MS_err  = SS_err  / df_err;

F = MS_cond / MS_err;
end