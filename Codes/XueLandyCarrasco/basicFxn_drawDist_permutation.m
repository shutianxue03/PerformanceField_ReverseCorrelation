
function basicFxn_drawDist_permutation(data_allIter_allSubj, ref, colors, str_namesCond, x_ticks, x_ticklabels, str_title, sz_fig, nIter, nSubj)

%% Define sizes
fsz_ticks = 15; % tunC: 20; BEHAV: xx; NOM: 15
fsz_title = 10; % font size of titlte
wd = 2; % line width of axis
wd_ref = wd; % line width of the reference line
wd_CI = wd;
fprintf('\n\n *** %s ***\n', str_title)

CI_level = .95; % CI range for stats only! Default for plotting is .68

%% 1) Validate / canonicalize shapes: enforce [nBoot x nSubj x nCond]
assert(ndims(data_allIter_allSubj) == 3, 'ALERT: data_allIter_allSubj must be 3D arrays.');

nCond = size(colors, 1);

data_allIter_allSubj = fxn_reshape(data_allIter_allSubj, nSubj, nCond, 'data_allIter_allSubj');

[nIter2, nSubj2, nCond2] = size(data_allIter_allSubj);
assert(nSubj2 == nSubj && nCond2 == nCond, ...
    'ALERT: Reshaped matrix has wrong nSubj/nCond.');
assert(nIter2 == nIter, ...
    'ALERT: Input nBoot=%d but inferred nBoot=%d from X.', nIter, nIter2);

%% Conduct ANOVA and pairwise comparison PER BOOT
% ANOVA
p_ANOVA_allIter = nan(nIter, 1);
F_allIter = nan(nIter, 1);

% Paired t-test per iteration
p_ttest_allIter = nan(nIter, 1);
t_allIter = nan(nIter, 1);
CohenD_allIter = nan(nIter, 1);

ANOVA_indCond = repmat(1:nCond, nSubj, 1);  % nSubj x nCond

for iBoot = 1:nIter
    data_perBoot = squeeze(data_allIter_allSubj(iBoot, :, :));

    % ANOVA
    [p, tbl] = anovan(data_perBoot(:), ANOVA_indCond(:), 'varnames', 'Condition', 'display', 'off');

    p_ANOVA_allIter(iBoot) = p;
    F_allIter(iBoot) = tbl{2, 6};

    % Pairwise t-test
    if nCond == 2
        [~, p, ~, stats] = ttest(data_perBoot(1, :), data_perBoot(2, :));
        CohenD = fxn_getES(data_perBoot(1, :), data_perBoot(2, :));
        p_ttest_allIter(iBoot) = p;
        t_allIter(iBoot) = stats.tstat;
        CohenD_allIter(iBoot) = CohenD;
    else % >2 levels: ONLY specified pairs, no exhaustive pairwise
        % think of this later
    end

end % iBoot

%% Get group means of all iterations (for comparing to ref)
[data_allIter_ave] = getCI(data_allIter_allSubj, 2, 2);
[~, data_ave_lb, data_ave_ub] = getCI(data_allIter_ave, 1, 1, CI_level);

%% Obtain median and CI
% [data_med, data_lb, data_ub] = getCI(data_allIter_allSubj, 1, 1);

% [ttest_p_med, ttest_p_lb, ttest_p_ub] = getCI(p_ttest_allIter, 1, 1, CI_level);
[ttest_t_med, ttest_t_lb, ttest_t_ub] = getCI(t_allIter, 1, 1, CI_level);
[ttest_d_med, ttest_d_lb, ttest_d_ub] = getCI(CohenD_allIter, 1, 1, CI_level);

% [p_ANOVA_med, p_ANOVA_lb, p_ANOVA_ub] = getCI(p_ANOVA_allIter, 1, 1, CI_level);
[F_med, F_lb, F_ub] = getCI(F_allIter, 1, 1, CI_level);
str_sig = ''; if F_lb*F_ub>0, str_sig = '*'; end

str_ANOVA = sprintf('ANOVA: F%s=%.2f [%.2f, %.2f]', str_sig, F_med, F_lb, F_ub);
str_ttest = sprintf('T-test: t=%.2f [%.2f, %.2f], CohenD=%.2f [%.2f, %.2f]', ttest_t_med, ttest_t_lb, ttest_t_ub, ttest_d_med, ttest_d_lb, ttest_d_ub);

%% Difference
diff_allIter_allSubj = squeeze(data_allIter_allSubj(:, :, 1) - data_allIter_allSubj(:, :, 2));
[diff_ave_allIter] = getCI(diff_allIter_allSubj, 2, 2);
[diff_ave_med, diff_ave_lb, diff_ave_ub] = getCI(diff_ave_allIter, 1, 1, CI_level);
str_sig = ''; if diff_ave_lb*diff_ave_ub>0, str_sig = '*'; end

str_diff = sprintf('Diff%s=%.2f [%.2f, %.2f]', str_sig, diff_ave_med, diff_ave_lb, diff_ave_ub);

%%
figure('Position', [0, 200, sz_fig])
hold on

%% Plot distribution of group means
for iCond = 1:nCond
    histogram(data_allIter_ave(:, iCond), 'FaceColor', colors(iCond, :), 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization','probability')
    xline(data_ave_lb(iCond), '-', 'Color', colors(iCond, :), 'LineWidth', wd_CI);
    xline(data_ave_ub(iCond), '-', 'Color', colors(iCond, :), 'LineWidth', wd_CI);
end % iCond

%% Draw ref and compare bars with ref (print in command window)
str_Comp2Ref = '';
if ~isnan(ref)
    if length(ref)>1; error('ALERT: there are more than one REF!!!'), end
    % yline(ref, 'color', ones(1,3)/2, 'linewidth', wd_ref, 'HandleVisibility', 'off');
    xline(ref, '--',  'color', ones(1,3)/2, 'linewidth', wd_ref, 'HandleVisibility', 'off');

    str_Comp2Ref = sprintf('Ref.=%.1f: ', ref);
    
    for iCond = 1:nCond
        str_sig = ''; if (data_ave_lb(iCond)-ref) * (data_ave_ub(iCond)-ref)>0, str_sig = '*'; end
        str_Comp2Ref = [str_Comp2Ref, sprintf('%s%s: [%.2f, %.2f] | ', str_namesCond{iCond}, str_sig, data_ave_lb(iCond), data_ave_ub(iCond))];
    end
end

%% ticks, limits and labels
if ~isnan(x_ticks), xticks(x_ticks), xlim(x_ticks([1, end])), end
if ~isnan(x_ticklabels), xticklabels(x_ticklabels), end

y_ticks = linspace(0, .4, 5);
yticks(y_ticks)
ylim(y_ticks([1, end]))

%% size
ax = gca;
ax.XAxis.FontSize = fsz_ticks;
ax.YAxis.FontSize = fsz_ticks;
ax.LineWidth = wd;

%% title
title(sprintf('%s\n%s\n%s\n%s\n%s', str_title, str_ANOVA, str_ttest, str_diff, str_Comp2Ref), 'fontsize', fsz_title)



