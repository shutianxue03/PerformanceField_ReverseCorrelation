

sz_marker = 20;
sz_ticks = 30;
wd_bar = .4;
wd = 5;
wd_line = 2;

get_ydiff = @(lb, ub) ub-(ub-lb)/6; %??

figure('Position', [0 0 6e2 4e2])

hold on
if im <= nmetrics, mm = squeeze(metrics_allSubj_(:, :, :, im));
elseif im == 9, mm = cs_allSubj_; % cs (1/cst)
else, mm = log(RT_allSubj_); % RT
end

% get idvd & group ave
[mm_med, ~, ~, mm_neg, mm_pos] = getCI(mm, 1, 2, 1, 1, CI_ratio);
mm_ave = mean(mm_med);
mm_sem = std(mm_med)/sqrt(nsubj);

%% BARS
if flag_plotScatter && nLoc==2, subplot(2,1,1), hold on, end

% bars
for iLoc = 1:nLoc
    bar(iLoc, mm_ave(iLoc), 'FaceColor', colors_comb_(iLoc, :), 'EdgeColor', colors_comb_(iLoc, :), 'Linewidth', 2, 'BarWidth', wd_bar, 'linewidth', wd_line)
    errorbar(iLoc, mm_ave(iLoc), mm_sem(iLoc), '.', 'CapSize', 0, 'color', colors_comb_(iLoc, :), 'Linewidth', 2, 'linewidth', wd_line)
end

% stats (1-way ANOVA)
[p, tbl] = anova1(mm_med, namesLocComb_, 'off');

% stats: t-test
if nLoc==2
    [~, p,~, stats] = ttest(mm_med(:, 1), mm_med(:, 2));
end

% compare with the ref
if ~isnan(yline_all{im})
    [h, p_ref, ~, stats_ref] = ttest(mm_med(:), yline_all{im});
    fprintf('%s vs. %.2f: t(%d) = %.3f, p = %.3f\n', namesMetrics_plus2{im}, yline_all{im}, stats_ref.df, stats_ref.tstat, p_ref)
end

%% grey lines: idvd
for isubj = 1:nsubj, plot([1.35, 1.65], mm_med(isubj, :), ['-', markers_allSubj{isubj}], 'color', ones(1,3)*.5,  'markerfacecolor', 'w', 'linewidth', wd_line, 'markersize', sz_marker), end

% errorbar of difference
if nLoc==2
    diff_y = 4.3;%get_ydiff(ticks_scatter_all{im}(1), ticks_scatter_all{im}(end));
    mm_diff = mm_med(:, 1) - mm_med(:, 2);
    [mm_diff_ave, ~, ~, mm_diff_SEM] = getCI(mm_diff, 2, 1);
    errorbar(1.5, diff_y, mm_diff_SEM, 'k.', 'CapSize', 0, 'linewidth', wd)
    plot([1,2], [diff_y, diff_y], 'k-', 'linewidth', wd)
end
% extra line
if ~isnan(yline_all{im}), yline(yline_all{im}, 'k-', 'linewidth', wd_line); end

%% axis setting
% axis square
% x-ticks and labels
xticks([])
% xticklabels(namesLocComb_)


% y-ticks and labels
switch im, case 10, yticklabels(round(exp(ticks_scatter_all{im})*100+500)), end % RT
yticks(ticks_scatter_all{im})
yticklabels(ticks_scatter_all{im})

% limit
xlim([1,nLoc] + [-.5, .5])
ylim(ticks_scatter_all{im}([1,end]))
yticklabels(ticks_scatter_all{im})
%%

    ax = gca;
    ax.XAxis.FontSize = sz_ticks;
    ax.YAxis.FontSize = sz_ticks;
    ax.LineWidth = wd;

%%
folderName = sprintf('%s/%s/behav/%s/', nameFigFolder, nameEnergySource, nameFileLoc);
folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end

saveas(gcf, sprintf('%sn%d_%s_grant.jpg', folderName, nsubj, namesMetrics_plus2{im}))
