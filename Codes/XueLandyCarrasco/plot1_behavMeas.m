
sz_marker_idvd = 20;
wd_bar = .5; % the width of the bar (not the bar edge!!)
wd_barEdge = 5;
wd_idvd = 2;
wd_ref = 3;
sz_title = 10;
sz_ticks = 50;

get_ydiff = @(lb, ub) ub-(ub-lb)/6; % get errorbar of difference
flag_plotScatter = 0;

if flag_plotScatter, figure('Position', [0 0 6e2 1e3])
else,% figure('Position', [0 0 3.5e2 4e2]),
    if nLoc==4, figure('Position', [0 0 1e3 7e2])
    else, figure('Position', [0 0 6e2 7e2])
    end
end
hold on
if im <= nmetrics, mm = squeeze(metrics_allSubj(indSubj, :, :, im));
elseif im == 9, mm = cs_allSubj_(indSubj, :, :); % cs (1/cst)
else, mm = log(RT_allSubj_); % RT
end
mm = pA_allSubj;

% get idvd & group ave
[mm_med, ~, ~, mm_neg, mm_pos] = getCI(mm, 1, 2, 1, 1, CI_ratio);
[nsubj, nLoc] = size(mm_med);
% if im==9, mm_med = 1./cs(indSubj, :); mm_med = [mean(mm_med(:, [2,4]) ,2), mean(mm_med(:, [5,3]), 2)]; end % delete
% if im==9, mm_med = 1./cs(indSubj, :); mm_med = mm_med(:, [5,3]); end % delete

mm_ave = mean(mm_med);
mm_sem = std(mm_med)/sqrt(nsubj);

%% BARS
if flag_plotScatter && nLoc==2, subplot(2,1,1), hold on, end

% bars
for iLoc = 1:nLoc
    bar(iLoc, mm_ave(iLoc), 'FaceColor', colors_comb(iLocComb_all(iLoc), :), 'EdgeColor',  colors_comb(iLocComb_all(iLoc), :), 'BarWidth', wd_bar, 'linewidth', wd_barEdge)
    errorbar(iLoc, mm_ave(iLoc), 0, mm_sem(iLoc), '.', 'CapSize', 0, 'color',  colors_comb(iLocComb_all(iLoc), :), 'linewidth', wd_barEdge)
    errorbar(iLoc, mm_ave(iLoc), mm_sem(iLoc), 0, '.', 'CapSize', 0, 'color', 'w', 'linewidth', wd_barEdge)
end

% stats (1-way ANOVA)
indLoc = repmat(1:nLoc, nsubj, 1);
text_ANOVA = print_nANOVA({'Loc'}, mm_med(:), {indLoc(:)}, nsubj, 1);

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
if flag_plotIDVD
    if nLoc==2, for isubj = 1:nsubj, plot([1.35, 1.65], mm_med(isubj, :), ['-', markers_allSubj{isubj}], 'color', ones(1,3)*.5,  'markerfacecolor', 'w', 'linewidth', wd_idvd, 'markersize', sz_marker_idvd), end
        %     if nLoc==2, for isubj = 1:nsubj, plot([1.35, 1.65], mm_med(isubj, :), '-', 'color', ones(1,3)*.7,  'markerfacecolor', 'w', 'linewidth', wd_idvd, 'markersize', sz_marker_idvd), end
    else
        for isubj = 1:nsubj
            plot(1:nLoc, mm_med(isubj, :),  ['-', markers_allSubj{isubj}], 'color', ones(1,3)*.8, 'linewidth', wd_ref, 'markersize', 20)
        end
    end
    
end

% errorbar of difference
if nLoc==2
    diff_y = get_ydiff(ticks_scatter_all{im}(1), ticks_scatter_all{im}(end));
    mm_diff = mm_med(:, 1) - mm_med(:, 2);
    [mm_diff_ave, ~, ~, mm_diff_SEM] = getCI(mm_diff, 2, 1);
    errorbar(1.5, diff_y, mm_diff_SEM, 'k.', 'CapSize', 0, 'linewidth', wd_barEdge)
    plot([1,2], [diff_y, diff_y], 'k-', 'linewidth', wd_barEdge)
end
% extra line
if ~isnan(yline_all{im}), yline(yline_all{im}, 'k-', 'linewidth', wd_ref); end

%% axis setting
% axis square
% x-ticks and labels
xticks(1:nLoc)
xticklabels(namesLocComb{iLocComb_all(iLoc)})

% y-ticks and labels
switch im, case 10, yticklabels(round(exp(ticks_scatter_all{im})*100+500)), end % RT
yticks(ticks_scatter_all{im})
yticklabels(ticks_scatter_all{im})

% limit
xlim([1,nLoc] + [-1, 1])
% ylim(ticks_scatter_all{im}([1,end]))
% yticklabels(ticks_scatter_all{im})

ax = gca;
ax.XAxis.FontSize = sz_ticks;
ax.YAxis.FontSize = sz_ticks;
ax.LineWidth = wd_barEdge;

%% scatter
if flag_plotScatter
    if nLoc==2
        subplot(2,1,2), hold on, box on
        % diagonal line
        plot(ticks_scatter_all{im}([1,end]), ticks_scatter_all{im}([1,end]), 'k-', 'linewidth', wd_barEdge)
        % IDVD data
        errorbar(mm_med(:,1), mm_med(:,2), ...
            mm_neg(:,2), mm_pos(:,2), ...
            mm_neg(:,1), mm_pos(:,1), ...
            '.', 'color', ones(1,3)*.5, 'CapSize', 0, 'linewidth', wd_idvd)
        for isubj = 1:nsubj
            plot(mm_med(isubj,1), mm_med(isubj,2), markers_allSubj{isubj}, ...
                'linestyle', 'none', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ones(1,3)*.5, 'linewidth', wd_idvd, 'MarkerSize', sz_marker_idvd)
        end
        
        xticks(ticks_scatter_all{im})
        yticks(ticks_scatter_all{im})
        
        ylim(ticks_scatter_all{im}([1, end]))
        xlim(ticks_scatter_all{im}([1, end]))
        
        xlabel(namesLocComb{iLocComb_all(iLoc)})
        ylabel(namesLocComb{iLocComb_all(iLoc)})
        axis square
    end % if iLoc==2
end

% set(findall(gcf, '-property', 'FontSize'), 'FontSize', 35)
% set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',3)

%% title
indLoc_ANOVA1 = repmat(1:nLoc, nsubj, 1);
text_ANOVA = print_nANOVA({'Loc'}, mm_med(:), indLoc_ANOVA1(:), nsubj);
if nLoc==2
    title(sprintf('%s\nt(%d)=%.2f, p=%.3f\n%s', namesMetrics_plus2{im}, stats.df, stats.tstat, p, text_ANOVA), 'fontsize', sz_title)
else
    title(sprintf('%s\n%s', namesMetrics_plus2{im}, text_ANOVA), 'fontsize', sz_title)
end

%%
folderName = sprintf('%s/behav/%s/', nameFigFolder, nameFileLoc);
folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end

% saveas(gcf, sprintf('%sn%d_%s.jpg', folderName, nsubj, namesMetrics_plus2{im}))
