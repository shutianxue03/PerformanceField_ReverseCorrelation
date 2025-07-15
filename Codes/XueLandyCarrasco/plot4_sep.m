
wd_border = 2;
sz_title = 40;
sz_label = 25; %40
sz_ticks = 35;
wd_all = 2;
sz_marker = 25;

[sep_med, ~, ~, sep_neg, sep_pos] = getCI(data_3allSubj(:, :, :, itype), 1, 2, 1,1);
[nsubj, nLoc] = size(sep_med);

figure('Position', [0 0 8e2 1.5e3])
hold on, box on
yline(.5, '--', 'linewidth', wd_all);

title_ = [];
for iLoc = 1:nLoc
    x = (1:nsubj) + .1 * iLoc;
    for isubj = 1:nsubj
        plot(x(isubj), sep_med(isubj, iLoc), markers_allSubj{isubj}, ...
            'MarkerEdgeColor', colors_comb(iLocComb_all(iLoc), :), 'MarkerFaceColor', 'w', 'MarkerSize', sz_marker, 'linewidth', 3)
    end
    title_ = [title_, sprintf('Loc#%d: %.2f (%.2f); ', iLoc, mean(sep_med(:, iLoc)), std(sep_med(:, iLoc))/sqrt(nsubj))];
    
end

ylabel('Correlation', 'FontSize', sz_label)
yticks([0,.5:.1:1]), ylim([0, 1])

xticks(1:nsubj), xticklabels(1:nsubj), xlim([0,nsubj+1])
xlabel('Observer #')

ax = gca; ax.XAxis.FontSize = sz_ticks; ax.YAxis.FontSize = sz_ticks; ax.LineWidth = wd_border;

%%% ANOVA
indLoc = repmat(1:nLoc, nsubj, 1);
text_ANOVA = print_nANOVA({'Loc'}, sep_med(:), {indLoc(:)}, nsubj);
%=== CI of ANOVA ===
for iB = 1:nB
    x=squeeze(data_3allSubj(:, iB, :, itype));
    [~, tbl] = print_nANOVA({'Loc'}, x(:), {indLoc(:)}, nsubj);
    p_allB(iB) = tbl{2,7};
end
[p_med, p_lb, p_ub] = getCI(p_allB, 1, 2); 
text_ANOVA = [text_ANOVA, sprintf('p=%.3f [%.3f, %.3f]\n', p_med, p_lb, p_ub)];
%===============
        
%%% ttest
if nLoc==2
    [h,p,ci,stats] = ttest(sep_med(:, 1), sep_med(:, 2));
    cohenD = fxn_getES(sep_med(:, 1), sep_med(:, 2));
end

%%%
folderName = sprintf('%s/%s/%s/', nameFigFolder, name_numFilters_Fitting, n);
folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
set(findall(gcf, '-property', 'FontSize'), 'FontSize', sz_ticks)

if nLoc==2
    title(sprintf('Loc#1: %.2f (%.2f); Loc#2: %.2f (%.2f)\nt(%d)=%.2f, p=%.2f, d=%.3f\n',...
        mean(sep_med(:, 1)), std(sep_med(:, 1)), mean(sep_med(:, 2)), std(sep_med(:, 2)), stats.df, stats.tstat, p, cohenD), 'FontSize', 30)
    saveas(gcf, sprintf('%sn%d_%s_%s.jpg', folderName, nsubj, nameFileLoc, namesType{itype}))
else
    title(sprintf('%s\n%s\n%s', n, title_, text_ANOVA), 'FontSize', 20)
    saveas(gcf, sprintf('%sn%d_%s_%s.jpg', folderName, nsubj, 'L653', namesType{itype}))
end



