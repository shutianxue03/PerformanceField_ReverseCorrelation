

figure('Position', [0 500 nparams*300 400])

for iparam = 1:nparams
   
    x = squeeze(params(:, itype, :, iparam)); % check the shape, each row should a subj
    if nsubj==1, x=x'; end
    
    %% dot plot
    subplot(2, nparams, iparam), hold on
    
    % plot idvd data
    for isubj_p = 1:nsubj, plot(1:nLoc, x(isubj_p, :), 'o-', 'MarkerFaceColor', colors(:, isubj_p), 'color', colors(:, isubj_p), 'MarkerEdgeColor', 'w'), end
    
    % plot ave data
    x_ave = nanmean(x,1);
    x_sem = nanstd(x, [], 1)/sqrt(nsubj);
    errorbar(1:nLoc, x_ave, x_sem, '-ok', 'MarkerFaceColor', 'k', 'linewidth', 2)
    
    if iparam == 1, legend(subjList), end
    if iparam == 3, yline(0); end
    xticks(1:nLoc), xticklabels(locNames), xtickangle(15)
    xlim([.5, nLoc+.5])
    title(titlesParams{iparam})
    
    %% polar plot
    polarAxesHandle = subplot(2, nparams, iparam+nparams); hold on
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    polaraxes('Units',polarAxesHandle.Units,'Position',polarAxesHandle.Position)
    hold on
    % plot idvd data
    for isubj_p = 1:nsubj, polarplot(polarAng, x(isubj_p, polarInd), 'o-', 'MarkerFaceColor', colors(:, isubj_p), 'color', colors(:, isubj_p), 'MarkerEdgeColor', 'w'), end
    %  plot average data
    polarplot(polarAng, x_ave(polarInd), '-ok', 'linewidth', 2)
    ax = gca;
    ax.ThetaTick = [];
    
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)

if nsubj==1, sgtitle(sprintf('%s- %s params', subjName, kernelNames{ikernel}), 'FontSize',25)
else, sgtitle(sprintf('%s params, tgt-%s (n=%d)', kernelNames{ikernel}, typeNames{itype}, nsubj), 'FontSize',25)
end






