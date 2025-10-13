
thresh_PMF_allLoc = nan(nLoc, nModels);
figure('Position', [0 0 2e3 2e3])
for iLoc = 1:nLocSingle

    subplot(3,3, iplots5(iLoc))
    
    hold on
    scaling=1/2;
    %-----------------------%
    fxn_plotPMF_singlePanel
    %-----------------------%
    title(sprintf('L%d', iLoc))
    
    thresh_PMF_allLoc(iLoc, :) = getCI(thresh_all{iLoc}, 1, 1);
    
end % iLoc_tgt

set(findall(gcf, '-property', 'fontsize'), 'fontsize',12)
set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
sgtitle(sprintf('%s %.0f%% [Bin%dFilter%d]', subjName, perfThresh_plot*100, flag_binData, flag_filterData))
saveas(gcf, sprintf('%s/%s_PMF_Bin%dFilter%d_%.0f.jpg', nameFolder_fig_thresh, subjName, flag_binData, flag_filterData, perfThresh_plot*100))

