

figure('Position', [0 200 1500 600])

for iLoc = 1:nLoc8
    
    subplot(2,4, iLoc)
    imagesc(filtersOri_all-90, filtersSF_all_log, squeeze(mean(stimPhase2D_allT_perComb{iLoc}, 1)))
    
    axis square
    
    cc = colorbar;
    caxis([-pi/8, pi/8])
    cc.YTick = [-pi/8, 0, pi/8];
    cc.YTickLabel = {'-\pi/8', '0', '\pi/8'};
    
    xline(0, 'r-', 'linewidth', 2);
    yline(1, 'r-', 'linewidth', 2);
    
    xticks(axisTicks_tuning{1})
    xlim([-90, 90])
    yticks(axisTicks_tuning{2})
    yticklabels(axisTL_tuning{2})
    ylim([0,2])
    ylabel('SF (cpd)'), xlabel('ORI (deg)')
    
    title(namesLocComb{iLoc})
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
sgtitle('Stim Phase', 'FontSize',20)

