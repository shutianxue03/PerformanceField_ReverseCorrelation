% close all
figure('Position', [200 500 1200 300])
for iType = 1:nTypes
    subplot(1,nTypes,iType), hold on
    plot(ebin_tgt(iType,:), pYES_tgt(iType,:), 'k-')
    ylim([0,1])
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
end