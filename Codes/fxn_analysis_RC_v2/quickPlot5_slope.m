% close all
figure('Position', [200 500 1200 300])
for itype = 1:ntypes
    subplot(1,ntypes,itype), hold on
    plot(ebin_tgt(itype,:), pYES_tgt(itype,:), 'k-')
    ylim([0,1])
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
end