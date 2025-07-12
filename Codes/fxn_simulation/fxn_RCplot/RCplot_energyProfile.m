

load(energyMatName, 'energy_allLoc', 'energy_norm_allLoc')

figure('Position', [0 800 1200 400])
for plotMode = [1,2]
    for iLoc = 1:nLoc
        subplot(2, nLoc, nLoc*(plotMode-1) + iLoc), hold on
        RCplot_energy(energy_norm_allLoc{iLoc}, data_both(:, iLoc, 1), filterSF_all, plotMode)
        if iLoc == 3,  if plotMode==1, title('plotted by SF'), else, title('plotted by trial#'),end, end
    end
end
sgtitle('energy-profile (y-axis is energy)', 'FontSize',15)
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)