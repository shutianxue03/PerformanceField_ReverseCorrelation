

eff_allSubj = squeeze(efficiency_allSubj(:,:,1));

figure('position', [0 0 800 300])

%% 1. bar plot
subplot(1,2,1),hold on

% idvd data
for isubj_p = 1:nsubj, plot((1:nLoc) + .1, eff_allSubj(isubj_p, :), [marks_allSubj{isubj_p}, 'k']), end

% ave
eff_ave = mean(eff_allSubj,1);
if nsubj>1, eff_err = std(eff_allSubj)/sqrt(nsubj); else, eff_err = zeros(1,nLoc); end

for iLoc = 1:nLoc
    bar(iLoc, eff_ave(iLoc), 'FaceColor', colors_comp(iLoc, :), 'facealpha', .3, 'handlevisibility', 'off', 'EdgeColor', 'w')
    errorbar(iLoc, eff_ave(iLoc), eff_err(iLoc), '.', 'color', colors_comp(iLoc, :), 'markersize', 10, 'linewidth', 2)
end

xticks(1:nLoc), xticklabels(locNames)
xlim([.5, nLoc+.5])
legend(subjList, 'Location', 'best','Orientation','horizontal')

%% 2. polar plot
subplot(1,2,2),hold on
polarAxesHandle = subplot(1,2,2);
ax = gca;
ax.XTick = [];
ax.YTick = [];
polaraxes('Units',polarAxesHandle.Units,'Position',polarAxesHandle.Position)
hold on
polarplot(polarAng, eff_ave(polarInd)/eff_ave(1), 'k')
ax = gca;
ax.ThetaTick = [];
% rlim([.8, 1])
title('relative to the center')

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
sgtitle('efficiency', 'FontSize',20)


