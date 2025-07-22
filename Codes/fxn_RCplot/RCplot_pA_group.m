figure('Position', [0 0 800 300]);

pA_ave = mean(pA_allSubj);
pA_SEM = std(pA_allSubj, [], 1)/sqrt(nsubj);

[pA_p, pA_table, pA_stats] = anova1(pA_allSubj, namesLoc2D, 'off');
fprintf('F(%d,%d) = %.3f, p = %.3f.\n', pA_table{2,3}, pA_table{3,3}, pA_table{2,5}, pA_p)

% 1. barplot
subplot(1,2,1), hold on

% idvd
for isubj_p = 1:nsubj, plot((1:nLoc) +.1, pA_allSubj(isubj_p, :), [marks_allSubj{isubj_p}, '-'], 'color', ones(1,3)*.5), end

% group
for iLoc = 1:nLoc
    bar(iLoc, pA_ave(iLoc), 'FaceColor', colors_comb(iLoc, :), 'facealpha', .3, 'handlevisibility', 'off', 'EdgeColor', 'w')
    errorbar(iLoc, pA_ave(iLoc), pA_SEM(iLoc), '.', 'color', colors_comb(iLoc, :), 'markersize', 10, 'linewidth', 2)
end

xlim([.5, 5.5]), xticklabels(namesLoc2D)
ylim([.5, .9])
ylabel('response consistency')
legend(subjList, 'Location', 'best','Orientation','horizontal')

% 2. polar plot
polarAxesHandle = subplot(1,2,2);
ax = gca;
ax.XTick = [];
ax.YTick = [];
polaraxes('Units',polarAxesHandle.Units,'Position',polarAxesHandle.Position)
hold on
polarplot(polarAng, pA_ave(polarInd)/pA_ave(1), 'k')
ax = gca;
ax.ThetaTick = [];
rlim([.8, 1])
title('relative to the center')



set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
sgtitle('Response consistency', 'FontSize',20)
