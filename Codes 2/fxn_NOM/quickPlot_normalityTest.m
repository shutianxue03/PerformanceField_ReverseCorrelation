figure('Position', [3e3 0 600 400])
bar(squeeze(mean(h_allB==0, 1)))
xticks(1:length(namesNormalityTest))
xticklabels(namesNormalityTest)
legend({'PRS', 'ABS'})
ylim([0,1])
title(subjName)
set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)


fprintf(' done\n')