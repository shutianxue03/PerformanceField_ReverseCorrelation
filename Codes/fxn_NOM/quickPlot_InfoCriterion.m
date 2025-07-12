
%% check information criterion
% IC_allB: ni x nIC
figure('Position', [3e3 800 600 200])

hold on
[IC_med, ~, ~, IC_neg, IC_pos] = getCI(IC_allB);

for ii = 1:ni, plot(1:nIC, IC_allB(ii, :), '-', 'color', ones(1,3)*.7), end
errorbar(1:nIC, IC_med, IC_neg, IC_pos, 'ok', 'CapSize', 0, 'Linewidth', 2)

xticks(1:nIC)
xticklabels(namesIC)
xlim([0,4])

title(sprintf('%s - %s - %s - %s', subjName, namesLocComb{iLocComb}, namesModelA{iModelA}, namesModelB{iModelB}))