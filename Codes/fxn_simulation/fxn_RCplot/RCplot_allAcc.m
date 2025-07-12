
grid on, hold on


for iLoc = 1:nLoc
    if iLoc == 4
        plot(1:nSess, pC_perSess_perLoc(:, iLoc), '--', 'color', colors_loc5(iLoc, :));
    else
        plot(1:nSess, pC_perSess_perLoc(:, iLoc), '-', 'color', colors_loc5(iLoc, :));
    end
end

yline(threshPerf, 'k-');
yticks(.5:.1:.8)
xlabel('Session # ')
ylabel('pC')

if nSess > 1, xlim([.5,nSess + .5])
    if nSess>2, xticks([1, round(nSess/2),nSess]), xticklabels([1, round(nSess/2),nSess]),else, xticks([1,2]), end
else, xlim([.5,2]), xticks(1)
end
legend(namesLoc2D, 'location', 'eastoutside')

%% inside panel; the average
axes('Position',[.15 .6 .1 .15])
box on, hold on
if nSess > 1
    for iLoc = 1:nLoc, errorbar(iLoc, mean(pC_perSess_perLoc(:, iLoc)), std(pC_perSess_perLoc(:, iLoc))/sqrt(nSess),'o', 'color', colors_loc5(iLoc, :)), end
end
yline(threshPerf, 'k-');
% xticks(1:nLoc), xticklabels(namesLoc2D), xtickangle(30)
xticklabels([])
xlim([.5, nLoc+.5])
% yticklabels([])
% yticks([.6, .7, .8])
yticks([])
ylim([.6, .8])
