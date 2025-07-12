titles = {'dprime', 'criterion'};
xticklabels_ = cell(nSimPairs,1);

figure('Position', [0 200 1200 300 ])

for isubplot = 1:length(titles)
    if isubplot == 1, x = dprime_allSubj; else, x = criterion_allSubj;end
    subplot(1,3,isubplot), hold on
    for ipair = 1:nSimPairs
        
        nsubj = simPairs(1,ipair);
        i_nsubj = find(nsubj == nsubj_all);
        
        ntrials = simPairs(2,ipair);
        i_ntrials = find(ntrials == ntrials_all);
        
        tgtCST = simPairs(3, ipair);
        iTgtCST = find(tgtCST == tgtCST_all);
        
        plot(ipair, x{i_nsubj, i_ntrials, iTgtCST}, 'o') % individual data
        x_ave = mean(x{i_nsubj, i_ntrials, iTgtCST});
        x_sem = std(x{i_nsubj, i_ntrials})/sqrt(nsubj);
        errorbar(ipair, x_ave, x_sem, 'ok')
        xticklabels_{ipair} = sprintf('tgt cst=%.2f #=%d', tgtCST, ntrials);
    end
    if isubplot == 2, plot([.5, nSimPairs+.5], [0,0], 'color', [.5 .5 .5]), ylim([-1, 1]), end
    xticks(1:nSimPairs), xticklabels(xticklabels_), xtickangle(45)
    xlim([.5, nSimPairs+.5])
    title(titles{isubplot})
end

%% linearity check
subplot(1,3,3), hold on
x_ = cell2mat(dprime_allSubj);
x = squeeze(mean(x_, [1,2]));
plot([0, tgtCST_all], [0, x'], 'o-')
xlim([0, tgtCST_all(end)]), xticks(tgtCST_all), xticklabels(tgtCST_all)
ylim([0, x(end)])
xlabel('target contrast')
ylabel('d prime')
title('linearity check')
%%
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)

