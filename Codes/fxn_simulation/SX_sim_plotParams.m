titles = {'center', 'sigma', 'alpha', 'baseline'};
nparams = length(titles);
figure('Position', [0 200 1200 400])

for isubplot = 1:nparams
    
    subplot(1,nparams,isubplot), hold on
    for ipair = 1:nSimPairs

        nsubj = simPairs(1,ipair);
        i_nsubj = find(nsubj == nsubj_all);
        
        ntrials = simPairs(2,ipair);
        i_ntrials = find(ntrials == ntrials_all);
        
        x = params_est_allPairs{i_nsubj, i_ntrials}(:,isubplot); 
        
        % individual 
        plot(ipair, x, 'o') 
        % ave 
        x_ave = mean(x);
        x_sem = std(x)/sqrt(nsubj);
        errorbar(ipair, x_ave, x_sem, 'ok')
        % real param
        realParam = params.stim.tuning_params(isubplot);
        plot([.5, nSimPairs+.5], [realParam,realParam], 'color', [.5 .5 .5])
        
        xticklabels_{ipair} = sprintf('n=%d #=%d', nsubj, ntrials);
    end
    
    xticks(1:nSimPairs), xticklabels(xticklabels_), xtickangle(45)
    xlim([.5, nSimPairs+.5])
    title(titles{isubplot})
end

sgtitle(sprintf('center=%.2f, sigma=%.2f, alpha=%.2f, b=%.2f\n', params.stim.tuning_params))
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)