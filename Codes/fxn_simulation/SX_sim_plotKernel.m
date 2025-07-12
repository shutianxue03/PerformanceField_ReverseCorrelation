figure('Position', [0 200 n_ntrials*300 nTgtCST*300 ])

colors = rand(3,max(nsubj_all));

for ipair = 1:nSimPairs
    nsubj = simPairs(1,ipair);
    i_nsubj = find(nsubj == nsubj_all);
    
    ntrials = simPairs(2,ipair);
    i_ntrials = find(ntrials == ntrials_all);
    ntrialsAll = ntrials*2;
    
    tgtCST = simPairs(3, ipair);
    iTgtCST = find(tgtCST == tgtCST_all);
    
    
    kernel_ave = mean(cat(1,kernel_allSubj{:, i_ntrials, iTgtCST}));
    kernel_sem = std(cat(1,kernel_allSubj{:, i_ntrials, iTgtCST}))/sqrt(nsubj);
    %     [kernel_ave_pred, params_ave_est, nLL_ave] = SX_sim08_fit(filterSF_all, kernel_ave, model, fitMode);
    
    subplot(nTgtCST, n_ntrials, n_ntrials*(iTgtCST-1)+i_ntrials), hold on
    % plot individual data
    for isubj = 1:nsubj
        plot(filterSF_all_log, kernel_allSubj{isubj, i_ntrials, iTgtCST}, 'o', 'color', colors(:, isubj))
        %         plot(filterSF_all_log, kernel_pred_allPairs{isubj, i_ntrials, iTgtCST}, '-', 'color', colors(:, isubj))
    end
    % plot ave data
    errorbar(filterSF_all_log, kernel_ave, kernel_sem, 'o-k')
    %     plot(filterSF_all_log, kernel_ave_pred, '-k', 'linewidth', 2)
    
    xticks(log([1,2,4])), xticklabels([1,2,4])
    ylimit = ylim;
    xlim(log([.8,4.8]))
    plot(log([.8,4.8]), [0 0], 'k', 'color',[.5,.5,.5])
    plot(log([2,2]), ylimit, 'k', 'color',[.5,.5,.5])
    
    if i_nsubj+i_ntrials == 2, xlabel('SF channel (cpd)'), ylabel('kernel (a.u.)'), end
    title(sprintf('tgt cst = %.2f, %d trials', tgtCST, ntrials*2))
    %     text(log(1), ylimit(1), sprintf('%.2f, %.2f, %.2f, %.2f\n', params_ave_est))
    %     text(log(1.5), ylimit(1), sprintf('center=%.2f\nsigma=%.2f\nalpha=%.2f\nb=%.2f\n', params_ave_est))
    
end

% sgtitle(sprintf('center=%.2f, sigma=%.2f, alpha=%.2f, b=%.2f\n', params.stim.tuning_params))
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
sgtitle(sprintf('kernel n=%d nCST = %.2f', nsubj, noise.noiseCST), 'fontsize', 20)

