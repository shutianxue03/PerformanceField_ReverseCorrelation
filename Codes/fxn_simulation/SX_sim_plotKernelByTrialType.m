

figure('Position', [0 200 n_ntrials*300 nGaborCST*300 ])

for ipair = 1:nSimPairs
    nsubj = simPairs(1,ipair);
    i_nsubj = find(nsubj == nsubj_all);
    
    ntrials = simPairs(2,ipair);
    i_ntrials = find(ntrials == ntrials_all);
    ntrialsAll = ntrials*2;
    
    stim.tgtCST = simPairs(3, ipair);
    iTgtCST = find(stim.tgtCST == tgtCST_all);
    
    kernel1_all = [];
    kernel2_all = [];
    kernel_both_all = [];
    
    for isubj = 1:nsubj
        energy_norm = energy_norm_allSubj{isubj, i_ntrials, iTgtCST};
        behav = behav_allSubj{isubj, i_ntrials, iTgtCST};
        [kernel1, ~] = SX_sim07_RC(nfilters, energy_norm(1:ntrials, :), behav(1:ntrials).');
% filtersSF_all, filtersOri_all, energy_both, y)
        [kernel2, ~] = SX_sim07_RC(nfilters, energy_norm(ntrials+1:end, :), behav(ntrials+1:end).');
        kernel_both = SX_sim07_RC(nfilters, energy_norm, behav');
%         kernel = kernel_allSubj{isubj, i_ntrials, iTgtCST};
        kernel1_all = cat(1, kernel1_all, kernel1);
        kernel2_all = cat(1, kernel2_all, kernel2);
        kernel_both_all = cat(1, kernel_both_all, kernel_both);
    end
    
    subplot(nTgtCST, n_ntrials, n_ntrials*(iTgtCST-1)+i_ntrials), hold on
    for n=1:3
        switch n 
            case 1, kk = kernel1_all; cc = '-r'; 
            case 2, kk = kernel2_all; cc = '-b'; 
            case 3, kk = kernel_both_all; cc = '-k'; 
        end
        errorbar(filterSF_all_log, mean(kk,1), std(kk, [],1)/sqrt(nsubj), cc)
    end
    
    plot(log([.8, 4.8]), [0,0], 'color', [.5, .5, .5])
    plot(log([2,2]), [-.2, 1], 'color',[.5,.5,.5])
    
    xticks(log([1,2,4])), xticklabels([1,2,4])
    xlim(log([.8,4.8]))
    ylim([-.2,1])
    title(sprintf('tgt cst = %.2f, %d trials', tgtCST_all(iTgtCST), ntrials*2))
    if ipair == 1, xlabel('SF channel (cpd)'), ylabel('kernel (a.u.)'), legend('tgt-prs trials', 'tgt-abs trials', 'all trials', 'location', 'best'), end
    
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
    sgtitle(sprintf('kernel by trial type n=%d nCST = %.2f', nsubj, noise.noiseCST), 'fontsize', 20)
end


