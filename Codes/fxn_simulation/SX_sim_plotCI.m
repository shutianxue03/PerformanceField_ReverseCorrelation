
nsubj = simPairs(1,ipair);
names = {'mean', 'var'};
for n = 1%:2 % CI mean and var
    if n==1, CI = CImean_allSubj; else, CI = CIvar_allSubj;end
    
    figure('Position', [0 200 n_ntrials*300 nTgtCST*300 ])
    for ipair = 1:nSimPairs
        
        ntrials = simPairs(2,ipair);
        i_ntrials = find(ntrials == ntrials_all);
        ntrialsAll = ntrials*2;
        
        tgtCST = simPairs(3, ipair);
        iTgtCST = find(tgtCST == tgtCST_all);
        
        subplot(nTgtCST, n_ntrials, n_ntrials*(iTgtCST-1)+i_ntrials), hold on
        imshow(mean(cat(3, CI{:, i_ntrials, iTgtCST}),3) + .1)
        title(sprintf('tgt cst = %.2f, %d trials', tgtCST, ntrials*2))
        
    end    
    
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
    sgtitle(sprintf('CI - %s nCST = %.2f', names{n}, noise.noiseCST), 'fontsize', 20)
end