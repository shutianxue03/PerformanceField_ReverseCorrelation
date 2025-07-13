

iLoc_all = 1:5;
nLoc = length(iLoc_all);
nm=3;

for isubj = 1:nsubj
    
    subjName = subjList{isubj};
    nblock = nblocks_allSubj(isubj);
    
    load(sprintf('%s/Data_OOD_%d%d/%s%d/%s_behavMeas.mat', nameFolder_Data, nORI, nSF, subjName, nblock, subjName), ...
        'cst_perSess_perLoc', 'pA3_perSess_perLoc', 'pC3_perSess_perLoc')
    
    %%
    figure('Position', [0 0 1e3 1e3])
    
    for im = 1:nm
        switch im
            case 1, m= cst_perSess_perLoc*100; y_label = 'Gabor CST (%)'; y_ticks = linspace(log10(.3), log10(.8), 5); y_ticklabels = round(10.^y_ticks*100);
                y_ticks = linspace(30,80, 5); y_ticklabels = y_ticks;
            case 2, m=squeeze(pA3_perSess_perLoc(:, :, 1))*100; y_label = 'pA (%)'; y_ticks = 40:15:100; y_ticklabels = y_ticks;
            case 3, m=squeeze(pC3_perSess_perLoc(:, :, 1))*100; y_label = 'pC (%)'; y_ticks = 50:10:90; y_ticklabels = y_ticks;
        end
        subplot(nm,1,im), hold on, grid on
        
        for iiLoc = 1:nLoc
            plot(m(:, iLoc_all(iiLoc)), 'color',colors_comb(iLoc_all(iiLoc), :)),
        end
        
        if im==2, yline(60, 'k'); yline(80, 'k'); end
        if im==3, yline(70, 'k'); end
        
        ylabel(y_label)
        yticks(y_ticks)
        yticklabels(y_ticklabels)
        ylim(y_ticks([1,end]))
        
        for iiLoc=1:nLoc
            iLoc = iLoc_all(iiLoc);
            meas_sd(isubj, im, iiLoc) = std(m(:, iLoc));
            meas_min(isubj, im, iiLoc) = min(m(:, iLoc));
            meas_max(isubj, im, iiLoc) = max(m(:, iLoc));
            meas_ave(isubj, im, iiLoc) = mean(m(:, iLoc));
        end
        
    end % im
    
    xlabel('Session #'),
    xlim([.5, nSess + .5])
    xticks(1:2:nSess)
    %     legend(namesLoc2D, 'location', 'best')
    
    sgtitle(sprintf('%s%d', subjName, nblock))
    
    set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
    set(findall(gcf, '-property', 'fontsize'), 'fontsize',15)
    
end % isubj

%%
clc

for im=3%:nm
    %     for iiLoc = 1:nLoc
    for in = 1:4
        switch in
            case 1, m = meas_sd; n='SD';
            case 2, m = meas_min; n='MIN';
            case 3, m = meas_max; n='MAX';
            case 4, m = meas_ave; n='AVE';
        end
        
        a=squeeze(m(:, im, :));
        indLoc = repmat(iLoc_all, nsubj, 1);
        text_ANOVA = print_nANOVA({'Loc'}, a(:), {indLoc(:)}, nsubj);
        
        [ave, ~, ~, SEM] = getCI(a(:), 2, 1);
        fprintf('%s: %.1f +- %.1f\n%s\n', n, ave, SEM, text_ANOVA)
    end
end
% CST: [29.7, 52.1%]