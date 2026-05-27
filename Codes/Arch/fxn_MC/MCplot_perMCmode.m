

indCands = 1:nCands;

figure('Position', [0,0,nn*800, 600]), hold on

for iIC = 1:nIC_
    iBest_unik = unique(iBest_allSubj(:, iIC));
    [freq, iBest_best] = groupcounts(iBest_allSubj(:, iIC));
    nChosen = length(freq);
    
    dev_groupAVE = squeeze(mean(dev_allSubj(:, :, iIC)));
    dev_groupSEM = squeeze(std(dev_allSubj(:, :, iIC)))/sqrt(nsubj);
    [~, irank_groupAVE] = sort(dev_groupAVE); % sort based on averaged dev (for CV) for IC
    iBEST_groupAVE = irank_groupAVE(1);
    dd_ = dev_groupAVE-min(dev_groupAVE);
    
    %% save the dev for comparing across MCmode and families
    if MCmode == 3
        dev_bestGroup_perModePerFamily{ifamily, MCmode, iIC} = dev_allSubj(:, iBEST_groupAVE, iIC);
        ibestGroup_perModePerFamily{ifamily, MCmode, iIC} = paramInd_all(iBEST_groupAVE, :);
    else
        dev_bestGroup_perModePerFamily{ifamily, MCmode, 1} = dev_allSubj(:, iBEST_groupAVE, 1);
        ibestGroup_perModePerFamily{ifamily, MCmode, 1} = paramInd_all(iBEST_groupAVE, :);
    end
    
    %% DEV
    if MCmode == 3
        subplot(2, nIC, iIC),
        title(sprintf('%s\n Model w/ lowest DEV/IC: M#%d [%s]\nDEV/IC = %.5f\n Model w/ highest freq: M#%d [%s]\nDEV/IC = %.5f', ...
            namesIC{iIC}, ...
            iBEST_groupAVE, num2str(paramInd_all(iBEST_groupAVE, :)), mean(dev_allSubj(:, iBEST_groupAVE, iIC)), ...
            iBest_best(end), num2str(paramInd_all(iBest_best(end), :)), mean(dev_allSubj(:, iBest_best(end), iIC))))
    else, subplot(2, 1, 1), title(namesYLabel{MCmode})
    end
    hold on
    
    bar(indCands, dd_(irank_groupAVE))
    errorbar(indCands, dd_(irank_groupAVE), dev_groupSEM, '.k', 'CapSize', 0)
    for isubj = 1:nsubj
        plot(indCands, squeeze(dev_allSubj(isubj, irank_groupAVE, iIC)) -min(dev_groupAVE) , '.-', 'color', ones(1,3)*.6)
    end
    xticks(indCands)
    xticklabels(indCands(irank_groupAVE))
    xlim([.5, nmodels+.5])
    ylabel(namesYLabel{MCmode})
    
    %% FREQ
    if MCmode == 3, subplot(2, nIC, iIC + nIC)
    else, subplot(2, 1, 2)
    end
    hold on
    for ichosen = 1:nChosen
        xPos = find(iBest_best(ichosen) == irank_groupAVE);
        bar(xPos, freq(ichosen), 'FaceColor', 'w')
    end
    xticks(indCands)
    xTL = cell(nmodels, 1); for ii = 1:nmodels, xTL{ii} = sprintf('%d [%s]', irank_groupAVE(ii), num2str(paramInd_all(irank_groupAVE(ii), :))); end
    xticklabels(xTL), xtickangle(-50)
    xlim([.5, nmodels+.5])
    ylim([0, nsubj/2])
    xlabel('Candidate model #')
    if iIC == 1, ylabel('Freq of being the best'), end
    
end % iIC

set(findall(gcf, '-property', 'FontSize'), 'FontSize',12)

if MCmode == 3
    sgtitle(sprintf('n = %d\n%s - Family #%d %s\n%s', ...
        nsubj, namesFeature{ifeature}, ifamily, namesFamily_all{ifamily}, namesMCmode{MCmode}))
else
    sgtitle(sprintf('n = %d - %s\n[%s - Family #%d %s]\n Model w/ lowest DEV/IC: M#%d [%s]\nDEV/IC = %.5f\n Model w/ highest freq: M#%d [%s]\nDEV/IC = %.5f', ...
        nsubj, namesMCmode{MCmode}, namesFeature{ifeature}, ifamily, namesFamily_all{ifamily}, ...
        iBEST_groupAVE, num2str(paramInd_all(iBEST_groupAVE, :)), mean(dev_allSubj(:, iBEST_groupAVE, 1)), ...
        iBest_best(end), num2str(paramInd_all(iBest_best(end), :)), mean(dev_allSubj(:, iBest_best(end), 1))))
end
