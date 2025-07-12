

for iIC = 1:nn
    iBest_unik = unique(iBest_allSubj(:, iIC));
    [freq, iBest_best] = groupcounts(iBest_allSubj(:, iIC));
    nChosen = length(freq);
    
    dev_groupAVE = squeeze(mean(dev_allSubj(:, :, iIC)));
    dev_groupSEM = squeeze(std(dev_allSubj(:, :, iIC)))/sqrt(nsubj);
    [~, irank_groupAVE] = sort(dev_groupAVE); % sort based on averaged dev (for CV) for IC
    iBest_group = irank_groupAVE(1);
    dd_ = dev_groupAVE-min(dev_groupAVE);
    
    figure('Position', [3e3,0, 2e3, 2e3]), hold on
    for isubj = 1:nsubj
        subjName = subjList{isubj};
        nblocks = nblocks_allSubj(isubj);
        nameKernel = sprintf('Data_OOD/%s%d/%s_kernel_%d_%d.mat', subjName, nblocks, subjName, nORI, nSF);
        load(nameKernel, 'margORI_perComb', 'margSF_perComb')
        
        iBest_idvd = iBest_allSubj(isubj);
        
        subplot(4, 3, isubj), hold on
        if iBest_idvd == iBest_group, box on, set(gca,'linewidth',2); else, box off; end
        switch ifeature, case 1, marg = margORI_perComb; axis_tuning_ln = axis_tuning{ifeature}; case 2, marg = margSF_perComb; axis_tuning_ln = 2.^axis_tuning{ifeature}; end
        for iline = 1:nLoc2
            iLocComb = iLocComb_all(iline);
            % data
            plot(axis_tuning{ifeature}, marg{iLocComb}(itype, :), 'o', 'color', colors2{iline})
            
            % pred (group best model)
            params_groupBest = params_est_allSubj{isubj, iBest_group};
            ypred_groupBest = MC_predKernel(paramInd_all(iBest_group, :), axis_tuning_ln, 2, nparams_full, params_groupBest, ifamily);
            plot(axis_tuning{ifeature}, ypred_groupBest(iline, :), '-', 'color', colors2{iline}, 'linewidth', 2)

            % pred (idvd best model)
            params_idvdBest = params_est_allSubj{isubj, iBest_idvd};
            ypred_groupBest = MC_predKernel(paramInd_all(iBest_idvd, :), axis_tuning_ln, 2, nparams_full, params_idvdBest, ifamily);
            plot(axis_tuning{ifeature}, ypred_groupBest(iline, :), 'k--', 'linewidth', 1)
        end % iline
        
        if isubj==1, legend({'Data', 'Pred (group best)', 'Pred (idvd best)'}, 'Location', 'best'), end
        title(sprintf('%s\nIDVD best M%d [%s], dev/IC = %.5f\nGroup best: dev/IC = %.5f', ...
            subjName,...
            iBest_idvd, num2str(paramInd_all(iBest_idvd, :)), dev_allSubj(isubj, iBest_idvd, iIC), ...
            dev_allSubj(isubj, iBest_group, iIC)))
    end % isubj
    if MCmode == 3
        sgtitle(sprintf('n = %d\n%s - Family #%d %s\n%s-%s\nBest group model M%d [%s]', ...
        nsubj, namesFeature{ifeature}, ifamily, namesFamily_all{ifamily}, namesMCmode{MCmode}, namesIC{iIC}, iBest_group, num2str(paramInd_all(iBest_group, :))))
    else
    sgtitle(sprintf('n = %d\n%s - Family #%d %s\n%s\nBest group model M%d [%s]', ...
        nsubj, namesFeature{ifeature}, ifamily, namesFamily_all{ifamily}, namesMCmode{MCmode}, iBest_group, num2str(paramInd_all(iBest_group, :))))
    end
end % iIC