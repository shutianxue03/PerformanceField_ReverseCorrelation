% estimated parameters
params_allSubj = [];

for iLoc = [1, 8, 6, 7, 5, 3]
    for isubj=1:nsubj
        
        subjName = subjList{isubj};
        nameFolder_NOM_trialWise = sprintf('%s/ORI%dSF%d/%s/L%d', nameFolder_NOM0, nORI, nSF, subjName, iLoc);
        nameFileModelIDVD_trialWise = sprintf('%s/n%d_A%dB%d', nameFolder_NOM_trialWise, ni, iModelA, iModelB);
        load(nameFileModelIDVD_trialWise, 'params_est_allB')
        
        params_allSubj(iLoc, isubj, :, :) = params_est_allB;
        nParams = size(params_est_allB,2);
    end
end

%%
for iLoc_ = {[1,8], [6,7], [5,3]}
    iLoc = iLoc_{1};
    figure('Position', [0 0 1e3 300])
    for iParam = 1:nParams
        params_allSubj_med = getCI(params_allSubj(iLoc, :, :, iParam), 1, 3);
        [params_ave, ~, ~, params_sem] = getCI(params_allSubj_med, 2, 2);

        subplot(1,4,iParam), hold on
        for iiLoc=1:2
            bar(iiLoc, params_ave(iiLoc), 'FaceColor', 'w', 'EdgeColor', colors_comb(iLoc(iiLoc), :))
            errorbar(iiLoc, params_ave(iiLoc), params_sem(iiLoc), 'CapSize', 0, 'color', colors_comb(iLoc(iiLoc), :))
        end
        for isubj=1:nsubj
            plot([1.2,1.8], params_allSubj_med(:, isubj), '-', 'color', ones(1,3)/2)
        end
        
        t_allB=[];
        p_allB=t_allB;
        for ii=1:ni
            a=params_allSubj(iLoc(1), :, ii, iParam);
            b=params_allSubj(iLoc(2), :, ii, iParam);
            [~, p, ~, stats] = ttest(a,b);
            t_allB(ii) = stats.tstat;
            p_allB(ii) = p;
        end
        
        title(sprintf('%s\nt=%.2f, p=%.3f', namesParamsModel_all{iModelB}{iParam}, median(t_allB), median(p_allB)))
    end
    sgtitle(sprintf('L%d%d', iLoc))
    set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
end
