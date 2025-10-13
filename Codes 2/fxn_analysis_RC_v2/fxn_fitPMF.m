
% Input: ccc_full [4 columns]: 1=Loc, 2=noise level, 3=cst, 4=correctness

cst_log_full = log10(ccc_full(:, 3));
corr_full = ccc_full(:, 4);
curveX_log = log10(fit.curveX);

%% get unique contrast levels and corresponding pC
if flag_binData
    nBins = fit.nBins;
    [N,edges,indBin] = histcounts(cst_log_full, nBins);
    cst_log_unik = []; nCorr = []; nData = [];
    for iBin = 1:nBins
        if N(iBin)>0
            cst_log_unik = [cst_log_unik, mean(cst_log_full(indBin==iBin))];
            nCorr = [nCorr, nansum(corr_full(indBin==iBin))];
            nData = [nData, nansum(indBin==iBin)]; assert(nansum(indBin==iBin) == N(iBin))
        end
    end
    assert(length(cst_log_unik) == sum(N>0))
else
    cst_log_unik = unique(cst_log_full);
    nCorr = nan(length(cst_log_unik), 1); nData=nCorr;
    for icst_unik = 1:length(cst_log_unik)
        indCST = cst_log_full==cst_log_unik(icst_unik);
        nCorr(icst_unik) = nansum(corr_full(indCST));
        nData(icst_unik) = nansum(indCST==1);
    end % icst_unik
    cst_log_unik = cst_log_unik';
    nCorr = nCorr';
    nData = nData';
end

pC = nCorr./nData;

%%
if flag_filterData
    % remove data with pC <.4 and when cst<2% pC>.65
    ind_noise = (pC<=.4) | (pC>.65 & cst_log_unik<log10(.01));
    ind=boolean(1-ind_noise);
    cst_log_unik = cst_log_unik(ind);
    nCorr = nCorr(ind);
    nData = nData(ind);
    pC = pC(ind);
    fprintf('[Removed ~%.0f%% trials] ...', 100*mean(1-ind))
    
    cst_log_unik = [-3, -2.5, cst_log_unik];
    nCorr = [20 20 nCorr];
    nData = [40 40 nData];
    pC = [.5, .5, pC];
    
end

%% fit PF
for iModel = 1:nModels
    fprintf('%s...', PMF_models{iModel})
    switch iModel
        case 1, fit.PF = @PAL_Logistic;
        case 2, fit.PF = @PAL_CumulativeNormal;
        case 3, fit.PF = @PAL_Gumbel;
        case 4, fit.PF = @PAL_Weibull;
    end
    
    % adjust the format of cst (linear/log)
    searchGrid = fit.searchGrid;
    if iModel==4, cst_fit = 10.^cst_log_unik; cst_fineG = 10.^curveX_log;  searchGrid.alpha = 10.^fit.searchGrid.alpha; % weibull fit to linear values
    else, cst_fit = cst_log_unik; cst_fineG = curveX_log; searchGrid.alpha = fit.searchGrid.alpha;
    end
    
    % estimate parameters
    [paramVals] = PAL_PFML_Fit(cst_fit, nCorr, nData, ...
        searchGrid, fit.paramsFree, fit.PF,'SearchOptions',fit.options,'LapseLimits',fit.lapseLimits,'GuessLimits',fit.guessLimits);
    
    % conduct bootstrapping
    [~, paramVals, LL, converged] = PAL_PFML_BootstrapParametric(cst_fit, nData, ...
        paramVals, fit.paramsFree, fit.nBoot, fit.PF,'LapseLimits',fit.lapseLimits,'GuessLimits',fit.guessLimits);
    
    fprintf(' (%.0f%%) ', mean(converged == 1)*100)
    
    % estimate thresh
    for iB = 1:fit.nBoot
        % save predicted pC
        if isnan(paramVals(iB, 2)), paramVals(iB, 2) = Inf; end
        yfit_allB(iB, iModel, :) = fit.PF(paramVals(iB, :), cst_fineG);
        
        % save estimated thresh
        for iPerf = 1:nPerf
            perfThresh = perfThresh_all(iPerf)/100;
            
            thresh = fit.PF(paramVals(iB, :), perfThresh, 'Inverse');
            if thresh == -Inf, thresh=nan; end
            if iModel==4, if thresh<0, error('ALERT: thresh is negative cannot take log10!!'); end, thresh=log10(thresh); end
            thresh_allB(iB, iModel, iPerf) = thresh;
        end
    end % iPerf
    
    LL_allB(:, iModel) = LL;
%     converged_allB(:, iModel) = converged;
%     slope_allB(:, iModel) = paramVals(:, 2) * ((1-.5)./sqrt(2.*pi)); % according to Strasburger 2001 eq.18
%     guess_allB(:, iModel) = paramVals(:, 3);
%     lapse_allB(:, iModel) = paramVals(:, 4);
%     estP_allB(:, iModel, :) = paramVals;
    
end % im=1:nModels
