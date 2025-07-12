
[nORI, nSF] = size(template);

%% histogram of IV
if flag_plotIVsDist, figure('Position', [0 0 2000 500]), end
for itest = 1:nNormTests
    switch itest % transform PRS and ABS together
        case 1, IV_trans = IV2; % raw
        case 2, IV_trans = exp(IV2);
        case 3, IV_trans = sqrt(IV2);
        case 4,  IV_trans = sqrt(IV2);
        case 5, IV_trans = 1./IV2;
    end
    
    IV_PRS_trans = IV_trans(1:nPRS);
    IV_ABS_trans = IV_trans(nPRS+1: end); assert(length(IV_ABS_trans) == nABS)
    
    [hPRS, pPRS] = swtest(IV_PRS_trans'); if hPRS == 0, s_PRS = 'succeeded'; else, s_PRS = 'FAILED'; end
    [hABS, pABS] = swtest(IV_ABS_trans); if hABS == 0, s_ABS = 'succeeded'; else, s_ABS = 'FAILED'; end
    
    h_allB(ii, itest, :) = [hPRS, hABS];
    f1 = fitdist(IV_PRS_trans, 'Normal');
    f2 = fitdist(IV_ABS_trans, 'Normal');
    
    if flag_plotIVsDist
        % IV histogram
        subplot(2, nNormTests*2, [(itest-1)*2+1, itest*2]), hold on
        h1 = histfit(IV_PRS_trans, 20); h1(1).FaceColor = 'r'; h1(1).LineWidth = .5; h1(1).FaceAlpha = .5; h1(2).Color = 'r'; h1(2).LineWidth = 2;
        h2 = histfit(IV_ABS_trans, 20); h2(1).FaceColor = 'b'; h2(1).LineWidth = .5; h2(1).FaceAlpha = .5; h2(2).Color = 'b'; h2(2).LineWidth = 2;
        
        title(sprintf('%s\nPRS~N(%.2f, %.3f) \nABS~N(%.2f, %.3f)\nFitted PRS~N(%.2f, %.3f)\nFitted ABS~N(%.2f, %.3f)', ...
            namesNormalityTest{itest}, median(IV_PRS_trans), std(IV_PRS_trans), median(IV_ABS_trans), std(IV_ABS_trans), ...
            f1.mu,  f1.sigma, f2.mu,  f2.sigma))
        
        % qq plot [PRS]
        subplot(2, nNormTests*2, nNormTests*2 + (itest-1)*2 +1),
        h = qqplot(IV_PRS_trans); set(h(1),'markeredgecolor', 'r');
        xlabel(''), ylabel(''), axis square, title(sprintf('%s\np = %.3f', s_PRS, pPRS))
        if ~hPRS, box on, set(gca,'linewidth',3), end
        
        % qq plot [ABS]
        subplot(2, nNormTests*2, nNormTests*2 + (itest-1)*2 +2)
        h = qqplot(IV_ABS_trans); set(h(1),'markeredgecolor', 'b');
        xlabel(''), ylabel(''), axis square, title(sprintf('%s\np = %.3f', s_ABS, pABS))
        if ~hABS, box on, set(gca,'linewidth',3); end
        
    end
end
if flag_plotIVsDist
    sgtitle(sprintf('%s-%s\n%s', subjName, namesLocComb{iLocComb}, namesModelA{iModelA}))
end

%% plot the IV per channel
if flag_plotIVsDist
    figure('Position', [0 0 500 400])
    % ORI
    subplot(2,2,1)
    histogram(maxPRS_allT(:, 1), axis_tuning{1})
    xticks(axisTicks_tuning{1}), xticklabels(axisTL_tuning{1})
    xline(0, 'r--', 'linewidth', 2);
    title('[PRS] ORI'), xline(axis_tuning{1}(ORI_bound(1)), 'r-', 'linewidth', 2); xline(axis_tuning{1}(ORI_bound(2)), 'r-', 'linewidth', 2);
    
    subplot(2,2,3)
    histogram(squeeze(maxABS_allT(:, 1)), axis_tuning{1}),
    xline(0, 'r--', 'linewidth', 2);
    xticks(axisTicks_tuning{1}), xticklabels(axisTL_tuning{1})
    title('[ABS] ORI'), xline(axis_tuning{1}(ORI_bound(1)), 'r-', 'linewidth', 2); xline(axis_tuning{1}(ORI_bound(2)), 'r-', 'linewidth', 2);
    
    % SF
    subplot(2,2,2)
    histogram(squeeze(maxPRS_allT(:, 2)), axis_tuning{2})
    xticks(axisTicks_tuning{2}), xticklabels(axisTL_tuning{2})
    xline(1, 'r--', 'linewidth', 2);
    title('[PRS] SF')
    
    subplot(2,2,4)
    histogram(squeeze(maxABS_allT(:, 2)), axis_tuning{2})
    xline(1, 'r--', 'linewidth', 2);
    xticks(axisTicks_tuning{2}), xticklabels(axisTL_tuning{2})
    title('[ABS] SF')
    
    sgtitle(sprintf('%s-%s-%s\nAt which ORI/SF channel the max IV is found', subjName, namesLocComb{iLocComb}, namesModelA{iModelA}))
    
end