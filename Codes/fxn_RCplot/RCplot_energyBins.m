

%%
% plot rHit/rFA vs.energy bins, at Gabor ori and SF (both when n=1 and n>1)

%%
ylabels = {'Hit rate', 'FA rate'};
varNames = {'Loc', 'Energy'};
if nsubj == 1, N = nB; errorType = 1;
else, N=nsubj; errorType = 2;
end
ind_Loc8 = repmat(1:nLoc8, [N, 1, nbins_e]);
ind_ebin = [];
for ib = 1:nbins_e
    ind_ebin_ = ones(N, nLoc8)*ib;
    ind_ebin = cat(3, ind_ebin, ind_ebin_);
end

%%
% for errorFlag = [1,0]
figure('Position', [0 300 ntypes*400 400])


for itype = 1:ntypes
    [r_ave, r_lb, r_ub, r_SEM_neg, r_SEM_pos] = getCI(r(:, itype, :, :), errorType, 1); 
    [e_ave, e_lb, e_ub, e_SEM_neg, e_SEM_pos] = getCI(e(:, itype, :, :), errorType, 1);
    
    slope_legends = cell(1, nLoc8);
    
    subplot(1,ntypes, itype), hold on, grid on
    
    % plot raw data and fitted line of each loc
    for iLoc = 1:nLoc8
        x = e_ave(iLoc, :) * 100;
        y = r_ave(iLoc, :);
        
        % raw data
        errorbar(x, y, r_SEM_neg(iLoc, :), r_SEM_pos(iLoc, :), e_SEM_neg(iLoc, :), e_SEM_pos(iLoc, :), '.', 'color', colors_comb(iLoc, :), 'linewidth', .05, 'CapSize', 0, 'handlevisibility', 'off')
        plot(x, y, 'o', 'MarkerFaceColor',  colors_comb(iLoc, :),  'MarkerEdgeColor', 'w')
        
        % fitting
        [pearsonr,p] = corr(x',y');
        ss = getString_starts(p);
        [beta, stats] = polyfit(x/100, y,1);
        slope = beta(1);
        yfit = polyval(beta, x/100);
        plot(x,yfit, '.-', 'color', colors_comb(iLoc, :), 'handlevisibility', 'off')
        
        % slope text
        yresid = y - yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(y)-1) * var(y);
        R2 = 1 - SSresid/SStotal;
        slope_legends{iLoc} = sprintf('[%s] r=%.2f%s', namesLocComb{iLoc}, pearsonr, ss);
        
    end
    legend(slope_legends, 'Location', 'best', 'NumColumns',2)
    
    % ANOVA2 (3 bins x 5 loc on rFA and rHit separately)
    if (nsubj > 1) || (nB > 1) % so that each dot stands for multiple data points
        r_anova = squeeze(r(:, itype, :, :)); % nB/nsubj x nLoc8 x nbins
        [p, tbl] = anovan(r_anova(:), {ind_Loc8(:), ind_ebin(:)},'model','full','varnames', varNames, 'display', 'off');
        anova_text = print_nANOVA(varNames, r_anova(:),  {ind_Loc8(:), ind_ebin(:)});
        title(anova_text)
    end
    ylim([0 1]), yticks(0:.25:1)
    xlabel('energy (%)')
    ylabel(namesRespType{itype})
    
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
make_sgtitle('Binned energy', subjName, nB, nsubj, nAllTrials)



