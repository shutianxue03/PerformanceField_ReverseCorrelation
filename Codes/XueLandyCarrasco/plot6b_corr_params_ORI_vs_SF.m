
% assess the corr between corresponding estP/tunC between ORI and SF

labels_estP_ORI = {'Gain', 'Width', 'Baseline'};
labels_estP_SF = {'Peak SF', 'Gain', 'Width', 'Baseline', 'Truncation'};
labels_tunC_ORI = {'PeakAmp', 'Bandwidth', 'Baseline'};
labels_tunC_SF = {'Peak SF', 'PeakAmp', 'Bandwidth', 'Baseline', 'Truncation'};

% assess the corr between corresponding paramsters/tuningC for ORI and SF
close all

for paramMode = 1:2
    if paramMode == 1 % Estimated params
        pORI = estP_ORI_allSubj; pSF = estP_SF_allSubj;
        xlabels = labels_estP_ORI; ylabels = labels_estP_SF; sgtitle_ = 'Estimated params'; figTitle = 'estP';
    else % Tuning Cs
        pORI = tuningC_ORI_allSubj; pSF = tuningC_SF_allSubj;
        xlabels = labels_tunC_ORI; ylabels = labels_tunC_SF; sgtitle_ = 'Tuning Cs'; figTitle = 'TunC';
    end
    
    indLoc_corr = repmat(1:nLoc, nsubj, 1);
    indSubj_corr = repmat((1:nsubj)', 1, nLoc);
    
    iplots = 1:nparams*nparams; iplots = reshape(iplots, [nparams, nparams])';
    
    for ip = 1:3 % gain, width, baseline
        % ORI
        [x, ~, ~, x_neg, x_pos] = getCI(pORI(:, :, :, itype, ip), 1, 2);
        % SF
        [y, ~, ~, y_neg, y_pos] = getCI(pSF(:, :, :, itype, ip+1), 1, 2);
        x_ = x-mean(x);
        y_ = y-mean(y);
        
        % partial corr
        [r, p] = partialcorr([x(:), y(:)], [indLoc_corr(:), indSubj_corr(:)]);
        r_partial = r(2,1); p_partial = p(2,1);
        titles = sprintf('Partial r = %.2f, p = %.3f\n', r_partial, p_partial);
        text_sig = '';
        if p_partial<.05, if r_partial>0, text_sig = 'pos'; else, text_sig = 'neg'; end, end
        
        % plot
        figure('Position', [3e3, 0, 500, 500]), hold on
        nsig = 0;
        for iLoc = 1:nLoc
            % idvd dots
            errorbar(x_(:, iLoc), y_(:, iLoc), x_neg(:, iLoc), x_pos(:, iLoc), 'horizontal', '.','color', colors_comb_(iLoc, :), 'Capsize', 0)
            errorbar(x_(:, iLoc), y_(:, iLoc), y_neg(:, iLoc), y_pos(:, iLoc), 'vertical', '.','color', colors_comb_(iLoc, :), 'Capsize', 0)
            for isubj = 1:nsubj
                plot(x_(isubj, iLoc), y_(isubj, iLoc), markers_allSubj{isubj}, 'markerfacecolor', 'w', 'markeredgecolor', colors_comb_(iLoc, :))
            end
            % corr per loc
            [r_, p_] = corr(x_(:, iLoc), y_(:, iLoc));
            if p_<.05
                titles = [titles, sprintf('%s%s', namesLocComb_{iLoc}, getString_starts(p_))];
                nsig = nsig+1;
            end
        end % iLoc
        
        % linear regression
        lm = polyfit(x_(:), y_(:), 1);
        x_lm2 = linspace(min(x_(:)) - std(x_(:)), max(x_(:)) + std(x_(:)), 2);
        yfit = polyval(lm, x_lm2);
        plot(x_lm2, yfit, 'k', 'linewidth', 2, 'handlevisibility', 'off');
        
        %
        xlabel(sprintf('Zero-meaned ORI P%d [%s]', ip, xlabels{ip}))
        ylabel(sprintf('Zero-meaned SF P%d [%s]', ip+1, ylabels{ip+1}))
        axis square
        title(sprintf('[%s] %s\n%s', figTitle, xlabels{ip}, titles))
        set(findall(gcf, '-property', 'fontsize'), 'fontsize', 18)
        set(findall(gcf, '-property', 'linewidth'), 'linewidth', 1.2)
        
        folderName = sprintf('VSS2023/fig/RC/%s/corr_params_ORI_vs_SF/all%d/', nameEnergySource, nLoc);
        folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
        saveas(gcf, sprintf('%sn%d_%s_%s_%s%d.jpg', folderName, nsubj, figTitle, xlabels{ip}, text_sig, nsig))
        
    end % ip
    
end % ii
% close all