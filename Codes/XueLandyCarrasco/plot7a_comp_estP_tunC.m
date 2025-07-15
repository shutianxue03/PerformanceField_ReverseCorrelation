
% compare the estimated params (estP) and the corresponding tuning characterstics (tunC)

labels_estP_ORI = {'Gain', 'Width', 'Baseline'};
labels_estP_SF = {'Peak SF', 'Gain', 'Width', 'Baseline', 'Truncation'};
labels_tunC_ORI = {'PeakAmp', 'Bandwidth', 'Baseline'};
labels_tunC_SF = {'Peak SF', 'PeakAmp', 'Bandwidth', 'Baseline', 'Truncation'};

for ifeature = 1:2
    fprintf('%s\n', namesFeature{ifeature})
    clear estP tuningC
    if ifeature == 1
        estP = estP_ORI_allSubj;
        tuningC = tuningC_ORI_allSubj;
        labels_estP = labels_estP_ORI;
        labels_tunC = labels_tunC_ORI;
        limits = {[0, .25], [5, 40], [-.05, .05]};
    else
        estP = estP_SF_allSubj;
        tuningC = tuningC_SF_allSubj;
        % peakSF, width full, left, right
        %     margParams_c(:, :, :, :, [1,3,6,7]) = 2.^margParams_c(:, :, :, :, [1,3,6,7]);
        %     margParams_c(:, :, :, :, 3) = margParams_c(:, :, :, :, 7);
        labels_estP = labels_estP_SF;
        labels_tunC = labels_tunC_SF;
        limits = {[1.5, 2.5], [0, .1], [0, 2], [-.05, .05], [-.01, .02]};
    end
    
    nparams = size(estP, 5);
    nsubj = size(estP, 1);
    
    %% 1. compare gain/peak, sigma/width or baseline/bottom
    figTitle = sprintf('comp_all%d', nLoc);
    
    IV_loc = nan(nsubj, nLoc, 2); IV_cond = IV_loc;
    for isubj = 1:nsubj
        IV_loc(isubj, :, :) = repmat((1:nLoc)', 1, 2);
        IV_cond(isubj, :, :) = repmat([1,2], nLoc, 1);
    end
    
    figure('Position', [3e3, 0, 1200, 800])
    for ip = 1:nparams
        subplot(2,3,ip), hold on, box on
        plot([limits{ip}(1), limits{ip}(2)], [limits{ip}(1), limits{ip}(2)], 'k-')
        
        % estimated params
        [a1, ~, ~, a1_neg, a1_pos] = getCI(estP(:, :, :, itype, ip), 1, 2);
        % tuning characteristics
        [a2, ~, ~, a2_neg, a2_pos] = getCI(tuningC(:, :, :, itype, ip), 1, 2);
        % plot
        
        for iLoc = 1:nLoc
            %         errorbar(a1(:, iLoc), a2(:, iLoc), a1_neg(:, iLoc), a1_pos(:, iLoc), 'horizontal', 'color', colors_comb_(iLoc, :), 'Capsize', 0)
            %         errorbar(a1(:, iLoc), a2(:, iLoc), a2_sem(:, iLoc), a2_pos(:, iLoc), 'vertical', 'color', colors_comb_(iLoc, :), 'Capsize', 0)
            for isubj = 1:nsubj
                plot(a1(isubj, iLoc), a2(isubj, iLoc), markers_allSubj{isubj}, 'markerfacecolor', 'w', 'markeredgecolor', colors_comb_(iLoc, :))
            end
        end
        xlabel(labels_estP{ip})
        ylabel(labels_tunC{ip})
        xlim(limits{ip})
        ylim(limits{ip})
        axis square
        % anova and ttest
        data = cat(3, a1, a2);
        anova_text = print_nANOVA({'Loc', 'Cond'}, data(:), {IV_loc(:), IV_cond(:)});
        title(anova_text)
    end
    set(findall(gcf, '-property', 'fontsize'), 'fontsize', 13)
    
    % save
    folderName = sprintf('VSS2023/fig/RC/%s/comp_estP_tunC/', nameEnergySource);
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    saveas(gcf, sprintf('%s%s_all%d_n%d.jpg', folderName, namesFeature{ifeature}, nLoc, nsubj))
    
end

