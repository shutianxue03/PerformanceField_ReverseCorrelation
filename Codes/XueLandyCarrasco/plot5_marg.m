clc
varNames = {'Loc', 'Channel'};

wd_border = 3;
sz_ticks = 20; % 30
sz_text = 30;
sz_stars = 20;

shift = [0,0];
for ifeature = 1:2
    
    if ifeature==1
        switch itype
            case 1, ymax = .1; ymin = -.05;
            case 2, ymax = .18; ymin = -.06;
            case 3, ymax = .1; ymin = -.05;
        end
    else
        switch itype
            case 1, ymax = .1; ymin = -.05;
            case 2, ymax = .1; ymin = -.05;
            case 3, ymax = .1; ymin = -.05;
        end
    end
    yrange = ymax - ymin;
    ytext_R2_all = {[ymax-yrange/12.5, ymax-yrange/6], [ymax-yrange/12.5, ymax-yrange/6]};
    ystars  = [ymin+yrange/4, ymin+yrange/4];
    
    
    xaxis = axis_tuning{ifeature};
    nfilters = length(xaxis);
    nameModelParams = namesParams_unit_perF{ifeature};
    nparams = length(nameModelParams);
    
    % for ANOVA
    ind_Loc = repmat(1:nLoc, nsubj, 1, nfilters);
    ind_channel = nan(nsubj, nLoc, nfilters);
    for isubj = 1:nsubj, ind_channel(isubj, :, :) = repmat(1:nfilters, nLoc, 1); end
    
    if ifeature == 1
        [marg_med, ~, ~, marg_neg, marg_pos] = getCI(margORI_allSubj_(:, :, :, itype, :), 1, 2);
        [margPred_med, margPred_CI_lb, margPred_CI_ub] = getCI(margPredORI_allSubj_(:, :, :, itype, :), 1, 2);
        [margR2_med, margR2_neg, margR2_pos] = getCI(margR2ORI_allSubj_(:, :, :, itype), 1, 2);
        
        %         bottom = marg_med(:, :, axis_tuning{ifeature} < -45 | axis_tuning{ifeature} > 45);
    else
        [marg_med, ~, ~, marg_neg, marg_pos] = getCI(margSF_allSubj_(:, :, :, itype, :), 1, 2);
        [margPred_med, margPred_CI_lb, margPred_CI_ub] = getCI(margPredSF_allSubj_(:, :, :, itype, :), 1, 2);
        [margR2_med, ~, ~, margR2_neg, margR2_pos]  = getCI(margR2SF_allSubj_(:, :, :, itype), 1, 2);
        
        %         bottom = marg_med(:, :, axis_tuning{ifeature} > 1.5);
    end
    
    %% whether the flanks are negative
    %     figure('Position', [1000 0 600 500]), hold on
    %     yyaxis left
    %     bottom_ = mean(bottom, 3);
    %     for iLoc = 1:nLoc
    %         [~, p, ~, stats] = ttest(bottom_(:));
    %         [bottom_ave , ~, ~, sem] = getCI(bottom_(:, iLoc), 2, 1);
    %         bar(iLoc, bottom_ave, 'edgecolor', colors_comb_(iLoc, :), 'facecolor', 'w', 'LineWidth', 1.5)
    %         errorbar(iLoc, bottom_ave, sem, 'color', colors_comb_(iLoc, :), 'CapSize', 0, 'LineWidth', 1.5)
    %         %         text(iLoc, bottom_ave+sem*1.5, getString_starts(p), 'HorizontalAlignment', 'center')
    %         text(iLoc, bottom_ave+sem*1.3, sprintf('%d%%', round(100*mean(bottom_(:,iLoc)<0))), 'HorizontalAlignment', 'center')
    %     end
    %
    %     %     if nLoc==2
    %     for isubj = 1:nsubj
    %         if nLoc == 2, plot([1.2, 1.8], bottom_(isubj, :), ['k-', markers_allSubj{isubj}], 'LineWidth', .5, 'color', ones(1,3)*.8)
    %         else, plot(1:nLoc, bottom_(isubj, :), ['k-', markers_allSubj{isubj}], 'LineWidth', .5, 'color', ones(1,3)*.8)
    %         end
    %     end
    %     %     end
    %
    %     xticks(1:nLoc)
    %     xticklabels(namesLocComb_)
    %     %     ylabel('Kernels (averaged across channels)', 'fontsize', sz_ticks)
    %
    %     ind1_Loc = repmat(1:nLoc, nsubj, 1, size(bottom, 3));
    %     ind2_channel = nan(nsubj, nLoc, size(bottom, 3));
    %     for isubj = 1:nsubj, ind2_channel(isubj, :, :) = repmat(1:size(bottom, 3), nLoc, 1); end
    %
    %     text_anova = print_nANOVA({'Loc', 'channel'}, bottom(:), {ind1_Loc(:), ind2_channel(:)}, nsubj);
    %     if ifeature==1
    %         title(sprintf('Whether the flanks are negative\nORI (<-45 and >45)\n\n%s', text_anova))
    %     else
    %         title(sprintf('Whether the flanks are negative\nSF (>2.83 cpd)\n\n%s', text_anova))
    %     end
    %     set(findall(gcf, '-property', 'fontsize'), 'fontsize',15)
    
    %% normalize margs (for each subj, divide both loc by the max of fovea)
    %     marg_med_norm = marg_med;
    %     margPred_med_norm = margPred_med;
    %     for isubj = 1:nsubj
    %     marg_med_norm(isubj, :, :) = marg_med(isubj, :, :)./max(marg_med(isubj, 1, :));
    %     margPred_med_norm(isubj, :, :) = margPred_med(isubj, :, :)./max(marg_med(isubj, 1, :));
    %     end
    %
    %     marg_med = marg_med_norm;
    %     margPred_med = margPred_med_norm;
    
    %% plot kernels & tuning fxn
    figure('Position', [0 0 1050 550])
    hold on, box on
    
    
    peak{ifeature} = nan(2, nsubj);
    
    %%
    if ifeature==1, x_band = {13:19}; else, x_band = {9:19, 24:26}; end
    nband = length(x_band);
    for iband=1:nband
        x_band_ = xaxis_(x_band{iband});
        nF = length(x_band_);
        patch([x_band_, flip(x_band_)], [ones(1, nF)*ymin, ones(1, nF)*ymax], 'k', 'FaceAlpha', .1, 'linestyle', 'none')
    end
    
    for iLoc = 1:nLoc
        
        
        
        
        %%
        color = colors_comb_(iLoc, :);
        
        if nsubj>1
            [marg_ave, ~, ~, marg_SEM] = getCI(marg_med(:, iLoc, :), 2, 1);
            [margPred_ave, margPred_lb, margPred_ub, ~, ~] = getCI(margPred_med(:, iLoc, :), 2, 1);
            [margR2_ave, ~, ~, margR2_SEM, ~] = getCI(margR2_med(:, iLoc), 2, 1);
        else
            marg_ave = marg_med(iLoc, :);
            marg_SEM = (marg_neg(iLoc, :) + marg_pos(iLoc, :))/2;
            margPred_ave = margPred_med(iLoc, :);
            margPred_lb = margPred_CI_lb(iLoc, :);
            margPred_ub = margPred_CI_ub(iLoc, :);
            margR2_ave = margR2_med(iLoc);
            margR2_SEM = (margR2_neg(iLoc, :) + margR2_pos(iLoc, :))/2;
        end
        % get the peak
        %         [~, maxInd] = max(squeeze(marg_med(:, iLoc, :)), [], 2);
        %         peak{ifeature}(iLoc, :) = xaxis(maxInd);
        
        % extra lines
        yline(0, 'handlevisibility', 'off', 'linewidth', 2, 'color', [.7, .7, .7]);
        xline(ifeature-1, 'handlevisibility', 'off', 'linewidth', 2, 'color', [.7, .7, .7]);
        
        % plot raw data
        xaxis_ = xaxis + shift(ifeature) * (-1) ^ iLoc;
        if nLoc==2
            errorbar(xaxis_, marg_ave, marg_SEM, 'color', color, 'CapSize',0, 'linestyle', 'none', 'linewidth', 2)
            plot(xaxis_, marg_ave, 'o', 'linestyle', 'none', 'markersize', 1/sqrt(nfilters)*50, 'markerfacecolor', color, 'markeredgecolor','w')
            
            %             if iLoc== 1, color_idvd = [.9, .8, 1]; else,  color_idvd = [.5, .9, 1]; end
            %             for isubj = 1:nsubj
            %             plot(xaxis_, squeeze(margPred_med(isubj, iLoc, :)), 'color', color_idvd, 'linewidth', 1)
            %             end
        else
            plot(xaxis_, marg_ave, 'o-', 'color', color_idvd, 'linewidth', 2)
        end
        
        
        % idvd data
        %         for isubj = 1:nsubj
        %             plot(xaxis, marg_med(isubj, :), 'o', 'color', ones(1,3)*.5)
        %             plot(xaxis, margPred_med(isubj, :), '-', 'color', ones(1,3)*.5)
        %             title(margR2_med(isubj))
        %         end
        
        %         plot fitting lines
        if nLoc == 2
            if nsubj>1
                plot(xaxis_, margPred_ave', '-', 'color', color, 'linewidth', wd_border)
                patch([xaxis_, flip(xaxis_)], [margPred_lb', flip(margPred_ub')], color, 'FaceAlpha', .3, 'linestyle', 'none')
            else
                plot(xaxis_, margPred_ave, '-', 'color', color, 'linewidth', wd_border)
                patch([xaxis_, flip(xaxis_)], [margPred_lb, flip(margPred_ub)], color, 'FaceAlpha', .3, 'linestyle', 'none')
            end
        end
        %  print R2
        if nLoc==2
            text(xaxis(1), ytext_R2_all{ifeature}(iLoc), sprintf('R^2 = %.2f (   %.2f)', margR2_ave, margR2_SEM), ...
                'color', color, 'fontsize', sz_text)
        end
    end % end of iiLoc
    
    % xaxis
    xlim(axisLim{ifeature})
    %     xlabel(xlabels_tuning{ifeature}, 'FontSize', sz_label)
    xticks(axisTicks_tuning{ifeature})
    xticklabels(axisTL_tuning{ifeature})
    
    yticks(linspace(ymin, ymax, 4))
    ylim([ymin, ymax])
    
    
    %% pairwise one-tailed t-tests at each channel
    %     if nLoc == 2 && (nsubj>1)
    %         for ifilter = 1:nfilters
    %             %             xline(xaxis(ifilter), 'k-');
    %             %             [h,p_] = ttest(squeeze(marg_med(:, 1, ifilter)), squeeze(marg_med(:, 2, ifilter)));
    %             [h,p_] = ttest(squeeze(marg_med(:, 1, ifilter)), squeeze(marg_med(:, 2, ifilter)), 'tail', 'right');
    %             p = p_*nfilters;
    %             if p<.05, plot(xaxis(ifilter), ystars(ifeature), 'k*', 'MarkerSize', sz_stars, 'linewidth', wd_border), end
    %         end
    %     end
    
    
    
    %% difference in kernels/tuning fxn between 2 loc
    if nsubj>1
        yyaxis right
        
        [marg_diff_ave, ~, ~, marg_diff_SEM] = getCI(marg_med(:, 1, :)-marg_med(:, 2, :), 2, 1);
        [margPred_diff_ave, margPred_diff_lb, margPred_diff_ub, ~, ~] = getCI(margPred_med(:, 1, :)-margPred_med(:, 2, :), 2, 1);
        
        % extra lines
        yline(0, 'handlevisibility', 'off', 'linewidth', 3, 'color', ones(1,3)*.8);
        xline(ifeature-1, 'handlevisibility', 'off', 'linewidth', 1, 'color', ones(1,3)*.8);
        
        % plot raw data
        errorbar(xaxis, marg_diff_ave, marg_diff_SEM, 'k', 'CapSize',0, 'linestyle', 'none', 'linewidth', 2)
        plot(xaxis, marg_diff_ave, 'o', 'linestyle', 'none', 'markersize', 1/sqrt(nfilters)*30, 'markerfacecolor', 'k', 'markeredgecolor','w')
        
        % plot fitting lines
        plot(xaxis, margPred_diff_ave', '-k', 'linewidth', 2)
        patch([xaxis, flip(xaxis)], [margPred_diff_lb', flip(margPred_diff_ub')], 'k', 'FaceAlpha', .3, 'linestyle', 'none')
        ylim([-0.04, .5])
        yticks(-.04:.04:.04)
        
    end
    
    %%
    ax = gca;
    ax.XAxis.FontSize = sz_ticks;
    ax.YAxis(1).FontSize = sz_ticks;
    ax.YAxis(1).Color= 'k';
    ax.LineWidth = wd_border/1.5;
    if nsubj>1
        ax.YAxis(2).FontSize = sz_ticks;
        ax.YAxis(2).Color= ones(1,3)*.5;
    end
    
    
    %% ANOVA2 (2 loc x 29 channels on kernels)
    r_anova = marg_med; % nsubj x nLoc5 x nORI/SF
    anova_text = print_nANOVA(varNames, r_anova(:),  {ind_Loc(:), ind_channel(:)}, nsubj);
    fprintf('[%s]\n%s', namesType{itype}, anova_text)
    %     title(sprintf('[%s]\n%s', namesType{itype}, anova_text))
    
    %%
    folderName = sprintf('%s/%s/marg_%s/', nameFigFolder, nameEnergySource, namesFeature{ifeature});
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    
    saveas(gcf, sprintf('%sn%d_%s_%s.jpg', folderName, nsubj, title_, namesType{itype}))
end % end of ifeature


