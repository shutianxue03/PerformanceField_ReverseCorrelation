


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% search for 'default' and change to default values!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varNames = {'Loc', 'Channel'};

wd_border = 4; % default 5
sz_ticks = 30;% default 35

shift = [0,0];

for ifeature = 1:2
    
    
    if ifeature==1,
        if iLocComb_all(1)==1, yticks_ = [-.03, 0, .05, .10, .15]; % fov vs. peri, higher ub
        else, yticks_ = [-.02, linspace(0, .12, 4)];
        end
    else,
        if iLocComb_all(1)==1, yticks_ = [-.01, linspace(0, .08, 4)];
        else, yticks_ = [-.01, linspace(0, .06, 4)]; %[-.02, 0, .02, .04, .06];
        end
    end
    
    ymax = max(yticks_);
    ymin = min(yticks_);
    %     yrange = ymax - ymin;
    %     ytext_R2_all = [ymax-yrange/10.5, ymax-yrange/6.5];
    %     ystars  = ymax-yrange/7.5;
    
    xaxis = axis_tuning{ifeature};
    nfilters = length(xaxis);
    
    % for ANOVA
    ind_Loc = repmat(1:nLoc, nsubj, 1, nfilters);
    ind_channel = nan(nsubj, nLoc, nfilters);
    for isubj = 1:nsubj, ind_channel(isubj, :, :) = repmat(1:nfilters, nLoc, 1); end
    
    if ifeature == 1
        [marg_med, marg_CI_lb, marg_CI_ub, marg_neg, marg_pos] = getCI(margORI_allSubj(:, :, :, itype, :), 1, 2);
        [margPred_med, margPred_CI_lb, margPred_CI_ub] = getCI(margPredORI_allSubj(:, :, :, itype, :), 1, 2);
        [margR2_med, margR2_neg, margR2_pos] = getCI(margR2ORI_allSubj(:, :, :, itype), 1, 2);
        
    else
        [marg_med, marg_CI_lb, marg_CI_ub, marg_neg, marg_pos] = getCI(margSF_allSubj(:, :, :, itype, :), 1, 2);
        [margPred_med, margPred_CI_lb, margPred_CI_ub] = getCI(margPredSF_allSubj(:, :, :, itype, :), 1, 2);
        [margR2_med, ~, ~, margR2_neg, margR2_pos]  = getCI(margR2SF_allSubj(:, :, :, itype), 1, 2);
        if ifamily_perF(2)==12
            peakSF_L_log_med = getCI(log2(margTuningC_SF_allSubj(:, :, :, itype, 1)), 1, 2);
            peakSF_R_log_med = getCI(log2(margTuningC_SF_allSubj(:, :, :, itype, 4)), 1, 2);
            [peakSF_L_log_ave, peakSF_L_log_lb, peakSF_L_log_ub] = getCI(peakSF_L_log_med, 2, 1);
            [peakSF_R_log_ave, peakSF_R_log_lb, peakSF_R_log_ub] = getCI(peakSF_R_log_med, 2, 1);
            
            peakAmp_L_med = getCI(margTuningC_SF_allSubj(:, :, :, itype, 2), 1, 2);
            peakAmp_R_med = getCI(margTuningC_SF_allSubj(:, :, :, itype, 5), 1, 2);
            peakAmp_L_ave = getCI(peakAmp_L_med, 2, 1);
            peakAmp_R_ave = getCI(peakAmp_R_med, 2, 1);
        elseif find(ifamily_perF(2)==[2,3])
            peakSF_log_med = getCI(log2(margTuningC_SF_allSubj(:, :, :, itype, 1)), 1, 2);
            [peakSF_log_ave, peakSF_log_lb, peakSF_log_ub] = getCI(peakSF_log_med, 2, 1);
            peakAmp_med = getCI(margTuningC_SF_allSubj(:, :, :, itype, 2), 1, 2);
            peakAmp_ave = getCI(peakAmp_med, 2, 1);
            
        end
    end
    
    %% plot kernels & tuning fxn
    figure('Position', [0 0 1.1e3 8e2]) % default 8e2
    %     figure('Position', [0 0 9e2 8e2]) % default 8e2
    hold on%, box on
    
    % ref lines
    yline(0, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);
    xline(ifeature-1, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);
    
    for iiLoc = 1:nLoc
        %         flag_LR = 1;
        color = colors_comb(iLocComb_all(iiLoc), :);
        
        if nsubj>1
            [marg_ave, ~, ~, marg_SEM_neg, marg_SEM_pos] = getCI(marg_med(:, iiLoc, :), 2, 1);
            [margPred_ave, margPred_lb, margPred_ub, ~, ~] = getCI(margPred_med(:, iiLoc, :), 2, 1);
            [margR2_ave, ~, ~, margR2_SEM, ~] = getCI(margR2_med(:, iiLoc), 2, 1);
            
            %%% NEW BOOTS %%%
            %             if ifeature==1
            %                 [marg_ave, ~, ~, marg_SEM_neg, marg_SEM_pos] = getCI(getCI(margORI_allSubj(:, :, iiLoc, itype, :), 2, 1), 1, 1); % nsubj x nBoot x nLoc x nTypes x nORI x nSF
            %                 [margPred_ave, margPred_lb, margPred_ub] = getCI(getCI(margPredORI_allSubj(:, :, iiLoc, itype, :), 2, 1), 1, 1);
            %                 [margR2_ave, margR2_lb, margR2_ub] = getCI(getCI(margR2ORI_allSubj(:, :, iiLoc, itype), 2, 1), 1, 2);
            %             else
            %                 [marg_ave, ~, ~, marg_SEM_neg, marg_SEM_pos] = getCI(getCI(margSF_allSubj(:, :, iiLoc, itype, :), 2, 1), 1, 1); % nsubj x nBoot x nLoc x nTypes x nORI x nSF
            %                 [margPred_ave, margPred_lb, margPred_ub] = getCI(getCI(margPredSF_allSubj(:, :, iiLoc, itype, :), 2, 1), 1, 1);
            %                 [margR2_ave, margR2_lb, margR2_ub] = getCI(getCI(margR2SF_allSubj(:, :, iiLoc, itype), 2, 1), 1, 2);
            %             end
            
        else
            marg_ave = marg_med(iiLoc, :);
            marg_SEM = (marg_neg(iiLoc, :) + marg_pos(iiLoc, :))/2;
            margPred_ave = margPred_med(iiLoc, :);
            margPred_lb = margPred_CI_lb(iiLoc, :);
            margPred_ub = margPred_CI_ub(iiLoc, :);
            margR2_ave = margR2_med(iiLoc);
            margR2_SEM = (margR2_neg(iiLoc, :) + margR2_pos(iiLoc, :))/2;
        end
        margR2_ave2(iiLoc) = margR2_ave;
        margR2_SEM2(iiLoc) = margR2_SEM;
        %         margR2_lb2(iiLoc) = margR2_lb;
        %         margR2_ub2(iiLoc) = margR2_ub;
        
        % plot raw data
        xaxis_ = xaxis + shift(ifeature) * (-1) ^ iiLoc;
        %         if flag_interpolate>1
        %             xaxis_itp =  linspace(axis_tuning{ifeature}(1), axis_tuning{ifeature}(end), flag_interpolate);
        %         else
        xaxis_itp = xaxis_;
        %         end
        
        if nLoc==2
            if flag_plotIDVD
                for isubj = 1:nsubj
                    plot(xaxis_, squeeze(margPred_med(isubj, iiLoc, :)), 'color', color, 'linewidth', 1)
                end
            end
            
            plot(xaxis_, marg_ave, 'o', 'linestyle', 'none', 'markersize', 1/sqrt(nfilters)*50, 'markerfacecolor', 'w', 'markeredgecolor',color, 'linewidth', 2)
            %             errorbar(xaxis_, marg_ave, marg_SEM, 'color', color, 'CapSize',0, 'linestyle', 'none', 'linewidth', 2)
            errorbar(xaxis_, marg_ave, marg_SEM_neg, marg_SEM_pos, 'color', color, 'CapSize',0, 'linestyle', 'none', 'linewidth', 2)
            
            
        else
            plot(xaxis_, marg_ave, 'o-', 'color', color_idvd, 'linewidth', 2)
        end
        
        %         plot fitting lines
        if nLoc == 2
            if nsubj>1
                plot(xaxis_itp, margPred_ave', '-', 'color', color, 'linewidth', wd_border)
                patch([xaxis_itp, flip(xaxis_itp)], [margPred_lb', flip(margPred_ub')], color, 'FaceAlpha', .3, 'linestyle', 'none')
            else
                plot(xaxis_itp, margPred_ave, '-', 'color', color, 'linewidth', wd_border)
                patch([xaxis_itp, flip(xaxis_itp)], [margPred_lb, flip(margPred_ub)], color, 'FaceAlpha', .3, 'linestyle', 'none')
            end
        end
        
        %% plot where the peak SF is (ifamily=12)
        if ifeature == 2 && flag_plotSFPeak
            plot([peakSF_log_ave(iiLoc), peakSF_log_ave(iiLoc)], [0, peakAmp_ave(iiLoc)], 'color', color, 'linewidth', 2)
            patch([peakSF_log_lb(iiLoc), peakSF_log_ub(iiLoc), peakSF_log_ub(iiLoc), peakSF_log_lb(iiLoc)], ...
                [0, 0, peakAmp_ave(iiLoc), peakAmp_ave(iiLoc)], color, 'FaceAlpha', .1, 'linestyle', 'none')%
            %                     switch flag_LR
            %                         case 1, % left peak
            %                             plot([peakSF_L_log_ave(iiLoc), peakSF_L_log_ave(iiLoc)], [0, peakAmp_L_ave(iiLoc)], 'color', color, 'linewidth', 2)
            %                             patch([peakSF_L_log_lb(iiLoc), peakSF_L_log_ub(iiLoc), peakSF_L_log_ub(iiLoc), peakSF_L_log_lb(iiLoc)], ...
            %                                 [0, 0, peakAmp_L_ave(iiLoc), peakAmp_L_ave(iiLoc)], color, 'FaceAlpha', .1, 'linestyle', 'none')
            %                         case 2, % right peak
            %                             plot([peakSF_R_log_ave(iiLoc), peakSF_R_log_ave(iiLoc)], [0, peakAmp_R_ave(iiLoc)], 'color', color, 'linewidth', 2)
            %                             patch([peakSF_R_log_lb(iiLoc), peakSF_R_log_ub(iiLoc), peakSF_R_log_ub(iiLoc), peakSF_R_log_lb(iiLoc)], ...
            %                                 [0, 0, peakAmp_R_ave(iiLoc), peakAmp_R_ave(iiLoc)], color, 'FaceAlpha', .1, 'linestyle', 'none')%
            %                     end
            for isubj=1:nsubj
                plot([peakSF_log_med(isubj, iiLoc),peakSF_log_med(isubj, iiLoc)], [-1e-3, peakAmp_med(isubj, iiLoc)], [markers_allSubj{isubj}, '-'], 'color', color, 'linewidth', 1, 'MarkerSize', 12),
                %                         switch flag_LR
                %                             case 1, plot([peakSF_L_log_med(isubj, iiLoc),peakSF_L_log_med(isubj, iiLoc)], [-1e-3, peakAmp_L_med(isubj, iiLoc)], [markers_allSubj{isubj}, '-'], 'color', color, 'linewidth', 1, 'MarkerSize', 12),
                %                             case 2, plot([peakSF_R_log_med(isubj, iiLoc), peakSF_R_log_med(isubj, iiLoc)], [-1e-3, peakAmp_R_med(isubj, iiLoc)], [markers_allSubj{isubj}, '--'], 'color', color, 'linewidth', 1, 'MarkerSize', 12)
                %                         end
            end
        end
        
    end % end of iiLoc
    
    % xaxis
    xlim(axisLim{ifeature})
    %     xlabel(xlabels_tuning{ifeature}, 'FontSize', sz_label)
    xticks(axisTicks_tuning{ifeature})
    xticklabels(axisTL_tuning{ifeature})
    if ifeature==2, xticklabels(round(axisTL_tuning{ifeature}, 1)), end
    
    yticks(yticks_)
    yticklabels(round(yticks_, 2))
    ylim(yticks_([1, end]))
    
    %% pairwise one-tailed t-tests at each channel
    for ifilter = 1:nfilters
        %             xline(xaxis(ifilter), 'k-');
        [h,p_] = ttest(squeeze(marg_med(:, 1, ifilter)), squeeze(marg_med(:, 2, ifilter)));
        %         [h,p_] = ttest(squeeze(marg_med(:, 1, ifilter)), squeeze(marg_med(:, 2, ifilter)), 'tail', 'right');
        p = p_*nfilters;
        if p<.05, ifilter, end%plot(xaxis(ifilter), ystars(ifeature), 'k*', 'MarkerSize', sz_stars, 'linewidth', wd_border), end
    end
    fprintf('\nPairwise comparison printed\n\n')
    %
    % %
    % manually enter below:
    if ifeature==1, x_band = {}; else, x_band = {}; end
    nband = length(x_band);
    for iband=1:nband
        x_band_ = xaxis_(x_band{iband});
        nF = length(x_band_);
        patch([x_band_, flip(x_band_)], [ones(1, nF)*ymin, ones(1, nF)*ymax], 'k', 'FaceAlpha', .1, 'linestyle', 'none')
    end
    
    %%
    ax = gca;
    ax.XAxis.FontSize = sz_ticks;
    ax.YAxis.FontSize = sz_ticks;
    ax.LineWidth = wd_border/1.5;
    
    %% ANOVA2 (2 loc x 29 channels on kernels)
    r_anova = marg_med; % nsubj x nLoc5 x nORI/SF
    anova_text = print_nANOVA(varNames, r_anova(:),  {ind_Loc(:), ind_channel(:)}, nsubj);
    %     fprintf('[%s]\n%s', namesType{itype}, anova_text)
    %     title(sprintf('[%s]\n%s\nR^2: %s: %.2f [%.2f, %.2f] // %s: %.2f [%.2f, %.2f]', namesType{itype}, anova_text, ...
    %         namesLocComb{iLocComb_all(1)}, margR2_ave2(1), margR2_lb2(1), margR2_ub2(1), ...
    %         namesLocComb{iLocComb_all(2)}, margR2_ave2(2), margR2_lb2(2), margR2_ub2(2)), 'fontsize', 20)
    
    title(sprintf('[%s]\n%s\nR^2: %s: %.2f [%.2f] // %s: %.2f [%.2f]', namesType{itype}, anova_text, ...
        namesLocComb{iLocComb_all(1)}, margR2_ave2(1), margR2_SEM2(1), ...
        namesLocComb{iLocComb_all(2)}, margR2_ave2(2), margR2_SEM2(2)), 'fontsize', 20)
    
    %% save
    folderName = sprintf('%s/%s/marg_%s/', nameFigFolder, name_numFilters_Fitting, namesFeature{ifeature});
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    saveas(gcf, sprintf('%sn%d_%s_%s.jpg', folderName, nsubj, title_, namesType{itype}))
    
end % end of ifeature


