


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% search for 'default' and change to default values!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varNames = {'Loc', 'Channel'};

wd_border = 4; % default 5
sz_ticks = 30;% default 35

shift = [0,0];

for iFeature = 1:nFeatures

    % Define directory to save the figure
    

    xaxis = axis_tuning{iFeature};
    nfilters = length(xaxis);

    % Set parameters for ANOVA
    ind_Loc = repmat(1:nLoc, nsubj, 1, nfilters);
    ind_channel = nan(nsubj, nLoc, nfilters);
    for isubj = 1:nsubj, ind_channel(isubj, :, :) = repmat(1:nfilters, nLoc, 1); end

    if iFeature == 1
        [marg_med, marg_CI_lb, marg_CI_ub, marg_neg, marg_pos] = getCI(margORI_allSubj(:, :, :, iType, :), 1, 2);
        [margPred_med, margPred_CI_lb, margPred_CI_ub] = getCI(margPredORI_allSubj(:, :, :, iType, :), 1, 2);
        [margR2_med, margR2_neg, margR2_pos] = getCI(margR2ORI_allSubj(:, :, :, iType), 1, 2);

    else
        [marg_med, marg_CI_lb, marg_CI_ub, marg_neg, marg_pos] = getCI(margSF_allSubj(:, :, :, iType, :), 1, 2);
        [margPred_med, margPred_CI_lb, margPred_CI_ub] = getCI(margPredSF_allSubj(:, :, :, iType, :), 1, 2);
        [margR2_med, ~, ~, margR2_neg, margR2_pos]  = getCI(margR2SF_allSubj(:, :, :, iType), 1, 2);
        if ifamily_perF(2)==12
            peakSF_L_log_med = getCI(log2(margTuningC_SF_allSubj(:, :, :, iType, 1)), 1, 2);
            peakSF_R_log_med = getCI(log2(margTuningC_SF_allSubj(:, :, :, iType, 4)), 1, 2);
            [peakSF_L_log_ave, peakSF_L_log_lb, peakSF_L_log_ub] = getCI(peakSF_L_log_med, 2, 1);
            [peakSF_R_log_ave, peakSF_R_log_lb, peakSF_R_log_ub] = getCI(peakSF_R_log_med, 2, 1);

            peakAmp_L_med = getCI(margTuningC_SF_allSubj(:, :, :, iType, 2), 1, 2);
            peakAmp_R_med = getCI(margTuningC_SF_allSubj(:, :, :, iType, 5), 1, 2);
            peakAmp_L_ave = getCI(peakAmp_L_med, 2, 1);
            peakAmp_R_ave = getCI(peakAmp_R_med, 2, 1);
        elseif find(ifamily_perF(2)==[2,3])
            peakSF_log_med = getCI(log2(margTuningC_SF_allSubj(:, :, :, iType, 1)), 1, 2);
            [peakSF_log_ave, peakSF_log_lb, peakSF_log_ub] = getCI(peakSF_log_med, 2, 1);
            peakAmp_med = getCI(margTuningC_SF_allSubj(:, :, :, iType, 2), 1, 2);
            peakAmp_ave = getCI(peakAmp_med, 2, 1);

        end
    end

    %% Plot kernels & tuning fxn
    figure('Position', [0 0 1.1e3 8e2]) % default 8e2
    %     figure('Position', [0 0 9e2 8e2]) % default 8e2
    hold on%, box on

    % Draw reference lines
    yline(0, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);
    xline(iFeature-1, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);

    % Loop through each location
    for iiLoc = 1:nLoc
        %         flag_LR = 1;
        color = colors_comb(iLocComb_all(iiLoc), :);

        if nsubj>1
            [marg_ave, ~, ~, marg_SEM_neg, marg_SEM_pos] = getCI(marg_med(:, iiLoc, :), 2, 1);
            [margPred_ave, margPred_lb, margPred_ub, ~, ~] = getCI(margPred_med(:, iiLoc, :), 2, 1);
            [margR2_ave, ~, ~, margR2_SEM, ~] = getCI(margR2_med(:, iiLoc), 2, 1);

            %%% NEW BOOTS %%%
            %             if ifeature==1
            %                 [marg_ave, ~, ~, marg_SEM_neg, marg_SEM_pos] = getCI(getCI(margORI_allSubj(:, :, iiLoc, iType, :), 2, 1), 1, 1); % nsubj x nBoot x nLoc x nTypes x nORI x nSF
            %                 [margPred_ave, margPred_lb, margPred_ub] = getCI(getCI(margPredORI_allSubj(:, :, iiLoc, iType, :), 2, 1), 1, 1);
            %                 [margR2_ave, margR2_lb, margR2_ub] = getCI(getCI(margR2ORI_allSubj(:, :, iiLoc, iType), 2, 1), 1, 2);
            %             else
            %                 [marg_ave, ~, ~, marg_SEM_neg, marg_SEM_pos] = getCI(getCI(margSF_allSubj(:, :, iiLoc, iType, :), 2, 1), 1, 1); % nsubj x nBoot x nLoc x nTypes x nORI x nSF
            %                 [margPred_ave, margPred_lb, margPred_ub] = getCI(getCI(margPredSF_allSubj(:, :, iiLoc, iType, :), 2, 1), 1, 1);
            %                 [margR2_ave, margR2_lb, margR2_ub] = getCI(getCI(margR2SF_allSubj(:, :, iiLoc, iType), 2, 1), 1, 2);
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
        xaxis_ = xaxis + shift(iFeature) * (-1) ^ iiLoc;
        %         if flag_interpolate>1
        %             xaxis_itp =  linspace(axis_tuning{ifeature}(1), axis_tuning{ifeature}(end), flag_interpolate);
        %         else
        xaxis_itp = xaxis_;
        %         end

        if nLoc == 2
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

        % Plot fitting lines
        if nLoc == 2
            if nsubj>1
                plot(xaxis_itp, margPred_ave', '-', 'color', color, 'linewidth', wd_border)
                patch([xaxis_itp, flip(xaxis_itp)], [margPred_lb', flip(margPred_ub')], color, 'FaceAlpha', .3, 'linestyle', 'none')
            else
                plot(xaxis_itp, margPred_ave, '-', 'color', color, 'linewidth', wd_border)
                patch([xaxis_itp, flip(xaxis_itp)], [margPred_lb, flip(margPred_ub)], color, 'FaceAlpha', .3, 'linestyle', 'none')
            end
        end
    end % end of iiLoc

    % xaxis
    xlim(axisLim{iFeature})
    %     xlabel(xlabels_tuning{ifeature}, 'FontSize', sz_label)
    xticks(axisTicks_tuning{iFeature})
    xticklabels(axisTL_tuning{iFeature})
    if iFeature==2, xticklabels(round(axisTL_tuning{iFeature}, 1)), end

    yticks(yticks_)
    yticklabels(round(yticks_, 2))
    ylim(yticks_([1, end]))

    %% pairwise one-tailed t-tests at each channel
    for iFilter = 1:nfilters
        %             xline(xaxis(ifilter), 'k-');
        [h,p_] = ttest(squeeze(marg_med(:, 1, iFilter)), squeeze(marg_med(:, 2, iFilter)));
        %         [h,p_] = ttest(squeeze(marg_med(:, 1, ifilter)), squeeze(marg_med(:, 2, ifilter)), 'tail', 'right');
        p = p_*nfilters;
        if p<.05, iFilter, end%plot(xaxis(ifilter), ystars(ifeature), 'k*', 'MarkerSize', sz_stars, 'linewidth', wd_border), end
    end
    fprintf('\nPairwise comparison printed\n\n')
    
    % manually enter below:
    if iFeature==1, x_band = {}; else, x_band = {}; end
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
    %     fprintf('[%s]\n%s', namesType{iType}, anova_text)
    %     title(sprintf('[%s]\n%s\nR^2: %s: %.2f [%.2f, %.2f] // %s: %.2f [%.2f, %.2f]', namesType{iType}, anova_text, ...
    %         namesLocComb{iLocComb_all(1)}, margR2_ave2(1), margR2_lb2(1), margR2_ub2(1), ...
    %         namesLocComb{iLocComb_all(2)}, margR2_ave2(2), margR2_lb2(2), margR2_ub2(2)), 'fontsize', 20)

    title(sprintf('[%s]\n%s\nR^2: %s: %.2f [%.2f] // %s: %.2f [%.2f]', namesType{iType}, anova_text, ...
        namesLocComb{iLocComb_all(1)}, margR2_ave2(1), margR2_SEM2(1), ...
        namesLocComb{iLocComb_all(2)}, margR2_ave2(2), margR2_SEM2(2)), 'fontsize', 20)

    %% Save figure
    saveas(gcf, sprintf('%sn%d_%s_%s.jpg', nameFolder_Fig_NOM_Tuning, nsubj, title_, namesType{iType}))

end % end of ifeature


