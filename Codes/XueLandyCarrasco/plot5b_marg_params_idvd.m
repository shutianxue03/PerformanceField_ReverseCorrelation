close all

% sz_stars = 5; sz_label=25; sz_text = 25; wd_border=1; sz_title=40;
sz_font = 10;
if flagPlotAllSubjInOnePlot==0, sz_text = 10;else, sz_text = 10; end

% ylim_diff  = [-.08, -.04];
if ifeature == 1
    ylim_tuning = [-.15, .2]; % left y-axis, limit
else
    ylim_tuning = [-.1, .15];
    %     ylim_diff = [-.05, ylim_tuning(2)*2];
end
ylim_diff = [-.05, ylim_tuning(2)*2]; % right y-axis, limit
ylim_diff_ = [-.05, .05]; % right y-axis, the window of y-ticks
ytext_R2 = ylim_tuning(2) - [.02, .04]; % left y-axis, the height of R2 for two loc
y_right_window = ylim_diff_(2)-ylim_diff_(1);
y_star = [-y_right_window/2.5, y_right_window/2.5]; % right y-axis, starts indicates multiple comparison of each channel

CI_ratio = .68;

if flagPlotAllSubjInOnePlot == 1, figure('Position', [0, 0, 2e3, 1.8e3]), hold on, end

for isubj = 1:nsubj
    subjName = subjList{isubj};
    
    xaxis = axis_tuning{ifeature};
    nfilters = length(xaxis);
    
    namesParams = {'Gain', 'Peak amp.', 'Bandwidth', 'Baseline', 'Peak SF', 'Truncation'};
    if ifeature == 1
        marg_allSubj = margORI_allSubj;
        margPred_allSubj = margPredORI_allSubj;
        margR2_allSubj = margR2ORI_allSubj;
        %         params_allSubj = cat(5, squeeze(estP_ORI_allSubj(:, :, :, :, 1)), tuningC_ORI_allSubj(:, :, :, :, 1:3));
    else
        marg_allSubj = margSF_allSubj;
        margPred_allSubj = margPredSF_allSubj;
        margR2_allSubj = margR2SF_allSubj;
        if ifamily_perF(2) == 12
%             peakSF1_log_med = getCI(log2(margTuningC_SF_allSubj(isubj, :, :, itype, 1)), 1, 2);
%             peakSF2_log_med = getCI(log2(margTuningC_SF_allSubj(isubj, :, :, itype, 4)), 1, 2);
            
            peakSF1_log_med = getCI(log2(margTuningC_SF_allSubj(isubj, :, :, itype, 1)), 1, 2);
            peakSF2_log_med = getCI(log2(margTuningC_SF_allSubj(isubj, :, :, itype, 4)), 1, 2);
            peakAmp1_med = getCI(margTuningC_SF_allSubj(isubj, :, :, itype, 2), 1, 2);
            peakAmp2_med = getCI(margTuningC_SF_allSubj(isubj, :, :, itype, 5), 1, 2);
        end
        
    end
    %     nparams = size(params_allSubj, 5);
    
    if flagPlotAllSubjInOnePlot == 0
        figure('Position', [3e3, 0, 2e3, 2e3]), hold on
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Fig 1: Kernels & Tuning fxns of two loc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flagPlotAllSubjInOnePlot == 1
        subplot(3,5, isubj), hold on
        yyaxis left
    else
        subplot(3,5, [1:3, 6:8]), hold on
    end
    marg_CI_2 = nan(nfilters, 2); % 1=lb of F/HM/LVM, 2=ub of P/VM/UVM
    
    for iiLoc = 1:nLoc
        color = colors_comb(iLocComb_all(iiLoc), :);
        [marg_med, ~, ~, marg_CI_neg, marg_CI_pos] = getCI(marg_allSubj(isubj, :, iiLoc, itype, :), 1, 2, 1,1, CI_ratio);
        [margPred_med, margPred_lb, margPred_ub, ~, ~] = getCI(margPred_allSubj(isubj, :, iiLoc, itype, :), 1, 2, 1,1, CI_ratio);
        [~, margR2_lb, margR2_ub] = getCI(margR2_allSubj(isubj, :, iiLoc, itype), 1, 2, 1, 1, CI_ratio);
        
        if iiLoc == 1, marg_CI_2(:, iiLoc) = marg_CI_neg;
        else, marg_CI_2(:, iiLoc) = marg_CI_pos;
        end
        
        % extra lines
        yline(0, 'handlevisibility', 'off', 'linewidth', 1, 'color', [.7, .7, .7]);
        xline(ifeature-1, 'handlevisibility', 'off', 'linewidth', 1, 'color', [.7, .7, .7]);
        
        if (~flag_mirrorMapping) && (ifeature==1)
            [~, ymax] =max(marg_med);
            xline(xaxis(ymax), '-', 'color', color, 'linewidth', 2);
        end
        % plot raw data
        if nLoc==2
            errorbar(xaxis, marg_med, marg_CI_neg, marg_CI_pos, 'color', color, 'CapSize',0, 'linestyle', 'none', 'linewidth', 2)
            plot(xaxis, marg_med, 'o', 'linestyle', 'none', 'markersize', 1/sqrt(nfilters)*30, 'markerfacecolor', color, 'markeredgecolor','w')
        else
            plot(xaxis, marg_med, 'o-', 'color', color, 'linewidth', 2)
        end
        
        % where the peaks are
        if ifeature == 2 && ifamily_perF(2)==12
            plot([peakSF1_log_med(iiLoc), peakSF1_log_med(iiLoc)], [-.1, peakAmp1_med(iiLoc)], '-', 'color', color, 'linewidth', 2)
            plot([peakSF2_log_med(iiLoc), peakSF2_log_med(iiLoc)], [-.1, peakAmp2_med(iiLoc)], '--', 'color', color, 'linewidth', 2)
        end
        
        % plot fitting lines
        %         if flag_interpolate>1
        %             xaxis_itp =  linspace(axis_tuning{ifeature}(1), axis_tuning{ifeature}(end), flag_interpolate);
        %         else
        xaxis_itp = xaxis;
        %         end
        if nLoc == 2
            plot(xaxis_itp, margPred_med', '-', 'color', color, 'linewidth', 2)
            patch([xaxis_itp, flip(xaxis_itp)], [margPred_lb', flip(margPred_ub')], color, 'FaceAlpha', .3, 'linestyle', 'none')
        end
        
        % R2
        if ifeature==1, xtext_R2 = xaxis(1); else, xtext_R2 = xaxis(20); end
        text(xtext_R2, ytext_R2(iiLoc), sprintf('R^2 = [%d%%, %d%%]', round(margR2_lb*100), round(margR2_ub*100)), ...
            'color', color, 'fontsize', sz_text)
        
    end % end of iiLoc
    
    % xaxis
    xlim(axisLim{ifeature})
    xticks(axisTicks_tuning{ifeature})
    xticklabels(axisTL_tuning{ifeature})
    
    % yaxis
    ylim(ylim_tuning), yticks(round(linspace(ylim_tuning(1), ylim_tuning(2), 5), 2))
    
    % xaxis
    xlim(axisLim{ifeature})
    xticks(axisTicks_tuning{ifeature})
    xticklabels(axisTL_tuning{ifeature})
    
    if isubj>1, xticklabels([]), yticklabels([]), end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Fig 2: difference in kernels/tuning fxn between 2 loc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if flagPlotAllSubjInOnePlot == 1, yyaxis right
    else, subplot(3,5, 11:13), hold on
    end
    
    [marg_diff_med, margPred_lb, margPred_ub, marg_diff_CI_neg, marg_diff_CI_pos] = getCI(marg_allSubj(isubj, :, 1, itype, :)-marg_allSubj(isubj, :, 2, itype, :), 1, 2, 1,1, CI_ratio);
    [margPred_diff_med, margPred_diff_lb, margPred_diff_ub, ~, ~] = getCI(margPred_allSubj(isubj, :, 1, itype, :)-margPred_allSubj(isubj, :, 2, itype, :), 1, 2, 1,1, CI_ratio);
    
    % extra lines
    yline(0, 'handlevisibility', 'off', 'linewidth', 3, 'color', ones(1,3)*.8);
    xline(ifeature-1, 'handlevisibility', 'off', 'linewidth', 1, 'color', ones(1,3)*.8);
    
    % plot raw data
    errorbar(xaxis, marg_diff_med, marg_diff_CI_neg, marg_diff_CI_pos, 'k', 'CapSize',0, 'linestyle', 'none', 'linewidth', 2)
    plot(xaxis, marg_diff_med, 'o', 'linestyle', 'none', 'markersize', 1/sqrt(nfilters)*30, 'markerfacecolor', 'k', 'markeredgecolor','w')
    
    % plot fitting lines
    plot(xaxis_itp, margPred_diff_med', '-k', 'linewidth', 2)
    patch([xaxis_itp, flip(xaxis_itp)], [margPred_diff_lb', flip(margPred_diff_ub')], 'k', 'FaceAlpha', .3, 'linestyle', 'none')
    
    if flagPlotAllSubjInOnePlot==1
        ylim(ylim_diff)
        yticks(linspace(ylim_diff_(1), ylim_diff_(2), 3))
        %         title(subjName)
        if ifeature==1, margTC = margTuningC_ORI_allSubj; else, margTC = margTuningC_SF_allSubj; end
        title(sprintf('S%d\nL%d: %s\nL%d: %s', isubj, ...
            iLocComb_all(1), num2str(round(getCI(margTC(isubj, :, 1, itype, :), 1, 2), 2).'), ...
            iLocComb_all(2), num2str(round(getCI(margTC(isubj, :, 2, itype, :), 1, 2), 2).')  ))
        sz_stars = 2;
    else
        %         ylim(ylim_diff(ifeature, :))
        %         sz_stars = 5;
        %         title(sprintf('%s minus %s', namesLocComb_{1}, namesLocComb_{2}))
    end
    
    % plot stars (on top of the plot: diff>0; at te bottom: diff < 0)
    %     for ifilter = 1:nfilters
    %         if margPred_ub(ifilter) > 0, plot(xaxis(ifilter), y_star(2), '*', 'color', ones(1,3)*.8, 'MarkerSize', sz_stars), end
    %         if margPred_lb(ifilter) < 0, plot(xaxis(ifilter), y_star(1), '*', 'color', ones(1,3)*.8, 'MarkerSize', sz_stars), end
    %     end
    
%     if isubj>1, xticklabels([]), yticklabels([]), end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Fig 3: params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flagPlotAllSubjInOnePlot == 0
        wd_errbar1 = 2;
        sz_marker1 = 15;
        isubplots = [4, 9, 14, 5, 10, 15];
        [params_med, params_lb, params_ub, params_CI_neg, params_CI_pos] = getCI(params_allSubj(isubj, :, :, itype, :), 1, 2, 1, 1, CI_ratio); % CI for errorbar in the scatter plot
        [params_diff_med, params_diff_lb, params_diff_ub, params_diff_CI_neg, params_diff_CI_pos] = getCI(params_allSubj(isubj, :, 1, itype, :) - params_allSubj(isubj, :, 2, itype, :), 1, 2, 1, 1, CI_ratio); % CI for errorbar in the scatter plot
        for iparam = 1:nparams
            subplot(3,5, isubplots(iparam)), hold on, box on
            
            yyaxis left
            for iiLoc = 1:nLoc
                errorbar(iiLoc, params_med(iiLoc, iparam), params_CI_neg(iiLoc, iparam), params_CI_pos(iiLoc, iparam), '.', ...
                    'color', colors_comb_(iiLoc, :), 'CapSize', 0, 'linewidth', wd_errbar1)
                plot(iiLoc, params_med(iiLoc, iparam), 'o', 'MarkerFaceColor' , colors_comb_(iiLoc, :), ...
                    'MarkerEdgeColor', 'w', 'MarkerSize', sz_marker1, 'linewidth', wd_errbar1)
            end % end of iiLoc
            xticks(1:nLoc)
            xlim([.3, nLoc+.7])
            xticklabels(namesLocComb_)
            ylim([min([params_lb(:, iparam); params_lb(:, iparam)]) - .01, max([params_ub(:, iparam); params_ub(:, iparam)]) + .01])
            
            yyaxis right
            errorbar(1.5, params_diff_med(iparam), params_diff_CI_neg(iparam), params_diff_CI_pos(iparam), 'k', 'CapSize', 0)
            plot(1.5, params_diff_med(iparam), 'ok')
            yline(0, 'color', ones(1,3)*.5);
            ylim([params_diff_lb(iparam) - .01, params_diff_ub(iparam) + .01])
            if params_diff_lb(iparam)>0, star = '[diff>0]';
            elseif params_diff_ub(iparam)<0, star = '[diff<0]';
            else , star = '';
            end
            
            title(sprintf('%s %s', namesParams{iparam}, star))
        end % iparam
        
    end
    
    ax = gca; ax.YAxis(1).Color= 'k'; ax.YAxis(2).Color= ones(1,3)*.5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    folderName = sprintf('%s/%s/idvd/tuningFxn/L%d%d/', nameFigFolder, name_numFilters_Fitting, iLocComb_all);
    folderDir = dir(folderName);
    if isempty(folderDir), mkdir(folderName), end
    
    if flagPlotAllSubjInOnePlot == 0
        sgtitle(sprintf('[%s] %s', namesType{itype}, subjName))
        set(findall(gcf, '-property', 'fontsize'), 'fontsize', sz_font)
        saveas(gcf, sprintf('%s%s_%s.jpg', folderName, subjName, namesType{itype}))
    end
end % isubj

if flagPlotAllSubjInOnePlot == 1
    sgtitle(namesType{itype})
    set(findall(gcf, '-property', 'fontsize'), 'fontsize', sz_font)
    saveas(gcf, sprintf('%sn%d_%s_%s.jpg', folderName, nsubj, namesFeature{ifeature}, namesType{itype}))
end

