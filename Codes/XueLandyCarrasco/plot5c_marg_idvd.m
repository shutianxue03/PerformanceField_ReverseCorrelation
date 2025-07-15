close all
nLoc3 = 3;
iLocComb_all3 = [6,5,3];
% sz_stars = 5; sz_label=25; sz_text = 25; wd_border=1; sz_title=40;
sz_font = 15;
% ylim_diff  = [-.08, -.04];
if ifeature == 1
    ylim_tuning = [-.05, .2]; % left y-axis, limit
else
    ylim_tuning = [-.05, .15];
    %     ylim_diff = [-.05, ylim_tuning(2)*2];
end
ylim_diff = [-.05, ylim_tuning(2)*2]; % right y-axis, limit
ylim_diff_ = [-.05, .05]; % right y-axis, the window of y-ticks
ytext_R2 = ylim_tuning(2) - [.02, .05]; % left y-axis, the height of R2 for two loc
y_right_window = ylim_diff_(2)-ylim_diff_(1);
y_star = [-y_right_window/2.5, y_right_window/2.5]; % right y-axis, starts indicates multiple comparison of each channel

CI_ratio = .68;

figure('Position', [0, 0, 2e3, 1.8e3]), hold on

for isubj = 1:nsubj
    subjName = subjList{isubj};
    
    xaxis = axis_tuning{ifeature};
    nfilters = length(xaxis);
    
    namesParams = {'Gain', 'Peak amp.', 'Bandwidth', 'Baseline', 'Peak SF', 'Truncation'};
    if ifeature == 1
        marg_3allSubj = margORI_3allSubj;
        margPred_3allSubj = margPredORI_3allSubj;
        margR2_3allSubj = margR2ORI_3allSubj;
    else
        marg_3allSubj = margSF_3allSubj;
        margPred_3allSubj = margPredSF_3allSubj;
        margR2_3allSubj = margR2SF_3allSubj;
    end
       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Fig 1: Kernels & Tuning fxns of two loc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(3,4, isubj), hold on
    
%     marg_CI_2 = nan(nfilters, 2); % 1=lb of F/HM/LVM, 2=ub of P/VM/UVM
    
    for iiLoc = nLoc3:-1:1
        color = colors_comb(iLocComb_all3(iiLoc), :);
        [marg_med, ~, ~, marg_CI_neg, marg_CI_pos] = getCI(marg_3allSubj(isubj, :, iiLoc, itype, :), 1, 2, 1,1, CI_ratio);
        [margPred_med, margPred_lb, margPred_ub, ~, ~] = getCI(margPred_3allSubj(isubj, :, iiLoc, itype, :), 1, 2, 1,1, CI_ratio);
        [~, margR2_lb, margR2_ub] = getCI(margR2_3allSubj(isubj, :, iiLoc, itype), 1, 2, 1, 1, CI_ratio);
        
%         if iiLoc == 1, marg_CI_2(:, iiLoc) = marg_CI_neg;
%         else, marg_CI_2(:, iiLoc) = marg_CI_pos;
%         end
        
        % extra lines
        yline(0, 'handlevisibility', 'off', 'linewidth', 1, 'color', [.7, .7, .7]);
        xline(ifeature-1, 'handlevisibility', 'off', 'linewidth', 1, 'color', [.7, .7, .7]);
        
        % plot raw data
        errorbar(xaxis, marg_med, marg_CI_neg, marg_CI_pos, '.', 'color', color, 'linewidth', 2, 'CapSize', 0)
        plot(xaxis, marg_med, 'o', 'markerfacecolor', color, 'MarkerEdgeColor', 'w', 'linewidth', 1)
        
        % plot fitting lines
        plot(xaxis, margPred_med', '-', 'color', color, 'linewidth', 2)
        %             patch([xaxis_itp, flip(xaxis_itp)], [margPred_lb', flip(margPred_ub')], color, 'FaceAlpha', .3, 'linestyle', 'none')

    end % end of iiLoc
    
    % xaxis
    xlim(axisLim{ifeature})
    xticks(axisTicks_tuning{ifeature})
    xticklabels(axisTL_tuning{ifeature})
    
    % yaxis
    ylim(ylim_tuning), yticks(round(linspace(ylim_tuning(1), ylim_tuning(2), 4), 2))
    
    % xaxis
    xlim(axisLim{ifeature})
    xticks(axisTicks_tuning{ifeature})
    xticklabels(axisTL_tuning{ifeature})
    
    
end % isubj

folderName = sprintf('%s/%s/idvd/tuningFxn/L%d%d%d/', nameFigFolder, name_numFilters_Fitting, iLocComb_all3);
folderDir = dir(folderName);
if isempty(folderDir), mkdir(folderName), end

sgtitle(namesType{itype})
set(findall(gcf, '-property', 'fontsize'), 'fontsize', sz_font)
saveas(gcf, sprintf('%sn%d_%s_%s.jpg', folderName, nsubj, namesFeature{ifeature}, namesType{itype}))


