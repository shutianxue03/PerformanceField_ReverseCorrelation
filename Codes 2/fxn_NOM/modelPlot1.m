

%% figure setting
wd_border = 4;
sz_title = 20;
sz_label = 50; %40
sz_ticks = 35; % 35
sz_marker = 30;
nMarkerMax = 11;

% dprime, criterion, pC, pHit, pFA, pA, pA_PRS, pA_ABS
ticks_allM = {linspace(0, 2.4, 5); ... % dprime
    linspace(-1.6, 1.6, 5);  ... % criterion
    round(linspace(55, 80, 5), 1);  ... % pC
    round(linspace(50, 100, 5), 1);  ... % pHit
    round(linspace(0, 100, 5), 1);...  % pFA
    round(linspace(50, 80, 5), 1);  ... % pA
    round(linspace(50, 100, 5), 1);  ...% pA_PRS
    round(linspace(50, 100, 5), 1)} ;% pA_ABS

flag_zeroMean = 0; 

for im = 1:8%1:nmetrics
    
    text_title = sprintf('[A%dB%d] %s', iModelA, iModelB, namesMetrics{im});
   
    if any(im==3:8)
        NOM_data_med_allSubj_ = NOM_data_med_allSubj*100;
        NOM_pred_med_allSubj_ = NOM_pred_med_allSubj*100;
    else
        NOM_data_med_allSubj_ = NOM_data_med_allSubj;
        NOM_pred_med_allSubj_ = NOM_pred_med_allSubj;
    end
    basicFxn_drawCorr(squeeze(NOM_data_med_allSubj_(:, :, im)), squeeze(NOM_pred_med_allSubj_(:, :, im)), ...
        colors, ticks_allM{im}, ticks_allM{im}, [], [], flag_zeroMean, 'pearson', 'both', text_title, markers_allSubj);
    plot([ticks_allM{im}(1), ticks_allM{im}(end)], [ticks_allM{im}(1), ticks_allM{im}(end)], 'k--', 'LineWidth', 1.5)
    % diagonal unit line (dashed)
    
%     
%     figure('Position', [0 200 1e3 1e3]); hold on, box on
%     for iLoc = 1:nLoc
%         if flag_plotCI
%             errorbar(NOM_data_med_allSubj(:, iLoc, im), NOM_pred_med_allSubj(:, iLoc, im), ...
%                 NOM_pred_neg_allSubj(:, iLoc, im), NOM_pred_pos_allSubj(:, iLoc, im), ...
%                 NOM_data_neg_allSubj(:, iLoc, im), NOM_data_pos_allSubj(:, iLoc, im),  '.','color', colors(iLoc, :), 'CapSize', 0)
%         end
%         for isubj = 1:nsubj
% %                         if isubj<=nMarkerMax,
%             facecolor = 'w';
% %                         else, facecolor = colors(iLoc, :); end
%             plot(NOM_data_med_allSubj(isubj, iLoc, im), NOM_pred_med_allSubj(isubj, iLoc, im), markers_allSubj{isubj}, ...
%                 'MarkerFaceColor', facecolor, 'MarkerEdgeColor', colors(iLoc, :), 'MarkerSize', sz_marker, 'linewidth', wd_border)
%         end % isubj
%     end % iLoc
%     
%     % diagonal line
%     plot([ticks_allM{im}(1), ticks_allM{im}(end)], [ticks_allM{im}(1), ticks_allM{im}(end)], 'k-', 'linewidth', wd_border)
%     
%     % anova (if pred vs. data has a main effect)
%     data_anova = [NOM_data_med_allSubj(:, :, im), NOM_pred_med_allSubj(:, :, im)];
%     % loc
%     iv1_anova = [repmat(1:nLoc, nsubj, 1), repmat(1:nLoc, nsubj, 1)];
%     % pred vs. data
%     iv2_anova = [ones(nsubj, nLoc), ones(nsubj, nLoc)*2];
%     text_anova = print_nANOVA({'Loc', 'Pred vs. Data'}, data_anova(:), {iv1_anova(:), iv2_anova(:)}, nsubj, 1);
%     
%     % correlation
%     a = NOM_data_med_allSubj(:, :, im);
%     b = NOM_pred_med_allSubj(:, :, im);
%     ANOVA_indLoc = repmat(1:nLoc, nsubj, 1);
%     [r_partial, p_partial] = partialcorr(a(:), b(:), ANOVA_indLoc(:));
%     text_corr = sprintf('Partial r=%.2f, p=%.3f\n', r_partial, p_partial);
%     
%     if nLoc==2
%         % ttest
%         ttest(NOM_data_med_allSubj(:, :, im), NOM_pred_med_allSubj(:, :, im))
%     end
%     %%
%     ax = gca;
%     ax.XAxis.FontSize = sz_ticks;
%     ax.YAxis.FontSize = sz_ticks;
%     ax.LineWidth = wd_border;
%     
%     %%%%
%     xlabel('Data')
%     ylabel('Prediction')
%     xlim(ticks_allM{im}([1,end]))
%     ylim(ticks_allM{im}([1,end]))
%     xticks(ticks_allM{im})
%     yticks(ticks_allM{im})
%     
%     axis square
%     title(sprintf('[A%dB%d] %s\n%s%s', iModelA, iModelB, namesMetrics{im}, text_anova, text_corr))
    
    %% save
    folderName = sprintf('%s/temp%d/data_vs_pred/%s/%s/', nameFolderFig, templateType, nameFileLoc, namesMetrics{im});
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    
    saveas(gcf, sprintf('%sn%d_A%dB%d.jpg', folderName, nsubj, iModelA, iModelB))
    
end % im

