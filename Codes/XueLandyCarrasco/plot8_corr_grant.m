
%% 
switch imetric
    case 1, x_med = getCI(cs_allSubj_, 1, 2); iNameMetric = 9;
    case 2, x_med = getCI(squeeze(metrics_allSubj(:, :, :, 6)), 1, 2); iNameMetric = 6;
    case 3, x_med = paramsNOM_med(:, :, 1); % nsubj x nLoc3 x nParamNOM [inducedN, constN, thresh]
end
x = x_med;


for imode = 1:nModes
    
%     pp_allSubj = params_allSubj{ifeature}(indSubj, :, :, itype, ip);
    pp_allSubj = params_allSubj{ifeature}(indSubj, :, :, itype, iparam);

    
    %% get y_med
    y = pp_allSubj;
    %     y_med = getCI(pp_allSubj, 1, 2);
    
    %% get y ticks, limits
    yticks_ = [];%linspace(yticks_all{ifeature, ip}(1), yticks_all{ifeature, ip}(2), 5);
    %     yticks_ = yticks_all{ifeature, ip}; assert(~isempty(yticks_))
    ylim_ = [];%yticks_([1, end]); %ylim_all(ip, :);
    %     if flag_plotINSET
    %     yticks_inset = yticks_inset_all{ifeature, ip};
    %     ylim_inset = yticks_inset([1, end]);%ylim_inset_all(ip, :);
    %     end
    fprintf('%s %s:\n', namesFeature{ifeature}, namesYaxis{iparam})
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = fxn_drawCorr(iparam, x, y, xticks_, yticks_, ANOVA_indLoc, iLocComb_all, flag_plotCI, flag_plotINSET, flag_plotConnect, markers_allSubj, colors_comb);    %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ylabel(sprintf('%s %s', namesFeature{ifeature}, namesYaxis{iparam}), 'fontsize', sz_label)
    
    %% save
    if p < 0.05, title_sig='_sig'; else, title_sig=''; end
%     if length(iLocComb_all)==2
%         folderName = sprintf('%s/%s/corr/%s/%s/', ...
%             nameFigFolder, name_numFilters_Fitting, nameFileLoc, namesMetrics_plus2{iNameMetric});
%     else
%         folderName = sprintf('%s/%s/corr/all%d/%s/', ...
%             nameFigFolder, name_numFilters_Fitting, length(iLocComb_all), namesMetrics_plus2{iNameMetric});
%     end
%     folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
%     saveas(gcf, sprintf('%sn%d_%s%d_%s%s.jpg', ...
%         folderName, nsubj, namesFeature{ifeature}, ip, namesParamsMode{paramMode}, title_sig))
    
%     if length(iLocComb_all)==2, saveas(gcf, sprintf('XueCarrasco/temp/corr_L%d%d_n%d_%s%d.jpg', iLocComb_all,  nsubj, namesFeature{ifeature}, ip)) , else, saveas(gcf, sprintf('XueCarrasco/temp/corr_Lall3_n%d_%s%d.jpg',  nsubj, namesFeature{ifeature}, ip)), end % delete
    
    
end %  end of ip
% end % iBehav


