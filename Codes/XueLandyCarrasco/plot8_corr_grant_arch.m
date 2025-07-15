
nLoc_corr = size(params_allSubj{1}, 3);
sz_label = 60; %40

%% for (partial) corr
ANOVA_indLoc = repmat(1:nLoc_corr, nsubj, 1); ANOVA_indLoc = ANOVA_indLoc(:);
% ANOVA_indSubj = repmat((1:nsubj)', 1, nLoc_corr); ANOVA_indSubj = ANOVA_indSubj(:);

namesYaxis = namesTunC_unit_perF{ifamily_perF(ifeature), paramMode};
nparams = length(namesYaxis);

% for estP, only plot gain
if paramMode == 1, if ifeature==1, pp_all = 1; else, pp_all = 2; end
else, pp_all = 1:nparams;
end

%% get med of x
switch imetric
    case 1, x_med = getCI(cs_allSubj_, 1, 2); iNameMetric = 9;
    case 2, x_med = getCI(squeeze(metrics_allSubj(:, :, :, 6)), 1, 2); iNameMetric = 6;
    case 3, x_med = paramsNOM_med(:, :, 1); % nsubj x nLoc3 x nParamNOM [inducedN, constN, thresh]
end
x = x_med;

for ip = pp_all
    
%     pp_allSubj = params_allSubj{ifeature}(indSubj, :, :, itype, ip);
    pp_allSubj = params_allSubj{ifeature}(indSubj, :, :, itype, ip);
    % specifics for SF
    if ifeature == 2
        if ifamily_perF(2) == 2
            % SF peak, plot log scale
            if ip == 1, pp_allSubj = log2(params_allSubj{ifeature}(indSubj, :, :, itype, ip)); end
            % SF bandwidth, plot log OR linear scale (saved as the last element)
            if (ip == 3) && (~flag_plotOctave), pp_allSubj = params_allSubj{ifeature}(indSubj, :, :, itype, end); end
        end
    end
    
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
    fprintf('%s %s:\n', namesFeature{ifeature}, namesYaxis{ip})
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = fxn_drawCorr(ip, x, y, xticks_, yticks_, ANOVA_indLoc, iLocComb_all, flag_plotCI, flag_plotINSET, flag_plotConnect, markers_allSubj, colors_comb);    %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ylabel(sprintf('%s %s', namesFeature{ifeature}, namesYaxis{ip}), 'fontsize', sz_label)
    
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


