

%%%%%%%%%%%%%%
% RCplot_setting_corrEEI
%%%%%%%%%%%%%%
sz_label = 55;

%%
x = cs_allSubj_;
y = squeeze(params_allSubj{ifeature}(indSubj, :, :, itype, :));

namesYaxis = namesTunC_unit_perF{ifamily_perF(ifeature), paramMode};
nparams = length(namesYaxis);

% for estP, only plot gain
if paramMode == 1, if ifeature==1, pp_all = 1; else, pp_all = 2; end
else, pp_all = 1:nparams;
end

for ip = pp_all
    fprintf('\n[n=%d L%d%d] %s%d: %s:\n', nsubj, iLocComb_all, namesFeature{ifeature}, ip, namesTunC_unit_perF{ifamily_perF(ifeature), 2}{ip})
    % specifics for SF
    if ifeature == 2
        if ifamily_perF(ifeature) == 2
        % SF peak, plot log scale
        if ip == 1, y = log2(y); end
        % bandwidth, plot log OR linear scale (saved as the last element)
        if (ip == 3) && ~flag_plotOctave, y = params_allSubj{ifeature}(:, :, :, itype, end); end
        end
    end
    
    yticks_ = [];%linspace(yticks_all{ifeature, ip}(1), yticks_all{ifeature, ip}(2), 5);
%     ylim_ = [min(yticks_), max(yticks_)];
    xticks_ = [];%linspace(xticks_all{ifeature, ip}(1), xticks_all{ifeature, ip}(2), 5);
%     ylim_ = [min(yticks_), max(yticks_)];
        
    figure('Position', [0 0 1200 1000]); hold on, box on
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = fxn_drawCorr_EEI(ip, x, y, xticks_, yticks_, namesType, namesLocComb, iLocComb_all, flag_plotCI, markers_allSubj, itype);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ylabel(sprintf('%s of %s %s', nameAsym, namesFeature{ifeature}, namesYaxis{ip}), 'fontsize', sz_label)
    if p<0.05, title_sig='_sig'; else, title_sig=''; end
    %%
    folderName = sprintf('%s/%s/corrEEI/%s/%s/CS/', nameFigFolder, name_numFilters_Fitting, nameFileLoc, namesType{itype});
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    
    saveas(gcf, sprintf('%sn%d_%s%d_%s%s.jpg', folderName, nsubj, namesFeature{ifeature}, ip, namesParamsMode{paramMode}, title_sig))
% switch iLocComb_all(1), case 5, saveas(gcf, sprintf('XueCarrasco/temp/corrVMA_n%d_%s%d.jpg',  nsubj, namesFeature{ifeature}, ip)), case 6, saveas(gcf, sprintf('XueCarrasco/temp/corrHVA_n%d_%s%d.jpg',  nsubj, namesFeature{ifeature}, ip)), end % delete


end % end of iparam
% end % iBehav

