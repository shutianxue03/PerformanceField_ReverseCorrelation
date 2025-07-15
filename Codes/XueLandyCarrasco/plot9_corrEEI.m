

%%%%%%%%%%%%%%
RCplot_setting_corrEEI
%%%%%%%%%%%%%%

%%
x = cs_allSubj_;
y = squeeze(params_allSubj{ifeature}(:, :, :, itype, :));
namesYaxis = namesParams_unit_perF{paramMode, ifeature};
nparams = length(namesYaxis);

% for estP, only plot gain
if paramMode == 1, if ifeature==1, pp_all = 1; else, pp_all = 2; end
else, pp_all = 1:nparams;
end

for ip = pp_all
    
    % specifics for SF
    if ifeature == 2
        % SF peak, plot log scale
        if ip == 1, y = log2(y); end
        % bandwidth, plot log OR linear scale (saved as the last element)
        if (ip == 3) && ~flag_plotOctave, y = params_allSubj{ifeature}(:, :, :, itype, end); end
    end
    
    figure('Position', [0 0 1200 1000]); hold on, box on
    
    % get EEI at each iteEEIn of boot
    EEI_x_allSubj = nan(nsubj, nB); EEI_y_allSubj = EEI_x_allSubj;
    yticks_ = yticks_all{ifeature, ip}; ylim_ = yticks_ ([1, end]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = fxn_drawCorr_EEI(ip, x, y, xticks_, yticks_, namesType, namesLocComb, iLocComb_all, flag_plotCI, markers_allSubj, itype);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if p<0.05, title_sig='_sig'; else, title_sig=''; end

    %%
    folderName = sprintf('%s/%s/corrEEI/%s/%s/CS/', nameFigFolder, nameEnergySource, nameFileLoc, namesType{itype});
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    
    saveas(gcf, sprintf('%sn%d_%s%d_%s%s.jpg', ...
        folderName, nsubj, namesFeature{ifeature}, ip, namesParamsMode{paramMode}, title_sig))
    
end % end of iparam
% end % iBehav

