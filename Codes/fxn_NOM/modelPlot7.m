
xlim_ = xticks_([1, end]);
xlim_inset = [2, 4]; xticks_inset = 2:4;
yticks_all = {-.6:.6:1.2;-.15:.1: .15;-.15:.15: .15};
yticks_inset_all = {0:.5:1; -.1:.2:.3; .2:.5:1.2};

namesYaxis = {'induced noise', 'constant noise', 'threshold'};

ANOVA_indLoc = repmat(1:nLoc, nsubj, 1); ANOVA_indLoc = ANOVA_indLoc(:);

y = squeeze(paramsNOM_med);

for ip = 1:2%:3
    % get y
    y_med = paramsNOM_med(:, :, ip);
    y_neg = paramsNOM_neg(:, :, ip);
    y_pos = paramsNOM_pos(:, :, ip);
    yticks_ = yticks_all{ip};
    yticks_inset = yticks_inset_all{ip};
    %     ylim_ = yticks_([1, end]);% ylim_all(ip, :);
    %     ylim_inset = yticks_inset([1, end]); %ylim_inset_all(ip, :);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = fxn_drawCorr(ip, x_allSubj, y_med, xticks_, yticks_, ANOVA_indLoc, iLocComb_all, flag_plotCI, flag_plotINSET, flag_plotConnect, markers_allSubj, colors_comb);    %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if p<.05, title_sig = 'sig'; else, title_sig = ''; end
    
    %% save
%     folderName = sprintf('%s/%s/corr/%s/%s/CS/NOM/%s', nameFigFolder, nameEnergySource, nameFileLoc, namesType{itype}, nameMetric);
%     
%     folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
%     saveas(gcf, sprintf('%s/n%d_A%dB%d_P%d.jpg', folderName, nsubj, iModelB, ip))
end
