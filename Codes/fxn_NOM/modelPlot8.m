
namesYaxis = {'Induced noise', 'Constant noise', 'Threshold'};
iModelA = 1;
iModelB=1;
x = cs_allSubj_;
y = squeeze(NOM_params_est_allSubj(:, :, iModelA, iModelB, :, :)); % nsubj x nLoc x ni
%%
yticks_ = -1.5:1:1.5;
ylim_ = yticks_ ([1, end]);
xticks_ = 0:.05:.2;
xlim_ = xticks_ ([1, end]);
    
for ip = 1:3 % induced noise, constant noise, threshold
    figure('Position', [0 0 1200 1000]); hold on, box on
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = fxn_drawCorr_EEI(ip, x, y, xticks_, yticks_, namesType, namesLocComb, iLocComb_all, flag_plotCI, markers_allSubj, itype);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if p<.05, title_sig = 'sig'; else, title_sig = ''; end

    
    %% save
%     folderName = sprintf('%s/%s/corrEEI/%s/', nameFigFolder, nameEnergySource, nameFileLoc);
%     folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
%     saveas(gcf, sprintf('%sn%d_A%dB%d_P%d.jpg', folderName, nsubj, iModelA, iModelB, ip))
end