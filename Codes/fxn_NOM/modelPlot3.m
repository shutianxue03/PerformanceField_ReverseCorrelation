
sz_fig = [250 440];

for iNOMp = 1:2%nNoise% induced and constant noise
    
    if flag_plotIDVD
        switch iNOMp            
            case 1, y_ticks = linspace(0, 2, 5);
            case 2, y_ticks = linspace(0, .28, 5);
            case 3, y_ticks = linspace(.7, 2.3, 5);
        end
    else
        switch iNOMp
            case 1, y_ticks = linspace(0, 1.2, 5);
                %                 case 1, y_ticks = linspace(0, .6, 5);
            case 2, y_ticks = linspace(0, .16, 5);
                %             case 2, y_ticks = linspace(0, .1, 5);
            case 3, y_ticks = linspace(.6, 1.8, 5);
        end
    end
    
    y_ticks = round(y_ticks, 3);
    
    %% save
    folderName = sprintf('%s/temp%d/comp_params/%s/', nameFolderFig, templateType, nameFileLoc);
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HM vs. VM
    basicFxn_drawBars([squeeze(NOM_params_med_allSubj(:, 1, iNOMp)),  squeeze(NOM_params_med_allSubj(:, 2, iNOMp))], [], ...
        colors_comb([6,7], :), namesLocComb([6,7]), ...
        y_ticks, y_ticks, flag_plotIDVD, 1, namesParamsNOM{iNOMp}, flag_pairwiseComp, sz_fig);
    saveas(gcf, sprintf('%sn%d_NOMp%d_A%dB%d_L67.jpg', folderName, nsubj, iNOMp, iModelA, iModelB))
    
    % LVM vs. UVM
    basicFxn_drawBars(squeeze(NOM_params_med_allSubj(:, 3:4, iNOMp)), [], ...
        colors_comb([5,3], :), namesLocComb([5,3]), ...
        y_ticks, y_ticks, flag_plotIDVD, 1, namesParamsNOM{iNOMp}, flag_pairwiseComp, sz_fig);
    saveas(gcf, sprintf('%sn%d_NOMp%d_A%dB%d_L53.jpg', folderName, nsubj, iNOMp, iModelA, iModelB))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end % ip
