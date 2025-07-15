close all, clc

% markers_allSubj = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o'};

if ~exist('indSubj', 'var'), indSubj = 1:nsubj; end
type_corr = 'Pearson'; % 'Pearson', 'Spearman'
type_tail = 'both'; % 'both', 'left', 'right' (corr assumed to be posi.)

iLocComb_all_all = {[6,5,3], [6,7], [5,3]};
iLocComb_all_all = {[6,5,3]};

for ii = iLocComb_all_all
    iLocComb_all = ii{1};
    if length(iLocComb_all)==3, nameFileLoc = 'all3'; else, nameFileLoc = sprintf('%d%d', iLocComb_all); end
    colors = colors_comb(iLocComb_all, :);
    
    folderNameFig = sprintf('%s/%s/corr/%s_vs_%s/L%s/', nameFigFolder, name_numFilters_Fitting, nameVarX, nameVarY, nameFileLoc);
    folderDir = dir(folderNameFig); if isempty(folderDir), mkdir(folderNameFig), end
    
    %% variable plot on X-axis
    switch nameVarX
        case 'CS'
            load(sprintf('%s/n%d_B1000_perf_L%s.mat', nameFolderCompile_BEHAV, nsubj, nameFileLoc))
            x_med_allSubj = getCI(cs_allSubj, 1, 2);
            x_ticks = linspace(-.8, .8, 5);
            if flag_zeroMean==2, x_ticks = linspace(-.4, .4, 5); end
            
        case 'pA'
            load(sprintf('%s/n%d_B1000_perf_L%s.mat', nameFolderCompile_BEHAV, nsubj, nameFileLoc))
            x_med_allSubj = getCI(pA_allSubj, 1, 2);
            x_ticks = linspace(-.1, .1, 5);
            if flag_zeroMean==2, x_ticks = linspace(-.08, .08, 5); end
            
        case 'NOMparam1' % induced noise
            load(sprintf('%s/%s/n%d_n%d_L%s_N_temp%d_B1_nA%d_conv1_IV1.mat', ...
                nameFolderCompile_NOM, fileName_NOM, nsubj, ni, nameFileLoc, templateType, nModelsA))
            x_med_allSubj = getCI(NOM_params_est_allSubj(:, :, iModelA, :, :, 1), 1, 5);
            x_ticks = linspace(-.7, .5, 5);
            if flag_zeroMean==2, x_ticks = linspace(-.7, .5, 5); end
            
        case 'NOMparam2' % constant noise
            load(sprintf('%s/%s/n%d_n%d_L%s_N_temp%d_B1_nA%d_conv1_IV1.mat', ...
                nameFolderCompile_NOM, fileName_NOM, nsubj, ni, nameFileLoc, templateType, nModelsA))
            x_med_allSubj = getCI(NOM_params_est_allSubj(:, :, iModelA, :, :, 2), 1, 5);
            x_ticks = linspace(-.24, .36, 5);
            if flag_zeroMean==2, x_ticks = linspace(-.24, .24, 5); end
            
        case 'NOMnLL' % constant noise
            load(sprintf('%s/%s/n%d_n%d_L%s_N_temp%d_B1_nA%d_conv1_IV1.mat', ...
                nameFolderCompile_NOM, fileName_NOM, nsubj, ni, nameFileLoc, templateType, nModelsA), 'NOM_nLL_allSubj')
            x_med_allSubj = getCI(NOM_nLL_allSubj(:, :, iModelA, :, :), 1, 5)/1e3;
            x_ticks = linspace(-.8, 1.2,5);
            if flag_zeroMean==2, x_ticks = linspace(6, 9, 5); end
            
            case 'NOMnLLdelta' % constant noise
            load('n13_NOM_nLL_L4', 'NOM_GoF_allSubj_allM_med')
            x_med_allSubj = NOM_GoF_allSubj_allM_med(:, [1,3,4], :);
            for iLoc = 1:3
            x_med_allSubj(:, iLoc, :) = x_med_allSubj(:, iLoc, :) - min(x_med_allSubj(:, iLoc, :), [], 3);
            end
            x_med_allSubj = squeeze(x_med_allSubj(:, :, 1)); % delta nLL of core model
            x_ticks = linspace(-.8, 1.2,5);
            if flag_zeroMean==2, x_ticks = linspace(6, 9, 5); end
    end
    
    %% variable plot on Y-axis
    switch nameVarY
        case 'tunC'
            load(sprintf('%s/%s/n%d_B1000_L%s_%d_%d_m_ORI%dSF%d.mat', nameFolderCompile_RC, fileName, nsubj, nameFileLoc, nORI, nSF, ifamily_perF))
            if length(iLocComb_all)==3
                margTuningC_ORI_allSubj = margTuningC_ORI_3allSubj;
                margTuningC_SF_allSubj = margTuningC_SF_3allSubj;
            end
            
            nIndY = length(namesTunC_unit_perF{ifamily_perF(1), paramMode}) + length(namesTunC_unit_perF{ ifamily_perF(2), paramMode});
            
        case 'pA'
            load(sprintf('%s/n%d_B1000_perf_L%s.mat', nameFolderCompile_BEHAV, nsubj, nameFileLoc))
            nIndY = 1;
            y_ticks = linspace(-.13,.11, 5);
            if flag_zeroMean==2, y_ticks = linspace(-.06, .06, 5); end
            
        case 'NOMparams'
            load(sprintf('%s/%s/n%d_n%d_L%s_N_temp%d_B1_nA%d_conv1_IV1.mat', ...
                nameFolderCompile_NOM, fileName_NOM, nsubj, ni, nameFileLoc, templateType, nModelsA))
            nIndY=2; % 1=induced; 2=constant internal noise
            
    end
    
    %%
    for indY = 1:nIndY
        
        switch nameVarY
            case 'tunC'
                if indY <= length(namesTunC_unit_perF{ifamily_perF(1), paramMode}), ifeature=1; iTunC = indY;
                else, ifeature=2; iTunC = indY - length(namesTunC_unit_perF{ifamily_perF(1), paramMode});
                end
                ifamily = ifamily_perF(ifeature);
                namesTunCs = namesTunC_unit_perF{ifamily, paramMode};
                nTunCs_full = length(namesTunCs);
                
                nameVarY_figTitle = sprintf('%s %s', namesFeature{ifeature}, namesTunC_unit_perF{ifamily_perF(ifeature), 2}{iTunC});
                nameVarY_fileTitle = sprintf('%s%d', namesFeature{ifeature}, iTunC);
                
                switch ifeature
                    case 1, y_med_allSubj = getCI(margTuningC_ORI_allSubj(:, :, :, itype, iTunC), 1, 2);  % nsubj x nB x nLoc3 x ntypes x nTunC
                    case 2, y_med_allSubj = getCI(margTuningC_SF_allSubj(:, :, :, itype, iTunC), 1, 2);  % nsubj x nB x nLoc3 x ntypes x nTunC
                end
                
                
                switch ifamily
                    case 1, y_ticks_allTunC_lb = -[10, .08, 25, .04]; y_ticks_allTunC_ub = [10, .08, 25, .04]; % pref ori//peak amp//band//baseline
                    case 8, y_ticks_allTunC_lb = -[10, .08, 44, .04, 34, .04]; y_ticks_allTunC_ub = [36, .08, 44, .04, 34, .04];% pref ori//peak amp//trough ori//tourgh mag.//band//baseline
                    case 2, y_ticks_allTunC_lb = -[1.8, .05, .6, .04]; y_ticks_allTunC_ub = [1.8, .07, .6, .04];
                    case 3, y_ticks_allTunC_lb = -[1.5, .04, .8, .04, .04];
                end
                
                if flag_zeroMean==2
                    switch ifamily
                        %                         case 1, y_ticks_allTunC_lb = -[10, .08, 25, .04]; y_ticks_allTunC_ub = [10, .08, 25, .04]; % pref ori//peak amp//band//baseline
                        case 8, y_ticks_allTunC_lb = -[15, .04, 60, .03, 20, .02]; y_ticks_allTunC_ub = [17, .04, 40, .02, 20, .02];% pref ori//peak amp//trough ori//tourgh mag.//band//baseline
                        case 2, y_ticks_allTunC_lb = -[1, .04, .6, .03]; y_ticks_allTunC_ub = [.6, .04, .6, .03];
                            %                         case 3, y_ticks_allTunC_lb = -[1.5, .04, .8, .04, .04];
                    end
                end
                y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);
                
                
            case 'NOMparams'
                y_med_allSubj = getCI(NOM_params_est_allSubj(:, :, iModelA, :, :, indY), 1, 5);
                nameVarY_figTitle = sprintf('%s', namesParamsModel_all{1}{indY});
                nameVarY_fileTitle = sprintf('NOMp%d', indY);
                if indY == 1, y_ticks = linspace(-.7, .5, 5);
                else, y_ticks = linspace(-.2, .4, 5);
                end
                
                if flag_zeroMean==2
                    if indY == 1, y_ticks = linspace(-.4, .4, 5);
                    else, y_ticks = linspace(-.3, .3, 5);
                    end
                end
                
            case 'pA'
                y_med_allSubj = getCI(pA_allSubj, 1, 2);
                nameVarY_figTitle = 'pA';
                nameVarY_fileTitle = 'pA';
        end
        
        if flag_zeroMean == 0, x_ticks=nan; y_ticks=nan; end
        x_ticklabels = x_ticks;
        y_ticklabels = y_ticks;
        
        text_title = sprintf('[L%s] %s vs. %s [%s (%s tail)] [zero mean = %d]', num2str(iLocComb_all), nameVarX, nameVarY_figTitle, type_corr, type_tail, flag_zeroMean);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        flag_sig = basicFxn_drawCorr(x_med_allSubj(indSubj,:), y_med_allSubj(indSubj, :), colors, ...
            x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, type_corr, type_tail, text_title, markers_allSubj(indSubj));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        saveas(gcf, sprintf('%sn%d_M%d_%s%s.jpg', folderNameFig, nsubj, flag_zeroMean, nameVarY_fileTitle, flag_sig))
    end % iIndY
    
    fprintf('\nL%s\n', num2str(iLocComb_all))
    %     close all
end % ii

fprintf('\n ======= DONE ======= \n\n')

