
if flag_noOL, indSubj = [1:8, 10:nsubj]; assert(strcmp(subjList{9}, 'FH'))
else, indSubj = 1:nsubj;
end

nameAsymX = 'HVA';
nameAsymY = 'VMA';
folderNameFig = sprintf('%s/%s/corrAsym/HVA_vs_VMA/',nameFigFolder, name_numFilters_Fitting);
folderDirFig = dir(folderNameFig); if isempty(folderDirFig), mkdir(folderNameFig), end

%% determine X/Y
nparams = 1;
indMedian = 2;
switch nameVar
    case 'tunC_ORI', nparams = length(namesTunC_unit_perF{ifamily_perF(1), 2});
    case 'tunC_SF', nparams = length(namesTunC_unit_perF{ifamily_perF(2), 2});
    case 'NOMparams', nparams = 2; indMedian = 3;
end

for indParam = 1:nparams
    switch nameVar
        case 'CS'
            load(sprintf('%s/n%d_B1000_perf_L67.mat', nameFolderCompile_BEHAV, nsubj))
            data67_allSubj = cs_allSubj;
            load(sprintf('%s/n%d_B1000_perf_L53.mat', nameFolderCompile_BEHAV, nsubj))
            data53_allSubj = cs_allSubj;
            
%             x_ticks = linspace(-.05, .15, 5); % HVA
%             y_ticks = linspace(0, .12, 5); % VMA (extent is smaller)
            x_ticks = linspace(-.14, .14, 5); % HVA
            y_ticks = x_ticks; % VMA (extent is smaller)
            
        case 'pA'
            load(sprintf('%s/n%d_B1000_perf_L67.mat', nameFolderCompile_BEHAV, nsubj))
            data67_allSubj = pA_allSubj;
            load(sprintf('%s/n%d_B1000_perf_L53.mat', nameFolderCompile_BEHAV, nsubj))
            data53_allSubj = pA_allSubj;
            x_ticks = linspace(-.1, .1, 5); % HVA
            y_ticks = linspace(-.1, .1, 5); % VMA (extent is smaller)
            
        case 'tunC_ORI'
            load(sprintf('%s/%s/n%d_B1000_L67_%d_%d_m_ORI%dSF%d.mat', nameFolderCompile_RC, fileName, nsubj, nORI, nSF, ifamily_perF), 'margTuningC_ORI_allSubj')
            data67_allSubj = margTuningC_ORI_allSubj(:, :, :, itype, indParam);
            load(sprintf('%s/%s/n%d_B1000_L53_%d_%d_m_ORI%dSF%d.mat', nameFolderCompile_RC, fileName, nsubj, nORI, nSF, ifamily_perF), 'margTuningC_ORI_allSubj')
            data53_allSubj = margTuningC_ORI_allSubj(:, :, :, itype, indParam);
            
            switch ifamily_perF(1)
                case 1, x_ticks_allTunC = [300, 60, 50, 320]/100; %y_ticks_allTunC = [5, .5, .6, 3];
                case 8, x_ticks_allTunC = [120, 50, 60, 320, 50, 320]/100; %y_ticks_allTunC = [5, .5, .6, 3];
            end
            y_ticks_allTunC = x_ticks_allTunC;
            x_ticks = linspace(-x_ticks_allTunC(indParam), x_ticks_allTunC(indParam), 5);
            y_ticks = linspace(-y_ticks_allTunC(indParam), y_ticks_allTunC(indParam), 5);
            
        case 'tunC_SF'
            load(sprintf('%s/%s/n%d_B1000_L67_%d_%d_m_ORI%dSF%d.mat', nameFolderCompile_RC, fileName, nsubj, nORI, nSF, ifamily_perF), 'margTuningC_SF_allSubj')
            data67_allSubj = margTuningC_SF_allSubj(:, :, :, itype, indParam);
            load(sprintf('%s/%s/n%d_B1000_L53_%d_%d_m_ORI%dSF%d.mat', nameFolderCompile_RC, fileName, nsubj, nORI, nSF, ifamily_perF), 'margTuningC_SF_allSubj')
            data53_allSubj = margTuningC_SF_allSubj(:, :, :, itype, indParam);
            
            switch ifamily_perF(2)
                case 2, x_ticks_allTunC = [50, 60, 70, 300]/100; %y_ticks_allTunC = [.3, .52, .6, 1.1];
                case 3, x_ticks_allTunC = [.5, .5, .6, 3.6, 1]; %y_ticks_allTunC = [.3, .52, .6, 1.1, 1.1];
            end
            
            y_ticks_allTunC = x_ticks_allTunC;
            x_ticks = linspace(-x_ticks_allTunC(indParam), x_ticks_allTunC(indParam), 5);
            y_ticks = linspace(-y_ticks_allTunC(indParam), y_ticks_allTunC(indParam), 5);
           
        case 'NOMparams'
            load(sprintf('%s/%s/n%d_n%d_L67_N_temp%d_B1_nA%d_conv1_IV1.mat', ...
                nameFolderCompile_NOM, fileName_NOM, nsubj, ni, templateType, nModelsA))
            data67_allSubj = squeeze(NOM_params_est_allSubj(:, :, iModelA, :, :, indParam)); %nsubj x nLoc x nModels x 1 x nB x nNOMparams
            load(sprintf('%s/%s/n%d_n%d_L53_N_temp%d_B1_nA%d_conv1_IV1.mat', ...
                nameFolderCompile_NOM, fileName_NOM, nsubj, ni, templateType, nModelsA))
            data53_allSubj = squeeze(NOM_params_est_allSubj(:, :, iModelA, :, :, indParam)); %nsubj x nLoc x nModels x 1 x nB x nNOMparams
            x_ticks = linspace(-120, 120, 5)/100;
            y_ticks = x_ticks;
    end
    
    if strcmp(nameVar, 'NOMparams')
        HVA_allSubj = (data67_allSubj(:, 1, :)-data67_allSubj(:, 2, :))./(data67_allSubj(:, 1, :)+data67_allSubj(:, 2, :));
        VMA_allSubj = (data53_allSubj(:, 1, :)-data53_allSubj(:, 2, :))./(data53_allSubj(:, 1, :)+data53_allSubj(:, 2, :));
    else
        HVA_allSubj = (data67_allSubj(:, :, 1)-data67_allSubj(:, :, 2))./(data67_allSubj(:, :, 1)+data67_allSubj(:, :, 2));
        VMA_allSubj = (data53_allSubj(:, :, 1)-data53_allSubj(:, :, 2))./(data53_allSubj(:, :, 1)+data53_allSubj(:, :, 2));
    end
    HVA_med_allSubj = getCI(HVA_allSubj, 1, indMedian);
    VMA_med_allSubj = getCI(VMA_allSubj, 1, indMedian);
    
    if nparams==1, text_title = sprintf('%s: HVA vs. VMA', nameVar); nameVarY_fileTitle = ['_', nameVar];
    else, text_title = sprintf('%s %d: HVA vs. VMA', nameVar, indParam); nameVarY_fileTitle = sprintf('_%s%d', nameVar, indParam);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    text_CI = ''; for iB=1:nB,[r, p] = corr(HVA_allSubj(:, iB), VMA_allSubj(:, iB), 'Type', 'Kendall'); p_allB(iB) = p; r_allB(iB) = r; end, [p_med, p_lb, p_ub] = getCI(p_allB, 1, 2); [r_med, r_lb, r_ub] = getCI(r_allB, 1, 2); text_CI = sprintf('r=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f]\n', r_med, r_lb, r_ub, p_med, p_lb, p_ub);
    text_title = sprintf('%s\n%s', text_title, text_CI);
    flag_sig = basicFxn_drawCorrAsym(HVA_med_allSubj(indSubj), VMA_med_allSubj(indSubj), ...
        x_ticks, y_ticks, round(x_ticks*100), round(y_ticks*100), text_title, markers_allSubj);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag_noOL, saveas(gcf, sprintf('%sn%d%s%s_noFH.jpg', folderNameFig, nsubj, nameVarY_fileTitle, flag_sig))
    else, saveas(gcf, sprintf('%sn%d%s%s.jpg', folderNameFig, nsubj, nameVarY_fileTitle, flag_sig))
    end
end % indparam

fprintf('\n ======= DONE ======= \n\n')
