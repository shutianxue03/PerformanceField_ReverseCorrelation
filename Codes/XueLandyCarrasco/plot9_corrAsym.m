% basicFxn_drawCorrAsym

indSubj = 1:nsubj;
iModelB = 1;

iLocComb_all = ii{1};

switch iLocComb_all(1)
    case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
    case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
end

nsubj, nameFileLoc = sprintf('%d%d', iLocComb_all);

folderNameFig = sprintf('%s/%s/corrAsym/%s%s_vs_%s%s/', nameFigFolder, name_numFilters_Fitting, nameVarX, nameAsymX, nameVarY, nameAsymY);
folderDir = dir(folderNameFig); if isempty(folderDir), mkdir(folderNameFig), end

%-------------------%
% decide asym on X %
%-------------------%
load(sprintf('%s/n%d_B1000_perf_L%s.mat', nameFolderCompile_BEHAV, nsubj, nameFileLoc))

switch nameVarX
    case 'CS'
        asymX_allSubj = (cs_allSubj(:, :, 1)-cs_allSubj(:, :, 2))./(cs_allSubj(:, :, 1)+cs_allSubj(:, :, 2));
        switch iLocComb_all(1)
            case 1, x_ticks = linspace(0, 20, 5); % EE
            case 6, x_ticks = linspace(-4, 12, 5); % HVA
            case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
        end
    case 'pA'
        asymX_allSubj = (pA_allSubj(:, :, 1)-pA_allSubj(:, :, 2))./(pA_allSubj(:, :, 1)+pA_allSubj(:, :, 2));
        switch iLocComb_all(1)
            case 6, x_ticks = linspace(0, 16, 5); % HVA
            case 5, x_ticks = linspace(0, 10, 5); % VMA (extent is smaller)
        end
end
asymX_med_allSubj = getCI(asymX_allSubj, 1, 2);
if iLocComb_all(1)==1, asymX_med_allSubj = mm_med; end
%-------------------%
% decide asym on Y %
%-------------------%
switch nameVarY
    case 'tunC'
        load(sprintf('%s/%s/n%d_B1000_L%s_%d_%d_m_ORI%dSF%d.mat', nameFolderCompile_RC, fileName, nsubj, nameFileLoc, nORI, nSF, ifamily_perF))
        
        for ifeature=1:2
            ifamily = ifamily_perF(ifeature);
            namesTunCs = namesTunC_unit_perF{ifamily, paramMode};
            nTunCs_full = length(namesTunCs);
            for iTunC = 1:nTunCs_full
                switch ifeature
                    case 1, tunC_allSubj = squeeze(margTuningC_ORI_allSubj(:, :, :, itype, iTunC));
                    case 2, tunC_allSubj = squeeze(margTuningC_SF_allSubj(:, :, :, itype, iTunC));  % nsubj x nB x nLoc x ntypes x nTunC
                end
                
                asymY_allSubj = (tunC_allSubj(:, :, 1)-tunC_allSubj(:, :, 2))./(tunC_allSubj(:, :, 1)+tunC_allSubj(:, :, 2));
                asymY_med_allSubj = getCI(asymY_allSubj, 1, 2);
                
                switch iLocComb_all(1)
                    case 1 % EE
                        y_ticks_allTunC_lb = 0;
                        y_ticks_allTunC_ub = 60;
                        
                    case 6 % HVA
                        switch ifamily
                            case 1, y_ticks_allTunC_lb = -[100, 50, 50, 100]; y_ticks_allTunC_ub = [100, 50, 50, 100];
                            case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                            case 2, y_ticks_allTunC_lb = -[50, 60, 50, 145]; y_ticks_allTunC_ub = [50, 60, 70, 265];
                                %                                 case 3, y_ticks_allTunC = [20, 50, 24, 360, 100];
                        end
                    case 5 % VMA
                        switch ifamily
                            case 1, y_ticks_allTunC_lb = -[100, 50, 50, 50]; y_ticks_allTunC_ub = [100, 50, 50, 50];
                            case 8, y_ticks_allTunC_lb = -[100, 40, 300, 320, 50, 320]; y_ticks_allTunC_ub = [100, 44, 300, 320, 50, 320];
                            case 2, y_ticks_allTunC_lb = -[50, 60, 60, 250]; y_ticks_allTunC_ub = [50, 60, 60, 250];
                                %                                 case 3, y_ticks_allTunC = [30, 50, 60, 110, 110];
                        end
                end
                y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);
                
                x_ticklabels = nan;
                y_ticklabels = nan;
                
                nameVarY_figTitle = sprintf('%s %s', namesFeature{ifeature}, namesTunC_unit_perF{ifamily_perF(ifeature), 2}{iTunC});
                nameVarY_fileTitle = sprintf('_%s%d', namesFeature{ifeature}, iTunC);
                
                text_title = sprintf('%s (%s) vs. %s (%s)', nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                text_CI = ''; for iB=1:nB,[r, p] = corr(asymX_allSubj(:, iB), asymY_allSubj(:, iB), 'Type', 'Kendall', 'Tail','right'); p_allB(iB) = p; r_allB(iB) = r; end, [p_med, p_lb, p_ub] = getCI(p_allB, 1, 2); [r_med, r_lb, r_ub] = getCI(r_allB, 1, 2); text_CI = sprintf('r=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f]', r_med, r_lb, r_ub, p_med, p_lb, p_ub);
                text_title = sprintf('%s\n%s', text_title, text_CI);
                flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj(indSubj)*100, asymY_med_allSubj(indSubj)*100, ...
                    x_ticks, y_ticks, x_ticklabels, y_ticklabels, text_title, markers_allSubj(indSubj));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if flag_noOL, saveas(gcf, sprintf('%sn%d%s%s_noOL.jpg', folderNameFig, nsubj, nameVarY_fileTitle, flag_sig))
                else, saveas(gcf, sprintf('%sn%d%s%s.jpg', folderNameFig, nsubj, nameVarY_fileTitle, flag_sig))
                end
            end % iTunC
        end % ifeature
        
    case 'pA'
        asymY_allSubj = (pA_allSubj(:, :, 1)-pA_allSubj(:, :, 2))./(pA_allSubj(:, :, 1)+pA_allSubj(:, :, 2));
        asymY_med_allSubj = getCI(asymY_allSubj, 1, 2);
        text_title = sprintf('%s (%s) vs. %s (%s)', nameVarX, nameAsymX, nameVarY, nameAsymY);
        nameVarY_fileTitle = '';
        y_ticks = linspace(-6, 6, 5);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        text_CI = ''; for iB=1:nB,[r, p] = corr(asymX_allSubj(:, iB), asymY_allSubj(:, iB), 'Type', 'Kendall', 'Tail','right'); p_allB(iB) = p; r_allB(iB) = r; end, [p_med, p_lb, p_ub] = getCI(p_allB, 1, 2); [r_med, r_lb, r_ub] = getCI(r_allB, 1, 2); text_CI = sprintf('r=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f]', r_med, r_lb, r_ub, p_med, p_lb, p_ub);
        text_title = sprintf('%s\n%s', text_title, text_CI);
        flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj(indSubj)*100, asymY_med_allSubj(indSubj)*100, ...
            x_ticks, y_ticks, [], [], text_title, markers_allSubj);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if flag_noOL, saveas(gcf, sprintf('%sn%d%s%s_noOL.jpg', folderNameFig, nsubj, nameVarY_fileTitle, flag_sig))
        else, saveas(gcf, sprintf('%sn%d%s%s.jpg', folderNameFig, nsubj, nameVarY_fileTitle, flag_sig))
        end
        
    case 'NOMparams'
        load(sprintf('%s/%s/n%d_n%d_L%s_N_temp%d_B1_nA%d_conv1_IV1.mat', ...
            nameFolderCompile_NOM, fileName_NOM, nsubj, ni, nameFileLoc, templateType, nModelsA))
        for iNOMp=1:2
            
            NOMparam_allSubj = squeeze(NOM_params_est_allSubj(:, :, iModelA, :, :, iNOMp));
            
            asymY_allSubj = (NOMparam_allSubj(:, 1, :)-NOMparam_allSubj(:, 2, :))./(NOMparam_allSubj(:, 1, :)+NOMparam_allSubj(:, 2, :));
            asymY_med_allSubj = getCI(asymY_allSubj, 1, 3);
            
            y_ticks = linspace(-120, 120, 5);
            x_ticklabels = nan;
            y_ticklabels = nan;
            
            nameVarY_figTitle = sprintf('%s', namesParamsModel_all{iModelB}{iNOMp});
            nameVarY_fileTitle = sprintf('_NOM%d', iNOMp);
            
            text_title = sprintf('%s (%s) vs. %s (%s)', nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            text_CI = ''; for iB=1:nB,[r, p] = corr(asymX_allSubj(:, iB), asymY_allSubj(:, iB), 'Type', 'Kendall', 'Tail','left'); p_allB(iB) = p; r_allB(iB) = r; end, [p_med, p_lb, p_ub] = getCI(p_allB, 1, 2); [r_med, r_lb, r_ub] = getCI(r_allB, 1, 2); text_CI = sprintf('r=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f]', r_med, r_lb, r_ub, p_med, p_lb, p_ub);
            text_title = sprintf('%s\n%s', text_title, text_CI);
            flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj(indSubj)*100, asymY_med_allSubj(indSubj)*100, ...
                x_ticks, y_ticks, x_ticklabels, y_ticklabels, text_title, markers_allSubj);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if flag_noOL, saveas(gcf, sprintf('%sn%d%s%s_noOL.jpg', folderNameFig, nsubj, nameVarY_fileTitle, flag_sig))
            else, saveas(gcf, sprintf('%sn%d%s%s.jpg', folderNameFig, nsubj, nameVarY_fileTitle, flag_sig))
            end
            
        end % for iNOMp
end

fprintf('\nL%s\n', num2str(iLocComb_all))


% fprintf('\n ======= DONE ======= \n\n')

% close all
