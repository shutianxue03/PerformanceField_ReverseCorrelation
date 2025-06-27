%%
addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('SX_toolbox'))

clear all, clc, close all, warning off

%% PMF fitting
%--------------%
fxn_PMF_setting
%--------------%
flag_binData = 1;
flag_filterData = 1;
flag_addData = 1;
nStairs = 2;
nLoc = 5; nLocSingle = nLoc;
nTrialsPerStair = 50;
nameFolder_fig_thresh = sprintf('XueCarrasco_JN/fig/THRESH'); if isempty(dir(nameFolder_fig_thresh)), mkdir(nameFolder_fig_thresh), end

%%
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'DT', 'DU'}; % Initials of 12 subjects
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 100];
subjList = {'SQ'};
nblocks_allSubj = [0];

nsubj = length(subjList);

thresh_PMF_allSubj = nan(nsubj, nLoc, nModels);
thresh_stair_allSubj = nan(nsubj, nLoc);

for isubj = 1:nsubj
    
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    fprintf('\n*** %s (%d/%d) ***', subjName, isubj, nsubj);
    
    %% compile ccc
    addpath(genpath(sprintf('Data/%s', subjName)))
    dirFile = dir(sprintf('Data/%s/%s_stair*', subjName, subjName));
    load(dirFile.name, 'record')
    
    ccc_allTrials = nan(nStairs*nTrialsPerStair, 4);
    itrial = 1;
    for iLoc = 1:nLoc
        itrialPerStair = 1;
        for iStair=1:nStairs
            for itrialPerStair = 1:nTrialsPerStair
                ccc_allTrials(itrial, 1) = iLoc;
                ccc_allTrials(itrial, 2) = iStair;
                ccc_allTrials(itrial, 3) = record.stairs_all{iLoc, iStair}.x(itrialPerStair);
                ccc_allTrials(itrial, 4) = record.stairs_all{iLoc, iStair}.response(itrialPerStair);
                itrial = itrial+1;
            end
        end
    end
    
    %% fit PMF
    nameFile_PMF = sprintf('Data_OOD/%s%d', subjName, nblocks);
    if isempty(dir(nameFile_PMF)), mkdir(nameFile_PMF), end
    nameFile_PMF = sprintf('%s/PMF.mat', nameFile_PMF);
    if isempty(dir(nameFile_PMF))
        
        cst_log_unik_all = cell(nLoc, 1);
        nCorr_all = cst_log_unik_all;
        nData_all = cst_log_unik_all;
        pC_all = cst_log_unik_all;
        yfit_all = cst_log_unik_all;
        thresh_all = cst_log_unik_all;
        LL_all = cst_log_unik_all;
        
        for iLoc = 1:nLoc
            fprintf('\n   Loc#%d/%d: ', iLoc, nLoc)
            ccc_full = ccc_allTrials(ccc_allTrials(:, 1)==iLoc, :);
            
            yfit_allB = nan(fit.nBoot, nModels, length(fit.curveX));
            thresh_allB = nan(fit.nBoot, nModels, nPerf);
            LL_allB = nan(fit.nBoot, nModels);
            
            %--------%
            fxn_fitPMF
            %--------%
            
            % save
            cst_log_unik_all{iLoc} = cst_log_unik;
            nCorr_all{iLoc} = nCorr;
            nData_all{iLoc} = nData;
            pC_all{iLoc} = pC;
            yfit_all{iLoc} = yfit_allB;
            thresh_all{iLoc} = thresh_allB;
            LL_all{iLoc} = thresh_allB;
        end % iLoc
        
        save(nameFile_PMF, '*_all')
        fprintf('\n   PMF SAVED\n')
    else
        load(nameFile_PMF)
        fprintf('\n   PMF LOADED\n')
    end
    
    %% PLOT
    %----------%
    fxn_plotPMF
    %----------%
    
    %----------%
    fxn_plotStairs
    %----------%
    
    thresh_PMF_allSubj(isubj, :, :) = thresh_PMF_allLoc; % nLoc x nModels
    thresh_stair_allSubj(isubj, :) = thresh_stair_single;
    close all
end % isubj

if nsubj>2
    save(sprintf('Data_compile/n%d_thresh', nsubj), 'thresh*')
end

%% get thresh finally in use (mean of the last 10 sessions)
thresh_final_allSubj = nan(nsubj, nLoc);

for isubj= 1:nsubj
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    load(sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjName, nblocks, subjName), 'cst_perSess_perLoc');
    
    thresh_final_allSubj(isubj, :) = mean(log10(cst_perSess_perLoc(end-6:end, :)));
end

%% plot thresh estimated from PMF against staircase
SX_RC1_setting

min_ = -1;
max_ = 0;

title_x_all = {'Stair','Final','Final'};
title_y_all = {'PMF','Stair','PMF'};
nCompModes = length(title_y_all);

iModel=2;

for iCompMode = 1:nCompModes
    figure('Position', [0 0 1e3 2e3])
    for iLoc = 1:nLoc
        switch iCompMode
            case 1 % stair vs. PMF
                data_x = thresh_stair_allSubj(:, iLoc);
                data_y = thresh_PMF_allSubj(:, iLoc, iModel);
            case 2 % final vs. stair
                
                data_x = thresh_final_allSubj(:, iLoc);
                data_y = thresh_stair_allSubj(:, iLoc);
            case 3 % final vs. PMF
                data_x = thresh_final_allSubj(:, iLoc);
                data_y = thresh_PMF_allSubj(:, iLoc, iModel);
        end
        subplot(3, 3, iplots5(iLoc)), hold on, box on, axis square, grid on
        
        for isubj=1:nsubj
            switch iCompMode
                case 1 % stair vs. PMF
                    plot(thresh_stair_allSubj(isubj, iLoc), thresh_PMF_allSubj(isubj, iLoc, iModel), markers_allSubj{isubj}, 'color', colors_allM{iModel}, 'MarkerSize', nblocks_allSubj(isubj)/20)
                case 2 % final vs. stair
                    plot(thresh_final_allSubj(isubj, iLoc), thresh_stair_allSubj(isubj, iLoc), markers_allSubj{isubj}, 'color', colors_allM{iModel}, 'MarkerSize', nblocks_allSubj(isubj)/20)
                case 3 % final vs. PMF
                    plot(thresh_final_allSubj(isubj, iLoc), thresh_PMF_allSubj(isubj, iLoc, iModel), markers_allSubj{isubj}, 'color', colors_allM{iModel}, 'MarkerSize', nblocks_allSubj(isubj)/20)
            end
        end
        
        plot([min_, max_], [min_, max_], 'k-')
        
        xlabel(sprintf('Est. thresh from %s', title_x_all{iCompMode}))
        ylabel(sprintf('Est. thresh from %s', title_y_all{iCompMode}))
        
        % corr
        [r, p] = corr(thresh_stair_allSubj(:, iLoc), thresh_final_allSubj(:, iLoc));
        
        % ANOVA
        data = [data_x; data_y];
        
        text_ANOVA = print_nANOVA({'method'}, data(:), {[ones(nsubj, 1); ones(nsubj, 1)*2]}, nsubj);
        
        
        title(sprintf('%s\nCorr: r=%.2f (p=%.3f)\n%s', namesLocComb{iLoc}, r, p, text_ANOVA))
    end % iLoc
    
    sgtitle(sprintf('[n=%d] %s vs. %s', nsubj, title_x_all{iCompMode}, title_y_all{iCompMode}))
    
    set(findall(gcf, '-property', 'fontsize'), 'fontsize',12)
    set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
end % iCompMode


%% ANOVA: assess the main effect of est. method
clc
for iModel = 1:nModels
    data = cat(3, thresh_stair_allSubj, squeeze(thresh_PMF_allSubj(:, :, iModel)));
    % ANOVA
    indLoc = repmat(1:nLoc, nsubj, 1, 2);
    indMethod = nan(nsubj, nLoc, 2);
    for isubj = 1:nsubj, indMethod(isubj, :, :) = repmat(1:2, nLoc, 1); end
    
    text_ANOVA = print_nANOVA({'Loc', 'method'}, data(:), {indLoc(:), indMethod(:)}, nsubj);
    fprintf('\n*** Model #%d ***\n%s', iModel, text_ANOVA)
end

%%%%%%
data = cat(3, thresh_stair_allSubj, thresh_final_allSubj);
% ANOVA
indLoc = repmat(1:nLoc, nsubj, 1, 2);
indMethod = nan(nsubj, nLoc, 2);
for isubj = 1:nsubj, indMethod(isubj, :, :) = repmat(1:2, nLoc, 1); end

text_ANOVA = print_nANOVA({'Loc', 'method'}, data(:), {indLoc(:), indMethod(:)}, nsubj);
fprintf('\n *** stair vs. final ***\n%s', text_ANOVA)


