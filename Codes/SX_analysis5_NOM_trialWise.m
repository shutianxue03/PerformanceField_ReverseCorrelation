% SX_analysis5_NOM_Trialwise.m
% Author: Shutian Xue
% Purpose: This script compiles, analyzes and plots boostrapped data for each observer saved Data_NOM_Trialwise

clear, clc, close all

addpath(genpath('PF_RC/Codes/')) % manually add to save time
addpath(genpath('PF_RC/Codes/fxn_analysis_RC_v2')) % manually add to save time

% Setting parameters for the analysis
%----------------
SX_RC1_setting
%----------------

set(0, 'DefaultFigureVisible', 'off') % avoid printing figures on the desktop

iModelA_sim_all = [1:2]; % 1=RC-derived template, 2=use IO template, 2=randomize template,
iModelB_sim_all = [1:4]; % see SX_RC1_setting for namesModelB
nIter = 200; % Number of iterations
nJob = 5;
nIterxJob = nIter*nJob;

iLocComb_all = [6,7, 5,3, 1,8]; % combination of locations
iLocGroups_all = {[6,7], [5,3], [1,8]} ; nGroups = length(iLocGroups_all);
iLocSingle_all = 1:5; nLocSingle = length(iLocSingle_all);
iLocSingle_allSets = {[6,5,3], [2,4,5,3], [1,6,5,3], [1,2,4,5,3]};

iFamily_ORI = 1; %1=scaled Gaussian
iFamily_SF = 2; % 2=log parabola
iFamily_perF = [iFamily_ORI, iFamily_SF];
namesMetrics_prob = {'pYES', 'pA'}; nMetrics_prob = length(namesMetrics_prob); namesMetrics_prob_full = {sprintf('Predicted detection prob.\nMeasured detection rate'), sprintf('Predicted consistency prob.\nMeasured resp. consistency')};
iDataset_plotRC = 1; %1=Template set; 2=Full set

iModelA_plot = 1; % just plot the RC-derived
iModelB_plot = 1; % 1=full model; 2=No shared variability; 3=No multiplicative variability
iModelB_plot_all = [1,3]; % generate NOMp-related plots
nNOMparams = length(namesModelBparams{iModelB_plot});
iDataset_plotNOM = 1; % 1=test set, 2=full; search for metrics_allCond

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NOM_mode = 2; % 1= aggregate model, 2 = trial-wise model
flag_subjIsHuman = 1; % 1=human subject, 0=IO data

iModelA_all = iModelA_sim_all;
iModelB_all = iModelB_sim_all;

for iRun=[1,4]
    % Define subject list
    if flag_subjIsHuman
        switch iRun
            case 1
                % n=15, everyone
                subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC', 'SR'}; nRows_subj = 3; nCols_subj = 5;
                nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195];
            case 2
                % n=14: no AS
                subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC', 'SR'}; nRows_subj = 3; nCols_subj = 5;
                nblocks_allSubj = [200, 240, 210, 240, 220, 220, 220, 205, 210, 205, 205, 205, 205, 195];

            case 3
                % n=13: no AS, CS
                subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'HL', 'FH', 'HA',  'DT', 'DU', 'RC', 'SR'}; nRows_subj = 3; nCols_subj = 5;
                nblocks_allSubj = [200, 240, 210, 240, 220, 220, 220, 205, 210, 205, 205, 205, 195];

            case 4
                % n=12: no AS, CS, RC
                subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'HL', 'FH', 'HA', 'DT', 'DU', 'SR'}; nRows_subj = 3; nCols_subj = 4;
                nblocks_allSubj = [200, 240, 210, 240, 220, 220, 220, 205, 210, 205, 205, 195];
        end
        assert(nRows_subj * nCols_subj >= length(subjList), 'ALERT: Not enough subplots for the number of subjects!');

    else
        % Find folders with names starting with 'IO' in Data/Data_NOM_Trialwise_2929
        % Filter for directories whose names start with 'IO'
        nameDir_IO_all = dir(sprintf('%s/IO*', nameFolder_Data_NOM_Trialwise));
        nIOs = length(nameDir_IO_all);
        subjList = cell(1, nIOs);
        % nRows_subj=4; nCols_subj=7;
        for iIO = 1:nIOs
            subjList{iIO} = nameDir_IO_all(iIO).name; % Store the names of IO folders
        end
        assert(nRows_subj * nCols_subj >= nIOs, 'Not enough subplots for the number of subjects!');
        iLocComb_all = 1; % For testing with IO data only
    end

    nsubj = length(subjList);

    markers_allSubj = markers_allSubj_full(1:nsubj);

    nameFolder_Data_SaveCompile = sprintf('%s/n%d_n%d', nameFolder_Data_NOM_Trialwise, nsubj, nIterxJob);

    % Names of the performance metrics being analyzed
    nModelsA = length(iModelA_all);
    nModelsB = length(iModelB_all);
    nLocComb8 = 8;
    iIC_plot = 3; % % Which information criterion to plot (1 = AIC, 2 = AICc, 3 = BIC)

    % Function handles for calculating information criteria (IC)
    getAIC_nLL = @(nLL, nParams) 2*nParams+2*nLL;
    getAICc_nLL = @(nLL, nParams, nData) getAIC_nLL(nLL, nParams) + (2*nParams*(nParams+1))/(nData-nParams-1);
    getBIC_nLL = @(nLL, nParams, nData) nParams*log(nData)+2*nLL;

    % Define folder to save figures (can't put in SX_RC1 because nsubj is not defined yet!)
    nameFolder_Fig_NOM_Trialwise = sprintf('%s/NOM_Trialwise_n%d', nameFolder_Figures, nsubj);
    if isempty(dir(nameFolder_Fig_NOM_Trialwise)), mkdir(nameFolder_Fig_NOM_Trialwise), end

    % Print a header to summarize the setting
    fprintf('NOM trial-wise analysis settings:\n')
    fprintf(' - nsubj = %d\n', nsubj)
    fprintf(' - nIter = %d x %d\n', nIter, nJob)
    fprintf(' - Models A: %s\n', strjoin(string(iModelA_all), ', '));
    fprintf(' - Models B: %s\n', strjoin(string(iModelB_all), ', '));
    fprintf(' - Combined locations:  %s\n', strjoin(string(iLocComb_all), ', '));
    fprintf(' - ORI tuning function:  %s\n', namesFamily_all{iFamily_ORI});
    fprintf(' - SF tuning function:  %s\n', namesFamily_all{iFamily_SF});
    fprintf(' - Number of Bins: %d\n', nBins);
    fprintf(' - Information Criterion to plot: %s\n\n', namesIC{iIC_plot});

    %%%%%% Check missing files %%%%%%
    flag_step_all = [1,2];
    %----------------%
    % checkMissingFiles
    %----------------%

    %% (Time-consuming!) Compile data from boostrapped data

    % % Preallocate arrays for storing results across all conditions
    % metrics_allCond = nan(nModelsA, nLocComb8, nsubj, nIter * nJob, nDatasets, nMetrics); % nDatasets: 1=Full; 2=Test; 11=number of metrics saved in OOD_xx_compIV.m; see "fxn_getMetrics"
    % template_tmpl_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nORI, nSF);
    % template_full_allCond = template_tmpl_allCond;
    % IV_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, nBins);
    % nTrials_allCond = IV_allCond;
    % metric_data_allCond = nan(nModelsA, nModelsB, nLocComb8, nMetrics_prob, nsubj, nIter * nJob, nBins); % "nMetrics" is predefined in this script
    % metric_pred_allCond = metric_data_allCond;
    % nLL_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob);
    % params_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, 3); % 3=Pre-allocate for max params
    % R2_NOM_allCond = nan(nModelsA, nModelsB, nLocComb8, nMetrics_prob, nsubj, nIter * nJob);
    % R2_w_NOM_allCond = R2_NOM_allCond;
    % 
    % % Fitting tuning curves
    % % nDatasets: 1=Full; 2=Template or test
    % sep_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, nDatasets);
    % margORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, nDatasets, nORI);
    % margPred_ORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, nDatasets, nORI);
    % margParams_ORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, nDatasets, length(namesParams_all{iFamily_ORI}));
    % margR2_ORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, nDatasets);
    % margTunC_ORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, nDatasets, length(namesTunC_unit_perF{iFamily_ORI, 2}));
    % 
    % margSF_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, nDatasets, nSF);
    % margPred_SF_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, nDatasets, nSF);
    % margParams_SF_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, nDatasets, length(namesParams_all{iFamily_SF}));
    % margR2_SF_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, nDatasets);
    % margTunC_SF_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIter * nJob, nDatasets, length(namesTunC_unit_perF{iFamily_SF, 2}));
    % 
    % % The main analysis loop for all conditions for each subject (Models x Loc x Metrics)
    % fprintf('\n ===== START COMPILING ===== \n\n')
    % for iJob = 1:nJob
    %     idxIter = (iJob-1)*nIter+1: iJob*nIter;
    %     fprintf('\n\n === Job %d/%d === ', iJob, nJob)
    % 
    %     % Loop through each model (1=RC-derived template; 2=ideal template)
    %     for iModelA = iModelA_all %1:nModelsA
    %         fprintf('\n\n   == Model A%d == ', iModelA)
    % 
    %         % Loop through each model variant ()
    %         for iModelB = iModelB_all %1:nModelsB
    %             fprintf('\n       Model B%d: ', iModelB)
    %             nParamsB = length(namesModelBparams{iModelB});
    % 
    %             % Loop through each SINGLE location
    %             for iLocSingle = iLocSingle_all
    %                 fprintf(' L%d ', iLocSingle)
    % 
    %                 % Loop through each subject
    %                 for isubj = 1:nsubj
    % 
    %                     subjName = subjList{isubj};
    % 
    %                     if flag_subjIsHuman
    %                         % nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, nameIO);
    %                         nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocSingle); % MUST be the same as OOD_NOM_Trialwise_compIV.m, Line 61
    %                     else
    %                         % nameFolder_NOM_save = sprintf('%s/ORI%dSF%d/%s/L%d', nameFolder_NOM0, nORI, nSF, subjName, iLocComb);
    %                         nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, subjName); % MUST be the same as OOD_NOM_Trialwise_compIV.m, Line 68
    %                     end
    % 
    %                     % Load IVs and derived templates
    %                     nameFile_compIV = sprintf('%s/n%d_J%d_A%d_compIV.mat', nameFolder_NOM_save, nIter, iJob, iModelA);
    %                     % if isempty(dir(nameFile_compIV)), error('ALERT: File %s does not exist!', nameFile_compIV); end
    %                     if isempty(dir(nameFile_compIV))
    %                         fprintf(' "%sL%dA%d" ', subjName, iLocSingle, iModelA)
    %                     else
    %                         load(nameFile_compIV, '*_allIter')
    %                     end
    % 
    %                     % Load predictions
    %                     nameFile_fitNOM = sprintf('%s/n%d_J%d_A%dB%d.mat', nameFolder_NOM_save, nIter, iJob, iModelA, iModelB);
    %                     % if isempty(dir(nameFile_fitNOM)), error('ALERT: File %s does not exist!', nameFile_fitNOM); end
    %                     if isempty(dir(nameFile_fitNOM))
    %                         fprintf(' "%sL%dA%dB%d" ', subjName, iLocSingle, iModelA, iModelB)
    %                     else
    %                         load(nameFile_fitNOM, '*_allIter')
    %                     end
    % 
    %                     % Loaded *measured* metrics of the full set and test set
    %                     % "metrics" are defined in OOD_xx_compIV: metrics = [metrics_full; metrics_tmpl; metrics_tmpl; metrics_test];
    %                     metrics_allCond(iModelA, iLocSingle, isubj, idxIter, :, :) = data_metrics_allIter(:, [4,1], :); % see OOD_xx_compIV.m search for "data_metrics_allIter"
    % 
    %                     % Pre-allocate temporary arrays for this subject and model
    %                     IV_allIter = nan(nIter, nBins); % IV for each bin
    %                     nTrials_allIter = IV_allIter; % Trial count for each bin
    %                     metric_data_allIter = nan(nMetrics_prob, nIter, nBins); % Measured data for each bin
    %                     metric_pred_allIter = metric_data_allIter; % Predicted data for each bin
    %                     margTuningC_ORI_allIter = nan(nIter, nDatasets, length(namesTunC_unit_perF{iFamily_ORI, 2}));
    %                     margTuningC_SF_allIter = nan(nIter, nDatasets, length(namesTunC_unit_perF{iFamily_SF, 2}));
    % 
    %                     % Loop through each iteration (within the job)
    %                     for iIter = 1:nIter
    % 
    %                         % Extract IV per bin
    %                         IV_allIter(iIter, :) = pred_metrics_allIter{iIter}.metrics.IV_allBins;
    % 
    %                         % Extract the number of trials per bin
    %                         nTrials_allIter(iIter, :) = pred_metrics_allIter{iIter}.metrics.nTrials_allBins;
    % 
    %                         % Extract metrics (data and predictions) per bin
    %                         pred_metrics = pred_metrics_allIter{iIter}.metrics; % Extract once for faster access
    %                         for iMetric_prob = 1:nMetrics_prob
    %                             switch namesMetrics_prob{iMetric_prob}
    %                                 case 'pYES'
    %                                     metric_data_allIter(iMetric_prob, iIter, :) = pred_metrics.pYES_data_allBins;
    %                                     metric_pred_allIter(iMetric_prob, iIter, :) = pred_metrics.pYES_pred_allBins;
    %                                 case 'pA'
    %                                     metric_data_allIter(iMetric_prob, iIter, :) = pred_metrics.pA_data_allBins;
    %                                     metric_pred_allIter(iMetric_prob, iIter, :) = pred_metrics.pA_pred_allBins;
    %                                 case 'pC'
    %                                     metric_data_allIter(iMetric_prob, iIter, :) = pred_metrics.pC_data_allBins;
    %                                     metric_pred_allIter(iMetric_prob, iIter, :) = pred_metrics.pC_pred_allBins;
    %                             end
    %                         end % iMetric
    %                         % Calculate information critertion based on nLL
    %                         nData = sum(nTrials_allIter(iIter, :));
    % 
    %                         % Compute tuning characteristics based on the fitted parameters
    %                         for iDataset=1:2 % needs to matchOOD_xx_compIV ("for iDataset = 1:2")
    %                             % ORI
    %                             iFeature= 1;
    %                             %======================%
    %                             margTuningC_ORI_allIter(iIter, iDataset, :) = fxn_getTuningC(axis_tuning{iFeature}, iFeature, iFamily_ORI, squeeze(margPred_ORI_allIter(iIter, iDataset, :)), squeeze(margParams_ORI_allIter(iIter, iDataset, :)));
    %                             %======================%
    % 
    %                             % SF
    %                             iFeature=2;
    %                             %======================%
    %                             margTuningC_SF_allIter(iIter, iDataset, :) = fxn_getTuningC(axis_tuning{iFeature}, iFeature, iFamily_SF, squeeze(margPred_SF_allIter(iIter, iDataset, :)), squeeze(margParams_SF_allIter(iIter, iDataset, :)));
    %                             %======================%
    %                         end
    %                     end % end of iIter
    % 
    %                     % Calculate R2 for NOM prediction (for pYES and pA)
    %                     R2_NOM_allIter = nan(nMetrics_prob, nIter);
    %                     R2_w_NOM_allIter = R2_NOM_allIter;
    %                     for iMetric_prob = 1:nMetrics_prob
    %                         for iIter = 1:nIter
    %                             nData = squeeze(nTrials_allIter(iIter, :));
    %                             metric_data  = squeeze(metric_data_allIter(iMetric_prob, iIter, :));  % ground truth
    %                             metric_pred = squeeze(metric_pred_allIter(iMetric_prob, iIter, :)); % prediction
    %                             % Remove NaNs if any
    %                             valid = ~(isnan(metric_data) | isnan(metric_pred));
    %                             metric_data = metric_data(valid);
    %                             metric_pred = metric_pred(valid);
    % 
    %                             % With weighting
    %                             metric_ave = sum(nData .* metric_data) / sum(nData);
    %                             SSres_w = sum(nData .* (metric_data - metric_pred).^2);
    %                             SStot_w = sum(nData .* (metric_data - metric_ave).^2);
    %                             if SStot_w == 0, R2_w = NaN; else, R2_w = 1 - SSres_w/SStot_w; end
    % 
    %                             % No weighting
    %                             SSres = sum((metric_data - metric_pred).^2);
    %                             SStot = sum((metric_data - mean(metric_data)).^2);
    %                             if SStot == 0, R2 = NaN; else, R2 = 1 - SSres/SStot; end
    % 
    %                             R2_NOM_allIter(iMetric_prob, iIter) = R2;
    %                             R2_w_NOM_allIter(iMetric_prob, iIter) = R2_w;
    %                         end
    %                     end
    % 
    %                     % Store results
    %                     % Templates
    %                     [template_tmpl_med, ~, ~, template_tmpl_sem] = getCI(template_tmpl_allIter, 1, 1);
    %                     [template_full_med, ~, ~, template_full_sem] = getCI(template_full_allIter, 1, 1);
    %                     template_tmpl_allCond(iModelA, iModelB, iLocSingle, isubj, :, :) = template_tmpl_med;
    %                     template_full_allCond(iModelA, iModelB, iLocSingle, isubj, :, :) = template_full_med;
    % 
    %                     % Tuning functions: marg, predictions, estimated parameters and separability
    %                     sep_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :) = sep_allIter;
    % 
    %                     margORI_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margORI_allIter; % directly loaded
    %                     margR2_ORI_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :) = margR2_ORI_allIter; % directly loaded
    %                     margPred_ORI_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margPred_ORI_allIter; % directly loaded
    %                     margParams_ORI_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margParams_ORI_allIter; % directly loaded
    %                     margTunC_ORI_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margTuningC_ORI_allIter; % derived above
    % 
    %                     margSF_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margSF_allIter; % directly loaded
    %                     margR2_SF_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :) = margR2_SF_allIter; % directly loaded
    %                     margPred_SF_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margPred_SF_allIter; % directly loaded
    %                     margParams_SF_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margParams_SF_allIter; % directly loaded
    %                     margTunC_SF_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margTuningC_SF_allIter; % derived above
    % 
    %                     % NOM: IVs, predictions, nLL and parameters
    %                     IV_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :) = IV_allIter;
    %                     nTrials_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :) = nTrials_allIter;
    %                     metric_data_allCond(iModelA, iModelB, iLocSingle, :, isubj, idxIter, :) = metric_data_allIter;
    %                     metric_pred_allCond(iModelA, iModelB, iLocSingle, :, isubj, idxIter, :) = metric_pred_allIter;
    %                     nLL_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter) = nLL_test_allIter;
    %                     % IC_nLL_allCond is compiled above
    %                     params_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, 1:nParamsB) = params_est_allIter;
    %                     R2_NOM_allCond(iModelA, iModelB, iLocSingle, :, isubj, idxIter) = R2_NOM_allIter;
    %                     R2_w_NOM_allCond(iModelA, iModelB, iLocSingle, :, isubj, idxIter) = R2_w_NOM_allIter;
    % 
    %                     clear *_allIter
    % 
    %                 end % isubj
    %             end % iLocComb
    %         end % iModelB
    %     end % iModelA
    % end % iJob
    % fprintf('\n\n ==== NOM outputs compiled ====\n\n')
    % 
    % %% Build combined locations by averaging single-location outputs
    % meanOverLoc = @(X, locDim, locOld) mean(X, locDim, 'omitnan');
    % 
    % for iMap = 6:8 % 6=HM, 7=VM, 8=Peri
    %     switch iMap
    %         case 6, locOld = [2 4];      % 6 (HM) = 2 (left) and 4 (right)
    %         case 7, locOld = [3 5];      % 7 (VM) = 3 (upper) and 5 (lower)
    %         case 8, locOld = 2:5;        % 8 (Perifovea) = 2 to 5
    %     end
    % 
    %     % locDim = 2
    %     metrics_allCond(:,iMap,:,:,:,:) = meanOverLoc(metrics_allCond(:,locOld,:,:,:,:), 2, locOld);
    % 
    %     % locDim = 3
    %     template_tmpl_allCond(:,:,iMap,:,:,:)    = meanOverLoc(template_tmpl_allCond(:,:,locOld,:,:,:), 3, locOld);
    %     template_full_allCond(:,:,iMap,:,:,:)    = meanOverLoc(template_full_allCond(:,:,locOld,:,:,:), 3, locOld);
    % 
    %     IV_allCond(:,:,iMap,:,:,:)               = meanOverLoc(IV_allCond(:,:,locOld,:,:,:), 3, locOld);
    %     nTrials_allCond(:,:,iMap,:,:,:)          = meanOverLoc(nTrials_allCond(:,:,locOld,:,:,:), 3, locOld);
    % 
    %     metric_data_allCond(:,:,iMap,:,:,:,:)    = meanOverLoc(metric_data_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    %     metric_pred_allCond(:,:,iMap,:,:,:,:)    = meanOverLoc(metric_pred_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    % 
    %     nLL_allCond(:,:,iMap,:,:)                = meanOverLoc(nLL_allCond(:,:,locOld,:,:), 3, locOld);
    %     params_allCond(:,:,iMap,:,:,:)           = meanOverLoc(params_allCond(:,:,locOld,:,:,:), 3, locOld);
    %     % IC_nLL_allCond(:,:,iMap,:,:,:)           = meanOverLoc(IC_nLL_allCond(:,:,locOld,:,:,:), 3, locOld);
    %     R2_NOM_allCond(:,:,iMap,:,:,:)           = meanOverLoc(R2_NOM_allCond(:,:,locOld,:,:,:), 3, locOld);
    %     R2_w_NOM_allCond(:,:,iMap,:,:,:)           = meanOverLoc(R2_w_NOM_allCond(:,:,locOld,:,:,:), 3, locOld);
    % 
    %     sep_allCond(:,:,iMap,:,:,:)              = meanOverLoc(sep_allCond(:,:,locOld,:,:,:), 3, locOld);
    % 
    %     margORI_allCond(:,:,iMap,:,:,:,:)        = meanOverLoc(margORI_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    %     margPred_ORI_allCond(:,:,iMap,:,:,:,:)   = meanOverLoc(margPred_ORI_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    %     margParams_ORI_allCond(:,:,iMap,:,:,:,:) = meanOverLoc(margParams_ORI_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    %     margR2_ORI_allCond(:,:,iMap,:,:,:)       = meanOverLoc(margR2_ORI_allCond(:,:,locOld,:,:,:), 3, locOld);
    %     margTunC_ORI_allCond(:,:,iMap,:,:,:,:)   = meanOverLoc(margTunC_ORI_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    % 
    %     margSF_allCond(:,:,iMap,:,:,:,:)         = meanOverLoc(margSF_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    %     margPred_SF_allCond(:,:,iMap,:,:,:,:)    = meanOverLoc(margPred_SF_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    %     margParams_SF_allCond(:,:,iMap,:,:,:,:)  = meanOverLoc(margParams_SF_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    %     margR2_SF_allCond(:,:,iMap,:,:,:)        = meanOverLoc(margR2_SF_allCond(:,:,locOld,:,:,:), 3, locOld);
    %     margTunC_SF_allCond(:,:,iMap,:,:,:,:)    = meanOverLoc(margTunC_SF_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    % end
    % 
    % fprintf('\n\n ==== Combined locations (6,7,8) created by averaging single locations ==== \n\n');
    % 
    % %%% Save the organized data for all subjects
    % save(nameFolder_Data_SaveCompile, '*_allCond')
    % clear *_allCond
    % fprintf('\n\n ==== *_allCond saved and cleared ==== \n\n');
    % 
    % %% Compile behavioral data (revise this part later, as the behav should also come from resampled data (metric_data_allCond), not from the raw data)
    % CS_allSubj = nan(nsubj, nLocComb8); % 8 is max number of combined locs,see namesLocComb
    % dprime_allSubj = CS_allSubj;
    % criterion_allSubj = CS_allSubj;
    % RT_allSubj = CS_allSubj;
    % pC_allSubj = CS_allSubj;
    % pA_allSubj = CS_allSubj;
    % 
    % for isubj = 1:nsubj
    %     subjName = subjList{isubj};
    %     nblocks = nblocks_allSubj(isubj);
    %     nameFile_behavMeas = sprintf('%s/%s%d/%s_behavMeas.mat', nameFolder_Data_OOD, subjName, nblocks, subjName);
    % 
    %     load(nameFile_behavMeas, '*_perSess_perLoc')
    % 
    %     for iLoc = 1:nLocComb8
    %         switch iLoc
    %             case 6, iLoc_allLoc = [2,4];
    %             case 7, iLoc_allLoc = [5,3];
    %             case 8, iLoc_allLoc = 2:5;
    %             otherwise
    %                 iLoc_allLoc = iLoc;
    %         end
    %         CS_allSubj(isubj, iLoc) = mean(1./cst_perSess_perLoc(:, iLoc_allLoc), 'all');
    %         dprime_allSubj(isubj, iLoc) = mean(dprime_perSess_perLoc(:, iLoc_allLoc), 'all');
    %         criterion_allSubj(isubj, iLoc) = mean(criterion_perSess_perLoc(:, iLoc_allLoc), 'all');
    %         RT_allSubj(isubj, iLoc) = median(RT_perSess_perLoc(:, iLoc_allLoc), 'all');
    %         pC_allSubj(isubj, iLoc) = mean(pC3_perSess_perLoc(:, iLoc_allLoc, 1), 'all');
    %         % pYES_allSubj(isubj, iLocComb) = mean(pYES_perSess_perLoc(:, iLoc_allLoc));
    %         pA_allSubj(isubj, iLoc) = mean(pA3_perSess_perLoc(:, iLoc_allLoc, 1), 'all');
    %     end % iLoc
    % 
    %     clear *_perSess_perLoc
    % 
    % end % isubj
    % 
    % % Save the organized data for all subjects
    % save(nameFolder_Data_SaveCompile, '*_allSubj', '-append')
    % 
    % fprintf('\n ==== Behav data compiled ==== \n\n')

    %% Behavioral metrics: 4 single locations
    clc, fprintf('\n\n 1/24 Plotting STARTED......\n\n')
    % Load data
    load(nameFolder_Data_SaveCompile, 'metrics_allCond')

    % Define folder for saving figures
    nameFolder_Fig_behav = sprintf('%s/Behav', nameFolder_Fig_NOM_Trialwise);
    if isempty(dir(nameFolder_Fig_behav)), mkdir(nameFolder_Fig_behav), end

    flag_plotIDVD = 1;
    flag_plotDiff = 1;
    namesMetrics_behav = {'CS', 'pA', 'dprime', 'criterion', 'pC', 'RT'}; nMetrics_behav = length(namesMetrics_behav);
    [d70,~] = SX_sim06_SDT(.7, .3);

    sz_wd_perBar = 200;

    for iSet = 1:numel(iLocSingle_allSets)

        iLocSingle_perSet = iLocSingle_allSets{iSet};

        nBars = numel(iLocSingle_perSet);
        sz_fig = [nBars*sz_wd_perBar, 500];

        for iMetric_prob = 1:nMetrics_behav
            switch iMetric_prob
                case 1, iMetric_vec = 10; x_ticks = linspace(1.6, 3.6, 5); flag_plotIDVD = 1; ref=nan;
                case 2, iMetric_vec = 6; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 1; ref=nan;
                case 3, iMetric_vec = 1; x_ticks = linspace(0, 1.6, 5); flag_plotIDVD = 0; ref=d70;
                case 4, iMetric_vec = 2; x_ticks = linspace(-1,1, 5); flag_plotIDVD = 0; ref=0;
                case 5, iMetric_vec = 3; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 0; ref=nan;
                case 6, iMetric_vec = 11; x_ticks = linspace(0, .2, 5); flag_plotIDVD = 0; ref=nan;
            end
            x_ticks = round(x_ticks, 2);

            data_allIter_allSubj = squeeze(metrics_allCond(iModelA_plot, iLocSingle_perSet, :, :, iDataset_plotNOM, iMetric_vec));

            str_title = sprintf('%s nIter=%d L%s', namesMetrics_behav{iMetric_prob}, nIterxJob, strjoin(string(iLocSingle_perSet), ''));
            %------------------------------%
            basicFxn_drawBars_permutation(data_allIter_allSubj, ref, colors_comb(iLocSingle_perSet, :), namesLocComb(iLocSingle_perSet), x_ticks, x_ticks, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj)
            % ------------------------------%
            ylabel(namesMetrics_behav{iMetric_prob})

            saveas(gcf, sprintf('%s/n%d_L%s_%s.png', nameFolder_Fig_behav, nsubj, strjoin(string(iLocSingle_perSet), ''), namesMetrics_behav{iMetric_prob}))
            close(gcf)
        end % iMetric
    end % iSet
    clear metrics_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% Behavioral metrics: paired locations
    clc, fprintf('\n\n 2/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'metrics_allCond')

    sz_wd_perBar = 200;
    nBars = 2;
    sz_fig = [nBars*sz_wd_perBar, 300];

    for iGroup = 1:nGroups
        iLocPair_all = iLocGroups_all{iGroup};
        for iMetric_prob = 1:nMetrics_behav
            switch iMetric_prob
                case 1, iMetric_vec = 10; x_ticks = linspace(1.6, 3.6, 5); flag_plotIDVD = 1; ref=nan;
                case 2, iMetric_vec = 6; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 1; ref=nan;
                case 3, iMetric_vec = 1; x_ticks = linspace(0, 1.6, 5); flag_plotIDVD = 0; ref=d70;
                case 4, iMetric_vec = 2; x_ticks = linspace(-1,1, 5); flag_plotIDVD = 0; ref=0;
                case 5, iMetric_vec = 3; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 0; ref=.7;
                case 6, iMetric_vec = 11; x_ticks = linspace(0, .2, 5); flag_plotIDVD = 0; ref=nan;
            end
            x_ticks = round(x_ticks, 2);

            data_allIter_allSubj = squeeze(metrics_allCond(iModelA_plot, iLocPair_all, :, :, iDataset_plotNOM, iMetric_vec));

            str_title = sprintf('%s nIter=%d L%d%d', namesMetrics_behav{iMetric_prob}, nIterxJob, iLocPair_all);
            %------------------------------%
            basicFxn_drawBars_permutation(data_allIter_allSubj, ref, colors_comb(iLocPair_all, :), namesLocComb(iLocPair_all), x_ticks, x_ticks, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj)
            %------------------------------%
            saveas(gcf, sprintf('%s/n%d_L%d%d_%s.png', nameFolder_Fig_behav, nsubj, iLocPair_all, namesMetrics_behav{iMetric_prob}))
            close(gcf)
        end % iMetric
    end % iGroup
    clear metrics_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% Plot templates for IDVD and group averages
    clc, fprintf('\n\n 3/24 Plotting STARTED......\n\n')
    % Load data
    load(nameFolder_Data_SaveCompile, 'template_*_allCond')

    pixMin = 0;
    pixMax = 0;
    % Define folder for saving figures
    nameFolder_Fig_NOM_Template = sprintf('%s/Template', nameFolder_Fig_NOM_Trialwise);
    if isempty(dir(nameFolder_Fig_NOM_Template)), mkdir(nameFolder_Fig_NOM_Template), end

    for iDataset = 1:nDatasets
        switch iDataset % needs to matchOOD_xx_compIV ("for iDataset = 1:2")
            case 1
                template_allSubj = squeeze(template_tmpl_allCond(iModelA_plot, iModelB_plot, :, :, :, :)); % ... nLoc x nSubj x nORI x nSF
            case 2
                template_allSubj = squeeze(template_full_allCond(iModelA_plot, iModelB_plot, :, :, :, :)); % ... nLoc x nSubj x nORI x nSF
        end
        str_dataset = namesDataset{iDataset};

        %%%% Group ave %%%%%
        for iLocComb = iLocComb_all
            figure('Position', [0 0 1e3 1e3]), hold on
            e2D_ave = getCI(template_allSubj(iLocComb, :, :, :), 2, 2);
            e2D_ave = e2D_ave';
            fprintf('\n%s L%d: Min = %.3f, Max=%.3f\n', str_dataset, iLocComb, min(e2D_ave(:)), max(e2D_ave(:)))

            if iLocComb==1
                cLim = [-.03, .26];
            else
                cLim = [-.02, .14];
                pixMin = min(pixMin, min(e2D_ave(:)));
                pixMax = max(pixMax, max(e2D_ave(:)));
            end
            RCplot_2Dkernel(e2D_ave, cLim)
            title(sprintf('n=%d nIter=%d L%d %s (%s)', nsubj, nIterxJob, iLocComb, namesLocComb{iLocComb}, str_dataset))
            saveas(gcf, sprintf('%s/n%d_group_L%d_A%d_%s.png', nameFolder_Fig_NOM_Template, nsubj, iLocComb, iModelA_plot, str_dataset))
            close(gcf)
        end % iiLoc

        %%%%% Idvd data in one figure, per loc %%%%%
        % for iiLoc = 1:nLocComb8
        %
        %     figure('Position', [0, 0, 2e3, 1.8e3])
        %     for isubj = 1:nsubj
        %
        %         e2D = squeeze(template_allSubj(iLocComb_all(iiLoc), isubj, :,:))';
        %
        %         subplot(nRows, nCols, isubj), hold on
        %         RCplot_2Dkernel(e2D)
        %         title(subjList{isubj})
        %
        %     end % isubj
        %     set(findall(gcf, '-property', 'FontSize'), 'FontSize',12)
        %     sgtitle(sprintf('L%d %s [A%dB%d] (%s)', iLocComb_all(iiLoc), namesLocComb{iLocComb_all(iiLoc)}, iModelA_plot, iModelB_plot, str_title), 'FontSize',20)
        %
        %     % save
        %     saveas(gcf, sprintf('%s/n%d_L%d_A%dB%d_%s.png', nameFolder_Fig_NOM_Template, nsubj, iLocComb_all(iiLoc), iModelA_plot, iModelB_plot, str_title))
        %     close(gcf)
        %
        % end % iiLoc
    end % iDataset

    % Print min and max pixel value
    fprintf('\n Summary: Min = %.3f, Max=%.3f\n', pixMin, pixMax)

    clear template_*_allCond
    fprintf('\n\n Plotting DONE\n\n')


    %% Plot separability
    clc; fprintf('\n\n 4/24 Plotting STARTED......\n\n')
    str_dataset = namesDataset{iDataset_plotRC};

    % Load compiled separability: sep_allCond
    load(nameFolder_Data_SaveCompile, 'sep_allCond')

    % Folder for saving figures
    nameFolder_Fig_Sep = sprintf('%s/Separability', nameFolder_Fig_NOM_Trialwise);
    if ~exist(nameFolder_Fig_Sep, 'dir'), mkdir(nameFolder_Fig_Sep); end

    sz_labelOffset = 0.05;
    sz_axOffset = 0.05; % extra breathing room

    % ---- Loop across sets of locations ----
    for iSet = 1:numel(iLocSingle_allSets)

        % Locations included in this figure
        iLocSingle_perSet = iLocSingle_allSets{iSet};
        nLoc = numel(iLocSingle_perSet);

        % Extract separability for this set
        % sep_allIter_allSubj: [nLoc x nSubj x nIter]
        sep_allIter_allSubj = squeeze(sep_allCond(iModelA_plot, iModelB_plot, iLocSingle_perSet, :, :, iDataset_plotRC));

        % Reformat to match basicFxn_drawBars_permutation:
        % data_allIter_allSubj: [nIter x nSubj x nCond], where nCond = nLoc
        data_allIter_allSubj = permute(sep_allIter_allSubj, [3 2 1]);  % [nIter x nSubj x nLoc]

        % Colors for each bar/condition (each location)
        colors = colors_comb(iLocSingle_perSet, :); % [nLoc x 3]

        % Condition labels (not displayed on x-axis inside the function, but used in title strings)
        x_ticks = cell(1, nLoc);
        for iiLoc = 1:nLoc
            x_ticks{iiLoc} = sprintf('L%d', iLocSingle_perSet(iiLoc));
        end

        % Plot settings for this figure
        ref = nan;                 % reference line used in your old plot
        y_ticks = [0.6 0.7 0.8 0.9 1.0];
        y_ticklabels = y_ticks;

        flag_plotIDVD = 1;         % show subject-level lines (recommended)
        flag_plotDiff = 0;         % only meaningful when nCond==2 (turn on if you want)
        sz_fig = [600 600];

        str_loc = strjoin(string(iLocSingle_perSet), '');
        str_title = sprintf('n%d nIter=%d [A%dB%d] L%s %s', nsubj, nIterxJob, iModelA_plot, iModelB_plot, str_loc, str_dataset);

        % -------------------%
        basicFxn_drawBars_permutation(data_allIter_allSubj, ref, colors, x_ticks, y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj);
        % -------------------%
        ylabel('Separability (Pearson''s r)')

        % Adjust distance between components
        ax = gca;
        ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
        ax.XLabel.Units = 'normalized';
        ax.YLabel.Units = 'normalized';

        ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
        ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left

        ax.Position = [ ...
            ti(1) + sz_axOffset, ...
            ti(2) + sz_axOffset, ...
            1 - ti(1) - ti(3) - 2*sz_axOffset, ...
            1 - ti(2) - ti(4) - 2*sz_axOffset];

        % ---- save ----
        outName = sprintf('%s/n%d_L%s_A%d_%s.png', nameFolder_Fig_Sep, nsubj, str_loc, iModelA_plot, str_dataset);
        saveas(gcf, outName);
        close(gcf);

    end % iSet

    clear sep_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% Tuning functions: group averages
    clc, fprintf('\n\n 5/24 Plotting STARTED......\n\n')
    % Load data
    load(nameFolder_Data_SaveCompile, 'marg*_allCond', 'margPred*_allCond', 'margR2*_allCond')

    % Define folder for saving figures
    nameFolder_Fig_NOM_Tuning = sprintf('%s/TuningFxns_group', nameFolder_Fig_NOM_Trialwise);
    if isempty(dir(nameFolder_Fig_NOM_Tuning)), mkdir(nameFolder_Fig_NOM_Tuning), end

    wd_border = 3.5; % default 5
    sz_ticks = 30;% default 35

    iplots = reshape((1:nGroups*nFeatures)', nGroups, nFeatures)';
    for iDataset=1:nDatasets
        switch iDataset
            case 1
                markerStyle = 'o';
                lineStyle = '-';
            case 2
                markerStyle = 's';
                lineStyle = '--';
        end

        for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
            iLocPair_all = iLocGroups_all{iGroup};

            for iFeature = 1:nFeatures

                xaxis = axis_tuning{iFeature};

                figure('Position', [0 0 1.1e3 6e2]) % default 8e2
                hold on
                for iLoc = iLocPair_all

                    % Set color
                    color_comb = colors_comb(iLoc, :);

                    % Set y-ticks
                    if iFeature==1 % ORI tuning fxn
                        marg_allCond = margORI_allCond;
                        margPred_allCond = margPred_ORI_allCond;
                        margR2_allCond = margR2_ORI_allCond;
                        if find(iLocPair_all==1), yticks_ = [-.03, linspace(0, .24, 4)]; % fov vs. peri, higher ub
                        else, yticks_ = [-.03, linspace(0, .12, 4)];
                        end
                    else % SF tuning fxn
                        marg_allCond = margSF_allCond;
                        margPred_allCond = margPred_SF_allCond;
                        margR2_allCond = margR2_SF_allCond;
                        if find(iLocPair_all==1), yticks_ = [-.01, linspace(0, .12, 4)];
                        else, yticks_ = [-.01, linspace(0, .06, 4)]; %[-.02, 0, .02, .04, .06];
                        end
                    end

                    ymax = max(yticks_);
                    ymin = min(yticks_);

                    % Obtain data and prediction
                    [marg_ave, ~, ~, marg_sem] = getCI(getCI(marg_allCond(iModelA_plot, iModelB_plot, iLoc, :, :, iDataset, :), 1, 5), 2, 1);
                    [margPred_ave, ~, ~, margPred_sem] = getCI(getCI(margPred_allCond(iModelA_plot, iModelB_plot, iLoc, :, :, iDataset, :), 1, 5), 2, 1);
                    [R2_ave, ~, ~, R2_sem] = getCI(getCI(margR2_allCond(iModelA_plot, iModelB_plot, iLoc, :, :, iDataset), 1, 5), 2, 1);

                    % Data (dots + errorbars)
                    errorbar(xaxis, marg_ave, marg_sem, markerStyle, 'Color', color_comb, 'CapSize',0)
                    plot(xaxis, marg_ave, 'o', 'color', color_comb, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'HandleVisibility','off')

                    % Prediction (lines + bands)
                    patch([xaxis, flip(xaxis)], [margPred_ave-margPred_sem, flip(margPred_ave+margPred_sem)], color_comb, 'FaceAlpha', .3, 'linestyle', 'none')
                    plot(xaxis, margPred_ave, '-', 'color', color_comb)

                    % Draw reference lines
                    yline(0, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);
                    xline(iFeature-1, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);

                    % Set ticks, labels and limits (for EACH loc, to print R2 at the right loc)
                    xlabel(namesFeature_axis{iFeature})
                    ylabel(namesFeature_axis_Tuning{iFeature})

                    xlim(axisLim{iFeature})
                    ylim([ymin, ymax])
                    yticks(yticks_)
                    %     xlabel(xlabels_tuning{ifeature}, 'FontSize', sz_label)
                    xticks(axisTicks_tuning{iFeature})
                    xticklabels(axisTL_tuning{iFeature})
                    if iFeature==2, xticklabels(round(axisTL_tuning{iFeature}, 2)), end

                    % Print R2 in the figure
                    xlims = axisLim{iFeature};
                    switch iFeature
                        case 1, x_R2 = xlims(1) + 0.12 * range(xlims);   % slightly right of left boundary
                        case 2, x_R2 = xlims(2) - 0.3 * range(xlims);   % slightly left of rightboundary
                    end
                    y_R2 = ymax - 0.01-0.08*(find(iLoc == iLocPair_all)-1) * (ymax-ymin);   % slightly below top boundary

                    text(x_R2, y_R2, ...
                        sprintf('$R^2 = %.2f \\pm %.2f$', R2_ave, R2_sem), ...
                        'Interpreter', 'latex', ...
                        'FontSize', 30, ...
                        'Color', color_comb, ...
                        'HorizontalAlignment', 'left', ...
                        'VerticalAlignment', 'top');
                end % iLoc

                ax = gca;
                ax.XAxis.FontSize = sz_ticks;
                ax.YAxis.FontSize = sz_ticks;
                ax.LineWidth = wd_border/1.5;

                [R2_ave, ~, ~, R2_sem] = getCI(getCI(margR2_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset), 1, 5), 2, 2);

                title(sprintf('n=%d L%d%d [A%dB%d] %s tuning (%s)', nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, namesDataset{iDataset}))

                set(findall(gcf, '-property', 'linewidth'), 'linewidth', 2)

                % Save the figure
                saveas(gcf, sprintf('%s/n%d_L%d%d_A%d_%s_%s.png', nameFolder_Fig_NOM_Tuning, nsubj, iLocPair_all, iModelA_plot, namesFeature{iFeature}, namesDataset{iDataset}))
                close(gcf)
            end % iFeature
        end % iGroup
    end % iDataset = 1:2
    clear marg*_allCond margPred*_allCond margR2*_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% Tuning functions: idvd (LARGE FIGURES!!)
    clc, fprintf('\n\n 6/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'marg*_allCond', 'margPred*_allCond', 'margParams*_allCond')

    % Define folder for saving figures
    nameFolder_Fig_NOM_Tuning = sprintf('%s/TuningFxns_IDVD', nameFolder_Fig_NOM_Trialwise);
    if isempty(dir(nameFolder_Fig_NOM_Tuning)), mkdir(nameFolder_Fig_NOM_Tuning), end

    wd_border = 4; % default 5
    sz_ticks = 30;% default 35

    for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
        iLocPair_all = iLocGroups_all{iGroup};

        for iFeature = 2%1:nFeatures

            xaxis = axis_tuning{iFeature};

            for iDataset = 1:nDatasets
                switch iDataset
                    case 1 % TmplSet
                        markerStyle = 'o';
                        lineStyle = '-';
                    case 2 % FullSet
                        markerStyle = 's';
                        lineStyle = '--';
                end

                figure('Position', [0 0 2e3 1.5e3])

                for isubj = 1:nsubj
                    subplot(nRows_subj, nCols_subj, isubj), hold on

                    subjName = subjList{isubj};

                    str_TunParams ='';

                    for iLoc = iLocPair_all
                        str_TunParams = sprintf('%s\nL%d', str_TunParams, iLoc);

                        % fprintf('\nL%d...', iLoc)
                        % Set color
                        color_comb = colors_comb(iLoc, :);

                        % Set y-ticks
                        if iFeature==1
                            marg_allCond = margORI_allCond;
                            margPred_allCond = margPred_ORI_allCond;
                            margParams_allCond = margParams_ORI_allCond;
                            if iLoc==1, yticks_ = [-.03, 0, .05, .10, .15]; % fov vs. peri, higher ub
                            else, yticks_ = [-.02, linspace(0, .12, 4)];
                            end
                        else
                            marg_allCond = margSF_allCond;
                            margPred_allCond = margPred_SF_allCond;
                            margParams_allCond = margParams_SF_allCond;
                            if iLoc==1, yticks_ = [-.01, linspace(0, .08, 4)];
                            else, yticks_ = [-.01, linspace(0, .06, 4)]; %[-.02, 0, .02, .04, .06];
                            end
                        end

                        ymax = max(yticks_);
                        ymin = min(yticks_);

                        % Obtain data, prediction, and params
                        [marg_med, ~, ~, marg_lb, marg_ub] = getCI(marg_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :), 1, 5);
                        [margPred_med, margPred_lb, margPred_ub] = getCI(margPred_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :), 1, 5);
                        [margParam_med, margParam_lb, margParam_ub] = getCI(margParams_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :), 1, 5);

                        % Data (dots + errorbars)
                        errorbar(xaxis, marg_med, marg_lb, marg_ub, markerStyle, 'Color', color_comb, 'CapSize',0)
                        plot(xaxis, marg_med, markerStyle, 'color', color_comb, 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'HandleVisibility','off')

                        % Prediction (lines + bands)
                        patch([xaxis, flip(xaxis)], [margPred_lb', flip(margPred_ub')], color_comb, 'FaceAlpha', .3, 'linestyle', 'none', 'linewidth', 2)
                        plot(xaxis, margPred_med, '-', 'color', color_comb)

                        % Print estimated parameters
                        % each loc has 3-4 params, think of how to arrange (maybe drop CI)
                        for iTunParam = 1:length(margParam_med)
                            str_TunParams = sprintf('%s | %.2f', str_TunParams, margParam_med(iTunParam));
                        end % iTunParam
                        % str_TunParams = sprintf('%s\n', str_TunParams);
                    end % iLoc

                    % Draw reference lines
                    yline(0, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);
                    xline(iFeature-1, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);

                    xlabel('Orientation (deg)')
                    ylabel('Amplitude')
                    title(sprintf('%s\n%s\n', subjName, str_TunParams))
                    % legend('Location', 'best')

                end % isubj
                sgtitle(sprintf('n=%d L%d%d [A%d] %s tuning (%s)', nsubj, iLocPair_all,iModelA_plot, namesFeature{iFeature}, namesDataset{iDataset}))
                set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10)

                % Save the figure
                saveas(gcf, sprintf('%s/n%d_L%d%d_A%d_%s_%s.png', nameFolder_Fig_NOM_Tuning, nsubj, iLocPair_all, iModelA_plot, namesFeature{iFeature}, namesDataset{iDataset}))
                close(gcf)
            end % i=1:2
        end % iFeature
    end % iGroup
    clear marg*_allCond margPred*_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% Tuning characteristics
    clc, fprintf('\n\n 7/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'margTunC_*_allCond')

    % Define folder for saving figures
    nameFolder_Fig_tunC = sprintf('%s/TuningCs', nameFolder_Fig_NOM_Trialwise);
    if isempty(dir(nameFolder_Fig_tunC)), mkdir(nameFolder_Fig_tunC), end

    flag_plotDist = 0;
    flag_plotIDVD = 1;
    flag_plotDiff = 1;
    paramMode = 2;
    sz_fig = [5e2 5e2]; % size of the figure canvas

    for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
        iLocPair_all = iLocGroups_all{iGroup};

        % Plotting settings
        colors = colors_comb(iLocPair_all, :);
        x_ticks = namesLocComb(iLocPair_all);

        % Loop through each feature
        for iFeature = 1:nFeatures

            % Define the x-axis and parameters for the current feature
            xaxis = axis_tuning{iFeature};
            nfilters = length(xaxis);
            iFamily = iFamily_perF(iFeature);
            namesTunC = namesTunC_unit_perF{iFamily, paramMode};
            nTunC_full = length(namesTunC);

            y_ticks_all = [];

            % Obtain y_ticks of ORI/SF domain
            switch flag_plotDist
                case 0
                    switch iFamily_perF(iFeature)
                        case 1 % ORI tuniningC | Scaled Gaussian
                            if flag_plotIDVD
                                if find(iLocPair_all==1)
                                    y_ticks_all{1} = linspace(0, .4, 5); % ORI peak amp
                                else
                                    y_ticks_all{1} = linspace(0, .2, 5); % ORI peak amp
                                end
                                y_ticks_all{2} = linspace(0, 80, 5); % ORI band
                                y_ticks_all{3} = linspace(-.05, .03, 5); % ORI baseline
                            else
                                if find(iLocPair_all==1)
                                    y_ticks_all{1} = linspace(0, .24, 5); % ORI peak amp
                                else
                                    y_ticks_all{1} = linspace(0, .12, 5); % ORI peak amp
                                end

                                y_ticks_all{2} = linspace(0, 100, 5); % ORI band
                                y_ticks_all{3} = linspace(-.08, .08, 5); % ORI baseline
                            end

                        case 2 % SF tuniningC | log parabola
                            if flag_plotIDVD

                                y_ticks_all{1} = linspace(0, 2, 5); % peak SF
                                if find(iLocPair_all==1)
                                    y_ticks_all{2} = linspace(.01, .2, 5); % SF peak amp
                                else
                                    y_ticks_all{2} = linspace(.0, .12, 5); % SF peak amp
                                end
                                y_ticks_all{3} = linspace(.5, 2.5, 5); % SF bandwidth
                                y_ticks_all{4} = linspace(-.1, .1, 5); % SF baseline
                            else
                                y_ticks_all{1} = linspace(0, 2, 5); % peak SF

                                if find(iLocPair_all==1)
                                    y_ticks_all{2} = linspace(.02, .12, 5); % SF peak amp
                                else
                                    y_ticks_all{2} = linspace(.02, .06, 5); % SF peak amp
                                end
                                y_ticks_all{3} = linspace(.3, 2, 5); % SF bandwidth
                                y_ticks_all{4} = linspace(-.08, 0, 5); % SF baseline
                            end
                    end
                case 1 % plot distribution of group averages, not bars
                    switch iFamily_perF(iFeature)
                        case 1 % ORI tuniningC | Scaled Gaussian
                            if find(iLocPair_all==1)
                                y_ticks_all{1} = linspace(0, .3, 5); % ORI peak amp (fovea)
                            else
                                y_ticks_all{1} = linspace(.06, .14, 5); % ORI peak amp
                            end
                            y_ticks_all{2} = linspace(10, 50, 5); % ORI band
                            y_ticks_all{3} = linspace(-.03, .01, 5); % ORI baseline

                        case 2 % SF tuniningC | log parabola
                            y_ticks_all{1} = linspace(log2(1.4), log2(2.8), 5); % peak SF
                            if find(iLocPair_all==1)
                                y_ticks_all{2} = linspace(.01, .13, 5); % SF peak amp (fovea)
                            else
                                y_ticks_all{2} = linspace(.03, .07, 5); % SF peak amp
                            end
                            y_ticks_all{3} = linspace(.5, .9, 5); % SF bandwidth
                            y_ticks_all{4} = linspace(-.04, .04, 5); % SF baseline

                    end
            end

            % Loop through each tuning characteristic
            for iTunC = 1:nTunC_full

                % Compute median for each observer
                % margTunC_ORI_allCond: nModelA x nModelB x nLoc_pair x nsubj x nIter x nDataset x nTunC
                % tunC_allSubj_allIter: nLoc_pair x nsubj x nIter
                switch iFeature
                    case 1
                        data_allIter_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                        data_obs_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, 1, 2, iTunC)); % full data, without resampling
                    case 2
                        data_allIter_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                        data_obs_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, 1, 2, iTunC));
                end

                % Obtain ytick labels and ref
                y_ticklabels = round(y_ticks_all{iTunC}, 2);
                ref = nan;
                switch iFamily_perF(iFeature)
                    case 1, if find(iTunC==3), ref = 0; end

                    case 2 % log parabola
                        switch iTunC
                            case 1, y_ticklabels = round(2.^y_ticks_all{iTunC}, 1); ref = log2(2);
                                data_allIter_allSubj = log2(data_allIter_allSubj); % pref SF
                            case 4, ref = 0; % baseline
                        end
                end

                str_title = sprintf('n=%d nIter=%d L%d%d [A%d] [%s] | %s %s', nsubj, nIterxJob, iLocPair_all, iModelA_plot, namesDataset{iDataset_plotRC}, namesFeature{iFeature}, namesTunC{iTunC});

                switch flag_plotDist
                    case 0
                        %------------------------------%
                        basicFxn_drawBars_permutation(data_allIter_allSubj, ref, colors, x_ticks, y_ticks_all{iTunC}, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj)
                        % ------------------------------%
                        ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily, 2}{iTunC}))
                    case 1
                        %------------------------------%
                        basicFxn_drawDist_permutation(data_allIter_allSubj, ref, colors, x_ticks, y_ticks_all{iTunC}, y_ticklabels, str_title, sz_fig, nIterxJob, nsubj)
                        %------------------------------%
                        xlabel(sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily, 2}{iTunC}))
                        ylabel('Probabillity')
                end
                % Save the figure
                saveas(gcf, sprintf('%s/n%d_L%d%d_A%d_%s%d.png', nameFolder_Fig_tunC, nsubj, iLocPair_all, iModelA_plot, namesFeature{iFeature}, iTunC))
                close(gcf)

            end % end of iTunC
        end % end of iFeature
    end % end of iGroup
    clear margTunC_*_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% Corr0: Corr between CS and tunParams (just to check bound-hitting)
    clc, fprintf('\n\n 8/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'margParams*_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'tunParam';
    sz_label = 55;
    sz_labelOffset = .05;
    sz_axOffset = 0.05; % extra breathing room

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

    flag_zeroMean = 0;
    flag_plotIdvdCI = 0;
    flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles

    for iSet = 1:length(iLocSingle_allSets)

        iLocCorr_all = iLocSingle_allSets{iSet};
        fprintf('   - L%s\n', strjoin(string(iLocCorr_all), ''))

        % Obtain CS (x-axis)
        X_allSubj = CS_allSubj(:, iLocCorr_all);
        x_ticks = linspace(1.5, 3.5, 5); % EE
        X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);

        for iFeature = 1:nFeatures
            iFamily = iFamily_perF(iFeature);
            namesTunCs = namesTunC_unit_perF{iFamily, 1};
            nTunCs_full = length(namesTunCs);

            for iTunC = 1:nTunCs_full
                switch iFeature
                    case 1, NOMp_allIter_allSubj = squeeze(margParams_ORI_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
                    case 2, NOMp_allIter_allSubj = squeeze(margParams_SF_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
                end

                switch iFeature
                    case 1, y_ticks_allTunC_lb = [0, 0, -.1]; y_ticks_allTunC_ub = [.36, 90, .1];% 3 values are ORI peak amplitude, width, baseline
                    case 2, y_ticks_allTunC_lb = [0, 0, 0, -.1]; y_ticks_allTunC_ub = [4, .24, 1, .1]; % 4 values are SF peak, peak amplitude, width, baseline
                end

                % y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);
                y_ticks = [];

                x_ticklabels = x_ticks;
                y_ticklabels = y_ticks;

                nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunCs{iTunC});
                nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

                str_title = sprintf('n=%d, nIter=%d %s vs. %s [L%s]', nsubj, nIterxJob, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

                %----------------------------%
                basicFxn_drawCorr_permutation(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, nIterxJob);
                %----------------------------%

                % plot ub and lb when fitting
                yline(ub_full_all{iFamily}(iTunC), 'k--');
                yline(lb_full_all{iFamily}(iTunC), 'k--');

                xlabel('Contrast sensitivity (1/contrast)', 'fontsize', sz_label)
                ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunCs{iTunC}), 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left

                ax.Position = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];

                saveas(gcf, sprintf('%s/n%d_L%s_%s.png', nameFolder_Fig_NOM_corr, nsubj, strjoin(string(iLocCorr_all), ''), nameVarY_fileTitle))
                close(gcf)

            end % iTunC
        end % iFeature
    end % iSet
    clear margParams*_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% Corr1: Corr between CS and tunC
    clc, fprintf('\n\n 9/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'margTunC*_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'tunC';
    sz_label = 50;
    sz_labelOffset = .03;
    sz_axOffset = 0.05; % extra breathing room

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

    flag_zeroMean = 0;
    flag_plotIdvdCI = 0;
    flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles

    for iSet = 1:numel(iLocSingle_allSets)

        iLocCorr_all = iLocSingle_allSets{iSet};
        fprintf('   - L%s\n', strjoin(string(iLocCorr_all), ''))

        % Obtain CS (x-axis)
        X_allSubj = CS_allSubj(:, iLocCorr_all);
        x_ticks = linspace(1.5, 3.5, 5); % EE
        X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);

        for iFeature = 1:nFeatures
            iFamily = iFamily_perF(iFeature);
            namesTunCs = namesTunC_unit_perF{iFamily, 2};
            nTunCs_full = length(namesTunCs);

            for iTunC = 1:nTunCs_full
                switch iFeature
                    case 1, NOMp_allIter_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
                    case 2, NOMp_allIter_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
                end

                switch iFeature
                    case 1, y_ticks_allTunC_lb = [0, 0, -.1]; y_ticks_allTunC_ub = [.2, 80, .1];% 3 values are ORI peak amplitude, width, baseline
                    case 2, y_ticks_allTunC_lb = [0, 0, 0, -.1]; y_ticks_allTunC_ub = [4, .12, 2.8, .1]; % 4 values are SF peak, peak amplitude, width, baseline
                end

                y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);

                x_ticklabels = x_ticks;
                y_ticklabels = y_ticks;

                nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily_perF(iFeature), 2}{iTunC});
                nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

                str_title = sprintf('n=%d, nIter=%d %s vs. %s [L%s]', nsubj, nIterxJob, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

                %======================%
                basicFxn_drawCorr_permutation(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, nIterxJob);
                %======================%

                xlabel('Contrast sensitivity (1/contrast)', 'fontsize', sz_label)
                ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily_perF(iFeature), 2}{iTunC}), 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left

                ax.Position = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];

                saveas(gcf, sprintf('%s/n%d_L%s_%s.png', nameFolder_Fig_NOM_corr, nsubj, strjoin(string(iLocCorr_all), ''), nameVarY_fileTitle))
                close(gcf)

            end % iTunC
        end % iFeature
    end % iSet
    clear margTunC*_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% Corr2: Corr between CS and pA
    clc, fprintf('\n\n 10/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'metric_data_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'pA'; iMetric_pA = 2;

    sz_label = 55;
    sz_labelOffset = .05;
    sz_axOffset = 0.05; % extra breathing room

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

    flag_zeroMean = 0;
    flag_plotIdvdCI = 0;
    flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles    

    for iSet = 1:numel(iLocSingle_allSets)

        iLocCorr_all = iLocSingle_allSets{iSet};
        fprintf('   - L%s\n', strjoin(string(iLocCorr_all), ''))

        % Obtain CS (x-axis)
        X_allSubj = CS_allSubj(:, iLocCorr_all);
        x_ticks = linspace(1.5, 3.5, 5); % EE
        X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);

        NOMp_allIter_allSubj = getCI(metric_data_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, iMetric_pA, :, :, :), 2, 7);

        y_ticks = round(linspace(.6, .9, 5), 2);

        x_ticklabels = x_ticks;
        y_ticklabels = y_ticks;

        nameVarY_figTitle = nameVarY;
        nameVarY_fileTitle = nameVarY_figTitle;

        str_title = sprintf('n=%d, nIter=%d%s vs. %s [L%s]', nsubj, nIterxJob, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

        basicFxn_drawCorr_permutation(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, nIterxJob);

        xlabel('Contrast sensitivity (1/contrast)', 'fontsize', sz_label)
        ylabel(nameVarY_figTitle, 'fontsize', sz_label)

        % Adjust distance between components
        ax = gca;
        ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
        ax.XLabel.Units = 'normalized';
        ax.YLabel.Units = 'normalized';

        ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
        ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left

        ax.Position = [ ...
            ti(1) + sz_axOffset, ...
            ti(2) + sz_axOffset, ...
            1 - ti(1) - ti(3) - 2*sz_axOffset, ...
            1 - ti(2) - ti(4) - 2*sz_axOffset];

        saveas(gcf, sprintf('%s/n%d_L%s.png', nameFolder_Fig_NOM_corr, nsubj, strjoin(string(iLocCorr_all), '')))
        close(gcf)

    end % iSet
    clear metric_data_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% CompAsym1: CS and tunC: compare binned estimates
    clc, fprintf('\n\n 11/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'margTunC_*_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'tunC';

    nameFolder_Fig_NOM_CompAsym = sprintf('%s/CompAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_CompAsym)); mkdir(nameFolder_Fig_NOM_CompAsym), end

    sz_label = 15;
    sz_labelOffset = 0.01;
    sz_axOffset = 0.05; % extra breathing room
    nBinsCompAsym = 2;
    sz_wd_perBar = 200;
    sz_fig = [1e3, 500];

    for iGroup = 1:nGroups
        iLocPair_all = iLocGroups_all{iGroup};
        fprintf('   - L%s\n', strjoin(string(iLocPair_all), ''))

        switch iLocPair_all(1)
            case 1, nameAsymX = 'Ecc. effect'; nameAsymY = 'Ecc. effect';
            case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
            case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
        end

        % X-axis
        asymX_allSubj = (CS_allSubj(:, iLocPair_all(1))-CS_allSubj(:, iLocPair_all(2)))./(CS_allSubj(:, iLocPair_all(1))+CS_allSubj(:, iLocPair_all(2)));
        switch iLocPair_all(1)
            case 1, x_ticks = linspace(0, 20, 5); % EE
            case 6, x_ticks = linspace(-5, 15, 5); % HVA
            case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
        end
        asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob); % [nsubj x nIterxJob]

        for iFeature = 1:nFeatures
            iFamily = iFamily_perF(iFeature);
            namesTunCs = namesTunC_unit_perF{iFamily, 2};
            nTunCs_full = length(namesTunCs);

            for iTunC = 1:nTunCs_full
                % Y-axis
                switch iFeature
                    case 1, NOMp_allIter_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                    case 2, NOMp_allIter_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                end

                % Y-value: [nsubj x nIterxJob]
                asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

                switch iLocPair_all(1)
                    case 1 % EE
                        switch iFamily
                            case 1, y_ticks_allTunC_lb = -[0, 32, 100]; y_ticks_allTunC_ub = [72, 20, 100];
                                % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                            case 2, y_ticks_allTunC_lb = -[50, 0, 50, 145]; y_ticks_allTunC_ub = [50, 80, 50, 265];
                        end
                    case 6 % HVA
                        switch iFamily
                            case 1, y_ticks_allTunC_lb = -[40, 30, 100]; y_ticks_allTunC_ub = [40, 50, 100];
                                % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                            case 2, y_ticks_allTunC_lb = -[50, 20, 60, 145]; y_ticks_allTunC_ub = [50, 60, 60, 265];
                        end
                    case 5 % VMA
                        switch iFamily
                            case 1, y_ticks_allTunC_lb = -[30, 70, 50]; y_ticks_allTunC_ub = [90, 30, 50];
                                % case 8, y_ticks_allTunC_lb = -[100, 40, 300, 320, 50, 320]; y_ticks_allTunC_ub = [100, 44, 300, 320, 50, 320];
                            case 2, y_ticks_allTunC_lb = -[50, 40, 60, 250]; y_ticks_allTunC_ub = [50, 100, 40, 250];
                        end
                end
                y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);

                x_ticklabels = nan;
                y_ticklabels = nan;

                nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily_perF(iFeature), 2}{iTunC});
                nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

                str_title = sprintf('n=%d, nIter=%d %s (%s) vs. %s (%s)', nsubj, nIterxJob, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

                % remove later!!s
                asymX_allIter_allSubj_ = repmat(asymX_allIter_allSubj, 5, 1);
                asymX_allIter_allSubj_ = asymX_allIter_allSubj_+randn(size(asymX_allIter_allSubj_))*mean(asymX_allIter_allSubj_(:))/20;
                asymY_allIter_allSubj_ = repmat(asymY_allIter_allSubj, 5, 1);
                asymY_allIter_allSubj_ = asymY_allIter_allSubj_+randn(size(asymY_allIter_allSubj_))*mean(asymY_allIter_allSubj_(:))/20;

                str_ylabel = sprintf('\\Delta %s %s (%%)', namesFeature{iFeature}, namesTunC_noUnit{iFamily_perF(iFeature), 2}{iTunC});
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                basicFxn_compAsym_permutation(asymX_allIter_allSubj_*100, asymY_allIter_allSubj_*100, nBinsCompAsym, y_ticks, sz_fig, str_title, str_ylabel)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % xlabel(sprintf('\\Delta contrast sensitivity (%%)'), 'FontSize', sz_label);
                % ylabel(str_ylabel, 'fontsize', sz_label)

                % Adjust distance between components
                % ax = gca;
                % ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
                %
                % ax.XLabel.Units = 'normalized';
                % ax.YLabel.Units = 'normalized';
                %
                % ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
                % ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left
                % ax.Position = [ ...
                %     ti(1) + sz_axOffset, ...
                %     ti(2) + sz_axOffset, ...
                %     1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                %     1 - ti(2) - ti(4) - 2*sz_axOffset];

                saveas(gcf, sprintf('%s/n%d_L%d%d_%s.png', nameFolder_Fig_NOM_CompAsym, nsubj, iLocPair_all, nameVarY_fileTitle))
                close(gcf)

            end % iTunC
        end % iFeature
    end % iGroup
    clear margTunC_*_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% CorrAsym1: (CS and tunC): corr between extents of EE/HVA/VMA
    flag_plotIdvdCI = 1;
    flag_plotUnikSymbol = 0;

    clc, fprintf('\n\n 12/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'margTunC_*_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'tunC';

    nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

    sz_label = 55;
    sz_labelOffset = 0.01;
    sz_axOffset = 0.05; % extra breathing room

    for iGroup = 1:nGroups
        iLocPair_all = iLocGroups_all{iGroup};
        fprintf('   - L%s\n', strjoin(string(iLocPair_all), ''))

        switch iLocPair_all(1)
            case 1, nameAsymX = 'Ecc. effect'; nameAsymY = 'Ecc. effect';
            case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
            case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
        end

        % X-axis
        asymX_allSubj = (CS_allSubj(:, iLocPair_all(1))-CS_allSubj(:, iLocPair_all(2)))./(CS_allSubj(:, iLocPair_all(1))+CS_allSubj(:, iLocPair_all(2)));
        switch iLocPair_all(1)
            case 1, x_ticks = linspace(0, 20, 5); % EE
            case 6, x_ticks = linspace(-5, 15, 5); % HVA
            case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
        end
        asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';

        for iFeature = 1:nFeatures
            iFamily = iFamily_perF(iFeature);
            namesTunCs = namesTunC_unit_perF{iFamily, 2};
            nTunCs_full = length(namesTunCs);

            for iTunC = 1:nTunCs_full
                % Y-axis
                switch iFeature
                    case 1, NOMp_allIter_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                    case 2, NOMp_allIter_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                end

                asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

                switch iLocPair_all(1)
                    case 1 % EE
                        switch iFamily
                            case 1, y_ticks_allTunC_lb = -[0, 40, 100]; y_ticks_allTunC_ub = [80, 40, 160];
                                % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                            case 2, y_ticks_allTunC_lb = -[10, -20, 50, 180]; y_ticks_allTunC_ub = [50, 70, 50, 200];
                        end
                    case 6 % HVA
                        switch iFamily
                            case 1, y_ticks_allTunC_lb = -[40, 40, 200]; y_ticks_allTunC_ub = [60, 60, 200];
                                % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                            case 2, y_ticks_allTunC_lb = -[50, 40, 80, 160]; y_ticks_allTunC_ub = [50, 60, 80, 200];
                        end
                    case 5 % VMA
                        switch iFamily
                            case 1, y_ticks_allTunC_lb = -[40, 70, 200]; y_ticks_allTunC_ub = [100, 50, 200];
                                % case 8, y_ticks_allTunC_lb = -[100, 40, 300, 320, 50, 320]; y_ticks_allTunC_ub = [100, 44, 300, 320, 50, 320];
                            case 2, y_ticks_allTunC_lb = -[50, 60, 100, 300]; y_ticks_allTunC_ub = [70, 100, 100, 300];
                        end
                end
                y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);

                x_ticklabels = nan;
                y_ticklabels = nan;

                nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily_perF(iFeature), 2}{iTunC});
                nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

                str_title = sprintf('n=%d, nIter=%d %s (%s) vs. %s (%s)', nsubj, nIterxJob, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                basicFxn_drawCorrAsym_permutation(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                xlabel(sprintf('\\Delta contrast sensitivity (%%)'), 'FontSize', sz_label);
                ylabel(sprintf('\\Delta %s %s (%%)', namesFeature{iFeature}, namesTunC_noUnit{iFamily_perF(iFeature), 2}{iTunC}), 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks

                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left
                ax.Position = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];

                saveas(gcf, sprintf('%s/n%d_L%d%d_%s.png', nameFolder_Fig_NOM_CorrAsym, nsubj, iLocPair_all, nameVarY_fileTitle))
                close(gcf)

            end % iTunC
        end % iFeature
    end % iGroup
    clear margTunC_*_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% CorrAsym2: (CS and pA): corr between extents of EE/HVA/VMA
    clc, fprintf('\n\n 13/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'metric_data_allCond', 'CS_allSubj', 'pA_allSubj')

    nameVarX = 'CS';
    nameVarY = 'pA';

    iMetric_pA = 2;

    nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

    sz_label = 60;
    sz_labelOffset = 0.05;
    sz_axOffset = 0.05; % extra breathing room

    for iGroup = 1:nGroups
        iLocPair_all = iLocGroups_all{iGroup};
        fprintf('   - L%s\n', strjoin(string(iLocPair_all), ''))

        switch iLocPair_all(1)
            case 1, nameAsymX = 'EE'; nameAsymY = 'EE';
            case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
            case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
        end

        asymX_allSubj = (CS_allSubj(:, iLocPair_all(1))-CS_allSubj(:, iLocPair_all(2)))./(CS_allSubj(:, iLocPair_all(1))+CS_allSubj(:, iLocPair_all(2)));
        switch iLocPair_all(1)
            case 1, x_ticks = linspace(0, 20, 5); % EE
            case 6, x_ticks = linspace(-5, 15, 5); % HVA
            case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
        end
        asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';

        NOMp_allIter_allSubj = getCI(metric_data_allCond(iModelA_plot, iModelB_plot, iLocPair_all, iMetric_pA, :, :, :), 2, 7);

        asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

        % switch iLocPair_all(1)
        %     case 1 % EE
        %         y_ticks_lb = -[10, 70, 12]; y_ticks_ub = [6, 40, 30];
        %     case 6 % HVA
        %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        %     case 5 % VMA
        %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        % end

        y_ticks = linspace(-10, 10, 5);

        x_ticklabels = nan;
        y_ticklabels = nan;

        nameVarY_figTitle = nameVarY;
        nameVarY_fileTitle = nameVarY_figTitle;

        str_title = sprintf('n=%d, nIter=%d %s (%s) vs. %s (%s)', nsubj, nIterxJob, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        basicFxn_drawCorrAsym_permutation(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        xlabel(sprintf('%s of contrast sensitivity (%%)', nameAsymX), 'fontsize', sz_label)
        ylabel(sprintf('%s of %s (%%)', nameAsymX, nameVarY_figTitle), 'fontsize', sz_label)

        % Adjust distance between components
        ax = gca;
        ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks

        ax.XLabel.Units = 'normalized';
        ax.YLabel.Units = 'normalized';

        ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
        ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left
        ax.Position = [ ...
            ti(1) + sz_axOffset, ...
            ti(2) + sz_axOffset, ...
            1 - ti(1) - ti(3) - 2*sz_axOffset, ...
            1 - ti(2) - ti(4) - 2*sz_axOffset];

        saveas(gcf, sprintf('%s/n%d_L%d%d.png', nameFolder_Fig_NOM_CorrAsym, nsubj, iLocPair_all))
        close(gcf)

    end % iGroup
    clear metric_data_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% [NOM] Metrics vs. DV; group averages
    clc, fprintf('\n\n 14/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'IV_allCond', 'nTrials_allCond', 'metric_*_allCond', 'R2_NOM_allCond')

    % Define folder for saving figures for each model A and model B
    nameFolder_Fig_NOM_metrics = sprintf('%s/NOMmetrics_A%d_group', nameFolder_Fig_NOM_Trialwise, iModelA_plot);
    if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end

    x_base   = 0.60;   % move a bit further left to make space for name column
    x_width  = 0.4;   % total width of the table block
    y_base   = 0.10;
    y_height = 0.25;

    sz_font = 15;
    sz_axis = 18;
    lineStyle_all = {'-', '--', ':', '-.'};
    namesPlotMode = {'Full', 'NoShared', 'NoMulti', 'NoAdd', 'AllModels', 'DataOnly'};
    namesModelB_plot = {'Full model', 'No $\sigma_{shared}$', 'No $\sigma_{mul}$', 'No $\sigma_{add}$'};
    szScaling = 10; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end
    szBase = 5;
    sz_font_num = 15;
    wd = 2;

    for iSet = 1:numel(iLocSingle_allSets)

        iLocSingle_perSet = iLocSingle_allSets{iSet};

        fprintf('\nL%s...', strjoin(string(iLocSingle_perSet), ''))

        for iPlotMode = 1:numel(namesPlotMode)

            % --------------------------
            % Decide what to plot
            % --------------------------
            scaleMode = 1;
            switch iPlotMode
                case {1,2,3,4} % Plot ONE model per figure (B1..B4)
                    iModelB_all_plot = iPlotMode;   % assumes your ModelB indices are 1..4
                    flag_plotPred = true;
                    flag_showR2table = true;
                    if iPlotMode==1, scaleMode = 1;
                    else, scaleMode=2;
                    end

                case 5 % Plot ALL models  in one figure
                    iModelB_all_plot = iModelB_all; % e.g., [1 2 3 4]
                    flag_plotPred = true;
                    flag_showR2table = true;

                case 6 % Plot DATA only (no predictions)
                    iModelB_all_plot = 1;           % use a single source for data
                    flag_plotPred = false;
                    flag_showR2table = false;       % set true if you still want table (would be 1-row)
                otherwise
                    error('Unknown iPlotMode=%d', iPlotMode);
            end

            nModelB_plot = numel(iModelB_all_plot);

            %%% --- PLOT ---
            for iMetric_prob = 1:nMetrics_prob
                figure('Position', [0 0 500/scaleMode 500/scaleMode]); hold on

                % R2 table: rows = plotted ModelB, cols = locations in this set
                R2_tab = nan(numel(iLocSingle_perSet), 2); % ave and sem

                for rModel = 1:nModelB_plot
                    iModelB = iModelB_all_plot(rModel);

                    for iiLoc = 1:numel(iLocSingle_perSet)
                        locID = iLocSingle_perSet(iiLoc);

                        % Compute medians and CIs
                        [IV_ave, ~, ~, IV_SEM] = getCI(getCI(IV_allCond(iModelA_plot, iModelB, locID, :, :, :), 1, 5), 2, 1);
                        [nTrials_ave, ~, ~, nTrials_SEM] = getCI(getCI(nTrials_allCond(iModelA_plot, iModelB, locID, :, :, :), 1, 5), 2, 1);
                        [data_ave, ~, ~, data_SEM] = getCI(getCI(metric_data_allCond(iModelA_plot, iModelB, locID, iMetric_prob, :, :, :), 1, 6), 2, 1);

                        % pred exists but may not be plotted
                        [pred_ave, ~, ~, pred_SEM] = getCI(getCI(metric_pred_allCond(iModelA_plot, iModelB, locID, iMetric_prob, :, :, :), 1, 6), 2, 1);
                        [R2_tab(iiLoc, 1), ~, ~, R2_tab(iiLoc, 2)] = getCI(getCI(R2_NOM_allCond(iModelA_plot, iModelB, locID, iMetric_prob, :, :), 1, 6), 2, 1);

                        % ---------
                        % Styling
                        % ---------
                        cLoc = colors_comb(locID, :);

                        % For "all models together", differentiate models by linestyle/linewidth
                        % For "single model", keep clean.
                        if numel(iModelB_all_plot) > 1
                            ls = lineStyle_all{iModelB};     % you already defined lineStyle_all
                            if iModelB == 1, lw = 3; else, lw = 1.5; end
                        else
                            ls = '-'; % plot solid line when only one model prediction is plotted
                            % ls = lineStyle_all{iModelB};     % you already defined lineStyle_all
                            lw = 3;
                        end
                        % lw = lw*scaleMode;

                        % --------------------------
                        % Plot prediction (optional)
                        % --------------------------
                        if flag_plotPred
                            patch([IV_ave; flip(IV_ave)], [pred_ave-pred_SEM; flip(pred_ave+pred_SEM)], cLoc, 'FaceAlpha', .15, 'LineStyle', 'none', 'HandleVisibility', 'off');
                            plot(IV_ave, pred_ave, 'LineStyle', ls, 'Color', cLoc, 'LineWidth', lw, 'HandleVisibility', 'on');
                        end

                        % --------------------------
                        % Plot measurement (always)
                        % --------------------------
                        errorbar(IV_ave, data_ave, data_SEM, '.', 'vertical', 'CapSize', 0, 'Color', cLoc, 'HandleVisibility', 'off', 'LineWidth', max(lw/1.5, 1));

                        % Plot averaged data of each bin (dot size indicates number of trials)
                        for iBin = 1:nBins
                            plot(IV_ave(iBin), data_ave(iBin), 'o', 'MarkerEdgeColor', cLoc, 'MarkerFaceColor', 'w', 'MarkerSize', nTrials_ave(iBin) / szScaling + szBase, 'LineWidth', max(lw/1.5, 1), 'LineStyle', 'none', 'HandleVisibility', 'off');
                        end
                    end % iiLoc
                end % rModel

                % reference line and y formatting
                yline(.5, '--', 'LineWidth', 2, 'Color', ones(1,3)/2);

                switch iMetric_prob
                    case 1, ylim([0, 1]);    yticks(0:.2:1)
                    case 2, ylim([.45, 1]);  yticks(.5:.1:1)
                end
                if iPlotMode==1
                    ylabel(namesMetrics_prob_full{iMetric_prob})
                end

                if any(iLocSingle_perSet == 1), x_ticks = linspace(0, 180, 5);
                else,                          x_ticks = linspace(0, 80, 5);
                end
                xticks(x_ticks); xlim(x_ticks([1, end]))
                xlabel('Binned decision variable')

                % --------------------------
                % Plot R2
                % --------------------------
                ax = gca;

                % nRows = nModelB_plot;
                nRows = numel(iLocSingle_perSet);

                % --- layout (more horizontal spacing) ---
                % Use a fixed-width "table box" in normalized coordinates.
                % x_width controls spacing between columns; increase it if you want more gap.

                % centers of each column/row
                x_cells = 0.9;
                y_cells = y_base  + ((nRows:-1:1) - 0.5) * (y_height    / nRows);    % 1 x nRows (top to bottom)

                % --- draw numbers ---
                for iRow = 1:nRows
                    str_cell = sprintf('%.0f%%\\pm%.0f%%', 100*R2_tab(iRow, :));
                    text(x_cells, y_cells(iRow), str_cell, ...
                        'Units','normalized', ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', ...
                        'FontSize', sz_font_num, ...
                        'FontWeight','normal', ...
                        'Interpreter','tex', ...   % supports \pm
                        'Color', colors_comb(iLocSingle_perSet(iRow), :));
                end

                % axis cosmetics
                ax = gca;
                ax.FontSize = sz_axis;
                ax.LineWidth = wd;

                title(sprintf('n=%d [A%d] [L%s] [nIter=%d] %s | %s', ...
                    nsubj, iModelA_plot, strjoin(string(iLocSingle_perSet), ''), nIterxJob, ...
                    namesMetrics_prob{iMetric_prob}, namesPlotMode{iPlotMode}), ...
                    'FontSize', 10);

                % Save
                saveas(gcf, sprintf('%s/n%d_L%s_%s_%s.png', ...
                    nameFolder_Fig_NOM_metrics, nsubj, strjoin(string(iLocSingle_perSet), ''), namesMetrics_prob{iMetric_prob}, namesPlotMode{iPlotMode}));
                close(gcf)

            end % iMetric_prob
        end % iPlotMode
    end % iSet

    fprintf('\n\n Plotting DONE\n\n')

    %% [NOM] Metrics vs. DV; per idvd
    clc, fprintf('\n\n 16/24 Plotting STARTED......\n\n')
    lineStyle_all = {'-', '-', '--', ':', '-.', '-', '--', ':', '-.'};

    % iModelA_plot = 1;%iModelA_all % just plot pred of the model using RC-derived template (happens to be yhe best model)
    % iModelA = iModelA_plot;
    iModelB_all_plot = 1;

    % Load data
    load(nameFolder_Data_SaveCompile, 'params_allCond', 'IV_allCond', 'nTrials_allCond', 'metric_data_allCond', 'metric_pred_allCond', 'R2_*NOM_allCond')

    % Define folder for saving figures for each model A and model B
    nameFolder_Fig_NOM_metrics = sprintf('%s/NOMmetrics_A%d_idvd', nameFolder_Fig_NOM_Trialwise, iModelA_plot);
    if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end

    for iLocSingle = 1:nLocComb8
        fprintf('\nL%d...', iLocSingle)
        szScaling = 80; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end

        %%%%%%%%%%%%
        for iMetric_prob = 1:nMetrics_prob
            figure('Position', [0 0 2e3 1.5e3])

            for isubj = 1:nsubj
                subplot(nRows_subj, nCols_subj, isubj), hold on

                subjName = subjList{isubj};

                % Create a string to store R2
                % str_R2 = 'R2: ';
                % str_NRMSE = 'NRMSE: ';
                str_NRMSE = '';

                % Create a string to store parameter estimates for display
                str_est = [];
                params_allCond(isnan(params_allCond))=0;
                for iModelB = iModelB_all_plot
                    vals = strjoin(string(round(getCI(params_allCond(iModelA_plot, iModelB, iLocSingle, isubj, :, :), 1, 5)',2)), ", ");
                    str_est = [str_est, sprintf('B%d: %s\n', iModelB, vals)];
                end

                for iModelB = iModelB_all_plot

                    % Compute medians and CIs
                    [IV_allBins, ~, ~, IV_allBins_neg, IV_allBins_pos] = getCI(IV_allCond(iModelA_plot, iModelB, iLocSingle, isubj, :, :), 1, 5);
                    [nTrials_allBins, nTrials_allBins_lb, nTrials_allBins_ub] = getCI(nTrials_allCond(iModelA_plot, iModelB, iLocSingle, isubj, :, :), 1, 5);
                    [data_allBins, ~, ~, data_allBins_neg, data_allBins_pos] = getCI(metric_data_allCond(iModelA_plot, iModelB, iLocSingle, iMetric_prob, isubj, :, :), 1, 6);
                    [pred_allBins, pred_allBins_lb, pred_allBins_ub] = getCI(metric_pred_allCond(iModelA_plot, iModelB, iLocSingle, iMetric_prob, isubj, :, :), 1, 6);
                    [params_allBins, params_allBins_lb, params_allBins_ub] = getCI(params_allCond(iModelA_plot, iModelB, iLocSingle, isubj, :, :), 1, 5);
                    [R2_NOM_med, R2_NOM_lb, R2_NOM_ub] = getCI(R2_NOM_allCond(iModelA_plot, iModelB, iLocSingle, iMetric_prob, isubj, :), 1, 6, .95);

                    % % Calculate the normalized RMSE
                    % % 0.1 is excellent, 0.3 is decent, >0.5 is poor
                    n = squeeze(nTrials_allCond(iModelA_plot, iModelB, iLocSingle, iMetric_prob, isubj, :, :));
                    p = squeeze(metric_pred_allCond(iModelA_plot, iModelB, iLocSingle, iMetric_prob, isubj, :, :));
                    d = squeeze(metric_data_allCond(iModelA_plot, iModelB, iLocSingle, iMetric_prob, isubj, :, :));

                    % Get weight
                    w = sqrt(n);
                    w=w/mean(w);
                    W = repmat(w(:)', nIterxJob, 1);
                    err2 = (d-p).^2;
                    rmse_w = sqrt(sum(W.*err2, 2, 'omitnan'))./sum(W, 2, 'omitnan');

                    % (2A) Normalize by weighted SD of data per iteration (statistical)
                    ybar_w = sum(W .* d, 2, 'omitnan') ./ sum(W, 2, 'omitnan');          % [nIter x 1]
                    sd_w   = sqrt( sum(W .* (d - ybar_w).^2, 2, 'omitnan') ./ sum(W, 2, 'omitnan') ); % [nIter x 1]
                    nrmse_w_sd = rmse_w ./ sd_w;     % [nIter x 1]
                    [NRMSE_allBins_med_sd, NRMSE_allBins_lb_sd, NRMSE_allBins_ub_sd] = getCI(nrmse_w_sd, 1, 1);

                    % (2B) Normalize by range of data per iteration (robust)
                    range_y = max(d, [], 2, 'omitnan') - min(d, [], 2, 'omitnan');    % [nIter x 1]
                    nrmse_w_range = rmse_w ./ range_y;  % [nIter x 1]
                    [NRMSE_allBins_med_range, NRMSE_allBins_lb_range, NRMSE_allBins_ub_range] = getCI(nrmse_w_range, 1, 1);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Plot prediction
                    if iModelB==1, color_pred = colors_comb(iLocSingle, :); lw = 3;
                    else, color_pred='k'; lw = 1;
                    end
                    patch([IV_allBins; flip(IV_allBins)], [pred_allBins_lb; flip(pred_allBins_ub)], ones(1, 3) / 2, 'FaceAlpha', .3, 'linestyle', 'none', 'handlevisibility', 'off')
                    plot(IV_allBins, pred_allBins, 'lineStyle', lineStyle_all{iModelB}, 'color', color_pred, 'linewidth', lw)

                    % Plot measurement
                    errorbar(IV_allBins, data_allBins, IV_allBins_neg, IV_allBins_pos, '.', 'horizontal', 'CapSize', 0, 'color', colors_comb(iLocSingle, :), 'handlevisibility', 'off', 'linewidth', lw)
                    errorbar(IV_allBins, data_allBins, data_allBins_neg, data_allBins_pos, '.', 'vertical', 'CapSize', 0, 'color', colors_comb(iLocSingle, :), 'handlevisibility', 'off', 'linewidth', lw)

                    yline(.5, 'k--');

                    % Plot averaged data of each bin
                    for iBin = 1:nBins
                        % if isubj>nMarkersMax, facecolor=colors_comb(iLocComb, :); else, facecolor='w'; end
                        facecolor = 'w';
                        plot(IV_allBins(iBin), data_allBins(iBin), 'o', 'markeredgecolor', colors_comb(iLocSingle, :), ...
                            'markerfacecolor', facecolor, 'MarkerSize', nTrials_allBins(iBin) / szScaling + 5, 'LineWidth', 1, 'LineStyle', 'none', 'handlevisibility', 'off')
                    end

                    % xlim([,150])
                    switch iMetric_prob
                        case 1, ylim([0, 1])
                        case 2, ylim([.5, 1])
                    end
                    ylabel(namesMetrics_prob{iMetric_prob})
                    xlabel('Binned decision variable')
                    % if isubj == 1, legend(namesModelB(iModelB_all), 'Location', 'best'), end

                    % Median and CI of GoF (R2/NRMSE)
                    str_NRMSE = sprintf('%s\n[B%d] NRMSE: sd=%.2f, range=%.2f\n', ...
                        str_NRMSE, iModelB, NRMSE_allBins_med_sd, NRMSE_allBins_med_range);
                    
                end % iModelB

                title(sprintf('[%s] %s', subjName, str_NRMSE), 'fontsize', 20)
            end % isubj

            sgtitle(sprintf('n=%d [A%d] [L%d] [nIter=%d] %s', nsubj, iModelA_plot, iLocSingle, nIterxJob, namesMetrics_prob{iMetric_prob}))
            % set(findall(gcf, '-property', 'fontsize'), 'fontsize', 15)

            % Save the figure
            saveas(gcf, sprintf('%s/n%d_%s_L%d_A%d.png', nameFolder_Fig_NOM_metrics, nsubj, namesMetrics_prob{iMetric_prob}, iLocSingle, iModelA_plot))
            close(gcf)
        end % iMetric

        fprintf('DONE\n')
    end % iLocComb

    fprintf('\n\n Plotting DONE\n\n')

    %% [NOM] ANOVA on nLL (ModelA x ModelB x Loc; with vars collapsed)
    % clc, fprintf('\n\n 17/24 Plotting STARTED......\n\n')
    % % Load data
    % load(nameFolder_Data_SaveCompile, 'nLL_allCond')
    %
    % % Define folder for saving GoF figures
    % nameFolder_Fig_NOM_nLL = fullfile(nameFolder_Fig_NOM_Trialwise, 'nLL');
    % if isempty(dir(nameFolder_Fig_NOM_nLL)), mkdir(nameFolder_Fig_NOM_nLL); end
    %
    % nPerm = 1e4;
    % CI_level = .95;
    %
    % % [NOM] Plot comparison when collapsing modelB or location
    % % iModelB_allSets = {[1,3], [1,3,2,4]};
    % iModelB_allSets = {[1,3,2,4]};
    % DimCollapse_all = {'ModelB', 'Loc'};
    % sz_label = 20;
    %
    % for iSetModelB = 1:numel(iModelB_allSets)
    %     iModelB_all = iModelB_allSets{iSetModelB};
    %
    %     for iSetLoc = 1:numel(iLocSingle_allSets)
    %
    %         iLocSingle_perSet = iLocSingle_allSets{iSetLoc};
    %
    %         % 3-way: modelA x moelB x single Loc
    %         iModelA_selected = [1,2];
    %         nLL_ANOVA = squeeze(nLL_allCond(iModelA_selected, iModelB_all, iLocSingle_perSet, :, :));
    %         stats_3way = rm3ANOVA_iter_perm(nLL_ANOVA, nPerm, CI_level);
    %
    %         % 2-way: modelB x single Loc when modelA=1 (data-derived temlate)
    %         iModelA_selected = 1;
    %         nLL_ANOVA = squeeze(nLL_allCond(iModelA_selected, iModelB_all, iLocSingle_perSet, :, :));
    %         [stats_2way, str_ANOVA2] = rm2ANOVA_A1(nLL_ANOVA, nPerm, CI_level);
    %
    %         figure('Position', [0 0 1e3 500])
    %         for iDimCollapse = 1:2
    %             subplot(1,2,iDimCollapse);
    %             DimCollapse = DimCollapse_all{iDimCollapse};
    %
    %             % nBars = numel(iLocSingle_all);
    %             % sz_fig = [nBars*sz_wd_perBar, 250];
    %
    %             % Extract
    %             nLL_collapse = squeeze(mean(squeeze(nLL_allCond(iModelA_selected, iModelB_all, iLocSingle_perSet, :, :)), iDimCollapse, 'omitnan')); % has the shape [nCond x nSubj x nIter]
    %             nLL_collapse = permute(nLL_collapse, [3,2,1]); % so that input has the shape [nIter x nSubj x nCond]
    %
    %             % For each subject and each iteration, subtract the minimum across models
    %             dnLL_collapse = nan(size(nLL_collapse));
    %             for isubj = 1:nsubj
    %                 parfor iIter = 1:nIterxJob
    %                     d = nLL_collapse(iIter, isubj, :);
    %                     d_min = min(d(:));
    %                     dnLL_collapse(iIter, isubj, :) = d - d_min;
    %                 end
    %             end
    %             dnLL_collapse = squeeze(dnLL_collapse);
    %
    %             str_title = sprintf('n=%d nIter=%d [L%s] [ModelB %s] %s collapsed', nsubj, nIterxJob, strjoin(string(iLocSingle_perSet), ''), strjoin(string(iModelB_all), ''), DimCollapse);
    %             ref = nan;
    %             switch DimCollapse
    %                 case 'ModelB'
    %                     colors = colors_comb(iLocSingle_perSet, :);
    %                     x_ticks = namesLocComb(iLocSingle_perSet);
    %                 case 'Loc'
    %                     colors = repmat(linspace(0, .5, numel(iModelB_all))', 1, 3);
    %                     x_ticks = namesModelB(iModelB_all);
    %             end
    %             y_ticks = nan;
    %             y_ticklabels = nan;
    %             sz_fig = nan; % set to be nan to not create a figure inside the fxn
    %             flag_plotIDVD=1;
    %             flag_plotDiff=1;
    %             %------------------------------%
    %             basicFxn_drawBars_permutation(dnLL_collapse, ref, colors, x_ticks, y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj)
    %             %------------------------------%
    %             % ylim(y_lim)
    %             ylabel(sprintf('\\Delta nLL'), 'FontSize', sz_label);
    %
    %         end % iDimCollapse
    %         sgtitle(sprintf('n=%d nIter=%d [L%s] [ModelB %s]\n%s\n%s', nsubj, nIterxJob, strjoin(string(iLocSingle_perSet), ''), strjoin(string(iModelB_all), ''), str_ANOVA2))
    %
    %         % save
    %         saveas(gcf, fullfile(nameFolder_Fig_NOM_nLL, sprintf('n%d_L%s_ModelB%s.png', nsubj, strjoin(string(iLocSingle_perSet), ''), strjoin(string(iModelB_all), ''))));
    %         close(gcf)
    %     end % iSetLoc
    % end % iSetModelB
    %
    % fprintf('\n\n Plotting DONE\n\n')

    %% [NOM] ANOVA on nLL (ModelB x Loc; NO vars collapsed)
    % Load data
    load(nameFolder_Data_SaveCompile, 'nLL_allCond')
    clc, fprintf('\n\n 18/24 Plotting STARTED......\n\n')

    % Define folder for saving GoF figures
    nameFolder_Fig_NOM_nLL = fullfile(nameFolder_Fig_NOM_Trialwise, 'nLL');
    if isempty(dir(nameFolder_Fig_NOM_nLL)), mkdir(nameFolder_Fig_NOM_nLL); end

    iModelA_selected = 1;

    % ModelB ordering + plotting flags
    iModelB_selected = [1,3,2,4];
    flag_plotIDVD = 0;
    flag_plotDiff = 1;
    wd = 2;

    y_ticks = linspace(0, 16, 5);
    sz_label = 20;

    % Loop over locations
    for iLocSingle = 1:nLocComb8

        nBars = numel(iModelB_selected);
        sz_fig = [nBars * 100, 300+nchoosek(nBars, 2)*30];

        % Extract raw nLL: [ModelB x Subj x Iter]
        nLL_allIter = squeeze(nLL_allCond(iModelA_selected, iModelB_selected, iLocSingle, :, :, :));

        % Sanity reshape to [nBars x nSubj x nIter]
        if ndims(nLL_allIter) == 2
            % If no iteration dimension survived, force nIter=1
            nLL_allIter = reshape(nLL_allIter, [nBars, nsubj, 1]);
        elseif ndims(nLL_allIter) == 3
            % assume [nBars x nSubj x nIter] already
        else
            error('Unexpected nLL dimensionality after squeeze: ndims=%d', ndims(nLL_allIter));
        end

        nLL_med = getCI(nLL_allIter, 1, 3);          % [nCond x nSubj]
        nLL_min_perSubj = min(nLL_med, [], 1);        % [1 x nSubj]
        dnLL_med = nLL_med - nLL_min_perSubj;         % [nBars x nSubj], >=0

        % ------------------------------------------------------------
        % basicFxn_drawBars_permutation expects [nIter x nSubj x nCond]
        % We now have only one "iteration" (the median-collapsed value), so set nIter=1.
        % ------------------------------------------------------------
        dnLL_allIter_allSubj = nan(1, nsubj, nBars);
        dnLL_allIter_allSubj(1,:,:) = dnLL_med.';     % transpose -> [nSubj x nBars]

        % Strings / plotting params
        str_title = sprintf('n=%d | Loc L%d | ModelB [%s] | nIter=%d', nsubj, iLocSingle, strjoin(string(iModelB_selected), ' '), nIterxJob);

        ref = nan;

        % ------------------------------%
        basicFxn_drawBars_permutation(dnLL_allIter_allSubj, ref, repmat(colors_comb(iLocSingle, :), nBars, 1), namesModelB(iModelB_selected), y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, 1, markers_allSubj);
        % ------------------------------%
        ylim(y_ticks([1, end]))
        
        ylabel('\Delta nLL', 'FontSize', sz_label);

        % [Plot] planned comparison brackets + CI of mean difference at midpoint
        pairs_bracket = [1 2; 1 3; 1 4];     % requested comparisons

        % dnLL_allIter_allSubj must be [nIter x nSubj x nCond]
        [~, nSubj, nCond] = size(dnLL_allIter_allSubj);

        % (1) Subject-level summary across iterations (median)
        data_med_allSubj = squeeze(median(dnLL_allIter_allSubj, 1, 'omitnan')); % [nSubj x nCond]

        % (2) Bootstrap group mean + CI for bars
        CI_plot = 0.68;         % 68% for plotting
        nBootPlot = 5000;

        ave_allBoot = nan(nBootPlot, nCond);
        for iBoot = 1:nBootPlot
            idx = randi(nSubj, [1 nSubj]);                 % resample subjects with replacement
            ave_allBoot(iBoot,:) = mean(data_med_allSubj(idx,:), 1, 'omitnan');
        end

        % getCI should return median/mean + CI bounds; here we just want point + lb/ub
        [data_ave, data_lb, data_ub, data_sem_neg, data_sem_pos] = getCI(ave_allBoot, 1, 1, CI_plot);

        % (3) Bootstrap CIs for mean differences for specific pairs
        pairs_all = nchoosek(1:nCond, 2);
        nPairs = size(pairs_all, 1);

        % Compute bootstrap for ALL pairs once (so you can reuse)
        diffCond_allBoot = nan(nBootPlot, nPairs);
        for iBoot = 1:nBootPlot
            idx = randi(nSubj, [1 nSubj]);
            Xb = data_med_allSubj(idx,:);  % [nSubj x nCond]
            for iPair = 1:nPairs
                iA = pairs_all(iPair,1);
                iB = pairs_all(iPair,2);
                diffCond_allBoot(iBoot,iPair) = mean(Xb(:,iA) - Xb(:,iB), 'omitnan');
            end % iPair
        end % iBoot
        [diffCond_med, diffCond_lb, diffCond_ub, diffCond_sem_neg, diffCond_sem_pos] = getCI(diffCond_allBoot, 1, 1, 0.95); % use 95% for reporting

        % Set limits
        yl = ylim;
        yMin = yl(1);
        yMax = yl(2);
        yRange = yMax - yMin;

        % Place brackets above the tallest bar+CI
        topData = y_ticks(end);
        if isnan(topData), topData = yMax; end
        yBase = yMin + 0.6*yRange;
        yStep = 0.15 * yRange;   % vertical spacing between brackets

        % Helper: find index in "pairs" for a given (iA,iB)
        getPairIdx = @(iA,iB) find(pairs_all(:,1)==min(iA,iB) & pairs_all(:,2)==max(iA,iB), 1, 'first');

        for iPair = 1:size(pairs_bracket,1)

            iA = pairs_bracket(iPair,1);
            iB = pairs_bracket(iPair,2);

            y = yBase + (iPair-1)*yStep;

            % --- bracket line ---
            plot([iA iB], [y y], 'k-', 'LineWidth', wd, 'HandleVisibility','off');

            % --- text label (left-aligned, above the bracket line) ---
            str_delta = sprintf('\\Delta=%.1f, [%.1f, %.1f]', diffCond_med(iPair), diffCond_lb(iPair), diffCond_ub(iPair));

            xText = iA;                    % left end of bracket
            yText = y + 0.02*yRange;       % a bit above the bracket line (tune 0.02)

            text(xText, yText, str_delta, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 15, ...
                'Color', 'k', ...
                'Interpreter', 'tex', ...
                'HandleVisibility', 'off');
        end % iPair

        saveas(gcf, fullfile(nameFolder_Fig_NOM_nLL, sprintf('n%d_L%d_ModelB%s.png', nsubj, iLocSingle, strjoin(string(iModelB_selected), ''))));
        close(gcf)

    end % iLocSingle

    clear nLL_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% [NOM] Compare NOM parameters across locations
    % Load data
    load(nameFolder_Data_SaveCompile, 'params_allCond')
    clc, fprintf('\n\n 19/24 Plotting STARTED......\n\n')

    % Define folder for saving GoF figures
    nameFolder_Fig_NOM_params = fullfile(nameFolder_Fig_NOM_Trialwise, 'NOMparams');
    if isempty(dir(nameFolder_Fig_NOM_params)), mkdir(nameFolder_Fig_NOM_params); end

    
    nBars = 2;
    sz_fig = [nBars * 200, 300 + nchoosek(nBars, 2)*30];
    flag_plotIDVD = 1;
    flag_plotDiff = 1;

    for iModelB_NOMplot = iModelB_plot_all

        % Number of parameters for this Model B
        nNOMparams = numel(namesModelBparams{iModelB_NOMplot});

        for ii = 1:2
            switch ii
                case 1, iLocG = iLocGroups_all;
                case 2, iLocG = iLocSingle_allSets;
            end
            for iSetLoc = 1:numel(iLocG)

                iLocPair_all = iLocG{iSetLoc};

                for iParam = 1:nNOMparams

                    NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_NOMplot, iLocPair_all, :, :, iParam));

                    x_ticks = namesLocComb(iLocPair_all);
                    switch nNOMparams
                        case 2
                            switch iParam
                                case 1, y_ticks = linspace(NOMp2_lb, 40, 5);
                                case 2, y_ticks = linspace(NOMp3_lb, 40, 5);
                            end
                        case 3
                            switch iParam
                                case 1, y_ticks = linspace(NOMp1_lb, .8, 5);
                                case 2, y_ticks = linspace(NOMp2_lb, 40, 5);
                                case 3, y_ticks = linspace(NOMp3_lb, 40, 5);
                            end
                    end
                    y_ticklabels = round(y_ticks, 2);

                    str_title = sprintf('n=%d nIter=%d L%s [A%dB%d] %s', nsubj, nIterxJob, strjoin(string(iLocPair_all), ''), iModelA_plot, iModelB_NOMplot, namesModelBparams{iModelB_NOMplot}{iParam});
                    %------------------------------%
                    basicFxn_drawBars_permutation(NOMp_allIter_allSubj, ref, colors_comb(iLocPair_all, :), x_ticks, y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj)
                    %------------------------------%
                    ylabel(sprintf('%s', namesModelBparams{iModelB_NOMplot}{iParam}))
                    % Save the figure
                    saveas(gcf, sprintf('%s/n%d_A%dB%d_L%s_NOMp%d.png', nameFolder_Fig_NOM_params, nsubj, iModelA_plot, iModelB_NOMplot, strjoin(string(iLocPair_all), ''), iParam))
                    close(gcf);
                end % iParam
            end % iGroup
        end % ii
    end % iModelB_NOMplot
    clear params_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% [NOM] Corr1: pA and NOM params
    clc, fprintf('\n\n 20/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'params_allCond', 'pA_allSubj')
    nameVarX = 'pA';
    nameVarY = 'NOMparams';

    sz_label = 55;
    sz_labelOffset = .05;
    sz_axOffset = 0.05; % extra breathing room

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

    flag_zeroMean = 0;
    flag_plotIdvdCI = 0;
    flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles

    

    for iModelB_NOMplot = iModelB_plot_all
        nNOMparams = length(namesModelBparams{iModelB_NOMplot});

        for iSet = 1:numel(iLocSingle_allSets)

            iLocCorr_all = iLocSingle_allSets{iSet};
            fprintf('   - L%s\n', strjoin(string(iLocCorr_all), ''))

            % Obtain x-axis values
            X_allSubj = pA_allSubj(:, iLocCorr_all);
            x_ticks = linspace(.6, .8, 5); % EE
            X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);

            for iNOMparam = 1:nNOMparams
                NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_NOMplot, iLocCorr_all, :, :, iNOMparam));

                switch iModelB_NOMplot
                    case 1, y_ticks_lb = [0, 5, 0]; y_ticks_ub = [.8, 45, 40];
                    case 3, y_ticks_lb = [5, 0]; y_ticks_ub = [45, 40];
                end

                y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

                x_ticklabels = x_ticks;
                y_ticklabels = y_ticks;

                nameVarY_figTitle = namesModelBparams{iModelB_NOMplot}{iNOMparam};
                nameVarY_fileTitle = nameVarY_figTitle;

                str_title = sprintf('n=%d, nIter=%d B%d %s vs. %s [L%s]', nsubj, nIterxJob, iModelB_NOMplot, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

                basicFxn_drawCorr_permutation(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, nIterxJob);

                xlabel(nameVarX, 'fontsize', sz_label)
                ylabel(nameVarY_figTitle, 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left

                ax.Position = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];

                saveas(gcf, sprintf('%s/n%d_B%d_L%s_NOMp%d.png', nameFolder_Fig_NOM_corr, nsubj, iModelB_NOMplot, strjoin(string(iLocCorr_all), ''), iNOMparam))
                close(gcf)
            end % iNOMparam
        end % iSet
    end % iModelB_NOMplot
    clear params_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% [NOM] CorrAsym1: pA and NOMparams: corr between extents of EE/HVA/VMA
    clc, fprintf('\n\n 21/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'params_allCond', 'pA_allSubj')

    nameVarX = 'pA';
    nameVarY = 'NOMparams';

    nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

    sz_label = 60;
    sz_labelOffset = 0.05;
    sz_axOffset = 0.05; % extra breathing room

    for iModelB_NOMplot = iModelB_plot_all
        nNOMparams = length(namesModelBparams{iModelB_NOMplot});

        for iGroup = 1:nGroups
            iLocPair_all = iLocGroups_all{iGroup};
            fprintf('   - L%s\n', strjoin(string(iLocPair_all), ''))

            switch iLocPair_all(1)
                case 1, nameAsymX = 'EE'; nameAsymY = 'EE';
                case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
                case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
            end

            % X-axis (pA)
            asymX_allSubj = (pA_allSubj(:, iLocPair_all(1))-pA_allSubj(:, iLocPair_all(2)))./(pA_allSubj(:, iLocPair_all(1))+pA_allSubj(:, iLocPair_all(2)));
            asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';
            % switch iLocPair_all(1)
            %     case 1, x_ticks = linspace(0, 20, 5); % EE
            %     case 6, x_ticks = linspace(-5, 15, 5); % HVA
            %     case 5, x_ticks = linspace(0, 16, 5); % VMA (extent is smaller)
            % end
            x_ticks = linspace(-6, 10, 5);

            for iNOMparam = 1:nNOMparams
                NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_NOMplot, iLocPair_all, :, :, iNOMparam));

                asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

                % switch iLocPair_all(1)
                %     case 1 % EE
                %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
                %     case 6 % HVA
                %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
                %     case 5 % VMA
                %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
                % end
                switch iModelB_NOMplot
                    case 1, y_ticks_lb = -[60, 40, 20]; y_ticks_ub = [40, 40, 60];
                    case 3, y_ticks_lb = -[40, 20]; y_ticks_ub = [40, 60];
                end
                y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

                x_ticklabels = nan;
                y_ticklabels = nan;

                nameVarY_figTitle = namesModelBparams{iModelB_NOMplot}{iNOMparam};
                nameVarY_fileTitle = nameVarY_figTitle;

                str_title = sprintf('n=%d, nIter=%d B%d %s (%s) vs. %s (%s)', nsubj, nIterxJob, iModelB_NOMplot, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                basicFxn_drawCorrAsym_permutation(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                xlabel(sprintf('%s of %s (%%)', nameAsymX, nameVarX), 'fontsize', sz_label)
                ylabel(sprintf('%s of %s (%%)', nameAsymX, nameVarY_figTitle), 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks

                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left
                ax.Position = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];

                saveas(gcf, sprintf('%s/n%d_B%d_L%s_NOMp%d.png', nameFolder_Fig_NOM_CorrAsym, nsubj, iModelB_NOMplot, strjoin(string(iLocPair_all), ''), iNOMparam))
                close(gcf)

            end % iNOMparam
        end % iGroup
    end % iModelB_NOMplot
    clear params_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% [NOM] Corr2: CS and NOM params
    clc, fprintf('\n\n 22/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'params_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'NOMparams';

    sz_label = 55;
    sz_labelOffset = .05;
    sz_axOffset = 0.05; % extra breathing room

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

    flag_zeroMean = 0;
    flag_plotIdvdCI = 0;
    flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles

    for iModelB_NOMplot = iModelB_plot_all
        nNOMparams = length(namesModelBparams{iModelB_NOMplot});

        for iSet = 1:numel(iLocSingle_allSets)

            iLocCorr_all = iLocSingle_allSets{iSet};
            fprintf('   - L%s\n', strjoin(string(iLocCorr_all), ''))

            % Obtain CS (x-axis)
            X_allSubj = CS_allSubj(:, iLocCorr_all);
            x_ticks = linspace(1.5, 3.5, 5); % EE
            X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);
            % X_med_allSubj = X_allSubj;

            for iNOMparam = 1:nNOMparams
                NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_NOMplot, iLocCorr_all, :, :, iNOMparam));

                Y_med_allSubj = getCI(NOMp_allIter_allSubj, 1, 3)'; % rotate to match the format needed by basicFxn_drawCorr

                switch iModelB_NOMplot
                    case 1, y_ticks_lb = [0, 0, 0]; y_ticks_ub = [.8, 40, 40];
                    case 3, y_ticks_lb = [0, 0]; y_ticks_ub = [40, 40];
                end

                y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

                x_ticklabels = x_ticks;
                y_ticklabels = y_ticks;

                nameVarY_figTitle = namesModelBparams{iModelB_NOMplot}{iNOMparam};
                nameVarY_fileTitle = nameVarY_figTitle;

                str_title = sprintf('n=%d, nIter=%d B%d %s vs. %s [L%s]', nsubj, nIterxJob, iModelB_NOMplot, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

                basicFxn_drawCorr_permutation(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, nIterxJob);

                xlabel('Contrast sensitivity (1/contrast)', 'fontsize', sz_label)
                ylabel(nameVarY_figTitle, 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left

                ax.Position = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];

                saveas(gcf, sprintf('%s/n%d_B%d_L%s_NOMp%d.png', nameFolder_Fig_NOM_corr, nsubj, iModelB_NOMplot, strjoin(string(iLocCorr_all), ''), iNOMparam))
                close(gcf)
            end % iNOMparam
        end % iSet
    end % iModelB_NOMplot
    clear params_allCond
    fprintf('\n\n Plotting DONE\n\n')


    %% [NOM] CorrAsym2: CS and NOMparams: corr between extents of EE/HVA/VMA
    clc, fprintf('\n\n 23/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'params_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'NOMparams';

    nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

    sz_label = 60;
    sz_labelOffset = 0.05;
    sz_axOffset = 0.05; % extra breathing room

    for iModelB_NOMplot = iModelB_plot_all
        nNOMparams = length(namesModelBparams{iModelB_NOMplot});

        for iGroup = 1:nGroups
            iLocPair_all = iLocGroups_all{iGroup};
            fprintf('   - L%s\n', strjoin(string(iLocPair_all), ''))

            switch iLocPair_all(1)
                case 1, nameAsymX = 'Ecc. effect'; nameAsymY = 'Ecc. effect';
                case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
                case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
            end

            asymX_allSubj = (CS_allSubj(:, iLocPair_all(1))-CS_allSubj(:, iLocPair_all(2)))./(CS_allSubj(:, iLocPair_all(1))+CS_allSubj(:, iLocPair_all(2)));
            switch iLocPair_all(1)
                case 1, x_ticks = linspace(0, 20, 5); % EE
                case 6, x_ticks = linspace(-5, 15, 5); % HVA
                case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
            end
            asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';

            for iNOMparam = 1:nNOMparams
                NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_NOMplot, iLocPair_all, :, :, iNOMparam));

                asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

                switch iModelB_NOMplot
                    case 1, y_ticks_lb = -[30, 20, 20]; y_ticks_ub = [50 40 40];
                    case 3, y_ticks_lb = -[20, 20]; y_ticks_ub = [40 40];
                end
                y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

                x_ticklabels = nan;
                y_ticklabels = nan;

                nameVarY_figTitle = namesModelBparams{iModelB_NOMplot}{iNOMparam};
                nameVarY_fileTitle = nameVarY_figTitle;

                str_title = sprintf('n=%d, nIter=%d B%d %s (%s) vs. %s (%s)', nsubj, nIterxJob, iModelB_NOMplot, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                basicFxn_drawCorrAsym_permutation(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                xlabel(sprintf('\\Delta contrast sensitivity (%%)'), 'FontSize', sz_label);
                ylabel(sprintf('\\Delta %s (%%)', nameVarY_figTitle), 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks

                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left
                ax.Position = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];

                saveas(gcf, sprintf('%s/n%d_B%d_L%s_NOMp%d.png', nameFolder_Fig_NOM_CorrAsym, nsubj, iModelB_NOMplot, strjoin(string(iLocPair_all), ''), iNOMparam))
                close(gcf)

            end % iNOMparam
        end % iGroup
    end % iModelB_NOMplot
    clear params_allCond
    fprintf('\n\n Plotting DONE\n\n')

    %% [NOM] CompAsym2: CS and NOMparams: compare binned estimates
    clc, fprintf('\n\n 24/24 Plotting STARTED......\n\n')

    % Load data
    load(nameFolder_Data_SaveCompile, 'params_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'NOMparams';

    nameFolder_Fig_NOM_CompAsym = sprintf('%s/CompAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_CompAsym)); mkdir(nameFolder_Fig_NOM_CompAsym), end

    sz_label = 15;
    sz_labelOffset = 0.05;
    sz_axOffset = 0.05; % extra breathing room

    for iModelB_NOMplot = iModelB_plot_all
        nNOMparams = length(namesModelBparams{iModelB_NOMplot});

        for iGroup = 1:nGroups
            iLocPair_all = iLocGroups_all{iGroup};
            fprintf('   - L%s\n', strjoin(string(iLocPair_all), ''))

            switch iLocPair_all(1)
                case 1, nameAsymX = 'Ecc. effect'; nameAsymY = 'Ecc. effect';
                case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
                case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
            end

            asymX_allSubj = (CS_allSubj(:, iLocPair_all(1))-CS_allSubj(:, iLocPair_all(2)))./(CS_allSubj(:, iLocPair_all(1))+CS_allSubj(:, iLocPair_all(2)));
            switch iLocPair_all(1)
                case 1, x_ticks = linspace(0, 20, 5); % EE
                case 6, x_ticks = linspace(-5, 15, 5); % HVA
                case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
            end
            asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';

            for iNOMparam = 1:nNOMparams
                NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_NOMplot, iLocPair_all, :, :, iNOMparam));

                asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

                switch iModelB_NOMplot
                    case 1, y_ticks_lb = -[30, 20, 20]; y_ticks_ub = [50 40 60];
                    case 3, y_ticks_lb = -[20, 20]; y_ticks_ub = [40 60];
                end
                y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

                x_ticklabels = nan;
                y_ticklabels = nan;

                nameVarY_figTitle = namesModelBparams{iModelB_NOMplot}{iNOMparam};
                nameVarY_fileTitle = nameVarY_figTitle;

                str_title = sprintf('n=%d, nIter=%d B%d %s (%s) vs. %s (%s)', nsubj, nIterxJob, iModelB_NOMplot, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

                % remove later!!s
                asymX_allIter_allSubj_ = repmat(asymX_allIter_allSubj', 5, 1);
                asymX_allIter_allSubj_ = asymX_allIter_allSubj_+randn(size(asymX_allIter_allSubj_))*mean(asymX_allIter_allSubj_(:))/20;
                asymY_allIter_allSubj_ = repmat(asymY_allIter_allSubj, 5, 1);
                asymY_allIter_allSubj_ = asymY_allIter_allSubj_+randn(size(asymY_allIter_allSubj_))*mean(asymY_allIter_allSubj_(:))/20;

                str_ylabel = sprintf('\\Delta %s (%%)', nameVarY_figTitle);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                basicFxn_compAsym_permutation(asymX_allIter_allSubj_*100, asymY_allIter_allSubj_*100, nBinsCompAsym, y_ticks, sz_fig, str_title, str_ylabel)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % xlabel(sprintf('\\Delta contrast sensitivity (%%)'), 'FontSize', sz_label);
                % ylabel(, 'fontsize', sz_label)

                % Adjust distance between components
                % ax = gca;
                % ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
                %
                % ax.XLabel.Units = 'normalized';
                % ax.YLabel.Units = 'normalized';
                %
                % ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
                % ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left
                % ax.Position = [ ...
                %     ti(1) + sz_axOffset, ...
                %     ti(2) + sz_axOffset, ...
                %     1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                %     1 - ti(2) - ti(4) - 2*sz_axOffset];

                saveas(gcf, sprintf('%s/n%d_B%d_L%s_NOMp%d.png', nameFolder_Fig_NOM_CompAsym, nsubj, iModelB_NOMplot, strjoin(string(iLocPair_all), ''), iNOMparam))
                close(gcf)

            end % iNOMparam
        end % iGroup
    end % iModelB_NOMplot
    clear params_allCond
    fprintf('\n\n Plotting DONE\n\n')

end% iRun
