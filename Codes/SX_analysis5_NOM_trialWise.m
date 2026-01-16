% SX_analysis5_NOM_Trialwise.m
% Author: Shutian Xue
% Purpose: This script compiles, analyzes and plots boostrapped data for each observer saved Data_NOM_Trialwise

clear all, clc, close all

addpath(genpath('PF_RC/Codes/')) % manually add to save time
addpath(genpath('PF_RC/Codes/fxn_analysis_RC_v2')) % manually add to save time

% Setting parameters for the analysis
%----------------
SX_RC1_setting
%----------------

%%%%%% Should copy from OOD_sim %%%%%%%%%%
% noiseCST_all = [0, .1, .2, .5]; % Noise contrast sensitivity thresholds
% gaborCST_all = [.1, .5]; % Gabor contrast sensitivity thresholds
% nTrials_all = [5000]; % Number of trials per condition
% noiseP_all = [0, 0.1, 0.2]; % Proportion of noise trials
iModelA_sim_all = [1:2]; % 1=RC-derived template, 2=use IO template, 2=randomize template,
iModelB_sim_all = [1:4]; % see SX_RC1_setting for namesModelB
nBoot = 100; % Number of bootstraps

iLocComb_all = [1,8, 6,7, 5,3]; % combination of locations
iLocGroups_all = {[1,8], [6,7], [5,3]} ; nGroups = length(iLocGroups_all);
iLocSingle_all = 1:5; nLocSingle = length(iLocSingle_all);
iLocSingle_allSets = {[1,2,4,5,3], [2,4,5,3], [1,6,5,3], [6,5,3]};

iFamily_ORI = 1; %1=scaled Gaussian
iFamily_SF = 2; % 2=log parabola
iFamily_perF = [iFamily_ORI, iFamily_SF];
namesMetrics_prob = {'pYES', 'pA'}; nMetrics_prob = length(namesMetrics_prob); namesMetrics_prob_full = {sprintf('Predicted detection prob.\nMeasured detection rate'), sprintf('Predicted consistency prob.\nMeasured resp. consistency')};
iDataset_plotRC = 1; %1=Template set; 2=Full set

iModelA_plot = 1; % just plot the RC-derived
iModelB_plot = 1; % 1=full model; 2=No Rho; 3=No induced noise
nNOMparams = length(namesModelBparams{iModelB_plot});
iDataset_plotNOM = 1; % 1=test set, 2=full; search for metrics_allCond

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NOM_mode = 2; % 1= aggregate model, 2 = trial-wise model
flag_subjIsHuman = 1; % 1=human subject, 0=IO data

iModelA_all = iModelA_sim_all;
iModelB_all = iModelB_sim_all;

% Define subject list
if flag_subjIsHuman
    % n=15, everyone
    subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC', 'SR'}; nRows_subj = 3; nCols_subj = 5;
    nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195];

    % n=12: no AS, CS, RC
    % subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'HL', 'FH', 'HA', 'DT', 'DU', 'SR'}; nRows_subj = 3; nCols_subj = 4;
    % nblocks_allSubj = [200, 240, 210, 240, 220, 220, 220, 205, 210, 205, 205, 195];

    % n=13: no AS, CS
    % subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'HL', 'FH', 'HA',  'DT', 'DU', 'RC', 'SR'}; nRows_subj = 3; nCols_subj = 5;
    % nblocks_allSubj = [200, 240, 210, 240, 220, 220, 220, 205, 210, 205, 205, 205, 195];

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

nameFolder_Data_SaveCompile = sprintf('%s/n%d_n%d', nameFolder_Data_NOM_Trialwise, nsubj, nBoot);

% Names of the performance metrics being analyzed
nModelsA = length(iModelA_all);
nModelsB = length(iModelB_all);
nLocComb8 = 8;
iIC_plot = 3; % % Which information criterion to plot (1 = AIC, 2 = AICc, 3 = BIC)

% Function handles for calculating information criteria (IC)
getAIC_nLL = @(nLL, nParams) 2*nParams+2*nLL;
getAICc_nLL = @(nLL, nParams, nData) getAIC_nLL(nLL, nParams) + (2*nParams*(nParams+1))/(nData-nParams-1);
getBIC_nLL = @(nLL, nParams, nData) nParams*log(nData)+2*nLL;

% Define folder to save figures
nameFolder_Fig_NOM_Trialwise = sprintf('%s/NOM_Trialwise_%d%d', nameFolder_Figures, nORI, nSF);
if isempty(dir(nameFolder_Fig_NOM_Trialwise)), mkdir(nameFolder_Fig_NOM_Trialwise), end

% Print a header to summarize the setting
fprintf('NOM trial-wise analysis settings:\n')
fprintf(' - nsubj = %d\n', nsubj)
fprintf(' - nBoot = %d\n', nBoot)
fprintf(' - Models A: %s\n', strjoin(string(iModelA_all), ', '));
fprintf(' - Models B: %s\n', strjoin(string(iModelB_all), ', '));
fprintf(' - Combined locations:  %s\n', strjoin(string(iLocComb_all), ', '));
fprintf(' - ORI tuning function:  %s\n', namesFamily_all{iFamily_ORI});
fprintf(' - SF tuning function:  %s\n', namesFamily_all{iFamily_SF});
fprintf(' - Number of Bins: %d\n', nBins);
fprintf(' - Information Criterion to plot: %s\n\n\n', namesIC{iIC_plot});

%% (Time-consuming!) Compile data from boostrapped data

% Preallocate arrays for storing results across all conditions
metrics_allCond = nan(nModelsA, nLocComb8, nsubj, nBoot, nDatasets, nMetrics); % nDatasets: 1=Full; 2=Test; 11=number of metrics saved in OOD_xx_compIV.m; see "fxn_getMetrics"
template_tmpl_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nORI, nSF);
template_full_allCond = template_tmpl_allCond;
IV_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nBins);
nTrials_allCond = IV_allCond;
metric_data_allCond = nan(nModelsA, nModelsB, nLocComb8, nMetrics_prob, nsubj, nBoot, nBins); % "nMetrics" is predefined in this script
metric_pred_allCond = metric_data_allCond;
nLL_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot);
params_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, 3); % 3=Pre-allocate for max params
IC_nLL_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, 3); % 3=AIC, AICc, BIC
R2_NOM_allCond= nan(nModelsA, nModelsB, nLocComb8, nMetrics_prob, nsubj, nBoot);

% Fitting tuning curves
% nDatasets: 1=Full; 2=Template or test
sep_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nDatasets);
margORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nDatasets, nORI);
margPred_ORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nDatasets, nORI);
margParams_ORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nDatasets, length(namesParams_all{iFamily_ORI}));
margR2_ORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nDatasets);
margTunC_ORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nDatasets, length(namesTunC_unit_perF{iFamily_ORI, 2}));

margSF_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nDatasets, nSF);
margPred_SF_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nDatasets, nSF);
margParams_SF_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nDatasets, length(namesParams_all{iFamily_SF}));
margR2_SF_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nDatasets);
margTunC_SF_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nDatasets, length(namesTunC_unit_perF{iFamily_SF, 2}));

% The main analysis loop for all conditions for each subject (Models x Loc x Metrics)
for iModelA = iModelA_all %1:nModelsA
    fprintf('\n\n === Model A%d === ', iModelA)

    for iModelB = iModelB_all %1:nModelsB
        fprintf('\n     Model B%d: ', iModelB)
        nParamsB = length(namesModelBparams{iModelB});

        for iLocSingle = iLocSingle_all
            % for iLocComb = iLocComb_all
            fprintf(' L%d ', iLocSingle)

            for isubj = 1:nsubj

                subjName = subjList{isubj};

                if flag_subjIsHuman
                    % nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, nameIO);
                    nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocSingle); % MUST be the same as OOD_NOM_Trialwise_compIV.m, Line 61
                else
                    % nameFolder_NOM_save = sprintf('%s/ORI%dSF%d/%s/L%d', nameFolder_NOM0, nORI, nSF, subjName, iLocComb);
                    nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, subjName); % MUST be the same as OOD_NOM_Trialwise_compIV.m, Line 68
                end

                % Load IVs and derived templates
                nameFile_compIV = sprintf('%s/n%d_A%d_compIV.mat', nameFolder_NOM_save, nBoot, iModelA);
                % if isempty(dir(nameFile_compIV)), error('ALERT: File %s does not exist!', nameFile_compIV); end
                if isempty(dir(nameFile_compIV))
                    fprintf(' "%sL%dA%d" ', subjName, iLocSingle, iModelA)
                else
                    load(nameFile_compIV, '*_allBoot')
                end

                % Load predictions
                nameFile_fitNOM = sprintf('%s/n%d_A%dB%d.mat', nameFolder_NOM_save, nBoot, iModelA, iModelB);
                % if isempty(dir(nameFile_fitNOM)), error('ALERT: File %s does not exist!', nameFile_fitNOM); end
                if isempty(dir(nameFile_fitNOM))
                    fprintf(' "%sL%dA%dB%d" ', subjName, iLocSingle, iModelA, iModelB)
                else
                    load(nameFile_fitNOM, '*_allBoot')
                end

                % Loaded *measured* metrics of the full set and test set
                % "metrics" is defined in OOD_xx_compIV: metrics = [metrics_full; metrics_tmpl; metrics_tmpl; metrics_test];
                metrics_allCond(iModelA, iLocSingle, isubj, :, :, :) = data_metrics_allBoot(:, [4,1], :); % see OOD_xx_compIV.m search for "data_metrics_allBoot"

                % Pre-allocate temporary arrays for this subject and model
                IV_allBoot = nan(nBoot, nBins); % IV for each bin
                nTrials_allBoot = IV_allBoot; % Trial count for each bin
                metric_data_allBoot = nan(nMetrics_prob, nBoot, nBins); % Measured data for each bin
                metric_pred_allBoot = metric_data_allBoot; % Predicted data for each bin
                margTuningC_ORI_allBoot = nan(nBoot, nDatasets, length(namesTunC_unit_perF{iFamily_ORI, 2}));
                margTuningC_SF_allBoot = nan(nBoot, nDatasets, length(namesTunC_unit_perF{iFamily_SF, 2}));

                for iBoot = 1:nBoot
                    % Extract IV per bin
                    IV_allBoot(iBoot, :) = pred_metrics_allBoot{iBoot}.metrics.IV_allBins;

                    % Extract the number of trials per bin
                    nTrials_allBoot(iBoot, :) = pred_metrics_allBoot{iBoot}.metrics.nTrials_allBins;

                    % Extract metrics (data and predictions) per bin
                    pred_metrics = pred_metrics_allBoot{iBoot}.metrics; % Extract once for faster access
                    for iMetric_prob = 1:nMetrics_prob
                        switch namesMetrics_prob{iMetric_prob}
                            case 'pYES'
                                metric_data_allBoot(iMetric_prob, iBoot, :) = pred_metrics.pYES_data_allBins;
                                metric_pred_allBoot(iMetric_prob, iBoot, :) = pred_metrics.pYES_pred_allBins;
                            case 'pA'
                                metric_data_allBoot(iMetric_prob, iBoot, :) = pred_metrics.pA_data_allBins;
                                metric_pred_allBoot(iMetric_prob, iBoot, :) = pred_metrics.pA_pred_allBins;
                            case 'pC'
                                metric_data_allBoot(iMetric_prob, iBoot, :) = pred_metrics.pC_data_allBins;
                                metric_pred_allBoot(iMetric_prob, iBoot, :) = pred_metrics.pC_pred_allBins;
                        end
                    end % iMetric
                    % Calculate information critertion based on nLL
                    nData = sum(nTrials_allBoot(iBoot, :));

                    % nLL_allBoot
                    % nLL = nLL_test_allBoot(iBoot);
                    % nLL_allCond(iModelA, iModelB, iLocSingle, isubj, iBoot) = IC_nLL;

                    % Derive ICs based on nLL
                    % IC_nLL = [getAIC_nLL(nLL, nParamsB), getAICc_nLL(nLL, nParamsB, nData), getBIC_nLL(nLL, nParamsB, nData)];
                    % IC_nLL_allCond(iModelA, iModelB, iLocSingle, isubj, iBoot, :) = IC_nLL;

                    % Compute tuning characteristics based on the fitted parameters
                    for iDataset=1:2 % needs to matchOOD_xx_compIV ("for iDataset = 1:2")
                        % ORI
                        iFeature= 1;
                        %======================%
                        margTuningC_ORI_allBoot(iBoot, iDataset, :) = fxn_getTuningC(axis_tuning{iFeature}, iFeature, iFamily_ORI, squeeze(margPred_ORI_allBoot(iBoot, iDataset, :)), squeeze(margParams_ORI_allBoot(iBoot, iDataset, :)));
                        %======================%
                        
                        % SF
                        iFeature=2;
                        %======================%
                        margTuningC_SF_allBoot(iBoot, iDataset, :) = fxn_getTuningC(axis_tuning{iFeature}, iFeature, iFamily_SF, squeeze(margPred_SF_allBoot(iBoot, iDataset, :)), squeeze(margParams_SF_allBoot(iBoot, iDataset, :)));
                        %======================%
                    end
                end % end of iBoot

                % Calculate R2 for NOM prediction (for pYES and pA)
                R2_NOM_allBoot = nan(nMetrics_prob, nBoot);
                for iMetric_prob = 1:nMetrics_prob
                    for iBoot = 1:nBoot
                        nData = squeeze(nTrials_allBoot(iBoot, :));
                        metric_data  = squeeze(metric_data_allBoot(iMetric_prob, iBoot, :));  % ground truth
                        metric_pred = squeeze(metric_pred_allBoot(iMetric_prob, iBoot, :)); % prediction
                        % Remove NaNs if any
                        valid = ~(isnan(metric_data) | isnan(metric_pred));
                        metric_data = metric_data(valid);
                        metric_pred = metric_pred(valid);

                        % With weighting
                        % metric_ave = sum(nData .* metric_data) / sum(nData);
                        % SSres = sum(nData .* (metric_data - metric_pred).^2);
                        % SStot = sum(nData .* (metric_data - metric_ave).^2);

                        % No weighting
                        SSres = sum((metric_data - metric_pred).^2);
                        SStot = sum((metric_data - mean(metric_data)).^2);

                        if SStot == 0, R2 = NaN; else, R2 = 1 - SSres/SStot; end
                        R2_NOM_allBoot(iMetric_prob, iBoot) = R2;
                    end
                end

                % Store results
                % Templates
                [template_tmpl_med, ~, ~, template_tmpl_sem] = getCI(template_tmpl_allBoot, 1, 1);
                [template_full_med, ~, ~, template_full_sem] = getCI(template_full_allBoot, 1, 1);
                template_tmpl_allCond(iModelA, iModelB, iLocSingle, isubj, :, :) = template_tmpl_med;
                template_full_allCond(iModelA, iModelB, iLocSingle, isubj, :, :) = template_full_med;

                % Tuning functions: marg, predictions, estimated parameters and separability
                sep_allCond(iModelA, iModelB, iLocSingle, isubj, :, :) = sep_allBoot;

                margORI_allCond(iModelA, iModelB, iLocSingle, isubj, :, :, :) = margORI_allBoot; % directly loaded
                margR2_ORI_allCond(iModelA, iModelB, iLocSingle, isubj, :, :) = margR2_ORI_allBoot; % directly loaded
                margPred_ORI_allCond(iModelA, iModelB, iLocSingle, isubj, :, :, :) = margPred_ORI_allBoot; % directly loaded
                margParams_ORI_allCond(iModelA, iModelB, iLocSingle, isubj, :, :, :) = margParams_ORI_allBoot; % directly loaded
                margTunC_ORI_allCond(iModelA, iModelB, iLocSingle, isubj, :, :, :) = margTuningC_ORI_allBoot;

                margSF_allCond(iModelA, iModelB, iLocSingle, isubj, :, :, :) = margSF_allBoot; % directly loaded
                margR2_SF_allCond(iModelA, iModelB, iLocSingle, isubj, :, :) = margR2_SF_allBoot; % directly loaded
                margPred_SF_allCond(iModelA, iModelB, iLocSingle, isubj, :, :, :) = margPred_SF_allBoot; % directly loaded
                margParams_SF_allCond(iModelA, iModelB, iLocSingle, isubj, :, :, :) = margParams_SF_allBoot; % directly loaded
                margTunC_SF_allCond(iModelA, iModelB, iLocSingle, isubj, :, :, :) = margTuningC_SF_allBoot;

                % NOM: IVs, predictions, nLL and parameters
                IV_allCond(iModelA, iModelB, iLocSingle, isubj, :, :) = IV_allBoot;
                nTrials_allCond(iModelA, iModelB, iLocSingle, isubj, :, :) = nTrials_allBoot;
                metric_data_allCond(iModelA, iModelB, iLocSingle, :, isubj, :, :) = metric_data_allBoot;
                metric_pred_allCond(iModelA, iModelB, iLocSingle, :, isubj, :, :) = metric_pred_allBoot;
                nLL_allCond(iModelA, iModelB, iLocSingle, isubj, :) = nLL_test_allBoot;
                % IC_nLL_allCond is compiled above
                params_allCond(iModelA, iModelB, iLocSingle, isubj, :, 1:nParamsB) = params_est_allBoot;
                R2_NOM_allCond(iModelA, iModelB, iLocSingle, :, isubj, :) = R2_NOM_allBoot;

                clear *allBoot

            end % isubj
        end % iLocComb
    end % iModelB
end % iModelA

fprintf('\n\n *** NOM outputs compiled ***\n\n')

%%% Build combined locations by averaging single-location outputs
meanOverLoc = @(X, locDim, locOld) mean(X, locDim, 'omitnan');

for iMap = 6:8
    switch iMap
        case 6, locOld = [2 4];      % 6 (HM) = 2 (left) and 4 (right)
        case 7, locOld = [3 5];      % 7 (VM) = 3 (upper) and 5 (lower)
        case 8, locOld = 2:5;        % 8 (Perifovea) = 2 to 5
    end

    % locDim = 2
    metrics_allCond(:,iMap,:,:,:,:) = meanOverLoc(metrics_allCond(:,locOld,:,:,:,:), 2, locOld);

    % locDim = 3
    template_tmpl_allCond(:,:,iMap,:,:,:)    = meanOverLoc(template_tmpl_allCond(:,:,locOld,:,:,:), 3, locOld);
    template_full_allCond(:,:,iMap,:,:,:)    = meanOverLoc(template_full_allCond(:,:,locOld,:,:,:), 3, locOld);

    IV_allCond(:,:,iMap,:,:,:)               = meanOverLoc(IV_allCond(:,:,locOld,:,:,:), 3, locOld);
    nTrials_allCond(:,:,iMap,:,:,:)          = meanOverLoc(nTrials_allCond(:,:,locOld,:,:,:), 3, locOld);

    metric_data_allCond(:,:,iMap,:,:,:,:)    = meanOverLoc(metric_data_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    metric_pred_allCond(:,:,iMap,:,:,:,:)    = meanOverLoc(metric_pred_allCond(:,:,locOld,:,:,:,:), 3, locOld);

    nLL_allCond(:,:,iMap,:,:)                = meanOverLoc(nLL_allCond(:,:,locOld,:,:), 3, locOld);
    params_allCond(:,:,iMap,:,:,:)           = meanOverLoc(params_allCond(:,:,locOld,:,:,:), 3, locOld);
    IC_nLL_allCond(:,:,iMap,:,:,:)           = meanOverLoc(IC_nLL_allCond(:,:,locOld,:,:,:), 3, locOld);
    R2_NOM_allCond(:,:,iMap,:,:,:)           = meanOverLoc(R2_NOM_allCond(:,:,locOld,:,:,:), 3, locOld);

    sep_allCond(:,:,iMap,:,:,:)              = meanOverLoc(sep_allCond(:,:,locOld,:,:,:), 3, locOld);

    margORI_allCond(:,:,iMap,:,:,:,:)        = meanOverLoc(margORI_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    margPred_ORI_allCond(:,:,iMap,:,:,:,:)   = meanOverLoc(margPred_ORI_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    margParams_ORI_allCond(:,:,iMap,:,:,:,:) = meanOverLoc(margParams_ORI_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    margR2_ORI_allCond(:,:,iMap,:,:,:)       = meanOverLoc(margR2_ORI_allCond(:,:,locOld,:,:,:), 3, locOld);
    margTunC_ORI_allCond(:,:,iMap,:,:,:,:)   = meanOverLoc(margTunC_ORI_allCond(:,:,locOld,:,:,:,:), 3, locOld);

    margSF_allCond(:,:,iMap,:,:,:,:)         = meanOverLoc(margSF_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    margPred_SF_allCond(:,:,iMap,:,:,:,:)    = meanOverLoc(margPred_SF_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    margParams_SF_allCond(:,:,iMap,:,:,:,:)  = meanOverLoc(margParams_SF_allCond(:,:,locOld,:,:,:,:), 3, locOld);
    margR2_SF_allCond(:,:,iMap,:,:,:)        = meanOverLoc(margR2_SF_allCond(:,:,locOld,:,:,:), 3, locOld);
    margTunC_SF_allCond(:,:,iMap,:,:,:,:)    = meanOverLoc(margTunC_SF_allCond(:,:,locOld,:,:,:,:), 3, locOld);
end

fprintf('\n\n *** Combined locations (6,7,8) created by averaging single locations *** \n\n');

%%% Save the organized data for all subjects
save(nameFolder_Data_SaveCompile, '*_allCond')
clear *_allCond
fprintf('\n\n *** *_allCond saved and cleared *** \n\n');

%% Compile behavioral data (revise this part later, as the behav should also come from bootstrrapped data (metric_data_allCond), not from the raw data)
CS_allSubj = nan(nsubj, nLocComb8); % 8 is max number of combined locs,see namesLocComb
dprime_allSubj = CS_allSubj;
criterion_allSubj = CS_allSubj;
RT_allSubj = CS_allSubj;
pC_allSubj = CS_allSubj;
pA_allSubj = CS_allSubj;

for isubj = 1:nsubj
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    nameFile_behavMeas = sprintf('%s/%s%d/%s_behavMeas.mat', nameFolder_Data_OOD, subjName, nblocks, subjName);

    load(nameFile_behavMeas, '*_perSess_perLoc')

    for iLoc = 1:nLocComb8
        switch iLoc
            case 6, iLoc_allLoc = [2,4];
            case 7, iLoc_allLoc = [5,3];
            case 8, iLoc_allLoc = 2:5;
            otherwise
                iLoc_allLoc = iLoc;
        end
        CS_allSubj(isubj, iLoc) = mean(1./cst_perSess_perLoc(:, iLoc_allLoc), 'all');
        dprime_allSubj(isubj, iLoc) = mean(dprime_perSess_perLoc(:, iLoc_allLoc), 'all');
        criterion_allSubj(isubj, iLoc) = mean(criterion_perSess_perLoc(:, iLoc_allLoc), 'all');
        RT_allSubj(isubj, iLoc) = median(RT_perSess_perLoc(:, iLoc_allLoc), 'all');
        pC_allSubj(isubj, iLoc) = mean(pC3_perSess_perLoc(:, iLoc_allLoc, 1), 'all');
        % pYES_allSubj(isubj, iLocComb) = mean(pYES_perSess_perLoc(:, iLoc_allLoc));
        pA_allSubj(isubj, iLoc) = mean(pA3_perSess_perLoc(:, iLoc_allLoc, 1), 'all');
    end % iLoc

    clear *_perSess_perLoc

end % isubj

% Save the organized data for all subjects
save(nameFolder_Data_SaveCompile, '*_allSubj', '-append')

fprintf('\nBehav data compiled\n')

%% Plot templates for IDVD and group averages
clc, fprintf('\n\n Plotting STARTED\n\n')
% Load data
load(nameFolder_Data_SaveCompile, 'template_*_allCond')

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
        
        if iLocComb==1, cLim = [-.03, .26];
        else, cLim = [-.03, .13];
        end
        RCplot_2Dkernel(e2D_ave, cLim)
        title(sprintf('n=%d L%d %s (%s)', nsubj, iLocComb, namesLocComb{iLocComb}, str_dataset))
        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_group_L%d_A%dB%d_%s.png', nameFolder_Fig_NOM_Template, nsubj, iLocComb, iModelA_plot, iModelB_plot, str_dataset))
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
    %     drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d_A%dB%d_%s.png', nameFolder_Fig_NOM_Template, nsubj, iLocComb_all(iiLoc), iModelA_plot, iModelB_plot, str_title))
    %     close(gcf)
    %
    % end % iiLoc
end % iDataset
clear template_*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Plot separability
clc, fprintf('\n\n Plotting STARTED\n\n')
str_dataset = namesDataset{iDataset_plotRC};

% Load data
load(nameFolder_Data_SaveCompile, 'sep_allCond')

% Define folder for saving figures
nameFolder_Fig_Sep = sprintf('%s/Separability', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_Sep)), mkdir(nameFolder_Fig_Sep), end

wd_border = 2;
sz_title = 10;
sz_label = 25; %40
sz_ticks = 35;
wd_all = 3;
sz_marker = 20;

for iSet = 1:numel(iLocSingle_allSets)

    iLocSingle_all = iLocSingle_allSets{iSet};

    [sep_med, ~, ~, sep_neg, sep_pos] = getCI(sep_allCond(iModelA_plot, iModelB_plot, iLocSingle_all, :, :, iDataset_plotRC), 1, 5);

    figure('Position', [0 0 1e3 600])
    hold on, box on
    yline(.5, '--', 'linewidth', wd_all);

    str_sepValues = [];
    for iiLoc = 1:length(iLocSingle_all)
        x = (1:nsubj) + .1 * iiLoc;
        for isubj = 1:nsubj
            errorbar(x(isubj), sep_med(iiLoc, isubj), sep_neg(iiLoc, isubj), sep_pos(iiLoc, isubj),'color', colors_comb(iLocSingle_all(iiLoc), :), 'CapSize', 0)
            plot(x(isubj), sep_med(iiLoc, isubj), markers_allSubj{isubj}, ...
                'MarkerFaceColor', colors_comb(iLocSingle_all(iiLoc), :), 'MarkerEdgeColor', colors_comb(iLocSingle_all(iiLoc), :), ...
                'MarkerSize', sz_marker, 'linewidth', 3)
        end
        str_sepValues = [str_sepValues, sprintf('L%d: %.2f (%.2f); ', iLocSingle_all(iiLoc), mean(sep_med(iiLoc, :)), std(sep_med(iiLoc, :))/sqrt(nsubj))];
    end

    ylabel('Correlation', 'FontSize', sz_label)
    yticks([0, .5:.1:1])
    ylim([.5, 1])

    xticks(1:nsubj)
    xticklabels(1:nsubj) % 
    % xticklabels(subjList) % mute later
    xlim([0,nsubj+1])
    xlabel('Observer #')

    ax = gca; ax.XAxis.FontSize = sz_ticks; ax.YAxis.FontSize = sz_ticks; ax.LineWidth = wd_border;

    %%% ANOVA
    text_ANOVA = '...ANOVA...';
    indLoc = repmat(1:length(iLocSingle_all), nsubj, 1)';
    text_ANOVA = print_nANOVA({'Loc'}, sep_med(:), {indLoc(:)}, nsubj);
    
    %=== CI of ANOVA ===
    p_allBoot = nan(nBoot, 1);
    for iBoot = 1:nBoot
        x = squeeze(sep_allCond(iModelA_plot, iModelB_plot, iLocSingle_all, :, iBoot, iDataset_plotRC));
        [~, tbl] = print_nANOVA({'Loc'}, x(:), {indLoc(:)}, nsubj);
        p_allBoot(iBoot) = tbl{2,7};
    end
    [p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 1);
    text_ANOVA = [text_ANOVA, sprintf('p=%.3f [%.3f, %.3f]\n', p_med, p_lb, p_ub)];
    %===============

    %%% ttest
    if length(iLocSingle_all)==2
        [h,p,ci,stats] = ttest(sep_med(:, 1), sep_med(:, 2));
        cohenD = fxn_getES(sep_med(:, 1), sep_med(:, 2));
    end

    set(findall(gcf, '-property', 'linewidth'), 'linewidth', 2)
    set(findall(gcf, '-property', 'fontsize'), 'fontsize', 30)

    title(sprintf('n%d [A%dB%d] L%s %s\n%s\n%s', nsubj, iModelA_plot, iModelB_plot, strjoin(string(iLocSingle_all), ''), str_dataset, str_sepValues, text_ANOVA), 'FontSize', sz_title)
    drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%s_A%dB%d_%s.png', nameFolder_Fig_Sep, nsubj, strjoin(string(iLocSingle_all), ''), iModelA_plot, iModelB_plot, str_dataset))
    close(gcf)
end

clear sep_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Tuning functions for group averages
clc, fprintf('\n\n Plotting STARTED\n\n')
% Load data
load(nameFolder_Data_SaveCompile, 'marg*_allCond', 'margPred*_allCond', 'margR2*_allCond')

% Define folder for saving figures
nameFolder_Fig_NOM_Tuning = sprintf('%s/Tuning', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_NOM_Tuning)), mkdir(nameFolder_Fig_NOM_Tuning), end

wd_border = 4; % default 5
sz_ticks = 35;% default 35

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

            figure('Position', [0 0 1.1e3 8e2]) % default 8e2
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
            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_group_L%d%d_A%dB%d_%s_%s.png', nameFolder_Fig_NOM_Tuning, nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, namesDataset{iDataset}))
            close(gcf)
        end % iFeature
    end % iGroup
end % iDataset = 1:2
clear marg*_allCond margPred*_allCond margR2*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Tuning functions for idvd (LARGE FIGURES!!)
clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'marg*_allCond', 'margPred*_allCond', 'margParams*_allCond')

% Define folder for saving figures
nameFolder_Fig_NOM_Tuning = sprintf('%s/Tuning', nameFolder_Fig_NOM_Trialwise);
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

            figure('Position', [0 0 2e3 2e3])

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
                    patch([xaxis, flip(xaxis)], [margPred_lb', flip(margPred_ub')], color_comb, 'FaceAlpha', .3, 'linestyle', 'none')
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
            sgtitle(sprintf('n=%d L%d%d [A%dB%d] %s tuning (%s)', nsubj, iLocPair_all,iModelA_plot, iModelB_plot, namesFeature{iFeature}, namesDataset{iDataset}))
            set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10)

            % Save the figure
            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_A%dB%d_%s_%s.png', nameFolder_Fig_NOM_Tuning, nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, namesDataset{iDataset}))
            close(gcf)
        end % i=1:2
    end % iFeature
end % iGroup
clear marg*_allCond margPred*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Tuning characteristics
clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'margTunC_*_allCond')

% Define folder for saving figures
nameFolder_Fig_tunC = sprintf('%s/TuningCs', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_tunC)), mkdir(nameFolder_Fig_tunC), end

flag_plotDist = 0;
flag_plotIDVD = 1;
flag_plotDiff = 1;
paramMode = 2;
sz_fig = [350 350]; % size of the figure canvas

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
                            y_ticks_all{2} = linspace(0, 100, 5); % ORI band
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
                                y_ticks_all{2} = linspace(.01, .09, 5); % SF peak amp
                            end
                            y_ticks_all{3} = linspace(0, 2, 5); % SF bandwidth
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
            % margTunC_ORI_allCond: nModelA x nModelB x nLoc_pair x nsubj x nBoot x nDataset x nTunC
            % tunC_allSubj_allBoot: nLoc_pair x nsubj x nBoot
            switch iFeature
                case 1
                    data_allBoot_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                    data_obs_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, 1, 2, iTunC)); % full data, without resampling
                case 2
                    data_allBoot_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
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
                            data_allBoot_allSubj = log2(data_allBoot_allSubj); % pref SF
                        case 4, ref = 0; % baseline
                    end
            end

            str_title = sprintf('n=%d L%d%d [A%dB%d] [%s] | %s %s', nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesDataset{iDataset_plotRC}, namesFeature{iFeature}, namesTunC{iTunC});

            switch flag_plotDist
                case 0
                    %------------------------------%
                    basicFxn_drawBars_permutation(data_allBoot_allSubj, ref, colors, x_ticks, y_ticks_all{iTunC}, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nBoot, markers_allSubj)
                    % ------------------------------%
                    ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily, 2}{iTunC}))
                case 1
                    %------------------------------%
                    basicFxn_drawDist_permutation(data_allBoot_allSubj, ref, colors, x_ticks, y_ticks_all{iTunC}, y_ticklabels, str_title, sz_fig, nBoot, nsubj)
                    %------------------------------%
                    xlabel(sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily, 2}{iTunC}))
                    ylabel('Probabillity')
            end
            % Save the figure
            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_A%dB%d_%s%d.png', nameFolder_Fig_tunC, nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, iTunC))
            close(gcf)

        end % end of iTunC
    end % end of iFeature
end % end of iGroup
clear margTunC_*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Behavioral metrics: single locations
fprintf('\n\n Plotting STARTED\n\n')
% Load data
load(nameFolder_Data_SaveCompile, 'metrics_allCond')

% Define folder for saving figures
nameFolder_Fig_behav = sprintf('%s/Behav', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_behav)), mkdir(nameFolder_Fig_behav), end

flag_plotIDVD = 1;
namesMetrics_behav = {'CS', 'pA', 'dprime', 'criterion', 'pC', 'RT'}; nMetrics_behav = length(namesMetrics_behav);
[d70,~] = SX_sim06_SDT(.7, .3);

sz_wd_perBar = 70;

for iSet = 1:numel(iLocSingle_allSets)

    iLocSingle_all = iLocSingle_allSets{iSet};

    nBars = numel(iLocSingle_all);
    sz_fig = [nBars*sz_wd_perBar, 350];

    for iMetric_prob=1:nMetrics_behav
        switch iMetric_prob
            case 1, iMetric_vec = 10; x_ticks = linspace(1.6, 3.6, 5); flag_plotIDVD = 1; ref=nan;
            case 2, iMetric_vec = 6; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 1; ref=nan;
            case 3, iMetric_vec = 1; x_ticks = linspace(0, 1.6, 5); flag_plotIDVD = 0; ref=d70;
            case 4, iMetric_vec = 2; x_ticks = linspace(-1,1, 5); flag_plotIDVD = 0; ref=0;
            case 5, iMetric_vec = 3; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 0; ref=nan;
            case 6, iMetric_vec = 11; x_ticks = linspace(0, .2, 5); flag_plotIDVD = 0; ref=nan;
        end
        x_ticks = round(x_ticks, 2);

        data_allBoot_allSubj = squeeze(metrics_allCond(iModelA_plot, iLocSingle_all, :, :, iDataset_plotNOM, iMetric_vec));

        str_title = sprintf('%s L%s', namesMetrics_behav{iMetric_prob}, strjoin(string(iLocSingle_all), ''));
        %------------------------------%
        basicFxn_drawBars_permutation(data_allBoot_allSubj, ref, colors_comb(iLocSingle_all, :), namesLocComb(iLocSingle_all), x_ticks, x_ticks, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nBoot, markers_allSubj)
        % ------------------------------%
        % ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily, 2}{iTunC}))
        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%s_%s.png', nameFolder_Fig_behav, nsubj, strjoin(string(iLocSingle_all), ''), namesMetrics_behav{iMetric_prob}))
        close(gcf)
    end % iMetric
end % iSet
clear metrics_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Behavioral metrics: paired locations
fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'metrics_allCond')

sz_wd_perBar = 100;
nBars = 2;
sz_fig = [nBars*sz_wd_perBar, 250];

for iGroup = 1:nGroups
    iLocPair_all = iLocGroups_all{iGroup};
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

        data_allBoot_allSubj = squeeze(metrics_allCond(iModelA_plot, iLocPair_all, :, :, iDataset_plotNOM, iMetric_vec));

        str_title = sprintf('%s L%d%d', namesMetrics_behav{iMetric_prob}, iLocPair_all);
        %------------------------------%
        basicFxn_drawBars_permutation(data_allBoot_allSubj, ref, colors_comb(iLocPair_all, :), namesLocComb(iLocPair_all), x_ticks, x_ticks, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nBoot, markers_allSubj)
        %------------------------------%
        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_%s.png', nameFolder_Fig_behav, nsubj, iLocPair_all, namesMetrics_behav{iMetric_prob}))
        close(gcf)
    end % iMetric
end % iGroup
clear metrics_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Corr0: Corr between CS and tun params (just to check bound-hitting)
clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'margParams*_allCond')

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
type_corr = 'pearson';
type_tail = 'both';

for iSet = 1

    iLocCorr_all = iLocSingle_allSets{iSet};

    % Obtain CS (x-axis)
    X_allSubj = CS_allSubj(:, iLocCorr_all);
    x_ticks = linspace(1.5, 3.5, 5); % EE
    X_allBoot_allSubj = repmat(X_allSubj, 1, 1, nBoot);

    for iFeature = 1:nFeatures
        iFamily = iFamily_perF(iFeature);
        namesTunCs = namesTunC_unit_perF{iFamily, 1};
        nTunCs_full = length(namesTunCs);
        
        for iTunC = 1:nTunCs_full
            switch iFeature
                case 1, NOMp_allBoot_allSubj = squeeze(margParams_ORI_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
                case 2, NOMp_allBoot_allSubj = squeeze(margParams_SF_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
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

            str_title = sprintf('n=%d, %s vs. %s [L%s]', nsubj, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

            basicFxn_drawCorr_permutation(X_allBoot_allSubj, NOMp_allBoot_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, type_corr, type_tail, str_title, markers_allSubj, nBoot);
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

            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_%s_L%s.png', nameFolder_Fig_NOM_corr, nsubj, nameVarY_fileTitle, strjoin(string(iLocCorr_all), '')))
            close(gcf)

        end % iTunC
    end % iFeature
end % iSet
clear margParams*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Corr1: Corr between CS and tunC
clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'margTunC*_allCond')

nameVarX = 'CS';
nameVarY = 'tunC';
sz_label = 55;
sz_labelOffset = .05;
sz_axOffset = 0.05; % extra breathing room

nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

flag_zeroMean = 0;
flag_plotIdvdCI = 0;
flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles
type_corr = 'pearson';
type_tail = 'both';

for iSet = 1:numel(iLocSingle_allSets)

    iLocCorr_all = iLocSingle_allSets{iSet};

    % Obtain CS (x-axis)
    X_allSubj = CS_allSubj(:, iLocCorr_all);
    x_ticks = linspace(1.5, 3.5, 5); % EE
    X_allBoot_allSubj = repmat(X_allSubj, 1, 1, nBoot);

    for iFeature = 1:nFeatures
        iFamily = iFamily_perF(iFeature);
        namesTunCs = namesTunC_unit_perF{iFamily, 2};
        nTunCs_full = length(namesTunCs);
        
        for iTunC = 1:nTunCs_full
            switch iFeature
                case 1, NOMp_allBoot_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
                case 2, NOMp_allBoot_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
            end

            switch iFeature
                case 1, y_ticks_allTunC_lb = [0, 0, -.1]; y_ticks_allTunC_ub = [.36, 90, .1];% 3 values are ORI peak amplitude, width, baseline
                case 2, y_ticks_allTunC_lb = [0, 0, 0, -.1]; y_ticks_allTunC_ub = [4, .24, 2.8, .1]; % 4 values are SF peak, peak amplitude, width, baseline
            end

            y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);

            x_ticklabels = x_ticks;
            y_ticklabels = y_ticks;

            nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily_perF(iFeature), 2}{iTunC});
            nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

            str_title = sprintf('n=%d, %s vs. %s [L%s]', nsubj, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

            %======================%
            basicFxn_drawCorr_permutation(X_allBoot_allSubj, NOMp_allBoot_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, type_corr, type_tail, str_title, markers_allSubj, nBoot);
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

            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_%s_L%s.png', nameFolder_Fig_NOM_corr, nsubj, nameVarY_fileTitle, strjoin(string(iLocCorr_all), '')))
            close(gcf)

        end % iTunC
    end % iFeature
end % iSet
clear margTunC*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Corr2: Corr between CS and NOM params
clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond')

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
type_corr = 'pearson';
type_tail = 'both';

for iSet = 1:numel(iLocSingle_allSets)

    iLocCorr_all = iLocSingle_allSets{iSet};

    % Obtain CS (x-axis)
    X_allSubj = CS_allSubj(:, iLocCorr_all);
    x_ticks = linspace(1.5, 3.5, 5); % EE
    X_allBoot_allSubj = repmat(X_allSubj, 1, 1, nBoot);
    % X_med_allSubj = X_allSubj;

    for iNOMparam = 1:nNOMparams
        NOMp_allBoot_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iNOMparam));

        Y_med_allSubj = getCI(NOMp_allBoot_allSubj, 1, 3)'; % rotate to match the format needed by basicFxn_drawCorr

        y_ticks_lb = [0, 0, 0]; y_ticks_ub = [.8, 40, 40];

        y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

        x_ticklabels = x_ticks;
        y_ticklabels = y_ticks;

        nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
        nameVarY_fileTitle = nameVarY_figTitle;

        str_title = sprintf('n=%d, %s vs. %s [L%s]', nsubj, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

        basicFxn_drawCorr_permutation(X_allBoot_allSubj, NOMp_allBoot_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, type_corr, type_tail, str_title, markers_allSubj, nBoot);

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

        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_NOM%d_L%s.png', nameFolder_Fig_NOM_corr, nsubj, iNOMparam, strjoin(string(iLocCorr_all), '')))
        close(gcf)
    end % iNOMparam
end % iSet
clear params_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Corr3: Corr between CS and pA
clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'metric_data_allCond')

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
type_corr = 'pearson';
type_tail = 'both';

for iSet = 1:numel(iLocSingle_allSets)

    iLocCorr_all = iLocSingle_allSets{iSet};

    % Obtain CS (x-axis)
    X_allSubj = CS_allSubj(:, iLocCorr_all);
    x_ticks = linspace(1.5, 3.5, 5); % EE
    X_allBoot_allSubj = repmat(X_allSubj, 1, 1, nBoot);

    NOMp_allBoot_allSubj = getCI(metric_data_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, iMetric_pA, :, :, :), 2, 7);

    y_ticks = round(linspace(.6, .9, 5), 2);

    x_ticklabels = x_ticks;
    y_ticklabels = y_ticks;

    nameVarY_figTitle = nameVarY;
    nameVarY_fileTitle = nameVarY_figTitle;

    str_title = sprintf('n=%d, %s vs. %s [L%s]', nsubj, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

    basicFxn_drawCorr_permutation(X_allBoot_allSubj, NOMp_allBoot_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, type_corr, type_tail, str_title, markers_allSubj, nBoot);

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

    drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%s.png', nameFolder_Fig_NOM_corr, nsubj, strjoin(string(iLocCorr_all), '')))
    close(gcf)

end % iSet
clear metric_data_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Corr4: Corr between pA and NOM params
clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond')
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
type_corr = 'pearson';
type_tail = 'both';

for iSet = 1:numel(iLocSingle_allSets)

    iLocCorr_all = iLocSingle_allSets{iSet};

    % Obtain x-axis values
    X_allSubj = pA_allSubj(:, iLocCorr_all);
    x_ticks = linspace(.6, .8, 5); % EE
    X_allBoot_allSubj = repmat(X_allSubj, 1, 1, nBoot);

    for iNOMparam = 1:nNOMparams
        NOMp_allBoot_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iNOMparam));

        y_ticks_lb = [0, 5, 0]; y_ticks_ub = [.8, 45, 40];

        y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

        x_ticklabels = x_ticks;
        y_ticklabels = y_ticks;

        nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
        nameVarY_fileTitle = nameVarY_figTitle;

        str_title = sprintf('n=%d, %s vs. %s [L%s]', nsubj, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

        basicFxn_drawCorr_permutation(X_allBoot_allSubj, NOMp_allBoot_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, type_corr, type_tail, str_title, markers_allSubj, nBoot);

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

        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_NOM%d_L%s.png', nameFolder_Fig_NOM_corr, nsubj, iNOMparam, strjoin(string(iLocCorr_all), '')))
        close(gcf)
    end % iNOMparam
end % iSet
clear params_allCond
fprintf('\n\n Plotting DONE\n\n')

%% CorrAsym1: (CS and tunC): corr between extents of EE/HVA/VMA
flag_plotIdvdCI = 1;
flag_plotUnikSymbol = 0;

clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'margTunC_*_allCond')

nameVarX = 'CS';
nameVarY = 'tunC';

nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

sz_label = 55;
sz_labelOffset = 0.01;
sz_axOffset = 0.05; % extra breathing room

for iGroup = 1:nGroups
    iLocPair_all = iLocGroups_all{iGroup};

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
    asymX_allBoot_allSubj = repmat(asymX_allSubj, 1, nBoot)';

    for iFeature = 1:nFeatures
        iFamily = iFamily_perF(iFeature);
        namesTunCs = namesTunC_unit_perF{iFamily, 2};
        nTunCs_full = length(namesTunCs);

        for iTunC = 1:nTunCs_full
            % Y-axis
            switch iFeature
                case 1, NOMp_allBoot_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                case 2, NOMp_allBoot_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
            end

            asymY_allBoot_allSubj = squeeze((NOMp_allBoot_allSubj(1, :, :)-NOMp_allBoot_allSubj(2, :, :))./(NOMp_allBoot_allSubj(1, :, :)+NOMp_allBoot_allSubj(2, :, :)));

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

            str_title = sprintf('n=%d, %s (%s) vs. %s (%s)', nsubj, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            basicFxn_drawCorrAsym_permutation(asymX_allBoot_allSubj*100, asymY_allBoot_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)
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

            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_%s.png', nameFolder_Fig_NOM_CorrAsym, nsubj, iLocPair_all, nameVarY_fileTitle))
            close(gcf)

        end % iTunC
    end % iFeature
end % iGroup
clear margTunC_*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% CorrAsym2: CorrAsym (CS and NOMparams): corr between extents of EE/HVA/VMA
clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond')

nameVarX = 'CS';
nameVarY = 'NOMparams';

nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

sz_label = 60;
sz_labelOffset = 0.05;
sz_axOffset = 0.05; % extra breathing room

for iGroup = 1:nGroups
    iLocPair_all = iLocGroups_all{iGroup};

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
    asymX_allBoot_allSubj = repmat(asymX_allSubj, 1, nBoot)';

    for iNOMparam = 1:nNOMparams
        NOMp_allBoot_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iNOMparam));

        asymY_allBoot_allSubj = squeeze((NOMp_allBoot_allSubj(1, :, :)-NOMp_allBoot_allSubj(2, :, :))./(NOMp_allBoot_allSubj(1, :, :)+NOMp_allBoot_allSubj(2, :, :)));

        y_ticks_lb = -[30, 20, 20]; y_ticks_ub = [50 40 40];
        y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

        x_ticklabels = nan;
        y_ticklabels = nan;

        nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
        nameVarY_fileTitle = nameVarY_figTitle;

        str_title = sprintf('n=%d, %s (%s) vs. %s (%s)', nsubj, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        basicFxn_drawCorrAsym_permutation(asymX_allBoot_allSubj*100, asymY_allBoot_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)
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

        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_NOM%d.png', nameFolder_Fig_NOM_CorrAsym, nsubj, iLocPair_all, iNOMparam))
        close(gcf)

    end % iNOMparam
end % iGroup
clear params_allCond
fprintf('\n\n Plotting DONE\n\n')

%% CorrAsym3: (CS and pA): corr between extents of EE/HVA/VMA
clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'metric_data_allCond')

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
    asymX_allBoot_allSubj = repmat(asymX_allSubj, 1, nBoot)';

    NOMp_allBoot_allSubj = getCI(metric_data_allCond(iModelA_plot, iModelB_plot, iLocPair_all, iMetric_pA, :, :, :), 2, 7);

    asymY_allBoot_allSubj = squeeze((NOMp_allBoot_allSubj(1, :, :)-NOMp_allBoot_allSubj(2, :, :))./(NOMp_allBoot_allSubj(1, :, :)+NOMp_allBoot_allSubj(2, :, :)));

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

    str_title = sprintf('n=%d, %s (%s) vs. %s (%s)', nsubj, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    basicFxn_drawCorrAsym_permutation(asymX_allBoot_allSubj*100, asymY_allBoot_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)
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

    drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d.png', nameFolder_Fig_NOM_CorrAsym, nsubj, iLocPair_all))
    close(gcf)

end % iGroup
clear metric_data_allCond
fprintf('\n\n Plotting DONE\n\n')

%% CorrAsym4: (pA and NOMparams): corr between extents of EE/HVA/VMA
clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond')

nameVarX = 'pA';
nameVarY = 'NOMparams';

nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

sz_label = 60;
sz_labelOffset = 0.05;
sz_axOffset = 0.05; % extra breathing room

for iGroup = 1:nGroups
    iLocPair_all = iLocGroups_all{iGroup};

    switch iLocPair_all(1)
        case 1, nameAsymX = 'EE'; nameAsymY = 'EE';
        case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
        case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
    end

    % X-axis (pA)
    asymX_allSubj = (pA_allSubj(:, iLocPair_all(1))-pA_allSubj(:, iLocPair_all(2)))./(pA_allSubj(:, iLocPair_all(1))+pA_allSubj(:, iLocPair_all(2)));
    asymX_allBoot_allSubj = repmat(asymX_allSubj, 1, nBoot)';
    % switch iLocPair_all(1)
    %     case 1, x_ticks = linspace(0, 20, 5); % EE
    %     case 6, x_ticks = linspace(-5, 15, 5); % HVA
    %     case 5, x_ticks = linspace(0, 16, 5); % VMA (extent is smaller)
    % end
    x_ticks = linspace(-6, 10, 5);

    for iNOMparam = 1:nNOMparams
        NOMp_allBoot_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iNOMparam));

        asymY_allBoot_allSubj = squeeze((NOMp_allBoot_allSubj(1, :, :)-NOMp_allBoot_allSubj(2, :, :))./(NOMp_allBoot_allSubj(1, :, :)+NOMp_allBoot_allSubj(2, :, :)));

        % switch iLocPair_all(1)
        %     case 1 % EE
        %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        %     case 6 % HVA
        %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        %     case 5 % VMA
        %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        % end
        y_ticks_lb = -[60, 40, 20]; y_ticks_ub = [40, 40, 60];
        y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

        x_ticklabels = nan;
        y_ticklabels = nan;

        nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
        nameVarY_fileTitle = nameVarY_figTitle;

        str_title = sprintf('n=%d, %s (%s) vs. %s (%s)', nsubj, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        basicFxn_drawCorrAsym_permutation(asymX_allBoot_allSubj*100, asymY_allBoot_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)
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

        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_NOM%d.png', nameFolder_Fig_NOM_CorrAsym, nsubj, iLocPair_all, iNOMparam))
        close(gcf)

    end % iNOMparam
end % iGroup
clear params_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Plot metrics vs. IV for group averages (single locations)
clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'IV_allCond', 'nTrials_allCond', 'metric_*_allCond', 'R2_NOM_allCond')

% Define folder for saving figures
nameFolder_Fig_NOM_Trialwise = sprintf('%s/NOM_Trialwise_%d%d', nameFolder_Figures, nORI, nSF);
if isempty(dir(nameFolder_Fig_NOM_Trialwise)), mkdir(nameFolder_Fig_NOM_Trialwise), end
namesModelB_plot = {'Full model', 'No $\sigma_{shared}$', 'No $\sigma_{mul}$', 'No $\sigma_{private}$'};

x_base   = 0.60;   % move a bit further left to make space for name column
x_width  = 0.4;   % total width of the table block
y_base   = 0.10;
y_height = 0.25;

sz_font = 15;
sz_axis = 18;
lineStyle_all = {'-', '-', '--', ':', ':', '-', '--', ':', '-.'};
namesPlotMode = {'NoPred', 'FullOnly', 'AllModels'};

for iPlotMode = 1:length(namesPlotMode)

    % Decide the number of models to plot
    if iPlotMode == 2
        iModelB_all_plot = 1;
    else
        % iModelB_all_plot = 1:4;
        iModelB_all_plot = iModelB_all;
    end
    nModelB_plot = length(iModelB_all_plot);

    for iModelA = 1%iModelA_all % just plot pred of the model using RC-derived template (happens to be yhe best model)

        % Define folder for saving figures for each model A and model B
        nameFolder_Fig_NOM_metrics = sprintf('%s/A%d', nameFolder_Fig_NOM_Trialwise, iModelA);
        if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end


        for iSet = 1:numel(iLocSingle_allSets)

            iLocSingle_all = iLocSingle_allSets{iSet};

            fprintf('\nL%s...', strjoin(string(iLocSingle_all), ''))
            scalingF = 80; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end

            figure('Position', [0 0 500 1e3])

            for iMetric_prob = 1:nMetrics_prob

                subplot(nMetrics_prob, 1, iMetric_prob), hold on

                R2_tab = nan(nModelB_plot, 2);   % [row = model, col = loc within pair]

                for iModelB = iModelB_all_plot

                    if iMetric_prob==1
                        fprintf('%s: %s | ', namesModelB{iModelB_all_plot(iModelB)}, lineStyle_all{iModelB})
                    end

                    for iiLoc = 1:numel(iLocSingle_all)
                        % Compute medians and CIs
                        [IV_ave, ~, ~, IV_SEM] = getCI(getCI(IV_allCond(iModelA, iModelB, iLocSingle_all(iiLoc), :, :, :), 1, 5), 2, 1);
                        [nTrials_ave, ~, ~, nTrials_SEM] = getCI(getCI(nTrials_allCond(iModelA, iModelB, iLocSingle_all(iiLoc), :, :, :), 1, 5), 2, 1);
                        [data_ave, ~, ~, data_SEM] = getCI(getCI(metric_data_allCond(iModelA, iModelB, iLocSingle_all(iiLoc), iMetric_prob, :, :, :), 1, 6), 2, 1);
                        [pred_ave, ~, ~, pred_SEM] = getCI(getCI(metric_pred_allCond(iModelA, iModelB, iLocSingle_all(iiLoc), iMetric_prob, :, :, :), 1, 6), 2, 1);
                        [R2_NOM_ave, ~, ~, R2_NOM_SEM] = getCI(getCI(R2_NOM_allCond(iModelA, iModelB, iLocSingle_all(iiLoc), iMetric_prob, :, :), 1, 6), 2, 1);

                        % ---- Store R2 mean into table ----
                        rowIdx = find(iModelB_all_plot == iModelB);  % which row for this model
                        colIdx = iiLoc;                              % col 1 or 2 for the pair
                        R2_tab(rowIdx, colIdx) = R2_NOM_ave;        % store mean R2 (not SEM)

                        if iModelB==1, color_pred = colors_comb(iLocSingle_all(iiLoc), :); lw = 3;
                        else, color_pred='k'; lw = 1.5;
                        end

                        if iPlotMode ~= 1
                            % Plot prediction
                            patch([IV_ave; flip(IV_ave)], [pred_ave-pred_SEM; flip(pred_ave+pred_SEM)], ones(1, 3) / 2, 'FaceAlpha', .3, 'linestyle', 'none', 'handlevisibility', 'off')
                            plot(IV_ave, pred_ave, 'lineStyle', lineStyle_all{iModelB}, 'color', colors_comb(iLocSingle_all(iiLoc), :), 'linewidth', lw)
                        end

                        % Plot measurement
                        % errorbar(IV_ave, data_ave, IV_SEM, '.', 'horizontal', 'CapSize', 0, 'color', colors_comb(iLoc_Sep_all(iiLoc), :), 'handlevisibility', 'off')
                        errorbar(IV_ave, data_ave, data_SEM, '.', 'vertical', 'CapSize', 0, 'color', colors_comb(iLocSingle_all(iiLoc), :), 'handlevisibility', 'off', 'linewidth', lw/1.5)

                        % Plot averaged data of each bin (dot size indicates number of trials)
                        for iBin = 1:nBins
                            facecolor = 'w';
                            plot(IV_ave(iBin), data_ave(iBin), 'o', 'markeredgecolor', colors_comb(iLocSingle_all(iiLoc), :), ...
                                'markerfacecolor', facecolor, 'MarkerSize', nTrials_ave(iBin) / scalingF + 5, 'LineWidth', lw/1.5, 'LineStyle', 'none', 'handlevisibility', 'off')
                        end
                    end % iiLoc

                    yline(.5, '--', 'linewidth', 2, 'color', ones(1,3)/2);
                    ylim([0, 1])
                    yticks(0:.2:1)
                    ylabel(namesMetrics_prob_full{iMetric_prob})

                    if find(iLocSingle_all==1), x_ticks = linspace(0, 180, 5);
                    else, x_ticks = linspace(0, 100, 5);
                    end
                    xticks(x_ticks)
                    xlim(x_ticks([1, end]))
                    xlabel('Binned decision variable')
                    if isubj == 1, legend(namesModelB(iModelB_all_plot), 'Location', 'best'), end

                end % iModelB

                % ==== Print R2 as a 4x3 table (left column = model name) ====
                ax = gca;

                nRows_subj = nModelB_plot;
                nCols= 3;   % now: [ModelName | Loc1 | Loc2]

                % Create normalized coordinates for each cell
                x_pos = linspace(x_base, x_base + x_width, nCols + 1);
                x_cells = x_pos(1:nCols) + diff(x_pos(1:2))/2;   % center of each column

                y_pos = linspace(y_base, y_base + y_height, nRows_subj + 1);
                y_cells = fliplr(y_pos(1:nRows_subj)) + diff(y_pos(1:2))/2;

                if iPlotMode ~= 1
                    % ---- Column 1: model names (B1, B2, ...) ----
                    for r = 1:nRows_subj
                        name_str = sprintf('%s', namesModelB_plot{iModelB_all_plot(r)});

                        text(x_cells(1), y_cells(r), name_str, 'Units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize', sz_font, 'FontWeight', 'bold', 'interpreter', 'latex');
                        % text(x_cells(1), y_cells(r), '$R^2$=', 'Units','normalized', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', sz_font, 'FontWeight', 'bold', 'interpreter', 'latex');
                    end

                    % ---- Columns 2–3: R2 values ----
                    for r = 1:nRows_subj
                        for c = 1:2   % two locations
                            R2_val = R2_tab(r, c);
                            if isnan(R2_val), continue; end

                            str_cell = sprintf('%.0f%%', R2_val * 100);

                            % text(x_cells(c+1), y_cells(r), str_cell, 'Units','normalized', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', sz_font);
                        end
                    end
                end % if iPlotMode==2
                ax = gca;
                ax.FontSize=sz_axis;
                ax.LineWidth = 1.5; % Adjust the value as desired
            end % iMetric

            sgtitle(sprintf('n=%d [A%d] [L%s] [nBoot=%d] %s', nsubj, iModelA, strjoin(string(iLocSingle_all), ''), nBoot, namesMetrics_prob{iMetric_prob}))

            % Save the figure
            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%s_A%d_group_%s.png', nameFolder_Fig_NOM_metrics, nsubj, strjoin(string(iLocSingle_all), ''), iModelA_plot, namesPlotMode{iPlotMode}))
            close(gcf)
            fprintf('DONE\n')
        end % iSet
    end % iModelA
end % iPlotMode
clear *_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Plot metrics vs. IV for group averages (paired locations)
clc, fprintf('\n\n Plotting STARTED\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'IV_allCond', 'nTrials_allCond', 'metric_*_allCond', 'R2_NOM_allCond')

% Define folder for saving figures
nameFolder_Fig_NOM_Trialwise = sprintf('%s/NOM_Trialwise_%d%d', nameFolder_Figures, nORI, nSF);
if isempty(dir(nameFolder_Fig_NOM_Trialwise)), mkdir(nameFolder_Fig_NOM_Trialwise), end
namesModelB_plot = {'Full model', 'No $\sigma_{shared}$', 'No $\sigma_{mul}$', 'No $\sigma_{private}$'};

x_base   = 0.60;   % move a bit further left to make space for name column
x_width  = 0.4;   % total width of the table block
y_base   = 0.10;
y_height = 0.25;

sz_font = 15;
sz_axis = 18;
lineStyle_all = {'-', '-', '--', ':', ':', '-', '--', ':', '-.'};
namesPlotMode = {'NoPred', 'FullOnly', 'AllModels'};

for iPlotMode = 1:length(namesPlotMode)

    % Decide the number of models to plot
    if iPlotMode == 2
        iModelB_all_plot = 1;
    else
        % iModelB_all_plot = 1:4;
        iModelB_all_plot = iModelB_all;
    end
    nModelB_plot = length(iModelB_all_plot);

    for iModelA = 1%iModelA_all % just plot pred of the model using RC-derived template (happens to be yhe best model)

        % Define folder for saving figures for each model A and model B
        nameFolder_Fig_NOM_metrics = sprintf('%s/A%d', nameFolder_Fig_NOM_Trialwise, iModelA);
        if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end

        % for iLocComb = iLocComb_all
        for iGroup = 1:nGroups
            iLocPair_all = iLocGroups_all{iGroup};

            fprintf('\nL%d%d...', iLocPair_all)
            scalingF = 80; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end

            figure('Position', [0 0 500 1e3])

            for iMetric_prob = 1:nMetrics_prob

                subplot(nMetrics_prob, 1, iMetric_prob), hold on

                R2_tab = nan(nModelB_plot, 2);   % [row = model, col = loc within pair]

                for iModelB = iModelB_all_plot

                    if iMetric_prob==1
                        fprintf('%s: %s | ', namesModelB{iModelB_all_plot(iModelB)}, lineStyle_all{iModelB})
                    end

                    for iiLoc = 1:2
                        % Compute medians and CIs
                        [IV_ave, ~, ~, IV_SEM] = getCI(getCI(IV_allCond(iModelA, iModelB, iLocPair_all(iiLoc), :, :, :), 1, 5), 2, 1);
                        [nTrials_ave, ~, ~, nTrials_SEM] = getCI(getCI(nTrials_allCond(iModelA, iModelB, iLocPair_all(iiLoc), :, :, :), 1, 5), 2, 1);
                        [data_ave, ~, ~, data_SEM] = getCI(getCI(metric_data_allCond(iModelA, iModelB, iLocPair_all(iiLoc), iMetric_prob, :, :, :), 1, 6), 2, 1);
                        [pred_ave, ~, ~, pred_SEM] = getCI(getCI(metric_pred_allCond(iModelA, iModelB, iLocPair_all(iiLoc), iMetric_prob, :, :, :), 1, 6), 2, 1);
                        [R2_NOM_ave, ~, ~, R2_NOM_SEM] = getCI(getCI(R2_NOM_allCond(iModelA, iModelB, iLocPair_all(iiLoc), iMetric_prob, :, :), 1, 6), 2, 1);

                        % ---- Store R2 mean into table ----
                        rowIdx = find(iModelB_all_plot == iModelB);  % which row for this model
                        colIdx = iiLoc;                              % col 1 or 2 for the pair
                        R2_tab(rowIdx, colIdx) = R2_NOM_ave;        % store mean R2 (not SEM)

                        if iModelB==1, color_pred = colors_comb(iLocPair_all(iiLoc), :); lw = 3;
                        else, color_pred='k'; lw = 1.5;
                        end

                        if iPlotMode ~= 1
                            % Plot prediction
                            patch([IV_ave; flip(IV_ave)], [pred_ave-pred_SEM; flip(pred_ave+pred_SEM)], ones(1, 3) / 2, 'FaceAlpha', .3, 'linestyle', 'none', 'handlevisibility', 'off')
                            plot(IV_ave, pred_ave, 'lineStyle', lineStyle_all{iModelB}, 'color', colors_comb(iLocPair_all(iiLoc), :), 'linewidth', lw)
                        end

                        % Plot measurement
                        % errorbar(IV_ave, data_ave, IV_SEM, '.', 'horizontal', 'CapSize', 0, 'color', colors_comb(iLocPair_all(iiLoc), :), 'handlevisibility', 'off')
                        errorbar(IV_ave, data_ave, data_SEM, '.', 'vertical', 'CapSize', 0, 'color', colors_comb(iLocPair_all(iiLoc), :), 'handlevisibility', 'off', 'linewidth', lw/1.5)

                        % Plot averaged data of each bin
                        for iBin = 1:nBins
                            facecolor = 'w';
                            plot(IV_ave(iBin), data_ave(iBin), 'o', 'markeredgecolor', colors_comb(iLocPair_all(iiLoc), :), ...
                                'markerfacecolor', facecolor, 'MarkerSize', nTrials_ave(iBin) / scalingF + 5, 'LineWidth', lw/1.5, 'LineStyle', 'none', 'handlevisibility', 'off')
                        end
                    end % iiLoc

                    yline(.5, '--', 'linewidth', 2);

                    % xlim([,150])
                    ylim([0, 1])
                    yticks(0:.2:1)
                    ylabel(namesMetrics_prob_full{iMetric_prob})

                    switch iLocPair_all(1)
                        case 1, xlim([0, 200])
                        case 6, xlim([0, 100])
                        case 5, xlim([0, 80])
                    end
                    xlabel('Binned decision variable')
                    if isubj == 1, legend(namesModelB(iModelB_all_plot), 'Location', 'best'), end

                end % iModelB

                % ==== Print R2 as a 4x3 table (left column = model name) ====
                ax = gca;

                nRows_subj = nModelB_plot;
                nCols= 3;   % now: [ModelName | Loc1 | Loc2]

                % Create normalized coordinates for each cell
                x_pos = linspace(x_base, x_base + x_width, nCols + 1);
                x_cells = x_pos(1:nCols) + diff(x_pos(1:2))/2;   % center of each column

                y_pos = linspace(y_base, y_base + y_height, nRows_subj + 1);
                y_cells = fliplr(y_pos(1:nRows_subj)) + diff(y_pos(1:2))/2;

                if iPlotMode ~= 1
                    % ---- Column 1: model names (B1, B2, ...) ----
                    for r = 1:nRows_subj
                        name_str = sprintf('%s', namesModelB_plot{iModelB_all_plot(r)});

                        text(x_cells(1), y_cells(r), name_str, 'Units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize', sz_font, 'FontWeight', 'bold', 'interpreter', 'latex');
                        % text(x_cells(1), y_cells(r), '$R^2$=', 'Units','normalized', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', sz_font, 'FontWeight', 'bold', 'interpreter', 'latex');
                    end

                    % ---- Columns 2–3: R2 values ----
                    for r = 1:nRows_subj
                        for c = 1:2   % two locations
                            R2_val = R2_tab(r, c);
                            if isnan(R2_val), continue; end

                            str_cell = sprintf('%.0f%%', R2_val * 100);

                            % text(x_cells(c+1), y_cells(r), str_cell, 'Units','normalized', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', sz_font);
                        end
                    end
                end % if iPlotMode==2
                ax = gca;
                ax.FontSize=sz_axis;
                ax.LineWidth = 1.5; % Adjust the value as desired
            end % iMetric

            sgtitle(sprintf('n=%d [A%d] [L%d%d] [nBoot=%d] %s', nsubj, iModelA, iLocPair_all, nBoot, namesMetrics_prob{iMetric_prob}))
            % set(findall(gcf, '-property', 'fontsize'), 'fontsize', 15)
            % set(findall(gcf, '-property', 'linewidth'), 'linewidth', 2)

            % Save the figure
            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_A%d_group_%s.png', nameFolder_Fig_NOM_metrics, nsubj, iLocPair_all, iModelA_plot, namesPlotMode{iPlotMode}))
            close(gcf)
            fprintf('DONE\n')
        end % iLocComb
    end % iModelA
end % iPlotMode

clear *_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Plot metrics vs. IV for each idvd
clc, fprintf('\n\n Plotting STARTED\n\n')
lineStyle_all = {'-', '-', '--', ':', '-.', '-', '--', ':', '-.'};

% iModelA_plot = 1;%iModelA_all % just plot pred of the model using RC-derived template (happens to be yhe best model)
% iModelA = iModelA_plot;
iModelB_all_plot = 1;

% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond', 'IV_allCond', 'nTrials_allCond', 'metric_data_allCond', 'metric_pred_allCond', 'R2_NOM_allCond')
% Define folder for saving figures for each model A and model B
nameFolder_Fig_NOM_metrics = sprintf('%s/A%d', nameFolder_Fig_NOM_Trialwise, iModelA_plot);
if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end

for iLocSingle = iLocComb_all
    fprintf('\nL%d...', iLocSingle)
    scalingF = 80; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end

    %%%%%%%%%%%%
    for iMetric_prob = 1:nMetrics_prob
        figure('Position', [0 0 2e3 2e3])

        for isubj = 1:nsubj
            subplot(nRows_subj, nCols_subj, isubj), hold on

            subjName = subjList{isubj};

            % Create a string to store R2
            str_R2 = '|';

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
                        'markerfacecolor', facecolor, 'MarkerSize', nTrials_allBins(iBin) / scalingF + 5, 'LineWidth', 1, 'LineStyle', 'none', 'handlevisibility', 'off')
                end

                % xlim([,150])
                ylim([0, 1])
                ylabel(namesMetrics_prob{iMetric_prob})
                xlabel('Binned decision variable')
                if isubj == 1, legend(namesModelB(iModelB_all), 'Location', 'best'), end

                % R2 median and CI
                str_R2 = sprintf('%s B%d: %.2f [%.2f, %.2f]', str_R2, iModelB, R2_NOM_med, R2_NOM_lb, R2_NOM_ub);
            end % iModelB

            % print estimates
            text(100, .2, str_est, 'FontSize', 6)
            title(sprintf('%s %s', subjName, str_R2))
        end % isubj

        sgtitle(sprintf('n=%d [A%d] [L%d] [nBoot=%d] %s', nsubj, iModelA_plot, iLocSingle, nBoot, namesMetrics_prob{iMetric_prob}))
        set(findall(gcf, '-property', 'fontsize'), 'fontsize', 12)

        % Save the figure
        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d_A%d_%s.png', nameFolder_Fig_NOM_metrics, nsubj, iLocSingle, iModelA_plot, namesMetrics_prob{iMetric_prob}))
        close(gcf)
    end % iMetric

    fprintf('DONE\n')
end % iLocComb

% end % iModelA

fprintf('\n\n Plotting DONE\n\n')

%% Plot GoF for ModelA x ModelB x Loc (use basicFxn_drawBars_permutation)
clc
clc, fprintf('\n\n Plotting STARTED\n\n')
% Load data
load(nameFolder_Data_SaveCompile, 'nLL_allCond')

% Define folder for saving GoF figures
nameFolder_Fig_NOM_GoF = fullfile(nameFolder_Fig_NOM_Trialwise, 'GoF');
if isempty(dir(nameFolder_Fig_NOM_GoF)), mkdir(nameFolder_Fig_NOM_GoF); end

% iModelB_selected = [1,3,2,4]; flag_plotIDVD = 0; flag_plotDiff = 0;
iModelB_selected = [1,2,3,4]; flag_plotIDVD = 0; flag_plotDiff = 0;
iModelB_selected = [1,3,4]; flag_plotIDVD = 1; flag_plotDiff = 1; %  {'FullModel'}  {'NoMultiN'}    {'NoPrivN'}  
% iModelB_selected = [1,3]; flag_plotIDVD = 1; flag_plotDiff = 1; %  {'FullModel'}  {'NoMultiN'}

% Goodness-of-fit (GoF) measures to plot
namesGoF = {'nLL', 'AIC-nLL', 'AICc-nLL', 'BIC-nLL'}; % no raw nLL
nGoFs = numel(namesGoF);

% Common y-limit for all GoF measures (already positive deltas)
y_lim = [0, 16];
sz_label = 20;
% Loop over GoF measures
for iGoF = 1% 1:nGoFs % just plot nLL

    for iModelA = 1%iModelA_all

        for iLocSingle = 1:nLocComb8

            sz_wd_perBar = 70;
            nBars = numel(iModelB_all(iModelB_selected));
            sz_fig = [nBars*sz_wd_perBar, 250];

            % data: [ModelA x ModelB x LocComb x Metric x Subj x Iteration]
            nLL_allBoot_allSubj = nLL_allCond(iModelA, iModelB_selected, iLocSingle, :, :);
            dnLL_allBoot_allSubj = nan(size(nLL_allBoot_allSubj));

            % For each subject and each bootstrap, subtract the minimum across models
            for isubj = 1:nsubj
                parfor iBoot = 1:nBoot
                    d = nLL_allBoot_allSubj(:, :, :, isubj, iBoot);
                    d_min = min(d(:));
                    dnLL_allBoot_allSubj(:, :, :, isubj, iBoot) = d - d_min;
                end
            end

            dnLL_allBoot_allSubj = squeeze(dnLL_allBoot_allSubj);

            str_title = sprintf('n=%d %s [L%d] ModelB [%s]', nsubj, namesGoF{iGoF}, iLocSingle, strjoin(string(iModelB_selected), ' '));
            x_ticks = nan;
            ref = nan;

            %------------------------------%
            basicFxn_drawBars_permutation(dnLL_allBoot_allSubj, ref, repmat(colors_comb(iLocSingle, :), nBars, 1), namesModelB(iModelB_all(iModelB_selected)), x_ticks, x_ticks, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nBoot, markers_allSubj)
            %------------------------------%
            % ylim(y_lim)
            ylabel(sprintf('\\Delta nLL'), 'FontSize', sz_label);

            drawnow, pause(1), saveas(gcf, fullfile(nameFolder_Fig_NOM_GoF, sprintf('n%d_%s_B%s_L%d.png', nsubj, namesGoF{iGoF}, strjoin(string(iModelB_selected), ''), iLocSingle)));
            close(gcf)

        end % iLocSingle
    end % iModelA
end % iGoF

clear nLL_allCond
fprintf('\n\n Plotting DONE\n\n')

%% ANOVA on models
clc
iModelB_ANOVA_all = [1,3];

% ANOVA
nModelB_ANOVA = length(iModelB_ANOVA_all);
GoF_perBoot = squeeze(GoF_allCond(1, iModelB_ANOVA_all, iLocComb_all, :, 1));
IV_NOMmodel = nan(size(GoF_perBoot));
IV_Loc = IV_NOMmodel;
for iModelB = 1:nModelB_ANOVA
    IV_NOMmodel(iModelB, :, :) = ones(length(iLocComb_all), nsubj)*iModelB_ANOVA_all(iModelB);
end
for iiLoc = 1:length(iLocComb_all)
    IV_Loc(:, iiLoc, :) = ones(nModelB_ANOVA, nsubj)*iLocComb_all(iiLoc);
end

% Placeholders for ANOVA
p_modelB_allBoot = nan(nBoot, 1);
p_loc_allBoot = p_modelB_allBoot;
p_interaction_allBoot = p_modelB_allBoot;

% Placeholders for t-test
t_allLoc_allBoot = nan(length(iLocComb_all), nBoot);
p_ttest_allLoc_allBoot = t_allLoc_allBoot;
diff_allLoc_allBoot = t_allLoc_allBoot;

for iBoot= 1:nBoot
    GoF_perBoot = squeeze(GoF_allCond(1, iModelB_ANOVA_all, iLocComb_all, :, iBoot));

    [~, tbl] = print_nANOVA({'NOMmodel', 'Loc'}, GoF_perBoot(:), {IV_NOMmodel(:), IV_Loc(:)}, nsubj);
    p_modelB_allBoot(iBoot) = tbl{2,6};
    p_loc_allBoot(iBoot) = tbl{3,6};
    p_interaction_allBoot(iBoot) = tbl{4,6};

    % conduct t-test
    if nModelB_ANOVA==2
        for iiLoc = 1:length(iLocComb_all)
            [~, p_ttest, ~, stats] = ttest(GoF_perBoot(1, iiLoc, :), GoF_perBoot(2, iiLoc, :));
            t_allLoc_allBoot(iiLoc, iBoot) = stats.tstat;
            p_ttest_allLoc_allBoot(iiLoc, iBoot) = p_ttest;
            diff_allLoc_allBoot(iiLoc, iBoot) = mean(GoF_perBoot(1, iiLoc, :)-GoF_perBoot(2, iiLoc, :));
        end
    end
end

% Print median and CI of p-values
clc
fprintf('\n\nModelB ANOVA: B%d-%s vs. B%d-%s\n', iModelB_ANOVA_all(1), namesModelB{iModelB_ANOVA_all(1)}, iModelB_ANOVA_all(2), namesModelB{iModelB_ANOVA_all(2)})
for iTest = 1:3
    switch iTest
        case 1
            p_allBoot = p_modelB_allBoot;
            str_name = 'ModelB';
        case 2
            p_allBoot = p_loc_allBoot;
            str_name = 'Loc';
        case 3
            p_allBoot = p_interaction_allBoot;
            str_name = 'Interaction';
    end
    [p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 1);
    fprintf('   %s effect: p=%.3f [CI=%.3f, %.3f]\n', str_name, p_med, p_lb, p_ub);
end

% Print t-test
if nModelB_ANOVA==2
    fprintf('\n\nModelB t-test: B%d-%s vs. B%d-%s\n', iModelB_ANOVA_all(1), namesModelB{iModelB_ANOVA_all(1)}, iModelB_ANOVA_all(2), namesModelB{iModelB_ANOVA_all(2)})
    if p_med<.05 % If interation is significant
        for iiLoc = 1:length(iLocComb_all)
            fprintf('*L%d*\n', iLocComb_all(iiLoc))
            [t_med, t_lb, t_ub] = getCI(t_allLoc_allBoot(iiLoc, :), 1, 2);
            [p_med, p_lb, p_ub] = getCI(p_ttest_allLoc_allBoot(iiLoc, :)*length(iLocComb_all), 1, 2);
            [diff_med, diff_lb, diff_ub] = getCI(diff_allLoc_allBoot(iiLoc, :), 1, 2);
            fprintf('   t = %.3f [CI=%.3f, %.3f]\n', t_med, t_lb, t_ub);
            fprintf('   p = %.3f [CI=%.3f, %.3f]\n', p_med, p_lb, p_ub);
            fprintf('   diff = %.3f [CI=%.3f, %.3f]\n', diff_med, diff_lb, diff_ub);
        end % iiLoc
    end
end

%% Compare NOM parameters across locations (use basicFxn)
clc, fprintf('\n\n Plotting STARTED\n\n')
% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond')

% Define folder for saving comparison figures
nameFolder_Fig_NOM_compParams = fullfile(nameFolder_Fig_NOM_Trialwise, 'CompNOMparams');
if isempty(dir(nameFolder_Fig_NOM_compParams))
    mkdir(nameFolder_Fig_NOM_compParams);
end

% nParamsMax = 3; % max three params
sz_wd_perBar = 120;
nBars = 2;
sz_fig = [nBars*sz_wd_perBar, 250];
flag_plotIDVD = 1;
flag_plotDiff = 1;

for iModelB = iModelB_plot %iModelB_all % just plot the best model

    % Number of parameters for this Model B
    nParamsB = numel(namesModelBparams{iModelB});

    for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
        iLocPair_all = iLocGroups_all{iGroup};

        for iParam = 1:nParamsB

            NOMp_allBoot_allSubj = squeeze(params_allCond(iModelA_plot, iModelB, iLocPair_all, :, :, iParam));

            x_ticks = namesLocComb(iLocPair_all);
            switch iModelB
                case 1
                    switch iParam
                        case 1, y_ticks = linspace(NOMp1_lb, .8, 5);
                        case 2, y_ticks = linspace(NOMp2_lb, 40, 5);
                        case 3, y_ticks = linspace(NOMp3_lb, 40, 5);
                    end
            end
            y_ticklabels = round(y_ticks, 2);

            str_title = sprintf('n=%d L%d%d [A%dB%d] %s', nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesModelBparams{iModelB}{iParam});
            %------------------------------%
            basicFxn_drawBars_permutation(NOMp_allBoot_allSubj, ref, colors_comb(iLocPair_all, :), x_ticks, y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nBoot, markers_allSubj)
            %------------------------------%
            ylabel(sprintf('%s', namesModelBparams{iModelB}{iParam}))
            % Save the figure
            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_A%dB%d_NOMp%d.png', nameFolder_Fig_NOM_compParams, nsubj, iLocPair_all, iModelA_plot, iModelB_plot, iParam))
            close(gcf);
        end % iParam
    end % iGroup
end % iModelB
clear params_allCond
fprintf('\n\n Plotting DONE\n\n')
