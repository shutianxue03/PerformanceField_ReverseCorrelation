
% SX_analysis5_NOM_Trialwise.m
% Author: Shutian Xue
% Purpose: This script compiles, analyzes and plots boostrapped data for each observer saved Data_NOM_Trialwise

clear, clc, close all
set(0, 'DefaultFigureVisible', 'off')

addpath(genpath('PF_RC/Codes/')) % manually add to save time
addpath(genpath('PF_RC/Codes/fxn_analysis_RC_v2')) % manually add to save time

% Setting parameters for the analysis
%----------------
SX_RC1_setting
%----------------

set(0, 'DefaultFigureVisible', 'off') % avoid printing figures on the desktop

% Define the settings for the analysis
nIter = 200; % Number of iterations
nJob = 5;
nIterxJob = nIter*nJob;
iModelA_fit_all = [1]; % 1=RC-derived template, 2=use IO template, 2=randomize template,
iModelB_fit_all = [1:7]; % see SX_RC1_setting for namesModelB
flag_plotPurpose = 'slide'; % 'slide' or 'paper'

% Define model families for tuning curves
iFamily_ORI = 10; %1=scaled Gaussian
iFamily_SF = 2; % 2=log parabola
iFamily_perF = [iFamily_ORI, iFamily_SF];

% Use compiled SF tuning characteristics. For family 14, store one combined width in octaves.
namesTunC_SF_unit = namesTunC_unit_perF{iFamily_SF, 2};
namesTunC_SF_noUnit = namesTunC_noUnit{iFamily_SF, 2};
if iFamily_SF == 14
    namesTunC_SF_unit = {'peak SF (cpd)', 'amplitude (a.u.)', 'width (octaves)', 'baseline (a.u.)'};
    namesTunC_SF_noUnit = {'peak SF', 'amplitude', 'width', 'baseline'};
end
nTunC_SF = length(namesTunC_SF_unit);

% Define confidence intervals for bootstrapping
CI95 = .95;
CI68 = 0.68; % shaded band
nPerm = 1e4; % the same as all fxn_xx
nBoot = 1e4; % the same as all fxn_xx
seedPerm = 1; % the same as all fxn_xx
seedBoot = 2; % the same as all fxn_xx
flag_UseRUseRho = 'useR'; % useR or useRho for correlation analysis (for both drawCorr and drawCorrAsym)

% Define locations and their combinations
iLocComb_all = [6,7, 5,3, 1,8]; % combination of locations
iLocGroups_all = {[6,7], [5,3], [1,8]} ; nGroups = length(iLocGroups_all);
% iLocGroups_all = {[6,3]} ; nGroups = length(iLocGroups_all);
iLocSingle_all = 1:5; nLocSingle = length(iLocSingle_all);
iLocSingle_allSets = {[1,6,5,3]};

% Define the number of bins for the NOM analysis
namesMetrics_prob = {'pYES', 'pA'}; nMetrics_prob = length(namesMetrics_prob); namesMetrics_prob_full = {sprintf('Predicted detection prob.\nMeasured detection rate'), sprintf('Predicted consistency prob.\nMeasured resp. consistency')};

% Define settings for plotting RC-derived templates
iDataset_plotRC = 1; %1=Template set; 2=Full set

% Define settings for plotting
iModelA_plot = 1; % just plot the RC-derived
iModelB_plot = 1; % 1=full model; 2=No shared variability; 3=No multiplicative variability
iModelB_plot_all = [1,2]; % generate NOMp-related plots
nNOMparams = length(namesModelBparams{iModelB_plot});
iDataset_plotNOM = 1; % 1=test set, 2=full; search for metrics_allCond

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NOM_mode = 2; % 1= aggregate model, 2 = trial-wise model
flag_subjIsHuman = 1; % 1=human subject, 0=IO data

iModelA_all = iModelA_fit_all;
iModelB_all = iModelB_fit_all;

%%
for iRun=4%[1,3,4]
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
                subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'HL', 'FH', 'HA', 'DT', 'DU', 'RC', 'SR'}; nRows_subj = 3; nCols_subj = 5;
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

    nSubj = length(subjList);

    markers_allSubj = markers_allSubj_full(1:nSubj);


    % Names of the performance metrics being analyzed
    nModelsA = length(iModelA_all);
    nModelsB = length(iModelB_all);
    nLocComb8 = 8;
    iIC_plot = 3; % % Which information criterion to plot (1 = AIC, 2 = AICc, 3 = BIC)

    % Function handles for calculating information criteria (IC)
    getAIC_nLL = @(nLL, nParams) 2*nParams+2*nLL;
    getAICc_nLL = @(nLL, nParams, nData) getAIC_nLL(nLL, nParams) + (2*nParams*(nParams+1))/(nData-nParams-1);
    getBIC_nLL = @(nLL, nParams, nData) nParams*log(nData)+2*nLL;

    % Define folder to save figures (can't put in SX_RC1 because nSubj is not defined yet!)
    nameFolder_Fig_NOM_Trialwise = sprintf('%s/NOM_Trialwise_n%d', nameFolder_Figures, nSubj);
    if isempty(dir(nameFolder_Fig_NOM_Trialwise)), mkdir(nameFolder_Fig_NOM_Trialwise), end

    % Define folder to save outputs
    nameFolder_Outputs_NOM_Trialwise = sprintf('%s/NOM_Trialwise_n%d', nameFolder_Outputs, nSubj);
    if isempty(dir(nameFolder_Outputs_NOM_Trialwise)), mkdir(nameFolder_Outputs_NOM_Trialwise), end

    % Define where to save the compiled data
    nameFolder_Output_SaveCompile = sprintf('%s/n%d_n%d', nameFolder_Outputs_NOM_Trialwise, nSubj, nIterxJob);

    % Print a header to summarize the setting
    fprintf('NOM trial-wise analysis settings:\n')
    fprintf(' - nSubj = %d\n', nSubj)
    fprintf(' - nIter = %d x %d\n', nIter, nJob)
    fprintf(' - Models A: %s\n', strjoin(string(iModelA_all), ', '));
    fprintf(' - Models B: %s\n', strjoin(string(iModelB_all), ', '));
    fprintf(' - Combined locations: %s\n', strjoin(string(iLocComb_all), ', '));
    fprintf(' - ORI tuning function: %s\n', namesFamily_all{iFamily_ORI});
    fprintf(' - SF tuning function: %s\n', namesFamily_all{iFamily_SF});
    fprintf(' - Number of Bins: %d\n', nBins);
    fprintf(' - Information Criterion to plot: %s\n\n', namesIC{iIC_plot});

    %% Compile/load data
    nameFile_compiledData = sprintf('%s.mat', nameFolder_Output_SaveCompile);
    if ~exist(nameFile_compiledData, 'file')

        % Preallocate arrays for storing results across all conditions
        metrics_allCond = nan(nModelsA, nLocComb8, nSubj, nIter * nJob, nDatasets, nMetrics); % nDatasets: 1=Full; 2=Test; metrics come from data_metrics_allIter in OOD_NOM_Trialwise_compDV_A12.m
        template_tmpl_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nORI, nSF);
        template_full_allCond = template_tmpl_allCond;
        DV_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, nBins);
        nTrials_allCond = DV_allCond;
        metric_data_allCond = nan(nModelsA, nModelsB, nLocComb8, nMetrics_prob, nSubj, nIter * nJob, nBins); % "nMetrics" is predefined in this script
        metric_pred_allCond = metric_data_allCond;
        nLL_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob);
        params_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, 3); % 3=Pre-allocate for max params
        R2_NOM_allCond = nan(nModelsA, nModelsB, nLocComb8, nMetrics_prob, nSubj, nIter * nJob);
        R2_w_NOM_allCond = R2_NOM_allCond;

        % Fitting tuning curves
        % nDatasets: 1=Full; 2=Template or test
        sep_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, nDatasets);
        margORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, nDatasets, nORI);
        margPred_ORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, nDatasets, nORI);
        margParams_ORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, nDatasets, length(namesParams_all{iFamily_ORI}));
        margR2_ORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, nDatasets);
        margTunC_ORI_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, nDatasets, length(namesTunC_unit_perF{iFamily_ORI, 2}));

        margSF_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, nDatasets, nSF);
        margPred_SF_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, nDatasets, nSF);
        margParams_SF_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, nDatasets, length(namesParams_all{iFamily_SF}));
        margR2_SF_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, nDatasets);
        margTunC_SF_allCond = nan(nModelsA, nModelsB, nLocComb8, nSubj, nIter * nJob, nDatasets, nTunC_SF);

        % The main analysis loop for all conditions for each subject (Models x Loc x Metrics)
        fprintf('\n ===== START COMPILING ===== \n\n')
        for iJob = 1:nJob
            idxIter = (iJob-1)*nIter+1: iJob*nIter;
            fprintf('\n\n === Job %d/%d === ', iJob, nJob)

            % Loop through each model (1=RC-derived template; 2=ideal template)
            for iModelA = iModelA_all %1:nModelsA
                fprintf('\n\n == Model A%d == ', iModelA)

                % Loop through each model variant ()
                for iModelB = iModelB_all %1:nModelsB
                    fprintf('\n Model B%d: ', iModelB)
                    nParamsB = length(namesModelBparams{iModelB});

                    % Loop through each SINGLE location
                    for iLocSingle = iLocSingle_all
                        fprintf(' L%d ', iLocSingle)

                        % Loop through each subject
                        for isubj = 1:nSubj

                            subjName = subjList{isubj};

                            if flag_subjIsHuman
                                % nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, nameIO);
                                nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocSingle); % MUST match OOD_NOM_Trialwise_compDV_A12.m
                            else
                                % nameFolder_NOM_save = sprintf('%s/ORI%dSF%d/%s/L%d', nameFolder_NOM0, nORI, nSF, subjName, iLocComb);
                                nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, subjName); % MUST match OOD_NOM_Trialwise_compDV_A12.m
                            end

                            % Load DVs / templates from stage A (compDV)
                            nameFile_compDV = sprintf('%s/n%d_J%d_A%d_compDV.mat', nameFolder_NOM_save, nIter, iJob, iModelA);
                            if isempty(dir(nameFile_compDV))
                                fprintf(' "%sL%dA%d" ', subjName, iLocSingle, iModelA)
                                continue
                            end
                            load(nameFile_compDV,  '*_allIter');

                            % Load predictions
                            nameFile_fitNOM = sprintf('%s/n%d_J%d_A%dB%d.mat', nameFolder_NOM_save, nIter, iJob, iModelA, iModelB);
                            if isempty(dir(nameFile_fitNOM))
                                fprintf(' "%sL%dA%dB%d" ', subjName, iLocSingle, iModelA, iModelB)
                                continue
                            end
                            load(nameFile_fitNOM, '*_allIter');

                            % Keep dependent arrays consistent with the current iJob.
                            if size(data_metrics_allIter, 1) ~= nIter
                                fprintf(' "BadMetricsSize:%sL%dA%d" ', subjName, iLocSingle, iModelA)
                                continue
                            end

                            % Loaded *measured* metrics of the full set and test set.
                            % In OOD_NOM_Trialwise_compDV_A12, rows are [full; tmpl; train; test].
                            metrics_allCond(iModelA, iLocSingle, isubj, idxIter, :, :) = data_metrics_allIter(:, [4,1], :); % keep [test, full] ordering used downstream here

                            % Pre-allocate temporary arrays for this subject and model
                            DV_allIter = nan(nIter, nBins); % DV for each bin
                            nTrials_allIter = DV_allIter; % Trial count for each bin
                            metric_data_allIter = nan(nMetrics_prob, nIter, nBins); % Measured data for each bin
                            metric_pred_allIter = metric_data_allIter; % Predicted data for each bin
                            margTuningC_ORI_allIter = nan(nIter, nDatasets, length(namesTunC_unit_perF{iFamily_ORI, 2}));
                            margTuningC_SF_allIter = nan(nIter, nDatasets, nTunC_SF);

                            % Loop through each iteration (within the job)
                            for iIter = 1:nIter

                                % Extract the binned decision variable.
                                DV_allIter(iIter, :) = pred_metrics_allIter{iIter}.metrics.DV_allBins;

                                % Extract the number of trials per bin
                                nTrials_allIter(iIter, :) = pred_metrics_allIter{iIter}.metrics.nTrials_allBins;

                                % Extract metrics (data and predictions) per bin
                                pred_metrics = pred_metrics_allIter{iIter}.metrics; % Extract once for faster access
                                for iMetric_prob = 1:nMetrics_prob
                                    switch namesMetrics_prob{iMetric_prob}
                                        case 'pYES'
                                            metric_data_allIter(iMetric_prob, iIter, :) = pred_metrics.pYES_data_allBins;
                                            metric_pred_allIter(iMetric_prob, iIter, :) = pred_metrics.pYES_pred_allBins;
                                        case 'pA'
                                            metric_data_allIter(iMetric_prob, iIter, :) = pred_metrics.pA_data_allBins;
                                            metric_pred_allIter(iMetric_prob, iIter, :) = pred_metrics.pA_pred_allBins;
                                        case 'pC'
                                            metric_data_allIter(iMetric_prob, iIter, :) = pred_metrics.pC_data_allBins;
                                            metric_pred_allIter(iMetric_prob, iIter, :) = pred_metrics.pC_pred_allBins;
                                    end
                                end % iMetric
                                % Calculate information critertion based on nLL
                                nData = sum(nTrials_allIter(iIter, :));

                                % Compute tuning characteristics based on the fitted parameters
                                for iDataset=1:2 % needs to match OOD_NOM_Trialwise_compDV_A12 (template/full tuning outputs)
                                    % ORI
                                    iFeature= 1;
                                    %======================%
                                    margTuningC_ORI_allIter(iIter, iDataset, :) = fxn_getTuningC( ...
                                        axis_tuning{iFeature}, iFeature, iFamily_ORI, squeeze(margPred_ORI_allIter(iIter, iDataset, :)), squeeze(margParams_ORI_allIter(iIter, iDataset, :)));
                                    % %======================%

                                    % SF
                                    iFeature=2;
                                    %======================%
                                    tuningC_SF = fxn_getTuningC(axis_tuning{iFeature}, iFeature, iFamily_SF, squeeze(margPred_SF_allIter(iIter, iDataset, :)), squeeze(margParams_SF_allIter(iIter, iDataset, :)));
                                    margTuningC_SF_allIter(iIter, iDataset, :) = tuningC_SF;
                                    %======================%
                                end
                            end % end of iIter

                            % Calculate R2 for NOM prediction (for pYES and pA)
                            R2_NOM_allIter = nan(nMetrics_prob, nIter);
                            R2_w_NOM_allIter = R2_NOM_allIter;
                            for iMetric_prob = 1:nMetrics_prob
                                for iIter = 1:nIter
                                    nData = squeeze(nTrials_allIter(iIter, :));
                                    metric_data = squeeze(metric_data_allIter(iMetric_prob, iIter, :)); % ground truth
                                    metric_pred = squeeze(metric_pred_allIter(iMetric_prob, iIter, :)); % prediction
                                    % Remove NaNs if any
                                    valid = ~(isnan(metric_data) | isnan(metric_pred));
                                    metric_data = metric_data(valid);
                                    metric_pred = metric_pred(valid);

                                    % With weighting
                                    metric_ave = sum(nData .* metric_data) / sum(nData);
                                    SSres_w = sum(nData .* (metric_data - metric_pred).^2);
                                    SStot_w = sum(nData .* (metric_data - metric_ave).^2);
                                    if SStot_w == 0, R2_w = NaN; else, R2_w = 1 - SSres_w/SStot_w; end

                                    % No weighting
                                    SSres = sum((metric_data - metric_pred).^2);
                                    SStot = sum((metric_data - mean(metric_data)).^2);
                                    if SStot == 0, R2 = NaN; else, R2 = 1 - SSres/SStot; end

                                    R2_NOM_allIter(iMetric_prob, iIter) = R2;
                                    R2_w_NOM_allIter(iMetric_prob, iIter) = R2_w;
                                end
                            end

                            % Store results
                            % Templates
                            [template_tmpl_med, ~, ~, template_tmpl_sem] = getCI(template_tmpl_allIter, 1, 1);
                            [template_full_med, ~, ~, template_full_sem] = getCI(template_full_allIter, 1, 1);
                            template_tmpl_allCond(iModelA, iModelB, iLocSingle, isubj, :, :) = template_tmpl_med;
                            template_full_allCond(iModelA, iModelB, iLocSingle, isubj, :, :) = template_full_med;

                            % Tuning functions: marg, predictions, estimated parameters and separability
                            sep_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :) = sep_allIter;

                            margORI_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margORI_allIter; % directly loaded
                            margR2_ORI_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :) = margR2_ORI_allIter; % directly loaded
                            margPred_ORI_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margPred_ORI_allIter; % directly loaded
                            margParams_ORI_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margParams_ORI_allIter; % directly loaded
                            margTunC_ORI_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margTuningC_ORI_allIter; % derived above

                            margSF_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margSF_allIter; % directly loaded
                            margR2_SF_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :) = margR2_SF_allIter; % directly loaded
                            margPred_SF_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margPred_SF_allIter; % directly loaded
                            margParams_SF_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margParams_SF_allIter; % directly loaded
                            margTunC_SF_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :, :) = margTuningC_SF_allIter; % derived above

                            % NOM: IVs, predictions, nLL and parameters
                            DV_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :) = DV_allIter;
                            nTrials_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, :) = nTrials_allIter;
                            metric_data_allCond(iModelA, iModelB, iLocSingle, :, isubj, idxIter, :) = metric_data_allIter;
                            metric_pred_allCond(iModelA, iModelB, iLocSingle, :, isubj, idxIter, :) = metric_pred_allIter;
                            nLL_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter) = nLL_test_allIter;
                            % IC_nLL_allCond is compiled above
                            params_allCond(iModelA, iModelB, iLocSingle, isubj, idxIter, 1:nParamsB) = params_est_allIter;
                            R2_NOM_allCond(iModelA, iModelB, iLocSingle, :, isubj, idxIter) = R2_NOM_allIter;
                            R2_w_NOM_allCond(iModelA, iModelB, iLocSingle, :, isubj, idxIter) = R2_w_NOM_allIter;

                            clear *_allIter

                        end % isubj
                    end % iLocComb
                end % iModelB
            end % iModelA
        end % iJob
        fprintf('\n\n ==== NOM outputs compiled ====\n\n')

        % Build combined locations by averaging single-location outputs
        meanOverLoc = @(X, locDim, locOld) mean(X, locDim, 'omitnan');

        for iMap = 6:8 % 6=HM, 7=VM, 8=Peri
            switch iMap
                case 6, locOld = [2 4]; % 6 (HM) = 2 (left) and 4 (right)
                case 7, locOld = [3 5]; % 7 (VM) = 3 (upper) and 5 (lower)
                case 8, locOld = 2:5; % 8 (Perifovea) = 2 to 5
            end

            % locDim = 2
            metrics_allCond(:,iMap,:,:,:,:) = meanOverLoc(metrics_allCond(:,locOld,:,:,:,:), 2, locOld);

            % locDim = 3
            template_tmpl_allCond(:,:,iMap,:,:,:) = meanOverLoc(template_tmpl_allCond(:,:,locOld,:,:,:), 3, locOld);
            template_full_allCond(:,:,iMap,:,:,:) = meanOverLoc(template_full_allCond(:,:,locOld,:,:,:), 3, locOld);

            DV_allCond(:,:,iMap,:,:,:) = meanOverLoc(DV_allCond(:,:,locOld,:,:,:), 3, locOld);
            nTrials_allCond(:,:,iMap,:,:,:) = meanOverLoc(nTrials_allCond(:,:,locOld,:,:,:), 3, locOld);

            metric_data_allCond(:,:,iMap,:,:,:,:) = meanOverLoc(metric_data_allCond(:,:,locOld,:,:,:,:), 3, locOld);
            metric_pred_allCond(:,:,iMap,:,:,:,:) = meanOverLoc(metric_pred_allCond(:,:,locOld,:,:,:,:), 3, locOld);

            nLL_allCond(:,:,iMap,:,:) = meanOverLoc(nLL_allCond(:,:,locOld,:,:), 3, locOld);
            params_allCond(:,:,iMap,:,:,:) = meanOverLoc(params_allCond(:,:,locOld,:,:,:), 3, locOld);
            % IC_nLL_allCond(:,:,iMap,:,:,:) = meanOverLoc(IC_nLL_allCond(:,:,locOld,:,:,:), 3, locOld);
            R2_NOM_allCond(:,:,iMap,:,:,:) = meanOverLoc(R2_NOM_allCond(:,:,locOld,:,:,:), 3, locOld);
            R2_w_NOM_allCond(:,:,iMap,:,:,:) = meanOverLoc(R2_w_NOM_allCond(:,:,locOld,:,:,:), 3, locOld);

            sep_allCond(:,:,iMap,:,:,:) = meanOverLoc(sep_allCond(:,:,locOld,:,:,:), 3, locOld);

            margORI_allCond(:,:,iMap,:,:,:,:) = meanOverLoc(margORI_allCond(:,:,locOld,:,:,:,:), 3, locOld);
            margPred_ORI_allCond(:,:,iMap,:,:,:,:) = meanOverLoc(margPred_ORI_allCond(:,:,locOld,:,:,:,:), 3, locOld);
            margParams_ORI_allCond(:,:,iMap,:,:,:,:) = meanOverLoc(margParams_ORI_allCond(:,:,locOld,:,:,:,:), 3, locOld);
            margR2_ORI_allCond(:,:,iMap,:,:,:) = meanOverLoc(margR2_ORI_allCond(:,:,locOld,:,:,:), 3, locOld);
            margTunC_ORI_allCond(:,:,iMap,:,:,:,:) = meanOverLoc(margTunC_ORI_allCond(:,:,locOld,:,:,:,:), 3, locOld);

            margSF_allCond(:,:,iMap,:,:,:,:) = meanOverLoc(margSF_allCond(:,:,locOld,:,:,:,:), 3, locOld);
            margPred_SF_allCond(:,:,iMap,:,:,:,:) = meanOverLoc(margPred_SF_allCond(:,:,locOld,:,:,:,:), 3, locOld);
            margParams_SF_allCond(:,:,iMap,:,:,:,:) = meanOverLoc(margParams_SF_allCond(:,:,locOld,:,:,:,:), 3, locOld);
            margR2_SF_allCond(:,:,iMap,:,:,:) = meanOverLoc(margR2_SF_allCond(:,:,locOld,:,:,:), 3, locOld);
            margTunC_SF_allCond(:,:,iMap,:,:,:,:) = meanOverLoc(margTunC_SF_allCond(:,:,locOld,:,:,:,:), 3, locOld);
        end

        fprintf('\n\n ==== Combined locations (6,7,8) created by averaging single locations ==== \n\n');

        %%% Save the organized data for all subjects
        save(nameFolder_Output_SaveCompile, '*_allCond')
        clear *_allCond
        fprintf('\n\n ==== *_allCond saved and cleared ==== \n\n');

        %% Compile behavioral data (revise this part later, as the behav should also come from resampled data (metric_data_allCond), not from the raw data)
        CS_allSubj = nan(nSubj, nLocComb8); % 8 is max number of combined locs,see namesLocComb
        dprime_allSubj = CS_allSubj;
        criterion_allSubj = CS_allSubj;
        RT_allSubj = CS_allSubj;
        pC_allSubj = CS_allSubj;
        pA_allSubj = CS_allSubj;

        for isubj = 1:nSubj
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
        save(nameFolder_Output_SaveCompile, '*_allSubj', '-append')

        fprintf('\n ==== Behav data compiled ==== \n\n')
    else
        fprintf('\n%s: Loading compiled data: %s ...', string(datetime('now')), nameFile_compiledData)
        load(nameFolder_Output_SaveCompile)
        fprintf('DONE \n\n')

    end

    %% 1. Behavioral metrics: paired locations
    clc; fprintf('\n%s: Fig 1/23: Behavioral metrics (paired locations) STARTED.\n', string(datetime('now')))

    % Load data
    load(nameFolder_Output_SaveCompile, 'metrics_allCond')

    % Define folder for saving figures
    nameFolder_Fig_behav = sprintf('%s/Behav_pairedLoc', nameFolder_Fig_NOM_Trialwise);
    if isempty(dir(nameFolder_Fig_behav)), mkdir(nameFolder_Fig_behav), end

    flag_plotIDVD = 1;
    flag_plotDiff = 1;
    namesMetrics_behav = {'CS', 'pA', 'dprime', 'criterion', 'pC', 'RT'}; nMetrics_behav = length(namesMetrics_behav);
    namesMetrics_behav_long = {'Contrast sensitivity', 'Consistency rate', 'Dprime', 'Criterion', 'Accuracy', 'Resp. time'};
    [d70,~] = SX_sim06_SDT(.7, .3);

    switch flag_plotPurpose
        case 'paper'
    sz_wd_perBar = 180;
        case 'slide'
            sz_wd_perBar = 100;
    end
    nBars = 2;
    sz_fig = [nBars*sz_wd_perBar, 300+nchoosek(nBars,2)*50];

    for iGroup = 1:nGroups
        iLocPair_all = iLocGroups_all{iGroup};
        for iMetric_prob = 1:nMetrics_behav
            switch iMetric_prob
                case 1, iMetric_vec = 10; x_ticks = linspace(2, 3.6, 5); flag_plotIDVD = 1; ref=nan; % CS
                case 2, iMetric_vec = 6; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 1; ref=nan; % pA
                case 3, iMetric_vec = 1; x_ticks = linspace(0, 1.6, 5); flag_plotIDVD = 0; ref=d70; % dprime
                case 4, iMetric_vec = 2; x_ticks = linspace(-1, 1, 5); flag_plotIDVD = 0; ref=0; % SDT criterion
                case 5, iMetric_vec = 3; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 0; ref=.7; % pC
                case 6, iMetric_vec = 11; x_ticks = linspace(0, .2, 5); flag_plotIDVD = 0; ref=nan; % RT
            end
            x_ticks = round(x_ticks, 2);

            data_allIter_allSubj = squeeze(metrics_allCond(iModelA_plot, iLocPair_all, :, :, iDataset_plotNOM, iMetric_vec));

            str_title = sprintf('n=%d %s nIter=%d L%d%d', nSubj, namesMetrics_behav{iMetric_prob}, nIterxJob, iLocPair_all);
            sz_text = 22;
            wd = 3;
            %------------------------------%
            fxn_drawBars(data_allIter_allSubj, ref, colors_comb(iLocPair_all, :), namesLocComb(iLocPair_all), x_ticks, x_ticks, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj, sz_text, wd, flag_plotPurpose);
            %------------------------------%
            fig = gcf;
            if strcmp(flag_plotPurpose, 'paper')
            ylabel(namesMetrics_behav_long{iMetric_prob})
            end
            saveAndCloseFigure(fig, sprintf('%s/n%d_L%d%d_%s.png', nameFolder_Fig_behav, nSubj, iLocPair_all, namesMetrics_behav{iMetric_prob}))
        end % iMetric
    end % iGroup
    clear metrics_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 2. Behavioral metrics: single locations
    clc; fprintf('\n%s: Fig 2/23: Behavioral metrics (single locations) STARTED.\n', string(datetime('now')))
    % Load data
    load(nameFolder_Output_SaveCompile, 'metrics_allCond')

    % Define folder for saving figures
    nameFolder_Fig_behav = sprintf('%s/Behav_singleLoc', nameFolder_Fig_NOM_Trialwise);
    if isempty(dir(nameFolder_Fig_behav)), mkdir(nameFolder_Fig_behav), end

    flag_plotIDVD = 1;
    flag_plotDiff = 1;
    namesMetrics_behav = {'CS', 'pA', 'dprime', 'criterion', 'pC', 'RT'}; nMetrics_behav = length(namesMetrics_behav);
    namesMetrics_behav_long = {'Contrast sensitivity', 'Consistency rate', 'Dprime', 'Criterion', 'Accuracy', 'Resp. time'};
    [d70,~] = SX_sim06_SDT(.7, .3);

    sz_wd_perBar = 150;

    for iSet = 1:numel(iLocSingle_allSets)

        iLocSingle_perSet = iLocSingle_allSets{iSet};

        nBars = numel(iLocSingle_perSet);
        sz_fig = [nBars*sz_wd_perBar, 300+nchoosek(nBars,2)*50];

        for iMetric_prob = 1:nMetrics_behav
            switch iMetric_prob
                case 1, iMetric_vec = 10; x_ticks = linspace(1.6, 3.6, 5); flag_plotIDVD = 1; ref=0; % CS
                case 2, iMetric_vec = 6; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 1; ref=0; % pA
                case 3, iMetric_vec = 1; x_ticks = linspace(0, 1.6, 5); flag_plotIDVD = 0; ref=d70; % dprime
                case 4, iMetric_vec = 2; x_ticks = linspace(-.5, .5, 5); flag_plotIDVD = 0; ref=0; % criterion
                case 5, iMetric_vec = 3; x_ticks = linspace(.6, .8, 5); flag_plotIDVD = 0; ref=.7; % pC
                case 6, iMetric_vec = 11; x_ticks = linspace(0, .2, 5); flag_plotIDVD = 0; ref=0; % RT
            end
            x_ticks = round(x_ticks, 2);

            data_allIter_allSubj = squeeze(metrics_allCond(iModelA_plot, iLocSingle_perSet, :, :, iDataset_plotNOM, iMetric_vec));

            str_title = sprintf('n=%d %s nIter=%d L%s', nSubj, namesMetrics_behav{iMetric_prob}, nIterxJob, strjoin(string(iLocSingle_perSet), ''));
            sz_text = 22;
            wd = 3;
            %------------------------------%
            fxn_drawBars(data_allIter_allSubj, ref, colors_comb(iLocSingle_perSet, :), namesLocComb(iLocSingle_perSet), x_ticks, x_ticks, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj, sz_text, wd, flag_plotPurpose);
            % ------------------------------%
            fig = gcf;
            ylabel(namesMetrics_behav_long{iMetric_prob})

            saveAndCloseFigure(fig, sprintf('%s/n%d_L%s_%s.png', nameFolder_Fig_behav, nSubj, strjoin(string(iLocSingle_perSet), ''), namesMetrics_behav{iMetric_prob}))
        end % iMetric
    end % iSet
    clear metrics_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 3. Templates for IDVD and group averages
    clc; fprintf('\n%s: Fig 3/23: Templates (group averages) STARTED.\n', string(datetime('now')))
    % Load data
    load(nameFolder_Output_SaveCompile, 'template_*_allCond')

    pixMin = 0;
    pixMax = 0;
    % Define folder for saving figures
    nameFolder_Fig_NOM_Template = sprintf('%s/Template', nameFolder_Fig_NOM_Trialwise);
    if isempty(dir(nameFolder_Fig_NOM_Template)), mkdir(nameFolder_Fig_NOM_Template), end

    for iDataset = 1%1:nDatasets
        switch iDataset % needs to match OOD_NOM_Trialwise_compDV_A12 (1=template, 2=full)
            case 1 % Template set
                template_med_allSubj = squeeze(template_tmpl_allCond(iModelA_plot, iModelB_plot, :, :, :, :)); % ... nLoc x nSubj x nORI x nSF
            case 2 % Full set
                template_med_allSubj = squeeze(template_full_allCond(iModelA_plot, iModelB_plot, :, :, :, :)); % ... nLoc x nSubj x nORI x nSF
        end
        str_dataset = namesDataset{iDataset};

        %%%% Group ave %%%%%
        for iLocComb = iLocComb_all
            fig = figure('Position', [0 0 1e3 1e3]); hold on
            e2D_ave = getCI(template_med_allSubj(iLocComb, :, :, :), 2, 2);
            e2D_ave = e2D_ave';
            fprintf('\n%s L%d: Min = %.3f, Max=%.3f\n', str_dataset, iLocComb, min(e2D_ave(:)), max(e2D_ave(:)))

            if iLocComb==1
                cLim = [-.03, .26];
            else
                cLim = [-1, 8]*1e-3;
                pixMin = min(pixMin, min(e2D_ave(:)));
                pixMax = max(pixMax, max(e2D_ave(:)));
            end
            %-----------------%
            RCplot_2Dkernel(e2D_ave, nan)
            %-----------------%
            title(sprintf('n=%d nIter=%d L%d %s (%s)', nSubj, nIterxJob, iLocComb, namesLocComb{iLocComb}, str_dataset))
            saveAndCloseFigure(fig, sprintf('%s/n%d_group_L%d_A%d_%s.png', nameFolder_Fig_NOM_Template, nSubj, iLocComb, iModelA_plot, str_dataset))
        end % iiLoc

        % %% Idvd data in one figure, per loc %%%%%
        for iiLoc = 1:nLocComb8

            fig = figure('Position', [0, 0, 2e3, 1.8e3]);
            for isubj = 1:nSubj
                e2D = squeeze(template_med_allSubj(iiLoc, isubj, :,:))';
                subplot(nRows_subj, nCols_subj, isubj), hold on
                %-----------------%
                RCplot_2Dkernel(e2D, nan)
                %-----------------%
                % function RCplot_2Dkernel(data2D, caxisLim, outline_pos, outline_neg)
                title(subjList{isubj})

            end % isubj
            set(findall(gcf, '-property', 'FontSize'), 'FontSize',12)
            sgtitle(sprintf('L%d %s [A%dB%d] (%s)', iiLoc, namesLocComb{iiLoc}, iModelA_plot, iModelB_plot, str_dataset), 'FontSize',20)

            % save
            saveAndCloseFigure(fig, sprintf('%s/n%d_L%d_A%dB%d_%s.png', nameFolder_Fig_NOM_Template, nSubj, iiLoc, iModelA_plot, iModelB_plot, str_dataset))

        end % iiLoc
    end % iDataset

    % Print min and max pixel value
    fprintf('\n Summary: Min = %.3f, Max=%.3f\n', pixMin, pixMax)

    clear template_*_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 4. Separability
    clc; fprintf('\n%s: Fig 4/23: Separability analysis STARTED.\n', string(datetime('now')))
    str_dataset = namesDataset{iDataset_plotRC};

    % Load compiled separability: sep_allCond
    load(nameFolder_Output_SaveCompile, 'sep_allCond')

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

        % Reformat to match fxn_drawBars:
        % data_allIter_allSubj: [nIter x nSubj x nCond], where nCond = nLoc
        data_allIter_allSubj = permute(sep_allIter_allSubj, [3 2 1]); % [nIter x nSubj x nLoc]

        % Colors for each bar/condition (each location)
        colors = colors_comb(iLocSingle_perSet, :); % [nLoc x 3]

        % Condition labels (not displayed on x-axis inside the function, but used in title strings)
        x_ticks = cell(1, nLoc);
        for iiLoc = 1:nLoc
            x_ticks{iiLoc} = sprintf('L%d', iLocSingle_perSet(iiLoc));
        end

        % Plot settings for this figure
        ref = nan; % reference line used in your old plot
        y_ticks = [0.6 0.7 0.8 0.9 1.0];
        y_ticklabels = y_ticks;

        flag_plotIDVD = 0; % show subject-level lines (recommended)
        flag_plotDiff = 0; % only meaningful when nCond==2 (turn on if you want)
        sz_wd_perBar = 150;
        nBars =nLoc;
        sz_fig = [nBars*sz_wd_perBar, 300+nchoosek(nBars,2)*50];

        str_loc = strjoin(string(iLocSingle_perSet), '');
        str_title = sprintf('n%d nIter=%d [A%dB%d] L%s %s', nSubj, nIterxJob, iModelA_plot, iModelB_plot, str_loc, str_dataset);
        sz_text = 22;
        wd = 3;

        % -------------------%
        fxn_drawBars(data_allIter_allSubj, ref, colors, x_ticks, y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj, sz_text, wd, flag_plotPurpose);
        % -------------------%
        fig = gcf;
        ylabel('Separability (Pearson''s r)')

        % Adjust distance between components
        ax = gca;
        drawnow; % Ensure text/tick extents are up to date before reading TightInset
        ti = ax.TightInset; % [left bottom right top] padding needed for labels/ticks
        ax.XLabel.Units = 'normalized';
        ax.YLabel.Units = 'normalized';

        ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset; % move label down
        ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset; % move label left

        newPos = [ ...
            ti(1) + sz_axOffset, ...
            ti(2) + sz_axOffset, ...
            1 - ti(1) - ti(3) - 2*sz_axOffset, ...
            1 - ti(2) - ti(4) - 2*sz_axOffset];
        if all(isfinite(newPos)) && newPos(3) >= 0.55 && newPos(4) >= 0.55
            ax.Position = newPos;
        end

        % ---- save ----
        outName = sprintf('%s/n%d_L%s_A%d_%s.png', nameFolder_Fig_Sep, nSubj, str_loc, iModelA_plot, str_dataset);
        saveAndCloseFigure(fig, outName);

    end % iSet

    clear sep_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 5. Tuning functions: group averages
    clc; fprintf('\n%s: Fig 5/23: Tuning functions (group averages) STARTED.\n', string(datetime('now')))
    % Load data
    load(nameFolder_Output_SaveCompile, 'marg*_allCond', 'margPred*_allCond', 'margR2*_allCond')

    % Define folder for saving figures
    nameFolder_Fig_NOM_Tuning = sprintf('%s/TuningFxns_group', nameFolder_Fig_NOM_Trialwise);
    if isempty(dir(nameFolder_Fig_NOM_Tuning)), mkdir(nameFolder_Fig_NOM_Tuning), end

    wd_border = 3.5; % default 5
    sz_ticks = 30;% default 35
    sz_label = sz_ticks;

    iplots = reshape((1:nGroups*nFeatures)', nGroups, nFeatures)';
    for iDataset=1%1:nDatasets
        switch iDataset
            case 1 % Template set
                markerStyle = 'o';
                lineStyle = '-';
            case 2 % Full set
                markerStyle = 's';
                lineStyle = '--';
        end

        for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
            iLocPair_all = iLocGroups_all{iGroup};

            for iFeature = 1:nFeatures

                xaxis = axis_tuning{iFeature};

                switch flag_plotPurpose
                    case 'slide'
                fig = figure('Position', [0 0 8e2 6e2]); % default 8e2
                    case 'paper'
                        fig = figure('Position', [0 0 1.1e3 6e2]); % default 8e2
                end
                hold on
                for iLoc = iLocPair_all

                    fprintf('\n   - %s L%s %s', namesDataset{iDataset}, strjoin(string(iLocPair_all), ''), namesFeature{iFeature})

                    % Set color
                    color_comb = colors_comb(iLoc, :);

                    % Set y-ticks
                    if iFeature==1 % ORI tuning fxn
                        marg_allCond = margORI_allCond;
                        margPred_allCond = margPred_ORI_allCond;
                        margR2_allCond = margR2_ORI_allCond;
                        if find(iLocPair_all==1), yticks_ = linspace(-2, 8, 5)*1e-3; % fov vs. peri, higher ub
                        else, yticks_ = linspace(-1, 5, 5)*1e-3; %yticks_ = [-.03, linspace(0, .12, 4)];
                        end
                    else % SF tuning fxn
                        marg_allCond = margSF_allCond;
                        margPred_allCond = margPred_SF_allCond;
                        margR2_allCond = margR2_SF_allCond;
                        if find(iLocPair_all==1), yticks_ = [-.5, linspace(0, 2, 4)]*1e-3;
                        else, yticks_ = linspace(-.5, 2, 5)*1e-3; %[-.02, 0, .02, .04, .06];
                        end
                    end

                    ymax = max(yticks_);
                    ymin = min(yticks_);

                    % Obtain idvd data (median of across-validated iterations)
                    marg_med = getCI(marg_allCond(iModelA_plot, iModelB_plot, iLoc, :, :, iDataset, :), 1, 5);
                    margPred_med = getCI(margPred_allCond(iModelA_plot, iModelB_plot, iLoc, :, :, iDataset, :), 1, 5);
                    R2_med = getCI(margR2_allCond(iModelA_plot, iModelB_plot, iLoc, :, :, iDataset), 1, 5);

                    % Bootstraping to derive median and CI of group averages
                    rng(seedBoot, 'twister');
                    indRand_allBoot = randi(nSubj, [nBoot, nSubj], 'uint16');

                    marg_groupAve_allBoot = nan(nBoot, size(marg_med, 2));
                    margPred_groupAve_allBoot = nan(nBoot, size(margPred_med, 2));
                    R2_groupAve_allBoot = nan(nBoot, 1);

                    parfor iBoot = 1:nBoot
                        indRandBoot = double(indRand_allBoot(iBoot,:));
                        marg_groupAve_allBoot(iBoot, :) = mean(marg_med(indRandBoot, :));
                        margPred_groupAve_allBoot(iBoot, :) = mean(margPred_med(indRandBoot, :));
                        R2_groupAve_allBoot(iBoot) = mean(R2_med(indRandBoot, :));
                    end

                    % Obtain median and 95% CI of boostrapped values
                    [marg_med, marg_lb, marg_ub, marg_sem_neg, marg_sem_pos] = getCI(marg_groupAve_allBoot, 1, 1, CI68);
                    [margPred_med, margPred_lb, margPred_ub] = getCI(margPred_groupAve_allBoot, 1, 1, CI68);
                    [R2_med, R2_lb, R2_ub] = getCI(R2_groupAve_allBoot, 1, 1, CI95);

                    % Prediction (lines + bands)
                    patch([xaxis, flip(xaxis)], [margPred_lb, flip(margPred_ub)], color_comb, 'FaceAlpha', .3, 'linestyle', 'none')
                    plot(xaxis, margPred_med, '-', 'color', color_comb, 'linewidth', wd_border*2)

                    % Data (dots + errorbars)
                    errorbar(xaxis, marg_med, marg_sem_neg, marg_sem_pos, markerStyle, 'Color', color_comb, 'CapSize',0)
                    plot(xaxis, marg_med, 's', 'color', color_comb, 'MarkerFaceColor', 'w', 'MarkerSize', 15, 'HandleVisibility','off')

                    % Draw reference lines
                    yline(0, '--', 'handlevisibility', 'off', 'linewidth', wd_border*2, 'color', [.7, .7, .7]);
                    xline(iFeature-1, '--', 'handlevisibility', 'off', 'linewidth', wd_border*2, 'color', [.7, .7, .7]);

                    % Set ticks, labels and limits (for EACH loc, to print R2 at the right loc)
                    xlim(axisLim{iFeature})
                    ylim([ymin, ymax])
                    yticks(yticks_)

                    xticks(axisTicks_tuning{iFeature})
                    xticklabels(axisTL_tuning{iFeature})
                    xlabel(namesFeature_axis{iFeature}, 'FontSize', sz_label)
                    ylabel(namesFeature_axis_Tuning{iFeature})

                    if iFeature==2, xticklabels(round(axisTL_tuning{iFeature}, 2)), end

                    % Print R2 in normalized axes coordinates so placement is stable across y-ranges
                    x_R2 = 0.98;
                    y_R2 = 0.95 - 0.08 * (find(iLoc == iLocPair_all) - 1);

                    text(x_R2, y_R2, ...
                        sprintf('R^2 = %.2f [%.2f, %.2f]', R2_med, R2_lb, R2_ub), ...
                        'Units', 'normalized', ...
                        'Interpreter', 'tex', ...
                        'FontSize', 35, ...
                        'Color', color_comb, ...
                        'HorizontalAlignment', 'right', ...
                        'VerticalAlignment', 'top', ...
                        'Interpreter','tex');
                end % iLoc

                ax = gca;
                ax.XAxis.FontSize = sz_ticks;
                ax.YAxis.FontSize = sz_ticks;
                ax.LineWidth = wd_border/1.5;

                title(sprintf('n=%d L%d%d [A%dB%d] %s tuning (%s)', nSubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, namesDataset{iDataset}))

                set(findall(gcf, '-property', 'linewidth'), 'linewidth', 2)

                % Save the figure
                saveAndCloseFigure(fig, sprintf('%s/n%d_L%d%d_A%d_%s_%s.png', nameFolder_Fig_NOM_Tuning, nSubj, iLocPair_all, iModelA_plot, namesFeature{iFeature}, namesDataset{iDataset}))
            end % iFeature
        end % iGroup
    end % iDataset = 1:2
    clear marg*_allCond margPred*_allCond margR2*_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 6. Tuning functions: idvd (LARGE FIGURES!!)
    clc; fprintf('\n%s: Fig 6/23: Tuning functions (individual) STARTED.\n', string(datetime('now')))

    % Load data
    load(nameFolder_Output_SaveCompile, 'marg*_allCond', 'margPred*_allCond', 'margParams*_allCond', 'margR2*_allCond')

    % Define folder for saving figures
    nameFolder_Fig_NOM_Tuning = sprintf('%s/TuningFxns_IDVD', nameFolder_Fig_NOM_Trialwise);
    if isempty(dir(nameFolder_Fig_NOM_Tuning)), mkdir(nameFolder_Fig_NOM_Tuning), end

    wd_border = 4; % default 5
    sz_ticks = 30;% default 35

    % Temporary diagnostics: per-iteration SF tuning plots (data + prediction + peak markers).
    flag_temp_plotIterSF = 1;
    temp_nIterPlotMax = 50;
    temp_iterPauseSec = 0.1;
    if flag_temp_plotIterSF
        nameFolder_Fig_NOM_Tuning_tmp = sprintf('%s/TuningFxns_IDVD_tmpIterSF', nameFolder_Fig_NOM_Trialwise);
        if isempty(dir(nameFolder_Fig_NOM_Tuning_tmp)), mkdir(nameFolder_Fig_NOM_Tuning_tmp), end
    end

    for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
        iLocPair_all = iLocGroups_all{iGroup};

        for iFeature = 2%1:nFeatures

            fprintf('\n   - L%s %s...', strjoin(string(iLocPair_all), ''), namesFeature{iFeature})

            xaxis = axis_tuning{iFeature};

            % Sets the tolerance for peak SF estimation in octaves.
            % This is used to calculate the percentage of iterations where the predicted peak SF is within this tolerance of the true peak SF (derived from fitted parameters). The tolerance is set to half the channel spacing in log2(SF), which reflects the resolution of the SF axis in the tuning function.
            if iFeature == 2
                peakSF_tol_oct = mean(diff(xaxis)) / 2; % tolerance=half channel spacing in log2(SF)
            end

            for iDataset = 1%:nDatasets
                switch iDataset
                    case 1 % TmplSet
                        markerStyle = 'o';
                        lineStyle = '-';
                    case 2 % FullSet
                        markerStyle = 's';
                        lineStyle = '--';
                end

                if nSubj==12
                    fig = figure('Position', [0 0 2e3 1.5e3]);
                else
                    fig = figure('Position', [0 0 2e3 1.2e3]);
                end

                for isubj = 1:nSubj
                    subplot(nRows_subj, nCols_subj, isubj), hold on

                    subjName = subjList{isubj};

                    str_TunParams ='';

                    % Draw reference lines
                    yline(0, 'handlevisibility', 'off', 'linewidth', wd_border/2, 'color', [.7, .7, .7]);
                    xline(iFeature-1, 'handlevisibility', 'off', 'linewidth', wd_border/2, 'color', [.7, .7, .7]);

                    for iLoc = iLocPair_all

                        % fprintf('\nL%d...', iLoc)
                        % Set color
                        color_comb = colors_comb(iLoc, :);

                        % Set y-ticks
                        if iFeature==1 % ORI
                            marg_allCond = margORI_allCond;
                            margPred_allCond = margPred_ORI_allCond;
                            margParams_allCond = margParams_ORI_allCond;
                            margR2_allCond = margR2_ORI_allCond;
                            if iLoc==1, yticks_ = linspace(-4,10,5)*1e-3; % fov vs. peri, higher ub
                            else, yticks_ = linspace(-4,6,5)*1e-3;
                            end

                        else % SF
                            marg_allCond = margSF_allCond;
                            margPred_allCond = margPred_SF_allCond;
                            margParams_allCond = margParams_SF_allCond;
                            margR2_allCond = margR2_SF_allCond;
                            if iLoc==1, yticks_ = linspace(-4,4,5)*1e-3;
                            else, yticks_ = linspace(-2,4,5)*1e-3;
                            end
                        end

                        ymax = max(yticks_);
                        ymin = min(yticks_);

                        % Obtain data, prediction, and params
                        [marg_med_perS, ~, ~, marg_sem_neg_perS, marg_sem_pos_perS] = getCI(marg_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :), 1, 5);
                        [margPred_med_perS, margPred_lb_perS, margPred_ub_perS] = getCI(margPred_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :), 1, 5);
                        [margParam_med_perS, margParam_lb_perS, margParam_ub_perS] = getCI(margParams_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :), 1, 5);
                        margR2_med_perS = getCI(margR2_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :), 1, 5);

                        % Per-iteration peakSF validation using already-compiled outputs.
                        if iFeature == 2
                            margPred_iter_perS = squeeze(margPred_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :));
                            margParam_iter_perS = squeeze(margParams_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :));
                            margPred_iter_perS = reshape(margPred_iter_perS, [], numel(xaxis));
                            margParam_iter_perS = reshape(margParam_iter_perS, [], size(margParam_iter_perS, ndims(margParam_iter_perS)));

                            xPeakParam_log2_allIter = log2(margParam_iter_perS(:, 1));
                            [~, idxPeakPred_allIter] = max(margPred_iter_perS, [], 2);
                            xPeakPred_log2_allIter = xaxis(idxPeakPred_allIter)';

                            peakErr_oct_allIter = abs(xPeakPred_log2_allIter - xPeakParam_log2_allIter);
                            peakErr_valid = peakErr_oct_allIter(isfinite(peakErr_oct_allIter));
                            peakPassRate_pct = 100 * mean(peakErr_valid <= peakSF_tol_oct);
                            peakErr_med_oct = median(peakErr_valid, 'omitnan');

                            gain_allIter = margParam_iter_perS(:, 2);
                            width_allIter = margParam_iter_perS(:, 3);
                            nInvalidShape = sum(~isfinite(gain_allIter) | ~isfinite(width_allIter) | gain_allIter <= 0 | width_allIter <= 0);
                        end

                        % Data (dots + errorbars)
                        errorbar(xaxis, marg_med_perS, marg_sem_neg_perS, marg_sem_pos_perS, markerStyle, 'Color', color_comb, 'CapSize',0, 'linewidth', wd_border/2)
                        plot(xaxis, marg_med_perS, markerStyle, 'color', color_comb, 'linewidth', wd_border/2, 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'HandleVisibility','off')

                        % Prediction (lines + bands)
                        patch([xaxis, flip(xaxis)], [margPred_lb_perS', flip(margPred_ub_perS')], color_comb, 'FaceAlpha', .3, 'linestyle', 'none')
                        plot(xaxis, margPred_med_perS, '-', 'color', color_comb, 'linewidth', wd_border)

                        % Print estimated parameters with scale-aware formatting.
                        paramStr = arrayfun(@formatMixedNumber, margParam_med_perS, 'UniformOutput', false);
                        str_TunParams = sprintf('%s\n[L%d] [R2=%.0f%%] %s', str_TunParams, iLoc, margR2_med_perS*100, strjoin(paramStr, ', '));
                        if iFeature == 2
                            % fprintf('\npeakChk: %.0f%%<=%.3f oct, medErr=%.3f, bad(g<=0|w<=0)=%d', peakPassRate_pct, peakSF_tol_oct, peakErr_med_oct, nInvalidShape);
                        end

                        % Draw peak SF
                        if iFeature==2
                            % Peak from the displayed median fitted curve (guaranteed to align visually).
                            [~, idxPeakPred] = max(margPred_med_perS);
                            xPeakPred_log2 = xaxis(idxPeakPred);
                            xline(xPeakPred_log2, 'Color', color_comb, 'LineWidth', wd_border);

                            % Parameter-median peakSF as reference (can differ from peak of median curve).
                            xPeakParam_log2 = log2(margParam_med_perS(1));
                            xline(xPeakParam_log2, '--', 'Color', color_comb, 'LineWidth', wd_border/1.5, 'HandleVisibility', 'off');
                            errorbar(xPeakParam_log2, 0, log2(margParam_med_perS(1))-log2(margParam_lb_perS(1)), log2(margParam_ub_perS(1))-log2(margParam_med_perS(1)), 'horizontal', 'Color', color_comb, 'LineWidth', wd_border)
                        end

                        % Temp: plot per-iteration SF tuning function (one figure per iteration).
                        if iFeature == 2 && flag_temp_plotIterSF
                            marg_iter_perS = squeeze(marg_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :));
                            margPred_iter_perS = squeeze(margPred_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :));
                            margParam_iter_perS = squeeze(margParams_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :));
                            margR2_iter_perS = squeeze(margR2_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :));

                            marg_iter_perS = reshape(marg_iter_perS, [], numel(xaxis));
                            margPred_iter_perS = reshape(margPred_iter_perS, [], numel(xaxis));
                            margParam_iter_perS = reshape(margParam_iter_perS, [], size(margParam_iter_perS, ndims(margParam_iter_perS)));
                            margR2_iter_perS = reshape(margR2_iter_perS, [], 1);

                            nIterAvail = size(margPred_iter_perS, 1);
                            nIterPlot = min(temp_nIterPlotMax, nIterAvail);
                            idxIterPlot = unique(round(linspace(1, nIterAvail, nIterPlot)));

                            for iIter = idxIterPlot
                                thisData = marg_iter_perS(iIter, :);
                                thisPred = margPred_iter_perS(iIter, :);
                                thisParam = margParam_iter_perS(iIter, :);
                                thisR2 = margR2_iter_perS(iIter);

                                xPeakParam_log2_iter = log2(thisParam(1));
                                [~, idxPeakPred_iter] = max(thisPred);
                                xPeakArgmax_log2_iter = xaxis(idxPeakPred_iter);

                                y_min_iter = min([thisData(:); thisPred(:)], [], 'omitnan');
                                y_max_iter = max([thisData(:); thisPred(:)], [], 'omitnan');
                                y_rng_iter = y_max_iter - y_min_iter;
                                if ~isfinite(y_rng_iter) || y_rng_iter == 0, y_rng_iter = 1; end
                                y_peakMarker = y_min_iter - 0.08 * y_rng_iter;

                                fig_tmp = figure('Position', [100 80 1000 700]);
                                ax_tmp = axes(fig_tmp); hold(ax_tmp, 'on')

                                h_data = plot(xaxis, thisData, 'o', 'Color', [0.55 0.55 0.55], 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.1);
                                h_pred = plot(xaxis, thisPred, '-', 'Color', [0.15 0.45 0.85], 'LineWidth', 2);
                                h_peakParam = plot(xPeakParam_log2_iter, y_peakMarker, 'x', 'Color', [0.85 0.2 0.2], 'MarkerSize', 9, 'LineWidth', 1.5);
                                h_peakArg = plot(xPeakArgmax_log2_iter, y_peakMarker, 'o', 'Color', [0.1 0.6 0.1], 'MarkerSize', 7, 'LineWidth', 1.2);
                                xline(xPeakParam_log2_iter, '--', 'Color', [0.85 0.2 0.2], 'LineWidth', 1.2, 'HandleVisibility', 'off');
                                xline(xPeakArgmax_log2_iter, '-', 'Color', [0.1 0.6 0.1], 'LineWidth', 1.2, 'HandleVisibility', 'off');

                                xlim(axisTicks_tuning{iFeature}([1, end]))
                                xticks(axisTicks_tuning{iFeature})
                                xticklabels(axisTL_tuning{iFeature})
                                ylim([y_peakMarker - 0.05 * y_rng_iter, y_max_iter + 0.08 * y_rng_iter])
                                xlabel(namesFeature_axis{iFeature})
                                ylabel('Marg. weights (a.u.)')

                                paramStr_iter = arrayfun(@formatMixedNumber, thisParam, 'UniformOutput', false);
                                title(sprintf('TEMP iter-SF | %s | L%d | %s | iter %d/%d | R2=%.1f%%\nparams=[%s]', ...
                                    subjName, iLoc, namesDataset{iDataset}, iIter, nIterAvail, thisR2*100, strjoin(paramStr_iter, ', ')))
                                legend([h_data, h_pred, h_peakParam, h_peakArg], {'data', 'prediction', 'peakSF param (log2)', 'argmax(pred)'}, 'Location', 'best')

                                drawnow
                                % if temp_iterPauseSec > 0, pause(temp_iterPauseSec), end
                                close(fig_tmp) 

                            end
                        end
                    end % iLoc

                    % y ticks
                    yticks(yticks_)
                    ylim([ymin, ymax])

                    % x ticks
                    xlim(axisTicks_tuning{iFeature}([1, end]))
                    xticks(axisTicks_tuning{iFeature})
                    xticklabels(axisTL_tuning{iFeature})
                    xlabel(namesFeature_axis{iFeature})
                    ylabel('Marg. weights (a.u.)')
                    title(sprintf('%s\n%s\n', subjName, str_TunParams));
                    % legend('Location', 'best')

                end % isubj
                sgtitle(sprintf('n=%d L%d%d [A%d] %s tuning (%s)\nParams (NOT tunC) are printed in the title: [%s]\n', ...
                    nSubj, iLocPair_all, iModelA_plot, namesFeature{iFeature}, namesDataset{iDataset}, ...
                    strjoin(namesParams_all{iFamily_perF(iFeature)}, ', ')))
                set(findall(gcf, '-property', 'fontsize'), 'fontsize', 15)

                % Save the figure
                saveAndCloseFigure(fig, sprintf('%s/n%d_L%d%d_A%d_%s_%s.png', nameFolder_Fig_NOM_Tuning, nSubj, iLocPair_all, iModelA_plot, namesFeature{iFeature}, namesDataset{iDataset}))
            end % iDataset
        end % iFeature
    end % iGroup
    clear marg*_allCond margPred*_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 7. Tuning characteristics
    clc; fprintf('\n%s: Fig 7/23: Tuning characteristics STARTED.\n', string(datetime('now')))

    % Load data
    load(nameFolder_Output_SaveCompile, 'margParams_*_allCond')

    % Define folder for saving figures
    nameFolder_Fig_tunC = sprintf('%s/TuningCs', nameFolder_Fig_NOM_Trialwise);
    if isempty(dir(nameFolder_Fig_tunC)), mkdir(nameFolder_Fig_tunC), end

    flag_plotDist = 0;
    flag_plotIDVD = 1;
    flag_plotDiff = 1;
    paramMode = 2;

    sz_text = 22;
    wd = 2;

    for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
        iLocPair_all = iLocGroups_all{iGroup};

        switch flag_plotPurpose
            case 'paper'
                sz_wd_perBar = 200;
            case 'slide'
                sz_wd_perBar = 250;
        end
        nBars = numel(iLocPair_all);
        sz_fig = [nBars*sz_wd_perBar, 400+nchoosek(nBars,2)*50];

        % Plotting settings
        colors = colors_comb(iLocPair_all, :);
        x_ticks = namesLocComb(iLocPair_all);

        % Loop through each feature
        for iFeature = 1:nFeatures

            % Define the x-axis and parameters for the current feature
            xaxis = axis_tuning{iFeature};
            nfilters = length(xaxis);
            iFamily = iFamily_perF(iFeature);
            str_family = sprintf('F%d %s', iFamily, namesFamily_all{iFamily});
            namesTunC = namesParams_all{iFamily};
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
                                y_ticks_all{2} = linspace(0, 60, 5); % ORI band
                                y_ticks_all{3} = linspace(-.05, .03, 5); % ORI baseline
                            else
                                if find(iLocPair_all==1)
                                    y_ticks_all{1} = linspace(0, .24, 5); % ORI peak amp
                                else
                                    y_ticks_all{1} = linspace(0, .12, 5); % ORI peak amp
                                end

                                y_ticks_all{2} = linspace(0, 60, 5); % ORI band
                                y_ticks_all{3} = linspace(-.08, .08, 5); % ORI baseline
                            end

                        case 10 % ORI tuningC | von Mises
                            if flag_plotIDVD
                                y_ticks_all{1} = linspace(0, 16, 5)*1e-3; % ORI peak amp
                                y_ticks_all{2} = linspace(0,6, 5); % ORI band
                                y_ticks_all{3} = linspace(-.05, .03, 5); % ORI baseline
                            else
                                if find(iLocPair_all==1)
                                    y_ticks_all{1} = linspace(0, .24, 5); % ORI peak amp
                                else
                                    y_ticks_all{1} = linspace(0, .12, 5); % ORI peak amp
                                end

                                y_ticks_all{2} = linspace(0, 6, 5); % ORI band
                                y_ticks_all{3} = linspace(-.08, .08, 5); % ORI baseline
                            end

                        case 2 % SF tuniningC | log parabola
                            if flag_plotIDVD

                                y_ticks_all{1} = linspace(-2, 2, 5); % peak SF
                                if find(iLocPair_all==1)
                                    y_ticks_all{2} = linspace(.01, .2, 5); % SF peak amp
                                else
                                    y_ticks_all{2} = linspace(.0, 4, 5)*1e-3; % SF peak amp
                                end
                                y_ticks_all{3} = linspace(1, 3, 5); % SF bandwidth
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

                        case 14 % SF tuningC | asymmetric Gaussian
                            if flag_plotIDVD

                                y_ticks_all{1} = linspace(0, 2, 5); % peak SF
                                if find(iLocPair_all==1)
                                    y_ticks_all{2} = linspace(.01, .2, 5); % SF peak amp
                                else
                                    y_ticks_all{2} = linspace(.0, .12, 5); % SF peak amp
                                end
                                y_ticks_all{3} = linspace(.6, 5, 5); % SF width (left+right), octaves
                                y_ticks_all{4} = linspace(-.1, .1, 5); % SF baseline
                            else
                                y_ticks_all{1} = linspace(0, 2, 5); % peak SF

                                if find(iLocPair_all==1)
                                    y_ticks_all{2} = linspace(.02, .12, 5); % SF peak amp
                                else
                                    y_ticks_all{2} = linspace(.02, .06, 5); % SF peak amp
                                end
                                y_ticks_all{3} = linspace(.6, 4, 5); % SF width (left+right), octaves
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

                        case 10 % ORI tuningC | von Mises
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

                        case 14 % SF tuningC | asymmetric Gaussian
                            y_ticks_all{1} = linspace(log2(1.4), log2(2.8), 5); % peak SF
                            if find(iLocPair_all==1)
                                y_ticks_all{2} = linspace(.01, .13, 5); % SF peak amp (fovea)
                            else
                                y_ticks_all{2} = linspace(.03, .07, 5); % SF peak amp
                            end
                            y_ticks_all{3} = linspace(1, 1.8, 5); % SF width (left+right), octaves
                            y_ticks_all{4} = linspace(-.04, .04, 5); % SF baseline

                    end
            end

            % Loop through each tuning characteristic
            for iTunC = 1:nTunC_full

                % Compute median for each observer
                % margParams_ORI_allCond: nModelA x nModelB x nLoc_pair x nSubj x nIter x nDataset x nParam
                % tunC_allSubj_allIter: nLoc_pair x nSubj x nIter
                switch iFeature
                    case 1
                        data_allIter_allSubj = squeeze(margParams_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                        data_obs_allSubj = squeeze(margParams_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, 1, 2, iTunC)); % full data, without resampling
                    case 2
                        data_allIter_allSubj = squeeze(margParams_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                        data_obs_allSubj = squeeze(margParams_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, 1, 2, iTunC));
                end

                % Let MATLAB format y-axis ticks automatically unless a custom unit remapping is needed.
                y_ticklabels = nan;
                ref = nan;
                switch iFamily_perF(iFeature)
                    case 1
                        if find(iTunC==3), ref = 0; end

                    case 10
                        if find(iTunC==3), ref = 0; end

                    case 2 % log parabola
                        switch iTunC
                            case 1
                                y_ticklabels = 2.^y_ticks_all{iTunC};
                                y_ticklabels = round(y_ticklabels, 2);
                                ref = log2(2);
                                data_allIter_allSubj = log2(data_allIter_allSubj); % pref SF
                            case 4, ref = 0; % baseline
                        end

                    case 14 % asymmetric Gaussian
                        switch iTunC
                            case 1, y_ticklabels = round(2.^y_ticks_all{iTunC}, 1); ref = log2(2);
                                data_allIter_allSubj = log2(data_allIter_allSubj); % pref SF
                            case 4, ref = 0; % baseline
                        end
                end

                str_title = sprintf('n=%d nIter=%d L%d%d [A%d] [%s] | %s %s | %s', ...
                    nSubj, nIterxJob, iLocPair_all, iModelA_plot, namesDataset{iDataset_plotRC}, namesFeature{iFeature}, namesTunC{iTunC}, str_family);

                switch flag_plotDist
                    case 0
                        %------------------------------%
                        fxn_drawBars(data_allIter_allSubj, ref, colors, x_ticks, y_ticks_all{iTunC}, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj, sz_text, wd, flag_plotPurpose);
                        % ------------------------------%
                        fig = gcf;
                        ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunC{iTunC}))
                    case 1
                        %------------------------------%
                        fxn_drawDist(data_allIter_allSubj, ref, colors, x_ticks, y_ticks_all{iTunC}, y_ticklabels, str_title, sz_fig, nIterxJob, nSubj)
                        %------------------------------%
                        fig = gcf;
                        xlabel(sprintf('%s %s [%s]', namesFeature{iFeature}, namesTunC{iTunC}, str_family))
                        ylabel('Probabillity')
                end

                % Constrain axis layout for Fig 7 to keep long y-labels visible.
                drawnow;
                ax = gca;
                if isgraphics(ax, 'axes')
                    oldUnits = ax.Units;
                    ax.Units = 'normalized';
                    ti = ax.TightInset; % [left bottom right top]

                    left = max(0.26, ti(1) + 0.04);
                    bottom = max(0.14, ti(2) + 0.03);
                    right = max(0.03, ti(3) + 0.02);
                    top = max(0.05, ti(4) + 0.02);

                    w = 1 - left - right;
                    h = 1 - bottom - top;
                    if isfinite(w) && isfinite(h) && w > 0.30 && h > 0.30
                        ax.Position = [left, bottom, w, h];
                    end
                    ax.Units = oldUnits;
                end

                % Save the figure
                saveAndCloseFigure(fig, sprintf('%s/n%d_L%d%d_A%d_%s%d.png', nameFolder_Fig_tunC, nSubj, iLocPair_all, iModelA_plot, namesFeature{iFeature}, iTunC))

            end % end of iTunC
        end % end of iFeature
    end % end of iGroup
    clear margTunC_*_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 8. Corr0: Corr between CS and tunParams (just to check bound-hitting)
    clc; fprintf('\n%s: Fig 8/23: Corr CS vs tuning params STARTED.\n', string(datetime('now')))

    % Load data
    load(nameFolder_Output_SaveCompile, 'margParams*_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'tunParam';
    sz_label = 55;
    sz_labelOffset = .05;
    sz_axOffset = 0.05; % extra breathing room

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

    nameFolder_Outputs_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Outputs_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Outputs_NOM_corr)); mkdir(nameFolder_Outputs_NOM_corr), end

    flag_zeroMean = 0;
    flag_plotIdvdCI = 0;
    flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles

    for iSet = 1:length(iLocSingle_allSets)

        iLocCorr_all = iLocSingle_allSets{iSet};
        fprintf(' - L%s\n', strjoin(string(iLocCorr_all), ''))

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
                    case 1, y_ticks_allTunC_lb = [0, 0, -.1]; y_ticks_allTunC_ub = [.36, 60, .1];% 3 values are ORI peak amplitude, width, baseline
                    case 2, y_ticks_allTunC_lb = [0, 0, 0, -.1]; y_ticks_allTunC_ub = [4, .24, 1, .1]; % 4 values are SF peak, peak amplitude, width, baseline
                end

                % y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);
                y_ticks = [];

                x_ticklabels = x_ticks;
                y_ticklabels = y_ticks;

                nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunCs{iTunC});
                nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

                str_title = sprintf('n=%d, nIter=%d %s vs. %s [L%s]', nSubj, nIterxJob, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

                %----------------------------%
                [PearsonR, SpearmanRho] = fxn_drawCorr(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), flag_UseRUseRho, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, nIterxJob, flag_plotPurpose);
                %----------------------------%
                fig = gcf;

                % Save CI for automatic CI range calcuation in CorrAsym
                save(sprintf('%s/n%d_L%s_%s.mat', nameFolder_Outputs_NOM_corr, nSubj, strjoin(string(iLocCorr_all), ''), nameVarY_fileTitle), 'PearsonR', 'SpearmanRho')

                % plot ub and lb when fitting
                yline(ub_full_all{iFamily}(iTunC), 'k--');
                yline(lb_full_all{iFamily}(iTunC), 'k--');

                xlabel('Contrast sensitivity (1/contrast)', 'fontsize', sz_label)
                ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunCs{iTunC}), 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                drawnow; % Ensure text/tick extents are up to date before reading TightInset
                ti = ax.TightInset; % [left bottom right top] padding needed for labels/ticks
                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset; % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset; % move label left

                newPos = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];
                if all(isfinite(newPos)) && newPos(3) >= 0.55 && newPos(4) >= 0.55
                    ax.Position = newPos;
                end

                saveAndCloseFigure(fig, sprintf('%s/n%d_L%s_%s.png', nameFolder_Fig_NOM_corr, nSubj, strjoin(string(iLocCorr_all), ''), nameVarY_fileTitle))

            end % iTunC
        end % iFeature
    end % iSet
    clear margParams*_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 9. Corr1: Corr between CS and tunC
    clc; fprintf('\n%s: Fig 9/23: Corr CS vs tuning characteristics STARTED.\n', string(datetime('now')))

    % Load data
    load(nameFolder_Output_SaveCompile, 'margParams*_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'tunC';
    sz_label = 50;
    sz_labelOffset = .03;
    sz_axOffset = 0.05; % extra breathing room

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

    nameFolder_Outputs_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Outputs_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Outputs_NOM_corr)); mkdir(nameFolder_Outputs_NOM_corr), end

    flag_zeroMean = 0;
    flag_plotIdvdCI = 0;
    flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles

    for iSet = 1:numel(iLocSingle_allSets)

        iLocCorr_all = iLocSingle_allSets{iSet};

        % Obtain CS (x-axis)
        X_allSubj = CS_allSubj(:, iLocCorr_all);
        x_ticks = linspace(1.5, 3.5, 5); % Asymmetry in CS
        X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);

        for iFeature = 1:nFeatures
            iFamily = iFamily_perF(iFeature);
            str_family = sprintf('F%d %s', iFamily, namesFamily_all{iFamily});
            namesTunCs = namesParams_all{iFamily};
            nTunCs_full = length(namesTunCs);

            for iTunC = 1:nTunCs_full
                fprintf(' - L%s %s%d [%s]\n', strjoin(string(iLocCorr_all), ''), namesFeature{iFeature}, iTunC, str_family)
                switch iFeature
                    case 1, NOMp_allIter_allSubj = squeeze(margParams_ORI_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
                    case 2, NOMp_allIter_allSubj = squeeze(margParams_SF_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
                end

                switch iFeature
                    case 1, y_ticks_allTunC_lb = [-.004, 0, -.1]; y_ticks_allTunC_ub = [.02, 6, .1];% 3 values are ORI peak amplitude, width, baseline
                    case 2, y_ticks_allTunC_lb = [-2, 0, .5, -.1]; y_ticks_allTunC_ub = [2, .004, 3.5, .1]; % 4 values are SF peak, peak amplitude, width, baseline
                end

                y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);

                x_ticklabels = x_ticks;
                y_ticklabels = y_ticks;
                if iFeature==2 && iTunC==1 % SF peak, needs to convert y tick label to linear scale
                    y_ticklabels = 2.^y_ticks;
                end

                nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunCs{iTunC});
                nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

                str_title = sprintf('n=%d, nIter=%d %s vs. %s [L%s]', nSubj, nIterxJob, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

                %----------------------------%
                [PearsonR, SpearmanRho] = fxn_drawCorr(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), flag_UseRUseRho, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, nIterxJob, flag_plotPurpose);
                %----------------------------%
                fig = gcf;
                % Save CI for automatic CI range calcuation in CorrAsym
                save(sprintf('%s/n%d_L%s_%s.mat', nameFolder_Outputs_NOM_corr, nSubj, strjoin(string(iLocCorr_all), ''), nameVarY_fileTitle), 'PearsonR', 'SpearmanRho')

                xlabel('Contrast sensitivity (1/contrast)', 'fontsize', sz_label)
                ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunCs{iTunC}), 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                drawnow; % Ensure text/tick extents are up to date before reading TightInset
                ti = ax.TightInset; % [left bottom right top] padding needed for labels/ticks
                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset; % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset; % move label left

                newPos = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];
                if all(isfinite(newPos)) && newPos(3) >= 0.55 && newPos(4) >= 0.55
                    ax.Position = newPos;
                end

                saveAndCloseFigure(fig, sprintf('%s/n%d_L%s_%s.png', nameFolder_Fig_NOM_corr, nSubj, strjoin(string(iLocCorr_all), ''), nameVarY_fileTitle))

            end % iTunC
        end % iFeature
    end % iSet
    clear margTunC*_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% Corr2: Corr between CS and pA
    clc; fprintf('\n%s: Fig 10/23: Corr CS vs pA STARTED.\n', string(datetime('now')))

    % Load data
    load(nameFolder_Output_SaveCompile, 'metric_data_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'pA'; iMetric_pA = 2;

    sz_label = 55;
    sz_labelOffset = .05;
    sz_axOffset = 0.05; % extra breathing room

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end
    nameFolder_Outputs_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Outputs_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Outputs_NOM_corr)); mkdir(nameFolder_Outputs_NOM_corr), end

    flag_zeroMean = 0;
    flag_plotIdvdCI = 0;
    flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles

    for iSet = 1:numel(iLocSingle_allSets)

        iLocCorr_all = iLocSingle_allSets{iSet};
        fprintf(' - L%s\n', strjoin(string(iLocCorr_all), ''))

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

        str_title = sprintf('n=%d, nIter=%d%s vs. %s [L%s]', nSubj, nIterxJob, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

        %----------------------------%
        [PearsonR, SpearmanRho] = fxn_drawCorr(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), flag_UseRUseRho, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, nIterxJob, flag_plotPurpose);
        %----------------------------%
        fig = gcf;
        % Save CI for automatic CI range calcuation in CorrAsym
        save(sprintf('%s/n%d_L%s_%s.mat', nameFolder_Outputs_NOM_corr, nSubj, strjoin(string(iLocCorr_all), ''), nameVarY_fileTitle), 'PearsonR', 'SpearmanRho')

        xlabel('Contrast sensitivity (1/contrast)', 'fontsize', sz_label)
        ylabel(nameVarY_figTitle, 'fontsize', sz_label)

        % Adjust distance between components
        ax = gca;
        drawnow; % Ensure text/tick extents are up to date before reading TightInset
        ti = ax.TightInset; % [left bottom right top] padding needed for labels/ticks
        ax.XLabel.Units = 'normalized';
        ax.YLabel.Units = 'normalized';

        ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset; % move label down
        ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset; % move label left

        newPos = [ ...
            ti(1) + sz_axOffset, ...
            ti(2) + sz_axOffset, ...
            1 - ti(1) - ti(3) - 2*sz_axOffset, ...
            1 - ti(2) - ti(4) - 2*sz_axOffset];
        if all(isfinite(newPos)) && newPos(3) >= 0.55 && newPos(4) >= 0.55
            ax.Position = newPos;
        end

        saveAndCloseFigure(fig, sprintf('%s/n%d_L%s.png', nameFolder_Fig_NOM_corr, nSubj, strjoin(string(iLocCorr_all), '')))

    end % iSet
    clear metric_data_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% CompAsym1: CS and tunC: compare binned estimates
    % clc; fprintf('\n%s: Fig 11/23: CompAsym CS vs tuning characteristics STARTED.\n', string(datetime('now')))
    %
    % % Load data
    % load(nameFolder_Data_SaveCompile, 'margTunC_*_allCond', 'CS_allSubj')
    %
    % nameVarX = 'CS';
    % nameVarY = 'tunC';
    %
    % nameFolder_Fig_NOM_CompAsym = sprintf('%s/CompAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    % if isempty(dir(nameFolder_Fig_NOM_CompAsym)); mkdir(nameFolder_Fig_NOM_CompAsym), end
    %
    % sz_label = 15;
    % sz_labelOffset = 0.01;
    % sz_axOffset = 0.05; % extra breathing room
    % nBinsCompAsym = 2;
    % sz_wd_perBar = 200;
    % sz_fig = [1e3, 500];
    %
    % for iGroup = 1:nGroups
    %     iLocPair_all = iLocGroups_all{iGroup};
    %     fprintf(' - L%s\n', strjoin(string(iLocPair_all), ''))
    %
    %     switch iLocPair_all(1)
    %         case 1, nameAsymX = 'Ecc. effect'; nameAsymY = 'Ecc. effect';
    %         case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
    %         case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
    %     end
    %
    %     % X-axis
    %     asymX_allSubj = (CS_allSubj(:, iLocPair_all(1))-CS_allSubj(:, iLocPair_all(2)))./(CS_allSubj(:, iLocPair_all(1))+CS_allSubj(:, iLocPair_all(2)));
    %     switch iLocPair_all(1)
    %         case 1, x_ticks = linspace(0, 20, 5); % EE
    %         case 6, x_ticks = linspace(-5, 15, 5); % HVA
    %         case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
    %     end
    %     asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob); % [nSubj x nIterxJob]
    %
    %     for iFeature = 1:nFeatures
    %         iFamily = iFamily_perF(iFeature);
    %         str_family = sprintf('F%d %s', iFamily, namesFamily_all{iFamily});
    %         if iFeature == 1
    %             namesTunCs = namesTunC_unit_perF{iFamily, 2};
    %             namesTunCs_noUnit = namesTunC_noUnit{iFamily, 2};
    %         else
    %             namesTunCs = namesTunC_SF_unit;
    %             namesTunCs_noUnit = namesTunC_SF_noUnit;
    %         end
    %         nTunCs_full = length(namesTunCs);
    %
    %         for iTunC = 1:nTunCs_full
    %             % Y-axis
    %             switch iFeature
    %                 case 1, NOMp_allIter_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
    %                 case 2, NOMp_allIter_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
    %             end
    %
    %             % Y-value: [nSubj x nIterxJob]
    %             asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));
    %
    %             switch iLocPair_all(1)
    %                 case 1 % EE
    %                     switch iFamily
    %                         case 1, y_ticks_allTunC_lb = -[0, 32, 100]; y_ticks_allTunC_ub = [72, 20, 100];
    %                             % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
    %                         case 2, y_ticks_allTunC_lb = -[50, 0, 50, 145]; y_ticks_allTunC_ub = [50, 80, 50, 265];
    %                         case 14, y_ticks_allTunC_lb = -[50, 0, 50, 145]; y_ticks_allTunC_ub = [50, 80, 50, 265];
    %                     end
    %                 case 6 % HVA
    %                     switch iFamily
    %                         case 1, y_ticks_allTunC_lb = -[40, 30, 100]; y_ticks_allTunC_ub = [40, 50, 100];
    %                             % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
    %                         case 2, y_ticks_allTunC_lb = -[50, 20, 60, 145]; y_ticks_allTunC_ub = [50, 60, 60, 265];
    %                         case 14, y_ticks_allTunC_lb = -[50, 20, 60, 145]; y_ticks_allTunC_ub = [50, 60, 60, 265];
    %                     end
    %                 case 5 % VMA
    %                     switch iFamily
    %                         case 1, y_ticks_allTunC_lb = -[30, 70, 50]; y_ticks_allTunC_ub = [90, 30, 50];
    %                             % case 8, y_ticks_allTunC_lb = -[100, 40, 300, 320, 50, 320]; y_ticks_allTunC_ub = [100, 44, 300, 320, 50, 320];
    %                         case 2, y_ticks_allTunC_lb = -[50, 40, 60, 250]; y_ticks_allTunC_ub = [50, 100, 40, 250];
    %                         case 14, y_ticks_allTunC_lb = -[50, 40, 60, 250]; y_ticks_allTunC_ub = [50, 100, 40, 250];
    %                     end
    %             end
    %             y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);
    %
    %             x_ticklabels = nan;
    %             y_ticklabels = nan;
    %
    %             nameVarY_figTitle = sprintf('%s %s [%s]', namesFeature{iFeature}, namesTunCs{iTunC}, str_family);
    %             nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);
    %
    %             str_title = sprintf('n=%d, nIter=%d %s (%s) vs. %s (%s)', nSubj, nIterxJob, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);
    %
    %             str_ylabel = sprintf('\\Delta %s %s (%%)', namesFeature{iFeature}, namesTunCs_noUnit{iTunC});
    %
    %             %----------------------------%
    %             fxn_compAsym(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, nBinsCompAsym, y_ticks, sz_fig, str_title, str_ylabel)
    %             %----------------------------%
    %
    %             % xlabel(sprintf('\\Delta contrast sensitivity (%%)'), 'FontSize', sz_label);
    %             % ylabel(str_ylabel, 'fontsize', sz_label)
    %
    %             % Adjust distance between components
    %             % ax = gca;
    %             % ti = ax.TightInset; % [left bottom right top] padding needed for labels/ticks
    %             %
    %             % ax.XLabel.Units = 'normalized';
    %             % ax.YLabel.Units = 'normalized';
    %             %
    %             % ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset; % move label down
    %             % ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset; % move label left
    %             % ax.Position = [ ...
    %             % ti(1) + sz_axOffset, ...
    %             % ti(2) + sz_axOffset, ...
    %             % 1 - ti(1) - ti(3) - 2*sz_axOffset, ...
    %             % 1 - ti(2) - ti(4) - 2*sz_axOffset];
    %
    %             saveas(gcf, sprintf('%s/n%d_L%d%d_%s.png', nameFolder_Fig_NOM_CompAsym, nSubj, iLocPair_all, nameVarY_fileTitle))
    %             close(gcf)
    %
    %         end % iTunC
    %     end % iFeature
    % end % iGroup
    % clear margTunC_*_allCond
    % fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 12. CorrAsym1: (CS and tunC): corr between extents of EE/HVA/VMA
    clc; fprintf('\n%s: Fig 12/23: CorrAsym CS vs tuning characteristics STARTED.\n', string(datetime('now')))

    flag_plotIdvdCI = 1;
    flag_plotUnikSymbol = 0;
    iLocCorr_all = [6,5,3];

    % Load data
    load(nameFolder_Output_SaveCompile, 'margParams_*_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'tunC';

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

    nameFolder_Outputs_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Outputs_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Outputs_NOM_corr)); mkdir(nameFolder_Outputs_NOM_corr), end

    sz_label = 55;
    sz_labelOffset = 0.01;
    sz_axOffset = 0.05; % extra breathing room

    for iGroup = 1:nGroups
        iLocPair_all = iLocGroups_all{iGroup};

        switch iLocPair_all(1)
            case 1, nameAsymX = 'EE'; nameAsymY = 'EE';
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
        x_ticks = linspace(-3, 17, 5); 
        asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';

        for iFeature = 1:nFeatures
            iFamily = iFamily_perF(iFeature);
            str_family = sprintf('F%d %s', iFamily, namesFamily_all{iFamily});
            if iFeature == 1
                namesTunCs = namesTunC_unit_perF{iFamily, 2};
                namesTunCs_noUnit = namesTunC_noUnit{iFamily, 2};
            else
                namesTunCs = namesTunC_SF_unit;
                namesTunCs_noUnit = namesTunC_SF_noUnit;
            end
            nTunCs_full = length(namesTunCs);

            for iTunC = 1:nTunCs_full
                fprintf(' - L%s %s%d\n', strjoin(string(iLocPair_all), ''), namesFeature{iFeature}, iTunC)

                % Y-axis
                switch iFeature
                    % case 1, NOMp_allIter_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                    % case 2, NOMp_allIter_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                        case 1, NOMp_allIter_allSubj = squeeze(margParams_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                    case 2, NOMp_allIter_allSubj = squeeze(margParams_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                end

                asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

                switch iLocPair_all(1)
                    case 1 % EE
                        switch iFamily
                            case 1, y_ticks_allTunC_lb = -[0, 40, 100]; y_ticks_allTunC_ub = [80, 40, 160];
                                % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                            case 2, y_ticks_allTunC_lb = -[40, -20, 50, 180]; y_ticks_allTunC_ub = [60, 70, 50, 200];
                            case 10, y_ticks_allTunC_lb = -[10, -20, 50]; y_ticks_allTunC_ub = [50, 70, 50];
                        end
                    case 6 % HVA
                        switch iFamily
                            % case 1, y_ticks_allTunC_lb = -[60, 60, 200]; y_ticks_allTunC_ub = [60, 60, 200]; % ORI—peak, width, baseline
                            case 2, y_ticks_allTunC_lb = -[50, 100, 60, 160]; y_ticks_allTunC_ub = [50, 100, 60, 200]; % SF—peak, peak amplitude, width, baseline
                            case 10, y_ticks_allTunC_lb = -[40, 18, 20]; y_ticks_allTunC_ub = [60, 10, 10]; % ORI—peak, width, baseline;
                        end
                    case 5 % VMA
                        switch iFamily
                            % case 1, y_ticks_allTunC_lb = -[100, 80, 200]; y_ticks_allTunC_ub = [100, 60, 200];
                                % case 8, y_ticks_allTunC_lb = -[100, 40, 300, 320, 50, 320]; y_ticks_allTunC_ub = [100, 44, 300, 320, 50, 320];
                            case 2, y_ticks_allTunC_lb = -[50, 120, 120, 300]; y_ticks_allTunC_ub = [70, 120, 120, 300];
                            case 10, y_ticks_allTunC_lb = -[80, 30, 100]; y_ticks_allTunC_ub = [100, 30, 100];
                        end
                end

                % Apply the same yticks to all panels
                switch iFamily
                    case 2, y_ticks_allTunC_lb = -[100, 100, 100, 100]; y_ticks_allTunC_ub = [100, 100, 100, 100];
                    case 10, y_ticks_allTunC_lb = -[100, 100, 100]; y_ticks_allTunC_ub = [100, 100, 100];
                end

                y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);
                x_ticklabels = nan;
                y_ticklabels = nan;

                nameVarY_figTitle = sprintf('%s %s [%s]', namesFeature{iFeature}, namesTunCs{iTunC}, str_family);
                nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

                str_title = sprintf('n=%d, nIter=%d %s (%s) vs. %s (%s)', nSubj, nIterxJob, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

                % Load corr analysis
                load(sprintf('%s/n%d_L%s_%s.mat', nameFolder_Outputs_NOM_corr, nSubj, strjoin(string(iLocCorr_all), ''), nameVarY_fileTitle), 'PearsonR', 'SpearmanRho');
                switch flag_UseRUseRho
                    case 'useR'
                        flag_CIrange = PearsonR(4) < 0.05; % use permutation p-value to determine whether one-tailed corr should be conducted
                        % flag_CIrange = PearsonR(2) * PearsonR(3)>0; % flag_CIrange=1 if r excludes 0, so one-tailed corr should be conducted, so CI range is 90%
                    case 'useRho'
                        flag_CIrange = SpearmanRho(4) < 0.05; % use permutation p-value to determine whether one-tailed corr should be conducted
                        % flag_CIrange = SpearmanRho(2) * SpearmanRho(3)>0; % flag_CIrange=1 if r excludes 0, so one-tailed corr should be conducted, so CI range is 95%
                end

                %----------------------------%
                fxn_drawCorrAsym(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, colors_comb(iLocPair_all, :), flag_UseRUseRho, flag_CIrange, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, flag_plotPurpose)
                %----------------------------%
                fig = gcf;

                xlabel(sprintf('%s of contrast sensitivity (%%)', nameAsymX), 'FontSize', sz_label);
                ylabel(sprintf('%s of %s %s (%%)', nameAsymY, namesFeature{iFeature}, namesTunCs_noUnit{iTunC}), 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                drawnow; % Ensure text/tick extents are up to date before reading TightInset
                ti = ax.TightInset; % [left bottom right top] padding needed for labels/ticks

                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset; % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset; % move label left
                newPos = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];
                if all(isfinite(newPos)) && newPos(3) >= 0.55 && newPos(4) >= 0.55
                    ax.Position = newPos;
                end

                saveAndCloseFigure(fig, sprintf('%s/n%d_L%d%d_%s.png', nameFolder_Fig_NOM_CorrAsym, nSubj, iLocPair_all, nameVarY_fileTitle))

            end % iTunC
        end % iFeature
    end % iGroup
    clear margTunC_*_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% CorrAsym2: (CS and pA): corr between extents of EE/HVA/VMA
    clc; fprintf('\n%s: Fig 13/23: CorrAsym CS vs pA STARTED.\n', string(datetime('now')))

    % Load data
    load(nameFolder_Output_SaveCompile, 'metric_data_allCond', 'CS_allSubj', 'pA_allSubj')

    nameVarX = 'CS';
    nameVarY = 'pA';

    iMetric_pA = 2;

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

    nameFolder_Outputs_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Outputs_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Outputs_NOM_corr)); mkdir(nameFolder_Outputs_NOM_corr), end

    sz_label = 60;
    sz_labelOffset = 0.02;
    sz_axOffset = 0.05; % extra breathing room

    for iGroup = 1:nGroups
        iLocPair_all = iLocGroups_all{iGroup};
        fprintf(' - L%s\n', strjoin(string(iLocPair_all), ''))

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
        % case 1 % EE
        % y_ticks_lb = -[10, 70, 12]; y_ticks_ub = [6, 40, 30];
        % case 6 % HVA
        % y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        % case 5 % VMA
        % y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        % end

        y_ticks = linspace(-10, 10, 5);

        x_ticklabels = nan;
        y_ticklabels = nan;

        nameVarY_figTitle = nameVarY;
        nameVarY_fileTitle = nameVarY_figTitle;

        str_title = sprintf('n=%d, nIter=%d %s (%s) vs. %s (%s)', nSubj, nIterxJob, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

        % Load corr analysis
        load(sprintf('%s/n%d_L%s_%s.mat', nameFolder_Outputs_NOM_corr, nSubj, strjoin(string(iLocCorr_all), ''), nameVarY_fileTitle), 'PearsonR', 'SpearmanRho');
        switch flag_UseRUseRho
            case 'useR'
                flag_CIrange = PearsonR(4) < 0.05; % use permutation p-value to determine whether one-tailed corr should be conducted
                % flag_CIrange = PearsonR(2) * PearsonR(3)>0; % flag_CIrange=1 if r excludes 0, so one-tailed corr should be conducted, so CI range is 90%
            case 'useRho'
                flag_CIrange = SpearmanRho(4) < 0.05; % use permutation p-value to determine whether one-tailed corr should be conducted
                % flag_CIrange = SpearmanRho(2) * SpearmanRho(3)>0; % flag_CIrange=1 if r excludes 0, so one-tailed corr should be conducted, so CI range is 95%
        end

        %----------------------------%
        fxn_drawCorrAsym(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, ...
            colors_comb(iLocPair_all, :), flag_UseRUseRho, flag_CIrange, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, flag_plotPurpose)
        %----------------------------%
        fig = gcf;

        xlabel(sprintf('%s of contrast sensitivity (%%)', nameAsymX), 'fontsize', sz_label)
        ylabel(sprintf('%s of %s (%%)', nameAsymX, nameVarY_figTitle), 'fontsize', sz_label)

        % Adjust distance between components
        ax = gca;
        drawnow; % Ensure text/tick extents are up to date before reading TightInset
        ti = ax.TightInset; % [left bottom right top] padding needed for labels/ticks

        ax.XLabel.Units = 'normalized';
        ax.YLabel.Units = 'normalized';

        ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset; % move label down
        ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset; % move label left
        newPos = [ ...
            ti(1) + sz_axOffset, ...
            ti(2) + sz_axOffset, ...
            1 - ti(1) - ti(3) - 2*sz_axOffset, ...
            1 - ti(2) - ti(4) - 2*sz_axOffset];
        if all(isfinite(newPos)) && newPos(3) >= 0.55 && newPos(4) >= 0.55
            ax.Position = newPos;
        end

        saveAndCloseFigure(fig, sprintf('%s/n%d_L%d%d.png', nameFolder_Fig_NOM_CorrAsym, nSubj, iLocPair_all))

    end % iGroup
    clear metric_data_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 14. [NOM] Metrics vs. DV; group averages
    clc; fprintf('\n%s: Fig 14/23: NOM metrics vs DV (group averages) STARTED.\n', string(datetime('now')))

    % Load data
    load(nameFolder_Output_SaveCompile);
    % Define folder for saving figures for each model A and model B
    nameFolder_Fig_NOM_metrics = sprintf('%s/NOMmetrics_A%d_group', nameFolder_Fig_NOM_Trialwise, iModelA_plot);
    if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end

    x_base = 0.60; % move a bit further left to make space for name column
    x_width = 0.4; % total width of the table block
    y_base = 0.1; % bottom of the R2 table; increase to move up; decrease to move down
    y_height = 0.3; % height of the R2 table; increase to make table taller, decrease to make it shorter

    lineStyle_all = {'-', '--', ':', '-.'};
    % Build plot modes from iModelB_plot_all:
    % 1) DataOnly, 2) AllModels, 3..end) one figure per model.
    nModes = numel(iModelB_plot_all) + 2;
    modeLabels = cell(1, nModes);
    modeModelB_all = cell(1, nModes);
    modePlotPred = false(1, nModes);
    modeScale = ones(1, nModes);

    modeLabels{1} = 'DataOnly';
    modeModelB_all{1} = iModelB_plot_all(1); % source for plotting data only
    modePlotPred(1) = false;

    modeLabels{2} = 'AllModels';
    modeModelB_all{2} = iModelB_plot_all; % only models requested in iModelB_plot_all
    modePlotPred(2) = true;

    for iModeSingle = 1:numel(iModelB_plot_all)
        iB = iModelB_plot_all(iModeSingle);
        modeLabels{iModeSingle + 2} = sprintf('B%d', iB);
        modeModelB_all{iModeSingle + 2} = iB;
        modePlotPred(iModeSingle + 2) = true;
    end

    szScaling = 10; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end
    szBase = 8;
    sz_font = 22; % full model: xx; reduced models: 15
    wd = 2;

    for iSet = 1:numel(iLocSingle_allSets)

        iLocSingle_perSet = iLocSingle_allSets{iSet};

        fprintf('\nL%s...', strjoin(string(iLocSingle_perSet), ''))

        for iPlotMode = 1:numel(modeLabels)

            % --------------------------
            % Decide what to plot
            % --------------------------
            str_plotMode = modeLabels{iPlotMode};
            fprintf(' %s ', str_plotMode)
            scaleMode = modeScale(iPlotMode);
            iModelB_all_plot = modeModelB_all{iPlotMode};
            flag_plotPred = modePlotPred(iPlotMode);
            flag_showR2table = flag_plotPred;

            nModelB_plot = numel(iModelB_all_plot);

            %%% --- PLOT ---
            for iMetric_prob = 1:nMetrics_prob
                fig = figure('Position', [0 0 500/scaleMode 500/scaleMode]); hold on

                % R2 table: rows = plotted ModelB, cols = locations in this set
                % R2_tab = nan(numel(iLocSingle_perSet), 2); % ave and sem

                for rModel = 1:nModelB_plot
                    iModelB = iModelB_all_plot(rModel);

                    % Preallocate
                    R2_w_ave = nan(numel(iLocSingle_perSet), 1);
                    R2_w_lb = R2_w_ave;
                    R2_w_ub = R2_w_ave;
                    NRMSE_w_ave = R2_w_ave;
                    NRMSE_w_lb = R2_w_ave;
                    NRMSE_w_ub = R2_w_ave;

                    for iiLoc = 1:numel(iLocSingle_perSet)
                        % ===== Bootstrap group averages (DATA + PRED) and plot median + 68% CI bands =====
                        locID = iLocSingle_perSet(iiLoc);

                        % --------------------------
                        % 1) Extract to [nBins x nSubj x nIter] (robust to extra singleton dims)
                        % --------------------------
                        DV_raw = squeeze(DV_allCond(iModelA_plot, iModelB, locID, :, :, :)); % -> [nBins x nSubj x nIter] (expected)
                        nTrials_raw = squeeze(nTrials_allCond(iModelA_plot, iModelB, locID, :, :, :));
                        data_raw = squeeze(metric_data_allCond(iModelA_plot, iModelB, locID, iMetric_prob, :, :, :));
                        pred_raw = squeeze(metric_pred_allCond(iModelA_plot, iModelB, locID, iMetric_prob, :, :, :));

                        % Force canonical shape: [nBins x nSubj x nIter]
                        canonize = @(X) local_canonize_to_bins_subj_iter(X, nBins, nSubj);
                        DV_raw = canonize(DV_raw);
                        nTrials_raw = canonize(nTrials_raw);
                        data_raw = canonize(data_raw);
                        pred_raw = canonize(pred_raw);

                        % --------------------------
                        % 2) Within-subject summary across iterations (median over iter) -> [nBins x nSubj]
                        % --------------------------
                        DV_medPerSubj = median(DV_raw, 3, 'omitnan');
                        nTrials_medPerSubj = median(nTrials_raw, 3, 'omitnan');
                        data_medPerSubj = median(data_raw, 3, 'omitnan');
                        pred_medPerSubj = median(pred_raw, 3, 'omitnan');

                        % --------------------------
                        % 3) Bootstrap across subjects: group average per bin -> boot arrays: [nBoot x nBins]
                        % --------------------------
                        DV_allBoot = nan(nBoot, nBins);
                        nTrials_allBoot = nan(nBoot, nBins);
                        data_allBoot = nan(nBoot, nBins);
                        pred_allBoot = nan(nBoot, nBins);
                        R2_w_allBoot = nan(nBoot, 1);
                        NRMSE_w_allBoot = nan(nBoot, 1);

                        % Pregenerate subj indices
                        rng(seedBoot, 'twister');
                        indRand_allBoot = randi(nSubj, [nBoot, nSubj], 'uint16');
                        parfor iBoot = 1:nBoot
                            % idindRandx = randi(nSubj, [1, nSubj]); % resample subjects with replacement
                            indRandBoot = double(indRand_allBoot(iBoot,:));
                            DV_allBoot(iBoot,:) = mean(DV_medPerSubj(:,indRandBoot), 2, 'omitnan').';
                            nTrials_allBoot(iBoot,:) = mean(nTrials_medPerSubj(:,indRandBoot), 2, 'omitnan').';
                            data_allBoot(iBoot,:) = mean(data_medPerSubj(:,indRandBoot), 2, 'omitnan').';
                            pred_allBoot(iBoot,:) = mean(pred_medPerSubj(:,indRandBoot), 2, 'omitnan').';

                            % Calculate weighted R2
                            y    = data_allBoot(iBoot,:);
                            yhat = pred_allBoot(iBoot,:);
                            w    = nTrials_allBoot(iBoot,:);
                            m = isfinite(y) & isfinite(yhat) & isfinite(w) & (w > 0);                             % keep valid entries
                            y = y(m); yhat = yhat(m); w = w(m);
                            ybar_w = sum(w .* y) / sum(w); % weighted mean
                            SSE_w = sum(w .* (y - yhat).^2); % weighted SST
                            SST_w = sum(w .* (y - ybar_w).^2);
                            R2_w_allBoot(iBoot) = 1 - (SSE_w / SST_w);

                            % Calculate weighted normalized RMSE
                            wRMSE = sqrt(SSE_w / sum(w));
                            yrng = max(y) - min(y);
                            NRMSE_w_allBoot(iBoot) = wRMSE / yrng;   % normalized RMSE (by data range)
                        end % iBoot

                        % --------------------------
                        % 4) Point + interval estimates from bootstrap distribution
                        % Use getCI() (median + CI) across BOOT dimension (dim=1)
                        % --------------------------
                        [DV_ave, DV_lb, DV_ub] = getCI(DV_allBoot, 1, 1, CI68); % each is [1 x nBins]
                        [nTrials_ave, ~, ~] = getCI(nTrials_allBoot, 1, 1, CI68); % [1 x nBins] (for marker size only)
                        [data_ave, data_lb, data_ub, data_sem_neg, data_sem_pos] = getCI(data_allBoot, 1, 1, CI68);
                        [pred_ave, pred_lb, pred_ub] = getCI(pred_allBoot, 1, 1, CI68);
                        [R2_w_ave(iiLoc), R2_w_lb(iiLoc), R2_w_ub(iiLoc)] = getCI(R2_w_allBoot, 1, 1, CI95);
                        [NRMSE_w_ave(iiLoc), NRMSE_w_lb(iiLoc), NRMSE_w_ub(iiLoc)] = getCI(NRMSE_w_allBoot, 1, 1, CI95);

                        % Make column vectors for patch/plot
                        DV_ave = DV_ave(:); DV_lb = DV_lb(:); DV_ub = DV_ub(:);
                        data_ave = data_ave(:); data_lb = data_lb(:); data_ub = data_ub(:);
                        pred_ave = pred_ave(:); pred_lb = pred_lb(:); pred_ub = pred_ub(:);
                        nTrials_ave = nTrials_ave(:);

                        if iMetric_prob==1 % detection rate
                            data_ave_detection = data_ave;
                            pred_ave_detection  = pred_ave;
                            pA_data_detection = data_ave_detection.^2+(1-data_ave_detection).^2;
                            pA_pred_detection = pred_ave_detection.^2+(1-pred_ave_detection).^2;
                        end

                        % Optional: enforce monotonic x for clean plotting
                        % [DV_ave, iSort] = sort(DV_ave);
                        % data_ave = data_ave(iSort); data_lb = data_lb(iSort); data_ub = data_ub(iSort);
                        % pred_ave = pred_ave(iSort); pred_lb = pred_lb(iSort); pred_ub = pred_ub(iSort);
                        % nTrials_ave = nTrials_ave(iSort);

                        % ---------
                        % Styling
                        % ---------
                        cLoc = colors_comb(locID, :);

                        if numel(iModelB_all_plot) > 1
                            ls = '-';%lineStyle_all{iModelB};
                            if iModelB == 1, lw = 3; else, lw = 1.5; end
                        else
                            ls = '-';
                            lw = 3;
                        end
                        lw = lw/scaleMode;

                        % --------------------------
                        % Plot prediction (optional): median + 68% CI shaded band
                        % --------------------------
                        if flag_plotPred
                            % Shaded band (y uncertainty only; x is DV_ave)
                            patch([DV_ave; flipud(DV_ave)], [pred_lb; flipud(pred_ub)], cLoc, 'FaceAlpha', .15, 'LineStyle', 'none', 'HandleVisibility', 'off');

                            % Median line
                            plot(DV_ave, pred_ave, 'LineStyle', ls, 'Color', cLoc, 'LineWidth', lw, 'HandleVisibility', 'on');
                        end

                        % --------------------------
                        % Plot measurement (always): median + 68% CI shaded band + dots
                        % --------------------------
                        % Median points (size ~ trials)
                        for iBin = 1:nBins
                            errorbar(DV_ave(iBin), data_ave(iBin), data_sem_neg(iBin), data_sem_pos(iBin), 'vertical', 'color', cLoc, 'capsize', 0, 'LineWidth', lw, 'HandleVisibility', 'off');

                            plot(DV_ave(iBin), data_ave(iBin), 's', ...
                                'MarkerEdgeColor', cLoc, ...
                                'MarkerFaceColor', 'w', ...
                                'MarkerSize', (nTrials_ave(iBin) / szScaling + szBase)/scaleMode, ...
                                'LineWidth', max(lw/1.5, 1), ...
                                'LineStyle', 'none', ...
                                'HandleVisibility', 'off');

                        end % iBin

                        % Plot emp detection rate converted to consistency
                        % if iMetric_prob==2
                        %     plot(DV_ave, pA_data_detection, 'o', 'color', cLoc)
                        %     plot(DV_ave, pA_pred_detection, '-', 'color', cLoc)
                        % end

                    end % iiLoc

                end % rModel

                % --------------------------
                % Add size legend (by plotting symbols with integer sizes)
                % --------------------------
                legend_nTrials_all = [10 25 50];
                ms = (legend_nTrials_all./szScaling + szBase) ./ scaleMode;
                hLegend = gobjects(numel(legend_nTrials_all),1);
                for iSz = 1:numel(legend_nTrials_all)
                    hLegend(iSz) = plot(nan, nan, 's', 'MarkerEdgeColor','k','MarkerFaceColor','w', 'MarkerSize', ms(iSz), 'LineWidth', wd, 'LineStyle', 'none');
                end
                lgd = legend(hLegend, compose('%d trials', legend_nTrials_all), 'Location','northwest', 'Box','off');
                lgd.Title.String = '# trials';

                % --------------------------
                % Draw reference line and y formatting
                % --------------------------
                yline(.5, '--', 'LineWidth', 2/scaleMode, 'Color', ones(1,3)/2, 'HandleVisibility', 'off');
                
                switch iMetric_prob
                    case 1, ylim([0, 1]); yticks(0:.2:1)
                    case 2, ylim([.45, 1]); yticks(.5:.1:1)
                end
                ylabel(namesMetrics_prob_full{iMetric_prob})

                x_ticks = 0:2:6;
                xlim([0,6])
                xticks(x_ticks);
                xticklabels(x_ticks);
                xlabel('Binned decision variable')

                % --------------------------
                % Plot weighted R2 summary
                % --------------------------
                if flag_plotPred
                    ax = gca;

                    % nRows = nModelB_plot;
                    nRows = numel(iLocSingle_perSet);

                    % --- layout (more horizontal spacing) ---
                    % centers of each column/row
                    x_cells = 0.55/scaleMode;
                    y_cells = y_base + ((nRows+1:-1:1) - 0.5*scaleMode) * (y_height / (nRows+1)); % 1 x nRows (top to bottom)

                    % --- draw numbers ---
                    for iRow = 0:nRows
                        if iRow==0
                            str_cell = sprintf('Weighted R^2:');
                            c='k';
                        else
                            str_cell = sprintf('%.0f%% [%.0f%%, %.0f%%]', 100*R2_w_ave(iRow), 100*R2_w_lb(iRow), 100*R2_w_ub(iRow));
                            c = colors_comb(iLocSingle_perSet(iRow), :);
                        end

                        text(x_cells, y_cells(iRow+1), str_cell, ...
                            'Units','normalized', ...
                            'HorizontalAlignment','left', ...
                            'VerticalAlignment','middle', ...
                            'FontSize', sz_font*scaleMode/1.5, ...
                            'FontWeight','normal', ...
                            'Interpreter','tex', ...
                            'Color', c);
                    end % iRow
                end
                % axis cosmetics
                ax = gca;
                ax.FontSize = sz_font/scaleMode;
                ax.LineWidth = wd/scaleMode/1.5;

                title(sprintf('n=%d [A%d] [L%s] [nIter=%d] %s | %s', ...
                    nSubj, iModelA_plot, strjoin(string(iLocSingle_perSet), ''), nIterxJob, ...
                    namesMetrics_prob{iMetric_prob}, str_plotMode), ...
                    'FontSize', 10/scaleMode);

                % Save
                saveAndCloseFigure(fig, sprintf('%s/n%d_L%s_%s_%s.png', ...
                    nameFolder_Fig_NOM_metrics, nSubj, strjoin(string(iLocSingle_perSet), ''), namesMetrics_prob{iMetric_prob}, str_plotMode));

            end % iMetric_prob
        end % iPlotMode
    end % iSet

    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 15. [NOM] Metrics vs. DV; per idvd
    clc; fprintf('\n%s: Fig 15/23: NOM metrics vs DV (individual) STARTED.\n', string(datetime('now')))
    lineStyle_all = {'-', '-', '--', ':', '-.', '-', '--', ':', '-.'};

    % iModelA_plot = 1;%iModelA_all % just plot pred of the model using RC-derived template (happens to be yhe best model)
    % iModelA = iModelA_plot;
    iModelB_all_plot = 1;

    % Load data
    load(nameFolder_Output_SaveCompile);

    % Define folder for saving figures for each model A and model B
    nameFolder_Fig_NOM_metrics = sprintf('%s/NOMmetrics_A%d_idvd', nameFolder_Fig_NOM_Trialwise, iModelA_plot);
    if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end

    for iLocSingle = 1:nLocComb8
        fprintf('\nL%d...', iLocSingle)
        szScaling = 80; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end

        %%%%%%%%%%%%
        for iMetric_prob = 1:nMetrics_prob
            fig = figure('Position', [0 0 2e3 1.5e3]);

            for isubj = 1:nSubj
                subplot(nRows_subj, nCols_subj, isubj), hold on

                subjName = subjList{isubj};

                % Create a string to store R2
                % str_R2 = 'R2: ';
                % str_NRMSE = 'NRMSE: ';
                str_title = '';

                % Create a string to store parameter estimates for display
                str_est = [];
                params_allCond(isnan(params_allCond))=0;
                for iModelB = iModelB_all_plot
                    vals = strjoin(string(round(getCI(params_allCond(iModelA_plot, iModelB, iLocSingle, isubj, :, :), 1, 5)',2)), ", ");
                    str_est = [str_est, sprintf('B%d: %s\n', iModelB, vals)];
                end

                for iModelB = iModelB_all_plot

                    % Compute medians and CIs
                    [DV_allBins, ~, ~, DV_allBins_neg, DV_allBins_pos] = getCI(DV_allCond(iModelA_plot, iModelB, iLocSingle, isubj, :, :), 1, 5);
                    [nTrials_allBins, nTrials_allBins_lb, nTrials_allBins_ub] = getCI(nTrials_allCond(iModelA_plot, iModelB, iLocSingle, isubj, :, :), 1, 5);
                    [data_allBins, ~, ~, data_allBins_neg, data_allBins_pos] = getCI(metric_data_allCond(iModelA_plot, iModelB, iLocSingle, iMetric_prob, isubj, :, :), 1, 6);
                    [pred_allBins, pred_allBins_lb, pred_allBins_ub] = getCI(metric_pred_allCond(iModelA_plot, iModelB, iLocSingle, iMetric_prob, isubj, :, :), 1, 6);
                    [params_allBins, params_allBins_lb, params_allBins_ub] = getCI(params_allCond(iModelA_plot, iModelB, iLocSingle, isubj, :, :), 1, 5);
                    [R2_w_NOM_med, R2_w_NOM_lb, R2_w_NOM_ub] = getCI(R2_w_NOM_allCond(iModelA_plot, iModelB, iLocSingle, iMetric_prob, isubj, :), 1, 6, .95);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Plot prediction
                    if iModelB==1, color_pred = colors_comb(iLocSingle, :); lw = 3;
                    else, color_pred='k'; lw = 1;
                    end
                    patch([DV_allBins; flip(DV_allBins)], [pred_allBins_lb; flip(pred_allBins_ub)], ones(1, 3) / 2, 'FaceAlpha', .3, 'linestyle', 'none', 'handlevisibility', 'off')
                    plot(DV_allBins, pred_allBins, 'lineStyle', lineStyle_all{iModelB}, 'color', color_pred, 'linewidth', lw)

                    % Plot measurement
                    errorbar(DV_allBins, data_allBins, DV_allBins_neg, DV_allBins_pos, '.', 'horizontal', 'CapSize', 0, 'color', colors_comb(iLocSingle, :), 'handlevisibility', 'off', 'linewidth', lw)
                    errorbar(DV_allBins, data_allBins, data_allBins_neg, data_allBins_pos, '.', 'vertical', 'CapSize', 0, 'color', colors_comb(iLocSingle, :), 'handlevisibility', 'off', 'linewidth', lw)

                    yline(.5, 'k--');

                    % Plot averaged data of each bin
                    for iBin = 1:nBins
                        % if isubj>nMarkersMax, facecolor=colors_comb(iLocComb, :); else, facecolor='w'; end
                        facecolor = 'w';
                        plot(DV_allBins(iBin), data_allBins(iBin), 'o', 'markeredgecolor', colors_comb(iLocSingle, :), ...
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

                    % Median and CI of GoF (trial-weighted R2)
                    str_title = sprintf('%s\n[B%d] wR^2: %.2f [%.2f, %.2f]\n', str_title, iModelB, R2_w_NOM_med, R2_w_NOM_lb, R2_w_NOM_ub);

                end % iModelB

                title(sprintf('[%s] %s', subjName, str_title), 'fontsize', 20)
            end % isubj

            sgtitle(sprintf('n=%d [A%d] [L%d] [nIter=%d] %s', nSubj, iModelA_plot, iLocSingle, nIterxJob, namesMetrics_prob{iMetric_prob}))
            % set(findall(gcf, '-property', 'fontsize'), 'fontsize', 15)

            % Save the figure
            saveAndCloseFigure(fig, sprintf('%s/n%d_%s_L%d_A%d.png', nameFolder_Fig_NOM_metrics, nSubj, namesMetrics_prob{iMetric_prob}, iLocSingle, iModelA_plot))
        end % iMetric

        fprintf('DONE\n')
    end % iLocComb

    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% [NOM] ANOVA on nLL (ModelA x ModelB x Loc; with vars collapsed)
    clc; fprintf('\n%s: Fig 16/23: ANOVA nLL (collapsed factors) STARTED.\n', string(datetime('now')))

    % Load data
    load(nameFolder_Output_SaveCompile, 'nLL_allCond')

    % Define folder for saving GoF figures
    nameFolder_Fig_NOM_nLL = fullfile(nameFolder_Fig_NOM_Trialwise, 'NOM_NLL');
    if isempty(dir(nameFolder_Fig_NOM_nLL)), mkdir(nameFolder_Fig_NOM_nLL); end

    % [NOM] Plot comparison when collapsing modelB or location
    iModelB_allSets = {[1:7]};
    DimCollapse_all = {'ModelB', 'Loc'};
    switch flag_plotPurpose
        case 'paper', sz_label = 20;
        case 'slide', sz_label = 30;
    end

    for iSetModelB = 1:numel(iModelB_allSets)
        iModelB_all = iModelB_allSets{iSetModelB};

        for iSetLoc = 1:numel(iLocSingle_allSets)

            iLocSingle_perSet = iLocSingle_allSets{iSetLoc};

            % --------------------------
            % 2-way: modelB x single Loc when modelA=1 (data-derived template)
            % --------------------------
            iModelA_selected = 1;
            nLL_ANOVA = squeeze(nLL_allCond(iModelA_selected, iModelB_all, iLocSingle_perSet, :, :)); % [nModelB x nLoc x nSubj x nIter]
            [stats_2way, str_ANOVA2] = rm2ANOVA_A1(nLL_ANOVA, nPerm, CI95);

            % --------------------------
            % [NEW] Bootstrap CIs for effect sizes (partial eta^2)
            % - analyses are run on the same summaries you analyze/plot: median across CV iterations
            % - bootstrap is over observers (rows) with replacement
            % --------------------------
            nModelB = numel(iModelB_all);
            nLoc = numel(iLocSingle_perSet);

            % Subject-level cell means (median across iterations)
            NLL_med = squeeze(median(nLL_ANOVA, 4, 'omitnan')); % [nModelB x nLoc x nSubj]
            NLL_med = permute(NLL_med, [3, 1, 2]); % [nSubj x nModelB x nLoc]

            % Observed partial eta^2 on the subject-level table
            [eta2p_obs_ModelB, eta2p_obs_Loc, eta2p_obs_Int] = fxn_eta2p_rm2(NLL_med);

            % Bootstrap distribution of partial eta^2
            eta2p_ModelB_allBoot = nan(nBoot, 1);
            eta2p_Loc_allBoot = nan(nBoot, 1);
            eta2p_Int_allBoot = nan(nBoot, 1);

            % Pregenerate subj indices
            rng(seedBoot, 'twister');
            indRand_allBoot = randi(nSubj, [nBoot, nSubj], 'uint16');
            parfor iBoot = 1:nBoot
                indRandBoot = double(indRand_allBoot(iBoot, :));
                nLL_b = NLL_med(indRandBoot, :, :); % [nSubj x nModelB x nLoc]
                [eta2p_ModelB_allBoot(iBoot), eta2p_Loc_allBoot(iBoot), eta2p_Int_allBoot(iBoot)] = fxn_eta2p_rm2(nLL_b);
            end

            % Point + interval estimates from bootstrap (keep your getCI convention)
            [eta2p_ModelB_med, eta2p_ModelB_lb, eta2p_ModelB_ub] = getCI(eta2p_ModelB_allBoot, 1, 1, CI95);
            [eta2p_Loc_med, eta2p_Loc_lb, eta2p_Loc_ub] = getCI(eta2p_Loc_allBoot, 1, 1, CI95);
            [eta2p_Int_med, eta2p_Int_lb, eta2p_Int_ub] = getCI(eta2p_Int_allBoot, 1, 1, CI95);

            % A compact string you can append to the panel title (effect size only)
            str_eta2p = sprintf('eta2p (boot %.0f%%CI): ModelB=%.2f [%.2f, %.2f], Loc=%.2f [%.2f, %.2f], Int=%.2f [%.2f, %.2f]', ...
                CI95*100, eta2p_ModelB_med, eta2p_ModelB_lb, eta2p_ModelB_ub, ...
                eta2p_Loc_med, eta2p_Loc_lb, eta2p_Loc_ub, ...
                eta2p_Int_med, eta2p_Int_lb, eta2p_Int_ub);

            fig = figure('Position', [0 0 1e3 500]);

            for iDimCollapse = 1:2
                subplot(1,2,iDimCollapse);
                DimCollapse = DimCollapse_all{iDimCollapse};

                % Extract (collapse dim iDimCollapse on [nModelB x nLoc x nSubj x nIter])
                nLL_collapse = squeeze(mean(squeeze(nLL_allCond(iModelA_selected, iModelB_all, iLocSingle_perSet, :, :)), iDimCollapse, 'omitnan')); % [nCond x nSubj x nIter]
                nLL_collapse = permute(nLL_collapse, [3,2,1]); % -> [nIter x nSubj x nCond]

                % For each subject and each iteration, subtract the minimum across conditions (models/locs)
                dnLL_collapse = nan(size(nLL_collapse));
                for isubj = 1:nSubj
                    parfor iIter = 1:nIterxJob
                        d = nLL_collapse(iIter, isubj, :);
                        d_min = min(d(:));
                        dnLL_collapse(iIter, isubj, :) = d - d_min;
                    end
                end
                dnLL_collapse = squeeze(dnLL_collapse);

                str_title = sprintf('n=%d nIter=%d [L%s] [ModelB %s] %s collapsed', ...
                    nSubj, nIterxJob, strjoin(string(iLocSingle_perSet), ''), strjoin(string(iModelB_all), ''), DimCollapse);

                ref = nan;
                switch DimCollapse
                    case 'ModelB'
                        colors = colors_comb(iLocSingle_perSet, :);
                        x_ticks = namesLocComb(iLocSingle_perSet);
                    case 'Loc'
                        colors = repmat(linspace(0, .5, numel(iModelB_all))', 1, 3);
                        x_ticks = namesModelB(iModelB_all);
                end

                y_ticks = nan;
                y_ticklabels = nan;
                sz_fig = nan; % set to be nan to not create a figure inside the fxn
                flag_plotIDVD = 1;
                flag_plotDiff = 1;
                sz_text = 35;
                wd = 2;

                %------------------------------%
                fxn_drawBars(dnLL_collapse, ref, colors, x_ticks, y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj, sz_text, wd, flag_plotPurpose);
                %------------------------------%

                ylabel(sprintf('\\Delta nLL'), 'FontSize', sz_label);

            end % iDimCollapse

            sgtitle(sprintf('n=%d nIter=%d [L%s] [ModelB %s]\n%s\n%s\n%s', ...
                nSubj, nIterxJob, strjoin(string(iLocSingle_perSet), ''), strjoin(string(iModelB_all), ''), ...
                str_ANOVA2, str_eta2p, ...
                sprintf('eta2p_obs (median-iter): ModelB=%.2f, Loc=%.2f, Int=%.2f', eta2p_obs_ModelB, eta2p_obs_Loc, eta2p_obs_Int)));

            % save
            saveAndCloseFigure(fig, fullfile(nameFolder_Fig_NOM_nLL, sprintf('n%d_L%s_ModelB%s.png', nSubj, strjoin(string(iLocSingle_perSet), ''), strjoin(string(iModelB_all), ''))));

        end % iSetLoc
    end % iSetModelB

    fprintf('\n%s: DONE.\n', string(datetime('now')))


    %% 17. [NOM] ANOVA on nLL (ModelB x Loc; NO vars collapsed)
    % Load data
    load(nameFolder_Output_SaveCompile, 'nLL_allCond')
    clc; fprintf('\n%s: Fig 17/23: ANOVA nLL (no collapse) STARTED.\n', string(datetime('now')))

    % Define folder for saving GoF figures
    nameFolder_Fig_NOM_nLL = fullfile(nameFolder_Fig_NOM_Trialwise, 'NOM_NLL');
    if isempty(dir(nameFolder_Fig_NOM_nLL)), mkdir(nameFolder_Fig_NOM_nLL); end

    iModelA_selected = 1;

    % ModelB ordering + plotting flags
    iModelB_selected = 1:7; y_ticks = linspace(0, 20, 5);
    iModelB_selected = 1:4; y_ticks = linspace(0, 6, 5);
    pairs_bracket = [1 2; 1 3; 1 4]; % requested comparisons
    pairs_bracket = [1 2]; % requested comparisons
    % pairs_bracket = nan; % no comparison needed
    flag_plotIDVD = 0;
    flag_plotDiff = 0;

    switch flag_plotPurpose
        case 'paper', sz_label = 20;
        case 'slide', sz_label = 30;
    end
    sz_text = sz_label;
    wd = 3;
    y_ticklabels = nan;

    % Loop over locations
    for iLocSingle = 1:nLocComb8

        nBars = numel(iModelB_selected);
        switch flag_plotPurpose
            case 'paper', sz_fig = [nBars * 180, 200+nchoosek(nBars, 2)*40];
            case 'slide', sz_fig = [nBars * 180, 300+nchoosek(nBars, 2)*40];
        end

        % Extract raw nLL: [ModelB x Subj x Iter]
        nLL_allIter = squeeze(nLL_allCond(iModelA_selected, iModelB_selected, iLocSingle, :, :, :));

        % Sanity reshape to [nBars x nSubj x nIter]
        if ndims(nLL_allIter) == 2
            % If no iteration dimension survived, force nIter=1
            nLL_allIter = reshape(nLL_allIter, [nBars, nSubj, 1]);
        elseif ndims(nLL_allIter) == 3
            % assume [nBars x nSubj x nIter] already
        else
            error('Unexpected nLL dimensionality after squeeze: ndims=%d', ndims(nLL_allIter));
        end

        NLL_med = getCI(nLL_allIter, 1, 3); % [nCond x nSubj]
        NLL_min_perSubj = min(NLL_med, [], 1); % [1 x nSubj]
        dNLL_med = NLL_med - NLL_min_perSubj; % [nBars x nSubj], >=0

        % ------------------------------------------------------------
        % fxn_drawBars expects [nIter x nSubj x nCond]
        % We now have only one "iteration" (the median-collapsed value), so set nIter=1.
        % ------------------------------------------------------------
        dNLL_allIter_allSubj = nan(1, nSubj, nBars);
        dNLL_allIter_allSubj(1,:,:) = dNLL_med.'; % transpose -> [nSubj x nBars]

        % Strings / plotting params
        str_title = sprintf('n=%d | Loc L%d | ModelB [%s] | nIter=%d', nSubj, iLocSingle, strjoin(string(iModelB_selected), ' '), nIterxJob);

        ref = nan;

        % ------------------------------%
        [pperm_allPairs , pairs] = fxn_drawBars(dNLL_allIter_allSubj, ref, repmat(colors_comb(iLocSingle, :), nBars, 1), namesModelB(iModelB_selected), y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, 1, markers_allSubj, sz_text, wd, flag_plotPurpose);
        % ------------------------------%
        figMain = gcf;
        % ylim(y_ticks([1, end]))

        ylabel('\Delta NLL', 'FontSize', sz_label);

        % dnLL_allIter_allSubj must be [nIter x nSubj x nCond]
        [~, nSubj, nCond] = size(dNLL_allIter_allSubj);

        % (1) Subject-level summary across iterations (median)
        data_med_allSubj = squeeze(median(dNLL_allIter_allSubj, 1, 'omitnan')); % [nSubj x nCond]

        % (2) Bootstrap group mean + CI for bars
        ave_allBoot = nan(nBoot, nCond);
        % Pregenerate subj indices
        rng(seedBoot, 'twister');
        indRand_allBoot = randi(nSubj, [nBoot, nSubj], 'uint16');
        parfor iBoot = 1:nBoot
            indRandBoot = double(indRand_allBoot(iBoot, :));
            ave_allBoot(iBoot,:) = mean(data_med_allSubj(indRandBoot,:), 1, 'omitnan');
        end

        % getCI should return median/mean + CI bounds; here we just want point + lb/ub
        [data_ave, data_lb, data_ub, data_sem_neg, data_sem_pos] = getCI(ave_allBoot, 1, 1, CI68);

        % (3) Bootstrap CIs for mean differences for specific pairs
        pairs_all = nchoosek(1:nCond, 2);
        nPairs = size(pairs_all, 1);

        % Compute bootstrap for ALL pairs once (so you can reuse)
        ave_allBoot = nan(nBoot, nCond);
        % Pregenerate subj indices
        rng(seedBoot, 'twister');
        indRand_allBoot = randi(nSubj, [nBoot, nSubj], 'uint16');
        parfor iBoot = 1:nBoot
            indRandBoot = double(indRand_allBoot(iBoot, :));
            ave_allBoot(iBoot,:) = mean(data_med_allSubj(indRandBoot,:), 1, 'omitnan');
        end

        diffCond_allBoot = nan(nBoot, nPairs);
        % Pregenerate subj indices
        rng(seedBoot, 'twister');
        indRand_allBoot = randi(nSubj, [nBoot, nSubj], 'uint16');
        parfor iBoot = 1:nBoot
            indRandBoot = double(indRand_allBoot(iBoot, :));
            Xb = data_med_allSubj(indRandBoot,:); % [nSubj x nCond]
            for iPair = 1:nPairs
                iA = pairs_all(iPair,1);
                iB = pairs_all(iPair,2);
                diffCond_allBoot(iBoot, iPair) = mean(Xb(:,iA) - Xb(:,iB), 'omitnan'); % M1-M2
            end % iPair
        end % iBoot

        [diffCond_med, diffCond_lb, diffCond_ub, diffCond_sem_neg, diffCond_sem_pos] = getCI(diffCond_allBoot, 1, 1, CI95); % use 95% for reporting

        % Set limits
        yl = ylim;
        yMin = yl(1);
        yMax = yl(2);
        yRange = yMax - yMin;

        % Place brackets above the tallest bar+CI
        topData = y_ticks(end);
        if isnan(topData), topData = yMax; end
        yBase = yMin + 0.2*yRange; % lower the scalar, lower yBase
        yStep = .1; % vertical spacing between brackets

        % Helper: find index in "pairs" for a given (iA,iB)
        getPairIdx = @(iA,iB) find(pairs_all(:,1)==min(iA,iB) & pairs_all(:,2)==max(iA,iB), 1, 'first');

        % [Plot] planned comparison brackets + CI of mean difference atmidpoint
        if ~isnan(pairs_bracket)
            for iPair = 1:size(pairs_bracket,1)

                iA = pairs_bracket(iPair,1);
                iB = pairs_bracket(iPair,2);

                y = yBase + (iPair-1)*yStep;

                % --- bracket line ---
                plot([iA iB], [y y], 'k-', 'LineWidth', wd, 'HandleVisibility','off');
                capHeight = 0.015 * yRange;
                plot([iA iA], [y-capHeight, y], 'k-', 'LineWidth', wd, 'HandleVisibility','off');
                plot([iB iB], [y-capHeight, y], 'k-', 'LineWidth', wd, 'HandleVisibility','off');

                % --- text label (left-aligned, above the bracket line) ---
                % "pperm_allPairs(1)" is hard-coded!! 1=comparing Full vs. NoMul
                str_delta = sprintf('$\\Delta_{CI95}=[%+.2f, %+.2f]$\n$\\mathit{p}=%.3f$', diffCond_lb(iPair), diffCond_ub(iPair), pperm_allPairs(1));

                xText = iA; % left end of bracket
                yText = y + 0.02*yRange; % a bit above the bracket line (tune 0.02)

                text(xText, yText, str_delta, ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', sz_text, ...
                    'Color', 'k', ...
                    'Interpreter', 'latex', ...
                    'HandleVisibility', 'off');
            end % iPair
        end
        % % ======== Inset ========
        % % Compare raw nLL_med (not delta nLL) for B1 vs B2 using fxn_drawBars
        % ax_inset = axes('Position', [0.08, 0.62, 0.33, 0.28]);
        % nLL_B12_allIter_allSubj = nan(1, nSubj, 2);
        % nLL_B12_allIter_allSubj(1,:,:) = nLL_med(1:2, :)';
        %
        % c1 = colors_comb(iLocSingle, :);
        % c2 = max(c1 - 0.35, 0);
        % colors_B12 = [c1; c2];
        %
        % flag_plotIDVD_inset = 1;
        % flag_plotDiff_inset = 0;
        % str_title_inset = sprintf('B1 vs B2 nLL_{med}');
        % fxn_drawBars(nLL_B12_allIter_allSubj, nan, colors_B12, namesModelB(iModelB_selected(1:2)), nan, nan, ...
        %     flag_plotIDVD_inset, flag_plotDiff_inset, str_title_inset, nan, 1, markers_allSubj);
        % ylabel('nLL', 'FontSize', 12);
        % title('B1 vs B2 (nLL_{med})', 'FontSize', 10);
        % set(ax_inset, 'LineWidth', 1.5);
        % %========================

        % xTL
        xticks(1:size(ave_allBoot, 2))
        xticklabels(namesModelB(iModelB_selected))
        xtickangle(45)

        % set axis font size
        ax = gca;
        ax.FontSize = sz_label * 0.8;
        ax.LineWidth = wd;

        title(str_title, 'FontSize', 10);

        saveas(figMain, fullfile(nameFolder_Fig_NOM_nLL, sprintf('n%d_L%d_ModelB%s.png', nSubj, iLocSingle, strjoin(string(iModelB_selected), ''))));
        % close(figMain); % do NOT close here, as we need to copy the axis geometry for the box-code figure; close at the end of the loop after copying geometry.

        % ======== Plot the box-code figure (only once) ========
        if iLocSingle==1

            axMain  = gca;

            % Copy main-axis geometry as normalized proportions
            % This preserves the left margin used for y-ticks/y-labels.
            oldUnitsFig = get(figMain, 'Units');
            oldUnitsAx  = get(axMain,  'Units');

            set(figMain, 'Units', 'pixels');
            set(axMain,  'Units', 'normalized');

            figMainPos = get(figMain, 'Position');   % only use width, not screen location
            axMainPos  = get(axMain,  'Position');   % [left bottom width height], normalized

            xLim_main = get(axMain, 'XLim');

            % Restore units
            set(figMain, 'Units', oldUnitsFig);
            set(axMain,  'Units', oldUnitsAx);

            % Box-code matrix
            % Rows: Nmul, Nadd, Nshared
            % Canonical columns (B1..B7): Full, NoMul, NoAdd, NoShared, OnlyMul, OnlyAdd, OnlyShared
            modelBox_all = [
                1 0 1 1 1 0 0;   % Nmul
                1 1 0 1 0 1 0;   % Nadd
                1 1 1 0 0 0 1    % Nshared
                ];
            modelBox = modelBox_all(:, iModelB_selected);

            [nRows, nCols] = size(modelBox);

            % Independent canvas
            figBoxW = figMainPos(3);   % same canvas width as main figure
            figBoxH = 180;             % independent, short canvas

            figBox = figure('Color', 'w', 'Units', 'pixels', 'Position', [100, 100, figBoxW, figBoxH], 'PaperPositionMode', 'auto');

            % Own axes, but with same left/right plotting geometry as main axes
            boxBottom = 0.20;
            boxHeight = 0.65;

            axBox = axes(figBox, 'Units', 'normalized', 'Position', [axMainPos(1), boxBottom, axMainPos(3), boxHeight]);

            hold(axBox, 'on');

            set(axBox, 'XLim', xLim_main, 'YLim', [0.5, nRows + 0.5], 'XTick', 1:nCols, 'YTick', [], 'Visible', 'off');

            markerSize = 600;
            edgeWidth  = 2.5;

            for iCol = 1:nCols
                for iRow = 1:nRows

                    x = iCol;
                    y = nRows - iRow + 1;  % first row appears at top

                    if modelBox(iRow, iCol) == 1
                        faceColor = 'k';
                    else
                        faceColor = 'w';
                    end

                    scatter(axBox, x, y, markerSize, 's', 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', 'k', 'LineWidth', edgeWidth);
                end % iRow
            end % iCol

            % Save the full independent canvas
            saveas(figBox, fullfile(nameFolder_Fig_NOM_nLL, sprintf('BoxCode_B%s.png', strjoin(string(iModelB_selected), ''))));
            close(figBox);
        end % if

        close(figMain) % MUST be here!! as the axis inherited for the box-code figure is from the main figure, so we can't close the main figure before copying the axis geometry.

    end % iLocSingle

    clear nLL_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 18. [NOM] Compare NOM parameters across locations
    % Load data
    load(nameFolder_Output_SaveCompile, 'params_allCond')
    clc; fprintf('\n%s: Fig 18/23: Compare NOM params across locations STARTED.\n', string(datetime('now')))

    % Define folder for saving GoF figures
    nameFolder_Fig_NOM_params = fullfile(nameFolder_Fig_NOM_Trialwise, 'NOMparams');
    if isempty(dir(nameFolder_Fig_NOM_params)), mkdir(nameFolder_Fig_NOM_params); end

    flag_plotIDVD = 1;
    flag_plotDiff = 1;
    nBars = 2;
    sz_fig = [nBars * 200, 400+nchoosek(nBars, 2)*30];
    sz_text = 22;
    wd = 2;

    for iModelB_NOMplot = 2; %iModelB_plot_all % Only plot NoMul

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
                    switch iModelB_NOMplot
                        case 1 % FullModel
                            switch iParam
                                case 1, y_ticks = linspace(NOMp2_lb, 1.2, 5);
                                case 2, y_ticks = linspace(NOMp3_lb, 3, 5);
                                case 3, y_ticks = linspace(NOMp3_lb, 3, 5);
                                case 4, y_ticks = linspace(NOMp3_lb, 40, 5);
                            end
                        case 2 % NoMul
                            switch iParam
                                case 1, y_ticks = linspace(0, 2.8, 5);
                                case 2, y_ticks = linspace(0, 2.8, 5);
                                case 3, y_ticks = linspace(NOMp3_lb, 40, 5);
                            end
                    end
                    y_ticklabels = round(y_ticks, 2);

                    str_title = sprintf('n=%d nIter=%d L%s [A%dB%d] %s', nSubj, nIterxJob, strjoin(string(iLocPair_all), ''), iModelA_plot, iModelB_NOMplot, namesModelBparams{iModelB_NOMplot}{iParam});

                    %------------------------------%
                    fxn_drawBars(NOMp_allIter_allSubj, ref, colors_comb(iLocPair_all, :), x_ticks, y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj, sz_text, wd, flag_plotPurpose);
                    %------------------------------%
                    fig = gcf;
                    ylabel(sprintf('%s', namesModelBparams{iModelB_NOMplot}{iParam}))
                    % Save the figure
                    saveAndCloseFigure(fig, sprintf('%s/n%d_A%dB%d_L%s_NOMp%d.png', nameFolder_Fig_NOM_params, nSubj, iModelA_plot, iModelB_NOMplot, strjoin(string(iLocPair_all), ''), iParam));
                end % iParam
            end % iGroup
        end % ii
    end % iModelB_NOMplot
    clear params_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 19. [NOM] Corr1: pA and NOM params
    clc; fprintf('\n%s: Fig 19/23: Corr pA vs NOM params STARTED.\n', string(datetime('now')))

    % Load data
    load(nameFolder_Output_SaveCompile, 'params_allCond', 'pA_allSubj')
    nameVarX = 'pA';
    nameVarY = 'NOMparams';

    sz_label = 55;
    sz_labelOffset = 0.02; % normalized units to move labels away from axes; tune as needed
    sz_axOffset = 0.05; % normalized units to move axes away from figure edges; tune as needed

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

    nameFolder_Outputs_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Outputs_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Outputs_NOM_corr)); mkdir(nameFolder_Outputs_NOM_corr), end

    flag_zeroMean = 0;
    flag_plotIdvdCI = 0;
    flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles

    for iModelB_NOMplot = 2%iModelB_plot_all % Just plot NoMul!!
        nNOMparams = length(namesModelBparams{iModelB_NOMplot});

        for iSet = 1:numel(iLocSingle_allSets)

            iLocCorr_all = iLocSingle_allSets{iSet};
            fprintf('\n - B%d L%s', iModelB_NOMplot, strjoin(string(iLocCorr_all), ''))

            % Obtain x-axis values
            X_allSubj = pA_allSubj(:, iLocCorr_all);
            x_ticks = linspace(.6, .8, 5); % EE
            X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);

            for iNOMparam = 1:nNOMparams
                fprintf(' NOMp%d ', iNOMparam);
                NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_NOMplot, iLocCorr_all, :, :, iNOMparam));

                switch iModelB_NOMplot
                    case 1, y_ticks_lb = [0, 0, 0, 0]; y_ticks_ub = [.8, 30, 30, 40];
                    case 2, y_ticks_lb = [0, 0, 0]; y_ticks_ub = [2.8, 2.8, 40]; %*
                    case 3, y_ticks_lb = [0, 0, 0]; y_ticks_ub = [.8, 30, 40];
                    case 4, y_ticks_lb = [0, 0, 0]; y_ticks_ub = [.8, 30, 40];
                    case 5, y_ticks_lb = [0, 0]; y_ticks_ub = [.8, 40];
                    case 6, y_ticks_lb = [0, 0]; y_ticks_ub = [30, 40];
                    case 7, y_ticks_lb = [0, 0]; y_ticks_ub = [30, 40];
                    otherwise, error('Unsupported iModelB_NOMplot=%d', iModelB_NOMplot);
                end

                y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

                x_ticklabels = x_ticks;
                y_ticklabels = y_ticks;

                nameVarY_figTitle = namesModelBparams{iModelB_NOMplot}{iNOMparam};
                nameVarY_fileTitle = nameVarY_figTitle;

                str_title = sprintf('n=%d, nIter=%d B%d %s vs. %s [L%s]', nSubj, nIterxJob, iModelB_NOMplot, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

                %----------------------------%
                [PearsonR, SpearmanRho] = fxn_drawCorr(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), flag_UseRUseRho, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, nIterxJob, flag_plotPurpose);
                %----------------------------%
                fig = gcf;
                % Save CI for automatic CI range calcuation in CorrAsym
                save(sprintf('%s/n%d_B%d_L%s_NOMp%d.mat', nameFolder_Outputs_NOM_corr, nSubj, iModelB_NOMplot, strjoin(string(iLocCorr_all), ''), iNOMparam), 'PearsonR', 'SpearmanRho')

                xlabel('Consistency rate', 'fontsize', sz_label)
                ylabel(nameVarY_figTitle, 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                drawnow; % Ensure text/tick extents are up to date before reading TightInset
                ti = ax.TightInset; % [left bottom right top] padding needed for labels/ticks
                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset; % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset; % move label left

                newPos = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];
                if all(isfinite(newPos)) && newPos(3) >= 0.55 && newPos(4) >= 0.55
                    ax.Position = newPos;
                end

                saveAndCloseFigure(fig, sprintf('%s/n%d_B%d_L%s_NOMp%d.png', nameFolder_Fig_NOM_corr, nSubj, iModelB_NOMplot, strjoin(string(iLocCorr_all), ''), iNOMparam))
            end % iNOMparam
        end % iSet
    end % iModelB_NOMplot
    clear params_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    % %% [NOM] CorrAsym1: pA and NOMparams: corr between extents of EE/HVA/VMA
    % clc, fprintf('\n\n 21/24 Plotting STARTED......\n\n')
    %
    % % Load data
    % load(nameFolder_Data_SaveCompile, 'params_allCond', 'pA_allSubj')
    %
    % nameVarX = 'pA';
    % nameVarY = 'NOMparams';
    %
    % nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    % if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end
    %
    % sz_label = 60;
    % sz_labelOffset = 0.05;
    % sz_axOffset = 0.05; % extra breathing room
    %
    % for iModelB_NOMplot = iModelB_plot_all
    %     nNOMparams = length(namesModelBparams{iModelB_NOMplot});
    %
    %     for iGroup = 1:nGroups
    %         iLocPair_all = iLocGroups_all{iGroup};
    %         fprintf(' - L%s\n', strjoin(string(iLocPair_all), ''))
    %
    %         switch iLocPair_all(1)
    %             case 1, nameAsymX = 'EE'; nameAsymY = 'EE';
    %             case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
    %             case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
    %         end
    %
    %         % X-axis (pA)
    %         asymX_allSubj = (pA_allSubj(:, iLocPair_all(1))-pA_allSubj(:, iLocPair_all(2)))./(pA_allSubj(:, iLocPair_all(1))+pA_allSubj(:, iLocPair_all(2)));
    %         asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';
    %         % switch iLocPair_all(1)
    %         % case 1, x_ticks = linspace(0, 20, 5); % EE
    %         % case 6, x_ticks = linspace(-5, 15, 5); % HVA
    %         % case 5, x_ticks = linspace(0, 16, 5); % VMA (extent is smaller)
    %         % end
    %         x_ticks = linspace(-6, 10, 5);
    %
    %         for iNOMparam = 1:nNOMparams
    %             NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_NOMplot, iLocPair_all, :, :, iNOMparam));
    %
    %             asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));
    %
    %             % switch iLocPair_all(1)
    %             % case 1 % EE
    %             % y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
    %             % case 6 % HVA
    %             % y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
    %             % case 5 % VMA
    %             % y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
    %             % end
    %             switch iModelB_NOMplot
    %                 case 1, y_ticks_lb = -[60, 40, 20]; y_ticks_ub = [40, 40, 60];
    %                 case 3, y_ticks_lb = -[40, 20]; y_ticks_ub = [40, 60];
    %             end
    %             y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);
    %
    %             x_ticklabels = nan;
    %             y_ticklabels = nan;
    %
    %             nameVarY_figTitle = namesModelBparams{iModelB_NOMplot}{iNOMparam};
    %             nameVarY_fileTitle = nameVarY_figTitle;
    %
    %             str_title = sprintf('n=%d, nIter=%d B%d %s (%s) vs. %s (%s)', nSubj, nIterxJob, iModelB_NOMplot, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);
    %
    %             %----------------------------%
    %             fxn_drawCorrAsym(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, flag_plotPurpose)
    %             %----------------------------%
    %
    %             xlabel(sprintf('%s of %s (%%)', nameAsymX, nameVarX), 'fontsize', sz_label)
    %             ylabel(sprintf('%s of %s (%%)', nameAsymX, nameVarY_figTitle), 'fontsize', sz_label)
    %
    %             % Adjust distance between components
    %             ax = gca;
    %             ti = ax.TightInset; % [left bottom right top] padding needed for labels/ticks
    %
    %             ax.XLabel.Units = 'normalized';
    %             ax.YLabel.Units = 'normalized';
    %
    %             ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset; % move label down
    %             ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset; % move label left
    %             ax.Position = [ ...
    %                 ti(1) + sz_axOffset, ...
    %                 ti(2) + sz_axOffset, ...
    %                 1 - ti(1) - ti(3) - 2*sz_axOffset, ...
    %                 1 - ti(2) - ti(4) - 2*sz_axOffset];
    %
    %             saveas(gcf, sprintf('%s/n%d_B%d_L%s_NOMp%d.png', nameFolder_Fig_NOM_CorrAsym, nSubj, iModelB_NOMplot, strjoin(string(iLocPair_all), ''), iNOMparam))
    %             close(gcf)
    %
    %         end % iNOMparam
    %     end % iGroup
    % end % iModelB_NOMplot
    % clear params_allCond
    % fprintf('\n\n Plotting DONE\n\n')

    %% 20. [NOM] Corr2: CS and NOM params
    clc; fprintf('\n%s: Fig 20/23: Corr CS vs NOM params STARTED.\n', string(datetime('now')))

    % Load data
    load(nameFolder_Output_SaveCompile, 'params_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'NOMparams';

    sz_label = 55;
    sz_labelOffset = 0.02;
    sz_axOffset = 0.05; % extra breathing room

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

    nameFolder_Outputs_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Outputs_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Outputs_NOM_corr)); mkdir(nameFolder_Outputs_NOM_corr), end

    flag_zeroMean = 0;
    flag_plotIdvdCI = 0;
    flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles

    for iModelB_NOMplot = 2;%iModelB_plot_all % Only plot NoMul!
        nNOMparams = length(namesModelBparams{iModelB_NOMplot});

        for iSet = 1:numel(iLocSingle_allSets)

            iLocCorr_all = iLocSingle_allSets{iSet};
            fprintf(' - B%d L%s\n', iModelB_NOMplot, strjoin(string(iLocCorr_all), ''))

            % Obtain CS (x-axis)
            X_allSubj = CS_allSubj(:, iLocCorr_all);
            x_ticks = linspace(1.5, 3.5, 5); % EE
            X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);
            % X_med_allSubj = X_allSubj;

            for iNOMparam = 1:nNOMparams
                NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_NOMplot, iLocCorr_all, :, :, iNOMparam));

                Y_med_allSubj = getCI(NOMp_allIter_allSubj, 1, 3)'; % rotate to match the format needed by fxn_drawCorr

                switch iModelB_NOMplot
                    case 1, y_ticks_lb = [0, 0, 0, 0]; y_ticks_ub = [.6, 32, 32, 40];
                    case 2, y_ticks_lb = [0, 0, 0]; y_ticks_ub = [2.8, 2.8, 40]; %NoMul*
                    case 3, y_ticks_lb = [0, 0, 0]; y_ticks_ub = [.6, 32, 40];
                    case 4, y_ticks_lb = [0, 0, 0]; y_ticks_ub = [.6, 32, 40];
                    case 5, y_ticks_lb = [0, 0]; y_ticks_ub = [.6, 40];
                    case 6, y_ticks_lb = [0, 0]; y_ticks_ub = [32, 40];
                    case 7, y_ticks_lb = [0, 0]; y_ticks_ub = [32, 40];
                    otherwise, error('Unsupported iModelB_NOMplot=%d', iModelB_NOMplot);
                end

                y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

                x_ticklabels = x_ticks;
                y_ticklabels = y_ticks;

                nameVarY_figTitle = namesModelBparams{iModelB_NOMplot}{iNOMparam};
                nameVarY_fileTitle = nameVarY_figTitle;

                str_title = sprintf('n=%d, nIter=%d B%d %s vs. %s [L%s]', nSubj, nIterxJob, iModelB_NOMplot, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

                %----------------------------%
                [PearsonR, SpearmanRho] = fxn_drawCorr(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), flag_UseRUseRho, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, nIterxJob, flag_plotPurpose);
                %----------------------------%
                fig = gcf;
                % Save CI for automatic CI range calcuation in CorrAsym
                save(sprintf('%s/n%d_B%d_L%s_NOMp%d.mat', nameFolder_Outputs_NOM_corr, nSubj, iModelB_NOMplot, strjoin(string(iLocCorr_all), ''), iNOMparam), 'PearsonR', 'SpearmanRho')

                xlabel('Contrast sensitivity (1/contrast)', 'fontsize', sz_label)
                ylabel(nameVarY_figTitle, 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                drawnow; % Ensure text/tick extents are up to date before reading TightInset
                ti = ax.TightInset; % [left bottom right top] padding needed for labels/ticks
                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset; % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset; % move label left

                newPos = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];
                if all(isfinite(newPos)) && newPos(3) >= 0.55 && newPos(4) >= 0.55
                    ax.Position = newPos;
                end

                saveAndCloseFigure(fig, sprintf('%s/n%d_B%d_L%s_NOMp%d.png', nameFolder_Fig_NOM_corr, nSubj, iModelB_NOMplot, strjoin(string(iLocCorr_all), ''), iNOMparam))
            end % iNOMparam
        end % iSet
    end % iModelB_NOMplot
    clear params_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% 21. [NOM] CorrAsym2: CS and NOMparams: corr between extents of EE/HVA/VMA
    clc; fprintf('\n%s: Fig 21/23: CorrAsym CS vs NOM params STARTED.\n', string(datetime('now')))

    % Load data
    load(nameFolder_Output_SaveCompile, 'params_allCond', 'CS_allSubj')

    nameVarX = 'CS';
    nameVarY = 'NOMparams';

    nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

    nameFolder_Outputs_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Outputs_NOM_Trialwise, nameVarX, nameVarY);
    if isempty(dir(nameFolder_Outputs_NOM_corr)); mkdir(nameFolder_Outputs_NOM_corr), end

    sz_label = 60;
    sz_labelOffset = 0.05;
    sz_axOffset = 0.05; % extra breathing room
    flag_plotIdvdCI = 1;
    flag_plotUnikSymbol = 0;
    % iLocCorr_all = [6,5,3];

    for iModelB_NOMplot = 2; %iModelB_plot_all % Only plot NoMul!
        nNOMparams = length(namesModelBparams{iModelB_NOMplot});

        for iGroup = 1:nGroups
            iLocPair_all = iLocGroups_all{iGroup};
            fprintf(' - B%d L%s\n', iModelB_NOMplot, strjoin(string(iLocPair_all), ''))

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
            x_ticks = linspace(-3, 17, 5);
            asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';

            for iNOMparam = 1:nNOMparams
                NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_NOMplot, iLocPair_all, :, :, iNOMparam));

                asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

                switch iModelB_NOMplot
                    case 1, y_ticks_lb = -[90, 50, 50, 40]; y_ticks_ub = [90, 50, 50, 40];
                    case 2, y_ticks_lb = -[50, 50, 40]; y_ticks_ub = [50, 50, 40]; %* NoMul
                    case 3, y_ticks_lb = -[90, 50, 40]; y_ticks_ub = [90, 50, 40];
                    case 4, y_ticks_lb = -[90, 50, 40]; y_ticks_ub = [90, 50, 40];
                    case 5, y_ticks_lb = -[90, 40]; y_ticks_ub = [90, 40];
                    case 6, y_ticks_lb = -[50, 40]; y_ticks_ub = [50, 40];
                    case 7, y_ticks_lb = -[50, 40]; y_ticks_ub = [50, 40];
                    otherwise, error('Unsupported iModelB_NOMplot=%d', iModelB_NOMplot);
                end
                y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

                x_ticklabels = nan;
                y_ticklabels = nan;

                nameVarY_figTitle = namesModelBparams{iModelB_NOMplot}{iNOMparam};
                nameVarY_fileTitle = nameVarY_figTitle;

                str_title = sprintf('n=%d, nIter=%d B%d %s (%s) vs. %s (%s)', nSubj, nIterxJob, iModelB_NOMplot, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

                % Load corr analysis
                load(sprintf('%s/n%d_B%d_L%s_NOMp%d.mat', nameFolder_Outputs_NOM_corr, nSubj, iModelB_NOMplot, strjoin(string(iLocCorr_all), ''), iNOMparam), 'PearsonR', 'SpearmanRho')
                switch flag_UseRUseRho
                    case 'useR'
                        flag_CIrange = PearsonR(4) < 0.05; % use permutation p-value to determine whether one-tailed corr should be conducted
                        % flag_CIrange = PearsonR(2) * PearsonR(3)>0; % flag_CIrange=1 if r excludes 0, so one-tailed corr should be conducted, so CI range is 90%
                    case 'useRho'
                        flag_CIrange = SpearmanRho(4) < 0.05; % use permutation p-value to determine whether one-tailed corr should be conducted
                        % flag_CIrange = SpearmanRho(2) * SpearmanRho(3)>0; % flag_CIrange=1 if r excludes 0, so one-tailed corr should be conducted, so CI range is 95%
                end

                %----------------------------%
                fxn_drawCorrAsym(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, colors_comb(iLocPair_all, :), flag_UseRUseRho, flag_CIrange, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj, flag_plotPurpose)
                %----------------------------%
                fig = gcf;

                xlabel(sprintf('%s of contrast sensitivity (%%)', nameAsymX), 'FontSize', sz_label);
                ylabel(sprintf('%s of %s (%%)', nameAsymY, nameVarY_figTitle), 'fontsize', sz_label)

                % Adjust distance between components
                ax = gca;
                drawnow; % Ensure text/tick extents are up to date before reading TightInset
                ti = ax.TightInset; % [left bottom right top] padding needed for labels/ticks

                ax.XLabel.Units = 'normalized';
                ax.YLabel.Units = 'normalized';

                ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset; % move label down
                ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset; % move label left
                newPos = [ ...
                    ti(1) + sz_axOffset, ...
                    ti(2) + sz_axOffset, ...
                    1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                    1 - ti(2) - ti(4) - 2*sz_axOffset];
                if all(isfinite(newPos)) && newPos(3) >= 0.55 && newPos(4) >= 0.55
                    ax.Position = newPos;
                end

                saveAndCloseFigure(fig, sprintf('%s/n%d_B%d_L%s_NOMp%d.png', nameFolder_Fig_NOM_CorrAsym, nSubj, iModelB_NOMplot, strjoin(string(iLocPair_all), ''), iNOMparam))

            end % iNOMparam
        end % iGroup
    end % iModelB_NOMplot
    clear params_allCond
    fprintf('\n%s: DONE.\n', string(datetime('now')))

    %% [NOM] CompAsym2: CS and NOMparams: compare binned estimates
    % clc; fprintf('\n%s: Fig 22/23: CompAsym CS vs NOM params STARTED.\n', string(datetime('now')))
    %
    % % Load data
    % load(nameFolder_Data_SaveCompile, 'params_allCond', 'CS_allSubj')
    %
    % nameVarX = 'CS';
    % nameVarY = 'NOMparams';
    %
    % nameFolder_Fig_NOM_CompAsym = sprintf('%s/CompAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
    % if isempty(dir(nameFolder_Fig_NOM_CompAsym)); mkdir(nameFolder_Fig_NOM_CompAsym), end
    %
    % sz_label = 15;
    % sz_labelOffset = 0.05;
    % sz_axOffset = 0.05; % extra breathing room
    %
    % for iModelB_NOMplot = iModelB_plot_all
    %     nNOMparams = length(namesModelBparams{iModelB_NOMplot});
    %
    %     for iGroup = 1:nGroups
    %         iLocPair_all = iLocGroups_all{iGroup};
    %         fprintf(' - L%s\n', strjoin(string(iLocPair_all), ''))
    %
    %         switch iLocPair_all(1)
    %             case 1, nameAsymX = 'Ecc. effect'; nameAsymY = 'Ecc. effect';
    %             case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
    %             case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
    %         end
    %
    %         asymX_allSubj = (CS_allSubj(:, iLocPair_all(1))-CS_allSubj(:, iLocPair_all(2)))./(CS_allSubj(:, iLocPair_all(1))+CS_allSubj(:, iLocPair_all(2)));
    %         switch iLocPair_all(1)
    %             case 1, x_ticks = linspace(0, 20, 5); % EE
    %             case 6, x_ticks = linspace(-5, 15, 5); % HVA
    %             case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
    %         end
    %         asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';
    %
    %         for iNOMparam = 1:nNOMparams
    %             NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_NOMplot, iLocPair_all, :, :, iNOMparam));
    %
    %             asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));
    %
    %             switch iModelB_NOMplot
    %                 case 1, y_ticks_lb = -[30, 20, 20, 20]; y_ticks_ub = [50, 40, 60, 60];
    %                 case 2, y_ticks_lb = -[20, 20, 20]; y_ticks_ub = [40, 60, 60];
    %                 case 3, y_ticks_lb = -[30, 20, 20]; y_ticks_ub = [50, 60, 60];
    %                 case 4, y_ticks_lb = -[30, 20, 20]; y_ticks_ub = [50, 40, 60];
    %                 case 5, y_ticks_lb = -[30, 20]; y_ticks_ub = [50, 60];
    %                 case 6, y_ticks_lb = -[20, 20]; y_ticks_ub = [40, 60];
    %                 case 7, y_ticks_lb = -[20, 20]; y_ticks_ub = [40, 60];
    %                 otherwise, error('Unsupported iModelB_NOMplot=%d', iModelB_NOMplot);
    %             end
    %             y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);
    %
    %             x_ticklabels = nan;
    %             y_ticklabels = nan;
    %
    %             nameVarY_figTitle = namesModelBparams{iModelB_NOMplot}{iNOMparam};
    %             nameVarY_fileTitle = nameVarY_figTitle;
    %
    %             str_title = sprintf('n=%d, nIter=%d B%d %s (%s) vs. %s (%s)', nSubj, nIterxJob, iModelB_NOMplot, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);
    %
    %             % remove later!!s
    %             asymX_allIter_allSubj_ = repmat(asymX_allIter_allSubj', 5, 1);
    %             asymX_allIter_allSubj_ = asymX_allIter_allSubj_+randn(size(asymX_allIter_allSubj_))*mean(asymX_allIter_allSubj_(:))/20;
    %             asymY_allIter_allSubj_ = repmat(asymY_allIter_allSubj, 5, 1);
    %             asymY_allIter_allSubj_ = asymY_allIter_allSubj_+randn(size(asymY_allIter_allSubj_))*mean(asymY_allIter_allSubj_(:))/20;
    %
    %             str_ylabel = sprintf('\\Delta %s (%%)', nameVarY_figTitle);
    %
    %             %----------------------------%
    %             fxn_compAsym(asymX_allIter_allSubj_*100, asymY_allIter_allSubj_*100, nBinsCompAsym, y_ticks, sz_fig, str_title, str_ylabel)
    %             %----------------------------%
    %
    %             % xlabel(sprintf('\\Delta contrast sensitivity (%%)'), 'FontSize', sz_label);
    %             % ylabel(, 'fontsize', sz_label)
    %
    %             % Adjust distance between components
    %             % ax = gca;
    %             % ti = ax.TightInset; % [left bottom right top] padding needed for labels/ticks
    %             %
    %             % ax.XLabel.Units = 'normalized';
    %             % ax.YLabel.Units = 'normalized';
    %             %
    %             % ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset; % move label down
    %             % ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset; % move label left
    %             % ax.Position = [ ...
    %             % ti(1) + sz_axOffset, ...
    %             % ti(2) + sz_axOffset, ...
    %             % 1 - ti(1) - ti(3) - 2*sz_axOffset, ...
    %             % 1 - ti(2) - ti(4) - 2*sz_axOffset];
    %
    %             saveas(gcf, sprintf('%s/n%d_B%d_L%s_NOMp%d.png', nameFolder_Fig_NOM_CompAsym, nSubj, iModelB_NOMplot, strjoin(string(iLocPair_all), ''), iNOMparam))
    %             close(gcf)
    %
    %         end % iNOMparam
    %     end % iGroup
    % end % iModelB_NOMplot
    % clear params_allCond
    % fprintf('\n%s: DONE.\n', string(datetime('now')))


end% iRun

%%
function saveAndCloseFigure(figHandle, outName)
if ~isgraphics(figHandle, 'figure')
    return
end

drawnow;
set(figHandle, 'PaperPositionMode', 'auto');

% Recover from occasional collapsed single-axis layout in batch/invisible mode.
axList = findall(figHandle, 'Type', 'axes');
if numel(axList) == 1
    ax = axList(1);
    oldUnits = ax.Units;
    ax.Units = 'normalized';

    % Keep large labels/ticks visible in exported images by reserving
    % margin from the current tight inset rather than using a fixed box.
    drawnow;
    ti = ax.TightInset; % [left bottom right top]

    % Figure 7 (TuningCs) often needs extra left margin for long y-labels.
    isTuningCsFig = ~isempty(strfind(outName, [filesep 'TuningCs' filesep]));
    if isTuningCsFig
        minLeft = 0.30;
        minBottom = 0.20;
        minTop = 0.12;
    else
        minLeft = 0.20;
        minBottom = 0.18;
        minTop = 0.04;
    end

    left = max(minLeft, ti(1) + 0.05);
    bottom = max(minBottom, ti(2) + 0.05);
    right = max(0.04, ti(3) + 0.02);
    top = max(minTop, ti(4) + 0.02);

    w = 1 - left - right;
    h = 1 - bottom - top;
    if isfinite(w) && isfinite(h) && w > 0.40 && h > 0.40
        ax.Position = [left, bottom, w, h];
    else
        ax.Position = [0.18, 0.16, 0.74, 0.74];
    end

    ax.Units = oldUnits;
end

try
    exportgraphics(figHandle, outName, 'Resolution', 300);
catch
    saveas(figHandle, outName);
end
close(figHandle);
end

function s = formatMixedNumber(v)
av = abs(v);
if isnan(v)
    s = 'NaN';
elseif av >= 0.1
    s = sprintf('%.2f', v);
elseif av >= 0.01
    s = sprintf('%.3f', v);
elseif av == 0
    s = '0';
else
    s = sprintf('%.1e', v);
end
end
