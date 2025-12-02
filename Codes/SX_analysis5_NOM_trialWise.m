% SX_analysis5_NOM_Trialwise.m
% Author: Shutian Xue
% Purpose: This script compiles, analyzes and plots boostrapped data for each observer saved Data_NOM_Trialwise

clc, clear all, close all

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
iModelB_sim_all = [1,2,3]; % see SX_RC1_setting for namesModelB
nBoot = 100; % Number of bootstraps (100 for raw IVs, 99 for transformed IVs)
iLocComb_all = [1,8, 6,7, 5,3]; % combination of locations
iLocGroups_all = {[1,8], [6,7], [5,3]} ; nGroups = length(iLocGroups_all);
iFamily_ORI = 1; %1=scaled Gaussian
iFamily_SF = 2; % 2=log parabola
iFamily_perF = [iFamily_ORI, iFamily_SF];
namesMetrics = {'pYES', 'pA'}; nMetrics = length(namesMetrics);
iDataset_plot = 2; %1=Training set; 2=Full set

iModelA_plot = 1; % just plot the RC-derived
iModelB_plot = 1; % 1=full model; 2=No Rho; 3=No induced noise

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
    % 
    % % n=13: no AS, CS
    % subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'HL', 'FH', 'HA',  'DT', 'DU', 'RC', 'SR'}; nRows_subj = 3; nCols_subj = 5;
    % nblocks_allSubj = [200, 240, 210, 240, 220, 220, 220, 205, 210, 205, 205, 205, 195];

    assert(nRows_subj * nCols_subj >= length(subjList), 'ALERT: Not enough subplots for the number of subjects!');

else
    % Find folders with names starting with 'IO' in Data/Data_NOM_Trialwise_2929
    % Filter for directories whose names start with 'IO'
    nameDir_IO_all = dir(sprintf('%s/IO*', nameFolder_Data_NOM_Trialwise));
    nIOs = length(nameDir_IO_all);
    subjList = cell(1, nIOs);
    nRows_subj=4; nCols_subj=7;
    for iIO = 1:nIOs
        subjList{iIO} = nameDir_IO_all(iIO).name; % Store the names of IO folders
    end
    assert(nRows_subj * nCols_subj >= nIOs, 'Not enough subplots for the number of subjects!');
    iLocComb_all = 1; % For testing with IO data only
end

nsubj = length(subjList);

% Names of the performance metrics being analyzed
nModelsA = length(iModelA_all);
nModelsB = length(iModelB_all);
nLocComb8 = length(iLocComb_all);
nBins = 10; % Number of bins for the analysis
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
fprintf(' - Information Criterion to plot: %s\n', namesIC{iIC_plot});

%% (Time-consuming!) Compile data from boostrapped data

% Preallocate arrays for storing results across all conditions
metrics_allCond = nan(nModelsA, nLocComb8, nsubj, nBoot, 2, 11); % 2: 1=Full; 2=Test; 11=number of metrics saved in OOD_xx_compIV.m; see "fxn_getMetrics"
% metrics_test = [dprime, criterion, [pC, pHit, pFA], nanmean(respC_test_rand), pYES, 1./mean(1./cst_full_rand_sel)];
template_train_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nORI, nSF);
template_full_allCond = template_train_allCond;
IV_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nBins);
nTrials_allCond = IV_allCond;
metric_data_allCond = nan(nModelsA, nModelsB, nLocComb8, nMetrics, nsubj, nBoot, nBins); % "nMetrics" is predefined in this script
metric_pred_allCond = metric_data_allCond;
nLL_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot);
params_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, 3); % Pre-allocate for max params
IC_nLL_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, 3); % AIC, AICc, BIC (3)
R2_NOM_allCond= nan(nModelsA, nModelsB, nLocComb8, nMetrics, nsubj, nBoot);

% Fitting tuning curves
nDatasets = 2; %1=template is derived from training set; 2=full dataset
sep_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nBoot, nDatasets); % 2 is for 1=training set and 2=full set

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
    fprintf('\n\nModel A%d: ', iModelA)

    for iModelB = iModelB_all %1:nModelsB
        fprintf('\n Model B%d: ', iModelB)
        nParamsB = length(namesModelBparams{iModelB});

        for iLocComb = iLocComb_all
            fprintf(' L%d ', iLocComb)

            for isubj = 1:nsubj

                subjName = subjList{isubj};

                if flag_subjIsHuman
                    % nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, nameIO);
                    nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocComb); % MUST be the same as OOD_NOM_Trialwise_compIV.m, Line 61
                else
                    % nameFolder_NOM_save = sprintf('%s/ORI%dSF%d/%s/L%d', nameFolder_NOM0, nORI, nSF, subjName, iLocComb);
                    nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, subjName); % MUST be the same as OOD_NOM_Trialwise_compIV.m, Line 68
                end

                % Load IVs and derived templates
                nameFile_compIV = sprintf('%s/n%d_A%d_compIV.mat', nameFolder_NOM_save, nBoot, iModelA);
                % if isempty(dir(nameFile_compIV)), error('ALERT: File %s does not exist!', nameFile_compIV); end
                if isempty(dir(nameFile_compIV))
                    fprintf(' "%sL%dA%d" ', subjName, iLocComb, iModelA)
                else
                    load(nameFile_compIV, '*_allBoot')
                end

                % Load predictions
                nameFile_fitNOM = sprintf('%s/n%d_A%dB%d.mat', nameFolder_NOM_save, nBoot, iModelA, iModelB);
                % if isempty(dir(nameFile_fitNOM)), error('ALERT: File %s does not exist!', nameFile_fitNOM); end
                if isempty(dir(nameFile_fitNOM))
                    fprintf(' "%sL%dA%dB%d" ', subjName, iLocComb, iModelA, iModelB)
                else
                    load(nameFile_fitNOM, '*_allBoot')
                end

                % Loaded "measured" metrics of the test set
                metrics_allCond(iModelA, iLocComb, isubj, :, :, :) = data_metrics_allBoot;

                % Pre-allocate temporary arrays for this subject and model
                IV_allBoot = nan(nBoot, nBins); % IV for each bin
                nTrials_allBoot = IV_allBoot; % Trial count for each bin
                metric_data_allBoot = nan(nMetrics, nBoot, nBins); % Measured data for each bin
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
                    for iMetric = 1:nMetrics
                        switch namesMetrics{iMetric}
                            case 'pYES'
                                metric_data_allBoot(iMetric, iBoot, :) = pred_metrics.pYES_data_allBins;
                                metric_pred_allBoot(iMetric, iBoot, :) = pred_metrics.pYES_pred_allBins;
                            case 'pA'
                                metric_data_allBoot(iMetric, iBoot, :) = pred_metrics.pA_data_allBins;
                                metric_pred_allBoot(iMetric, iBoot, :) = pred_metrics.pA_pred_allBins;
                            case 'pC'
                                metric_data_allBoot(iMetric, iBoot, :) = pred_metrics.pC_data_allBins;
                                metric_pred_allBoot(iMetric, iBoot, :) = pred_metrics.pC_pred_allBins;
                        end
                    end % iMetric
                    % Calculate information critertion based on nLL
                    nData = sum(nTrials_allBoot(iBoot, :));

                    % nLL_allBoot
                    nLL = nLL_allBoot(iBoot);
                    IC_nLL = [getAIC_nLL(nLL, nParamsB), getAICc_nLL(nLL, nParamsB, nData), getBIC_nLL(nLL, nParamsB, nData)];
                    IC_nLL_allCond(iModelA, iModelB, iLocComb, isubj, iBoot, :) = IC_nLL;

                    % Compute tuning characteristics based on the fitted parameters
                    for iDataset=1:2
                        % ORI
                        iFeature= 1;
                        margTuningC_ORI_allBoot(iBoot, iDataset, :) = fxn_getTuningC(axis_tuning{iFeature}, iFeature, iFamily_ORI, squeeze(margPred_ORI_allBoot(iBoot, iDataset, :)), squeeze(margParams_ORI_allBoot(iBoot, iDataset, :)));
                        % fxn_getTuningC(axis_tuning{iFeature}, iFeature, iFamily_ORI, squeeze(margPred_ORI_allBoot(iBoot, iDataset, :)), squeeze(margParams_ORI_allBoot(iBoot, iDataset, :)))

                        % SF
                        iFeature=2;
                        margTuningC_SF_allBoot(iBoot, iDataset, :) = fxn_getTuningC(axis_tuning{iFeature}, iFeature, iFamily_SF, squeeze(margPred_SF_allBoot(iBoot, iDataset, :)), squeeze(margParams_SF_allBoot(iBoot, iDataset, :)));
                    end
                end % end of iBoot

                % Calculate R2 for NOM prediction (for pYES and pA)
                R2_NOM_allBoot = nan(nMetrics, nBoot);
                for iMetric = 1:nMetrics
                    for iBoot = 1:nBoot
                        nData = squeeze(nTrials_allBoot(iBoot, :));
                        metric_data  = squeeze(metric_data_allBoot(iMetric, iBoot, :));  % ground truth
                        metric_pred = squeeze(metric_pred_allBoot(iMetric, iBoot, :)); % prediction
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
                        R2_NOM_allBoot(iMetric, iBoot) = R2;
                    end 
                end

                % Store results
                % Templates
                [template_train_med, ~, ~, template_train_sem] = getCI(template_train_allBoot, 1, 1);
                [template_full_med, ~, ~, template_full_sem] = getCI(template_full_allBoot, 1, 1);
                template_train_allCond(iModelA, iModelB, iLocComb, isubj, :, :) = template_train_med;
                template_full_allCond(iModelA, iModelB, iLocComb, isubj, :, :) = template_full_med;

                % Tuning functions: marg, predictions, estimated parameters and separability
                sep_allCond(iModelA, iModelB, iLocComb, isubj, :, :) = sep_allBoot;

                margORI_allCond(iModelA, iModelB, iLocComb, isubj, :, :, :) = margORI_allBoot; % directly loaded
                margR2_ORI_allCond(iModelA, iModelB, iLocComb, isubj, :, :) = margR2_ORI_allBoot; % directly loaded
                margPred_ORI_allCond(iModelA, iModelB, iLocComb, isubj, :, :, :) = margPred_ORI_allBoot; % directly loaded
                margParams_ORI_allCond(iModelA, iModelB, iLocComb, isubj, :, :, :) = margParams_ORI_allBoot; % directly loaded
                margTunC_ORI_allCond(iModelA, iModelB, iLocComb, isubj, :, :, :) = margTuningC_ORI_allBoot;

                margSF_allCond(iModelA, iModelB, iLocComb, isubj, :, :, :) = margSF_allBoot; % directly loaded
                margR2_SF_allCond(iModelA, iModelB, iLocComb, isubj, :, :) = margR2_SF_allBoot; % directly loaded
                margPred_SF_allCond(iModelA, iModelB, iLocComb, isubj, :, :, :) = margPred_SF_allBoot; % directly loaded
                margParams_SF_allCond(iModelA, iModelB, iLocComb, isubj, :, :, :) = margParams_SF_allBoot; % directly loaded
                margTunC_SF_allCond(iModelA, iModelB, iLocComb, isubj, :, :, :) = margTuningC_SF_allBoot;

                % NOM: IVs, predictions, nLL and parameters
                IV_allCond(iModelA, iModelB, iLocComb, isubj, :, :) = IV_allBoot;
                nTrials_allCond(iModelA, iModelB, iLocComb, isubj, :, :) = nTrials_allBoot;
                metric_data_allCond(iModelA, iModelB, iLocComb, :, isubj, :, :) = metric_data_allBoot;
                metric_pred_allCond(iModelA, iModelB, iLocComb, :, isubj, :, :) = metric_pred_allBoot;
                nLL_allCond(iModelA, iModelB, iLocComb, isubj, :) = nLL_allBoot;
                % IC_nLL_allCond is compiled above
                params_allCond(iModelA, iModelB, iLocComb, isubj, :, 1:nParamsB) = params_est_allBoot;
                R2_NOM_allCond(iModelA, iModelB, iLocComb, :, isubj, :) = R2_NOM_allBoot;

                % clear *allBoot

            end % isubj
        end % iLocComb
    end % iModelB
end % iModelA

fprintf('\nNOM outputs compiled\n')

% Save the organized data for all subjects
save(sprintf('%s/n%d_n%d', nameFolder_Data_NOM_Trialwise, nsubj, nBoot), '*_allCond')

%% Compile behavioral data (revise this part later, as the behav should also come from bootstrrapped data (metric_data_allCond), not from the raw data)
CS_allSubj = nan(nsubj, 8); % 8 is max number of combined locs,see namesLocComb
dprime_allSubj = CS_allSubj;
criterion_allSubj = CS_allSubj;
RT_allSubj = CS_allSubj;
pC_allSubj = CS_allSubj;
% pYES_allSubj = CS_allSubj;
pA_allSubj = CS_allSubj;

for isubj=1:nsubj
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    nameFile_behavMeas = sprintf('%s/%s%d/%s_behavMeas.mat', nameFolder_Data_OOD, subjName, nblocks, subjName);

    load(nameFile_behavMeas)

    for iLocComb = iLocComb_all
        switch iLocComb
            case 6, iLoc_allLoc = [2,4];
            case 7, iLoc_allLoc = [5,3];
            case 8, iLoc_allLoc = 2:5;
            otherwise
                iLoc_allLoc = iLocComb;
        end
        CS_allSubj(isubj, iLocComb) = 1./mean(cst_perSess_perLoc(:, iLoc_allLoc), 'all');
        dprime_allSubj(isubj, iLocComb) = mean(dprime_perSess_perLoc(:, iLoc_allLoc), 'all');
        criterion_allSubj(isubj, iLocComb) = mean(criterion_perSess_perLoc(:, iLoc_allLoc), 'all');
        RT_allSubj(isubj, iLocComb) = median(RT_perSess_perLoc(:, iLoc_allLoc), 'all');
        pC_allSubj(isubj, iLocComb) = mean(pC3_perSess_perLoc(:, iLoc_allLoc, 1), 'all');
        % pYES_allSubj(isubj, iLocComb) = mean(pYES_perSess_perLoc(:, iLoc_allLoc));
        pA_allSubj(isubj, iLocComb) = mean(pA3_perSess_perLoc(:, iLoc_allLoc, 1), 'all');
    end

    delete *_allT *_perSess *_perLoc dataMatrix
end

% Save the organized data for all subjects
save(sprintf('%s/n%d_n%d', nameFolder_Data_NOM_Trialwise, nsubj, nBoot), '*_allSubj', '-append')

fprintf('\nBehav data compiled\n')

%% Plot templates for IDVD and group averages
clc, fprintf('\n\n Plotting STARTING\n\n')
% Define folder for saving figures
nameFolder_Fig_NOM_Template = sprintf('%s/Template', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_NOM_Template)), mkdir(nameFolder_Fig_NOM_Template), end

for iDataset = 1:nDatasets
    switch iDataset
        case 1
            template_allSubj = squeeze(template_train_allCond(iModelA_plot, iModelB_plot, :, :, :, :)); % ... nLoc x nSubj x nORI x nSF
        case 2
            template_allSubj = squeeze(template_full_allCond(iModelA_plot, iModelB_plot, :, :, :, :)); % ... nLoc x nSubj x nORI x nSF
    end
    str_dataset = namesDataset{iDataset};
    %%%% Group ave %%%%%

    for iiLoc = 1:nLocComb8
        figure('Position', [0 0 1e3 1e3])
        e2D_ave = getCI(template_allSubj(iLocComb_all(iiLoc), :, :, :), 2, 2);
        e2D_ave = e2D_ave';
        fprintf('\n%s L%d: Min = %.3f, Max=%.3f\n', str_dataset, iLocComb_all(iiLoc), min(e2D_ave(:)), max(e2D_ave(:)))
        % subplot(1, nLocComb8, iiLoc),
        hold on
        if iLocComb_all(iiLoc)==1, cLim = [-.03, .26];
        else, cLim = [-.03, .13];
        end
        RCplot_2Dkernel(e2D_ave, cLim)
        title(sprintf('n=%d L%d %s (%s)', nsubj, iLocComb_all(iiLoc), namesLocComb{iLocComb_all(iiLoc)}, str_dataset))
        % set(findall(gcf, '-property', 'FontSize'), 'FontSize',12)
        % sgtitle(sprintf('n=%d [A%dB%d] (%s)', nsubj, iModelA_plot, iModelB_plot, str_title), 'FontSize',20)

        % title(namesLocComb{iLocComb_all(1)}, 'FontSize', sz_title)
        saveas(gcf, sprintf('%s/n%d_group_L%d_A%dB%d_%s.jpg', nameFolder_Fig_NOM_Template, nsubj, iLocComb_all(iiLoc), iModelA_plot, iModelB_plot, str_dataset))

    end % iiLoc
    % save
    % saveas(gcf, sprintf('%s/n%d_group_A%dB%d_%s.jpg', nameFolder_Fig_NOM_Template, nsubj, iModelA_plot, iModelB_plot, str_title))
    close all

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
    %     saveas(gcf, sprintf('%s/n%d_L%d_A%dB%d_%s.jpg', nameFolder_Fig_NOM_Template, nsubj, iLocComb_all(iiLoc), iModelA_plot, iModelB_plot, str_title))
    %     close all
    %
    % end % iiLoc
end % i
fprintf('\n\n Plotting DONE\n\n')

%% Plot separability
clc, fprintf('\n\n Plotting STARTING\n\n')
str_dataset = namesDataset{iDataset_plot};
% Define folder for saving figures
nameFolder_Fig_Sep = sprintf('%s/Separability', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_Sep)), mkdir(nameFolder_Fig_Sep), end

iLoc_Sep_all = [1,6,5,3];
wd_border = 2;
sz_title = 10;
sz_label = 25; %40
sz_ticks = 35;
wd_all = 3;
sz_marker = 20;

[sep_med, ~, ~, sep_neg, sep_pos] = getCI(sep_allCond(iModelA_plot, iModelB_plot, iLoc_Sep_all, :, :, iDataset_plot), 1, 5);

figure('Position', [0 0 1e3 600])
hold on, box on
yline(.5, '--', 'linewidth', wd_all);

str_sepValues = [];
for iiLoc = 1:length(iLoc_Sep_all)
    x = (1:nsubj) + .1 * iiLoc;
    for isubj = 1:nsubj
        errorbar(x(isubj), sep_med(iiLoc, isubj), sep_neg(iiLoc, isubj), sep_pos(iiLoc, isubj), 'CapSize', 0)
        plot(x(isubj), sep_med(iiLoc, isubj), markers_allSubj{isubj}, ...
            'MarkerEdgeColor', colors_comb(iLocComb_all(iLoc_Sep_all(iiLoc)), :), 'MarkerFaceColor', 'w', 'MarkerSize', sz_marker, 'linewidth', 3)
    end
    str_sepValues = [str_sepValues, sprintf('L%d: %.2f (%.2f); ', iLoc_Sep_all(iiLoc), mean(sep_med(iiLoc, :)), std(sep_med(iiLoc, :))/sqrt(nsubj))];

end

ylabel('Correlation', 'FontSize', sz_label)
yticks([0, .5:.1:1])
ylim([.5, 1])

xticks(1:nsubj)
xticklabels(1:nsubj)
xlim([0,nsubj+1])
xlabel('Observer #')

ax = gca; ax.XAxis.FontSize = sz_ticks; ax.YAxis.FontSize = sz_ticks; ax.LineWidth = wd_border;

%%% ANOVA
text_ANOVA = '...ANOVA...';
indLoc = repmat(1:length(iLoc_Sep_all), nsubj, 1)';
text_ANOVA = print_nANOVA({'Loc'}, sep_med(:), {indLoc(:)}, nsubj);
%=== CI of ANOVA ===
p_allBoot = nan(nBoot, 1);
for iBoot = 1:nBoot
    x = squeeze(sep_allCond(iModelA_plot, iModelB_plot, iLoc_Sep_all, :, iBoot, iDataset_plot));
    [~, tbl] = print_nANOVA({'Loc'}, x(:), {indLoc(:)}, nsubj);
    p_allBoot(iBoot) = tbl{2,7};
end
[p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 1);
text_ANOVA = [text_ANOVA, sprintf('p=%.3f [%.3f, %.3f]\n', p_med, p_lb, p_ub)];
%===============

%%% ttest
if length(iLoc_Sep_all)==2
    [h,p,ci,stats] = ttest(sep_med(:, 1), sep_med(:, 2));
    cohenD = fxn_getES(sep_med(:, 1), sep_med(:, 2));
end

set(findall(gcf, '-property', 'linewidth'), 'linewidth', 2)
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 30)

title(sprintf('n%d [A%dB%d] %s\n%s\n%s', nsubj, iModelA_plot, iModelB_plot, str_dataset, str_sepValues, text_ANOVA), 'FontSize', sz_title)
saveas(gcf, sprintf('%s/n%d_L%s_A%dB%d_%s.jpg', nameFolder_Fig_Sep, nsubj, strjoin(string(iLoc_Sep_all), ''), iModelA_plot, iModelB_plot, str_dataset))

%% Plot tuning functions for idvd
clc, fprintf('\n\n Plotting STARTING\n\n')
% Define folder for saving figures
nameFolder_Fig_NOM_Tuning = sprintf('%s/Tuning', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_NOM_Tuning)), mkdir(nameFolder_Fig_NOM_Tuning), end

wd_border = 4; % default 5
sz_ticks = 30;% default 35

for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
    iLocPair_all = iLocGroups_all{iGroup};

    for iFeature = 1:nFeatures

        xaxis = axis_tuning{iFeature};

        for iDataset=1:nDatasets
            switch iDataset
                case 1
                    markerStyle = 'o';
                    lineStyle = '-';
                    str_dataset = 'TrainingSet';
                case 2
                    markerStyle = 's';
                    lineStyle = '--';
                    str_dataset = 'FullSet';
            end

            figure('Position', [0 0 2e3 2e3])

            for isubj = 1:nsubj
                subplot(nRows_subj, nCols_subj, isubj), hold on

                subjName = subjList{isubj};

                for iLoc = iLocPair_all
                    % fprintf('\nL%d...', iLoc)
                    % Set color
                    color_comb = colors_comb(iLoc, :);
                    % Set y-ticks
                    if iFeature==1
                        marg_allCond = margORI_allCond;
                        margPred_allCond = margPred_ORI_allCond;

                        if iLoc==1, yticks_ = [-.03, 0, .05, .10, .15]; % fov vs. peri, higher ub
                        else, yticks_ = [-.02, linspace(0, .12, 4)];
                        end
                    else
                        marg_allCond = margSF_allCond;
                        margPred_allCond = margPred_SF_allCond;

                        if iLoc==1, yticks_ = [-.01, linspace(0, .08, 4)];
                        else, yticks_ = [-.01, linspace(0, .06, 4)]; %[-.02, 0, .02, .04, .06];
                        end
                    end

                    ymax = max(yticks_);
                    ymin = min(yticks_);

                    % Obtain data and prediction
                    [marg_med, ~, ~, marg_lb, marg_ub] = getCI(marg_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :), 1, 5);
                    [margPred_med, margPred_lb, margPred_ub] = getCI(margPred_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :), 1, 5);

                    % Data (dots + errorbars)
                    errorbar(xaxis, marg_med, marg_lb, marg_ub, markerStyle, 'Color', color_comb, 'CapSize',0)
                    plot(xaxis, marg_med, markerStyle, 'color', color_comb, 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'HandleVisibility','off')

                    % Prediction (lines + bands)
                    patch([xaxis, flip(xaxis)], [margPred_lb', flip(margPred_ub')], color_comb, 'FaceAlpha', .3, 'linestyle', 'none')
                    plot(xaxis, margPred_med, '-', 'color', color_comb)

                end % iLoc

                % Draw reference lines
                yline(0, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);
                xline(iFeature-1, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);

                xlabel('Orientation (deg)')
                ylabel('Amplitude')
                title(sprintf('%s', subjName))
                % legend('Location', 'best')

            end % isubj
            sgtitle(sprintf('n=%d L%d%d [A%dB%d] %s tuning (%s)', nsubj, iLocPair_all,iModelA_plot, iModelB_plot, namesFeature{iFeature}, str_dataset))
            set(findall(gcf, '-property', 'fontsize'), 'fontsize', 12)

            % Save the figure
            saveas(gcf, sprintf('%s/n%d_L%d%d_A%dB%d_%s_%s.jpg', nameFolder_Fig_NOM_Tuning, nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, str_dataset))
            close all
        end % i=1:2
    end % iFeature
end % iGroup

fprintf('\n\n Plotting DONE\n\n')

%% Plot tuning functions for group averages
clc, fprintf('\n\n Plotting STARTING\n\n')
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
            str_dataset = 'TrainingSet';
        case 2
            markerStyle = 's';
            lineStyle = '--';
            str_dataset = 'FullSet';
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

                % Data (dots + errorbars)
                errorbar(xaxis, marg_ave, marg_sem, markerStyle, 'Color', color_comb, 'CapSize',0)
                plot(xaxis, marg_ave, 'o', 'color', color_comb, 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'HandleVisibility','off')

                % Prediction (lines + bands)
                patch([xaxis, flip(xaxis)], [margPred_ave-margPred_sem, flip(margPred_ave+margPred_sem)], color_comb, 'FaceAlpha', .3, 'linestyle', 'none')
                plot(xaxis, margPred_ave, '-', 'color', color_comb)

            end % iLoc

            % Draw reference lines
            yline(0, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);
            xline(iFeature-1, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);

            xlabel(namesFeature_axis{iFeature})
            ylabel(namesFeature_axis_Tuning{iFeature})

            xlim(axisLim{iFeature})
            %     xlabel(xlabels_tuning{ifeature}, 'FontSize', sz_label)
            xticks(axisTicks_tuning{iFeature})
            xticklabels(axisTL_tuning{iFeature})
            if iFeature==2, xticklabels(round(axisTL_tuning{iFeature}, 2)), end

            yticks(yticks_)
            yticklabels(round(yticks_, 2))
            ylim(yticks_([1, end]))

            ax = gca;
            ax.XAxis.FontSize = sz_ticks;
            ax.YAxis.FontSize = sz_ticks;
            ax.LineWidth = wd_border/1.5;

            [R2_ave, ~, ~, R2_sem] = getCI(getCI(margR2_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset), 1, 5), 2, 2);

            title(sprintf('n=%d L%d%d [A%dB%d] %s tuning (%s)\nR2: L%d: %.2f (+-%.2f) L%d: %.2f (+-%.2f) ', ...
                nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, str_dataset, ...
                iLocPair_all(1), R2_ave(1), R2_sem(1), iLocPair_all(2), R2_ave(2), R2_sem(2)))
            % set(findall(gcf, '-property', 'fontsize'), 'fontsize', 15)
            set(findall(gcf, '-property', 'linewidth'), 'linewidth', 2)

            % Save the figure
            saveas(gcf, sprintf('%s/n%d_group_L%d%d_A%dB%d_%s_%s.jpg', nameFolder_Fig_NOM_Tuning, nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, str_dataset))
            close all
        end % iFeature
    end % iGroup


end % i=1:2
fprintf('\n\n Plotting DONE\n\n')

%% Tuning characteristics
nameFolder_Fig_tunC = sprintf('%s/TuningCs', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_tunC)), mkdir(nameFolder_Fig_tunC), end

flag_plotIDVD = 1;
flag_plotDiff = 1;
paramMode = 2;
sz_fig = [350 300]; % size of the figure canvas

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

        % Obtain y_ticks of ORI/SF domain (3 values)
        switch iFamily_perF(iFeature)
            case 1 % ORI tuniningC | Scaled Gaussian
                if flag_plotIDVD
                    if find(iLocPair_all==1)
                        y_ticks_all{1} = linspace(0, .3, 5); % ORI peak amp
                    else
                        y_ticks_all{1} = linspace(0, .2, 5); % ORI peak amp
                    end
                    y_ticks_all{2} = linspace(10, 50, 5); % ORI band
                    y_ticks_all{3} = linspace(-.05, .03, 5); % ORI baseline
                else
                    if find(iLocPair_all==1)
                        y_ticks_all{1} = linspace(0, .24, 5); % ORI peak amp
                    else
                        y_ticks_all{1} = linspace(0, .12, 5); % ORI peak amp
                    end

                    y_ticks_all{2} = linspace(10, 50, 5); % ORI band
                    y_ticks_all{3} = linspace(-.05, .03, 5); % ORI baseline
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
                    y_ticks_all{4} = linspace(-.08, 0, 5); % SF baseline
                else
                    y_ticks_all{1} = linspace(0, 2, 5); % peak SF

                    if find(iLocPair_all==1)
                        y_ticks_all{2} = linspace(.02, .12, 5); % SF peak amp
                    else
                        y_ticks_all{2} = linspace(.02, .06, 5); % SF peak amp
                    end
                    y_ticks_all{3} = linspace(.3, 1.5, 5); % SF bandwidth
                    y_ticks_all{4} = linspace(-.08, 0, 5); % SF baseline
                end

        end

        % Loop through each tuning characteristic
        for iTunC = 1:nTunC_full

            text_title = sprintf('%s %s', namesFeature{iFeature}, namesTunC{iTunC});

            % Compute median for each observer
            % margTunC_ORI_allCond: nModelA x nModelB x nLoc_pair x nsubj x nBoot x nDataset x nTunC
            % tunC_allSubj_allBoot: nLoc_pair x nsubj x nBoot
            switch iFeature
                case 1, Y_allSubj_allBoot = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plot, iTunC));
                case 2, Y_allSubj_allBoot = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plot, iTunC));
            end

            tunC_med_allSubj = getCI(Y_allSubj_allBoot, 1, 3)'; % rotate to match the format of basicFxn_drawBars()

            % Obtain ytick labels and ref
            y_ticklabels = round(y_ticks_all{iTunC}, 2);
            ref = nan;
            switch iFamily_perF(iFeature)
                case 1, if find(iTunC==3), ref = 0; end

                case 2 % log parabola
                    switch iTunC
                        case 1, y_ticklabels = round(2.^y_ticks_all{iTunC}, 1); ref = log2(2); tunC_med_allSubj = log2(tunC_med_allSubj); % pref SF
                        case 4, ref = 0; % baseline
                    end
            end

            % Obtain CI of t-test statistics (to print in the title)
            p_allBoot = nan(nBoot, 1);
            t_allBoot = p_allBoot;
            cohenD_allBoot = p_allBoot;
            for iBoot = 1:nBoot
                x=squeeze(Y_allSubj_allBoot(:, :, iBoot));
                [~, p, ~, stats] = ttest(x(1, :), x(2, :));
                cohenD = fxn_getES(x(1, :), x(2, :));
                p_allBoot(iBoot) = p;
                t_allBoot(iBoot) = stats.tstat;
                cohenD_allBoot(iBoot) = cohenD;
            end
            [p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 1);
            [t_med, t_lb, t_ub] = getCI(t_allBoot, 1, 1);
            [d_med, d_lb, d_ub] = getCI(cohenD_allBoot, 1, 1);

            % Define the title for the plot
            title_CI = sprintf('t=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f], d=%.3f [%.3f, %.3f]', t_med, t_lb, t_ub, p_med, p_lb, p_ub, d_med, d_lb, d_ub);

            % Plot the data
            text_title = sprintf('n=%d L%d%d [A%dB%d] [%s]\n%s\n%s', nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesDataset{iDataset_plot}, text_title, title_CI);
            %------------------------------%
            flag_sig = basicFxn_drawBars(tunC_med_allSubj, ref, colors, x_ticks, y_ticks_all{iTunC}, y_ticklabels, flag_plotIDVD, flag_plotDiff, text_title, 0, sz_fig);
            %------------------------------%
            ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily, 2}{iTunC}))
            % Save the figure
            saveas(gcf, sprintf('%s/n%d_L%d%d_A%dB%d_%s%d%s.jpg', nameFolder_Fig_tunC, nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, iTunC, flag_sig))

        end % end of iTunC
    end % end of iFeature

    close all
end % end of iGroup

%% Behavioral metrics in grouped locations
nameFolder_Fig_behav = sprintf('%s/behav', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_behav)), mkdir(nameFolder_Fig_behav), end

sz_wd_perBar = 100;
nBars = 2;
sz_fig = [nBars*sz_wd_perBar, 250];
flag_plotIDVD = 1;
namesMetrics_behav = {'CS', 'pA', 'dprime', 'criterion', 'pC', 'RT'}; nMetrics_behav = length(namesMetrics_behav);
[d70,~] = SX_sim06_SDT(.7, .3);
iDataset_vec = 1; % 1=full, 2=test; see OOD_xx_compIV.m search for "data_metrics_allBoot"
for iGroup=1:nGroups
    iLocPair_all = iLocGroups_all{iGroup};
    for iMetric=1:nMetrics_behav
        switch iMetric
            case 1, iMetric_vec = 10; x_ticks = linspace(1.6, 3.6, 5); flag_plotIDVD = 1; ref=nan;
            case 2, iMetric_vec = 6; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 1; ref=nan;
            case 3, iMetric_vec = 1; x_ticks = linspace(0, 1.6, 5); flag_plotIDVD = 0; ref=d70;
            case 4, iMetric_vec = 2; x_ticks = linspace(-1,1, 5); flag_plotIDVD = 0; ref=0;
            case 5, iMetric_vec = 3; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 0; ref=nan;
            case 6, iMetric_vec = 11; x_ticks = linspace(0, .2, 5); flag_plotIDVD = 0; ref=nan;
        end
        x_ticks = round(x_ticks, 2);
        % metrics_allCond = nan(nModelsA, nLocComb8, nsubj, nBoot, 2, 11); % 2: 1=Full; 2=Test; 11=number of metrics saved in OOD_xx_compIV.m; see "data_metrics_allBoot"
        data_med_allSubj = getCI(metrics_allCond(iModelA_plot, iLocPair_all, :, :, iDataset_vec, iMetric_vec), 1, 4)';
        % data_allSubj = data_allSubj(:, iLocPair_all);

        basicFxn_drawBars(data_med_allSubj, ref, colors_comb(iLocPair_all, :), namesLocComb(iLocPair_all), x_ticks, x_ticks, flag_plotIDVD, 1, namesMetrics_behav{iMetric}, 0, sz_fig);
        saveas(gcf, sprintf('%s/n%d_L%d%d_%s.jpg', nameFolder_Fig_behav, nsubj, iLocPair_all, namesMetrics_behav{iMetric}))
    end % iMetric
    close all
end % iGroup

%% [1.1] Corr between CS and tunC
nameVarX = 'CS';
nameVarY = 'tunC';
sz_label = 55;
sz_labelOffset = .05;
sz_axOffset = 0.05; % extra breathing room

nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

iLocCorr_all = [1,6,5,3];
flag_zeroMean = 0;
type_corr='pearson';
type_tail = 'both';

% Obtain CS (x-axis)
X_allSubj = CS_allSubj(:, iLocCorr_all);
x_ticks = linspace(1.5, 3.5, 5); % EE
X_med_allSubj = X_allSubj;

for iFeature = 1:nFeatures
    iFamily = iFamily_perF(iFeature);
    namesTunCs = namesTunC_unit_perF{iFamily, 2};
    nTunCs_full = length(namesTunCs);
    for iTunC = 1:nTunCs_full
        switch iFeature
            case 1, Y_allSubj_allBoot = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plot, iTunC));
            case 2, Y_allSubj_allBoot = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plot, iTunC));
        end

        Y_med_allSubj = getCI(Y_allSubj_allBoot, 1, 3)'; % rotate to match the format needed by basicFxn_drawCorr

        switch iFeature
            case 1, y_ticks_allTunC_lb = [0, 10, -.1]; y_ticks_allTunC_ub = [.36, 50, .1];
            case 2, y_ticks_allTunC_lb = [0, 0, 0, -.1]; y_ticks_allTunC_ub = [4, .24, 1.2, .1];
        end

        y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);

        x_ticklabels = x_ticks;
        y_ticklabels = y_ticks;

        nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily_perF(iFeature), 2}{iTunC});
        nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

        text_title = sprintf('%s vs. %s', nameVarX, nameVarY_figTitle);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        text_CI = ''; r_allBoot = nan(nBoot, 1); p_allBoot = r_allBoot;
        % for iBoot=1:nBoot,[r, p] = corr(X_allSubj, Y_allSubj_allBoot(:, iBoot), 'Type', 'Kendall', 'Tail','right'); p_allBoot(iBoot) = p; r_allBoot(iBoot) = r; end
        % [p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 1); [r_med, r_lb, r_ub] = getCI(r_allBoot, 1, 1);
        % text_CI = sprintf('r=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f]', r_med, r_lb, r_ub, p_med, p_lb, p_ub);
        text_title = sprintf('%s\n%s', text_title, text_CI);

        flag_sig = basicFxn_drawCorr(X_med_allSubj, Y_med_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, type_corr, type_tail, text_title, markers_allSubj);

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

        saveas(gcf, sprintf('%s/n%d_L%s_%s.jpg', nameFolder_Fig_NOM_corr, nsubj, strjoin(string(iLocCorr_all), ''), nameVarY_fileTitle))

    end % iTunC
    close all
end % iFeature


%% [1.2] Corr between extents of EE/HVA/VMA (CS and tunC)
nameVarX = 'CS';
nameVarY = 'tunC';

nameFolder_Fig_NOM_corrAsym = sprintf('%s/corrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corrAsym)); mkdir(nameFolder_Fig_NOM_corrAsym), end

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
        case 5, x_ticks = linspace(0, 16, 5); % VMA (extent is smaller)
    end
    asymX_med_allSubj = asymX_allSubj;

    for iFeature = 1:nFeatures
        iFamily = iFamily_perF(iFeature);
        namesTunCs = namesTunC_unit_perF{iFamily, 2};
        nTunCs_full = length(namesTunCs);
        for iTunC = 1:nTunCs_full
            switch iFeature
                case 1, Y_allSubj_allBoot = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plot, iTunC));
                case 2, Y_allSubj_allBoot = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plot, iTunC));
            end

            asymY_allSubj = squeeze((Y_allSubj_allBoot(1, :, :)-Y_allSubj_allBoot(2, :, :))./(Y_allSubj_allBoot(1, :, :)+Y_allSubj_allBoot(2, :, :)));
            asymY_med_allSubj = getCI(asymY_allSubj, 1, 2);

            switch iLocPair_all(1)
                case 1 % EE
                    switch iFamily
                        case 1, y_ticks_allTunC_lb = -[0, 20, 100]; y_ticks_allTunC_ub = [72, 20, 100];
                            % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                        case 2, y_ticks_allTunC_lb = -[50, 0, 50, 145]; y_ticks_allTunC_ub = [50, 80, 50, 265];
                    end
                case 6 % HVA
                    switch iFamily
                        case 1, y_ticks_allTunC_lb = -[40, 40, 100]; y_ticks_allTunC_ub = [40, 40, 100];
                            % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                        case 2, y_ticks_allTunC_lb = -[50, 20, 60, 145]; y_ticks_allTunC_ub = [50, 60, 60, 265];
                    end
                case 5 % VMA
                    switch iFamily
                        case 1, y_ticks_allTunC_lb = -[30, 50, 50]; y_ticks_allTunC_ub = [90, 30, 50];
                            % case 8, y_ticks_allTunC_lb = -[100, 40, 300, 320, 50, 320]; y_ticks_allTunC_ub = [100, 44, 300, 320, 50, 320];
                        case 2, y_ticks_allTunC_lb = -[50, 40, 60, 250]; y_ticks_allTunC_ub = [50, 100, 40, 250];
                    end
            end
            y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);

            x_ticklabels = nan;
            y_ticklabels = nan;

            nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily_perF(iFeature), 2}{iTunC});
            nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

            text_title = sprintf('%s (%s) vs. %s (%s)', nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            text_CI = ''; r_allBoot = nan(nBoot, 1); p_allBoot = r_allBoot;
            for iBoot=1:nBoot,[r, p] = corr(asymX_allSubj, asymY_allSubj(:, iBoot), 'Type', 'Kendall', 'Tail','right'); p_allBoot(iBoot) = p; r_allBoot(iBoot) = r; end
            [p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 1); [r_med, r_lb, r_ub] = getCI(r_allBoot, 1, 1);
            text_CI = sprintf('r=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f]', r_med, r_lb, r_ub, p_med, p_lb, p_ub);
            text_title = sprintf('%s\n%s', text_title, text_CI);

            flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj*100, asymY_med_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, text_title, markers_allSubj);
            % flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj*100, asymY_med_allSubj*100, [], [], [], [], text_title, markers_allSubj);
            xlabel(sprintf('%s of contrast sensitivity (%%)', nameAsymX), 'fontsize', sz_label)
            ylabel(sprintf('%s of %s %s (%%)', nameAsymX, namesFeature{iFeature}, namesTunC_noUnit{iFamily_perF(iFeature), 2}{iTunC}), 'fontsize', sz_label)

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

            saveas(gcf, sprintf('%s/n%d_L%d%d_%s.jpg', nameFolder_Fig_NOM_corrAsym, nsubj, iLocPair_all, nameVarY_fileTitle))

            close all

        end % iTunC

    end % iFeature

end % iGroup

%% [2.1] Corr between CS and NOM params
nameVarX = 'CS';
nameVarY = 'NOMparams';
sz_label = 55;
sz_labelOffset = .05;
sz_axOffset = 0.05; % extra breathing room

nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

iLocCorr_all = [1,6,5,3];
flag_zeroMean = 0;
type_corr='pearson';
type_tail = 'both';

% Obtain CS (x-axis)
X_allSubj = CS_allSubj(:, iLocCorr_all);
x_ticks = linspace(1.5, 3.5, 5); % EE
X_med_allSubj = X_allSubj;

nNOMparams = length(namesModelBparams{iModelB_plot});

for iNOMparam = 1:nNOMparams
    Y_allSubj_allBoot = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iNOMparam));

    Y_med_allSubj = getCI(Y_allSubj_allBoot, 1, 3)'; % rotate to match the format needed by basicFxn_drawCorr

    y_ticks_lb = [0, 5, 0]; y_ticks_ub = [.8, 45, .8];

    y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

    x_ticklabels = x_ticks;
    y_ticklabels = y_ticks;

    nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
    nameVarY_fileTitle = nameVarY_figTitle;

    text_title = sprintf('%s vs. %s', nameVarX, nameVarY_figTitle);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    text_CI = ''; r_allBoot = nan(nBoot, 1); p_allBoot = r_allBoot;
    % for iBoot=1:nBoot,[r, p] = corr(X_allSubj, Y_allSubj_allBoot(:, iBoot), 'Type', 'Kendall', 'Tail','right'); p_allBoot(iBoot) = p; r_allBoot(iBoot) = r; end
    % [p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 1); [r_med, r_lb, r_ub] = getCI(r_allBoot, 1, 1);
    % text_CI = sprintf('r=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f]', r_med, r_lb, r_ub, p_med, p_lb, p_ub);
    text_title = sprintf('%s\n%s', text_title, text_CI);

    flag_sig = basicFxn_drawCorr(X_med_allSubj, Y_med_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, type_corr, type_tail, text_title, markers_allSubj);
    % flag_sig = basicFxn_drawCorr(X_med_allSubj, Y_med_allSubj, colors_comb(iLocCorr_all, :), [], [] ,[] ,[], flag_zeroMean, type_corr, type_tail, text_title, markers_allSubj);

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
    saveas(gcf, sprintf('%s/n%d_L%s_NOM%d.jpg', nameFolder_Fig_NOM_corr, nsubj, strjoin(string(iLocCorr_all), ''), iNOMparam))

    close all
end % iNOMparam

%% [2.2] Corr between extents of EE/HVA/VMA (CS and NOM params)
nameVarX = 'CS';
nameVarY = 'NOMparams';

nameFolder_Fig_NOM_corrAsym = sprintf('%s/corrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corrAsym)); mkdir(nameFolder_Fig_NOM_corrAsym), end

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
        case 5, x_ticks = linspace(0, 16, 5); % VMA (extent is smaller)
    end
    asymX_med_allSubj = asymX_allSubj;

    for iNOMparam = 1:nNOMparams
        Y_allSubj_allBoot = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iNOMparam));

        asymY_allSubj = squeeze((Y_allSubj_allBoot(1, :, :)-Y_allSubj_allBoot(2, :, :))./(Y_allSubj_allBoot(1, :, :)+Y_allSubj_allBoot(2, :, :)));
        asymY_med_allSubj = getCI(asymY_allSubj, 1, 2);

        switch iLocPair_all(1)
            case 1 % EE
                y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
            case 6 % HVA
                y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
            case 5 % VMA
                y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        end
        y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

        x_ticklabels = nan;
        y_ticklabels = nan;

        nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
        nameVarY_fileTitle = nameVarY_figTitle;

        text_title = sprintf('%s (%s) vs. %s (%s)', nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        text_CI = ''; r_allBoot = nan(nBoot, 1); p_allBoot = r_allBoot;
        for iBoot=1:nBoot,[r, p] = corr(asymX_allSubj, asymY_allSubj(:, iBoot), 'Type', 'Kendall', 'Tail','right'); p_allBoot(iBoot) = p; r_allBoot(iBoot) = r; end
        [p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 1); [r_med, r_lb, r_ub] = getCI(r_allBoot, 1, 1);
        text_CI = sprintf('r=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f]', r_med, r_lb, r_ub, p_med, p_lb, p_ub);
        text_title = sprintf('%s\n%s', text_title, text_CI);

        flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj*100, asymY_med_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, text_title, markers_allSubj);
        % flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj*100, asymY_med_allSubj*100, [], [], [], [], text_title, markers_allSubj);
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

        saveas(gcf, sprintf('%s/n%d_L%d%d_NOM%d.jpg', nameFolder_Fig_NOM_corrAsym, nsubj, iLocPair_all, iNOMparam))

        close all

    end % iNOMparam

end % iGroup

%% [3.1] Corr between CS and pA
nameVarX = 'CS';
nameVarY = 'pA'; iMetric_pA = 2;
sz_label = 55;
sz_labelOffset = .05;
sz_axOffset = 0.05; % extra breathing room

nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

iLocCorr_all = [1,6,5,3];
flag_zeroMean = 0;
type_corr='pearson';
type_tail = 'both';

% Obtain CS (x-axis)
X_allSubj = CS_allSubj(:, iLocCorr_all);
x_ticks = linspace(1.5, 3.5, 5); % EE
X_med_allSubj = X_allSubj;

Y_allSubj_allBoot = getCI(metric_data_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, iMetric_pA, :, :, :), 2, 7);

Y_med_allSubj = getCI(Y_allSubj_allBoot, 1, 3)'; % rotate to match the format needed by basicFxn_drawCorr


y_ticks = linspace(.5, 1, 5);

x_ticklabels = x_ticks;
y_ticklabels = y_ticks;

nameVarY_figTitle = nameVarY;
nameVarY_fileTitle = nameVarY_figTitle;

text_title = sprintf('%s vs. %s', nameVarX, nameVarY_figTitle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text_CI = ''; r_allBoot = nan(nBoot, 1); p_allBoot = r_allBoot;
% for iBoot=1:nBoot,[r, p] = corr(X_allSubj, Y_allSubj_allBoot(:, iBoot), 'Type', 'Kendall', 'Tail','right'); p_allBoot(iBoot) = p; r_allBoot(iBoot) = r; end
% [p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 1); [r_med, r_lb, r_ub] = getCI(r_allBoot, 1, 1);
% text_CI = sprintf('r=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f]', r_med, r_lb, r_ub, p_med, p_lb, p_ub);
text_title = sprintf('%s\n%s', text_title, text_CI);

flag_sig = basicFxn_drawCorr(X_med_allSubj, Y_med_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, type_corr, type_tail, text_title, markers_allSubj);
% flag_sig = basicFxn_drawCorr(X_med_allSubj, Y_med_allSubj, colors_comb(iLocCorr_all, :), [], [] ,[] ,[], flag_zeroMean, type_corr, type_tail, text_title, markers_allSubj);

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
saveas(gcf, sprintf('%s/n%d_L%s.jpg', nameFolder_Fig_NOM_corr, nsubj, strjoin(string(iLocCorr_all), '')))

close all
% end % iNOMparam

%% [3.2] Corr between extents of EE/HVA/VMA (CS and pA)
nameVarX = 'CS';
nameVarY = 'pA';

nameFolder_Fig_NOM_corrAsym = sprintf('%s/corrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corrAsym)); mkdir(nameFolder_Fig_NOM_corrAsym), end

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
        case 5, x_ticks = linspace(0, 16, 5); % VMA (extent is smaller)
    end
    asymX_med_allSubj = asymX_allSubj;

    % for iNOMparam = 1:nNOMparams
    Y_allSubj_allBoot = getCI(metric_data_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, iMetric_pA, :, :, :), 2, 7);

    asymY_allSubj = squeeze((Y_allSubj_allBoot(1, :, :)-Y_allSubj_allBoot(2, :, :))./(Y_allSubj_allBoot(1, :, :)+Y_allSubj_allBoot(2, :, :)));
    asymY_med_allSubj = getCI(asymY_allSubj, 1, 2);

    switch iLocPair_all(1)
        case 1 % EE
            y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        case 6 % HVA
            y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        case 5 % VMA
            y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
    end
    y_ticks = linspace(0,100, 5);

    x_ticklabels = nan;
    y_ticklabels = nan;

    nameVarY_figTitle = nameVarY;
    nameVarY_fileTitle = nameVarY_figTitle;

    text_title = sprintf('%s (%s) vs. %s (%s)', nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    text_CI = ''; r_allBoot = nan(nBoot, 1); p_allBoot = r_allBoot;
    for iBoot=1:nBoot,[r, p] = corr(asymX_allSubj, asymY_allSubj(:, iBoot), 'Type', 'Kendall', 'Tail','right'); p_allBoot(iBoot) = p; r_allBoot(iBoot) = r; end
    [p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 1); [r_med, r_lb, r_ub] = getCI(r_allBoot, 1, 1);
    text_CI = sprintf('r=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f]', r_med, r_lb, r_ub, p_med, p_lb, p_ub);
    text_title = sprintf('%s\n%s', text_title, text_CI);

    % flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj*100, asymY_med_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, text_title, markers_allSubj);
    flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj*100, asymY_med_allSubj*100, [], [], [], [], text_title, markers_allSubj);
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

    saveas(gcf, sprintf('%s/n%d_L%d%d.jpg', nameFolder_Fig_NOM_corrAsym, nsubj, iLocPair_all))

    close all

    % end % iNOMparam

end % iGroup

%% [4.1] Corr between CS and NOM params
nameVarX = 'pA';
nameVarY = 'NOMparams';
sz_label = 55;
sz_labelOffset = .05;
sz_axOffset = 0.05; % extra breathing room

nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

iLocCorr_all = [1,6,5,3];
flag_zeroMean = 0;
type_corr='pearson';
type_tail = 'both';

% Obtain CS (x-axis)
X_allSubj = pA_allSubj(:, iLocCorr_all);
x_ticks = linspace(.5, .9, 5); % EE
X_med_allSubj = X_allSubj;

nNOMparams = length(namesModelBparams{iModelB_plot});

for iNOMparam = 1:nNOMparams
    Y_allSubj_allBoot = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iNOMparam));

    Y_med_allSubj = getCI(Y_allSubj_allBoot, 1, 3)'; % rotate to match the format needed by basicFxn_drawCorr

    y_ticks_lb = [0, 5, 0]; y_ticks_ub = [.8, 45, .8];

    y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

    x_ticklabels = x_ticks;
    y_ticklabels = y_ticks;

    nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
    nameVarY_fileTitle = nameVarY_figTitle;

    text_title = sprintf('%s vs. %s', nameVarX, nameVarY_figTitle);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    text_CI = ''; r_allBoot = nan(nBoot, 1); p_allBoot = r_allBoot;
    % for iBoot=1:nBoot,[r, p] = corr(X_allSubj, Y_allSubj_allBoot(:, iBoot), 'Type', 'Kendall', 'Tail','right'); p_allBoot(iBoot) = p; r_allBoot(iBoot) = r; end
    % [p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 1); [r_med, r_lb, r_ub] = getCI(r_allBoot, 1, 1);
    % text_CI = sprintf('r=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f]', r_med, r_lb, r_ub, p_med, p_lb, p_ub);
    text_title = sprintf('%s\n%s', text_title, text_CI);

    flag_sig = basicFxn_drawCorr(X_med_allSubj, Y_med_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, type_corr, type_tail, text_title, markers_allSubj);
    % flag_sig = basicFxn_drawCorr(X_med_allSubj, Y_med_allSubj, colors_comb(iLocCorr_all, :), [], [] ,[] ,[], flag_zeroMean, type_corr, type_tail, text_title, markers_allSubj);

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
    saveas(gcf, sprintf('%s/n%d_L%s_NOM%d.jpg', nameFolder_Fig_NOM_corr, nsubj, strjoin(string(iLocCorr_all), ''), iNOMparam))

    close all
end % iNOMparam

%% [4.2] Corr between extents of EE/HVA/VMA (CS and NOM params)
nameVarX = 'pA';
nameVarY = 'NOMparams';

nameFolder_Fig_NOM_corrAsym = sprintf('%s/corrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corrAsym)); mkdir(nameFolder_Fig_NOM_corrAsym), end

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

    asymX_allSubj = (pA_allSubj(:, iLocPair_all(1))-pA_allSubj(:, iLocPair_all(2)))./(pA_allSubj(:, iLocPair_all(1))+pA_allSubj(:, iLocPair_all(2)));
    switch iLocPair_all(1)
        case 1, x_ticks = linspace(0, 20, 5); % EE
        case 6, x_ticks = linspace(-5, 15, 5); % HVA
        case 5, x_ticks = linspace(0, 16, 5); % VMA (extent is smaller)
    end
    asymX_med_allSubj = asymX_allSubj;

    for iNOMparam = 1:nNOMparams
        Y_allSubj_allBoot = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iNOMparam));

        asymY_allSubj = squeeze((Y_allSubj_allBoot(1, :, :)-Y_allSubj_allBoot(2, :, :))./(Y_allSubj_allBoot(1, :, :)+Y_allSubj_allBoot(2, :, :)));
        asymY_med_allSubj = getCI(asymY_allSubj, 1, 2);

        switch iLocPair_all(1)
            case 1 % EE
                y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
            case 6 % HVA
                y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
            case 5 % VMA
                y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        end
        y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

        x_ticklabels = nan;
        y_ticklabels = nan;

        nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
        nameVarY_fileTitle = nameVarY_figTitle;

        text_title = sprintf('%s (%s) vs. %s (%s)', nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        text_CI = ''; r_allBoot = nan(nBoot, 1); p_allBoot = r_allBoot;
        for iBoot=1:nBoot,[r, p] = corr(asymX_allSubj, asymY_allSubj(:, iBoot), 'Type', 'Kendall', 'Tail','right'); p_allBoot(iBoot) = p; r_allBoot(iBoot) = r; end
        [p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 1); [r_med, r_lb, r_ub] = getCI(r_allBoot, 1, 1);
        text_CI = sprintf('r=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f]', r_med, r_lb, r_ub, p_med, p_lb, p_ub);
        text_title = sprintf('%s\n%s', text_title, text_CI);

        % flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj*100, asymY_med_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, text_title, markers_allSubj);
        flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj*100, asymY_med_allSubj*100, [], [], [], [], text_title, markers_allSubj);
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

        saveas(gcf, sprintf('%s/n%d_L%d%d_NOM%d.jpg', nameFolder_Fig_NOM_corrAsym, nsubj, iLocPair_all, iNOMparam))

        close all

    end % iNOMparam

end % iGroup

%% Plot metrics vs. IV for group averages
clc, fprintf('\n\n Plotting STARTING\n\n')
lineStyle_all = {'-', '-', '--', '-', ':', '-', '--', ':', '-.'};

% Define folder for saving figures
nameFolder_Fig_NOM_Trialwise = sprintf('%s/NOM_Trialwise_%d%d', nameFolder_Figures, nORI, nSF);
if isempty(dir(nameFolder_Fig_NOM_Trialwise)), mkdir(nameFolder_Fig_NOM_Trialwise), end
iModelB_all_plot = [1,2,3];
namesModelB_plot = {'Full model', 'No $\rho$', 'No $N$', 'xx', 'No $\sigma_{constant}$', 'xx'};
% namesModelB_plot = {'Full model', 'No $\rho$', 'No $N$', 'No $\sigma_{constant}$', 'xx'};
nModelB_plot = length(iModelB_all_plot);

x_base   = 0.60;   % move a bit further left to make space for name column
x_width  = 0.35;   % total width of the table block
y_base   = 0.10;
y_height = 0.25;
sz_font = 7;

for iModelA = 1%iModelA_all % just plot pred of the model using RC-derived template (happens to be yhe best model)

    % Define folder for saving figures for each model A and model B
    nameFolder_Fig_NOM_metrics = sprintf('%s/A%d', nameFolder_Fig_NOM_Trialwise, iModelA);
    if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end

    % for iLocComb = iLocComb_all
    for iGroup = 1:nGroups
        iLocPair_all = iLocGroups_all{iGroup};

        fprintf('\nL%d%d...', iLocPair_all)
        scalingF = 80; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end

        figure('Position', [0 0 400 1e3])

        for iMetric = 1:nMetrics

            subplot(nMetrics, 1, iMetric), hold on

            R2_tab = nan(nModelB_plot, 2);   % [row = model, col = loc within pair]

            for iModelB = iModelB_all_plot
                
                for iiLoc = 1:2
                    % Compute medians and CIs
                    [IV_ave, ~, ~, IV_SEM] = getCI(getCI(IV_allCond(iModelA, iModelB, iLocPair_all(iiLoc), :, :, :), 1, 5), 2, 1);
                    [nTrials_ave, ~, ~, nTrials_SEM] = getCI(getCI(nTrials_allCond(iModelA, iModelB, iLocPair_all(iiLoc), :, :, :), 1, 5), 2, 1);
                    [data_ave, ~, ~, data_SEM] = getCI(getCI(metric_data_allCond(iModelA, iModelB, iLocPair_all(iiLoc), iMetric, :, :, :), 1, 6), 2, 1);
                    [pred_ave, ~, ~, pred_SEM] = getCI(getCI(metric_pred_allCond(iModelA, iModelB, iLocPair_all(iiLoc), iMetric, :, :, :), 1, 6), 2, 1);
                    [R2_NOM_ave, ~, ~, R2_NOM_SEM] = getCI(getCI(R2_NOM_allCond(iModelA, iModelB, iLocPair_all(iiLoc), iMetric, :, :), 1, 6), 2, 1);

                    % ---- Store R2 mean into table ----
                    rowIdx = find(iModelB_all_plot == iModelB);  % which row for this model
                    colIdx = iiLoc;                              % col 1 or 2 for the pair
                    R2_tab(rowIdx, colIdx) = R2_NOM_ave;        % store mean R2 (not SEM)

                    % Plot prediction
                    if iModelB==1, color_pred = colors_comb(iLocPair_all(iiLoc), :); lw = 3;
                    else, color_pred='k'; lw = 1;
                    end
                    patch([IV_ave; flip(IV_ave)], [pred_ave-pred_SEM; flip(pred_ave+pred_SEM)], ones(1, 3) / 2, 'FaceAlpha', .3, 'linestyle', 'none', 'handlevisibility', 'off')
                    plot(IV_ave, pred_ave, 'lineStyle', lineStyle_all{iModelB}, 'color', colors_comb(iLocPair_all(iiLoc), :), 'linewidth', lw)
                    
                    % Plot measurement
                    errorbar(IV_ave, data_ave, IV_SEM, '.', 'horizontal', 'CapSize', 0, 'color', colors_comb(iLocPair_all(iiLoc), :), 'handlevisibility', 'off')
                    errorbar(IV_ave, data_ave, data_SEM, '.', 'vertical', 'CapSize', 0, 'color', colors_comb(iLocPair_all(iiLoc), :), 'handlevisibility', 'off')

                    % Plot averaged data of each bin
                    for iBin = 1:nBins
                        facecolor = 'w';
                        plot(IV_ave(iBin), data_ave(iBin), 'o', 'markeredgecolor', colors_comb(iLocPair_all(iiLoc), :), ...
                            'markerfacecolor', facecolor, 'MarkerSize', nTrials_ave(iBin) / scalingF + 5, 'LineWidth', 1, 'LineStyle', 'none', 'handlevisibility', 'off')
                    end
                end % iiLoc

                % str_R2_table = str_R2_table + row_str + "\newline";

                yline(.5, '--', 'linewidth', 2);

                % xlim([,150])
                ylim([0, 1])
                yticks(0:.2:1)
                ylabel(namesMetrics{iMetric})
                xlabel('Internal variable')
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

            % ---- Column 1: model names (B1, B2, ...) ----
            for r = 1:nRows_subj
                % modelID = iModelB_all_plot(r);
                name_str = sprintf('%s', namesModelB_plot{iModelB_all_plot(r)});

                text(x_cells(1), y_cells(r), name_str, ...
                    'Units','normalized', ...
                    'HorizontalAlignment','right', ...
                    'VerticalAlignment','middle', ...
                    'FontSize', sz_font*1.5, 'FontWeight', 'bold', 'interpreter', 'latex');
            end

            % ---- Columns 2–3: R2 values ----
            for r = 1:nRows_subj
                for c = 1:2   % two locations
                    R2_val = R2_tab(r, c);
                    if isnan(R2_val), continue; end

                    str_cell = sprintf('%.0f%%', R2_val * 100);

                    text(x_cells(c+1), y_cells(r), str_cell, ...
                        'Units','normalized', ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', ...
                        'FontSize', sz_font);
                end
            end

            title(namesMetrics{iMetric})

        end % iMetric

        sgtitle(sprintf('n=%d [A%d] [L%d%d] [nBoot=%d] %s', nsubj, iModelA, iLocPair_all, nBoot, namesMetrics{iMetric}))
        set(findall(gcf, '-property', 'fontsize'), 'fontsize', 15)
        % set(findall(gcf, '-property', 'linewidth'), 'linewidth', 2)

        % Save the figure
        saveas(gcf, sprintf('%s/n%d_L%d%d_A%d_group.jpg', nameFolder_Fig_NOM_metrics, nsubj, iLocPair_all, iModelA_plot))
        close all
        fprintf('DONE\n')
    end % iLocComb

end % iModelA

fprintf('\n\n Plotting DONE\n\n')


%% Plot metrics vs. IV for each idvd
clc, fprintf('\n\n Plotting STARTING\n\n')
lineStyle_all = {'-', '-', '--', ':', '-.', '-', '--', ':', '-.'};

% iModelA_plot = 1;%iModelA_all % just plot pred of the model using RC-derived template (happens to be yhe best model)
% iModelA = iModelA_plot;
iModelB_all_plot = [1:3];
% Define folder for saving figures for each model A and model B
nameFolder_Fig_NOM_metrics = sprintf('%s/A%d', nameFolder_Fig_NOM_Trialwise, iModelA_plot);
if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end

for iLocComb = iLocComb_all
    fprintf('\nL%d...', iLocComb)
    scalingF = 80; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end

    %%%%%%%%%%%%
    for iMetric = 1:nMetrics
        figure('Position', [0 0 2e3 2e3])

        for isubj = 1:nsubj
            subplot(nRows_subj, nCols_subj, isubj), hold on

            subjName = subjList{isubj};

            % Create a string to store parameter estimates for display
            str_est = [];
            params_allCond(isnan(params_allCond))=0;
            for iModelB = iModelB_all_plot
                vals = strjoin(string(round(getCI(params_allCond(iModelA_plot, iModelB, iLocComb, isubj, :, :), 1, 5)',2)), ", ");
                str_est = [str_est, sprintf('B%d: %s\n', iModelB, vals)];
            end

            for iModelB = iModelB_all_plot

                % % Define file names for loading data
                % nameFile_compIV = sprintf('%s/n%d_A%dB%d', nameFolder_NOM_save, nBoot, iModelA, iModelB); % MUST be the same as OOD_NOM_Trialwise_compIV.m, Line 77

                % Compute medians and CIs
                [IV_allBins, ~, ~, IV_allBins_neg, IV_allBins_pos] = getCI(IV_allCond(iModelA_plot, iModelB, iLocComb, isubj, :, :), 1, 5);
                [nTrials_allBins, nTrials_allBins_lb, nTrials_allBins_ub] = getCI(nTrials_allCond(iModelA_plot, iModelB, iLocComb, isubj, :, :), 1, 5);
                [data_allBins, ~, ~, data_allBins_neg, data_allBins_pos] = getCI(metric_data_allCond(iModelA_plot, iModelB, iLocComb, iMetric, isubj, :, :), 1, 6);
                [pred_allBins, pred_allBins_lb, pred_allBins_ub] = getCI(metric_pred_allCond(iModelA_plot, iModelB, iLocComb, iMetric, isubj, :, :), 1, 6);
                [params_allBins, params_allBins_lb, params_allBins_ub] = getCI(params_allCond(iModelA_plot, iModelB, iLocComb, isubj, :, :), 1, 5);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Plot prediction
                if iModelB==1, color_pred = colors_comb(iLocComb, :); lw = 3;
                else, color_pred='k'; lw = 1;
                end
                patch([IV_allBins; flip(IV_allBins)], [pred_allBins_lb; flip(pred_allBins_ub)], ones(1, 3) / 2, 'FaceAlpha', .3, 'linestyle', 'none', 'handlevisibility', 'off')
                plot(IV_allBins, pred_allBins, 'lineStyle', lineStyle_all{iModelB}, 'color', color_pred, 'linewidth', lw)
                % measurement
                errorbar(IV_allBins, data_allBins, IV_allBins_neg, IV_allBins_pos, '.', 'horizontal', 'CapSize', 0, 'color', colors_comb(iLocComb, :), 'handlevisibility', 'off')
                errorbar(IV_allBins, data_allBins, data_allBins_neg, data_allBins_pos, '.', 'vertical', 'CapSize', 0, 'color', colors_comb(iLocComb, :), 'handlevisibility', 'off')

                yline(.5, 'k--');

                % Plot averaged data of each bin
                for iBin = 1:nBins
                    % if isubj>nMarkersMax, facecolor=colors_comb(iLocComb, :); else, facecolor='w'; end
                    facecolor = 'w';
                    plot(IV_allBins(iBin), data_allBins(iBin), 'o', 'markeredgecolor', colors_comb(iLocComb, :), ...
                        'markerfacecolor', facecolor, 'MarkerSize', nTrials_allBins(iBin) / scalingF + 5, 'LineWidth', 1, 'LineStyle', 'none', 'handlevisibility', 'off')
                end

                % xlim([,150])
                ylim([0, 1])
                ylabel(namesMetrics{iMetric})
                xlabel('Internal variable')
                if isubj == 1, legend(namesModelB(iModelB_all), 'Location', 'best'), end

            end % iModelB

            % print estimates
            text(100, .2, str_est, 'FontSize', 6)
            title(subjName)
        end % isubj

        sgtitle(sprintf('n=%d [A%d] [L%d] [nBoot=%d] %s', nsubj, iModelA_plot, iLocComb, nBoot, namesMetrics{iMetric}))
        set(findall(gcf, '-property', 'fontsize'), 'fontsize', 12)

        % Save the figure
        saveas(gcf, sprintf('%s/n%d_L%d_A%d_%s.jpg', nameFolder_Fig_NOM_metrics, nsubj, iLocComb, iModelA_plot, namesMetrics{iMetric}))
        close all
    end % iMetric

    fprintf('DONE\n')
end % iLocComb

% end % iModelA

fprintf('\n\n Plotting DONE\n\n')

%% Plot GoF for ModelA x ModelB x Loc
clc
clc, fprintf('\n\n Plotting STARTING\n\n')
% Define folder for saving GoF figures
nameFolder_Fig_NOM_GoF = fullfile(nameFolder_Fig_NOM_Trialwise, 'GoF');
if isempty(dir(nameFolder_Fig_NOM_GoF)), mkdir(nameFolder_Fig_NOM_GoF); end

flag_plotIDVD = 0;

% Goodness-of-fit (GoF) measures to plot
namesGoF = {'AIC-nLL', 'AICc-nLL', 'BIC-nLL'}; % no raw nLL
nGoFs = numel(namesGoF);

% Common y-limit for all GoF measures (already positive deltas)
y_lim = [0, 2e3];

% Loop over GoF measures
for iGoF = 1:nGoFs

    % Extract the relevant GoF tensor (AIC, AICc, BIC based on iGoF)
    GoF_allCond = squeeze(IC_nLL_allCond(:, :, :, :, :, iGoF));

    figure('Position', [0 200 2e3 600])

    % Loop over location combinations
    for iLocComb = iLocComb_all

        % data: [ModelA x ModelB x LocComb x Metric x Subj x Iteration]
        GoF = GoF_allCond(:, iModelB_all, iLocComb, :, :);
        data_delta = nan(size(GoF));

        % For each subject and each bootstrap, subtract the minimum across models
        for isubj = 1:nsubj
            parfor iBoot = 1:nBoot
                d = GoF(:, :, :, isubj, iBoot);
                d_min = min(d(:));
                data_delta(:, :, :, isubj, iBoot) = d - d_min;
            end
        end

        % Subplot index for this location
        subplotIdx = find(iLocComb == iLocComb_all);
        subplot(1, numel(iLocComb_all), subplotIdx); hold on

        % First: median across bootstraps within each subject
        GoF_delta_med_allSubj = getCI(data_delta, 1, 5);

        % Then: mean and SEM across subjects
        [GoF_delta_ave, ~, ~, GoF_delta_sem] = getCI(GoF_delta_med_allSubj, 2, 3);

        % Group-level summary (ModelA x ModelB)
        for iModelA = iModelA_all

            switch iModelA
                case 1 % RC-derived template: bar
                    bar(1:nModelsB, GoF_delta_ave(iModelA, :), 'FaceColor', 'w', 'EdgeColor', colors_comb(iLocComb, :), 'LineWidth', 2);

                case 2 % Ideal template: circle markers
                    plot(1:nModelsB, GoF_delta_ave(iModelA, :), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors_comb(iLocComb, :), 'LineWidth', 2);
            end

            % Error bars (SEM across subjects)
            errorbar(1:nModelsB, GoF_delta_ave(iModelA, :), GoF_delta_sem(iModelA, :), '.', 'CapSize', 0, 'Color', colors_comb(iLocComb, :), 'LineWidth', 2);
        end

        % Individual-subject data (for ModelA = 1 only)
        if flag_plotIDVD
            iModelA_idvd = 1;
            for isubj = 1:nsubj

                if isubj <= nMarkersMax, faceColor = 'w'; else, faceColor = ones(1, 3) / 2; end

                plot(1:nModelsB, GoF_delta_med_allSubj(iModelA_idvd, :, isubj), ['-', markers_allSubj{isubj}], 'Color', ones(1, 3) / 2, 'MarkerFaceColor', faceColor);
            end
        end
        % Axes formatting
        xticks(1:nModelsB);
        xticklabels(namesModelB(iModelB_all));
        xtickangle(45);

        ylim(y_lim);
        title(sprintf('L%d', iLocComb));
    end % iLocComb

    % Figure-level formatting and save
    set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 2);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

    sgtitle(sprintf('n=%d %s', nsubj, namesGoF{iGoF}));

    saveas(gcf, fullfile(nameFolder_Fig_NOM_GoF, sprintf('n%d_%s.jpg', nsubj, namesGoF{iGoF})));

end % iGoF

close all

fprintf('\n\n Plotting DONE\n\n')

%% ANOVA on models
iModelB_ANOVA_all = [1,3];
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

p_modelB_allBoot = nan(nBoot, 1);
p_loc_allBoot = p_modelB_allBoot;
p_interaction_allBoot = p_modelB_allBoot;
t_allBoot = p_modelB_allBoot;
p_ttest_allBoot = p_modelB_allBoot;
diff_allBoot = p_modelB_allBoot;
for iBoot= 1:nBoot
    GoF_perBoot = squeeze(GoF_allCond(1, iModelB_ANOVA_all, iLocComb_all, :, iBoot));

    [~, tbl] = print_nANOVA({'NOMmodel', 'Loc'}, GoF_perBoot(:), {IV_NOMmodel(:), IV_Loc(:)}, nsubj);
    p_modelB_allBoot(iBoot) = tbl{2,6};
    p_loc_allBoot(iBoot) = tbl{3,6};
    p_interaction_allBoot(iBoot) = tbl{4,6};
    if nModelB_ANOVA==2
        % conduct t-test
        [~, p_ttest, ~, stats] = ttest(GoF_perBoot(1, :), GoF_perBoot(2, :));
        t_allBoot(iBoot) = stats.tstat;
        p_ttest_allBoot(iBoot) = p_ttest;
        diff_allBoot(iBoot) = mean(GoF_perBoot(1, :)-GoF_perBoot(2, :));
    end
end

% Print median and CI of p-values
clc
fprintf('\n\nModelB ANOVA: B%d-%s vs. B%d-%s\n', iModelB_ANOVA_all(1), namesModelB{iModelB_ANOVA_all(1)}, iModelB_ANOVA_all(2), namesModelB{iModelB_ANOVA_all(2)})
for iDataset=1:3
    switch iDataset
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
    [t_med, t_lb, t_ub] = getCI(t_allBoot, 1, 1);
    [p_med, p_lb, p_ub] = getCI(p_ttest_allBoot, 1, 1);
    [diff_med, diff_lb, diff_ub] = getCI(diff_allBoot, 1, 1);
    fprintf('   t=%.3f [CI=%.3f, %.3f]\n', t_med, t_lb, t_ub);
    fprintf('   p=%.3f [CI=%.3f, %.3f]\n', p_med, p_lb, p_ub);
    fprintf('   diff=%.3f [CI=%.3f, %.3f]\n', diff_med, diff_lb, diff_ub);
end


%% Compare parameters across locations
% Define folder for saving comparison figures
nameFolder_Fig_NOM_compParams = fullfile(nameFolder_Fig_NOM_Trialwise, 'compParams');
if isempty(dir(nameFolder_Fig_NOM_compParams))
    mkdir(nameFolder_Fig_NOM_compParams);
end

nParamsMax = 3; % max three params

% Layout for subplots (row = location-pair index, column = parameter index)
isubplots = reshape(1:12, 4, 3)'; % 3x4 mapping

for iModelB = iModelB_plot %iModelB_all % just plot the best model

    % Number of parameters for this Model B
    nParamsB = numel(namesModelBparams{iModelB});

    figure('Position', [0 0 2e3 2e3]);
    % iiLocComb = 1;

    % Compare parameters between specific locations (here: [6, 7])
    % To change locations, modify the cell contents below.

    iplots = reshape(1:nGroups*nParamsMax, [nParamsMax, nGroups])';

    for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
        iLocPair_all = iLocGroups_all{iGroup};
        % iLoc = iLoc_{1}; % numeric vector of location indices (e.g., [6 7])

        for iParam = 1:nParamsB
            subplot(nGroups, nParamsMax, iplots(iGroup, iParam)); hold on;

            % params_allSubj: [nLoc (here 2) x nSubj x nBoot]
            params_allSubj = squeeze(params_allCond(iModelA_plot, iModelB, iLocPair_all, :, :, iParam));
            params_allSubj_med = getCI(params_allSubj, 1, 3); % median across bootstraps
            [params_ave, ~, ~, params_sem] = getCI(params_allSubj_med, 2, 2); % mean and SEM across subjects

            % Group-level bars with error bars
            for iiLoc = 1:length(iLocPair_all)
                bar(iiLoc, params_ave(iiLoc), ...
                    'FaceColor', 'w', ...
                    'EdgeColor', colors_comb(iLocPair_all(iiLoc), :), ...
                    'BarWidth', 0.5);
                errorbar(iiLoc, params_ave(iiLoc), params_sem(iiLoc), ...
                    'CapSize', 0, ...
                    'Color', colors_comb(iLocPair_all(iiLoc), :));
            end

            % Individual subject traces (paired between the two locations)
            for isubj = 1:nsubj
                if isubj <= nMarkersMax
                    faceColor = 'w';
                else
                    faceColor = ones(1, 3) / 2;
                end
                % Assumes exactly 2 locations → x = [1.2, 1.8]
                plot([1.2, 1.8], params_allSubj_med(:, isubj), ...
                    ['-', markers_allSubj{isubj}], ...
                    'Color', ones(1, 3) / 2, ...
                    'MarkerFaceColor', faceColor, ...
                    'MarkerSize', 10);
            end

            % Paired t-tests across bootstraps
            t_allBoot = nan(1, nBoot);
            p_allBoot = nan(1, nBoot);

            for iBoot = 1:nBoot
                a = params_allSubj(1, :, iBoot); % location 1
                b = params_allSubj(2, :, iBoot); % location 2

                [~, p, ~, stats] = ttest(a, b);
                t_allBoot(iBoot) = stats.tstat;
                p_allBoot(iBoot) = p;
            end

            xticks(1:numel(iLocPair_all));
            xticklabels(namesLocComb(iLocPair_all));

            % CI of t and p across bootstraps
            [t_med, t_lb, t_ub] = getCI(t_allBoot, 1, 2);
            [p_med, p_lb, p_ub] = getCI(p_allBoot, 1, 2);

            title(sprintf('%s\nt = %.2f [%.2f, %.2f], p = %.3f [%.3f, %.3f]', ...
                namesModelBparams{iModelB}{iParam}, ...
                t_med, t_lb, t_ub, ...
                p_med, p_lb, p_ub));

        end % iParam

    end % iLoc_

    sgtitle(sprintf('A%d [%s] B%d [%s]\n[n = %d] [ni = %d]', iModelA_plot, namesModelA{iModelA_plot}, iModelB, namesModelB{iModelB}, nsubj, nBoot));
    set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 1.5);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);

    % Save the figure
    saveas(gcf, fullfile(nameFolder_Fig_NOM_compParams, sprintf('n%d_A%dB%d.jpg', nsubj, iModelA_plot, iModelB)));

end % iModelB

close all;

%% Template for ANOVA
% for ii=1:ni
% GoF_all_perIte = squeeze(GoF_all(:, :, :, :, :, ii, :)); % select data of each iteraction of bootstrap
% % nModelsA x nModelsB x nLoc x nMetrics x nsubj x nGoF
%
% % 5-way ANOVA: GoF x metric x model A x model B
% % get index
% indModelA = nan(size(GoF_all_perIte));
% indModelB = indModelA;
% indLocComb = indModelA;
% indMetric = indModelA;
% indGoF = indModelA;
%
% for iModelA=1:nModelsA
% for iModelB=1:nModelsB
% for iLocComb=iLocComb_all
% for iMetric=1:nMetrics
% for iGoF=1:nGoFs
% indModelA(iModelA, iModelB, iLocComb, iMetric, :, iGoF) = ones(nsubj, 1) * iModelA;
% indModelB(iModelA, iModelB, iLocComb, iMetric, :, iGoF) = ones(nsubj, 1) * iModelB;
% indLocComb(iModelA, iModelB, iLocComb, iMetric, :, iGoF) = ones(nsubj, 1) * iLocComb;
% indMetric(iModelA, iModelB, iLocComb, iMetric, :, iGoF) = ones(nsubj, 1) * iMetric;
% indGoF(iModelA, iModelB, iLocComb, iMetric, :, iGoF) = ones(nsubj, 1) * iGoF;
% end
% end
% end
% end
% end
%
% GoF_all_perIte = GoF_all_perIte(:, :, iLocComb_all, :, :, :, :);
% indModelA = indModelA(:, :, iLocComb_all, :, :, :, :);
% indModelB = indModelB(:, :, iLocComb_all, :, :, :, :);
% indLocComb = indLocComb(:, :, iLocComb_all, :, :, :, :);
% indMetric = indMetric(:, :, iLocComb_all, :, :, :, :);
% indGoF = indGoF(:, :, iLocComb_all, :, :, :, :);
%
% % conduct ANOVA
% str_ANOVA = print_nANOVA({'ModelA', 'ModelB', 'LocComb', 'Metric', 'GoF'}, GoF_all_perIte(:), {indModelA(:), indModelB(:), indLocComb(:), indMetric(:), indGoF(:)}, nsubj)
% end

%% Compare across ModelAs x ModelBs
% % IC_nLL_allCond: nModelsA x nModelsB x nLoc8 x nMetrics x nsubj x nBoot x nGoF
% namesGoF = {'AIC-nLL', 'AICc-nLL', 'BIC-nLL'}; % no raw nLL
% nGoFs = numel(namesGoF);
%
% % 2-way ANOVA: Model A x Model B; for each Loc x Metric x GoF
% namesANOVA_Var = {'A', 'B', 'AxB'}; % 1 = main effect of ModelA, 2 = main effect of ModelB, 3 = interaction
% nANOVA_Var = numel(namesANOVA_Var);
%
% % Preallocate:
% % p_allB_allCond(loc, metric, gof, bootstrap, anovaVar)
% p_allB_allCond = nan(nLocComb8, nMetrics, nGoFs, nBoot, nANOVA_Var);
%
% for iGoF = [1,3]%1:nGoFs
%
%     for iLocComb = iLocComb_all
%
%         for iMetric = 1:nMetrics
%
%             parfor iBoot = 1:nBoot
%
%                 % Select data for this GoF, location, metric, and bootstrap
%                 % GoF_all_perIte: [nModelsA x nModelsB x nsubj]
%                 GoF_all_perIte = squeeze(IC_nLL_allCond(:, :, iLocComb, :, iBoot, iGoF));
%
%                 % Build factor indices for ANOVA (same size as GoF_all_perIte)
%                 indModelA = nan(size(GoF_all_perIte));
%                 indModelB = nan(size(GoF_all_perIte));
%
%                 for iModelA = iModelA_all
%                     for iModelB = iModelB_all
%                         indModelA(iModelA, iModelB, :) = iModelA * ones(nsubj, 1);
%                         indModelB(iModelA, iModelB, :) = iModelB * ones(nsubj, 1);
%                     end
%                 end
%
%                 % Conduct 2-way ANOVA: factors = ModelA, ModelB
%                 [~, tbl] = print_nANOVA( ...
%                     {'ModelA', 'ModelB'}, ...
%                     GoF_all_perIte(:), ...
%                     {indModelA(:), indModelB(:)}, ...
%                     nsubj, ...
%                     0);
%
%                 % p-values: row 2 = A, row 3 = B, row 4 = A×B, col 7 = p
%                 p_allB_allCond(iLocComb, iMetric, iGoF, iBoot, :) = ...
%                     [tbl{2, 7}, tbl{3, 7}, tbl{4, 7}];
%
%             end % iBoot
%
%         end % iMetric
%
%     end % iLocComb
%
% end % iGoF
%
% % Plotting
% nameFolder_Fig_NOM_ANOVA_AxB = sprintf('%s/ANOVA/ANOVA_AxB', nameFolder_Fig_NOM_Trialwise);
% if isempty(dir(nameFolder_Fig_NOM_ANOVA_AxB)), mkdir(nameFolder_Fig_NOM_ANOVA_AxB), end
%
% for iGoF = 1:nGoFs
%     figure('Position', [0 200 numel(iLocComb_all) * 333 1e3])
%
%     for iLocComb = iLocComb_all
%
%         % p_allB_allCond(loc, metric, gof, bootstrap, anovaVar)
%         [p_med, p_lb, p_ub] = getCI(p_allB_allCond(iLocComb, iMetric, iGoF, :, :), 1, 4);
%
%         subplot(nGoFs, numel(iLocComb_all), find(iLocComb == iLocComb_all)); hold on
%
%         bar(1:nANOVA_Var, p_med, ...
%             'EdgeColor', colors_comb(iLocComb, :), ...
%             'FaceColor', 'w', ...
%             'LineWidth', 1.5);
%
%         errorbar(1:nANOVA_Var, p_med, p_lb, p_ub, ...
%             '.', ...
%             'Color', colors_comb(iLocComb, :), ...
%             'CapSize', 0, ...
%             'LineWidth', 1.5);
%
%         xticks(1:nANOVA_Var);
%         xticklabels(namesANOVA_Var);
%         ylabel('p-value')
%
%         yline(0.05, 'k--');
%         xlim([0, nANOVA_Var + 1]);
%         ylim([0, max(p_ub)]);
%
%         % text(0, 0.09, sprintf('%s - L%d', namesGoF{iGoF}, iLocComb));
%
%         title(sprintf('L%d', iLocComb));
%     end % iLocComb
%
%     % Figure-level formatting and save
%     set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 2);
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);
%
%     sgtitle(sprintf('%s', namesGoF{iGoF}));
%
%     % saveas(gcf, fullfile(nameFolder_Fig_NOM_ANOVA_AxB, ...
%     %     sprintf('n%d_%s.jpg', nsubj, namesGoF{iGoF})));
% end
% close all

%% Plot pred vs. meas metrics (scatter plots)
% wd_errorbar = 2;
% close all
%
% for iModelA = 1%iModelA_all
%
%     for iModelB = iModelB_all
%
%         % Define folder for saving figures for each model A and model B
%         nameFolder_Fig_NOM_metrics = sprintf('%s/A%d', nameFolder_Fig_NOM_Trialwise, iModelA);
%         if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end
%
%         for iLocComb = iLocComb_all
%
%             figure('Position', [0 200 2e3 1e3])
%
%             for iMetric = 1:nMetrics
%                 data_allSubj = squeeze(metric_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :));
%                 pred_allSubj = squeeze(pred_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :));
%
%                 data_allSubj_med = getCI(data_allSubj, 1, 2);
%                 pred_allSubj_med = getCI(pred_allSubj, 1, 2);
%
%                 if ~flag_subjIsHuman
%                     data_allSubj_med = getCI(data_allSubj, 1, 1);
%                     pred_allSubj_med = getCI(pred_allSubj, 1, 1);
%                 end
%
%                 [data_ave, ~, ~, data_sem] = getCI(data_allSubj_med, 2, 2);
%                 [pred_ave, ~, ~, pred_sem] = getCI(pred_allSubj_med, 2, 2);
%
%                 if ~flag_subjIsHuman
%                     data_ave = data_allSubj_med;
%                     data_sem = zeros(size(data_ave));
%                     pred_ave = pred_allSubj_med;
%                     pred_sem = data_sem;
%                 end
%
%                 % correlation between measurement and prediction & t-tests
%                 ind_subj = repmat((1:nsubj)', 1, nBins);
%                 r_allBoot = nan(nBoot, 1);
%                 p_corr_allBoot = r_allBoot;
%                 t_allBoot = r_allBoot;
%                 p_ttest_allBoot = r_allBoot;
%
%                 for iBoot = 1:nBoot
%
%                     if flag_subjIsHuman
%                         x = squeeze(data_allSubj(:, iBoot, :)); y = squeeze(pred_allSubj(:, iBoot, :));
%                     else
%                         % x = squeeze(data_allSubj(:, iBoot, :)); y = squeeze(pred_allSubj(:, iBoot, :));
%                     end
%
%                     [r, p] = partialcorr(x(:), y(:), ind_subj(:));
%                     r_allBoot(iBoot) = r;
%                     p_corr_allBoot(iBoot) = p;
%
%                     [~, p, ~, stats] = ttest(x(:), y(:));
%                     t_allBoot(iBoot) = stats.tstat;
%                     p_ttest_allBoot(iBoot) = p;
%                 end
%
%                 subplot(1, nMetrics, iMetric), hold on
%                 axis square
%                 for isubj = 1:nsubj
%                     plot(data_allSubj_med(isubj, :), pred_allSubj_med(isubj, :), 'o-', 'color', ones(1, 3) / 2, 'HandleVisibility', 'off')
%                     if isubj <= nMarkersMax, facecolor = 'w'; else, facecolor = ones(1, 3) / 2; end
%                     plot(data_allSubj_med(isubj, 1), pred_allSubj_med(isubj, 1), markers_allSubj{isubj}, 'markeredgecolor', 'k', 'markerfacecolor', facecolor)
%                     plot(data_allSubj_med(isubj, end), pred_allSubj_med(isubj, end), markers_allSubj{isubj}, 'markeredgecolor', 'k', 'markerfacecolor', facecolor, 'HandleVisibility', 'off')
%                 end
%
%                 errorbar(data_ave, pred_ave, data_sem, 'horizontal', 'color', colors_comb(iLocComb, :), 'CapSize', 0, 'LineWidth', wd_errorbar, 'HandleVisibility', 'off')
%                 errorbar(data_ave, pred_ave, pred_sem, 'vertical', 'color', colors_comb(iLocComb, :), 'CapSize', 0, 'LineWidth', wd_errorbar, 'HandleVisibility', 'off')
%                 plot(data_ave, pred_ave, 'o', 'markerfacecolor', colors_comb(iLocComb, :), 'markeredgecolor', 'w', 'HandleVisibility', 'off')
%                 plot([0, 1], [0, 1], 'k-', 'HandleVisibility', 'off')
%                 legend(subjList, 'Location', 'best', 'NumColumns', 3)
%                 xticks(linspace(0, 1, 5)), yticks(linspace(0, 1, 5))
%                 xlabel(sprintf('Measured %s', namesMetrics{iMetric}))
%                 ylabel(sprintf('Predicted %s', namesMetrics{iMetric}))
%
%                 % Group-averaged GoF
%                 % [R2_ave, ~, ~, R2_sem] = getCI(getCI(R2_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6), 2, 1);
%                 % [R2_w_ave, ~, ~, R2_w_sem] = getCI(getCI(R2_w_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6), 2, 1);
%                 % [GoF_delta_ave, ~, ~, GoF_delta_sem] = getCI(getCI(SSE_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6), 2, 1);
%                 % [SSE_w_ave, ~, ~, SSE_w_sem] = getCI(getCI(SSE_w_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6), 2, 1);
%                 [IC_nLL_ave, ~, ~, IC_nLL_sem] = getCI(getCI(IC_nLL_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :), 1, 6), 2, 1);
%                 % IC_SSE_ave = getCI(getCI(IC_SSE_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :), 1, 6), 2, 1);
%                 % IC_SSE_w_ave = getCI(getCI(IC_SSE_w_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :), 1, 6), 2, 1);
%
%                 % if ~flag_subjIsHuman
%                 % [R2_ave, ~, ~, R2_sem] = getCI(R2_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6);
%                 % [R2_w_ave, ~, ~, R2_w_sem] = getCI(R2_w_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6);
%                 % [GoF_delta_ave, ~, ~, GoF_delta_sem] = getCI(SSE_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6);
%                 % [SSE_w_ave, ~, ~, SSE_w_sem] = getCI(SSE_w_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6);
%                 % IC_SSE_ave = getCI(IC_SSE_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :), 1, 6);
%                 % IC_SSE_w_ave = getCI(IC_SSE_w_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :), 1, 6);
%                 % end
%
%                 % title(sprintf('IC_nLL (%s):%.2f (%.2f) || R2: %.2f (%.2f) || %.2f (%.2f)\nSSE: %.2f (%.2f) || %.2f (%.2f)\n%s: %.0f || %.0f\nPartial corr: r = %.2f, p=%.3f\nttest: t = %.2f, p=%.3f', ...
%                 % namesIC{iIC_plot}, IC_nLL_ave(iIC_plot), IC_nLL_sem(iIC_plot), ...
%                 % R2_w_ave, R2_w_sem, R2_ave, R2_sem, ...
%                 % SSE_w_ave, SSE_w_sem, GoF_delta_ave, GoF_delta_sem, ...
%                 % namesIC{iIC_plot}, IC_SSE_w_ave(iIC_plot), IC_SSE_ave(iIC_plot), ...
%                 % nanmedian(r_allB), nanmedian(p_corr_allB), ...
%                 % nanmedian(t_allB), nanmedian(p_ttest_allB)))
%
%                 title(sprintf('Partial corr: r = %.2f, p=%.3f\nttest: t = %.2f, p=%.3f', ...
%                     nanmedian(r_allBoot), nanmedian(p_corr_allBoot), ...
%                     nanmedian(t_allBoot), nanmedian(p_ttest_allBoot)))
%
%             end % iMetric
%
%             sgtitle(sprintf('[A%dB%d] [L%d] [ni=%d] nsubj=%d', iModelA, iModelB, iLocComb, nBoot, nsubj))
%             set(findall(gcf, '-property', 'linewidth'), 'linewidth', 1.5)
%             set(findall(gcf, '-property', 'fontsize'), 'fontsize', 20)
%
%             % Save the figure
%             saveas(gcf, sprintf('%s/n%d_L%d_B%d_pred_vs_meas.jpg', nameFolder_Fig_NOM_metrics, nsubj, iLocComb, iModelB))
%             %%%%%%%%%%
%         end % iLocComb
%
%         close all
%     end % iModelB
%
% end % iModelA
