% SX_analysis5_NOM_Trialwise.m
% Author: Shutian Xue
% Purpose: This script compiles, analyzes and plots boostrapped data for each observer saved Data_NOM_Trialwise

clc, clear all, close all

% addpath(genpath('fxn_analysis_RC_v2')) % manually add to save time

% Setting parameters for the analysis
%----------------
SX_RC1_setting
%----------------

%%%%%% Should copy from OOD_sim %%%%%%%%%%
noiseCST_all = [0, .1, .2, .5]; % Noise contrast sensitivity thresholds
gaborCST_all = [.1, .5]; % Gabor contrast sensitivity thresholds
nTrials_all = [5000]; % Number of trials per condition
noiseP_all = [0, 0.1, 0.2]; % Proportion of noise trials
iModelA_sim_all = [1,3]; % 1=core model, 2=randomize template, 3=use IO template
iModelB_sim_all = [4];
%               1 = estimate lapse rate, additive noise (SDadd) and criterion;
%               2 = estimate lapse rate, multiplicative noise (Nmul) and criterion;
%               3 = estimate lapse rate and criterion;
%               4 = estimate additive noise (SDadd) and criterion;
%               5 = estimate multiplicative noise (Nmul) and criterion;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nIterations = 100; % Number of iterations (100 for raw IVs, 99 for transformed IVs)
NOM_mode = 2; % 1= aggregate model, 2 = trial-wise model
flag_subjIsHuman = 0;

nRows=4; nCols=6;

iModelA_all = iModelA_sim_all;
iModelB_all = iModelB_sim_all;
% Define subject list
if flag_subjIsHuman
    subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC', 'SR'};
    iLocComb_all = [1, 8, 6, 7, 5, 3]; % combination of locations

else
    % Find folders with names starting with 'IO' in Data/Data_NOM_Trialwise_2929
    % Filter for directories whose names start with 'IO'
    nameDir_IO_all = dir(sprintf('%s/IO*', nameFolder_Data_NOM_Trialwise));
    nIOs = length(nameDir_IO_all);
    subjList = cell(1, nIOs);

    for iIO = 1:nIOs
        subjList{iIO} = nameDir_IO_all(iIO).name; % Store the names of IO folders
    end
    assert(nRows * nCols >= nIOs, 'Not enough subplots for the number of subjects!');
    iLocComb_all = 1;, % For testing with IO data only
end

nsubj = length(subjList);

% Define the folder structure based on the mode (NOM_mode)
% switch NOM_mode
%     case 1, nameFolder_Data_NOM_Trialwise = 'Data_NOM_aggregate';
%     case 2, nameFolder_Data_NOM_Trialwise = 'Data_NOM_Trialwise';
% end

% Names of the performance metrics being analyzed
namesMetrics = {'pC', 'pYES', 'pA'};
nMetrics = length(namesMetrics);
nModelsA = 3; % do not use length()!!
nModelsB = 5; % do not use length()!!
nLocComb8 = 8; % do not use length()!!
nBins = 10; % Number of bins for the analysis
NOM_mode = 2; % 1 for aggregate data, 2 for trial-wise data
iIC_plot = 3; % % Which information criterion to plot (1 = AIC, 2 = AICc, 3 = BIC)

% Function handles for calculating information criteria (IC)
getAIC_SSE = @(SSE, nParams, nData) nData * log(SSE / nData) + 2 * nParams;
getAICc_SSE = @(SSE, nParams, nData) nData * log(SSE / nData) + 2 * nParams + 2 * nParams * (nParams + 1) / (nData - nParams - 1);
getBIC_SSE = @(SSE, nParams, nData) nData * log(SSE / nData) + nParams * log(nData);

%% organize data from boostrapped data
% Preallocate arrays for storing results across all conditions
IV_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIterations, nBins);
nTrials_allCond = IV_allCond;
data_allCond = nan(nModelsA, nModelsB, nLocComb8, nMetrics, nsubj, nIterations, nBins);
pred_allCond = data_allCond;
nLL_allCond = nan(nModelsA, nModelsB, nLocComb8, nMetrics, nsubj, nIterations);
params_allCond = nan(nModelsA, nModelsB, nLocComb8, nsubj, nIterations, 4); % Pre-allocate for max params

R2_allCond = nan(nModelsA, nModelsB, nLocComb8, nMetrics, nsubj, nIterations);
R2_w_allCond = R2_allCond;
SSE_allCond = R2_allCond;
SSE_w_allCond = SSE_allCond;
IC_SSE_allCond = nan(nModelsA, nModelsB, nLocComb8, nMetrics, nsubj, nIterations, 3); % AIC, AICc, BIC (3)
IC_SSE_w_allCond = IC_SSE_allCond;

% The main analysis loop for all conditions for each subject (Models x Loc x Metrics)
for iModelA = iModelA_all %1:nModelsA
    fprintf('\nModelA #%d: ', iModelA)

    for iModelB = iModelB_all %1:nModelsB
        fprintf('\n  ModelB #%d: ', iModelB)
        nParams = length(namesParamsModel_all{iModelB});

        for iLocComb = iLocComb_all
            fprintf('\n    L%d: ', iLocComb)

            for iMetric = 1:nMetrics
                fprintf(' Metric%d: ', iMetric)

                for isubj = 1:nsubj

                    subjName = subjList{isubj};

                    if flag_subjIsHuman
                        % nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, nameIO);
                        nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocComb); % MUST be the same as OOD_NOM_Trialwise_Est.m, Line 61
                    else
                        % nameFolder_NOM_save = sprintf('%s/ORI%dSF%d/%s/L%d', nameFolder_NOM0, nORI, nSF, subjName, iLocComb);
                        nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, subjName); % MUST be the same as OOD_NOM_Trialwise_Est.m, Line 68
                    end

                    % nameFile_Est = sprintf('%s/n%d_A%dB%d', nameFolder_NOM_Trialwise, ni, iModelA, iModelB);
                    nameFile_Est = sprintf('%s/n%d_A%dB%d', nameFolder_NOM_save, nIterations, iModelA, iModelB); % MUST be the same as OOD_NOM_Trialwise_Est.m, Line 77
                    load(nameFile_Est, 'pred_metrics_allB', 'params_est_allB', 'nLL_allB')

                    % Pre-allocate temporary arrays for this subject and model
                    IV_allB = nan(nIterations, nBins); % IV for each bin
                    nTrials_allB = IV_allB; % Trial count for each bin
                    data_allB = IV_allB; % Observed data for each bin
                    pred_allB = IV_allB; % Predicted data for each bin

                    for iIteration = 1:nIterations
                        % Extract IV per bin
                        IV_allB(iIteration, :) = pred_metrics_allB{iIteration}.metrics.IV_allBins;

                        % Extract the number of trials per bin
                        nTrials_allB(iIteration, :) = pred_metrics_allB{iIteration}.metrics.nTrials_allBins;

                        % Extract metrics (data and predictions) per bin
                        pred_metrics = pred_metrics_allB{iIteration}.metrics; % Extract once for faster access

                        switch namesMetrics{iMetric}
                            case 'pC'
                                data_allB(iIteration, :) = pred_metrics.pC_data_allBins;
                                pred_allB(iIteration, :) = pred_metrics.pC_pred_allBins;
                            case 'pA'
                                data_allB(iIteration, :) = pred_metrics.pA_data_allBins;
                                pred_allB(iIteration, :) = pred_metrics.pA_pred_allBins;
                            case 'pYES'
                                data_allB(iIteration, :) = pred_metrics.pYES_data_allBins;
                                pred_allB(iIteration, :) = pred_metrics.pYES_pred_allBins;
                        end

                        % Calculate goodness-of-fit and information criteria for both weighted and unweighted cases
                        for flag_weighted = [0, 1]

                            if flag_weighted, weights = nTrials_allB(iIteration, :);
                            else , weights = nan; % in the function, vector of ones will be created
                            end

                            % Calculate SSE and R2
                            [R2, SSE] = getR2(data_allB(iIteration, :), pred_allB(iIteration, :), weights);

                            % Calculate information critertion based on SSE
                            nData = sum(nTrials_allB(iIteration, :));
                            IC_SSE = [getAIC_SSE(SSE, nParams, nData), getAICc_SSE(SSE, nParams, nData), getBIC_SSE(SSE, nParams, nData)];

                            % Store results in the preallocated arrays
                            if flag_weighted
                                R2_w_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, iIteration) = R2;
                                SSE_w_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, iIteration) = SSE;
                                IC_SSE_w_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, iIteration, :) = IC_SSE;
                            else
                                R2_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, iIteration) = R2;
                                SSE_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, iIteration) = SSE;
                                IC_SSE_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, iIteration, :) = IC_SSE;
                            end

                        end % for flag_weighted=[0,1]

                    end % end of iTeration

                    % Store results
                    IV_allCond(iModelA, iModelB, iLocComb, isubj, :, :) = IV_allB;
                    nTrials_allCond(iModelA, iModelB, iLocComb, isubj, :, :) = nTrials_allB;
                    data_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, :, :) = data_allB;
                    pred_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, :, :) = pred_allB;
                    nLL_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, :) = nLL_allB;
                    params_allCond(iModelA, iModelB, iLocComb, isubj, :, 1:nParams) = params_est_allB;
                end % isubj

            end % iMetric

        end % iLocComb

    end % iModelB

end % iModelA

fprintf('\n\n========= ORGANIZATION DONE =============\n\n')

% Save the organized data for all subjects
save(sprintf('%s/n%d_n%d', nameFolder_Data_NOM_Trialwise, nsubj, nIterations), '*_allCond')

%% Plot metrics vs. IV for each idvd
lineStyle_all = {'-', '--', ':', '-.', '-', '--', ':', '-.'};

% Define folder for saving figures
nameFolder_Fig_NOM_Trialwise = sprintf('%s/NOM_Trialwise_%d%d', nameFolder_Figures, nORI, nSF);
if isempty(dir(nameFolder_Fig_NOM_Trialwise)), mkdir(nameFolder_Fig_NOM_Trialwise), end

for iModelA = iModelA_all

    % Define folder for saving figures for each model A and model B
    nameFolder_Fig_NOM_AB = sprintf('%s/A%d', nameFolder_Fig_NOM_Trialwise, iModelA);
    if isempty(dir(nameFolder_Fig_NOM_AB)), mkdir(nameFolder_Fig_NOM_AB), end

    for iLocComb = iLocComb_all
        scalingF = 80; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end

        %%%%%%%%%%%%
        for iMetric = 1:nMetrics
            figure('Position', [0 0 2e3 2e3])

            for isubj = 1:nsubj
                subplot(nRows, nCols, isubj), hold on

                subjName = subjList{isubj};

                % Define folder for loading NOM file
                % nameFolder_NOM_save = sprintf('%s/ORI%dSF%d/%s/L%d', nameFolder_Data_NOM_Trialwise, nORI, nSF, subjName, iLocComb);
                if flag_subjIsHuman
                    nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocComb); % To Save results
                else
                    nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, subjName); % To Save results
                end

                str_est = [];

                for iModelB = iModelB_all

                    % Define file names for loading data
                    nameFile_Est = sprintf('%s/n%d_A%dB%d', nameFolder_NOM_save, nIterations, iModelA, iModelB); % MUST be the same as OOD_NOM_Trialwise_Est.m, Line 77

                    % Compute medians and CIs
                    [IV_allBins, ~, ~, IV_allBins_neg, IV_allBins_pos] = getCI(IV_allCond(iModelA, iModelB, iLocComb, isubj, :, :), 1, 5);
                    [nTrials_allBins, nTrials_allBins_lb, nTrials_allBins_ub] = getCI(nTrials_allCond(iModelA, iModelB, iLocComb, isubj, :, :), 1, 5);
                    [data_allBins, ~, ~, data_allBins_neg, data_allBins_pos] = getCI(data_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, :, :), 1, 6);
                    [pred_allBins, pred_allBins_lb, pred_allBins_ub] = getCI(pred_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, :, :), 1, 6);
                    [params_allBins, params_allBins_lb, params_allBins_ub] = getCI(params_allCond(iModelA, iModelB, iLocComb, isubj, :, :), 1, 5);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Plot prediction
                    plot(IV_allBins, pred_allBins, [lineStyle_all{iModelB}, 'k'])
                    patch([IV_allBins; flip(IV_allBins)], [pred_allBins_lb; flip(pred_allBins_ub)], ones(1, 3) / 2, 'FaceAlpha', .3, 'linestyle', 'none', 'handlevisibility', 'off')
                    % measurement
                    errorbar(IV_allBins, data_allBins, IV_allBins_neg, IV_allBins_pos, '.', 'horizontal', 'CapSize', 0, 'color', colors_comb(iLocComb, :), 'handlevisibility', 'off')
                    errorbar(IV_allBins, data_allBins, data_allBins_neg, data_allBins_pos, '.', 'vertical', 'CapSize', 0, 'color', colors_comb(iLocComb, :), 'handlevisibility', 'off')

                    yline(.5, 'k--');

                    % Plot averaged data of each bin
                    for iBin = 1:nBins
                        %                         if isubj>nMarkersMax, facecolor=colors_comb(iLocComb, :); else, facecolor='w'; end
                        facecolor = 'w';
                        plot(IV_allBins(iBin), data_allBins(iBin), 'o', 'markeredgecolor', colors_comb(iLocComb, :), ...
                            'markerfacecolor', facecolor, 'MarkerSize', nTrials_allBins(iBin) / scalingF + 5, 'LineWidth', 1, 'LineStyle', 'none', 'handlevisibility', 'off')
                    end

                    %                     xlim([-1,4])
                    ylim([0, 1])
                    ylabel(namesMetrics{iMetric})
                    xlabel('Internal variable')
                    if isubj == 1, legend(namesModelB, 'Location', 'best'), end

                    str_est = [str_est, sprintf('B%d: %s\n', iModelB, num2str(round(params_allBins', 4)))];

                    % print GoF indices
                    %                     [R2_w_med, R2_w_lb, R2_w_ub] = getCI(R2_w_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, :), 1, 6);
                    %                     [R2_med, R2_lb, R2_ub] = getCI(R2_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, :), 1, 6);
                    %                     [SSE_w_med, SSE_w_lb, SSE_w_ub] = getCI(SSE_w_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, :), 1, 6);
                    %                     [SSE_med, SSE_lb, SSE_ub] = getCI(SSE_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, :), 1, 6);
                    %                     IC_SSE_w_med = getCI(IC_SSE_w_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, :, :), 1, 6);
                    %                     IC_SSE_med = getCI(IC_SSE_allCond(iModelA, iModelB, iLocComb, iMetric, isubj, :, :), 1, 6);
                    %
                    %                     title(sprintf('%s\nR2: %.2f [%.2f, %.2f] || %.2f [%.2f, %.2f]\nSSE: %.2f [%.2f, %.2f] || %.2f [%.2f, %.2f]\n%s: %.0f  || %.0f \nEst.: %s', ...
                    %                         subjName, ...
                    %                         R2_w_med, R2_w_lb, R2_w_ub,      R2_med, R2_lb, R2_ub, ...
                    %                         SSE_w_med, SSE_w_lb, SSE_w_ub,      SSE_med, SSE_lb, SSE_ub, ...
                    %                         namesIC{iIC_plot}, IC_SSE_w_med(iIC_plot), IC_SSE_med(iIC_plot), ...
                    %                         num2str(round(params_allBins', 4))))
                    %                     title(sprintf('%s\nEst.: %s', subjName, num2str(round(params_allBins', 4))))
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end % iModelB

                % print estimates
                text(.01, .1, str_est, 'FontSize', 8, 'Units', 'normalized')
                title(subjName)
            end % isubj

            sgtitle(sprintf('[A%d] [L%d] [nBoot=%d] %s', iModelA, iLocComb, nIterations, namesMetrics{iMetric}))
            set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10)

            % Save the figure
            saveas(gcf, sprintf('%s/n%d_L%d_metric_vs_IV_%s.jpg', nameFolder_Fig_NOM_AB, nsubj, iLocComb, namesMetrics{iMetric}))

        end % iMetric

    end % iLocComb

    close all
end % iModelA

%% Plot pred vs. meas metrics
wd_errorbar = 2;
close all

for iModelA = iModelA_all

    for iModelB = iModelB_all

        % Define folder for saving figures for each model A and model B
        nameFolder_Fig_NOM_AB = sprintf('%s/A%d', nameFolder_Fig_NOM_Trialwise, iModelA);
        if isempty(dir(nameFolder_Fig_NOM_AB)), mkdir(nameFolder_Fig_NOM_AB), end

        for iLocComb = iLocComb_all

            figure('Position', [0 200 2e3 400])

            for iMetric = 1:nMetrics
                data_allSubj = squeeze(data_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :));
                pred_allSubj = squeeze(pred_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :));

                data_allSubj_med = getCI(data_allSubj, 1, 2);
                pred_allSubj_med = getCI(pred_allSubj, 1, 2);

                if ~flag_subjIsHuman
                    data_allSubj_med = getCI(data_allSubj, 1, 1);
                    pred_allSubj_med = getCI(pred_allSubj, 1, 1);
                end

                [data_ave, ~, ~, data_sem] = getCI(data_allSubj_med, 2, 1);
                [pred_ave, ~, ~, pred_sem] = getCI(pred_allSubj_med, 2, 1);

                if ~flag_subjIsHuman
                    data_ave = data_allSubj_med;
                    data_sem = zeros(size(data_ave));
                    pred_ave = pred_allSubj_med;
                    pred_sem = data_sem;
                end

                % correlation between measurement and prediction & t-tests
                ind_subj = repmat((1:nsubj)', 1, nBins);
                r_allB = nan(nIterations, 1);
                p_corr_allB = r_allB;
                t_allB = r_allB;
                p_ttest_allB = r_allB;

                for iIteration = 1:nIterations

                    if ~~flag_subjIsHuman
                        x = squeeze(data_allSubj(iIteration, :)); y = squeeze(pred_allSubj(iIteration, :));
                    else , x = squeeze(data_allSubj(:, iIteration, :)); y = squeeze(pred_allSubj(:, iIteration, :));
                    end

                    [r, p] = partialcorr(x(:), y(:), ind_subj(:));
                    r_allB(iIteration) = r;
                    p_corr_allB(iIteration) = p;

                    [~, p, ~, stats] = ttest(x(:), y(:));
                    t_allB(iIteration) = stats.tstat;
                    p_ttest_allB(iIteration) = p;
                end

                subplot(1, nMetrics, iMetric), hold on

                for isubj = 1:nsubj
                    plot(data_allSubj_med(isubj, :), pred_allSubj_med(isubj, :), '-', 'color', ones(1, 3) / 2, 'HandleVisibility', 'off')
                    if isubj <= nMarkersMax, facecolor = 'w'; else, facecolor = ones(1, 3) / 2; end
                    plot(data_allSubj_med(isubj, 1), pred_allSubj_med(isubj, 1), markers_allSubj{isubj}, 'markeredgecolor', 'k', 'markerfacecolor', facecolor)
                    plot(data_allSubj_med(isubj, end), pred_allSubj_med(isubj, end), markers_allSubj{isubj}, 'markeredgecolor', 'k', 'markerfacecolor', facecolor, 'HandleVisibility', 'off')
                end

                errorbar(data_ave, pred_ave, data_sem, 'horizontal', 'color', colors_comb(iLocComb, :), 'CapSize', 0, 'LineWidth', wd_errorbar, 'HandleVisibility', 'off')
                errorbar(data_ave, pred_ave, pred_sem, 'vertical', 'color', colors_comb(iLocComb, :), 'CapSize', 0, 'LineWidth', wd_errorbar, 'HandleVisibility', 'off')
                plot(data_ave, pred_ave, 'o', 'markerfacecolor', colors_comb(iLocComb, :), 'markeredgecolor', 'w', 'HandleVisibility', 'off')
                plot([0, 1], [0, 1], 'k-', 'HandleVisibility', 'off')
                legend(subjList, 'Location', 'best', 'NumColumns', 3)
                xticks(linspace(0, 1, 5)), yticks(linspace(0, 1, 5))
                xlabel(sprintf('Measured %s', namesMetrics{iMetric}))
                ylabel(sprintf('Predicted %s', namesMetrics{iMetric}))

                % Group-averaged GoF
                [R2_ave, ~, ~, R2_sem] = getCI(getCI(R2_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6), 2, 1);
                [R2_w_ave, ~, ~, R2_w_sem] = getCI(getCI(R2_w_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6), 2, 1);
                [GoF_delta_ave, ~, ~, GoF_delta_sem] = getCI(getCI(SSE_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6), 2, 1);
                [SSE_w_ave, ~, ~, SSE_w_sem] = getCI(getCI(SSE_w_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6), 2, 1);
                IC_SSE_ave = getCI(getCI(IC_SSE_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :), 1, 6), 2, 1);
                IC_SSE_w_ave = getCI(getCI(IC_SSE_w_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :), 1, 6), 2, 1);

                if ~flag_subjIsHuman
                    [R2_ave, ~, ~, R2_sem] = getCI(R2_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6);
                    [R2_w_ave, ~, ~, R2_w_sem] = getCI(R2_w_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6);
                    [GoF_delta_ave, ~, ~, GoF_delta_sem] = getCI(SSE_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6);
                    [SSE_w_ave, ~, ~, SSE_w_sem] = getCI(SSE_w_allCond(iModelA, iModelB, iLocComb, iMetric, :, :), 1, 6);
                    IC_SSE_ave = getCI(IC_SSE_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :), 1, 6);
                    IC_SSE_w_ave = getCI(IC_SSE_w_allCond(iModelA, iModelB, iLocComb, iMetric, :, :, :), 1, 6);
                end

                title(sprintf('R2: %.2f (%.2f) || %.2f (%.2f)\nSSE: %.2f (%.2f) || %.2f (%.2f)\n%s: %.0f || %.0f\nPartial corr: r = %.2f, p=%.3f\nttest: t = %.2f, p=%.3f', ...
                    R2_w_ave, R2_w_sem, R2_ave, R2_sem, ...
                    SSE_w_ave, SSE_w_sem, GoF_delta_ave, GoF_delta_sem, ...
                    namesIC{iIC_plot}, IC_SSE_w_ave(iIC_plot), IC_SSE_ave(iIC_plot), ...
                    nanmedian(r_allB), nanmedian(p_corr_allB), ...
                    nanmedian(t_allB), nanmedian(p_ttest_allB)))

            end % iMetric

            set(findall(gcf, '-property', 'fontsize'), 'fontsize', 12)
            sgtitle(sprintf('[A%dB%d] [L%d] [ni=%d] nsubj=%d', iModelA, iModelB, iLocComb, nIterations, nsubj))

            % Save the figure
            saveas(gcf, sprintf('%s/n%d_L%d_B%d_pred_vs_meas.jpg', nameFolder_Fig_NOM_AB, nsubj, iLocComb, iModelB))
            %%%%%%%%%%
        end % iLocComb

        close all
    end % iModelB

end % iModelA

%% Plot GoF

clc

% Define folder for saving GoF figures
nameFolder_Fig_NOM_GoF = sprintf('%s/GoF', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_NOM_GoF)), mkdir(nameFolder_Fig_NOM_GoF), end

% namesGoF = {'R2', 'R2w', 'SSE', 'SSEw', 'AIC', 'AICw', 'AICc', 'AICcw', 'BIC', 'BICw', 'nLL'};
namesGoF = {'R2', 'R2w', 'SSE', 'SSEw', 'AIC', 'AICw', 'AICc', 'AICcw', 'BIC', 'BICw'}; % no nLL
nGoFs = length(namesGoF);
iGoF_delta = 3:11; % do not calculate delta values for R2, which does not make sense
x_ticks = 1:nModelsA * nModelsB;

% Create x ticks labels
x_ticklabels = {};

for iModelA = iModelA_all

    for iModelB = iModelB_all
        %         x_ticklabels=[x_ticklabels, sprintf('%s-%s', namesModelA{iModelA}, namesModelB{iModelB})];
        x_ticklabels = [x_ticklabels, sprintf('%s', namesModelB{iModelB})];
    end

end

%
for iGoF = 5:nGoFs

    switch iGoF
        case 1, GoF_allCond = R2_allCond; y_lim = [-3, 1];
        case 2, GoF_allCond = R2_w_allCond; y_lim = [-3, 1];
        case 3, GoF_allCond = SSE_allCond; y_lim = [0, 2.4];
        case 4, GoF_allCond = SSE_w_allCond; y_lim = [0, 600];
        case 5, GoF_allCond = squeeze(IC_SSE_allCond(:, :, :, :, :, :, 1)); y_lim = [-6e4, 0]; y_lim = [0, 1.2e4];
        case 6, GoF_allCond = squeeze(IC_SSE_w_allCond(:, :, :, :, :, :, 1)); y_lim = [-4e4, 0]; y_lim = [0, 3e4];
        case 7, GoF_allCond = squeeze(IC_SSE_allCond(:, :, :, :, :, :, 2)); y_lim = [-6e4, 0]; y_lim = [0, 2.4e4];
        case 8, GoF_allCond = squeeze(IC_SSE_w_allCond(:, :, :, :, :, :, 2)); y_lim = [-4e4, 0]; y_lim = [0, 3e4];
        case 9, GoF_allCond = squeeze(IC_SSE_allCond(:, :, :, :, :, :, 3)); y_lim = [-6e4, 0]; y_lim = [0, 2e4];
        case 10, GoF_allCond = squeeze(IC_SSE_w_allCond(:, :, :, :, :, :, 3)); y_lim = [-4e4, 0]; y_lim = [0, 2.4e4];
        case 11, GoF_allCond = nLL_allCond; y_lim = [0, 8e3]; y_lim = [0, 4e3];
    end

    for iMetric = 1:nMetrics
        figure('Position', [0 200 length(iLocComb_all) * 333 1e3])
        isubplot = 1;

        for iLocComb = iLocComb_all

            data = GoF_allCond(:, :, iLocComb, iMetric, :, :);

            if any(iGoF == iGoF_delta)
                data_delta = nan(size(data));
                % for each subj at each iteraction, delete the min across 8 model variations (3A x 3B)
                for isubj = 1:nsubj

                    parfor iIteration = 1:nIterations
                        d = data(:, :, :, :, isubj, iIteration);
                        d_min = min(d(:));
                        data_delta(:, :, :, :, isubj, iIteration) = data(:, :, :, :, isubj, iIteration) - d_min;
                    end

                end

            end

            %             subplot(nMetrics, length(iLocComb_all), isubplot), hold on % one plot for one GoF
            subplot(3, 2, isubplot), hold on % one plot for one GoF and one metric
            GoF_delta_med_allSubj = getCI(data_delta, 1, 6);
            [GoF_delta_ave, ~, ~, GoF_delta_sem] = getCI(GoF_delta_med_allSubj, 2, 3);

            % group average
            x_all = [];

            for iModelA = iModelA_all
                x = (1:nModelsB) + (iModelA - 1) * (nModelsB + .5);
                bar(x, GoF_delta_ave(iModelA, :), 'FaceColor', 'w', 'EdgeColor', colors_comb(iLocComb, :), 'LineWidth', 2)
                errorbar(x, GoF_delta_ave(iModelA, :), GoF_delta_sem(iModelA, :), '.', 'CapSize', 0, 'Color', colors_comb(iLocComb, :), 'LineWidth', 2)
                x_all = [x_all, x];
            end % iModelA

            % idvd data
            str_all = '';

            for isubj = 1:nsubj
                y_all = [];

                for iModelA = iModelA_all
                    y = GoF_delta_med_allSubj(iModelA, :, isubj);
                    y_all = [y_all, y];
                end

                if isubj <= nMarkersMax, facecolor = 'w'; else, facecolor = ones(1, 3) / 2; end
                plot(x_all, y_all, ['-', markers_allSubj{isubj}], 'color', ones(1, 3) / 2, 'markerfacecolor', facecolor)
                iModelBest = find(y_all == min(y_all)); if length(iModelBest) > 1, fprintf('More than two models are the best\n'), end

                if any(iModelBest(1) == [2, 3, 5, 6:15])
                    str = sprintf('\n%s %s L%d %s%s: [%d] %.0f (%.0f)', namesGoF{iGoF}, namesMetrics{iMetric}, iLocComb, subjList{isubj}, markers_allSubj{isubj}, find(y_all == min(y_all)), y_all(y_all == min(y_all)), y_all(1));
                    str_all = [str_all, str];
                end

            end

            isubplot = isubplot + 1;
            xticks(x_all), xticklabels(x_ticklabels), xtickangle(45)
            y_lim_ = y_lim;

            if any(iLocComb == [1, 3, 5]), y_lim_(2) = y_lim(2) / 4;
            elseif any(iLocComb == [6, 7]), y_lim_(2) = y_lim(2) / 2;
            end

            y_ticks = linspace(y_lim_(1), y_lim_(2), 5);
            ylim(y_lim_), yticks(y_ticks)
            text(0, y_ticks(end), str_all)
            title(sprintf('L%d %s', iLocComb))
        end % iLocComb

        if any(iGoF == iGoF_delta), sgtitle(sprintf('Delta %s %s', namesGoF{iGoF}, namesMetrics{iMetric}))
        else , sgtitle(sprintf('%s %s', namesGoF{iGoF}, namesMetrics{iMetric}))
        end

        % Save the figure
        saveas(gcf, sprintf('%s/n%d_%s_%s.jpg', nameFolder_Fig_NOM_GoF, nsubj, namesGoF{iGoF}, namesMetrics{iMetric}))

    end % iMetric

    %     if any(iGoF == iGoF_delta), sgtitle(sprintf('Delta %s', namesGoF{iGoF}))
    %     else, sgtitle(sprintf('%s', namesGoF{iGoF}))
    %     end
    %     saveas(gcf, sprintf('%s/n%d_%s.jpg', nameFolder_Fig_NOM_GoF, nsubj, namesGoF{iGoF}))
end % iGoF

close all

%% Analyze GoF
clc
GoF_all = cat(7, R2_allCond, R2_w_allCond, IC_SSE_allCond, IC_SSE_w_allCond); % nModelsA x nModelsB x nLoc8 x nMetrics x nsubj x ni x nGoF
namesGoF = {'R2', 'R2w', 'AIC', 'AICw', 'AICc', 'AICcw', 'BIC', 'BICw'}; % no nLL
nGoFs = size(GoF_all, 7); assert(length(namesGoF) == nGoFs) % number of GoF measures included

% 2-way ANOVA: model A x model B; for each Loc x Metric x GoF
namesANOVA_Var = {'A', 'B', 'AxB'};
nANOVA_Var = length(namesANOVA_Var); %1=main effect of ModelA; 2=main effect of ModelB; 3=interaction
p_allB_allCond = nan(nGoFs, nLocComb8, nMetrics, nIterations, nANOVA_Var);

for iGoF = 1:nGoFs

    for iLocComb = iLocComb_all

        for iMetric = 1:nMetrics

            parfor iIteration = 1:nIterations
                GoF_all_perIte = squeeze(GoF_all(:, :, iLocComb, iMetric, :, iIteration, iGoF)); % select data of each iteraction of bootstrap

                % get index
                indModelA = nan(size(GoF_all_perIte));
                indModelB = indModelA;

                for iModelA = iModelA_all

                    for iModelB = iModelB_all
                        indModelA(iModelA, iModelB, :) = ones(nsubj, 1) * iModelA;
                        indModelB(iModelA, iModelB, :) = ones(nsubj, 1) * iModelB;
                    end

                end

                % conduct ANOVA
                [~, tbl] = print_nANOVA({'ModelA', 'ModelB'}, GoF_all_perIte(:), {indModelA(:), indModelB(:)}, nsubj, 0);
                p_allB_allCond(iLocComb, iMetric, iGoF, iIteration, :) = [tbl{2, 7}, tbl{3, 7}, tbl{4, 7}];
            end % ii

        end % iMetric

    end % iLocComb

end % iGoF

% plot
nameFolder_Fig_NOM_ANOVA = sprintf('%s/ANOVA/ANOVA_AB', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_NOM_ANOVA)), mkdir(nameFolder_Fig_NOM_ANOVA), end

for iMetric = 1:nMetrics
    figure('Position', [0 200 length(iLocComb_all) * 333 1e3])
    isubplot = 1;

    for iGoF = 1:nGoFs

        for iLocComb = iLocComb_all
            [p_med, p_lb, p_ub] = getCI(p_allB_allCond(iLocComb, iMetric, iGoF, :, :), 1, 4);

            subplot(nGoFs, length(iLocComb_all), isubplot), hold on
            bar(1:nANOVA_Var, p_med, 'edgecolor', colors_comb(iLocComb, :), 'FaceColor', 'w', 'LineWidth', 1.5)
            errorbar(1:nANOVA_Var, p_med, p_lb, p_ub, '.', 'color', colors_comb(iLocComb, :), 'CapSize', 0, 'LineWidth', 1.5)
            xticks(1:nANOVA_Var)
            xticklabels(namesANOVA_Var)
            yline(.05, 'k--');
            xlim([0, nANOVA_Var + 1])
            %             ylim([0, .1])
            ylim([0, max(p_ub)])
            text(0, .09, sprintf('%s-L%d', namesGoF{iGoF}, iLocComb))
            isubplot = isubplot + 1;
        end % iLocComb

    end % iGoF

    sgtitle(sprintf('n=%d %s', nsubj, namesMetrics{iMetric}))

    % Save the figure
    saveas(gcf, sprintf('%s/n%d_%s.jpg', nameFolder_Fig_NOM_ANOVA, nsubj, namesMetrics{iMetric}))
end % iMetric

close all

%% 2-way ANOVA: Model (with randTemp models aveaged) x Metric; for each Loc x GoF
namesANOVA_Var = {'Model', 'Metric', 'x'};
nANOVA_Var = length(namesANOVA_Var);
p_allB_allCond = nan(nGoFs, nLocComb8, nIterations, nANOVA_Var);

for iGoF = 1:nGoFs

    for iLocComb = iLocComb_all

        parfor iIteration = 1:nIterations
            GoF_all_perIte = squeeze(GoF_all(:, :, iLocComb, :, :, iIteration, iGoF)); % select data of each iteraction of bootstrap
            % nModelsA x nModelsB x (nLoc8) x nMetrics x nsubj x (ni) x (nGoF)

            nModels = nModelsB + 1;
            temp = nan(nModels, nMetrics, nsubj); % ModelA=1 (1-4), and averaged of model A=2
            temp(1:nModelsB, :, :) = GoF_all_perIte(1, :, :, :);
            temp(nModelsB + 1, :, :) = mean(GoF_all_perIte(2, :, :, :), 2);
            GoF_all_perIte = temp;

            % get index
            indModel = nan(size(GoF_all_perIte));
            indMetric = indModel;

            for iModel = 1:nModels

                for iMetric = 1:nMetrics
                    indModel(iModel, iMetric, :) = ones(nsubj, 1) * iModel;
                    indMetric(iModel, iMetric, :) = ones(nsubj, 1) * iMetric;
                end

            end

            % conduct ANOVA
            [~, tbl] = print_nANOVA({'Model', 'Metric'}, GoF_all_perIte(:), {indModel(:), indMetric(:)}, nsubj, 0);
            p_allB_allCond(iLocComb, iGoF, iIteration, :) = [tbl{2, 7}, tbl{3, 7}, tbl{4, 7}];
        end % ii

    end % iLocComb

end % iGoF

% Define folder for saving ANOVA figures
nameFolder_Fig_NOM_ANOVA = sprintf('%s/ANOVA/ANOVA_Model_Metric', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_NOM_ANOVA)), mkdir(nameFolder_Fig_NOM_ANOVA), end

figure('Position', [0 200 length(iLocComb_all) * 333 1e3])
isubplot = 1;

for iGoF = 1:nGoFs

    for iLocComb = iLocComb_all
        [p_med, p_lb, p_ub] = getCI(p_allB_allCond(iLocComb, iGoF, :, :), 1, 3);

        subplot(nGoFs, length(iLocComb_all), isubplot), hold on

        bar(1:nANOVA_Var, p_med, 'edgecolor', colors_comb(iLocComb, :), 'FaceColor', 'w', 'LineWidth', 1.5)
        errorbar(1:nANOVA_Var, p_med, p_lb, p_ub, '.', 'color', colors_comb(iLocComb, :), 'CapSize', 0, 'LineWidth', 1.5)

        xticks(1:nANOVA_Var)
        xticklabels(namesANOVA_Var)
        yline(.05, 'k--');
        xlim([0, nANOVA_Var + 1])
        ylim([0, .1])

        text(0, .09, sprintf('%s-L%d', namesGoF{iGoF}, iLocComb))
        isubplot = isubplot + 1;

    end % iLocComb

end % iGoF

sgtitle(sprintf('n=%d', nsubj))

% Save the figure
saveas(gcf, sprintf('%s/n%d.jpg', nameFolder_Fig_NOM_ANOVA, nsubj))

%% Assess corr between pA and IN parameters
% Define folder for saving figures
nameFolder_Fig_NOM_corr = sprintf('%s/corr_IN_pA', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_NOM_corr)), mkdir(nameFolder_Fig_NOM_corr), end

iModelA = 1; % models with normal templates
iMetric = 3; % pA

iLocComb_corr = [6, 5, 3];
% iLocComb_corr = [1,8];
% iLocComb_corr = 1;

for iModelB = [1, 2, 4, 5]

    switch iModelB
        case 1, iParams_all = 2; % constant noise
        case 2, iParams_all = 2; % induced noise
        case 4, iParams_all = 1; % induced noise
        case 5, iParams_all = 1; % induced noise
    end

    figure('Position', [0 100 1e3 400])

    for iParam = iParams_all
        pA_allLoc = getCI(data_allCond(iModelA, iModelB, :, iMetric, :, :, :), 2, 7); % averaged across all bins
        paramIN_allLoc = squeeze(params_allCond(iModelA, iModelB, :, :, :, iParam));
        %         paramIN_allLoc = log10(paramIN_allLoc);
        subplot(1, 2, find(iParam == iParams_all)), hold on

        for iLocComb = iLocComb_corr

            if ~flag_subjIsHuman
                pA_med = getCI(pA_allLoc(iLocComb, :, :), 1, 2);
                paramIN_med = getCI(paramIN_allLoc(iLocComb, :, :), 1, 2);
            else
                pA_med = getCI(pA_allLoc(iLocComb, :, :), 1, 3);
                paramIN_med = getCI(paramIN_allLoc(iLocComb, :, :), 1, 3);
            end

            % idvd data
            for isubj = 1:nsubj
                edgecolor = colors_comb(iLocComb, :);
                if isubj <= nMarkersMax, facecolor = 'w'; else, facecolor = ones(1, 3) / 2; end
                plot(pA_med(isubj), paramIN_med(isubj), markers_allSubj{isubj}, 'markeredgecolor', edgecolor', 'markerfacecolor', facecolor, 'LineWidth', 2, 'MarkerSize', 10)
            end

            % group ave per loc
            %             [pA_ave, ~, ~, pA_sem] = getCI(pA_med(:), 2, 1);
            %             [paramIN_ave, ~, ~, paramIN_sem] = getCI(paramIN_med(:), 2, 1);
            %             errorbar(pA_ave, paramIN_ave, pA_sem, 'horizontal', '.', 'color', edgecolor, 'CapSize', 0, 'LineWidth', 2)
            %             errorbar(pA_ave, paramIN_ave, paramIN_sem, 'vertical', '.', 'color', edgecolor, 'CapSize', 0, 'LineWidth', 2)
        end

        % group ave per subj
        %         for isubj=1:nsubj
        %             pA_med = getCI(pA_allLoc(:, isubj, :), 1, 3);
        %             paramIN_med = getCI(paramIN_allLoc(:, isubj, :), 1, 3);
        %
        %             [pA_ave, ~, ~, pA_sem] = getCI(pA_med(:), 2, 1);
        %             [paramIN_ave, ~, ~, paramIN_sem] = getCI(paramIN_med(:), 2, 1);
        %             errorbar(pA_ave, paramIN_ave, pA_sem, 'horizontal', '.', 'color', ones(1,3)/2, 'CapSize', 0, 'LineWidth', 1.5)
        %             errorbar(pA_ave, paramIN_ave, paramIN_sem, 'vertical', '.', 'color', ones(1,3)/2, 'CapSize', 0, 'LineWidth', 1.5)
        %             if isubj<=nMarkersMax, facecolor='w'; else, facecolor=ones(1,3)/2; end
        %             plot(pA_ave, paramIN_ave, markers_allSubj{isubj}, 'markeredgecolor', ones(1,3)/2', 'markerfacecolor', facecolor, 'LineWidth', 1.5)
        %         end

        r_allB = [];
        p_allB = [];

        for iIteration = 1:nIterations

            if ~~flag_subjIsHuman
                a = squeeze(pA_allLoc(iLocComb_corr, iIteration));
                b = squeeze(paramIN_allLoc(iLocComb_corr, iIteration));
            else
                a = squeeze(pA_allLoc(iLocComb_corr, :, iIteration));
                b = squeeze(paramIN_allLoc(iLocComb_corr, :, iIteration));
            end

            %             a=a-mean(a,1);
            %             b=b-mean(b,1);
            %             [r, p] = corr(a(:),b(:),'type', 'Kendall');
            [r, p] = corr(a(:), b(:));
            r_allB(iIteration) = r;
            p_allB(iIteration) = p;
        end % ii

        [r_med, r_lb, r_ub] = getCI(r_allB(:), 1, 1);
        [p_med, p_lb, p_ub] = getCI(p_allB(:), 1, 1);

        legend(subjList, 'NumColumns', 3, 'Location', 'best')
        xlabel(sprintf('%s', namesMetrics{iMetric}))
        ylabel(sprintf('[A%dB%d] %s', iModelA, iModelB, namesParamsModel_all{iModelB}{iParam}))

        title(sprintf('%s\nr=%.2f [%.2f, %.2f], p=%.3f [%.3f, %.3f]', namesParamsModel_all{iModelB}{iParam}, r_med, r_lb, r_ub, p_med, p_lb, p_ub))
    end % iParam

    sgtitle(sprintf('[A%dB%d] n=%d [Location: %s]', iModelA, iModelB, nsubj, num2str(iLocComb_corr)))

    % Save the figure
    saveas(gcf, sprintf('%s/n%d_B%d.jpg', nameFolder_Fig_NOM_corr, nsubj, iModelB))
end % iModelB

%% Compare params across locations
clc, close all

% Define folder for saving comparison figures
nameFolder_Fig_NOM_compParams = sprintf('%s/compParams', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_NOM_compParams)), mkdir(nameFolder_Fig_NOM_compParams), end

iModelA = 1;
isubplots = reshape(1:12, 4, 3)';

for iModelB = iModelB_all
    nParams = length(namesParamsModel_all{iModelB});
    figure('Position', [0 0 2e3 2e3])
    iiLocComb = 1;

    for iLoc_ = {[1, 8], [6, 7], [5, 3]}
        %     for iLoc_ = {[1,8]}

        iLoc = iLoc_{1};

        for iParam = 1:nParams
            subplot(3, 4, isubplots(iiLocComb, iParam)), hold on
            params_allSubj = squeeze(params_allCond(iModelA, iModelB, iLoc, :, :, iParam));
            params_allSubj_med = getCI(params_allSubj, 1, 3);
            [params_ave, ~, ~, params_sem] = getCI(params_allSubj_med, 2, 2);

            for iiLoc = 1:length(iLoc)
                bar(iiLoc, params_ave(iiLoc), 'FaceColor', 'w', 'EdgeColor', colors_comb(iLoc(iiLoc), :), 'BarWidth', .5)
                errorbar(iiLoc, params_ave(iiLoc), params_sem(iiLoc), 'CapSize', 0, 'color', colors_comb(iLoc(iiLoc), :))
            end

            for isubj = 1:nsubj
                if isubj <= nMarkersMax, facecolor = 'w'; else, facecolor = ones(1, 3) / 2; end
                plot([1.2, 1.8], params_allSubj_med(:, isubj), ['-', markers_allSubj{isubj}], 'color', ones(1, 3) / 2, 'markerfacecolor', facecolor, 'MarkerSize', 10)
            end

            t_allB = [];
            p_allB = t_allB;

            for iIteration = 1:nIterations
                a = params_allSubj(1, :, iIteration);
                b = params_allSubj(2, :, iIteration);
                [~, p, ~, stats] = ttest(a, b);
                t_allB(iIteration) = stats.tstat;
                p_allB(iIteration) = p;
            end

            xticks(1:length(iLoc))
            xticklabels(namesLocComb(iLoc))

            [t_med, t_lb, t_ub] = getCI(t_allB, 1, 2);
            [p_med, p_lb, p_ub] = getCI(p_allB, 1, 2);

            title(sprintf('%s\nt=%.2f [%.2f, %.2f], p=%.3f [%.3f, %.3f]', namesParamsModel_all{iModelB}{iParam}, t_med, t_lb, t_ub, p_med, p_lb, p_ub))

        end % iParam

        iiLocComb = iiLocComb + 1;
    end % iLoc_

    sgtitle(sprintf('[A%dB%d] [n=%d] [ni=%d]', iModelA, iModelB, nsubj, nIterations))
    set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 1.5)

    % Save the figure
    saveas(gcf, sprintf('%s/n%d_A%dB%d.jpg', nameFolder_Fig_NOM_compParams, nsubj, iModelA, iModelB))
end % iModelB

% close all

%% template for ANOVA
% for ii=1:ni
%     GoF_all_perIte = squeeze(GoF_all(:, :, :, :, :, ii, :)); % select data of each iteraction of bootstrap
%     % nModelsA x nModelsB x nLoc x nMetrics x nsubj x nGoF
%
%     % 5-way ANOVA: GoF x metric x model A x model B
%     % get index
%     indModelA = nan(size(GoF_all_perIte));
%     indModelB = indModelA;
%     indLocComb = indModelA;
%     indMetric = indModelA;
%     indGoF = indModelA;
%
%     for iModelA=1:nModelsA
%         for iModelB=1:nModelsB
%             for iLocComb=iLocComb_all
%                 for iMetric=1:nMetrics
%                     for iGoF=1:nGoFs
%                         indModelA(iModelA, iModelB, iLocComb, iMetric, :, iGoF) = ones(nsubj, 1) * iModelA;
%                         indModelB(iModelA, iModelB, iLocComb, iMetric, :, iGoF) = ones(nsubj, 1) * iModelB;
%                         indLocComb(iModelA, iModelB, iLocComb, iMetric, :, iGoF) = ones(nsubj, 1) * iLocComb;
%                         indMetric(iModelA, iModelB, iLocComb, iMetric, :, iGoF) = ones(nsubj, 1) * iMetric;
%                         indGoF(iModelA, iModelB, iLocComb, iMetric, :, iGoF) = ones(nsubj, 1) * iGoF;
%                     end
%                 end
%             end
%         end
%     end
%
%     GoF_all_perIte = GoF_all_perIte(:, :, iLocComb_all, :, :, :, :);
%     indModelA = indModelA(:, :, iLocComb_all, :, :, :, :);
%     indModelB = indModelB(:, :, iLocComb_all, :, :, :, :);
%     indLocComb = indLocComb(:, :, iLocComb_all, :, :, :, :);
%     indMetric = indMetric(:, :, iLocComb_all, :, :, :, :);
%     indGoF = indGoF(:, :, iLocComb_all, :, :, :, :);
%
%     % conduct ANOVA
%     str_ANOVA = print_nANOVA({'ModelA', 'ModelB', 'LocComb', 'Metric', 'GoF'}, GoF_all_perIte(:), {indModelA(:), indModelB(:), indLocComb(:), indMetric(:), indGoF(:)}, nsubj)
% end
