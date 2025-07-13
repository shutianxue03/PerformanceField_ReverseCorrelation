% Created by Shutian Xue on October 22, 2024

% This simulation aims to explain the underestimation of response consistency
% in the trial-wise noisy observer model.

% Hypothesis:
% The underestimation occurs because the template derived from reverse
% correlation (RC) is not accurate. As a result, when this template is used
% to compute the internal variable (IV) (i.e., the cross-correlation
% between the RC-derived template and the energy profile), the estimated IV
% deviates from the true IV (calculated based on the true template). Hence,
% at the most ambiguous trials (i.e., pA is at chance level, 0.5) the
% predicted pA will become a value higher than 0.5. Averaging over these
% trials leads to an expected pA that is greater than 0.5.

% Objective:
% The purpose of this simulation is to test this hypothesis. If the derived
% template deviates from the true template, we expect the predicted pA to
% be systematically underestimated. We will compute the confidence interval
% (CI) of this underestimation and see where human observers fall into this CI.

% Steps:
% 1. Define a true template and internal noise.
% 2. Simulate stimuli and responses for a subject, referred to as 'IO'.
% 3. Derive the template using RC and compare it to the true template.
% 4. Fit the trial-wise model to the simulated stimuli and data.
% 5. Compare predicted vs. measured metrics (pC, pYES, and pA), and save the
%    confidence interval of the predicted pA to evaluate how human observers
%    perform relative to these predictions.

% Setup (needs to be the same as OOD_sim)
clear all, clc, close all
addpath(genpath('Codes/fxn_exp'))  % Path to experimental functions
addpath(genpath('Codes/fxn_analysis_RC_v2'))  % Path to analysis functions
addpath(genpath('Codes/SX_toolbox'))  % Path to analysis functions

nameFolder_NOM0 = 'Data_NOM_trialWise';

% The data directory in the Carrasco Lab server
nameFolder_Data = '/Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Shutian_server/PF_RC/Data';

%=========================
% MUST be the same as shell_all_runSim
nORI = 29;
fxn_getSigma = @(SF) 3*sqrt(2*log(2))/(2*pi*SF);
noiseCST_allCond=.2;
gaborCST_allCond=.2;%.2:.1:.5;
nTrials_allCond = 5e3;%[5e3, 1e4];
iModelB_sim_allCond=[4,5];
noiseP_allCond = [0, .5];%[0, .1, .5, 1];
%--------------%
SX_RC1_setting
%--------------%
%=========================

ni = 100; % MUST be the same as OOD_sim!!
lapseRate=0;
iModelA = 1;%[1, 5];  % 1: use the derived template, 5: use the ideal template
iModelB_fit_all = iModelB_sim_allCond; % make them consistent for now; later try differ modelB for model recovery
iLocComb=1;
patchMode = 1; % 1=energy calculated from target-patches; 2=from noise patches

% Names of the performance metrics being analyzed
namesMetrics = {'pC', 'pYES', 'pA'};
nMetrics = length(namesMetrics);
nModelsA = 2; % do not use length()!! 1=legit template, 2=permuted template
nModelsB = 5;% do not use length()!!
nLocComb8 = 8;% do not use length()!!clc
nBins = 10;  % Number of bins for the analysis
iIC_plot = 3; % % Which information criterion to plot (1 = AIC, 2 = AICc, 3 = BIC)

% Function handles for calculating information criteria (IC)
getAIC_SSE = @(SSE, nParams, nData) nData * log(SSE/nData) + 2*nParams;
getAICc_SSE = @(SSE, nParams, nData) nData * log(SSE/nData) + 2*nParams + 2*nParams*(nParams+1)/(nData-nParams-1);
getBIC_SSE = @(SSE, nParams, nData) nData * log(SSE/nData) + nParams * log(nData);

% Define the noise sampling range based on pre-set values
noise.SF_low_sampling = noise.SF_low;
noise.SF_high_sampling = noise.SF_high;

% Initialize variables for the condition
cst_ln_template = 1.0;  % Constant for linear template
noise.SF_low = 1;  % Lower limit of SF sampling
noise.SF_high = 4;  % Upper limit of SF sampling
gaborSD = 0.8;  % Standard deviation of Gabor patch (default: 0.8)

% Define parameters related to stimulus
ppd = 32;  % Pixels per degree
sz_dva = 3;  % Size of stimulus in degrees of visual angle
sz_pix = sz_dva * ppd;  % Convert size to pixels
signalORI = 90;  % Orientation of signal (degrees)
gaborSF = 2;  % Spatial frequency of Gabor patches (cycles per degree)
stim.psz = sz_dva;
stim.aper_psz = sz_pix;
stim.targetOri = signalORI;
stim.phase = 0;
stim.gaborSF = gaborSF;
stim.gaborSD = gaborSD;
stim.gabor_sz = sz_dva;
stim.mask = exp_CreateCircularApertureSin(stim);

% Define noise parameters
noise.ppd = ppd;
noise.psz = sz_pix;
noise.fix_contrast = 1;
noise.ratio_gaborInTgt = 0.5;  % Ratio of Gabor signal in target
noise.ratio_base = 0.5;  % Base ratio for noise

limit0to1 = @(x) min(max(x, 0), 1);
nTrialsPerSess = 100;  % Number of trials per session

% Define folder to save results
nameFolder_Fig_NOM = sprintf('Figures/%s/ORI%dSF%d/FIGURES_IO_n%d_A%d', nameFolder_NOM0, nORI, nSF, ni, iModelA);
if isempty(dir(nameFolder_Fig_NOM)), mkdir(nameFolder_Fig_NOM), end

% Create a pool of Gabor filter
%------------------------%
[filter_sin, filter_cos] = SX_sim02_setFilters(stim, noise.filtersSF_all, filtersOri_all, fxn_getSigma, 0);
%------------------------%
fprintf('\nPool of filters created: nORI=%d, nSF=%d\n', length(filtersOri_all), length(noise.filtersSF_all))

% Define the TRUE template (ORI x SF)
stim.gaborSD = fxn_getSigma(stim.gaborSF);
template_true = exp_CreateGabor(stim, cst_ln_template);
template_true = SX_RC4_Energy_parfor(stim.mask, {template_true}, filter_sin, filter_cos);
template_true = squeeze(template_true);  %  Remove singleton dimension
new_max =.4; new_min = -.05;
template_true = (template_true * (new_max - new_min)) + new_min;

margORI_true = mean(template_true, 2);
margSF_true = mean(template_true, 1);

%% COMPILE
clc
iModelA=1; % simulating data using the legit template
nParams=2; % for iModelB=4 and 5

% Preallocate
ncriterion_allCond = 1; icriterion=1;
noiseCST = noiseCST_allCond;
kernels_allCond = nan(length(gaborCST_allCond), length(nTrials_allCond), length(noiseP_allCond), length(iModelB_sim_allCond), ncriterion_allCond, nTypes, nORI, nSF);
IV_allCond = nan(length(gaborCST_allCond), length(nTrials_allCond), length(noiseP_allCond), length(iModelB_sim_allCond), ncriterion_allCond, ni, nBins);
nTrials_perBin_allCond = IV_allCond;
data_allCond = nan(length(gaborCST_allCond), length(nTrials_allCond), length(noiseP_allCond), length(iModelB_sim_allCond), ncriterion_allCond, nMetrics, ni, nBins);
pred_allCond = data_allCond;
nLL_allCond = nan(length(gaborCST_allCond), length(nTrials_allCond), length(noiseP_allCond), length(iModelB_sim_allCond), ncriterion_allCond, ni);
params_allCond = nan(length(gaborCST_allCond), length(nTrials_allCond), length(noiseP_allCond), length(iModelB_sim_allCond), ncriterion_allCond, ni, 3); % Pre-allocate for max params
criterion_allCond = nan(length(gaborCST_allCond), length(nTrials_allCond), length(noiseP_allCond), length(iModelB_sim_allCond)); % Pre-allocate for max params
IC_nLL_allCond = nan(length(gaborCST_allCond), length(nTrials_allCond), length(noiseP_allCond), length(iModelB_sim_allCond), ncriterion_allCond, ni, 3); % AIC, AICc, BIC (3)
R2_allCond = nan(length(gaborCST_allCond), length(nTrials_allCond), length(noiseP_allCond), length(iModelB_sim_allCond), ncriterion_allCond, nMetrics, ni);
R2_w_allCond = R2_allCond;
SSE_allCond = R2_allCond;
SSE_w_allCond = SSE_allCond;
IC_SSE_allCond = nan(length(gaborCST_allCond), length(nTrials_allCond), length(noiseP_allCond), length(iModelB_sim_allCond), ncriterion_allCond, nMetrics, ni, 3); % AIC, AICc, BIC (3)
IC_SSE_w_allCond = IC_SSE_allCond;

for gaborCST = gaborCST_allCond
    igaborCST = find(gaborCST == gaborCST_allCond);
    
    for nTrials = nTrials_allCond
        inTrials = find(nTrials == nTrials_allCond);
        
        for noiseP = noiseP_allCond
            inoiseP = find(noiseP == noiseP_allCond);
            
            for iModelB_sim = iModelB_sim_allCond
                iiModelB_sim = find(iModelB_sim == iModelB_sim_allCond);
                
                nameCond_OOD = sprintf('%s/Data_OOD/ORI%dSF%d/IO/IO_nC%.0f_gC%.0f_nT%s_N%.3f_B%d', ...
                    nameFolder_Data, noiseCST*100, gaborCST*100, format_num2exp(nTrials), noiseP, iModelB_sim);
                nameCond = sprintf('IO_nC%.0f_gC%.0f_nT%s_N%.3f_B%d', noiseCST*100, gaborCST*100, format_num2exp(nTrials), noiseP, iModelB_sim);
                % LOAD KERNELS (from Data_OOD)
                nameFolder_kernel = sprintf('%s/Data_OOD/ORI%dSF%d/IO/%s/kernels*.mat', nameFolder_Data, nORI, nSF, nameCond);
                dirFolder_kernel = dir(nameFolder_kernel);
                load(sprintf('%s/%s', dirFolder_kernel.folder, dirFolder_kernel.name), 'kernels2D', 'criterion_true') % nTrials x 11
                
                str_c = sprintf('c%.1f', criterion_true);
                kernels_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, :, :, :) = kernels2D;
                clear kernels2D
                
                % LOAD NOM fitting (from Data_NOM_trialWise)
                nameFileModelIDVD_trialWise = sprintf('%s/ORI%dSF%d/IO/%s_%s/n%d_A%dB%d', nameFolder_NOM0, nORI, nSF, nameCond, str_c, ni, iModelA, iModelB_sim);
                load(nameFileModelIDVD_trialWise, 'pred_metrics_allB', 'params_est_allB', 'nLL_allB')
                
                % Calculate ICs based on nLL
                IC_nLL = [getAIC_SSE(nLL_allB, nParams, nTrials), getAICc_SSE(nLL_allB, nParams, nTrials), getBIC_SSE(nLL_allB, nParams, nTrials)];
                IC_nLL_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, :, :) = IC_nLL;
                
                for iMetric = 1:nMetrics
                    % Pre-allocate temporary arrays for this subject and model
                    IV_allB = nan(ni, nBins);  % IV for each bin
                    nTrials_allB = IV_allB;    % Trial count for each bin
                    data_allB = IV_allB;       % Observed data for each bin
                    pred_allB = IV_allB;       % Predicted data for each bin
                    
                    for ii=1:ni
                        % Extract IV per bin
                        IV_allB(ii, :) = pred_metrics_allB{ii}.metrics.IV_allBins;
                        
                        % Extract the number of trials per bin
                        nTrials_allB(ii, :) = pred_metrics_allB{ii}.metrics.nTrials_allBins;
                        
                        % Extract metrics (data and predictions) per bin
                        pred_metrics = pred_metrics_allB{ii}.metrics; % Extract once for faster access
                        switch namesMetrics{iMetric}
                            case 'pC'
                                data_allB(ii, :) = pred_metrics.pC_data_allBins;
                                pred_allB(ii, :) = pred_metrics.pC_pred_allBins;
                            case 'pA'
                                data_allB(ii, :) = pred_metrics.pA_data_allBins;
                                pred_allB(ii, :) = pred_metrics.pA_pred_allBins;
                            case 'pYES'
                                data_allB(ii, :) = pred_metrics.pYES_data_allBins;
                                pred_allB(ii, :) = pred_metrics.pYES_pred_allBins;
                        end
                        
                        % Calculate goodness-of-fit and information criteria for both weighted and unweighted cases
                        for flag_weighted=[0,1]
                            if flag_weighted, weights = nTrials_allB(ii, :);
                            else, weights = nan; % in the function, vector of ones will be created
                            end
                            % calculate SSE and R2
                            [R2, SSE] = getR2(data_allB(ii, :), pred_allB(ii, :), weights);
                            
                            % calculate information critertion based on SSE
                            nData = sum(nTrials_allB(ii, :));
                            IC_SSE = [getAIC_SSE(SSE, nParams, nData), getAICc_SSE(SSE, nParams, nData), getBIC_SSE(SSE, nParams, nData)];
                            
                            % organize
                            if flag_weighted
                                R2_w_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, iMetric, ii) = R2;
                                SSE_w_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, iMetric, ii) = SSE;
                                IC_SSE_w_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, iMetric, ii, :) = IC_SSE;
                            else
                                R2_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, iMetric, ii) = R2;
                                SSE_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, iMetric, ii) = SSE;
                                IC_SSE_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, iMetric, ii, :) = IC_SSE;
                            end % if flag_weighted
                        end % for flag_weighted=[0,1]
                        
                    end % ii
                    data_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, iMetric, :, :) = data_allB;
                    pred_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, iMetric, :, :) = pred_allB;
                end % iMetric
                
                IV_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, :, :) = IV_allB;
                nTrials_perBin_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, :, :) = nTrials_allB;
                nLL_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, :) = nLL_allB;
                params_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, :, 1:nParams) = params_est_allB;
                criterion_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim) = criterion_true;
                %                 end % criterion
            end % iModelB_sim
        end % noiseP
    end % nTrials
end % gaborCST

fprintf('\n=========== Compile finished ===========\n')

%% Template Recovery
% Compare RC-derived template to the true template
nameFolder_Fig_tempRecov = sprintf('%s/TemplateRecovery', nameFolder_Fig_NOM);
if isempty(dir(nameFolder_Fig_tempRecov)), mkdir(nameFolder_Fig_tempRecov), end

iPass = 1; %1=all trials, 2=PassA, 3=passB
iType_all = 2;
nRows = 3;%1=2D, 2=ORI marg, 3=SF marg

% plot the true template
figure, hold on % each row is all trials | pass A | pass B
imagesc(axis_tuning{1}, axis_tuning{2}, template_true')
[iORI_max, iSF_max] = find(template_true == max(template_true(:)));
plot(axis_tuning{1}(iORI_max), axis_tuning{2}(iSF_max), 'r*', 'LineWidth', 1.5)
colorbar, axis square, %caxis([-.1, 1]),
xline(0, 'r'); yline(log2(stim.gaborSF), 'r');
xticks(axisTicks_tuning{1}), xticklabels(axisTL_tuning{1}), xlim([min(axisTicks_tuning{1}), max(axisTicks_tuning{1})])
yticks(axisTicks_tuning{2}), yticklabels(axisTL_tuning{2}), ylim([min(axis_tuning{2}), max(axis_tuning{2})])
title(sprintf('True template (nORI=%d nSF=%d)', nORI, nSF))
set(findall(gcf, '-property', 'linewidth'), 'linewidth',2)
set(findall(gcf, '-property', 'fontsize'), 'fontsize',15)
saveas(gcf, sprintf('%s/true_template_%d_%d.jpg', nameFolder_Fig_tempRecov, nORI, nSF))

% for icriterion = 1:ncriterion_allCond
for iiModelB_sim = 1:length(iModelB_sim_allCond)
    for inTrials = 1:length(nTrials_allCond)
        figure('Position', [0 200 length(gaborCST_allCond)*400 1e3]) % if plot all three types
        
        isubplots = reshape(1:length(gaborCST_allCond)*nRows, length(gaborCST_allCond), nRows)';
        for igaborCST = 1:length(gaborCST_allCond)
            
            for iType = iType_all
                kernels2D = squeeze(kernels_allCond(igaborCST, inTrials, :, iiModelB_sim, icriterion, iType, :, :));
                
                % mirror the template
                kernels2D_mir = nan(size(kernels2D));
                indMir = (nORI-1)/2;
                kk_left = kernels2D(:, 1:indMir, :);
                kk_right = kernels2D(:, indMir+2:end, :);
                kk_mid = kernels2D(:, indMir+1,:);
                kk_ave = (kk_left + flip(kk_right, 2))/2;
                kernels2D_plot = cat(2, kk_ave, kk_mid, flip(kk_ave,2));
                
                % [1] plot the derived 2D kernels (averaged across noise levels)
                kernels2D_ave = getCI(kernels2D_plot, 2, 1);
                subplot(nRows, length(gaborCST_allCond), isubplots(1, igaborCST)), hold on
                [iORI_max, iSF_max] = find(kernels2D_ave == max(kernels2D_ave(:)));
                imagesc(axis_tuning{1}, axis_tuning{2}, kernels2D_ave')
                plot(axis_tuning{1}(iORI_max), axis_tuning{2}(iSF_max), 'r*', 'LineWidth', 1.5)
                colorbar, axis square, %caxis([-.1, 1]),
                xline(0, 'r', 'LineWidth', 2); yline(log2(stim.gaborSF), 'r', 'LineWidth', 2);
                xticks(axisTicks_tuning{1}), xticklabels(axisTL_tuning{1}), xlim([min(axisTicks_tuning{1}), max(axisTicks_tuning{1})])
                yticks(axisTicks_tuning{2}), yticklabels(axisTL_tuning{2}), ylim([min(axis_tuning{2}), max(axis_tuning{2})])
                [r, p] = corr(kernels2D_ave(:), template_true(:));
                title(sprintf('Gabor CST=%.0f%%\ncorr with the true temp: r=%.2f, p=%.3f', gaborCST_allCond(igaborCST)*100, r, p))
                
                % [2-3] plot marginalized kernels
                margORI = mean(kernels2D_plot, 3);
                margSF = mean(kernels2D_plot, 2);
                for iFeature=1:2
                    subplot(nRows, length(gaborCST_allCond), isubplots(iFeature+1, igaborCST)), hold on
                    switch iFeature
                        case 1, marg = margORI; marg_true = margORI_true; ref=0;y_lim = [-.15, .5];
                        case 2, marg = margSF; marg_true = margSF_true; ref=log2(stim.gaborSF); y_lim = [-.05, .2];
                    end
                    
                    str_title = cell(length(noiseP_allCond),1);
                    % plot estimated template
                    for inoiseP = 1:length(noiseP_allCond)
                        pC = getCI(getCI(data_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, :, :), 1, 6), 2,1);
                        plot(axis_tuning{iFeature}, marg(inoiseP, :), '-', 'color', ones(1,3)*(inoiseP-1)/length(noiseP_allCond), 'LineWidth', 1.5)
                        str_title{inoiseP} = sprintf('%.3f (%.0f%%)', noiseP_allCond(inoiseP), pC*100);
                    end
                    % plot the true margs
                    plot(axis_tuning{iFeature}, marg_true, 'r--', 'LineWidth', 2)
                    
                    xline(ref, 'k--'); yline(0, 'k--'); 
%                     ylim(y_lim)
                    xlabel(namesFeature{iFeature})
                    xticks(axisTicks_tuning{iFeature}), xticklabels(axisTL_tuning{iFeature}), xlim([min(axis_tuning{iFeature}), max(axis_tuning{iFeature})])
                    legend(str_title, 'Location', 'best')
                    
                end % iFeature
            end % iType
        end % igaborCST
        set(findall(gcf, '-property', 'fontsize'), 'fontsize',12)
        sgtitle(sprintf('ORI%dSF%d Signal SF=%d [%.2f, %.2f] ModelB%d, nTrials=%d', ...
            nORI, nSF, stim.gaborSF, noise.SF_low_sampling, noise.SF_high_sampling, iModelB_sim_allCond(iiModelB_sim), nTrials_allCond(inTrials)))
        
        saveas(gcf, sprintf('%s/B%d_nT%s.jpg', nameFolder_Fig_tempRecov, iModelB_sim_allCond(iiModelB_sim), format_num2exp(nTrials_allCond(inTrials))))
    end % inTrials
end % iiModelB_sim
% end % icriterion

% close all

%% Metric Recovery
% Compare measured vs. predicted metrics (as a fxn of binned IV)
nameFolder_Fig_MetricsRecov = sprintf('%s/MetricsRecovery', nameFolder_Fig_NOM);
if isempty(dir(nameFolder_Fig_MetricsRecov)), mkdir(nameFolder_Fig_MetricsRecov), end

scalingF=80; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end

for iiModelB_sim = 1:length(iModelB_sim_allCond)
    for inTrials = 1:length(nTrials_allCond)
        figure('Position', [0 200 1.5e3 1e3]) % if plot all three types
        
        isubplots = reshape(1:length(gaborCST_allCond)*nMetrics, length(gaborCST_allCond), nMetrics)';
        for igaborCST = 1:length(gaborCST_allCond)
            % extract
            [IV_allBins, ~, ~, IV_allBins_neg, IV_allBins_pos] = getCI(IV_allCond(igaborCST, inTrials, :, iiModelB_sim, icriterion, :, :), 1, 6);
            [nTrials_perBin_allBins, nTrials_perBin_allBins_lb, nTrials_perBin_allBins_ub] = getCI(nTrials_perBin_allCond(igaborCST, inTrials, :, iiModelB_sim, icriterion, :, :), 1, 6);
            [nLL_allBins, nLL_allBins_lb, nLL_allBins_ub] =           getCI(nLL_allCond(igaborCST, inTrials, :, iiModelB_sim, icriterion, :), 1, 6);
            [param_allBins, param_allBins_lb, param_allBins_ub] = getCI(params_allCond(igaborCST, inTrials, :, iiModelB_sim, icriterion, :, 1:nParams), 1, 6);
            
            for iMetric = 1:nMetrics
                % extract
                [data_allBins, ~, ~, data_allBins_neg, data_allBins_pos] = getCI(data_allCond(igaborCST, inTrials, :, iiModelB_sim, icriterion, iMetric, :, :), 1, 7);
                [pred_allBins, pred_allBins_lb, pred_allBins_ub] =            getCI(pred_allCond(igaborCST, inTrials, :, iiModelB_sim, icriterion, iMetric, :, :), 1, 7);
                
                subplot(nMetrics, length(gaborCST_allCond), isubplots(iMetric, igaborCST)), hold on
                yline(.5, 'k--', 'handlevisibility', 'off');
                ylim([0,1])
                ylabel(namesMetrics{iMetric})
                xlim([-1, 5])
                xlabel('Internal variable')
                if isubplots(iMetric, igaborCST)<=length(gaborCST_allCond), title(sprintf('Gabor CST=%.0f%%', gaborCST_allCond(igaborCST)*100)), end
                str_title = cell(length(noiseP_allCond), 1);
                for inoiseP = 1:length(noiseP_allCond)
                    str_title{inoiseP} = sprintf('%.3f', noiseP_allCond(inoiseP));
                    color = ones(1,3)*(inoiseP-1)/length(noiseP_allCond);
                    % Prediction
                    plot(IV_allBins(inoiseP, :), pred_allBins(inoiseP, :), '-', 'color', color)
                    %                         patch([IV_allBins(inoiseP, :)'; flip(IV_allBins(inoiseP, :)')], [pred_allBins_lb(inoiseP, :)'; flip(pred_allBins_ub(inoiseP, :)')], color, 'FaceAlpha', 0.3, 'LineStyle', 'none', 'HandleVisibility', 'off');
                    % Errobars of the measurement
                    %                         errorbar(IV_allBins(inoiseP, :), data_allBins(inoiseP, :), IV_allBins_neg(inoiseP, :), IV_allBins_pos(inoiseP, :), '.','horizontal', 'CapSize', 0, 'color', color, 'handlevisibility', 'off')
                    %                         errorbar(IV_allBins(inoiseP, :), data_allBins(inoiseP, :), data_allBins_neg(inoiseP, :), data_allBins_pos(inoiseP, :), '.','vertical', 'CapSize', 0, 'color', color, 'handlevisibility', 'off')
                    % Ave of the measurment
                    for iBin=1:nBins
                        facecolor='w';
                        plot(IV_allBins(inoiseP, iBin), data_allBins(inoiseP, iBin), '-o', 'markeredgecolor', color, 'markerfacecolor', facecolor, 'MarkerSize', nTrials_perBin_allBins(inoiseP, iBin)/scalingF+5, 'LineWidth', 1,  'handlevisibility', 'off')
                    end % iBin
%                     pause
                end % inoiseP
                %                     if iCond==1, legend({'Predicted metric', 'Measured metric'}, 'Location', 'best'), end
                %                     str_est = [str_est, sprintf('B%d: %s\n', iModelB_fit, num2str(round(param_allBins', 4)))];
                legend(str_title, 'Location', 'best')
            end % iMetric
        end % igaborCST
        sgtitle(sprintf('ModelB%d, nTrials=%d', iModelB_sim_allCond(iiModelB_sim), nTrials_allCond(inTrials)))
        saveas(gcf, sprintf('%s/B%d_nT%s.jpg', nameFolder_Fig_MetricsRecov, iModelB_sim_allCond(iiModelB_sim), format_num2exp(nTrials_allCond(inTrials))))
    end % inTrials
end % iiModelB_sim

close all

%% Parameter Recovery
% Compare true vs. estimated parameters
nameFolder_Fig_ParamRecov = sprintf('%s/ParamsRecovery', nameFolder_Fig_NOM);
if isempty(dir(nameFolder_Fig_ParamRecov)), mkdir(nameFolder_Fig_ParamRecov), end

for iiModelB_sim = 1:length(iModelB_sim_allCond)
    for inTrials = 1:length(nTrials_allCond)
        figure('Position', [0 0 1.5e3 400])
        nRows=2; %1=noise, 2=criterion
        
        isubplots = reshape(1:length(gaborCST_allCond)*nRows, length(gaborCST_allCond), nRows)';
        for igaborCST = 1:length(gaborCST_allCond)
            
            % 1. NOISE parameter
            iparam=1;
            subplot(nRows, length(gaborCST_allCond),isubplots(iparam, igaborCST)), hold on % noiseP
            str_legend = {};
           
            [param_med, ~, ~, p_sem] =getCI(params_allCond(igaborCST, inTrials, :, iiModelB_sim, :, :, iparam), 1, 6);
%             plot(noiseP_allCond, param_med, 'ko-')
%           
            for inoiseP=1:length(noiseP_allCond)
                plot(noiseP_allCond(inoiseP), param_med(inoiseP), ...
                    'o-','MarkerFaceColor', ones(1,3)*(inoiseP-1)/length(noiseP_allCond), 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
            end
            
            noiseP_min = min(noiseP_allCond);
            noiseP_max = max(noiseP_allCond);
            plot([noiseP_min, noiseP_max], [noiseP_min, noiseP_max], '--r'), %xlim([0,1]), ylim([0,1])
            xlabel('True param'), ylabel('Est. param')
            title(sprintf('Gabor CST=%.0f%%\n%s', gaborCST_allCond(igaborCST)*100, namesParamsModel_all{iModelB_sim_allCond(iiModelB_sim)}{iparam}))
            axis square
            
            % 2. CRITERION
            iparam=2;
            subplot(nRows, length(gaborCST_allCond), isubplots(iparam, igaborCST)), hold on % criterion
            for inoiseP=1:length(noiseP_allCond)
                [param_med, ~, ~, p_sem] =getCI(params_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, :, :, iparam), 1, 6);
                plot(criterion_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim), param_med, ...
                    'o-','MarkerFaceColor', ones(1,3)*(inoiseP-1)/length(noiseP_allCond), 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
            end
            c_min = min(criterion_allCond(:));
            c_max = max(criterion_allCond(:));
            plot([c_min, c_max], [c_min, c_max], '--r'),%xlim([4,6]), ylim([4,6])
%             legend(str_legend, 'Location', 'bestoutside')
            xlabel('True param'), ylabel('Est. param')
            title(namesParamsModel_all{iModelB_sim_allCond(iiModelB_sim)}(iparam))
            axis square
            
        end % igaborCST
        sgtitle(sprintf('ModelB%d nTrials=%d', iModelB_sim_allCond(iiModelB_sim), nTrials_allCond(inTrials)))
        saveas(gcf, sprintf('%s/B%d_nT%s.jpg', nameFolder_Fig_ParamRecov, iModelB_sim_allCond(iiModelB_sim), format_num2exp(nTrials_allCond(inTrials))))
    end % inTrials
end % iiModelB_sim
close all

%% Correlation between estimated IN and measured pA
nameFolder_Fig_corrPA = sprintf('%s/Corr_pA', nameFolder_Fig_NOM);
if isempty(dir(nameFolder_Fig_corrPA)), mkdir(nameFolder_Fig_corrPA), end

iMetric_pA = 3;
icriterion = 1;
nParams=1; %1=noise, 2=criterion

for iiModelB_sim = 1:length(iModelB_sim_allCond)
    for inTrials = 1:length(nTrials_allCond)
        figure('Position', [0 0 1.5e3 400])
        
        isubplots = reshape(1:length(gaborCST_allCond)*nParams, length(gaborCST_allCond), nParams)';
        for igaborCST = 1:length(gaborCST_allCond)
            
            % NOISE
            iparam=1;
            subplot(nParams, length(gaborCST_allCond),isubplots(iparam, igaborCST)), hold on % noiseP
            str_legend = {};
            for inoiseP=1:length(noiseP_allCond)
                [param_med, ~, ~, p_sem] =getCI(params_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, :, iparam), 1, 6);
%                 [pA_med, ~, ~, pA_sem] =getCI(getCI(data_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, iMetric_pA, :, :), 1, 7), 1, 2);
                [pA_med, ~, ~, pA_sem] =getCI(getCI(data_allCond(igaborCST, inTrials, inoiseP, iiModelB_sim, icriterion, iMetric_pA, :, :), 1, 7), 2, 1);
%                 plot(pA_med, param_med, 'o-','color', ones(1,3)*(icriterion-1)/ncriterion_allCond);
                plot(pA_med, param_med, ...
                    'o-','MarkerFaceColor', ones(1,3)*(inoiseP-1)/length(noiseP_allCond), 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
                
%                 str_legend{icriterion}=sprintf('%.1f', criterion_allCond(icriterion));
            end
            %             plot([0, 1], [0, 1], '--k'), %xlim([0,1]), ylim([0,1])
%             legend(str_legend, 'Location', 'bestoutside')
            xlabel('Sim. pA'), ylabel('Est. param')
            title(sprintf('Gabor CST=%.0f%%\n%s', gaborCST_allCond(igaborCST)*100, namesParamsModel_all{iModelB_sim_allCond(iiModelB_sim)}{iparam}))
            axis square
            
            % CRITERION
            %             iparam=2;
            %             subplot(nParams, length(gaborCST_allCond), isubplots(iparam, igaborCST)), hold on % criterion
            %             str_legend = {};
            %             for inoiseP=1:length(noiseP_allCond)
            %                 [param_med, ~, ~, p_sem] =getCI(params_allCond(igaborCST, inTrials, inoiseP, :, iiModelB_sim, :, iparam), 1, 6);
            %                 [pA_med, ~, ~, pA_sem] =getCI(getCI(data_allCond(igaborCST, inTrials, inoiseP, :, iiModelB_sim, iMetric_pA, :, :), 2, 8), 1, 2);
            %                 plot(pA_med, param_med, 'o-','color', ones(1,3)*(inoiseP-1)/length(noiseP_allCond));
            %                 str_legend{inoiseP}=sprintf('%.1f', noiseP_allCond(inoiseP));
            %             end
            %             %             plot([4, 6], [4, 6], '--k'),%xlim([4,6]), ylim([4,6])
            %             legend(str_legend, 'Location', 'bestoutside')
            %             xlabel('Sim. pA'), ylabel('Est. param')
            %             title(namesParamsModel_all{iModelB_sim}(iparam))
            %             axis square
            
        end % igaborCST
        sgtitle(sprintf('ModelB%d nTrials=%d', iModelB_sim_allCond(iiModelB_sim), nTrials_allCond(inTrials)))
        saveas(gcf, sprintf('%s/B%d_nT%s.jpg', nameFolder_Fig_corrPA, iModelB_sim_allCond(iiModelB_sim), format_num2exp(nTrials_allCond(inTrials))))
    end % inTrials
end % iiModelB_sim
close all

%% BIC
nameFolder_Fig_BIC = sprintf('%s/BIC', nameFolder_Fig_NOM);
if isempty(dir(nameFolder_Fig_BIC)), mkdir(nameFolder_Fig_BIC), end
iBIC=3;
linestyle_allB = {'-', '--'};

figure('Position', [0 0 400 1e3])
isubplot=1;
for igaborCST = 1:length(gaborCST_allCond)
    %     for icriterion = 1:ncriterion_allCond
    subplot(length(gaborCST_allCond), ncriterion_allCond, isubplot), hold on
    isubplot = isubplot+1;
    ii=1;
    str_legend = cell(length(iModelB_sim_allCond)*length(nTrials_allCond), 1);
    for iiModelB_sim = 1:length(iModelB_sim_allCond)
        for inTrials = 1:length(nTrials_allCond)
            BIC_nLL_med = getCI(IC_nLL_allCond(igaborCST, inTrials, :, iiModelB_sim, icriterion, :, iBIC), 1, 6);
            plot(noiseP_allCond, BIC_nLL_med, '-o', 'color', ones(1,3)*(inTrials-1)/nTrials_allCond(inTrials), 'LineStyle', linestyle_allB{iiModelB_sim})
            str_legend{ii} = sprintf('B%d, n=%d', iModelB_sim_allCond(iiModelB_sim), nTrials_allCond(inTrials));
            ii=ii+1;
        end % inTrials
    end % iiModelB_sim
    
    if igaborCST+icriterion==2, legend(str_legend, 'Location', 'best'), end
    %         ylim([-5e3, 2e3])
    xlabel('True noise'), ylabel('BIC')
    title(sprintf('Gabor CST=%.0f%% Criterion=%.1f', gaborCST_allCond(igaborCST)*100, criterion_allCond(icriterion)))
    
    %     end % iiModelB_sim
    
end % igaborCST
sgtitle('BIC (nLL)')
saveas(gcf, sprintf('%s/BIC.jpg', nameFolder_Fig_BIC))
close all

