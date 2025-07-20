%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Script name: OOD_sim.m
% Script type: script
% Author: Shutian Xue
% Date created: 
% Last updated: 07/15/2025

% Description: 
%       This script simulates observer responses based on a trialwise noisy observer model (NOM), and fits the model to the simulated data.

% Example:
%       noiseCST=.2, gaborCST=.5, nTrials=5e3, noiseP=0, iModelB_sim=5
%       OOD_sim(noiseCST, gaborCST, nTrials, noiseP, iModelB_sim)

% Inputs:
%       noiseCST: contrast of noise patches (0-1)
%       gaborCST: contrast of Gabor patches (0-1)
%       nTrials: Total number of trials (>1000)
%       noiseP: Noise parameter (can be both SDadd and Nmul)
%       iModelB_sim: Model index for simulation 
%               1 = estimate lapse rate, additive noise (SDadd) and criterion;
%               2 = estimate lapse rate, multiplicative noise (Nmul) and criterion;
%               3 = estimate lapse rate and criterion;
%               4 = estimate additive noise (SDadd) and criterion;
%               5 = estimate multiplicative noise (Nmul) and criterion;

% Outputs:
%       Simulated data and energy profiles saved in the specified directories.

% Dependencies: none
% Notes: none
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function OOD_sim(noiseCST, gaborCST, nTrials, noiseP, iModelB_sim)
noiseCST_all = [0, .1, .2, .5]; % Noise contrast sensitivity thresholds
gaborCST_all = [.1, .5]; % Gabor contrast sensitivity thresholds
nTrials_all = [5000]; % Number of trials per condition
noiseP_all = [0, 0.1, 0.2]; % Proportion of noise trials
iModelA_sim_all = [1,3]; % 1=core model, 2=randomize template, 3=use IO template
iModelB_sim_all = [4];
% noiseCST=.2, gaborCST=.5, nTrials=5e3, noiseP=0, iModelB_sim=4
% for noiseCST = noiseCST_all
%     for gaborCST = gaborCST_all
%         for nTrials = nTrials_all
%             for noiseP = noiseP_all
%                 for iModelB_sim = iModelB_sim_all
%                     OOD_sim
%                 end
%             end
%         end
%     end
% end
clc, close all
time_start = datetime('now')

% Add paths for custom functions
addpath(genpath('fxn_exp'))  % Path to experimental functions
addpath(genpath('fxn_analysis_RC_v2'))  % Path to analysis functions
addpath(genpath('SX_toolbox'))  % Path to analysis functions

%%
%--------------%
SX_RC1_setting
%--------------%

fprintf('\nnoiseCST=%.1f, gaborCST=%.1f, nTrials=%d, noiseP=%.3f, iModelB_sim=%d, nORI=%d\n', ...
    noiseCST, gaborCST, nTrials, noiseP, iModelB_sim, nORI)

% Define the IO's name based on stim and condition parameters
nameIO = sprintf('IO_nC%.0f_gC%.0f_nT%s_N%.3f_B%d', noiseCST*100, gaborCST*100, format_num2exp(nTrials), noiseP, iModelB_sim);

% Define folder to save results in Data_OOD
nameFolder_Data_OOD_IO = sprintf('%s/%s', nameFolder_Data_OOD, nameIO);
if isempty(dir(nameFolder_Data_OOD_IO)), mkdir(nameFolder_Data_OOD_IO), end

% Define folder to save results in Data_NOM
nameFolder_Data_NOM_IO = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, nameIO);
if isempty(dir(nameFolder_Data_NOM_IO)), mkdir(nameFolder_Data_NOM_IO), end

%% Parameters for the simulation
lapseRate=0;
pC_titrate = .7; % the accuracy at which threshold is measured
ni = 100; %
iModelA_fit_all = iModelA_sim_all;% 1=core model, 2=randomize template, 3=use IO template
iModelB_fit_all = iModelB_sim_all; % make them consistent for now; later try differ modelB for model recovery
iLocComb = 1;
patchMode = 'N'; % 'T'=energy calculated from target-patches; 'N'=from noise patches

% Define the noise sampling range based on pre-set values
noise.SF_low_sampling = noise.SF_low;
noise.SF_high_sampling = noise.SF_high;

% Initialize variables for the condition
cst_ln_template = 1.0;  % Constant for linear template
noise.noiseCST = noiseCST;
noise.SF_low = 1;  % Lower limit of SF sampling
noise.SF_high = 4;  % Upper limit of SF sampling
gaborSD = 0.8;  % Standard deviation of Gabor patch (default: 0.8)

% Define parameters related to stimulus
ppd = 32;  % Pixels per degree
sz_dva = 3;  % Size of stimulus in degrees of visual angle
sz_pix = sz_dva * ppd;  % Convert size to pixels
signalORI = 90;  % Orientation of signal (degrees)
gaborSF = 2;  % Spatial frequency of Gabor patches (cycles per degree)
stim.gaborCST = gaborCST;
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

%% Create a pool of Gabor filters
%------------------------%
[filter_sin, filter_cos] = SX_sim02_setFilters(stim, noise.filtersSF_all, filtersOri_all, fxn_getSigma_SPdomain, 0);
%------------------------%
fprintf('\nPool of filters created: nORI=%d, nSF=%d\n', length(filtersOri_all), length(noise.filtersSF_all))

% ==== Mute below when running on HPC ====
% % Plot the energy profile of these Gabor filters
% ind_plot = 1:nSF;
% for iORI=15 % Just select one orientation, as change in ori is relfected as vertical shift in energy profiles
%     % for iORI=1:nORI
%     figure('Position', [0 0 2e3 2e3])
%     isubplot = 1;
%     for iSF=ind_plot
%         subplot(5,6, isubplot)
%
%         imagesc(filter_sin{iORI, iSF})
%
%         template_true = SX_RC4_Energy_parfor(stim.mask, filter_sin(iORI, iSF), filter_sin, filter_cos);
%         template_true = squeeze(template_true);
%
%         imagesc(axis_tuning{2}, axis_tuning{1}, template_true)
%         xline(axis_tuning{2}(iSF), 'r-', 'linewidth', 2);
%         yline(axis_tuning{1}(iORI), 'r-', 'linewidth', 2);
%         axis square
%         xticks([]), yticks([])
%         title(sprintf('ORI=%d SF=%.2f [x%.1f]', axis_tuning{1}(iORI), 2.^axis_tuning{2}(iSF), fxn_getSigma_SPdomain(axis_tuning{2}(iSF))))
%         isubplot=isubplot+1;
%     end % iSF
%     sgtitle(sprintf('Energy profile of filters at ORI=%dº\n Value in [] = Filter bandwidth, normalized by SF', axis_tuning{1}(iORI)))
% end % iORI
% ====================================

%% Define the TRUE template (ORI x SF)
stim.gaborSD = fxn_getSigma_SPdomain(stim.gaborSF);
template_true = exp_CreateGabor(stim, cst_ln_template);
template_true = SX_RC4_Energy_parfor(stim.mask, {template_true}, filter_sin, filter_cos);
template_true = squeeze(template_true);  %  Remove singleton dimension
new_max =.4; new_min = -.05;
template_true = (template_true * (new_max - new_min)) + new_min; % scale the signal energy to roughly match the range of templates derived from subbj data
save(sprintf('%s/signalEnergy', nameFolder_Data_OOD), 'template_true')

%% Simulate stimuli and responses
nPairs = nTrials / 2;  % Half the trials are pairs
nSess = nTrials / nTrialsPerSess;  % Number of sessions

% Define iPRS, iPass and iPair of each pair in two passess
iPRS_allT_OnePass = [ones(nPairs/2, 1); zeros(nPairs/2, 1)];
iPRS_allT = [iPRS_allT_OnePass; iPRS_allT_OnePass]; %1=PRS, 0=ABS
iPair_allT = [1:nPairs,1:nPairs]';
iPass_allT = [ones(nPairs, 1); ones(nPairs, 1)*2]; %1=PassA, 2=PassB
iSess_allT = repmat(1:nSess, 1, nTrialsPerSess)';

% Preallocate data structures for performance and energy
dataMatrix = nan(nTrials, 11);  % Stores trial data
e3D_target_allT = nan(nTrials, nORI, nSF);  % Nnergy for each trial
e3D_noise_allT = e3D_target_allT;
IV_target_sim_allT = nan(nTrials, 1);  % Internal variable for each trial
IV_noise_sim_allT = IV_target_sim_allT;

% Parallel loop to simulate responses for each pair of trials
mask = stim.mask;
ratio_base = noise.ratio_base;
ratio_gaborInTgt = noise.ratio_gaborInTgt;

% Loop through each pair of trials
fprintf('\nRunning simulation (%d pairs): \n', nPairs)
parfor iPair = 1:nPairs

    % Create the signal Gabor patch
    gabor = exp_CreateGabor(stim, stim.gaborCST, 0);

    % Generate filtered noise
    filtered_noise = exp_CreateFilteredNoise(noise);

    % Create stimulus patches (signal-present and absent)
    stim_PRS = limit0to1(ratio_base + gabor * ratio_gaborInTgt + filtered_noise) .* mask + ratio_base * (1 - mask);
    stim_ABS = limit0to1(ratio_base + filtered_noise) .* mask + ratio_base * (1 - mask);

    % Select the target patch based on the PRS/ABS condition
    if iPRS_allT(iPair)==1, patch_target = stim_PRS;
    else, patch_target = stim_ABS;
    end
    patch_noise = stim_ABS;  % Noise patch is always the same for both conditions

    % Calculate energy profile of the stimulus
    e3D_target = SX_RC4_Energy_parfor(mask, {patch_target}, filter_sin, filter_cos);
    e3D_target = squeeze(e3D_target);  % Remove singleton dimensions
    e3D_noise = SX_RC4_Energy_parfor(mask, {patch_noise}, filter_sin, filter_cos);
    e3D_noise = squeeze(e3D_noise);  % Remove singleton dimensions

    % Compute internal variable (cross-correlation with template)
    IV_target = sum(e3D_target .* template_true, 'all');
    IV_target_sim_allT(iPair) = IV_target;
    IV_noise = sum(e3D_noise .* template_true, 'all');
    IV_noise_sim_allT(iPair) = IV_noise;

    % Store data in the matrix
    if iPRS_allT(iPair)==1, dataMatrix(iPair, :) = [nan, nan, nan, nan, iLocComb, nan, nan, nan, nan, nan, stim.gaborCST];
    else, dataMatrix(iPair, :) = [nan, nan, nan, nan, iLocComb, nan, nan, nan, nan, nan, 0];
    end

    % Store energy profile
    e3D_target_allT(iPair, :, :) = e3D_target;
    e3D_noise_allT(iPair, :, :) = e3D_noise;

    % Display progress
    if mod(iPair, nPairs/10) == 0, fprintf('%d/%d \n', iPair, nPairs), end

end % end of iPair
fprintf('\nDONE\n')

%% Copy information from pass A to pass B
% Since each pair of trials is simulated twice (once in each pass), copy results from pass A to pass B
dataMatrix(nPairs+1:end, :, :) = dataMatrix(1:nPairs, :, :);
IV_target_sim_allT(nPairs+1:end, :) = IV_target_sim_allT(1:nPairs, :);
IV_noise_sim_allT(nPairs+1:end, :) = IV_noise_sim_allT(1:nPairs, :);
e3D_noise_allT(nPairs+1:end, :, :) = e3D_noise_allT(1:nPairs, :, :);  % Copy noise energy
e3D_target_allT(nPairs+1:end, :, :) = e3D_target_allT(1:nPairs, :, :);  % Copy target energy

% Assert that pass A and pass B results are identical (for validation)
assert(sum(e3D_target_allT(nPairs+1:end, :, :) - e3D_target_allT(1:nPairs, :, :), 'all') == 0)
assert(sum(e3D_noise_allT(nPairs+1:end, :, :) - e3D_noise_allT(1:nPairs, :, :), 'all') == 0)

%% Save the energy profile for this IO (behav data saved later)
save(sprintf('%s/energy_T_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF), 'e3D_target_allT')
save(sprintf('%s/energy_N_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF), 'e3D_noise_allT')

fprintf('\n========== Saved simulated energy  (%d trials) ========== \n', nTrials)

%% Sample internal noise and calculate response for each trial
switch patchMode
    case 'T', IV_sim_allT = IV_target_sim_allT;
    case 'N', IV_sim_allT = IV_noise_sim_allT;
end

% Add noise to the internal variable (IV) to simulate trial-by-trial variability
switch iModelB_sim
    case 1, Nmul=0; SDadd=noiseP;
    case 2, Nmul=noiseP; SDadd=0;
    case 3, Nmul=0; SDadd=0;
    case 4, Nmul=0; SDadd=noiseP; lapseRate = 0;
    case 5, Nmul=noiseP; SDadd=0; lapseRate = 0;
end

% Sample and add internal noise to IV
IV_noisy_sim_allT = IV_sim_allT + randn(size(IV_sim_allT)) * SDadd + randn(size(IV_sim_allT)).* (IV_sim_allT*Nmul);

% Do a simple fitting to determine the criterion_true so that accuracy matches a certain level
fxn_loss_pC_ = @(criterion_potential) fxn_loss_pC(criterion_potential, pC_titrate, IV_noisy_sim_allT, iPRS_allT);
criterion_true_ = fmincon(fxn_loss_pC_, median(IV_noisy_sim_allT), [],[],[],[], min(IV_noisy_sim_allT), max(IV_noisy_sim_allT)) ;
criterion_true = round(criterion_true_, 1);

% Determine the binary response for each trial (1 = "yes", 0 = "no")
resp_allT = IV_noisy_sim_allT > criterion_true;

% Implement the lapse rate by randomly flipping responses
lapse_mask = rand(size(IV_noisy_sim_allT)) < lapseRate;
resp_allT(lapse_mask) = ~resp_allT(lapse_mask);

% Compute the predicted "yes" probability for each trial based on the noisy IV
sigma_pred = sqrt((IV_sim_allT*Nmul).^2 + SDadd^2);
pYES_pred_allT = 1 - normcdf(criterion_true, IV_sim_allT, sigma_pred);

% Organize data in the matrix
dataMatrix(:, 1) = 1:nTrials;  % Trial index
dataMatrix(:, 2) = iSess_allT;  % Session index
dataMatrix(:, 6) = iPRS_allT;  % PRS/ABS condition
dataMatrix(:, 7) = iPass_allT;  % Pass index (A or B)
dataMatrix(:, 8) = iPair_allT;  % Pair index
dataMatrix(:, 9) = resp_allT;  % Response for each trial

% Behavioral metrics
pHit = sum((iPRS_allT == 1) & (resp_allT == 1)) / sum(iPRS_allT == 1);
pFA = sum((iPRS_allT == 0) & (resp_allT == 1)) / sum(iPRS_allT == 0);
pC = (pHit + 1 - pFA) / 2;
[dprime, c] = SX_sim06_SDT(pHit, pFA);  % Sensitivity (d') and criterion (c)
respC = nan(nPairs, 1);  % Preallocate response consistency
for iPairUnik = 1:nPairs
    respAB = resp_allT(iPair_allT == iPairUnik);
    respC(iPairUnik) = respAB(1) == respAB(2);  % Response consistency between the two passes
end
pA = mean(respC);  % Average response consistency across pairs
metrics_sim = [dprime, c, pC, pHit, pFA, pA];  % Collect all behavioral metrics

% Save the simulated data and metrics
% save(sprintf('%s/behavMeas_c%.1f.mat', nameFolder_Data_OOD, criterion_true), 'dataMatrix', 'metrics_sim')
save(sprintf('%s/behavMeas.mat', nameFolder_Data_OOD_IO), 'dataMatrix', 'metrics_sim')

fprintf('\n========== Saved simulated behav data  (%d trials) ========== \n', nTrials)

%% Plot performance (optional)
% Define figure folder for histogram
% nameFolder_Fig_hist = sprintf('%s/IO/Fig', nameFolder_Data_OOD);
% if isempty(dir(nameFolder_Fig_hist)), mkdir(nameFolder_Fig_hist), end

% figure('Position', [0 200 600 1e3])
% subplot(2,1,1), hold on
% histogram(IV_noisy_sim_allT(iPRS_allT == 1), 'FaceColor', 'r')
% histogram(IV_noisy_sim_allT(iPRS_allT == 0), 'FaceColor', 'b')
% xline(criterion_true, 'linewidth', 2);
% xlim([min(IV_noisy_sim_allT), max(IV_noisy_sim_allT)])

% subplot(2,1,2), hold on
% plot(IV_sim_allT, resp_allT, 'ro')
% plot(IV_sim_allT, pYES_pred_allT, 'k+')
% xline(criterion_true, 'linewidth', 2); xlabel('IV space'), ylabel('Binary Resp')
% xlim([min(IV_noisy_sim_allT), max(IV_noisy_sim_allT)])
% metrics_sim_ = metrics_sim; metrics_sim_(3:6) = metrics_sim_(3:6)*100;
% sgtitle(sprintf('IV 95%% CI [%.1f, %.1f] Median = %.1f\n[TRUE] GaborCST=%.0f%%, criterion=%.1f, %s=%.2f\n[MEASURED] criterion=%.1f, pC=%.0f%%, pHit=%.0f%%, pFA=%.0f%%, pA=%.0f%%', ...
%     round(quantile(IV_noisy_sim_allT, [.05, .95, .5]), 1), ...
%     gaborCST*100, criterion_true, namesParamsModel_all{iModelB_sim}{1}, noiseP, ...
%     metrics_sim_(2:end)))

% Save the figure
% saveas(gcf, sprintf('%s/%s.jpg', nameFolder_Fig_hist, nameCond_NOM))

% Standardize the energy to normalize across trials
switch patchMode
    case 'T', e3D_allT = e3D_target_allT;
    case 'N', e3D_allT = e3D_noise_allT;
end
e3D_norm = normEnergy(e3D_allT, dataMatrix(:, 11), iPRS_allT);

% Perform reverse correlation to derive kernels (templates) for signal-present and signal-absent trials
kernels2D = SX_RC6_kernel_parfor(e3D_norm, dataMatrix, filtersSF_all, filtersOri_all);

% Save the derived kernels for future use
% save(sprintf('%s/kernels_c%.1f.mat', nameFolder_NOM_IO, criterion_true), 'kernels2D', 'criterion_true')
save(sprintf('%s/kernel_%s.mat', nameFolder_Data_OOD_IO, patchMode), 'kernels2D', 'criterion_true')
fprintf('\n========== Kernels saved ========== \n\n\n')

%% Fit trial-wise model to simulated data

fprintf('\n\n============= Fit Fit trial-wise model to simulated data =============\n')
for iModelA_fit = iModelA_fit_all
    OOD_NOM_Trialwise_beforeEst({nameIO, criterion_true}, iLocComb, iModelA_fit, ni);
    for iModelB_fit = iModelB_fit_all
        OOD_NOM_Trialwise_Est({nameIO, criterion_true}, iLocComb, iModelA_fit, iModelB_fit, ni)
    end
end
fprintf('\n\n============= CRITERION=%.1f DONE =============\n\n', criterion_true)
% end % criterion

%% Clean up variables for the next iteration
clear *allT e2D* e3D* kernels2D* dataMatrix
% delete(sprintf('%s/energy_T_%d_%d.mat', nameFolder_NOM_IO, nORI, nSF)) % I don't know why this is needed

time_end = datetime('now')

time_end - time_start


function loss = fxn_loss_pC(criterion_potential, pC_titrate, IV_noisy_sim_allT, iPRS_allT)
resp_allT = IV_noisy_sim_allT > criterion_potential;
pHit = sum((iPRS_allT == 1) & (resp_allT == 1)) / sum(iPRS_allT == 1);
pFA = sum((iPRS_allT == 0) & (resp_allT == 1)) / sum(iPRS_allT == 0);
pC_fit = (pHit + 1 - pFA) / 2;
loss = (pC_fit-pC_titrate)^2;
end
