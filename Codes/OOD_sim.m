%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Script name: OOD_sim.m
% Adapted by Shutian Xue on 08/27/2025
% This script simulates observer responses based on a trialwise noisy observer model (NOM), and fits the model to the simulated data.
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
%               6 = estimate multiplicative noise (Nmul) and SDadd;
%               7 = estimate multiplicative noise (Nmul) and SDadd, and a corr param reflecting shared noise across passes;
% Outputs:
%       Simulated data and energy profiles saved in the specified directories.

% Notes: The Gabor SD was normalized by SF only when creating the filters, not when creating the Gabor patches.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all, clc, close all

function OOD_sim(noiseCST, gaborCST, nTrials, noiseP_true, iModelB_sim)


% noiseCST_allCond = [0, .1, .2, .5]; % Noise contrast sensitivity thresholds
% gaborCST_allCond = [.1, .5]; % Gabor contrast sensitivity thresholds
% noiseP_allCond = [0, 0.1, 0.2]; % Proportion of noise trials
% nTrials_allCond = [1e3, 5e3, 1e4, 5e4]; % Number of trials per condition
iModelA_sim_allCond = [1]; % 1=core model, 2=randomize template, 3=use IO template
iModelB_sim_allCond = iModelB_sim; % 1=estimate lapse rate, 2=estimate additive noise (SDadd), 3=estimate multiplicative noise (Nmul)
%                             % 4=estimate additive noise (SDadd);
%                             % 5=estimate multiplicative noise (Nmul);
%                             % 6=estimate both multiplicative and additive noise;
%                             % 7=estimate both multiplicative and additive noise, and a shared noise across passes;
templateType_true = 1; % (1) raw (2) reconstructed kernel (3) mirrored template
IVType_true = 1; % 1=sum of the dot product/convolution; 2=max; 3=normalized
itype_template_true = 2; % 1=estimate template from PRS trials, ABS trials, or BOTH trials

%%
% noiseCST=.2; gaborCST=.2; nTrials=4e3; noiseP_true=.3; iModelB_sim=5;
convolveType_true = 1; % 1=dot product, 2=convolution; %used in fxn_getIV_v3
flag_PatchMode_true = 1; if flag_PatchMode_true == 1, patchMode = 'T'; else, patchMode = 'N'; end % T'=energy calculated from target-patches; 'N'=from noise patches

modelA_true = 3; % as long as modelA_true is not 2; matters for fxn_getIV_v3
nIterations = 10; % Number of iterations for simulating the data & fitting the model
flag_fminconORbads = 2; % 1=use fmincon when fitting NOM to data, faster; 2=bads, slower but better
% flag_logIV = 0;

clc, close all
time_start = datetime('now')

% Set rng seed for reproducibility
rng(1)

% Add paths for custom functions
addpath(genpath('fxn_exp'))  % Path to experimental functions
addpath(genpath('fxn_NOM'))  % Path to experimental functions
addpath(genpath('fxn_analysis_RC_v2'))  % Path to analysis functions
addpath(genpath('SX_toolbox'))  % Path to analysis functions

%%
%--------------%
SX_RC1_setting
%--------------%
nParams = length(noiseP_true);

switch iModelB_sim
    case 6
        fprintf('\nnoiseCST=%.1f, gaborCST=%.1f, nTrials=%d, Nmul=%.3f, SDadd = %.3f, iModelB_sim=%d, nORI=%d\n', ...
            noiseCST, gaborCST, nTrials, noiseP_true, iModelB_sim, nORI)
    case 7
        fprintf('\nnoiseCST=%.1f, gaborCST=%.1f, nTrials=%d, Nmul=%.3f, SDadd = %.3f, sharedN = %.3f, iModelB_sim=%d, nORI=%d\n', ...
            noiseCST, gaborCST, nTrials, noiseP_true, iModelB_sim, nORI)
    otherwise
        fprintf('\nnoiseCST=%.1f, gaborCST=%.1f, nTrials=%d, noiseP=%.3f, iModelB_sim=%d, nORI=%d\n', ...
            noiseCST, gaborCST, nTrials, noiseP_true, iModelB_sim, nORI)
end

% Define the IO's name based on stim and condition parameters
switch iModelB_sim
    case 6
        nameIO = sprintf('IO_nC%.0f_gC%.0f_nT%s_Nmul%.3f_SDadd%.3f_B%d', noiseCST*100, gaborCST*100, format_num2exp(nTrials), noiseP_true, iModelB_sim);
    case 7
        nameIO = sprintf('IO_nC%.0f_gC%.0f_nT%s_Nmul%.3f_SDadd%.3fShared%.3f_B%d', noiseCST*100, gaborCST*100, format_num2exp(nTrials), noiseP_true, iModelB_sim);
    otherwise
        nameIO = sprintf('IO_nC%.0f_gC%.0f_nT%s_N%.3f_B%d', noiseCST*100, gaborCST*100, format_num2exp(nTrials), noiseP_true, iModelB_sim);
end
% Define folder to save results in Data_OOD
nameFolder_Data_OOD_IO = sprintf('%s/%s', nameFolder_Data_OOD, nameIO); % nameFolder_Data_OOD is created in SX_RC1_setting
if isempty(dir(nameFolder_Data_OOD_IO)), mkdir(nameFolder_Data_OOD_IO), end

% Define folder to save results in Data_NOM
nameFolder_Data_NOM_IO = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, nameIO);
if isempty(dir(nameFolder_Data_NOM_IO)), mkdir(nameFolder_Data_NOM_IO), end

%% Parameters for the simulation
lapseRate_true=0;
pC_titrate = .7; % the accuracy at which threshold is measured
iModelA_fit_all = iModelA_sim_allCond;% 1=core model, 2=randomize template, 3=use IO template
iModelB_fit_all = iModelB_sim_allCond; % make them consistent for now; later try differ modelB for model recovery
iLocComb = 1;

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
fprintf('\n Created a pool of Gabor filters: nORI=%d, nSF=%d\n', length(filtersOri_all), length(noise.filtersSF_all))

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
% stim.gaborSD = fxn_getSigma_SPdomain(stim.gaborSF); % do NOT normalize by SF here
template_gabor_true = exp_CreateGabor(stim, stim.gaborCST);
template_true = SX_RC4_Energy_parfor(stim.mask, {template_gabor_true}, filter_sin, filter_cos);
template_true = squeeze(template_true);  %  Remove singleton dimension
template_true = fxn_getTemplate(template_true, templateType_true, 0);

% Normalize template (to match the value scale of derived template in OOD_xx_beforeEst)
template_true = template_true / max(template_true(:)) * 0.2; % scale the signal energy to roughly match the range of templates derived from subj data
% The same as OOD_NOM_xx_beforeTest
save(sprintf('%s/truth.mat', nameFolder_Data_OOD_IO), '*_true')
fprintf('\n Defined the true template\n')

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

%% Loop through each pair of trials
fprintf('\n Running simulation (%d pairs): \n', nPairs)
parfor iPair = 1:nPairs

    % Create the signal Gabor patch
    % gabor = exp_CreateGabor(stim, stim.gaborCST, 0);
    gabor = template_gabor_true;
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

    % Compute internal variable (cross-correlation with the true template)
    % IV_target = sum(e3D_target .* template_true, 'all');
    IV_target = fxn_getIV_v3(modelA_true, e3D_target, convolveType_true, IVType_true, template_true, [1,29]);
    assert(~isnan(IV_target), 'IV_target is NaN');
    IV_target_sim_allT(iPair) = IV_target;

    % IV_noise = sum(e3D_noise .* template_true, 'all');
    IV_noise = fxn_getIV_v3(modelA_true, e3D_noise, convolveType_true, IVType_true, template_true, [1,29]);
    assert(~isnan(IV_noise), 'IV_noise is NaN');
    IV_noise_sim_allT(iPair) = IV_noise;

    % %Visualize patches and energy profiles (optional)
    % figure('Position', [0 0 2e3 1e3])
    % subplot(2, 2, 1), imagesc(patch_target), axis square, title('Target Patch'), colormap gray
    % subplot(2, 2, 2), imagesc(patch_noise), axis square, title('Noise Patch'), colormap gray
    % subplot(2, 2, 3), imagesc(e3D_target), axis square, title('Energy Profile - Target'), colorbar, caxis([0, 6]), xline(15, 'r-'); yline(15, 'r-');
    % subplot(2, 2, 4), imagesc(e3D_noise), axis square, title('Energy Profile - Noise'), colorbar, caxis([0, 6]),xline(15, 'r-'); yline(15, 'r-');
    % sgtitle(sprintf('noise cst = %.2f, gabor cst = %.2f\nPair %d/%d: PRS=%d, Pass=%d\nPress ANY key to close this fig and proceed', ...
    %     noiseCST, gaborCST, iPair, nPairs, iPRS_allT(iPair), iPass_allT(iPair)))
    % pause, close all

    % Store data in the matrix
    if iPRS_allT(iPair)==1, dataMatrix(iPair, :) = [nan, nan, nan, nan, iLocComb, nan, nan, nan, nan, nan, stim.gaborCST];
    else, dataMatrix(iPair, :) = [nan, nan, nan, nan, iLocComb, nan, nan, nan, nan, nan, 0];
    end

    % Store energy profiles
    e3D_target_allT(iPair, :, :) = e3D_target;
    e3D_noise_allT(iPair, :, :) = e3D_noise;

    % Display progress
    if mod(iPair, nPairs/10) == 0, fprintf('    %d/%d \n', iPair, nPairs), end

end % end of iPair
fprintf('\n    All pairs simulated.\n')

%% Copy information from pass A to pass B
% Since each pair of trials is simulated twice (once in each pass), copy results from pass A to pass B

% Copy behav data
dataMatrix(nPairs+1:end, :, :) = dataMatrix(1:nPairs, :, :);
% Copy 3D energy
e3D_noise_allT(nPairs+1:end, :, :) = e3D_noise_allT(1:nPairs, :, :);
e3D_target_allT(nPairs+1:end, :, :) = e3D_target_allT(1:nPairs, :, :);
% Copy IV
IV_target_sim_allT(nPairs+1:end, :) = IV_target_sim_allT(1:nPairs, :);
IV_noise_sim_allT(nPairs+1:end, :) = IV_noise_sim_allT(1:nPairs, :);

% Assert that pass A and pass B results are identical (for validation)
assert(sum(e3D_target_allT(nPairs+1:end, :, :) - e3D_target_allT(1:nPairs, :, :), 'all') == 0)
assert(sum(e3D_noise_allT(nPairs+1:end, :, :) - e3D_noise_allT(1:nPairs, :, :), 'all') == 0)

fprintf('\n Copied trial info of pass A to B\n')

%% Normalize IV
% IV_target_sim_allT = (IV_target_sim_allT - mean(IV_target_sim_allT)) / std(IV_target_sim_allT);
% IV_noise_sim_allT = (IV_noise_sim_allT - mean(IV_noise_sim_allT)) / std(IV_noise_sim_allT);
% fprintf('\n Normalized IVs\n')nth

%% Save the energy profile for this IO (behav data saved later)
save(sprintf('%s/energy_T_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF), 'e3D_target_allT')
save(sprintf('%s/energy_N_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF), 'e3D_noise_allT')

fprintf('\n Saved simulated 3D energy (%d trials) \n', nTrials)

%% Sample internal noise and create noisy IVs
% switch patchMode
%     case 'T', IV_sim_allT = IV_target_sim_allT;
%     case 'N', IV_sim_allT = IV_noise_sim_allT;
% end

% I should be using the energy with the signal to predict response
IV_sim_allT = IV_target_sim_allT;

% Add noise to the internal variable (IV) to simulate trial-by-trial variability
switch iModelB_sim
    case 1, Nmul_true=0; SDadd_true=noiseP_true;
    case 2, Nmul_true=noiseP_true; SDadd_true=0;
    case 3, Nmul_true=0; SDadd_true=0;
    case 4, Nmul_true=0; SDadd_true=noiseP_true; lapseRate_true = 0;
    case 5, Nmul_true=noiseP_true; SDadd_true=0; lapseRate_true = 0;
    case 6, Nmul_true=noiseP_true(1); SDadd_true=noiseP_true(2); lapseRate_true = 0;
    case 7, Nmul_true=noiseP_true(1); SDadd_true=noiseP_true(2); shared = noiseP_true(3); lapseRate_true = 0;
end

% Sample and add internal noise to IV
noisyIV_sim_allT =IV_sim_allT + randn(size(IV_sim_allT)) .* sqrt(SDadd_true^2+ (IV_sim_allT.*Nmul_true).^2);

%% Derive response given noisy IVs
% Do a simple fitting to determine the criterion_true so that accuracy matches a certain level
fxn_loss_pC_ = @(criterion_potential) fxn_loss_pC(criterion_potential, pC_titrate, noisyIV_sim_allT, iPRS_allT);
criterion_true = fmincon(fxn_loss_pC_, median(noisyIV_sim_allT), [],[],[],[], min(noisyIV_sim_allT), max(noisyIV_sim_allT)) ;
% "criterion_true" is in unit of IV

% Determine the binary response for each trial (1 = "yes", 0 = "no")
resp_allT = noisyIV_sim_allT > criterion_true;

% Implement the lapse rate by randomly flipping responses
lapse_mask = rand(size(noisyIV_sim_allT)) < lapseRate_true;
resp_allT(lapse_mask) = ~resp_allT(lapse_mask);

% Compute the predicted "yes" probability for each trial based on simulated IVs
sigma_pred_allT = sqrt((IV_sim_allT*Nmul_true).^2 + SDadd_true^2);
pYES_pred_allT = 1 - normcdf(criterion_true, IV_sim_allT, sigma_pred_allT);

% Organize data in the matrix
dataMatrix(:, 1) = 1:nTrials;  % Trial index
dataMatrix(:, 2) = iSess_allT;  % Session index
dataMatrix(:, 6) = iPRS_allT;  % PRS/ABS condition
dataMatrix(:, 7) = iPass_allT;  % Pass index (A or B)
dataMatrix(:, 8) = iPair_allT;  % Pair index
dataMatrix(:, 9) = resp_allT;  % Response for each trial

% Behavioral metrics
pYES = mean(resp_allT == 1);
pHit = sum((iPRS_allT == 1) & (resp_allT == 1)) / sum(iPRS_allT == 1);
pFA = sum((iPRS_allT == 0) & (resp_allT == 1)) / sum(iPRS_allT == 0);
% pC = (pHit + 1 - pFA) / 2;
pC = mean( (iPRS_allT==1 & resp_allT==1) | (iPRS_allT==0 & resp_allT==0) );
[dprime, c_zscore] = SX_sim06_SDT(pHit, pFA);  % Sensitivity (d') and criterion (c)
respC = nan(nPairs, 1);  % Preallocate response consistency
for iPairUnik = 1:nPairs
    respAB = resp_allT(iPair_allT == iPairUnik);
    respC(iPairUnik) = respAB(1) == respAB(2);  % Response consistency between the two passes
end
pA = mean(respC);  % Average response consistency across pairs
metrics_sim = [dprime, c_zscore, pC, pHit, pFA, pA, pYES];  % Collect all behavioral metrics

% Compare the true criterion (in IV unit) and converted criterion (also in IV unit)
c_IV_all = median(IV_sim_allT) + c_zscore * sigma_pred_allT;
% figure, hold on
% histogram(c_IV_all)
% xline(median(c_IV_all), 'k--', 'DisplayName', 'Median');
% xline(mean(c_IV_all), 'k.', 'DisplayName', 'Mean');
% xline(criterion_true, 'r-', 'DisplayName', 'True criterion');
% xlabel('Converted criterion (in IV unit)')
% legend('show', 'location', 'best')

fprintf('\n\n true criterion (in IV unit)  = %.3f; converted criterion (also in IV unit) = %.3f\n\n', criterion_true, median(c_IV_all))

% Save the simulated data and metrics
% save(sprintf('%s/behavMeas_c%.1f.mat', nameFolder_Data_OOD, criterion_true), 'dataMatrix', 'metrics_sim')
save(sprintf('%s/behavMeas.mat', nameFolder_Data_OOD_IO), 'dataMatrix', 'metrics_sim')

fprintf('\n Saved simulated behav data  (%d trials) \n', nTrials)

%% Plot performance (optional)
% Define figure folder for histogram
nameFolder_Fig_IO = sprintf('Figures/');
if isempty(dir(nameFolder_Fig_IO)), mkdir(nameFolder_Fig_IO), end

figure('Position', [0 200 600 1e3])
subplot(2,1,1), hold on
histogram(noisyIV_sim_allT(iPRS_allT == 1), 'FaceColor', 'r', 'displayname', 'Signal Present')
histogram(noisyIV_sim_allT(iPRS_allT == 0), 'FaceColor', 'b', 'displayname', 'Signal Absent')
xline(criterion_true, 'linewidth', 2, 'displayname', 'True criterion');
xlim([min(noisyIV_sim_allT), max(noisyIV_sim_allT)])
% Add legend
legend('show', 'location', 'best')

subplot(2,1,2), hold on
plot(IV_sim_allT, resp_allT, 'ro', 'displayname', 'Binary response (0 or 1)')
plot(IV_sim_allT, pYES_pred_allT, 'k+', 'displayname', 'Pred pYES')
xline(criterion_true, 'linewidth', 2, 'displayname', 'True criterion');
xlabel('IV space'), ylabel('pYES')
yline(.5, 'k--')
xlim([min(noisyIV_sim_allT), max(noisyIV_sim_allT)])
metrics_sim_ = metrics_sim; metrics_sim_(3:end) = metrics_sim_(3:end)*100;
% add legend
legend('show', 'location', 'best')
if iModelB_sim==6
    sgtitle(sprintf('IV 95%% CI [%.1f, %.1f] Median = %.1f\n[TRUE] GaborCST=%.0f%%, criterion=%.1f, %s=%.2f, %s=%.2f\n[MEASURED] criterion=%.1f, pC=%.0f%%, pHit=%.0f%%, pFA=%.0f%%, pA=%.0f%%, pYES=%.0f%%', ...
        round(quantile(noisyIV_sim_allT, [.05, .95, .5]), 1), ...
        gaborCST*100, criterion_true, ...
        namesParamsModel_all{iModelB_sim}{1}, noiseP_true(1), ...
        namesParamsModel_all{iModelB_sim}{2}, noiseP_true(2), ...
        metrics_sim_(2:end)))
else
    sgtitle(sprintf('IV 95%% CI [%.1f, %.1f] Median = %.1f\n[TRUE] GaborCST=%.0f%%, criterion=%.1f, %s=%.2f\n[MEASURED] criterion=%.1f, pC=%.0f%%, pHit=%.0f%%, pFA=%.0f%%, pA=%.0f%%, pYES=%.0f%%', ...
        round(quantile(noisyIV_sim_allT, [.05, .95, .5]), 1), ...
        gaborCST*100, criterion_true, namesParamsModel_all{iModelB_sim}{1}, noiseP_true, ...
        metrics_sim_(2:end)))
end

% Save the figure
saveas(gcf, sprintf('%s/SimulatedPerfHist.jpg', nameFolder_Fig_IO))
close all

%% Standardize the energy and conduct RC
switch flag_PatchMode_true
    case 1, e3D_allT = e3D_target_allT;
    case 2, e3D_allT = e3D_noise_allT;
end
e3D_norm = normEnergy(e3D_allT, dataMatrix(:, 11), iPRS_allT);
% size of e3D_norm is nTrials x nORI x nSF (some-k x 29 x 29)
% dataMatrix(:, 11) is gabor contrast (if PRS, is a positive number; if ABS, is 0)

% Perform reverse correlation to derive kernels (templates) for signal-present and signal-absent trials
kernels2D = SX_RC6_kernel_parfor(e3D_norm, dataMatrix, filtersSF_all, filtersOri_all);
% size of kernels2D: nTypes x nORI x nSF
kernels2D = kernels2D / max(kernels2D(:)) * 0.2;  % scale the signal energy to roughly match the range of templates derived from subj data

% Plot kernels
figure('Position', [100, 100, 2e3, 2e3])
for iType=1:nTypes

    % The recovered template for each type of trials
    subplot(3,nTypes,iType), hold on
    imagesc(axis_tuning{2}, axis_tuning{1}, squeeze(kernels2D(iType, :, :)))
    axis square
    colorbar, clim([0, .2])
    xline(1, 'r-'); % log gabor SF
    yline(0, 'r-'); % Gabor ori
    xlabel('Spatial Frequency'), ylabel('Orientation')
    xticks(axisTicks_tuning{2}), xticklabels(axisTL_tuning{2})
    yticks(axisTicks_tuning{1}), yticklabels(axisTL_tuning{1})
    title(namesType{iType})

    % 2D plot of the difference between true and recovered template
    subplot(3,nTypes,iType+nTypes), hold on
    imagesc(axis_tuning{2}, axis_tuning{1}, template_true-squeeze(kernels2D(iType, :, :)))
    axis square
    colorbar, clim([-.05, .05])
    xline(1, 'r-'); % log gabor SF
    yline(0, 'r-'); % Gabor ori
    xlabel('Spatial Frequency'), ylabel('Orientation')
    xticks(axisTicks_tuning{2}), xticklabels(axisTL_tuning{2})
    yticks(axisTicks_tuning{1}), yticklabels(axisTL_tuning{1})

    % Scatter plot of the difference
    subplot(3,nTypes,iType+nTypes*2), hold on
    t = squeeze(kernels2D(iType, :, :));
    plot(t(:), template_true(:), '.')
    plot([0, .2], [0, .2], '-', 'color', [.5, .5, .5])
    xlim([0, .2]), ylim([0, .2])
    axis square
    xlabel('Recovered template'), ylabel('True template')

    % Calculate Pearson's corr
    kA = kernels2D(iType, :, :);
    kA = kA(:);
    kB = template_true(:);
    [r, p] = corr(kA, kB);

    title(sprintf('Pearson''s corr: r=%.2f, p=%.3f', r, p))

end
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 12)
set(findall(gcf, '-property', 'linewidth'), 'linewidth',2)
sgtitle(sprintf('Row 1: recovered template from all trials (before training-testing split) \n Row 2 and 3: Difference from the true template'))
saveas(gcf, sprintf('%s/RecoveredTemplate_NoSplit.jpg', nameFolder_Fig_IO))

% select the kernel for fitting
kernels2D = squeeze(kernels2D(itype_template_true, :, :));
kernels2D = fxn_getTemplate(kernels2D, templateType_true, 0);

% Save the derived kernels for future use
% save(sprintf('%s/kernel_%s.mat', nameFolder_Data_OOD_IO, patchMode), 'kernels2D', 'criterion_true')
% fprintf('\n========== Kernels saved ========== \n\n\n')

%% Fit trial-wise model to simulated data
fprintf('\n\n============= Fit the trial-wise model to simulated data =============\n')

for iModelA_fit = iModelA_fit_all

    % Simulate data and derive template via RC using the training set
    OOD_NOM_Trialwise_beforeEst({nameIO, criterion_true}, iLocComb, iModelA_fit, IVType_true, templateType_true, flag_PatchMode_true, itype_template_true, nIterations);

    for iModelB_fit = iModelB_fit_all
        % Predict response of the testing set, using the derived template
        OOD_NOM_Trialwise_Est({nameIO, criterion_true}, iLocComb, iModelA_fit, iModelB_fit, nIterations, flag_fminconORbads)
    end
end

fprintf('\n\n============= CRITERION=%.1f DONE =============\n\n', criterion_true)
% end % criterion

%% Clean up variables for the next iteration
clear *allT e2D* e3D* kernels2D* dataMatrix
% delete(sprintf('%s/energy_T_%d_%d.mat', nameFolder_NOM_IO, nORI, nSF)) % I don't know why this is needed

%% Delete the saved energy files
delete(sprintf('%s/energy_T_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF))
delete(sprintf('%s/energy_N_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF))

time_end = datetime('now')

time_end - time_start
end

