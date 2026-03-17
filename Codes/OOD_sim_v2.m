% function OOD_sim_v2(noiseCST, gaborCST, nTrials, noiseP_true, iModelB_sim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script name: OOD_sim.m
% Adapted by Shutian Xue on 08/27/2025
%
% This script simulates observer responses based on a trial-wise
% noisy observer model (NOM), and then runs the updated
% OOD_NOM_Trialwise_compIV / OOD_NOM_Trialwise_fitNOM pipeline
% on the simulated data.
%
% Inputs:
% noiseCST : contrast of noise patches (0–1)
% gaborCST : contrast of Gabor patches (0–1)
% nTrials : total number of trials (> 1000; must be even)
% noiseP_true : [Nmul_true, SDadd_true, rho_true]
% Nmul_true = multiplicative noise coefficient
% SDadd_true = additive noise (IV units)
% rho_true = shared-noise correlation across passes
% iModelB_sim : Model index (for labeling / recovery):
% 1 = Nmul + Nadd + Nshared (full model)
% 2 = no Nmul
% 3 = No Nadd
% 4 = No Nshared
%
% Outputs:
% - In Data_OOD/nameIO:
% behavMeas.mat: dataMatrix, metrics_sim
% energy_T_%d_%d.mat: e3D_target_allT, e3D_noise_allT,
% noise, filtersOri_all, stim, nBins
% truth.mat: template_true and related "true" values
% - In Data_NOM_Trialwise/nameIO:
% n<nIter>_A<iModelA>_compIV.mat
% n<nIter>_A<iModelA>B<iModelB>.mat
%
% Notes:
% - Gabor SD is normalized by SF only when creating the filters,
% not when creating the Gabor patches.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---------------------------- Basic setup ----------------------------- %
clc; close all;
time_start = datetime('now');
fprintf('\n===== OOD_sim START (%s) =====\n', char(time_start));

% Set RNG for reproducibility
rng(1);

% Add paths for custom functions
addpath(genpath('fxn_exp')); % experimental
addpath(genpath('fxn_NOM')); % NOM-related
addpath(genpath('fxn_analysis_RC_v2'));
addpath(genpath('SX_toolbox'));

% Global RC / NOM settings (defines nORI, nSF, nBins, folders, etc.)
%--------------%
SX_RC1_setting;
%--------------%

% Number of bootstraps for compIV / fitNOM
nIter = 10; % <-- adjust as needed
nJob = 1;
iJob = 1;

% Model A/B indices for fitting
iModelA_sim_allCond = [1]; % 1 = RC-derived template (Model A), 2=idealtemplate; 3=permuted template
iModelB_sim_allCond = iModelB_sim; % use same B for sim and fit

templateType_true = 1; % 1 = raw;
IVType_true = 1; % 1 = sum of dot product;
convolveType_true = 1; % 1 = dot product; 2 = convolution (for fxn_getIV_v3)
flag_PatchMode_true = 1; % 1 = use target energy ('T')
flag_permT = 0; %1=permute the input template per trial
IVType = 1; % 1=sum of the dot product/convolution; 2=max; 3=normalized
templateType = 1; % (1) raw (2) reconstructed kernel (3) mirrored template
itype_template = 2; % 1=estimate template from PRS trials, ABS trials, or BOTH trials
flag_PatchMode = 1; % if flag_PatchMode == 1, patchMode = 'T'; else, patchMode = 'N'; end
flag_plot_template = 0;
nRep = 20;
iFamily_ORI = 1; % 1=scaled gaussian, 8=DoG
iFamily_SF = 2; % 2=log parabola
problem_setting = MultiStart('StartPointsToRun', 'bounds','UseParallel', 1, 'Display', 'off');
flag_plot_tuning = 0;

% Internal-noise true parameters (from input)
Nmul_true = noiseP_true(1);
Nadd_true = noiseP_true(2);
Nshared_true = noiseP_true(3);

pC_titrate = 0.7; % criteron titration target accuracy
iLocComb = 1; % use single-location index (e.g., fovea) for IO

lapseRate_true = 0; % no lapses in this sim

% ---------------------------- Print info ------------------------------ %
nParams = length(noiseP_true);

fprintf('\nSim settings:\n');
fprintf(' - noiseCST = %.2f\n', noiseCST);
fprintf(' - gaborCST = %.2f\n', gaborCST);
fprintf(' - nTrials = %d\n', nTrials);
fprintf(' - Nmul_true = %.1f, Nadd_true = %.0f, Nshared_true = %.1f (nParams=%d)\n', Nmul_true, Nadd_true, Nshared_true, nParams);
fprintf(' - iModelB_sim = %d (see updated Model B definitions)\n', iModelB_sim);
fprintf(' - nIter (for compIV/fitNOM) = %d\n\n', nIter);

% ---------------------- Define IO name & folders ---------------------- %

nameIO = sprintf('IO_cN%.0f_cG%.0f_nT%s_Nm%.2f_Na%.0f_Ns%.1f_%d%d_B%d', ...
    noiseCST*100, gaborCST*100, format_num2exp(nTrials), Nmul_true, Nadd_true, Nshared_true, nORI, nSF, iModelB_sim);

% Folder to save IO data (energy + behav)
nameFolder_Data_OOD_IO = sprintf('%s/%s', nameFolder_Data_OOD, nameIO);
if isempty(dir(nameFolder_Data_OOD_IO)), mkdir(nameFolder_Data_OOD_IO); end

% Folder to save NOM trial-wise fits
nameFolder_Data_NOM_IO = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, nameIO);
if isempty(dir(nameFolder_Data_NOM_IO)), mkdir(nameFolder_Data_NOM_IO); end

% Folder to save figures
nameFolder_Figures_perSubj = sprintf('%s/IO/%s', nameFolder_Figures, nameIO);

if isempty(dir(nameFolder_Figures_perSubj))
    mkdir(nameFolder_Figures_perSubj);
end

fprintf('\nIO name: %s\n', nameIO);
fprintf('Data_OOD folder: %s\n', nameFolder_Data_OOD_IO);
fprintf('Data_NOM folder: %s\n\n', nameFolder_Data_NOM_IO);

%% Stimulus / noise parameters --------------------- %
% Some values are taken from SX_RC1_setting (nORI, nSF, noise.filtersSF_all, etc.)

% Define noise structure
noise.noiseCST = noiseCST;
noise.SF_low = 1;
noise.SF_high = 4;
noise.SF_low_sampling = noise.SF_low;
noise.SF_high_sampling= noise.SF_high;
noise.ppd = 32; % pixels per deg
noise.psz = 3 * noise.ppd;
noise.fix_contrast = 1;
noise.ratio_gaborInTgt= 0.5;
noise.ratio_base = 0.5;

% Stim parameters
ppd = noise.ppd;
sz_dva = 3;
sz_pix = sz_dva * ppd;
signalORI= 90;
gaborSF = 2;
gaborSD = 0.8;

stim.gaborCST = gaborCST;
stim.psz = sz_dva;
stim.aper_psz = sz_pix;
stim.targetOri= signalORI;
stim.phase = 0;
stim.gaborSF = gaborSF;
stim.gaborSD = gaborSD;
stim.gabor_sz = sz_dva;
stim.mask = exp_CreateCircularApertureSin(stim);

limit0to1 = @(x) min(max(x, 0), 1);
nTrialsPerSess = 100;
nPairs = nTrials / 2;
nSess = nTrials / nTrialsPerSess;

if mod(nTrials, 2) ~= 0
    error('nTrials must be even (because nPairs = nTrials/2).');
end

%% Create Gabor filters used for energy computation --------------------- %
[filter_sin, filter_cos] = SX_sim02_setFilters(stim, noise.filtersSF_all, filtersOri_all, fxn_getSigma_SPdomain, 0);
fprintf('\nCreated filter pool: nORI=%d, nSF=%d\n', length(filtersOri_all), length(noise.filtersSF_all));

%% True 2D template (ORI×SF) --------------------- %
template_gabor_true = exp_CreateGabor(stim, stim.gaborCST);
%---------------------------%
template_true = SX_RC4_Energy_parfor(stim.mask, {template_gabor_true}, filter_sin, filter_cos);
%---------------------------%
template_true = squeeze(template_true);
%---------------------------%
template_true = fxn_getTemplate(template_true, templateType_true, 0);
%---------------------------%

% Normalize template to roughly match scale of subject-derived templates
template_true = template_true / max(template_true(:)) * 0.2;

% Save "truth" (for IO template in Model A = 2)
save(sprintf('%s/truth.mat', nameFolder_Data_OOD_IO), '*_true');
fprintf('\nDefined and saved the true template\n');

%% Preallocate sim arrays -------------------------- %
nMetrics = 11;
iPRS_allT_OnePass = [ones(nPairs/2, 1); zeros(nPairs/2, 1)]; % 1=PRS, 0=ABS
iPRS_allT = [iPRS_allT_OnePass; iPRS_allT_OnePass];
iPair_allT = [1:nPairs, 1:nPairs]';
iPass_allT = [ones(nPairs,1); ones(nPairs,1)*2];
iSess_allT = repmat(1:nSess, 1, nTrialsPerSess)';

dataMatrix = nan(nTrials, nMetrics);
e3D_target_allT = nan(nTrials, nORI, nSF);
% e3D_noise_allT = nan(nTrials, nORI, nSF);
IV_target_sim_allT= nan(nTrials, 1);
% IV_noise_sim_allT = nan(nTrials, 1);

mask = stim.mask;
ratio_base = noise.ratio_base;
ratio_gaborInTgt = noise.ratio_gaborInTgt;

%% ---------------------- Simulate stimuli & IV ------------------------- %
fprintf('\nRunning simulation (%d pairs):\n', nPairs);

for iPair = 1:nPairs

    % Filtered noise
    %----------------------------------%
    filtered_noise = exp_CreateFilteredNoise(noise);
    %----------------------------------%

    % Signal-present and signal-absent patches
    if iPRS_allT(iPair) == 1
        % Use the true Gabor as "signal"
        gabor = template_gabor_true;
        patch_target = limit0to1(ratio_base + gabor * ratio_gaborInTgt + filtered_noise) .* mask + ratio_base * (1 - mask);
    else
        patch_target = limit0to1(ratio_base + filtered_noise) .* mask + ratio_base * (1 - mask);
    end
    % patch_noise = stim_ABS;

    % Energy profiles
    %---------------------------%
    e3D_target = SX_RC4_Energy_parfor(mask, {patch_target}, filter_sin, filter_cos);
    %---------------------------%
    e3D_target = squeeze(e3D_target);
    % e3D_noise = SX_RC4_Energy_parfor(mask, {patch_noise}, filter_sin, filter_cos);
    % e3D_noise = squeeze(e3D_noise);

    % Internal variable from energy × template_true
    %---------------------------%
    IV_target = fxn_getIV_v3(e3D_target, template_true, convolveType_true, IVType_true, flag_permT, [1, nORI]);
    %---------------------------%
    assert(~isnan(IV_target), 'IV_target is NaN');
    IV_target_sim_allT(iPair) = IV_target;

    % Data matrix: col 11 = gabor contrast
    if iPRS_allT(iPair) == 1
        dataMatrix(iPair, :) = [nan, nan, nan, nan, iLocComb, nan, nan, nan, nan, nan, stim.gaborCST];
    else
        dataMatrix(iPair, :) = [nan, nan, nan, nan, iLocComb, nan, nan, nan, nan, nan, 0];
    end

    % Store energy
    e3D_target_allT(iPair, :, :) = e3D_target;
    % e3D_noise_allT(iPair, :, :) = e3D_noise;

end

fprintf('\nAll pairs (pass A) simulated.\n');

% Copy pass A -> pass B
dataMatrix(nPairs+1:end, :, :) = dataMatrix(1:nPairs, :, :);
e3D_target_allT(nPairs+1:end,:,:) = e3D_target_allT(1:nPairs, :, :);
% e3D_noise_allT(nPairs+1:end,:,:) = e3D_noise_allT(1:nPairs, :, :);
IV_target_sim_allT(nPairs+1:end,:) = IV_target_sim_allT(1:nPairs, :);
% IV_noise_sim_allT(nPairs+1:end,:) = IV_noise_sim_allT(1:nPairs, :);

assert(sum(e3D_target_allT(nPairs+1:end,:,:) - e3D_target_allT(1:nPairs,:,:), 'all') == 0);
% assert(sum(e3D_noise_allT(nPairs+1:end,:,:) - e3D_noise_allT(1:nPairs,:,:), 'all') == 0);

fprintf('\nCopied trial info from pass A to pass B.\n');

%% --------------------- Save energy for compIV ------------------------- %
% Use patchMode = 'T' in compIV; this file must contain both target + noise
% energy as well as noise, filters, stim, nBins.
save(sprintf('%s/energy_T_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF), 'e3D_target_allT', 'noise', 'filtersOri_all', 'stim', 'nBins');

fprintf('\nSaved simulated 3D energy (%d trials).\n', nTrials);

%% ---------------------- Internal noise & responses -------------------- %
% For now, use IV from target-only energy to drive responses
IV_sim_allT = IV_target_sim_allT;

% Sample internal noise (full model: Nadd + Nmul*IV)
% sigma_pred_allT = sqrt((IV_sim_allT .* Nmul_true).^2 + Nadd_true^2 + Nshared_true^2);

% Independent (pass-specific) SD per trial
sigma_priv_allT = sqrt((IV_sim_allT .* Nmul_true).^2 + Nadd_true^2);

% Draw independent z per trial (per pass)
z_priv_allT = randn(size(IV_sim_allT));

% Shared SD per trial (often constant)
sigma_shared_allT  = Nshared_true * ones(size(IV_sim_allT));
pairIDs = unique(iPair_allT); % Draw one shared z per PAIR, reuse for both passes
z_sh_perPair = randn(numel(pairIDs), 1);
[~, idxPair] = ismember(iPair_allT, pairIDs); % Map pair -> z_sh for each trial
z_sh_allT = z_sh_perPair(idxPair);

% Construct noisy DV / noisyIV
noisyIV_sim_allT = IV_sim_allT + sigma_shared_allT .* z_sh_allT + sigma_priv_allT .* z_priv_allT;

% Total predicted SD if you need it for anything
sigma_pred_allT = sqrt(sigma_priv_allT.^2 + sigma_shared_allT.^2);

% noisyIV_sim_allT = IV_sim_allT + randn(size(IV_sim_allT)) .* sigma_pred_allT;

% Titrate a criterion to reach target accuracy pC_titrate
%---------------------------%
fxn_loss_pC_ = @(criterion_potential) fxn_loss_pC(criterion_potential, pC_titrate, noisyIV_sim_allT, iPRS_allT);
%---------------------------%
criterion_true = fmincon(fxn_loss_pC_, median(noisyIV_sim_allT), [], [], [], [], min(noisyIV_sim_allT), max(noisyIV_sim_allT));

% Binary responses
resp_allT = noisyIV_sim_allT > criterion_true;

% Lapses (if any)
lapse_mask = rand(size(noisyIV_sim_allT)) < lapseRate_true;
resp_allT(lapse_mask) = ~resp_allT(lapse_mask);

% pYES per trial from SDT mapping
pYES_pred_allT = 1 - normcdf(criterion_true, IV_sim_allT, sigma_pred_allT);

% Fill behavior columns in dataMatrix
dataMatrix(:, 1) = (1:nTrials)'; % trial index
dataMatrix(:, 2) = iSess_allT; % session index
dataMatrix(:, 5) = iLocComb; % location index
dataMatrix(:, 6) = iPRS_allT; % PRS/ABS
dataMatrix(:, 7) = iPass_allT; % pass index
dataMatrix(:, 8) = iPair_allT; % pair index
dataMatrix(:, 9) = resp_allT; % response

%% ---------------------- Behavioral metrics ---------------------------- %
pYES_sim = mean(resp_allT == 1);
pHit_sim = sum((iPRS_allT == 1) & (resp_allT == 1)) / sum(iPRS_allT == 1);
pFA_sim = sum((iPRS_allT == 0) & (resp_allT == 1)) / sum(iPRS_allT == 0);
pC_sim = mean( (iPRS_allT==1 & resp_allT==1) | (iPRS_allT==0 & resp_allT==0) );

%--------------------------------------%
[dprime_sim, c_zscore] = SX_sim06_SDT(pHit_sim, pFA_sim);
%--------------------------------------%

respC_sim = nan(nPairs, 1);
for iPairUnik = 1:nPairs
    respAB = resp_allT(iPair_allT == iPairUnik);
    respC_sim(iPairUnik) = (respAB(1) == respAB(2));
end
pA_sim = mean(respC_sim);

metrics_sim = [dprime_sim, c_zscore, pC_sim, pHit_sim, pFA_sim, pA_sim, pYES_sim];

% Compare true criterion (IV space) vs SDT-derived criterion (converted)
c_IV_all = median(IV_sim_allT) + c_zscore .* sigma_pred_allT;
fprintf('\nTrue criterion (IV units) = %.3f\n', criterion_true);
fprintf('Converted criterion (IV) = %.3f (median over trials)\n\n', median(c_IV_all));

% Save behavioral measures in behavMeas.mat (as expected by compIV)
save(sprintf('%s/behavMeas.mat', nameFolder_Data_OOD_IO), 'dataMatrix', 'metrics_sim');
fprintf('\nSaved simulated behavioral data (%d trials).\n', nTrials);

%% -------------------------- Plotting DV dist and behav metrics------------------------- %
figure('Position', [0 200 600 1e3]);
subplot(2,1,1); hold on;
histogram(noisyIV_sim_allT(iPRS_allT == 1), 'FaceColor', 'r', 'DisplayName', 'Signal Present', 'normalization', 'probability');
histogram(noisyIV_sim_allT(iPRS_allT == 0), 'FaceColor', 'b', 'DisplayName', 'Signal Absent', 'normalization', 'probability');
xline(criterion_true, 'LineWidth', 2, 'DisplayName', 'True criterion');
xlim([min(noisyIV_sim_allT), max(noisyIV_sim_allT)]);
ylabel('Proportion');
legend('show', 'location', 'best');

subplot(2,1,2); hold on;
plot(IV_sim_allT, resp_allT, 'ro', 'DisplayName', 'Binary response');
plot(IV_sim_allT, pYES_pred_allT, 'k+', 'DisplayName', 'Pred pYES');
xline(criterion_true, 'LineWidth', 2, 'DisplayName', 'True criterion');
xlabel('IV'); ylabel('pYES');
yline(0.5, 'k--');
xlim([min(noisyIV_sim_allT), max(noisyIV_sim_allT)]);
metrics_sim_ = metrics_sim; metrics_sim_(3:end) = metrics_sim_(3:end)*100;

sgtitle(sprintf('CI_{95}=[%.1f, %.1f], Median=%.1f\n[TRUE] gN=%.0f%%, gC=%.0f%%, crit=%.1f, Nmul=%.1f, Nadd=%.0f, Nshared=%.0f\n[MEASURED] c_z=%.1f, pC=%.0f%%, pHit=%.0f%%, pFA=%.0f%%\npA=%.0f%%, pYES=%.0f%%', ...
    round(quantile(noisyIV_sim_allT, [.05, .95, .5]), 1), noiseCST*100, gaborCST*100, criterion_true, Nmul_true, Nadd_true, Nshared_true, metrics_sim_(2:end)));

saveas(gcf, sprintf('%s/0PerfHist.jpg', nameFolder_Figures_perSubj));
close all;

%% ---------------------- Criterion sweep + plotting ---------------------- %
nCriterion = 5;

% Criteria evenly spaced between 10th and 90th percentile of noisy IV
criterion_bounds = quantile(noisyIV_sim_allT, [0.2, 0.8]);
criterion_all = linspace(criterion_bounds(1), criterion_bounds(2), nCriterion);

% Preallocate metrics: [criterion, pYES, pHit, pFA, pC, dprime, c_zscore, pA]
metrics_byCriterion = nan(nCriterion, 8);

str_title = '';
figure('Position', [100 50 800 300*nCriterion]);

for iC = 1:nCriterion
    criterion = criterion_all(iC);

    % Binary responses
    resp_allT = noisyIV_sim_allT > criterion;

    % pYES per trial from SDT mapping
    pYES_pred_allT = 1 - normcdf(criterion, IV_sim_allT, sigma_pred_allT);

    % ---------------------- Behavioral metrics ------------------------- %
    pYES_sim = mean(resp_allT == 1);
    pHit_sim = sum((iPRS_allT == 1) & (resp_allT == 1)) / sum(iPRS_allT == 1);
    pFA_sim  = sum((iPRS_allT == 0) & (resp_allT == 1)) / sum(iPRS_allT == 0);
    pC_sim   = mean((iPRS_allT == 1 & resp_allT == 1) | (iPRS_allT == 0 & resp_allT == 0));

    [dprime_sim, c_zscore] = SX_sim06_SDT(pHit_sim, pFA_sim);

    respC_sim = nan(nPairs, 1);
    for iPairUnik = 1:nPairs
        respAB = resp_allT(iPair_allT == iPairUnik);
        respC_sim(iPairUnik) = (respAB(1) == respAB(2));
    end
    pA_sim = mean(respC_sim);

    % Store metrics
    metrics_byCriterion(iC, :) = [criterion, pYES_sim, pHit_sim, pFA_sim, pC_sim, dprime_sim, c_zscore, pA_sim];

    % -------------------------- Plotting ------------------------------- %
    % Top: histogram
    subplot(nCriterion, 2, 2*iC - 1); hold on;
    histogram(noisyIV_sim_allT(iPRS_allT == 1), ...
        'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', 'Signal Present', 'normalization', 'probability');
    histogram(noisyIV_sim_allT(iPRS_allT == 0), ...
        'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', 'Signal Absent', 'normalization', 'probability');
    xline(criterion, 'k-', 'LineWidth', 2, 'DisplayName', 'Criterion');
    xlim([min(noisyIV_sim_allT), max(noisyIV_sim_allT)]);
    ylabel('Proportion');
    % title(sprintf('Criterion = %.2f | pC = %.1f%% | pHit = %.1f%% | pFA = %.1f%%', ...
    %     criterion, 100*pC_sim, 100*pHit_sim, 100*pFA_sim));

    if iC == 1
        legend('show', 'Location', 'best');
    end

    % Bottom/right: scatter + predicted pYES
    subplot(nCriterion, 2, 2*iC); hold on;
    plot(IV_sim_allT, resp_allT, 'ro', 'MarkerSize', 4, 'DisplayName', 'Binary response');
    plot(IV_sim_allT, pYES_pred_allT, 'k+', 'DisplayName', 'Pred pYES');
    xline(criterion, 'k-', 'LineWidth', 2, 'DisplayName', 'Criterion');
    yline(0.5, 'k--');
    xlim([min(noisyIV_sim_allT), max(noisyIV_sim_allT)]);
    ylim([-0.05, 1.05]);
    xlabel('IV');
    ylabel('pYES / Response');

    if iC == 1
        legend('show', 'Location', 'best');
    end
    str_title = sprintf('%s\ncriterion=%.2f, pC=%.0f%%, dprime=%.1f, pHit=%.0f%%, pFA=%.0f%%, criterion(z)=%.1f', ...
        str_title, criterion, 100*pC_sim, dprime_sim, 100*pHit_sim, 100*pFA_sim, c_zscore);
end

sgtitle(sprintf('Behavior across %d criterion values\n%s', nCriterion, str_title));

saveas(gcf, sprintf('%s/0PerfHist_%dcriterion.jpg', nameFolder_Figures_perSubj, nCriterion));
close all

%% ---------------- Run compIV and fitNOM on this IO --------------------- %
for iModelA_fit = iModelA_sim_allCond
    % Step 1: compute IVs, templates, and test-set metrics
    %------------------------------------%
    OOD_NOM_Trialwise_compIV({nameIO, criterion_true}, iLocComb, iModelA_fit, nIter, nJob, iJob)
    %------------------------------------%

    for iModelB_fit = iModelB_sim_allCond
        % Step 2: fit NOM parameters and predict metrics
        %------------------------------------%
        OOD_NOM_Trialwise_fitNOM({nameIO, criterion_true}, iLocComb, iModelA_fit, iModelB_fit, nIter, nJob, iJob)
        %------------------------------------%
    end
end

%% --------------------- Clean-up temporary files ------------------------ %
clear *allT e2D* e3D* dataMatrix;

% Remove the energy file to save space (compIV has already loaded it)
delete(sprintf('%s/energy_T_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF));

time_end = datetime('now');
fprintf('\n===== OOD_sim END (%s), elapsed = %s =====\n\n', char(time_end), char(time_end - time_start));

% end
