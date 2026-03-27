function OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, Cz_true, lambda_whiten, iModelA_sim, iModelB_sim)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script name: OOD_sim.m
% Adapted by Shutian Xue on 03/27/2026
%
% This script simulates observer responses based on a trial-wise
% noisy observer model (NOM), and then runs the updated
% OOD_NOM_Trialwise_compIV / OOD_NOM_Trialwise_fitNOM pipeline
% on the simulated data.
%
% Notes:
% - Gabor SD is normalized by SF only when creating the filters,
% not when creating the Gabor patches.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---------------------------- Basic setup ----------------------------- %
clc; close all;
warning off; % (You may want to remove this once things are stable.)
format compact;

time_start = datetime('now');
fprintf('\n\n%s: Simulation starts\n\n', datetime('now'))

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
nORI = nORI; % nORI is defined in SX_RC1_setting; repeated here just for parfor loop to work
nSF = nORI;

% Number of bootstraps for compIV / fitNOM
nIter = 100; % <-- adjust as needed
nJob = 1;
iJob = 1;

% Model A/B indices for fitting
iModelA_sim_allCond = iModelA_sim; % 1 = RC-derived template (Model A), 2=idealtemplate; 3=permuted template
iModelB_sim_allCond = iModelB_sim; % use same B for sim and fit

templateType_true = 1; % 1 = raw;
IVType_true = 1; % 1 = sum of dot product;
convolveType_true = 1; % 1 = dot product; 2 = convolution (for fxn_getIV_v3)
flag_permT = 0; %1=permute the input template per trial
flag_plotDist = 0;

% Internal-noise true parameters (from input)
switch iModelB_sim
    case 1
        % Nmul_true = noiseP_true(1);
        % Nadd_true = noiseP_true(2);
        % Nshared_true = noiseP_true(3);
        % Cz_true = noiseP_true(4);
        str_IN = sprintf('Nmul_true = %.1f, Nadd_true = %.0f, Nshared_true = %.1f,', Nmul_true, Nadd_true, Nshared_true);
    case 2 % no Nshared
        % Nmul_true = noiseP_true(1);
        % Nadd_true = noiseP_true(2);
        Nshared_true = 0;
        % Cz_true = noiseP_true(4);
        str_IN = sprintf('Nmul_true = %.1f, Nadd_true = %.0f,', Nmul_true, Nadd_true);
    case 3 % no Nmul
        Nmul_true = 0;
        % Nadd_true = noiseP_true(2);
        % Nshared_true = noiseP_true(3);
        % Cz_true = noiseP_true(4);
        str_IN = sprintf('Nadd_true = %.0f, Nshared_true = %.1f,', Nadd_true, Nshared_true);
    case 4 % no Nadd
        % Nmul_true = noiseP_true(1);
        Nadd_true = 0;
        % Nshared_true = noiseP_true(3);
        % Cz_true = noiseP_true(4);
        str_IN = sprintf('Nmul_true = %.1f, Nshared_true = %.1f,', Nmul_true, Nshared_true);
    otherwise
        error('ALERT: no matching ModelB for simulation!!')
end

iLocComb = 1; % use single-location index (e.g., fovea) for IO

% ---------------------------- Print info ------------------------------ %
nParams = numel(namesModelBparams{iModelB_sim});

fprintf(' - nIter (for compIV/fitNOM) = %d\n', nIter);
fprintf(' - noiseCST = %.2f\n', noiseCST);
fprintf(' - gaborCST = %.2f\n', gaborCST);
fprintf(' - nTrials = %d\n', nTrials);
fprintf(' - iModelB_sim = %d (see updated Model B definitions)\n', iModelB_sim);
fprintf(' - %d params: %s\n', nParams, str_IN);
fprintf(' - Criterion (z unit)=%.1f \n', Cz_true);
fprintf(' - Whitening strength (lambda): %.1f \n\n', lambda_whiten);

% ---------------------- Define IO name & folders ---------------------- %

nameIO = sprintf('IO_cN%.0f_cG%.0f_nT%s_Nm%s_Na%s_Ns%s_Cz%.1f_whiten%.1f_%d%d_B%d', ...
    noiseCST*100, gaborCST*100, format_num2exp(nTrials), ...
    format_num2exp(Nmul_true), format_num2exp(Nadd_true), format_num2exp(Nshared_true), ...
    Cz_true, lambda_whiten, nORI, nSF, iModelB_sim);

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

fprintf('%s: Filter banks (nORI=%d, nSF=%d) created and saved.\n\n', datetime('now'), length(filtersOri_all), length(noise.filtersSF_all))

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
% template_true = template_true / max(template_true(:)) * 0.2;

% Save "truth" (for IO template in Model A = 2)
save(sprintf('%s/truth.mat', nameFolder_Data_OOD_IO), '*_true');
fprintf('\n\n%s: Truth (including the ideal template) created and saved.\n\n', datetime('now'))

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
DV_target_sim_allT= nan(nTrials, 1);
% IV_noise_sim_allT = nan(nTrials, 1);

mask = stim.mask;
ratio_base = noise.ratio_base;
ratio_gaborInTgt = noise.ratio_gaborInTgt;

%% ---------------------- Simulate stimuli & IV ------------------------- %
fprintf('%s: started running %d pairs.\n\n', datetime('now'), nPairs)

parfor iPair = 1:nPairs

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

    % Internal variable from energy × template_true
    %---------------------------%
    DV_target = fxn_getIV_v3(e3D_target, template_true, convolveType_true, IVType_true, flag_permT, [1, nORI]);
    %---------------------------%
    assert(~isnan(DV_target), 'IV_target is NaN');
    DV_target_sim_allT(iPair) = DV_target;

    % Data matrix: col 11 = gabor contrast
    if iPRS_allT(iPair) == 1
        dataMatrix(iPair, :) = [nan, nan, nan, nan, iLocComb, nan, nan, nan, nan, nan, stim.gaborCST];
    else
        dataMatrix(iPair, :) = [nan, nan, nan, nan, iLocComb, nan, nan, nan, nan, nan, 0];
    end

    % Store energy
    e3D_target_allT(iPair, :, :) = e3D_target;
    % e3D_noise_allT(iPair, :, :) = e3D_noise;

end % iPair

fprintf('\n\n%s: Simulation of all pairs done.\n\n', datetime('now'))

% Copy pass A -> pass B
dataMatrix(nPairs+1:end, :, :) = dataMatrix(1:nPairs, :, :);
e3D_target_allT(nPairs+1:end,:,:) = e3D_target_allT(1:nPairs, :, :);
DV_target_sim_allT(nPairs+1:end,:) = DV_target_sim_allT(1:nPairs, :);

assert(sum(e3D_target_allT(nPairs+1:end,:,:) - e3D_target_allT(1:nPairs,:,:), 'all') == 0);
% assert(sum(e3D_noise_allT(nPairs+1:end,:,:) - e3D_noise_allT(1:nPairs,:,:), 'all') == 0);

fprintf('%s: Pass A copied to pass B done.\n\n', datetime('now'))

%% --------------------- Save energy for compIV ------------------------- %
% Use patchMode = 'T' in compIV; this file must contain both target + noise
% energy as well as noise, filters, stim, nBins.
save(sprintf('%s/energy_T_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF), 'e3D_target_allT', 'noise', 'filtersOri_all', 'stim', 'nBins');

fprintf('%s: Simulated 3D energy saved (%d trials).\n\n', datetime('now'), nTrials)

%% ---------------------- Internal noise & responses -------------------- %
% For now, use DV from target-only energy to drive responses
DVclean_sim_allT = DV_target_sim_allT;

% Independent (pass-specific) SD per trial
sigma_priv_allT = sqrt((DVclean_sim_allT .* Nmul_true).^2 + Nadd_true^2);

% Draw independent z per trial (per pass)
z_priv_allT = randn(size(DVclean_sim_allT));

% Shared SD per trial (often constant)
sigma_shared_allT  = Nshared_true * ones(size(DVclean_sim_allT));
pairIDs = unique(iPair_allT); % Draw one shared z per PAIR, reuse for both passes
z_sh_perPair = randn(numel(pairIDs), 1);
[~, idxPair] = ismember(iPair_allT, pairIDs); % Map pair -> z_sh for each trial
z_sh_allT = z_sh_perPair(idxPair);

% Construct noisy DV / DV_noisy
DVnoisy_sim_allT = DVclean_sim_allT + sigma_shared_allT .* z_sh_allT + sigma_priv_allT .* z_priv_allT;

% Total predicted SD if you need it for anything
sigma_pred_allT = sqrt(sigma_priv_allT.^2 + sigma_shared_allT.^2);

% Correlation between the noisy DV of two passes
% (assuming first half and second half of trials are the same trials)
corrPass = corr(DVnoisy_sim_allT(1:numel(DVnoisy_sim_allT)/2), DVnoisy_sim_allT(numel(DVnoisy_sim_allT)/2+1:numel(DVnoisy_sim_allT)));

fprintf('%s: Internal noise sampled and added.\n\n', datetime('now'))

%% Set a true criterion given the true criterion in z-score
% Old way: Titrate a criterion to reach target accuracy pC_titrate
%---------------------------%
% fxn_loss_pC_ = @(criterion_potential) fxn_loss_pC(criterion_potential, pC_titrate, DVnoisy_sim_allT, iPRS_allT);
% %---------------------------%
% criterion_true = fmincon(fxn_loss_pC_, median(DVnoisy_sim_allT), [], [], [], [], min(DVnoisy_sim_allT), max(DVnoisy_sim_allT));

% Candidate criteria in DV space:
DV_sorted = sort(unique(DVnoisy_sim_allT(:)));
criterion_grid = [-Inf; (DV_sorted(1:end-1) + DV_sorted(2:end))/2; Inf];

loss_grid = nan(size(criterion_grid));
c_z_grid = nan(size(criterion_grid));

for j = 1:numel(criterion_grid)
    c_try = criterion_grid(j);
    [loss_grid(j), c_z_grid(j)] = criterion_loss_cz(c_try, DVnoisy_sim_allT, iPRS_allT, Cz_true);
end
% Choose the criterion in DV unit
[~, idx_best] = min(loss_grid);
criterion_DV_true = criterion_grid(idx_best);

% Binary responses
resp_allT = DVnoisy_sim_allT > criterion_DV_true;

% % Lapses (if any)
% lapse_mask = rand(size(DVnoisy_sim_allT)) < lapseRate_true;
% resp_allT(lapse_mask) = ~resp_allT(lapse_mask);

% pYES per trial from SDT mapping
pYES_pred_allT = 1 - normcdf(criterion_DV_true, DVclean_sim_allT, sigma_pred_allT);

% Fill behavior columns in dataMatrix
dataMatrix(:, 1) = (1:nTrials)'; % trial index
dataMatrix(:, 2) = iSess_allT; % session index
dataMatrix(:, 5) = iLocComb; % location index
dataMatrix(:, 6) = iPRS_allT; % PRS/ABS
dataMatrix(:, 7) = iPass_allT; % pass index
dataMatrix(:, 8) = iPair_allT; % pair index
dataMatrix(:, 9) = resp_allT; % response

fprintf('%s: Criterion calculated and responses generated.\n\n', datetime('now'))

%% ---------------------- Behavioral metrics ---------------------------- %
pYES_sim = mean(resp_allT == 1);
pHit_sim = sum((iPRS_allT == 1) & (resp_allT == 1)) / sum(iPRS_allT == 1);
pFA_sim = sum((iPRS_allT == 0) & (resp_allT == 1)) / sum(iPRS_allT == 0);
pC_sim = mean( (iPRS_allT==1 & resp_allT==1) | (iPRS_allT==0 & resp_allT==0) );

%--------------------------------------%
[dprime_sim, criterion_z_sim] = SX_sim06_SDT(pHit_sim, pFA_sim);
%--------------------------------------%

respC_sim = nan(nPairs, 1);
for iPairUnik = 1:nPairs
    respAB = resp_allT(iPair_allT == iPairUnik);
    respC_sim(iPairUnik) = (respAB(1) == respAB(2));
end
pA_sim = mean(respC_sim);

metrics_sim = [dprime_sim, criterion_z_sim, pC_sim, pHit_sim, pFA_sim, pA_sim, pYES_sim];

% Compare true criterion (IV space) vs SDT-derived criterion (converted)
% c_IV_all = median(DVclean_sim_allT) + criterion_z .* sigma_pred_allT;
fprintf('   True criterion (z unit) = %.1f, (DV unit)=%.1f\n', Cz_true, criterion_DV_true);
fprintf('   Simulated criterion (z unit) = %.1f\n', criterion_z_sim);

% Save behavioral measures in behavMeas.mat (as expected by compIV)
save(sprintf('%s/behavMeas.mat', nameFolder_Data_OOD_IO), 'dataMatrix', 'metrics_sim');

fprintf('\n%s: Behavioral data saved.\n\nReady for template generation\n\n', datetime('now'))

%% -------------------------- Plotting DV dist and behav metrics------------------------- %
if flag_plotDist
figure('Position', [0 200 600 1e3]);
subplot(2,1,1); hold on;
histogram(DVnoisy_sim_allT(iPRS_allT == 1), 'FaceColor', 'r', 'DisplayName', 'Signal Present', 'normalization', 'probability');
histogram(DVnoisy_sim_allT(iPRS_allT == 0), 'FaceColor', 'b', 'DisplayName', 'Signal Absent', 'normalization', 'probability');
xline(criterion_DV_true, 'LineWidth', 2, 'DisplayName', 'True criterion');
xlim([min(DVnoisy_sim_allT), max(DVnoisy_sim_allT)]);
ylabel('Proportion');
legend('show', 'location', 'best');

subplot(2,1,2); hold on;
plot(DVclean_sim_allT, resp_allT, 'ro', 'DisplayName', 'Binary response');
plot(DVclean_sim_allT, pYES_pred_allT, 'k+', 'DisplayName', 'Pred pYES');
xline(criterion_DV_true, 'LineWidth', 2, 'DisplayName', 'True criterion');
xlabel('IV'); ylabel('pYES');
yline(0.5, 'k--');
xlim([min(DVnoisy_sim_allT), max(DVnoisy_sim_allT)]);
metrics_sim_ = metrics_sim; metrics_sim_(3:end) = metrics_sim_(3:end)*100;

sgtitle(sprintf(['CI_{95}=[%.1f, %.1f], Median=%.1f\n[TRUE] gN=%.2f, gC=%.2f, Nmul=%s, Nadd=%s, Nshared=%s, crit=%.1f (%.1f)' ...
    '\n[MEASURED] pYES=%.2f, pC=%.2f, pHit=%.2f, pFA=%.2f, pA=%.2f, corrPass=%.2f'], ...
    round(quantile(DVnoisy_sim_allT, [.05, .95, .5]), 1), noiseCST, gaborCST, ...
    format_num2exp(Nmul_true), format_num2exp(Nadd_true), format_num2exp(Nshared_true), Cz_true, criterion_z_sim, ...
    pYES_sim, pC_sim, pHit_sim, pFA_sim, pA_sim, corrPass), ...
    'fontsize', 12);

saveas(gcf, sprintf('%s/0PerfHist.jpg', nameFolder_Figures_perSubj));
close all;
fprintf('%s: Distribution of DV drawn.\n\n', datetime('now'))
end

%% ---------------- Run compIV and fitNOM on this IO --------------------- %
for iModelA_fit = iModelA_sim_allCond
    % Step 1: compute IVs, templates, and test-set metrics
    %------------------------------------%
    OOD_NOM_Trialwise_compIV({nameIO, criterion_DV_true}, iLocComb, lambda_whiten, iModelA_fit, nIter, nJob, iJob)
    %------------------------------------%

    for iModelB_fit = iModelB_sim_allCond
        % Step 2: fit NOM parameters and predict metrics
        %------------------------------------%
        OOD_NOM_Trialwise_fitNOM({nameIO, criterion_DV_true}, iLocComb, iModelA_fit, iModelB_fit, nIter, nJob, iJob)
        %------------------------------------%
    end
end

% Play sound to indicate end of analysis
fs = 44100;              % sampling rate
dur = 0.25;              % duration of each beep (seconds)
pauseDur = 0.1;          % silence between beeps
t = 0:1/fs:dur;
beep1 = sin(2*pi*400*t);
beep2 = sin(2*pi*600*t);
beep3 = sin(2*pi*800*t);
silence = zeros(1, round(fs*pauseDur));
% sound([beep1 silence beep2 silence beep3], fs)

%% --------------------- Clean-up temporary files ------------------------ %
clear *allT e2D* e3D* dataMatrix;

% Remove the energy file to save space (compIV has already loaded it)
delete(sprintf('%s/energy_T_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF));
fprintf('%s: Energy file removed to save space.\n\n', datetime('now'))

%% -------------------- End timing -------------------- %%
time_end = datetime('now');
fprintf('%s: Simulation done.\n\n', time_end)
elapsed = time_end - time_start;
fprintf('Time used: %s\n\n\n\n', char(elapsed));

end

%% helper
function [loss, c_z, pHit, pFA] = criterion_loss_cz(k, DV_noisy, iPRS_allT, target_cz)

nPRS = sum(iPRS_allT == 1);
nABS = sum(iPRS_allT == 0);
eps_ = 0.5 / min(nPRS, nABS);

resp = double(DV_noisy > k);

pHit = sum(resp == 1 & iPRS_allT == 1) / nPRS;
pFA  = sum(resp == 1 & iPRS_allT == 0) / nABS;

% clip rates to avoid infinities
pHit = min(max(pHit, eps_), 1 - eps_);
pFA  = min(max(pFA,  eps_), 1 - eps_);

c_z = -0.5 * (norminv(pHit) + norminv(pFA));

loss = (c_z - target_cz).^2;
end