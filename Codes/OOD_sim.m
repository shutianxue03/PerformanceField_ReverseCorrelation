function OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, iModelB_sim, nIter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script name: OOD_sim.m
% Adapted by Shutian Xue on 05/27/2026
%
% This script simulates observer responses based on a trial-wise
% noisy observer model (NOM), and then runs the updated
% OOD_NOM_Trialwise_compDV / OOD_NOM_Trialwise_fitNOM pipeline
% on the simulated data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation input parameters
iModelA_fit_all = [1]; % DO NOT CHANGE! 1 = RC-derived template (Model A), 2=ideal template; 3=permuted template
iModelB_fit_all = 1:7; %1=full, 2=No Nmul, 3=No Nadd, 4=No Nshared, 5=Nmul-only, 6=Nadd-only, 7=Nshared-only, 8=criterion-only

% Set hyperparameters
hyperparams = [];

hyperparams.nBasisORI = 3;%3:2:9;
hyperparams.nBasisSF = 3;%3:2:9;
hyperparams.basisWidthORI = .6;%.3:.2:.9;
hyperparams.basisWidthSF = .4;%.3:.2:.9;
hyperparams.asymSF_rightLeftRatio = 0.9;%1.1:.2:1.5;
hyperparams.Ridge = 100;%[1,10,100];

% Enforce reduced-model ground truth by zeroing excluded IN terms.
switch iModelB_sim
    case 1  % Full: Nmul + Nadd + Nshared
        % keep all as provided
    case 2  % No Nmul
        Nmul_true = 0;
    case 3  % No Nadd
        Nadd_true = 0;
    case 4  % No Nshared
        Nshared_true = 0;
    case 5  % Nmul only
        Nadd_true = 0;
        Nshared_true = 0;
    case 6  % Nadd only
        Nmul_true = 0;
        Nshared_true = 0;
    case 7  % Nshared only
        Nmul_true = 0;
        Nadd_true = 0;
    otherwise
        error('Unknown iModelB_sim = %d', iModelB_sim);
end

% Basic setup -%
clc; close all;
warning off;
set(0, 'DefaultFigureVisible', 'off')

time_start = datetime('now');
fprintf('\n\n%s: Simulation starts\n\n', datetime('now'))

% Set RNG for reproducibility
rng(1);

% Add paths for custom functions using the function's own location.
code_root = fileparts(mfilename('fullpath'));
addpath(genpath(code_root));

% Global RC / NOM settings (defines nORI, nSF, nBins, folders, etc.)
SX_RC1_setting;
nORI = nORI; % nORI is defined in SX_RC1_setting; repeated here just for parfor loop to work
nSF = nORI;

% Number of bootstraps for compDV / fitNOM
nJob = 1;
iJob = 1;
pC_filter = [.6, .8]; % Only proceed with compDV/fitNOM if simulated pC falls within this range; otherwise, discard this simulation and try again with different random seed or parameters.
ORI_bound = [5, 14]; % orientation window passed to fxn_getDV_v3

str_templateType_true = 'raw'; % see 'fxn_getTemplate' for other options
str_DVType_true = 'sum'; % see 'fxn_getDV_v3' for other options
str_convolveType_true = 'dot'; % see 'fxn_getDV_v3' for other options
flag_normDV = 1; % 1=normalize simulated DV by global SD before adding internal noise
flag_permT = 0; %1=permute the input template per trial
flag_plotDist = 1;
if strcmp(str_envir, 'HPC'), flag_plotDist = 0; end % don't plot when running on HPC
iLocComb = 1; % use single-location index (e.g., fovea) for IO

% Print info --%
fprintf(' - nIter (for compDV/fitNOM) = %d\n', nIter);
fprintf(' - noiseCST = %.2f\n', noiseCST);
fprintf(' - gaborCST = %.2f\n', gaborCST);
fprintf(' - nTrials = %d\n', nTrials);
fprintf(' - Simulated with ModelB = %d \n', iModelB_sim);
fprintf(' - True IN params: Nmul=%.3g, Nadd=%.3g, Nshared=%.3g\n', Nmul_true, Nadd_true, Nshared_true);
fprintf(' - Normalize simulated DV (global SD): %d\n', flag_normDV);
fprintf(' - Fitted with ModelA = %s (1=Data-derived template; 2=ideal template)\n', strjoin(string(iModelA_fit_all), ' '));
fprintf(' - Fitted with ModelB = %s (1=full, 2=No Nmul, 3=No Nadd, 4=No Nshared, 5=Nmul-only, 6=Nadd-only, 7=Nshared-only, 8=criterion-only)\n', strjoin(string(iModelB_fit_all), ' '));
fprintf(' - SDT Criterion=%.1f \n', cSDT_true);
fprintf(' - Contribution of criterion loss =%.1f \n', C_contribution);
fprintf(' - Regression type (1=Univariate; 2=Multi+smoothing): %d \n\n', flag_regressType);

% Define IO name & folders %
% Define the IO name
if isempty(hyperparams) % if basis functions' hyperparameters are NOT specified
    nameIO = sprintf('IO_Bsim%d_cN%.0f_cG%.0f_nT%s_Nm%.1f_Na%.1f_Ns%.1f_cSDT%.1f', ...
        iModelB_sim, ...
        noiseCST*100, gaborCST*100, format_num2exp(nTrials), ...
        Nmul_true, Nadd_true, Nshared_true, cSDT_true);
else % if hyperparameters are specified, include them in the IO name
    nameIO = sprintf('IO_Bsim%d_Bori%d_%.1f_Bsf%d_%.1f_RL%.1f_ridge%d_cN%.0f_cG%.0f_nT%s_Nm%.1f_Na%.1f_Ns%.1f_cSDT%.1f', ...
        iModelB_sim, ...
        hyperparams.nBasisORI, hyperparams.basisWidthORI, hyperparams.nBasisSF, hyperparams.basisWidthSF, hyperparams.asymSF_rightLeftRatio, hyperparams.Ridge, ...
        noiseCST*100, gaborCST*100, format_num2exp(nTrials), ...
        Nmul_true, Nadd_true, Nshared_true, cSDT_true);
end

% Folder to save IO data (energy + behav)
nameFolder_Data_OOD_IO = sprintf('%s/%s', nameFolder_Data_OOD, nameIO);

% Folder to save NOM trial-wise fits
nameFolder_Data_NOM_IO = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, nameIO);

% Folder to save figures
nameFolder_Figures_perSubj = sprintf('%s/IO/%s', nameFolder_Figures, nameIO);

fprintf('\nIO name: %s\n', nameIO);
fprintf('Data_OOD folder: %s\n', nameFolder_Data_OOD_IO);
fprintf('Data_NOM folder: %s\n\n', nameFolder_Data_NOM_IO);

%% Stimulus / noise parameters
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

%% Create Gabor filters used for energy computation
[filter_sin, filter_cos] = SX_sim02_setFilters(stim, noise.filtersSF_all, filtersOri_all, fxn_getSigma_SPdomain, 0);
fprintf('%s: Filter banks (nORI=%d, nSF=%d) created and saved.\n\n', datetime('now'), length(filtersOri_all), length(noise.filtersSF_all))

%% Create the ground-truth template (ORI x SF)
% Make sure this part is exactly the same as how ideal template is created
% in OOD_NOM_compDV_A12
template_gabor_true = exp_CreateGabor(stim, stim.gaborCST);
%--------------------------------------------%
template_true = SX_RC4_Energy_parfor(stim.mask, {template_gabor_true}, filter_sin, filter_cos);
%--------------------------------------------%
template_true = squeeze(template_true); % remove singleton dim
%--------------------------------------------%
template_true = fxn_getTemplate(template_true, str_templateType_true, 0);
%--------------------------------------------%

% Do L2 normalization
% ENSURE this step is consistent with OOD_xx_compDV_A12 when creating the ideal template
template_true = template_true / norm(template_true(:));

%% Preallocate sim arrays
nMetrics = 11;
iPRS_allT_OnePass = [ones(nPairs/2, 1); zeros(nPairs/2, 1)]; % 1=PRS, 0=ABS
iPRS_allT = [iPRS_allT_OnePass; iPRS_allT_OnePass];
iPair_allT = [1:nPairs, 1:nPairs]';
iPass_allT = [ones(nPairs,1); ones(nPairs,1)*2];
iSess_allT = repmat(1:nSess, 1, nTrialsPerSess)';

dataMatrix = nan(nTrials, nMetrics);
e3D_target_allT = nan(nTrials, nORI, nSF);
% DV_target_sim_allT= nan(nTrials, 1);

mask = stim.mask;
ratio_base = noise.ratio_base;
ratio_gaborInTgt = noise.ratio_gaborInTgt;

%% Simulate stimuli & DV
fprintf('%s: started running %d pairs.\n\n', datetime('now'), nPairs)

parfor iPair = 1:nPairs

    % Filtered noise
    filtered_noise = exp_CreateFilteredNoise(noise);

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
    e3D_target = SX_RC4_Energy_parfor(mask, {patch_target}, filter_sin, filter_cos);
    e3D_target = squeeze(e3D_target);

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
% DV_target_sim_allT(nPairs+1:end,:) = DV_target_sim_allT(1:nPairs, :);

assert(sum(e3D_target_allT(nPairs+1:end,:,:) - e3D_target_allT(1:nPairs,:,:), 'all') == 0);
% assert(sum(e3D_noise_allT(nPairs+1:end,:,:) - e3D_noise_allT(1:nPairs,:,:), 'all') == 0);

fprintf('%s: Pass A copied to pass B done.\n\n', datetime('now'))

%% Compute simulated DV
DV_target_sim_allT = fxn_getDV_v3(e3D_target_allT, template_true, str_convolveType_true, str_DVType_true, flag_permT, ORI_bound);

assert(~any(isnan(DV_target_sim_allT)), 'DV_target contains NaN');

%% Normalize simulated DV before internal-noise sampling
if flag_normDV
    DV_scaleFactor_sim = std(DV_target_sim_allT);
    if ~isfinite(DV_scaleFactor_sim) || DV_scaleFactor_sim < 1e-8
        DV_scaleFactor_sim = max(abs(DV_target_sim_allT));
    end
    if ~isfinite(DV_scaleFactor_sim) || DV_scaleFactor_sim < 1e-8
        DV_scaleFactor_sim = 1;
    end
else
    DV_scaleFactor_sim = 1;
end

fprintf('%s: Sim DV normalization scale = %.6g; true params used on normalized DV Nmul=%.3g, Nadd=%.3g, Nshared=%.3g\n\n', ...
    datetime('now'), DV_scaleFactor_sim, Nmul_true, Nadd_true, Nshared_true);

% Empty placeholders kept ONLY for backward compatibility with downstream
% code that may still expect these variable names. Do NOT use them in
% any computation.
mu_cov_DV       = [];
sigma_cov_DV    = [];
W_white_DV      = [];
beta_true_basis = [];
basisOpts_true  = [];

fprintf('%s: DV computed in the raw-energy space.\n\n', datetime('now'))

%% Add internal variability and generate responses
% For now, use DV from target-only energy to drive responses
DVclean_sim_allT = DV_target_sim_allT ./ DV_scaleFactor_sim;

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

fprintf('%s: Internal noise sampled and added.\n\n', datetime('now'))

%% Set a true criterion given the true SDT criterion
% Candidate criteria in DV space:
DV_sorted = sort(unique(DVnoisy_sim_allT(:)));
cDV_grid = [-Inf; (DV_sorted(1:end-1) + DV_sorted(2:end))/2; Inf];
loss_grid = nan(size(cDV_grid));
cSDT_grid = nan(size(cDV_grid));
parfor iC = 1:numel(cDV_grid)
    c_try = cDV_grid(iC);
    [loss_grid(iC), cSDT_grid(iC)] = fxn_loss_cSDT(c_try, DVnoisy_sim_allT, iPRS_allT, cSDT_true);
end

% Choose the criterion in DV unit
[~, idx_best] = min(loss_grid);
criterion_DV_true = cDV_grid(idx_best);

fprintf('%s: Defined a true criterion in DV unit given the true SDT criterion.\n\n', datetime('now'))

%% Simulate responses and calculate pYES
% Binary responses
resp_allT = DVnoisy_sim_allT > criterion_DV_true;

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

fprintf('%s: Responses generated.\n\n', datetime('now'))

%% Behavioral metrics
pYES_sim = mean(resp_allT == 1);
pHit_sim = sum((iPRS_allT == 1) & (resp_allT == 1)) / sum(iPRS_allT == 1);
pFA_sim = sum((iPRS_allT == 0) & (resp_allT == 1)) / sum(iPRS_allT == 0);
pC_sim = mean( (iPRS_allT==1 & resp_allT==1) | (iPRS_allT==0 & resp_allT==0) );

[dprime_sim, cSDT_sim] = SX_sim06_SDT(pHit_sim, pFA_sim);

respC_sim = nan(nPairs, 1);
for iPairUnik = 1:nPairs
    respAB = resp_allT(iPair_allT == iPairUnik);
    respC_sim(iPairUnik) = (respAB(1) == respAB(2));
end
pA_sim = mean(respC_sim);

metrics_sim = [dprime_sim, cSDT_sim, pC_sim, pHit_sim, pFA_sim, pA_sim, pYES_sim];

%% Plotting DV dist and behav metrics
% Create figure folder only when plotting is enabled (never on HPC).
if flag_plotDist
    mkdir(nameFolder_Figures_perSubj);
    NOMplot_dist
end

%%
% Stop if simulated pC is too low or too high
if (pC_sim >= pC_filter(1)) && (pC_sim <= pC_filter(2))
    fprintf('\n\n ** Simulated pC=%.2f, within the range [%.2f, %.2f] ** \n\n', pC_sim, pC_filter)

    % Create output folders only after the accuracy gate is passed.
    if isempty(dir(nameFolder_Data_OOD_IO)), mkdir(nameFolder_Data_OOD_IO); end
    if isempty(dir(nameFolder_Data_NOM_IO)), mkdir(nameFolder_Data_NOM_IO); end

    % Save only variables used by simPlot4_VaryOneDim.m in one atomic write.
    save(sprintf('%s/truth.mat', nameFolder_Data_OOD_IO), '*_true', 'mu_cov_DV', 'sigma_cov_DV', 'W_white_DV', 'beta_true_basis', 'flag_normDV', 'DV_scaleFactor_sim');
    fprintf('\n\n%s: truth.mat saved (template_true, criterion_DV_true).\n\n', datetime('now'))

    % Save the target-energy tensor
    save(sprintf('%s/energy_T_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF), 'e3D_target_allT', 'noise', 'filtersOri_all', 'stim', 'nBins');
    fprintf('%s: Simulated 3D energy saved (%d trials).\n\n', datetime('now'), nTrials)

    % Save behavioral measures in behavMeas.mat (as expected by compDV)
    save(sprintf('%s/behavMeas.mat', nameFolder_Data_OOD_IO), 'dataMatrix', 'metrics_sim', 'iPRS_allT', 'iPass_allT', 'iPair_allT', 'resp_allT', 'pC_filter');
    fprintf('\n%s: Behavioral data saved.\n\nReady for template generation\n\n', datetime('now'))

    %% Step 1: Estimate the template and compute DV
    % Step 1: compute DVs, templates, and test-set metrics
    %----------------------------%
    OOD_NOM_Trialwise_compDV_A12({nameIO, criterion_DV_true}, iLocComb, nIter, nJob, iJob, flag_normDV, hyperparams)
    %----------------------------%

    %% Step 2: fit NOM parameters and predict metrics
    for iModelA_fit = iModelA_fit_all

        for iModelB_fit = iModelB_fit_all
            %----------------------------%
            OOD_NOM_Trialwise_fitNOM({nameIO, criterion_DV_true}, iLocComb, iModelA_fit, iModelB_fit, nIter, nJob, iJob)
            %----------------------------%

        end % for iModelB_fit

    end % for iModelA_fit
else
    fprintf('\n\n ** Simulated pC=%.2f, OUT OF the range [%.2f, %.2f] ** \n\n', pC_sim, pC_filter)
end % if


%% Clean-up temporary files
clear *allT e2D* e3D* dataMatrix;

% Remove the energy file to save space (compDV has already loaded it)
nameFile_energy = sprintf('%s/energy_T_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF);
if exist(nameFile_energy, 'file')
    delete(nameFile_energy);
    fprintf('%s: Energy file removed to save space.\n\n', datetime('now'))
end

%%  End timing
time_end = datetime('now');
fprintf('%s: Simulation done.\n\n', time_end)
elapsed = time_end - time_start;
fprintf('Time used: %s\n\n\n\n', char(elapsed));

end % end of the OOD_sim function

%% helper
function [loss, cSDT, pHit, pFA] = fxn_loss_cSDT(k, DV_noisy, iPRS_allT, target_cSDT)
% This function computes the loss for a given criterion k in DV space, based on the simulated noisy DV and the true labels (PRS/ABS). The loss is defined as the squared error between the computed cSDT and the target cSDT. The function also returns the computed cSDT, pHit, and pFA for that criterion.
nPRS = sum(iPRS_allT == 1);
nABS = sum(iPRS_allT == 0);
eps_ = 0.5 / min(nPRS, nABS);

resp = double(DV_noisy > k);

pHit = sum(resp == 1 & iPRS_allT == 1) / nPRS;
pFA  = sum(resp == 1 & iPRS_allT == 0) / nABS;

% clip rates to avoid infinities
pHit = min(max(pHit, eps_), 1 - eps_);
pFA  = min(max(pFA,  eps_), 1 - eps_);

cSDT = -0.5 * (norminv(pHit) + norminv(pFA));

loss = (cSDT- target_cSDT).^2;
end
