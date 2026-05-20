% function OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, lambda_whiten, iModelB_sim, nIter, nBasisORI, basisWidthORI, nBasisSF, basisWidthSF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script name: OOD_sim.m
% Adapted by Shutian Xue on 03/27/2026
%
% This script simulates observer responses based on a trial-wise
% noisy observer model (NOM), and then runs the updated
% OOD_NOM_Trialwise_compDV / OOD_NOM_Trialwise_fitNOM pipeline
% on the simulated data.
%
% Notes:
% - Gabor SD is normalized by SF only when creating the filters,
% not when creating the Gabor patches.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation input parameters
iModelA_fit_all = [1:2]; % DO NOT CHANGE! 1 = RC-derived template (Model A), 2=ideal template; 3=permuted template
iModelB_fit_all = iModelB_sim; %1=full, 2=No Nmul, 3=No Nadd, 4=No Nshared, 5=Nmul-only, 6=Nadd-only, 7=Nshared-only, 8=criterion-only

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

templateType_true = 1; % 1 = raw;
DVType_true = 1; % 1 = sum of dot product;
convolveType_true = 1; % 1 = dot product; 2 = convolution (for fxn_getDV_v3)
flag_permT = 0; %1=permute the input template per trial
eps_whiten = 1e-3; % floor for whitening eigenvalues
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
fprintf(' - Fitted with ModelA = %s (1=Data-derived template; 2=ideal template)\n', strjoin(string(iModelA_fit_all), ' '));
fprintf(' - Fitted with ModelB = %s (1=full, 2=No Nmul, 3=No Nadd, 4=No Nshared, 5=Nmul-only, 6=Nadd-only, 7=Nshared-only, 8=criterion-only)\n', strjoin(string(iModelB_fit_all), ' '));
fprintf(' - SDT Criterion=%.1f \n', cSDT_true);
fprintf(' - Contribution of criterion loss =%.1f \n', C_contribution);
fprintf(' - Whitening strength (lambda): %.1f \n', lambda_whiten);
fprintf(' - Regression type (1=Univariate; 2=Multi+smoothing): %d \n\n', flag_regressType);

% Define IO name & folders %
% Define the IO name
nameIO = sprintf('IO_Bsim%d_cN%.0f_cG%.0f_nT%s_Nm%.1f_Na%.1f_Ns%.1f_cSDT%.1f_whiten%.1f_bORI%d_%.1f_bSF%d_%.1f', ...
    iModelB_sim, ...
    noiseCST*100, gaborCST*100, format_num2exp(nTrials), ...
    Nmul_true, Nadd_true, Nshared_true, cSDT_true, ...
    lambda_whiten, nBasisORI, basisWidthORI, nBasisSF, basisWidthSF);

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

%% True 2D template (ORI×SF)
template_gabor_true = exp_CreateGabor(stim, stim.gaborCST);
template_true = SX_RC4_Energy_parfor(stim.mask, {template_gabor_true}, filter_sin, filter_cos);
template_true = squeeze(template_true);
template_true = fxn_getTemplate(template_true, templateType_true, 0);

% Normalize template (Unit L2-norm)
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
DV_target_sim_allT= nan(nTrials, 1);

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

%% Compute DV in a fixed transformed space
% Fixed transform estimated from ABS trials:
% raw energy -> global ABS z-score (one mean/SD per channel) -> optional whitening.

% (1) Build fixed transform from ABS trials
% to obtain whitening matrix (Q) + normalization stats (mu, sigma) for each channel.
Tfix = fxn_buildFixedTransformFromABS(e3D_target_allT, iPRS_allT, lambda_whiten, eps_whiten);

% (2) Apply fixed transform to energy of all trials, to prepare for DV computation
[nTrials, nORI_local, nSF_local] = size(e3D_target_allT);
e3D_allT = reshape(e3D_target_allT, [nTrials, nORI_local * nSF_local]);
% To z-score
e3D_allT_z = (e3D_allT - Tfix.mu_abs) ./ Tfix.sigma_abs;
% To whiten
e3D_allT_whitened = e3D_allT_z * Tfix.Q;
e3D_forDV = reshape(e3D_allT_whitened, [nTrials, nORI_local, nSF_local]);

%% Sanity check: fixed transform correctness
idxABS_check = (iPRS_allT == 0);
X_abs = reshape(e3D_forDV(idxABS_check, :, :), [sum(idxABS_check), nORI_local * nSF_local]);
X_prs = reshape(e3D_forDV(~idxABS_check, :, :), [sum(~idxABS_check), nORI_local * nSF_local]);

figure('Position', [100 100 1000 400]);

% --- Check 2: ABS-trial covariance ≈ identity ---
C_abs = cov(X_abs);                          % should be ~I after whitening
diag_vals    = diag(C_abs);
offdiag_vals = C_abs(~logical(eye(size(C_abs, 1))));

subplot(1, 2, 1);
histogram(diag_vals, 20, 'FaceColor', [0.2 0.5 0.8]);
xline(1, 'r--', 'LineWidth', 1.5);
xlabel('Diagonal value'); ylabel('Count');
title(sprintf('Diag: mean=%.3f, sd=%.3f\nShould concentrate at 1', mean(diag_vals), std(diag_vals)));

subplot(1, 2, 2);
histogram(offdiag_vals, 50, 'FaceColor', [0.7 0.4 0.2]);
xline(0, 'r--', 'LineWidth', 1.5);
xlabel('Off-diagonal value'); ylabel('Count');
title(sprintf('Off-diag: mean=%.3f, sd=%.3f\nShould concentrate at 0', mean(offdiag_vals), std(offdiag_vals)));

sgtitle('ABS-trial Covariance Check', 'Interpreter', 'none');

%% Sanity check: Compute eigenvalues of actual ABS covariance in whitened space
Xz_abs = reshape(e3D_z(iPRS_allT==0,:,:), [sum(iPRS_allT==0), nORI_local*nSF_local]);
Xz_abs_c = bsxfun(@minus, Xz_abs, Tfix.mu_cov_abs);
[V_check, D_check] = eig(cov(Xz_abs_c, 1));
d_check = diag(D_check);

% After whitening, variance in each eigendirection = d_i / max(d_i, eps_whiten)
post_whiten_var = d_check ./ max(d_check, eps_whiten);

figure;
subplot(1,2,1);
histogram(d_check, 50); xline(eps_whiten, 'r--');
xlabel('Eigenvalue of Sigma'); title('Many below eps\_whiten floor (red)');

subplot(1,2,2);
histogram(post_whiten_var, 50); xline(1, 'r--');
xlabel('Post-whitening variance per eigendirection');
title('Should cluster at 1 (floored dirs < 1)');


%%
% (3) Convert template into the same transformed space.
sigma_vec = Tfix.sigma_abs(:);
sigma_vec(~isfinite(sigma_vec) | sigma_vec < 1e-8) = 1e-8;
template_z = template_true(:) .* sigma_vec;
if Tfix.useWhiten
    t_trans = Tfix.Q \ template_z;
else
    t_trans = template_z;
end
template_forDV = reshape(t_trans, size(template_true));

%% Visualize template in the original vs transformed space
figure,
% Scale the transformed template back to the same max as the original template for better visualization; this scaling does not affect the DV since it's just a constant factor.
template_forDV_scaled = template_forDV / max(template_forDV(:)) * max(template_true(:));
margORI_true = mean(template_true, 2);
margORI_forDV = mean(template_forDV, 2);
margORI_forDV_scaled = margORI_forDV / max(margORI_forDV(:)) * max(margORI_true(:));
margSF_true = mean(template_true, 1);
margSF_forDV = mean(template_forDV, 1);
margSF_forDV_scaled = margSF_forDV / max(margSF_forDV(:)) * max(margSF_true(:));

% Row 1: 2D template;
subplot(3,3,1), imagesc(template_true), colorbar, title('Ground-truth template (no scaling)')
subplot(3,3,2), imagesc(template_forDV_scaled), colorbar, title('Transformed template (no further scaling)')
subplot(3,3,3), imagesc(template_true - template_forDV_scaled), colorbar, title('Difference')
% Row 2: ORI tuning curve (averaged across SF)
subplot(3,3,5), hold on
plot(axis_tuning{1}, margORI_true, '-r')
plot(axis_tuning{1}, margORI_forDV_scaled, 'k--')
xticks(axisTicks_tuning{1})
title('Mean(true template, 2) and scaled template for DV')
subplot(3,3,6), plot(axis_tuning{1}, margORI_true-margORI_forDV_scaled), xticks(axisTicks_tuning{1}), ylim([-.1, .1])
title('Difference')
% Row 3: SF tuning curve (averaged across SF)
subplot(3,3,8), hold on
plot(axis_tuning{2}, margSF_true, '-r'),
plot(axis_tuning{2}, margSF_forDV_scaled, 'k--')
xticks(axisTicks_tuning{2})
subplot(3,3,9), plot(axis_tuning{2}, margSF_true-margSF_forDV_scaled), xticks(axisTicks_tuning{2}), ylim([-.1, .1])

% Results: whitened template is slightly wider than the raw template

%% Project whitened energy and template onto the basis functions

    basisCfg_true = SX_RC_getBasisSettings(struct( ...
        'nBasisORI', nBasisORI, ...
        'nBasisSF', nBasisSF, ...
        'basisWidthScaleORI', basisWidthORI, ...
        'basisWidthScaleSF', basisWidthSF));

    basisOpts_true = struct();
    basisOpts_true.nBasisORI = basisCfg_true.nBasisORI;
    basisOpts_true.nBasisSF = basisCfg_true.nBasisSF;
    basisOpts_true.basisWidthScaleORI = basisCfg_true.basisWidthScaleORI;
    basisOpts_true.basisWidthScaleSF = basisCfg_true.basisWidthScaleSF;
    basisOpts_true.basisFamilyORI = basisCfg_true.basisFamilyORI;
    basisOpts_true.basisFamilySF = basisCfg_true.basisFamilySF;
    basisOpts_true.asymSF_rightLeftRatio = basisCfg_true.asymSF_rightLeftRatio;
    basisOpts_true.oriPeriod_deg = basisCfg_true.oriPeriod_deg;

    %------------------%
    Z_true_basis = SX_RC_basisProject(e3D_forDV, filtersOri_all - 90, noise.filtersSF_all_log, basisOpts_true);
    %------------------%
    %------------------%
    beta_true_basis = SX_RC_basisProject(reshape(template_forDV, [1, nORI, nSF]), filtersOri_all - 90, noise.filtersSF_all_log, basisOpts_true);
    %------------------%

    %% Compute DV from the basis representation
    DV_target_sim_allT = Z_true_basis * beta_true_basis(:);
    dv_space_used = 'zscore_absContrast_whiten_basis';

    % Store Q, mu and signal for both
mu_cov_DV = Tfix.mu_cov_abs;
sigma_cov_DV = Tfix.normStats.sigma_global_2D;
W_white_DV = Tfix.Q;
template_space_true = 'fixed_transformed';

assert(~any(isnan(DV_target_sim_allT)), 'DV_target contains NaN');

%% Add internal variability and generate responses
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

fprintf('%s: Internal noise sampled and added.\n\n', datetime('now'))

%% Set a true criterion given the true SDT criterion
% Candidate criteria in DV space:
DV_sorted = sort(unique(DVnoisy_sim_allT(:)));
cDV_grid = [-Inf; (DV_sorted(1:end-1) + DV_sorted(2:end))/2; Inf];
loss_grid = nan(size(cDV_grid));
cSDT_grid = nan(size(cDV_grid));
for iC = 1:numel(cDV_grid)
    c_try = cDV_grid(iC);
    [loss_grid(iC), cSDT_grid(iC)] = fxn_loss_cSDT(c_try, DVnoisy_sim_allT, iPRS_allT, cSDT_true);
end
% Choose the criterion in DV unit
[~, idx_best] = min(loss_grid);
criterion_DV_true = cDV_grid(idx_best);

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

fprintf('%s: Criterion calculated and responses generated.\n\n', datetime('now'))

%% Behavioral metrics
pYES_sim = mean(resp_allT == 1);
pHit_sim = sum((iPRS_allT == 1) & (resp_allT == 1)) / sum(iPRS_allT == 1);
pFA_sim = sum((iPRS_allT == 0) & (resp_allT == 1)) / sum(iPRS_allT == 0);
pC_sim = mean( (iPRS_allT==1 & resp_allT==1) | (iPRS_allT==0 & resp_allT==0) );

% Stop if simulated pC is too low or too high
if (pC_sim >= pC_filter(1)) && (pC_sim <= pC_filter(2))
    fprintf('\n\n ** Simulated pC=%.2f, withtin the range [%.2f, %.2f] ** \n\n', pC_sim, pC_filter)

    % Create output folders only after the accuracy gate is passed.
    if isempty(dir(nameFolder_Data_OOD_IO)), mkdir(nameFolder_Data_OOD_IO); end
    if isempty(dir(nameFolder_Data_NOM_IO)), mkdir(nameFolder_Data_NOM_IO); end

    % Save only variables used by simPlot4_VaryOneDim.m in one atomic write.
    save(sprintf('%s/truth.mat', nameFolder_Data_OOD_IO), ...
        'template_true', 'template_forDV', 'template_space_true', 'dv_space_used', ...
        'mu_cov_DV', 'sigma_cov_DV', 'W_white_DV', 'lambda_whiten', ...
        'beta_true_basis', 'basisOpts_true', ...
        'criterion_DV_true', 'N*_true');
    fprintf('\n\n%s: truth.mat saved (template_true, criterion_DV_true).\n\n', datetime('now'))

    % Save energy for compDV after the accuracy gate.
    % Use patchMode = ''T'' in compDV; this file must contain both target + noise
    % energy as well as noise, filters, stim, nBins.
    save(sprintf('%s/energy_T_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF), 'e3D_target_allT', 'noise', 'filtersOri_all', 'stim', 'nBins');
    fprintf('%s: Simulated 3D energy saved (%d trials).\n\n', datetime('now'), nTrials)

    [dprime_sim, cSDT_sim] = SX_sim06_SDT(pHit_sim, pFA_sim);

    respC_sim = nan(nPairs, 1);
    for iPairUnik = 1:nPairs
        respAB = resp_allT(iPair_allT == iPairUnik);
        respC_sim(iPairUnik) = (respAB(1) == respAB(2));
    end
    pA_sim = mean(respC_sim);

    metrics_sim = [dprime_sim, cSDT_sim, pC_sim, pHit_sim, pFA_sim, pA_sim, pYES_sim];

    % Save behavioral measures in behavMeas.mat (as expected by compDV)
    save(sprintf('%s/behavMeas.mat', nameFolder_Data_OOD_IO), 'dataMatrix', 'metrics_sim', 'iPRS_allT', 'iPass_allT', 'iPair_allT', 'resp_allT', 'pC_filter');

    fprintf('\n%s: Behavioral data saved.\n\nReady for template generation\n\n', datetime('now'))

    % Create figure folder only when plotting is enabled (never on HPC).
    if flag_plotDist && isempty(dir(nameFolder_Figures_perSubj))
        mkdir(nameFolder_Figures_perSubj);
    end

    %% Plotting DV dist and behav metrics
    if flag_plotDist
        NOMplot_dist
    end

    %% Run compDV and fitNOM on this IO
    % Step 1: compute DVs, templates, and test-set metrics
    %----------------------------%
    OOD_NOM_Trialwise_compDV_A12({nameIO, criterion_DV_true}, iLocComb, lambda_whiten, nIter, nJob, iJob)
    %----------------------------%

    %% Step 2: fit NOM parameters and predict metrics
    for iModelA_fit = iModelA_fit_all

        for iModelB_fit = iModelB_fit_all
            %----------------------------%
            OOD_NOM_Trialwise_fitNOM({nameIO, criterion_DV_true}, iLocComb, iModelA_fit, iModelB_fit, nIter, nJob, iJob)
            %----------------------------%

            if iModelB_fit == iModelB_sim
                % For the simulated model, also plot DV vs DV scatter and parameter recovery

                if flag_plotDist
                    % plot_CorrBasisSetting(nameFolder_Data_NOM_IO, nameFolder_Figures_perSubj, nameIO, ...
                    %     nIter, iJob, iLocComb, iModelA_fit, iModelB_fit_all, namesModelBparams_short, ...
                    %     Nmul_true, Nadd_true, Nshared_true, criterion_DV_true, template_true);
                end
            end % if iModelB_fit == iModelB_sim
        end % for iModelB_fit

    end % for iModelA_fit
    %%
    % Summarize model comparison after all A/B fits finish.
    if flag_plotDist
        % plot_fit_model_comparison( ...
        %     nameFolder_Data_NOM_IO, nameFolder_Figures_perSubj, nameIO, ...
        %     iModelA_fit_all, iModelB_fit_all, namesModelA, namesModelB, namesModelBparams_short, ...
        %     nIter, iJob, iLocComb, iModelB_sim, ...
        %     Nmul_true, Nadd_true, Nshared_true, criterion_DV_true);
    end
else
    fprintf('\n\n ** Simulated pC=%.2f, OUT OF the range [%.2f, %.2f] ** \n\n', pC_sim, pC_filter)
end % if

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

% end % end of the OOD_sim function

%% helper
% function plot_fit_model_comparison(nameFolder_Data_NOM_IO, nameFolder_Figures_perSubj, nameIO, ...
%     iModelA_fit_all, iModelB_fit_all, namesModelA, namesModelB, namesModelBparams_short, ...
%     nIter, iJob, iLocComb, iModelB_sim, ...
%     Nmul_true, Nadd_true, Nshared_true, criterion_DV_true)
% % This helpder function summarizes the model comparison results across all fitted A/B models, and plots the summary.
% metricNames = {'pYES', 'pC', 'pA'};
% metricLabels = {'pYES', 'pC', 'pA'};
% paramNames = {'Nmul', 'Nadd', 'Nshared', 'criterion_DV'};
% paramLabels = {'Nmul', 'Nadd', 'Nshared', 'criterion DV'};

% nA = numel(iModelA_fit_all);
% nB = numel(iModelB_fit_all);

% summaryA = repmat(struct( ...
%     'winRate', nan(1, nB), ...
%     'deltaNLL', nan(1, nB), ...
%     'metricRMSE', nan(nB, numel(metricNames)), ...
%     'paramRMSE', nan(nB, numel(paramNames)), ...
%     'nValidIter', 0), nA, 1);

% trueParamVals = [Nmul_true, Nadd_true, Nshared_true, criterion_DV_true];

% for iA = 1:nA
%     iModelA_fit = iModelA_fit_all(iA);
%     nLL_test_byModel = nan(nIter, nB);

%     for iB = 1:nB
%         iModelB_fit = iModelB_fit_all(iB);
%         nameFile_fitNOM = sprintf('%s/n%d_J%d_A%dB%d.mat', ...
%             nameFolder_Data_NOM_IO, nIter, iJob, iModelA_fit, iModelB_fit);

%         if ~exist(nameFile_fitNOM, 'file')
%             continue
%         end

%         S_fit = load(nameFile_fitNOM, 'nLL_test_allIter', 'pred_metrics_allIter', 'params_est_allIter');

%         if isfield(S_fit, 'nLL_test_allIter') && ~isempty(S_fit.nLL_test_allIter)
%             nCopy = min(nIter, numel(S_fit.nLL_test_allIter));
%             nLL_test_byModel(1:nCopy, iB) = S_fit.nLL_test_allIter(1:nCopy);
%         end

%         if isfield(S_fit, 'pred_metrics_allIter') && ~isempty(S_fit.pred_metrics_allIter)
%             for iMetric = 1:numel(metricNames)
%                 rmse_iter = metric_recovery_iter(S_fit.pred_metrics_allIter, metricNames{iMetric});
%                 summaryA(iA).metricRMSE(iB, iMetric) = median(rmse_iter(isfinite(rmse_iter)), 'omitnan');
%             end
%         end

%         if isfield(S_fit, 'params_est_allIter') && ~isempty(S_fit.params_est_allIter)
%             paramRMSE_thisModel = compute_param_rmse_by_model( ...
%                 S_fit.params_est_allIter, namesModelBparams_short{iModelB_fit}, paramNames, trueParamVals);
%             summaryA(iA).paramRMSE(iB, :) = paramRMSE_thisModel;
%         end
%     end

%     [summaryA(iA).winRate, summaryA(iA).deltaNLL, summaryA(iA).nValidIter] = summarize_nll_compare(nLL_test_byModel);
% end

% figPos = [100, 100, 1300, 1100];
% figure('Position', figPos);
% tiledlayout(4, nA, 'TileSpacing', 'compact', 'Padding', 'compact');

% metricColors = [0.25 0.45 0.75; 0.20 0.65 0.35; 0.80 0.35 0.20];
% paramColors = [0.20 0.45 0.75; 0.20 0.70 0.55; 0.80 0.55 0.20; 0.55 0.35 0.75];
% xTicks = 1:nB;
% xLabels = cellfun(@(x) char(string(x)), namesModelB(iModelB_fit_all), 'UniformOutput', false);

% for iA = 1:nA
%     S = summaryA(iA);

%     nexttile(iA); hold on;
%     bar(xTicks, S.winRate, 0.75, 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
%     ylim([0 1]);
%     xlim([0.4 nB + 0.6]);
%     set(gca, 'XTick', xTicks, 'XTickLabel', xLabels, 'XTickLabelRotation', 0);
%     ylabel('Win rate');
%     title(sprintf('A%d %s', iModelA_fit_all(iA), namesModelA{iModelA_fit_all(iA)}), 'Interpreter', 'none');
%     box on;

%     nexttile(nA + iA); hold on;
%     bar(xTicks, S.deltaNLL, 0.75, 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
%     xlim([0.4 nB + 0.6]);
%     set(gca, 'XTick', xTicks, 'XTickLabel', xLabels, 'XTickLabelRotation', 0);
%     ylabel('\Delta nLL from best');
%     box on;

%     nexttile(2*nA + iA); hold on;
%     hb = bar(S.metricRMSE, 'grouped');
%     for iMetric = 1:numel(hb)
%         hb(iMetric).FaceColor = metricColors(iMetric, :);
%         hb(iMetric).EdgeColor = 'none';
%     end
%     xlim([0.4 nB + 0.6]);
%     set(gca, 'XTick', xTicks, 'XTickLabel', xLabels, 'XTickLabelRotation', 0);
%     ylabel('Metric RMSE');
%     box on;
%     if iA == nA
%         legend(metricLabels, 'Location', 'best', 'Box', 'off');
%     end

%     nexttile(3*nA + iA); hold on;
%     hb = bar(S.paramRMSE, 'grouped');
%     for iParam = 1:numel(hb)
%         hb(iParam).FaceColor = paramColors(iParam, :);
%         hb(iParam).EdgeColor = 'none';
%     end
%     xlim([0.4 nB + 0.6]);
%     set(gca, 'XTick', xTicks, 'XTickLabel', xLabels, 'XTickLabelRotation', 0);
%     ylabel('Parameter RMSE');
%     xlabel('Fitting model');
%     box on;
%     if iA == nA
%         legend(paramLabels, 'Location', 'best', 'Box', 'off');
%     end
% end

% sgtitle(sprintf(['Model comparison across fitted B models\n' ...
%     '%s | Bsim=%d | L%d | nIter=%d'], ...
%     nameIO, iModelB_sim, iLocComb, nIter), 'Interpreter', 'none');

% saveas(gcf, sprintf('%s/3Summary.jpg', nameFolder_Figures_perSubj));
% close(gcf);
% end

%% helper
function [winRate, deltaNLL, nValidIter] = summarize_nll_compare(nLL_test_byModel)
% This function summarizes the model comparison based on test-set nLL across all iterations, and returns the win rate and median delta nLL for each model.
nB = size(nLL_test_byModel, 2);
winCount = zeros(1, nB);
deltaVals = cell(1, nB);
nValidIter = 0;

for iIter = 1:size(nLL_test_byModel, 1)
    thisNLL = nLL_test_byModel(iIter, :);
    finiteMask = isfinite(thisNLL);
    if ~any(finiteMask)
        continue
    end

    bestNLL = min(thisNLL(finiteMask));
    idxFinite = find(finiteMask);
    idxBestLocal = find(thisNLL(finiteMask) == bestNLL, 1, 'first');
    idxBest = idxFinite(idxBestLocal);

    winCount(idxBest) = winCount(idxBest) + 1;
    nValidIter = nValidIter + 1;

    for iB = idxFinite
        deltaVals{iB}(end+1, 1) = thisNLL(iB) - bestNLL; %#ok<AGROW>
    end
end

winRate = nan(1, nB);
deltaNLL = nan(1, nB);

if nValidIter > 0
    winRate = winCount / nValidIter;
end

for iB = 1:nB
    if ~isempty(deltaVals{iB})
        deltaNLL(iB) = median(deltaVals{iB}, 'omitnan');
    end
end
end

%% helper
function rmse_iter = metric_recovery_iter(pred_metrics_allIter, metricName)
% This function computes the RMSE of a specific metric (e.g., pYES) for each iteration, based on the predicted vs data metric values across all bins.
nIter = numel(pred_metrics_allIter);
rmse_iter = nan(nIter, 1);

for iIter = 1:nIter
    yData = pred_metrics_allIter{iIter}.metrics.([metricName '_data_allBins'])(:);
    yPred = pred_metrics_allIter{iIter}.metrics.([metricName '_pred_allBins'])(:);
    nTrials = pred_metrics_allIter{iIter}.metrics.nTrials_allBins(:);

    good = isfinite(yData) & isfinite(yPred) & isfinite(nTrials) & (nTrials > 0);
    yData = yData(good);
    yPred = yPred(good);
    w = nTrials(good);

    if isempty(yData)
        continue
    end

    rmse_iter(iIter) = sqrt(sum(w .* (yData - yPred).^2) / sum(w));
end
end

%% helper
function paramRMSE = compute_param_rmse_by_model(params_est_allIter, paramNamesThisModel, paramNamesAll, trueParamVals)
% This function computes the RMSE of each parameter for a specific model, based on the estimated vs true parameter values across all iterations. The input params_est_allIter can be either a cell array (nIter x 1) or a numeric array (nIter x nParamsThisModel).
paramRMSE = nan(1, numel(paramNamesAll));

if iscell(params_est_allIter)
    P = cell2mat(cellfun(@(x) x(:)', params_est_allIter, 'UniformOutput', false));
else
    P = params_est_allIter;
end

for iParam = 1:numel(paramNamesThisModel)
    idxParamAll = find(strcmp(paramNamesAll, paramNamesThisModel{iParam}), 1, 'first');
    if isempty(idxParamAll)
        continue
    end

    est_iter = P(:, iParam);
    err_iter = abs(est_iter - trueParamVals(idxParamAll));
    paramRMSE(idxParamAll) = median(err_iter(isfinite(err_iter)), 'omitnan');
end
end

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

% %% helper
% %% helper
% function [mu_cov, sigma_cov_2D, e3D_centered] = recenter_e3D(e3D_in, eps_sigma)
% % This function recenters a 3D array e3D_in by subtracting the mean across trials.
% % It also computes the standard deviation for each element, ensuring a minimum value of eps_sigma.
% [nTrials, nORI, nSF] = size(e3D_in);
% X = reshape(e3D_in, [nTrials, nORI * nSF]);
% mu_cov = mean(X, 1);
% X_centered = X - mu_cov;

% sigma_cov = std(X, [], 1);
% sigma_cov(~isfinite(sigma_cov) | sigma_cov < eps_sigma) = eps_sigma;
% sigma_cov_2D = reshape(sigma_cov, [nORI, nSF]);

% e3D_centered = reshape(X_centered, [nTrials, nORI, nSF]);
% end

%% helper
% function W_white = build_whitening_matrix(e3D_centered, lambda_whiten, eps_whiten)
% % This function builds a whitening matrix W_white based on the covariance of the centered 3D data.
% % It applies shrinkage to the covariance matrix to ensure numerical stability, and then computes the whitening transformation.
% [nTrials, nORI, nSF] = size(e3D_centered);
% X_centered = reshape(e3D_centered, [nTrials, nORI * nSF]);

% Sigma = cov(X_centered, 1);
% Sigma_shrink = (1 - lambda_whiten) * Sigma + lambda_whiten * mean(diag(Sigma)) * eye(size(Sigma, 1));
% Sigma_sym = (Sigma_shrink + Sigma_shrink') / 2;

% [V, D] = eig(Sigma_sym);
% d = diag(D);
% d(~isfinite(d) | d < eps_whiten) = eps_whiten;
% W_white = V * diag(1 ./ sqrt(d)) * V';
% end

% %% helper
% function e3D_out = apply_linear_transform_e3D(e3D_centered, W)
% % This function applies a linear transformation W to the centered 3D data e3D_centered.
% % The transformation is applied to the trial dimension, and the output is reshaped back to the original 3D format.
% [nTrials, nORI, nSF] = size(e3D_centered);
% X = reshape(e3D_centered, [nTrials, nORI * nSF]);
% X_out = X * W;
% e3D_out = reshape(X_out, [nTrials, nORI, nSF]);
% end

% %% helper
% function template_out = convert_template_space(template_in, fromSpace, toSpace, sigma_cov_2D, W_white, eps_sigma)
% % This function converts a template from one space to another.
% % Supported spaces: 'raw_centered', 'z', 'white'.
% template_vec = template_in(:);

% if strcmp(fromSpace, toSpace)
%     template_out = template_in;
%     return;
% end

% switch fromSpace
%     case 'raw_centered'
%         template_raw = template_vec;
%     case 'z'
%         sigma_vec = sigma_cov_2D(:);
%         sigma_vec(~isfinite(sigma_vec) | sigma_vec < eps_sigma) = eps_sigma;
%         template_raw = template_vec ./ sigma_vec;
%     case 'white'
%         if isempty(W_white)
%             error('W_white is required to convert from white to raw_centered.');
%         end
%         template_raw = W_white * template_vec;
%     otherwise
%         error('Unknown template space: %s', fromSpace);
% end

% switch toSpace
%     case 'raw_centered'
%         template_vec_out = template_raw;
%     case 'z'
%         sigma_vec = sigma_cov_2D(:);
%         sigma_vec(~isfinite(sigma_vec) | sigma_vec < eps_sigma) = eps_sigma;
%         template_vec_out = template_raw .* sigma_vec;
%     case 'white'
%         if isempty(W_white)
%             error('W_white is required to convert from raw_centered to white.');
%         end
%         template_vec_out = W_white \ template_raw;
%     otherwise
%         error('Unknown target template space: %s', toSpace);
% end

% template_out = reshape(template_vec_out, size(template_in));
% end
