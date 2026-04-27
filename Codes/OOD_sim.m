% function OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, Cz_true, ...
%     lambda_whiten, flag_regressType, flag_incluCrit, C_contribution, iModelB_sim, nIter, nBasisORI, nBasisSF)

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

% Simulation input parameters
noiseCST = .2;
gaborCST = .5;
nTrials = 4e3;
Nmul_true = .5;
Nadd_true = 5;
Nshared_true = 5;
Cz_true = 0;
lambda_whiten = 0;
flag_regressType = 2;
flag_incluCrit = 1;
C_contribution = .5;
iModelB_sim = 1;
nIter = 5;
nBasisORI = 6;
nBasisSF = 5;

% Basic setup -%
clc; close all;
warning off;
set(0, 'DefaultFigureVisible', 'off')

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
%--%
SX_RC1_setting;
%--%

nORI = nORI; % nORI is defined in SX_RC1_setting; repeated here just for parfor loop to work
nSF = nORI;

% Number of bootstraps for compIV / fitNOM
nJob = 1;
iJob = 1;

% Model A/B indices for fitting
iModelA_fit_all = [1,2]; % DO NOT CHANGE! 1 = RC-derived template (Model A), 2=ideal template; 3=permuted template
iModelB_fit_all = 1:4;
pC_filter = [.6, .8];

templateType_true = 1; % 1 = raw;
IVType_true = 1; % 1 = sum of dot product;
convolveType_true = 1; % 1 = dot product; 2 = convolution (for fxn_getIV_v3)
flag_permT = 0; %1=permute the input template per trial
flag_plotDist = 1;
if strcmp(str_envir, 'HPC'), flag_plotDist = 0; end % don't plot when running on HPC
iLocComb = 1; % use single-location index (e.g., fovea) for IO

% Print info --%

fprintf(' - nIter (for compIV/fitNOM) = %d\n', nIter);
fprintf(' - noiseCST = %.2f\n', noiseCST);
fprintf(' - gaborCST = %.2f\n', gaborCST);
fprintf(' - nTrials = %d\n', nTrials);
fprintf(' - Simulated with ModelB = %d \n', iModelB_sim);
fprintf(' - Fitted with ModelA = %s (1=Data-derived template; 2=ideal template)\n', strjoin(string(iModelA_fit_all), ' '));
fprintf(' - Fitted with ModelB = %s (1=full, 2=No Nshared; 3=No Nmul; 4=No Nadd)\n', strjoin(string(iModelB_fit_all), ' '));
fprintf(' - Criterion (z unit)=%.1f \n', Cz_true);
fprintf(' - Contribution of criterion loss =%.1f \n', C_contribution);
fprintf(' - Whitening strength (lambda): %.1f \n', lambda_whiten);
fprintf(' - Regression type (1=Univariate; 2=Multi+smoothing): %d \n\n', flag_regressType);

% Define IO name & folders %
% Define the IO name
nameIO = sprintf('IO_cN%.0f_cG%.0f_nT%s_Nm%s_Na%s_Ns%s_Cz%.1f_cont%.1f_whiten%.1f_R%d_%d%d_B%d', ...
    noiseCST*100, gaborCST*100, format_num2exp(nTrials), ...
    format_num2exp(Nmul_true), format_num2exp(Nadd_true), format_num2exp(Nshared_true), ...
    Cz_true, C_contribution, lambda_whiten, flag_regressType, nORI, nSF, iModelB_sim);

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
%-%
template_true = SX_RC4_Energy_parfor(stim.mask, {template_gabor_true}, filter_sin, filter_cos);
%-%
template_true = squeeze(template_true);
%-%
template_true = fxn_getTemplate(template_true, templateType_true, 0);
%-%

% Normalize template (Unit L2-norm)
template_true = template_true / norm(template_true(:));

% Save "truth" (for IO template in Model A = 2)
save(sprintf('%s/truth.mat', nameFolder_Data_OOD_IO), ...
    'nameIO', '*_true', 'noiseCST', 'gaborCST', 'nTrials', 'lambda_whiten', 'flag_regressType', 'flag_incluCrit', 'C_contribution', 'iModelB_sim', 'iModelA_fit_all', 'iModelB_fit_all', 'nIter');
fprintf('\n\n%s: Truth (including the ideal template) created and saved.\n\n', datetime('now'))

%% Preallocate sim arrays
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

%% Simulate stimuli & IV
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

    % Decision variable from energy × template_true
    DV_target = fxn_getIV_v3(e3D_target, template_true, convolveType_true, IVType_true, flag_permT, [1, nORI]);
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

%% Save energy for compIV
% Use patchMode = 'T' in compIV; this file must contain both target + noise
% energy as well as noise, filters, stim, nBins.
save(sprintf('%s/energy_T_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF), 'e3D_target_allT', 'noise', 'filtersOri_all', 'stim', 'nBins');

fprintf('%s: Simulated 3D energy saved (%d trials).\n\n', datetime('now'), nTrials)

%% Internal noise & responses
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

%% Set a true criterion given the true criterion in z-score
% Old way: Titrate a criterion to reach target accuracy pC_titrate
%-%
% fxn_loss_pC_ = @(criterion_potential) fxn_loss_pC(criterion_potential, pC_titrate, DVnoisy_sim_allT, iPRS_allT);
% %-%
% criterion_true = fmincon(fxn_loss_pC_, median(DVnoisy_sim_allT), [], [], [], [], min(DVnoisy_sim_allT), max(DVnoisy_sim_allT));

% Candidate criteria in DV space:
DV_sorted = sort(unique(DVnoisy_sim_allT(:)));
criterion_grid = [-Inf; (DV_sorted(1:end-1) + DV_sorted(2:end))/2; Inf];

loss_grid = nan(size(criterion_grid));
c_z_grid = nan(size(criterion_grid));

for iC = 1:numel(criterion_grid)
    c_try = criterion_grid(iC);
    [loss_grid(iC), c_z_grid(iC)] = criterion_loss_cz(c_try, DVnoisy_sim_allT, iPRS_allT, Cz_true);
end
% Choose the criterion in DV unit
[~, idx_best] = min(loss_grid);
criterion_DV_true = criterion_grid(idx_best);

save(sprintf('%s/truth.mat', nameFolder_Data_OOD_IO), 'criterion_DV_true', '-append');

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
    [dprime_sim, Cz_sim] = SX_sim06_SDT(pHit_sim, pFA_sim);

    respC_sim = nan(nPairs, 1);
    for iPairUnik = 1:nPairs
        respAB = resp_allT(iPair_allT == iPairUnik);
        respC_sim(iPairUnik) = (respAB(1) == respAB(2));
    end
    pA_sim = mean(respC_sim);

    metrics_sim = [dprime_sim, Cz_sim, pC_sim, pHit_sim, pFA_sim, pA_sim, pYES_sim];

    % Save behavioral measures in behavMeas.mat (as expected by compIV)
    save(sprintf('%s/behavMeas.mat', nameFolder_Data_OOD_IO), 'dataMatrix', 'metrics_sim', 'iPRS_allT', 'iPass_allT', 'iPair_allT', 'resp_allT', 'pC_filter');

    fprintf('\n%s: Behavioral data saved.\n\nReady for template generation\n\n', datetime('now'))

    %% Plotting DV dist and behav metrics
    if flag_plotDist
        NOMplot_dist
    end

    %% Run compIV and fitNOM on this IO
    % Step 1: compute IVs, templates, and test-set metrics
    %----------------------------%
    OOD_NOM_Trialwise_compIV_A12(nBasisORI, nBasisSF, {nameIO, criterion_DV_true}, iLocComb, lambda_whiten, flag_regressType, nIter, nJob, iJob)
    %----------------------------%

    %% Step 2: fit NOM parameters and predict metrics
    for iModelA_fit=1:2
        
        for iModelB_fit = iModelB_fit_all
            %----------------------------%
            OOD_NOM_Trialwise_fitNOM(nBasisORI, nBasisSF, {nameIO, criterion_DV_true}, iLocComb, flag_incluCrit, C_contribution, iModelA_fit, iModelB_fit, nIter, nJob, iJob)
            %----------------------------%

            if iModelB_fit == iModelB_sim
                % For the simulated model, also plot IV vs DV scatter and parameter recovery
                
                plot_CorrBasisSetting(nameFolder_Data_NOM_IO, nameFolder_Figures_perSubj, nameIO, ...
                    nIter, iJob, iLocComb, iModelA_fit, iModelB_fit_all, namesModelBparams_short, ...
                    Nmul_true, Nadd_true, Nshared_true, criterion_DV_true, template_true);
                
            end % if iModelB_fit == iModelB_sim
        end % for iModelB_fit
        
    end % for iModelA_fit
%%
    % Summarize model comparison after all A/B fits finish.
    plot_fit_model_comparison( ...
        nameFolder_Data_NOM_IO, nameFolder_Figures_perSubj, nameIO, ...
        iModelA_fit_all, iModelB_fit_all, namesModelA, namesModelB, namesModelBparams_short, ...
        nIter, iJob, iLocComb, iModelB_sim, ...
        Nmul_true, Nadd_true, Nshared_true, criterion_DV_true);

else
    fprintf('\n\n ** Simulated pC=%.2f, outof the range [%.2f, %.2f] ** \n\n', pC_sim, pC_filter)
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

% Remove the energy file to save space (compIV has already loaded it)
delete(sprintf('%s/energy_T_%d_%d.mat', nameFolder_Data_OOD_IO, nORI, nSF));
fprintf('%s: Energy file removed to save space.\n\n', datetime('now'))

%%  End timing 
time_end = datetime('now');
fprintf('%s: Simulation done.\n\n', time_end)
elapsed = time_end - time_start;
fprintf('Time used: %s\n\n\n\n', char(elapsed));

% end

%% helper
function plot_fit_model_comparison(nameFolder_Data_NOM_IO, nameFolder_Figures_perSubj, nameIO, ...
    iModelA_fit_all, iModelB_fit_all, namesModelA, namesModelB, namesModelBparams_short, ...
    nIter, iJob, iLocComb, iModelB_sim, ...
    Nmul_true, Nadd_true, Nshared_true, criterion_DV_true)

metricNames = {'pYES', 'pC', 'pA'};
metricLabels = {'pYES', 'pC', 'pA'};
paramNames = {'Nmul', 'Nadd', 'Nshared', 'criterion_DV'};
paramLabels = {'Nmul', 'Nadd', 'Nshared', 'criterion DV'};

nA = numel(iModelA_fit_all);
nB = numel(iModelB_fit_all);

summaryA = repmat(struct( ...
    'winRate', nan(1, nB), ...
    'deltaNLL', nan(1, nB), ...
    'metricRMSE', nan(nB, numel(metricNames)), ...
    'paramRMSE', nan(nB, numel(paramNames)), ...
    'nValidIter', 0), nA, 1);

trueParamVals = [Nmul_true, Nadd_true, Nshared_true, criterion_DV_true];

for iA = 1:nA
    iModelA_fit = iModelA_fit_all(iA);
    nLL_test_byModel = nan(nIter, nB);

    for iB = 1:nB
        iModelB_fit = iModelB_fit_all(iB);
        nameFile_fitNOM = sprintf('%s/n%d_J%d_A%dB%d.mat', ...
            nameFolder_Data_NOM_IO, nIter, iJob, iModelA_fit, iModelB_fit);

        if ~exist(nameFile_fitNOM, 'file')
            continue
        end

        S_fit = load(nameFile_fitNOM, 'nLL_test_allIter', 'pred_metrics_allIter', 'params_est_allIter');

        if isfield(S_fit, 'nLL_test_allIter') && ~isempty(S_fit.nLL_test_allIter)
            nCopy = min(nIter, numel(S_fit.nLL_test_allIter));
            nLL_test_byModel(1:nCopy, iB) = S_fit.nLL_test_allIter(1:nCopy);
        end

        if isfield(S_fit, 'pred_metrics_allIter') && ~isempty(S_fit.pred_metrics_allIter)
            for iMetric = 1:numel(metricNames)
                rmse_iter = metric_recovery_iter(S_fit.pred_metrics_allIter, metricNames{iMetric});
                summaryA(iA).metricRMSE(iB, iMetric) = median(rmse_iter(isfinite(rmse_iter)), 'omitnan');
            end
        end

        if isfield(S_fit, 'params_est_allIter') && ~isempty(S_fit.params_est_allIter)
            paramRMSE_thisModel = compute_param_rmse_by_model( ...
                S_fit.params_est_allIter, namesModelBparams_short{iModelB_fit}, paramNames, trueParamVals);
            summaryA(iA).paramRMSE(iB, :) = paramRMSE_thisModel;
        end
    end

    [summaryA(iA).winRate, summaryA(iA).deltaNLL, summaryA(iA).nValidIter] = summarize_nll_compare(nLL_test_byModel);
end

figPos = [100, 100, 1300, 1100];
figure('Position', figPos);
tiledlayout(4, nA, 'TileSpacing', 'compact', 'Padding', 'compact');

metricColors = [0.25 0.45 0.75; 0.20 0.65 0.35; 0.80 0.35 0.20];
paramColors = [0.20 0.45 0.75; 0.20 0.70 0.55; 0.80 0.55 0.20; 0.55 0.35 0.75];
xTicks = 1:nB;
xLabels = cellfun(@(x) char(string(x)), namesModelB(iModelB_fit_all), 'UniformOutput', false);

for iA = 1:nA
    S = summaryA(iA);

    nexttile(iA); hold on;
    bar(xTicks, S.winRate, 0.75, 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
    ylim([0 1]);
    xlim([0.4 nB + 0.6]);
    set(gca, 'XTick', xTicks, 'XTickLabel', xLabels, 'XTickLabelRotation', 0);
    ylabel('Win rate');
    title(sprintf('A%d %s', iModelA_fit_all(iA), namesModelA{iModelA_fit_all(iA)}), 'Interpreter', 'none');
    box on;

    nexttile(nA + iA); hold on;
    bar(xTicks, S.deltaNLL, 0.75, 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
    xlim([0.4 nB + 0.6]);
    set(gca, 'XTick', xTicks, 'XTickLabel', xLabels, 'XTickLabelRotation', 0);
    ylabel('\Delta nLL from best');
    box on;

    nexttile(2*nA + iA); hold on;
    hb = bar(S.metricRMSE, 'grouped');
    for iMetric = 1:numel(hb)
        hb(iMetric).FaceColor = metricColors(iMetric, :);
        hb(iMetric).EdgeColor = 'none';
    end
    xlim([0.4 nB + 0.6]);
    set(gca, 'XTick', xTicks, 'XTickLabel', xLabels, 'XTickLabelRotation', 0);
    ylabel('Metric RMSE');
    box on;
    if iA == nA
        legend(metricLabels, 'Location', 'best', 'Box', 'off');
    end

    nexttile(3*nA + iA); hold on;
    hb = bar(S.paramRMSE, 'grouped');
    for iParam = 1:numel(hb)
        hb(iParam).FaceColor = paramColors(iParam, :);
        hb(iParam).EdgeColor = 'none';
    end
    xlim([0.4 nB + 0.6]);
    set(gca, 'XTick', xTicks, 'XTickLabel', xLabels, 'XTickLabelRotation', 0);
    ylabel('Parameter RMSE');
    xlabel('Fitting model');
    box on;
    if iA == nA
        legend(paramLabels, 'Location', 'best', 'Box', 'off');
    end
end

sgtitle(sprintf(['Model comparison across fitted B models\n' ...
    '%s | Bsim=%d | L%d | nIter=%d'], ...
    nameIO, iModelB_sim, iLocComb, nIter), 'Interpreter', 'none');

saveas(gcf, sprintf('%s/3Summary.jpg', nameFolder_Figures_perSubj));
close(gcf);
end

%% helper
function [winRate, deltaNLL, nValidIter] = summarize_nll_compare(nLL_test_byModel)

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

%% helper
function plot_CorrBasisSetting(nameFolder_Data_NOM_IO, nameFolder_Figures_perSubj, nameIO, ...
    nIter, iJob, iLocComb, iModelA_fit, iModelB_fit_all, namesModelBparams_short, ...
    Nmul_true, Nadd_true, Nshared_true, criterion_DV_true, template_true)

sz_font = 15;
nameFile_compIV = sprintf('%s/n%d_J%d_A%d_compIV.mat', nameFolder_Data_NOM_IO, nIter, iJob, iModelA_fit);
if ~exist(nameFile_compIV, 'file')
    nameFile_compIV = sprintf('%s/n%d_J%d_A%d_compIV', nameFolder_Data_NOM_IO, nIter, iJob, iModelA_fit);
end
if ~exist(nameFile_compIV, 'file')
    fprintf('WARNING: compIV file for A%d not found. Skip IVxDV scatter.\n', iModelA_fit);
    return
end

S = load(nameFile_compIV, ...
    'template_full_allIter', 'margORI_allIter', 'margSF_allIter', ...
    'margParams_ORI_allIter', 'margParams_SF_allIter', ...
    'nBasisORI_tmpl_allIter', 'nBasisSF_tmpl_allIter', ...
    'basisFxnORI_tmpl_allIter', 'basisFxnSF_tmpl_allIter', 'ridge_tmpl_allIter', ...
    'subjName', 'namesModelA');

if ~isfield(S, 'template_full_allIter') || isempty(S.template_full_allIter)
    fprintf('WARNING: template_full_allIter missing in A%d compIV file. Skip IVxDV scatter.\n', iModelA_fit);
    return
end

% IV list (x-axis in each row)
IV_vals = {
    S.nBasisORI_tmpl_allIter, ...
    S.nBasisSF_tmpl_allIter, ...
    S.basisFxnORI_tmpl_allIter, ...
    S.basisFxnSF_tmpl_allIter, ...
    S.ridge_tmpl_allIter ...
    };
IV_names = {'nBasisORI', 'nBasisSF', 'basisFamilyORI', 'basisFamilySF', 'ridge'};

% Reference template and reference marginals
template_ref = template_true;
template_ref = template_ref / max(norm(template_ref(:)), eps);
marg_ref_ORI = squeeze(mean(template_ref, 2))';
marg_ref_SF = squeeze(mean(template_ref, 1));
marg_ref_ORI = normalize_by_maxabs_local(marg_ref_ORI);
marg_ref_SF = normalize_by_maxabs_local(marg_ref_SF);

% A2 may not include marginalized profiles/params in some pipelines.
% If missing, construct profiles from template and keep param-DVs as NaN.
if ~isfield(S, 'margORI_allIter') || isempty(S.margORI_allIter)
    S.margORI_allIter = nan(nIter, 2, size(S.template_full_allIter, 2));
    for iIter = 1:min(nIter, size(S.template_full_allIter, 1))
        tmpl = squeeze(S.template_full_allIter(iIter, :, :));
        S.margORI_allIter(iIter, 2, :) = mean(tmpl, 2);
    end
end
if ~isfield(S, 'margSF_allIter') || isempty(S.margSF_allIter)
    S.margSF_allIter = nan(nIter, 2, size(S.template_full_allIter, 3));
    for iIter = 1:min(nIter, size(S.template_full_allIter, 1))
        tmpl = squeeze(S.template_full_allIter(iIter, :, :));
        S.margSF_allIter(iIter, 2, :) = mean(tmpl, 1);
    end
end
if ~isfield(S, 'margParams_ORI_allIter') || isempty(S.margParams_ORI_allIter)
    S.margParams_ORI_allIter = nan(nIter, 2, 0);
end
if ~isfield(S, 'margParams_SF_allIter') || isempty(S.margParams_SF_allIter)
    S.margParams_SF_allIter = nan(nIter, 2, 0);
end

if isfield(S, 'subjName')
    subjTag = S.subjName;
else
    subjTag = nameIO;
end
if isfield(S, 'namesModelA') && numel(S.namesModelA) >= iModelA_fit
    modelA_name = S.namesModelA{iModelA_fit};
else
    modelA_name = sprintf('A%d', iModelA_fit);
end

% Figure 1 per A: core RMSE (3 columns)
DV_core = nan(nIter, 3);
for iIter = 1:nIter
    template_i = squeeze(S.template_full_allIter(iIter, :, :));
    template_i = template_i / max(norm(template_i(:)), eps);
    DV_core(iIter, 1) = rmse_simple_local(template_i(:), template_ref(:));

    margORI_i = squeeze(S.margORI_allIter(iIter, 2, :))';
    margORI_i = normalize_by_maxabs_local(margORI_i);
    DV_core(iIter, 2) = rmse_simple_local(margORI_i(:), marg_ref_ORI(:));

    margSF_i = squeeze(S.margSF_allIter(iIter, 2, :))';
    margSF_i = normalize_by_maxabs_local(margSF_i);
    DV_core(iIter, 3) = rmse_simple_local(margSF_i(:), marg_ref_SF(:));
end

dvNames_core = {'Template RMSE', 'ORI tuning RMSE', 'SF tuning RMSE'};
figure('Position', [20 20 1600 1600]);
rng(1);
for iIV = 1:5
    [x_raw, x_tick, x_ticklabel, isCategorical] = encode_iv_for_scatter_local(IV_vals{iIV});
    for iDV = 1:3
        subplot(5, 3, (iIV - 1) * 3 + iDV); hold on;
        y = DV_core(:, iDV);
        valid = isfinite(x_raw) & isfinite(y);
        if ~any(valid)
            text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
            axis off;
            continue
        end
        xv = x_raw(valid);
        yv = y(valid);
        x_for_corr = xv;

        if iIV == 5
            % ridge row: plot in log10 steps and show unique ridge values as ticks
            pos = xv > 0;
            xv = xv(pos);
            yv = yv(pos);
            x_for_corr = xv;
            if isempty(xv)
                text(0.5, 0.5, 'No positive ridge', 'HorizontalAlignment', 'center');
                axis off;
                continue
            end
            xplot = log10(xv);
            scatter(xplot, yv, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45);
            uR = unique(xv);
            xticks(log10(uR));
            xticklabels(compose('%.3g', uR));
            if numel(uR) > 1
                xlim([min(log10(uR)) - 0.2, max(log10(uR)) + 0.2]);
            else
                xlim(log10(uR) + [-0.5, 0.5]);
            end
            % xtickangle(45);
            if numel(unique(xplot)) > 1
                pfit = polyfit(xplot, yv, 1);
                xr = linspace(min(xplot), max(xplot), 100);
                plot(xr, polyval(pfit, xr), 'r-', 'LineWidth', 1.5);
            end
        elseif isCategorical
            x_plot = xv + 0.2 * (rand(sum(valid), 1) - 0.5);
            scatter(x_plot, yv, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45);
            xticks(x_tick); xticklabels(x_ticklabel); %xtickangle(45);
            xlim([0.5, numel(x_tick) + 0.5]); % buffer on both ends for categorical rows
        else
            scatter(xv, yv, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45);
            if iIV == 1 || iIV == 2
                uInt = unique(round(xv));
                xticks(uInt);
                xlim([min(uInt) - 0.5, max(uInt) + 0.5]);
            end
            if numel(unique(xv)) > 1
                pfit = polyfit(xv, yv, 1);
                xr = linspace(min(xv), max(xv), 100);
                plot(xr, polyval(pfit, xr), 'r-', 'LineWidth', 1.5);
            end
        end
        if numel(x_for_corr) >= 3
            [rho, pval] = corr(x_for_corr, yv, 'Type', 'Spearman', 'Rows', 'complete');
            title(sprintf('\\rho=%.2f, p=%.3f', rho, pval));
        else
            title('n<3');
        end
        if iDV == 1, xlabel(IV_names{iIV}); end
        if iIV == 1, ylabel(dvNames_core{iDV}); end
        set(gca, 'XTickLabelRotation', 0);
        box on; grid on;
    end
end
sgtitle(sprintf('Correlation between basis settings and template recovery\n [A%d, L%d, %d iter]\n%s', ...
    iModelA_fit, iLocComb, nIter, subjTag));
set(findall(gcf, '-property', 'fontsize'), 'fontsize', sz_font);
saveas(gcf, sprintf('%s/4CorrBasisSetting_TempRMSE_A%d.jpg', nameFolder_Figures_perSubj, iModelA_fit));
close(gcf);

% Additional nBfit figures per A: parameter RMSE by B-fit model
paramNamesAll = {'Nmul', 'Nadd', 'Nshared', 'criterion_DV'};
trueParamVals = [Nmul_true, Nadd_true, Nshared_true, criterion_DV_true];

for iiB = 1:numel(iModelB_fit_all)
    iModelB_fit = iModelB_fit_all(iiB);
    nameFile_fitNOM = sprintf('%s/n%d_J%d_A%dB%d.mat', nameFolder_Data_NOM_IO, nIter, iJob, iModelA_fit, iModelB_fit);
    if ~exist(nameFile_fitNOM, 'file')
        fprintf('WARNING: fitNOM file not found for A%dB%d. Skip parameter-RMSE figure.\n', iModelA_fit, iModelB_fit);
        continue
    end

    Sfit = load(nameFile_fitNOM, 'params_est_allIter');
    if ~isfield(Sfit, 'params_est_allIter') || isempty(Sfit.params_est_allIter)
        fprintf('WARNING: params_est_allIter missing for A%dB%d. Skip parameter-RMSE figure.\n', iModelA_fit, iModelB_fit);
        continue
    end

    if iscell(Sfit.params_est_allIter)
        P = cell2mat(cellfun(@(x) x(:)', Sfit.params_est_allIter, 'UniformOutput', false));
    else
        P = Sfit.params_est_allIter;
    end

    paramNamesThisModel = namesModelBparams_short{iModelB_fit};
    nParamThis = numel(paramNamesThisModel);
    if nParamThis == 0
        continue
    end

    DV_param = nan(nIter, nParamThis);
    for iP = 1:nParamThis
        idxAll = find(strcmp(paramNamesAll, paramNamesThisModel{iP}), 1, 'first');
        if isempty(idxAll) || size(P, 2) < iP
            continue
        end
        nRows = min(nIter, size(P, 1));
        DV_param(1:nRows, iP) = abs(P(1:nRows, iP) - trueParamVals(idxAll));
    end

    figure('Position', [20 20 max(1200, 350*nParamThis) 1600]);
    rng(1);
    for iIV = 1:5
        [x_raw, x_tick, x_ticklabel, isCategorical] = encode_iv_for_scatter_local(IV_vals{iIV});
        for iDV = 1:nParamThis
            subplot(5, nParamThis, (iIV - 1) * nParamThis + iDV); hold on;
            y = DV_param(:, iDV);
            valid = isfinite(x_raw) & isfinite(y);
            if ~any(valid)
                text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
                axis off;
                continue
            end
            xv = x_raw(valid);
            yv = y(valid);
            x_for_corr = xv;

            if iIV == 5
                pos = xv > 0;
                xv = xv(pos);
                yv = yv(pos);
                x_for_corr = xv;
                if isempty(xv)
                    text(0.5, 0.5, 'No positive ridge', 'HorizontalAlignment', 'center');
                    axis off;
                    continue
                end
                xplot = log10(xv);
                scatter(xplot, yv, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45);
                uR = unique(xv);
                xticks(log10(uR));
                xticklabels(compose('%.3g', uR));
                if numel(uR) > 1
                    xlim([min(log10(uR)) - 0.2, max(log10(uR)) + 0.2]);
                else
                    xlim(log10(uR) + [-0.5, 0.5]);
                end
                % xtickangle(45);
                if numel(unique(xplot)) > 1
                    pfit = polyfit(xplot, yv, 1);
                    xr = linspace(min(xplot), max(xplot), 100);
                    plot(xr, polyval(pfit, xr), 'r-', 'LineWidth', 1.5);
                end
            elseif isCategorical
                x_plot = xv + 0.2 * (rand(sum(valid), 1) - 0.5);
                scatter(x_plot, yv, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45);
                xticks(x_tick); xticklabels(x_ticklabel); %xtickangle(45);
                xlim([0.5, numel(x_tick) + 0.5]);
            else
                scatter(xv, yv, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45);
                if iIV == 1 || iIV == 2
                    uInt = unique(round(xv));
                    xticks(uInt);
                    xlim([min(uInt) - 0.5, max(uInt) + 0.5]);
                end
                if numel(unique(xv)) > 1
                    pfit = polyfit(xv, yv, 1);
                    xr = linspace(min(xv), max(xv), 100);
                    plot(xr, polyval(pfit, xr), 'r-', 'LineWidth', 1.5);
                end
            end
            if numel(x_for_corr) >= 3
                [rho, pval] = corr(x_for_corr, yv, 'Type', 'Spearman', 'Rows', 'complete');
                title(sprintf('\\rho=%.2f, p=%.3f', rho, pval));
            else
                title('n<3');
            end
            if iDV == 1, xlabel(IV_names{iIV}); end
            if iIV == 1, ylabel(sprintf('RMSE %s', paramNamesThisModel{iDV})); end
            set(gca, 'XTickLabelRotation', 0);
            box on; grid on;
        end
    end

    sgtitle(sprintf('Correlation between basis settings and params recovery\n [A%dB%d, L%d, %d iter]\n%s', ...
        iModelA_fit, iModelB_fit, iLocComb, nIter, subjTag));
    set(findall(gcf, '-property', 'fontsize'), 'fontsize', sz_font);
    saveas(gcf, sprintf('%s/4CorrBasisSetting_ParamsRMSE_A%dB%d.jpg', nameFolder_Figures_perSubj, iModelA_fit, iModelB_fit));
    close(gcf);
end

fprintf('%s: Correlation between basis settings and params recovery saved for A%d.\n\n', datetime('now'), iModelA_fit);
end

%% helper
function out = rmse_simple_local(a, b)
a = a(:); b = b(:);
ok = isfinite(a) & isfinite(b);
if ~any(ok)
    out = NaN;
else
    d = a(ok) - b(ok);
    out = sqrt(mean(d.^2));
end
end

%% helper
function v = normalize_by_maxabs_local(v)
v = v(:)';
m = max(abs(v));
if isfinite(m) && m > 0
    v = v / m;
end
end

%% helper
function [x, xt, xtl, isCategorical] = encode_iv_for_scatter_local(vals)
if isnumeric(vals) || islogical(vals)
    x = vals(:);
    isCategorical = false;
    xt = [];
    xtl = {};
    return
end

if iscell(vals) || isstring(vals) || iscategorical(vals)
    labels = string(vals(:));
    bad = strlength(labels) == 0;
    labels(bad) = "<empty>";
    [u, ~, ic] = unique(labels, 'stable');
    x = double(ic);
    xt = 1:numel(u);
    xtl = cellstr(u);
    isCategorical = true;
    return
end

error('Unsupported IV type for scatter encoding.');
end

%% helper
function [pORI, pSF] = fit_marginal_params_from_template_local(template_ref, axis_tuning, iFamily_ORI, iFamily_SF, lb_full_all, ub_full_all, options_fmin)
pORI = [NaN NaN];
pSF = [NaN NaN];

try
    nRep_local = 20;
    problem_setting_local = MultiStart('StartPointsToRun', 'bounds', 'UseParallel', 1, 'Display', 'off');

    margORI_ref = mean(template_ref, 2)';
    margSF_ref = mean(template_ref, 1);

    xORI = axis_tuning{1};
    fxn_tuningLoss_ORI = @(param_est) sum((margORI_ref - predSFkernel(xORI, iFamily_ORI, param_est, 0)).^2);
    problem_ORI = createOptimProblem('fmincon', ...
        'objective', fxn_tuningLoss_ORI, ...
        'x0', (ub_full_all{iFamily_ORI} + lb_full_all{iFamily_ORI}) / 2, ...
        'lb', lb_full_all{iFamily_ORI}, ...
        'ub', ub_full_all{iFamily_ORI}, ...
        'options', options_fmin);
    pORI_fit = run(problem_setting_local, problem_ORI, nRep_local);

    xSF = axis_tuning{2};
    xSF_ln = 2.^xSF;
    fxn_tuningLoss_SF = @(param_est) sum((margSF_ref - predSFkernel(xSF_ln, iFamily_SF, param_est, 0)).^2);
    problem_SF = createOptimProblem('fmincon', ...
        'objective', fxn_tuningLoss_SF, ...
        'x0', (ub_full_all{iFamily_SF} + lb_full_all{iFamily_SF}) / 2, ...
        'lb', lb_full_all{iFamily_SF}, ...
        'ub', ub_full_all{iFamily_SF}, ...
        'options', options_fmin);
    pSF_fit = run(problem_setting_local, problem_SF, nRep_local);

    pORI(1:min(2, numel(pORI_fit))) = pORI_fit(1:min(2, numel(pORI_fit)));
    pSF(1:min(2, numel(pSF_fit))) = pSF_fit(1:min(2, numel(pSF_fit)));
catch
    % keep NaN so plotting can continue for other DVs
end
end

%% helper
function name = get_param_name_local(nameList, idx)
if iscell(nameList) && numel(nameList) >= idx
    name = nameList{idx};
else
    name = sprintf('param%d', idx);
end
end