function OOD_NOM_Trialwise_compIV_A12(nBasisORI, nBasisSF, isubj, iLocComb, lambda_whiten, flag_regressType, flag_whitenDV, nIter, nJob, iJob)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OOD_NOM_Trialwise_compIV.m
%
% Trial-wise noisy observer model – PRE-ESTIMATION STAGE
%
% Created by Shutian Xue on August 27, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all;
warning off; % (You may want to remove this once things are stable.)
format compact;

time_start = datetime('now');
fprintf('\n=======================================\n')
fprintf('Part 1: Compute DVs and derive templates')
fprintf('\n=======================================\n')

% fprintf('%s: Step 1 started.\n\n', time_start)

addpath(genpath('fxn_exp'));
addpath(genpath('fxn_NOM'));
addpath(genpath('fxn_RCplot'));
addpath(genpath('fxn_analysis_RC_v2'));
addpath(genpath('SX_toolbox/bads-master'));

%% Global settings
% Define directories
% --------------%
SX_RC1_setting; % defines nameFolder_*, nORI, nSF, namesLocComb, namesModelA, etc.
% --------------%

% nameFolder_Data = sprintf('%s/Data_Part2_nBasis%d%d_noMirror', nameFolder_server, nBasisORI, nBasisSF) ;
% nameFolder_Data_OOD = sprintf('%s/Data_OOD_%d%d', nameFolder_Data, nORI, nSF);  % Folder to save data
% nameFolder_Data_NOM_Trialwise = sprintf('%s/Data_NOM_Trialwise_%d%d', nameFolder_Data, nORI, nSF);  % Folder to save data

%% Deterministic RNG (grand seed + per-iteration substreams)
S_seed = GetGrandSeed(nIter, iJob, nJob, nameFolder_Data);

% One RNG stream for the whole job; each iteration uses its own Substream
stream = RandStream('Threefry', 'Seed', S_seed.grandSeed);
RandStream.setGlobalStream(stream);

% fprintf('%s: seed determined.\n\n', datetime('now'))

%% General parameters
IVType = 1; % 1=sum of the dot product/convolution; 2=max; 3=normalized
templateType = 1; % (1) raw (2) reconstructed kernel (3) mirrored template
itype_template = 2; % 1=estimate template from PRS trials, ABS trials, or BOTH trials
flag_PatchMode = 1; % if flag_PatchMode == 1, patchMode = 'T'; else, patchMode = 'N'; end
flag_plot_template = 0;

% Class weighting for template estimation (signal-absent emphasized)
wABS = 0.5;
wPRS = 1 - wABS;
weightScale_RC = 10; % integer replication scale for weighted regression
assert(wABS >= 0 && wABS <= 1, 'wABS must be in [0,1].');
assert(abs(wABS + wPRS - 1) < 1e-10, 'wABS + wPRS must equal 1.');

ratio_split = [.8, .15, .05]; % proportion of trials in template set (for RC), training set (for estimating parameters) and test set (for metric predictions)
assert(abs(sum(ratio_split)-1)<1e-10)
ORI_bound = [5, 14]; % orientation window, passed to fxn_getIV_v3

% whiten features
eps_whiten = 1e-3;  % floor for eigenvalues

% Setting for multivariate regression with smoothing
% 1=Univariate; 2=MultiSmooth and Univariate
if flag_regressType == 2
    opts = struct();
    % candidateORI = 3:9;
    % candidateSF = 3:9;
    candidateORI = 3:6;
    candidateSF = 3:6;
    % candidateSF = 6:8;
    candidateBasisFamilyORI = {'circ_gaussian', 'vonmises'};
    candidateBasisFamilyORI = {'vonmises'};
    candidateBasisFamilySF  = {'gaussianLog2', 'asymGaussianLog2', 'LogParabola', 'asymLogParabola'};
    % candidateBasisFamilySF  = {'asymGaussianLog2', 'asymLogParabola'};
    candidateRidge  = [0, 10.^(-2:3)];
    % candidateRidge  = 100;
    opts.basisWidthScaleORI = [.5:.2:.9]; % set vector to search ORI width scales
    opts.basisWidthScaleSF  = [.5:.2:.9]; % set vector to search SF width scales
    opts. asymSF_rightLeftRatio = 1.2;
    opts.nFolds = 5;
    opts.link = 'probit';
    opts.nBasisORI = nBasisORI;
    opts.nBasisSF = nBasisSF;
    opts.sigmaORI_deg = [];
    opts.sigmaSF_log2 = [];
    opts.oriPeriod_deg = 180;
    opts.zscorePredictor = true;
    opts.maxIter = 100;
    opts.tol = 1e-6;
end

% Settings for fitting tuning functions
nRep = 20;
iFamily_ORI = 1; % 1=scaled gaussian, 8=DoG
iFamily_SF = 2; % 2=log parabola
problem_setting = MultiStart('StartPointsToRun', 'bounds','UseParallel', 1, 'Display', 'off');

flag_plot_tuning = 0;

iSess_start = 1; % first session included
convolveType = 1; % IV from 1=cross-correlation; 2=convolution (fxn_getIV_v3)
flag_standEnergy = 1; % 1=z-score energy before RC
flag_plot_compIV = 1; % plot IV distributions and kernels at the end
if strcmp(str_envir,'HPC'), flag_plot_compIV = 0; end % don't plot when running on HPC

if flag_PatchMode == 1
    namePatchMode = 'T'; % target patch
else
    namePatchMode = 'N'; % noise patch
end

%% Subject / IO and paths
if isnumeric(isubj)
    % Human subjects
    subjList = {'YK','SP','SX','LS','RE','MD','AS','HL','FH','HA','CS','DT','DU','RC','SR'};
    nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195, 0];

    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);

    % Folder to load behav + energy
    nameFolder_OOD_load = sprintf('%s/%s%d', nameFolder_Data_OOD, subjName, nblocks);

    % Folder to save NOM trial-wise results
    nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocComb);

    % Folder to save figures
    nameFolder_Figures_perSubj = sprintf('%s/Human/%s', nameFolder_Figures, subjName);

    % Behavioral measurements
    load(sprintf('%s/%s%d/%s_behavMeas.mat', nameFolder_Data_OOD, subjName, nblocks, subjName), 'dataMatrix'); % ensure dataMatrix is loaded

    % Energy
    load(sprintf('%s/%s_energy_%s_%d_%d.mat', nameFolder_OOD_load, subjName, namePatchMode, nORI, nSF), ...
        'e3D_target_allT', 'e3D_noise_allT', 'noise', 'filtersOri_all', 'stim', 'nBins');

else % Ideal observer (IO)
    subjName = isubj{1}; % "nameIO" from OOD_sim
    nblocks = 0;

    % Folder to load behav + energy
    nameFolder_OOD_load = sprintf('%s/%s', nameFolder_Data_OOD, subjName);

    % Folder to save NOM trial-wise results
    nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, subjName);

    % Folder to save figures
    nameFolder_Figures_perSubj = sprintf('%s/IO/%s', nameFolder_Figures, subjName);

    % Behavioral measurements
    load(sprintf('%s/behavMeas.mat', nameFolder_OOD_load), 'dataMatrix');

    % Energy
    load(sprintf('%s/energy_%s_%d_%d.mat', nameFolder_OOD_load, namePatchMode, nORI, nSF), ...
        'e3D_target_allT', 'e3D_noise_allT', 'noise', 'filtersOri_all', 'stim', 'nBins');
end

% Create the directory if absent
if isempty(dir(nameFolder_NOM_save))
    mkdir(nameFolder_NOM_save);
end

if isempty(dir(nameFolder_Figures_perSubj))
    mkdir(nameFolder_Figures_perSubj);
end

%% Print run info
fprintf(['\nSubject/IO name: %s ' ...
    '\n - L%d [%s]', ...
    '\n - Number of iterations = %d', ...
    '\n - Number of blocks = %d', ...
    '\n - nORI=%d | nSF=%d', ...
    '\n - Template derived from %s trials', ...
    '\n - RC class weights: wABS=%.2f, wPRS=%.2f (scale=%d)', ...
    '\n - Energy source: %d (1=TARGET, 2=NOISE)', ...
    '\n - Convolve type: %s', ...
    '\n - IV type: %s\n\n'], ...
    subjName, ...
    iLocComb, namesLocComb{iLocComb}, ...
    nIter, ...
    nblocks, nORI, nSF, ...
    namesType{itype_template}, wABS, wPRS, weightScale_RC, flag_PatchMode, ...
    namesConvolveType{convolveType}, namesIVType{IVType});

%% Choose energy source for NOM
switch flag_PatchMode
    case 1 % target patch energy
        e3D_allT = e3D_target_allT;
    case 2 % noise patch energy
        e3D_allT = e3D_noise_allT;
    otherwise
        error('flag_PatchMode must be 1 (target) or 2 (noise).');
end

%% Sanity check: iPRS consistent within pairs
% Col#1 = itrial_allT; Col#6 = iPRS; Col#8 = iPair;
for iPair = 1:max(dataMatrix(:, 8))
    iIter = find(dataMatrix(:, 8) == iPair);
    iPRS_A = dataMatrix(dataMatrix(:, 1) == iIter(1), 6);
    iPRS_B = dataMatrix(dataMatrix(:, 1) == iIter(2), 6);
    if iPRS_A ~= iPRS_B
        fprintf('WARNING: when ipair=%d, iPRS_A=%d, iPRS_B=%d\n', iPair, iPRS_A, iPRS_B);
    end
end

%% Session and location selection
nSess = length(unique(dataMatrix(:, 2)));

% Sessions included (can be subsetted here)
iSess_select = iSess_start:nSess;

% Special case: subject 'RE' missing some sessions
if strcmp(subjName, 'RE')
    iSess_select = [1:15, 17:23, 25:nSess];
end

% Map iLocComb to actual location indices
if iLocComb < 6
    iLoc_all = iLocComb;
else
    switch iLocComb
        case 6, iLoc_all = [2, 4]; % horizontal meridian
        case 7, iLoc_all = [3, 5]; % vertical meridian
        case 8, iLoc_all = 2:5; % all perifoveal locations
        otherwise
            error('Unknown iLocComb = %d', iLocComb);
    end
end

% Precompute non-random pools
e3D_nonrand = [];
cst_nonrand = [];

for iSS = 1:length(iSess_select)
    inSess = (dataMatrix(:, 2) == iSess_select(iSS));
    inLoc = ismember(dataMatrix(:, 5), iLoc_all); % supports single or multiple locations
    indLoc = inSess & inLoc;
    e3D_nonrand = cat(1, e3D_nonrand, e3D_allT(indLoc, :, :));
    cst_nonrand = [cst_nonrand; dataMatrix(indLoc,11)];
end

% The gabor contrs
gaborCST = mean(cst_nonrand);
templateType_true = 1; % doesn't matter much for the IO Gabor template

stim.gaborCST = gaborCST;

% Create a pool of Gabor filters
[filter_sin, filter_cos] = SX_sim02_setFilters(stim, noise.filtersSF_all, filtersOri_all, fxn_getSigma_SPdomain, 0);

%% Define the IDEAL template (ORI x SF)
template_gabor_ideal = exp_CreateGabor(stim, stim.gaborCST);
%--------------------------------------------%
template_ideal = SX_RC4_Energy_parfor(stim.mask, {template_gabor_ideal}, filter_sin, filter_cos);
%--------------------------------------------%
template_ideal = squeeze(template_ideal); % remove singleton dim
%--------------------------------------------%
template_ideal = fxn_getTemplate(template_ideal, templateType_true, 0);
%--------------------------------------------%

% Normalize template (Unit L2-norm)
template_ideal = template_ideal / norm(template_ideal(:));

fprintf('%s: Ideal template created.\n\n', datetime('now'))

%% Stage overview
% A1: estimate trial-wise templates from data and compute DVs.
% A2: keep trial resampling/whitening/DV pipeline the same, but replace the
% template with a fixed "true" template scaled to A1's mean derived amplitude.
% Execution order matters because A2 depends on A1 outputs.

%% A1
% Output file name (before estimation)
iModelA_fit = 1;
nameFile_compIV_A1 = sprintf('%s/n%d_J%d_A%d_compIV', nameFolder_NOM_save, nIter, iJob, iModelA_fit);

%% A1—Main loop over iterations
% fprintf('%s: Creating empty placeholders for running iterations.\n\n', datetime('now'))
data_metrics_allIter = nan(nIter, nDatasets_full, nMetrics); % see fxn_getMetrics for all 11 metrics
data_train_allIter = cell(nIter, 1);
data_test_allIter = data_train_allIter;
template_tmpl_allIter = nan(nIter, nORI, nSF);
template_full_allIter = template_tmpl_allIter;
nBasisORI_tmpl_allIter = nan(nIter, 1);
nBasisSF_tmpl_allIter = nan(nIter, 1);
basisFxnORI_tmpl_allIter = data_train_allIter;
basisFxnSF_tmpl_allIter = data_train_allIter;
ridge_tmpl_allIter = nBasisSF_tmpl_allIter;
basisWidthScaleORI_tmpl_allIter = nBasisSF_tmpl_allIter;
basisWidthScaleSF_tmpl_allIter  = nBasisSF_tmpl_allIter;

sep_allIter = nan(nIter, 2); % 1=template set, 2=full set
margORI_allIter = nan(nIter, 2, nORI);
margPred_ORI_allIter = nan(nIter, 2, nORI);
margParams_ORI_allIter = nan(nIter, 2, length(namesParams_all{iFamily_ORI}));
margR2_ORI_allIter = nan(nIter, 2);
margSF_allIter = nan(nIter, 2, nSF);
margPred_SF_allIter = nan(nIter, 2, nSF);
margParams_SF_allIter = nan(nIter, 2, length(namesParams_all{iFamily_SF}));
margR2_SF_allIter = nan(nIter, 2);

nameFile_progress = [nameFile_compIV_A1, '.mat'];

fprintf('%s: A1 Started running %d iterations.\n\n', datetime('now'), nIter)

for iIter = 1:nIter
    fprintf('\n%s: %d...', datetime('now'), iIter);

    % Deterministic randomness for THIS iteration (global index across jobs)
    stream.Substream = S_seed.iterIdxList(iIter);

    %% 1. Resample trials into FULL or TEMPLATE/TRAIN / TEST
    %----------------%
    fxn_resampleTrials;
    %----------------%

    %% 2. Select which trials to use to estimate the template (TEMPLATE set)
    [useIdx_full, useIdx_tmpl] = selectTemplateIndices(itype_template, iPRS_full_rand, iPRS_tmpl_rand);

    % sel for "selected"
    e3D_tmpl_rand_sel = e3D_tmpl_rand(useIdx_tmpl, :, :);
    cst_tmpl_rand_sel = cst_tmpl_rand(useIdx_tmpl);
    resp_tmpl_rand_sel = resp_tmpl_rand(useIdx_tmpl);
    iPRS_tmpl_rand_sel = iPRS_tmpl_rand(useIdx_tmpl);

    e3D_full_rand_sel = e3D_full_rand(useIdx_full, :, :);
    cst_full_rand_sel = cst_full_rand(useIdx_full);
    resp_full_rand_sel = resp_full_rand(useIdx_full);
    iPRS_full_rand_sel = iPRS_full_rand(useIdx_full);

    % Standardize energy if requested
    if flag_standEnergy
        %--------------------------------------------%
        e3D_tmpl_rand_norm = normEnergy(e3D_tmpl_rand_sel, cst_tmpl_rand_sel, iPRS_tmpl_rand_sel);
        e3D_full_rand_norm = normEnergy(e3D_full_rand_sel, cst_full_rand_sel, iPRS_full_rand_sel);
        %--------------------------------------------%
    else
        e3D_tmpl_rand_norm = e3D_tmpl_rand_sel;
        e3D_full_rand_norm = e3D_full_rand_sel;
    end

    %% Deal with channel correlation
    % Two alternative strategies (lambda_whiten = NaN  --> Way 2)
    % Way 1: whiten the feature/channel energy before template estimation
    % Way 2: estimate the template in the original feature space, then apply
    %   - this is closer in spirit to AE1999 / AE2002 covariance correction

    % -------------------------------------------------------------------------
    % 1. Estimate covariance in channel-energy space separately for TEMPLATE and FULL
    % -------------------------------------------------------------------------
    [mu_cov_tmpl, ~, W_white_tmpl] = computeWhiteningParams(e3D_tmpl_rand_norm, lambda_whiten, eps_whiten);
    [mu_cov_full, ~, W_white_full] = computeWhiteningParams(e3D_full_rand_norm, lambda_whiten, eps_whiten);

    % -------------------------------------------------------------------------
    % 2. Choose covariance-handling strategy
    % -------------------------------------------------------------------------
    if isnan(lambda_whiten)
        error('lambda_whiten = NaN path is not implemented in this script variant (A12).');

    else
        % =====================================================================
        % Way 1: whiten the channel-energy predictors before template estimation
        % =====================================================================

        % lambda_whiten = 0 --> use the full empirical covariance
        % lambda_whiten = 1 --> use a scaled identity covariance (ignore channel correlations)

        % Apply whitening to each dataset using its own mean and covariance
        e3D_tmpl_use = whiten_e3D(e3D_tmpl_rand_norm, mu_cov_tmpl, W_white_tmpl);
        e3D_full_use = whiten_e3D(e3D_full_rand_norm, mu_cov_full, W_white_full);

        % Safe defaults for branch-specific outputs (used only for case 2)
        nBasisORI_tmpl = nan;
        nBasisSF_tmpl = nan;
        basisFxnORI_tmpl = '';
        basisFxnSF_tmpl = '';
        basisWidthScaleORI_tmpl = nan;
        basisWidthScaleSF_tmpl = nan;
        ridge_tmpl = nan;

        % Apply class weighting before template estimation (ABS emphasized).
        [e3D_tmpl_forRC, resp_tmpl_forRC] = applyClassWeightsForRC(e3D_tmpl_use, resp_tmpl_rand_sel, iPRS_tmpl_rand_sel, wABS, wPRS, weightScale_RC);
        [e3D_full_forRC, resp_full_forRC] = applyClassWeightsForRC(e3D_full_use, resp_full_rand_sel, iPRS_full_rand_sel, wABS, wPRS, weightScale_RC);

        % Estimate template in whitened feature space
        switch flag_regressType
            case 1  % Univariate
                template_tmpl_raw = SX_sim07_RC(e3D_tmpl_forRC, resp_tmpl_forRC);
                template_full_raw = SX_sim07_RC(e3D_full_forRC, resp_full_forRC);

            case 2  % Multivariate + smoothing
                opts.template_ideal = template_ideal;
                out = SX_RC_selectBasis_cv(e3D_tmpl_forRC, resp_tmpl_forRC, axis_tuning{1}, axis_tuning{2}, candidateORI, candidateSF, candidateBasisFamilyORI, candidateBasisFamilySF, candidateRidge, opts);
                template_tmpl_raw = out.template2D;
                nBasisORI_tmpl = out.nBasisORI;
                nBasisSF_tmpl = out.nBasisSF;
                basisFxnORI_tmpl = out.basisFamilyORI;
                basisFxnSF_tmpl = out.basisFamilySF;
                basisWidthScaleORI_tmpl = out.basisWidthScaleORI;
                basisWidthScaleSF_tmpl = out.basisWidthScaleSF;
                ridge_tmpl = out.ridge;

                out = SX_RC_selectBasis_cv(e3D_full_forRC, resp_full_forRC, axis_tuning{1}, axis_tuning{2}, candidateORI, candidateSF, candidateBasisFamilyORI, candidateBasisFamilySF, candidateRidge, opts);
                template_full_raw = out.template2D;

        end % switch

        % No back-transformation here:
        % template_tmpl_raw and template_full_raw now live in their respective
        % whitened feature spaces and should stay there for subsequent DV computation.
    end % isnan()

    %% Regularize the derived template
    template_full = fxn_getTemplate(template_full_raw, templateType, flag_plot_template);
    template_tmpl = fxn_getTemplate(template_tmpl_raw, templateType, flag_plot_template);

    %% Check how channel correlation affects template estimation
    % -----PCA-based template recovery =====
    if flag_plot_compIV
        % NOMplot_PCA
    end

    % ----- Debug: check whether channels are too correlated thus redundant
    % !!!! Generate one figure per iteration !!!!
    if flag_plot_compIV
        % NOMplot_checkChannelCorr
    end

    % ----- Compare the univariate regression and multivariation regression with smoothing
    % One of the two may not be generated; need to pause and generate manually
    if flag_plot_compIV
        % template_smooth = template_full_raw;
        % template_noSmoothing = SX_sim07_RC(e3D_full_rand_norm, resp_full_rand_sel);
        % RCplot_compTempSmoothing
    end

    %% 4. Marginalization and fit tuning functions
    template_notNormed_tmpl = template_tmpl;
    template_notNormed_full = template_full;

    for iDataset = 1:2
        switch iDataset
            case 1
                template_notNormed = template_notNormed_tmpl;
            case 2
                template_notNormed = template_notNormed_full;
        end

        % Calculate separability
        templateType_recon = 2; % 2=reconstructed template
        %--------------------------------------------%
        template_recon = fxn_getTemplate(template_notNormed, templateType_recon, flag_plot_template);
        %--------------------------------------------%
        sep = corr(template_notNormed(:), template_recon(:));

        % Marginalize the template into 1D ORI and SF profiles
        margORI = mean(template_notNormed, 2)';
        margSF = mean(template_notNormed, 1);

        % Fitting
        if any(isnan(margORI)), error('ALERT: NaN in margORI!'), end

        xORI = axis_tuning{1};
        fxn_tuningLoss_ORI = @(param_est) sum((margORI - predSFkernel(xORI, iFamily_ORI, param_est, flag_plot_tuning)).^2);
        problem_ORI = createOptimProblem('fmincon','objective', fxn_tuningLoss_ORI,'x0', (ub_full_all{iFamily_ORI}+lb_full_all{iFamily_ORI})/2,'lb',lb_full_all{iFamily_ORI},'ub',ub_full_all{iFamily_ORI},'options',options_fmin);
        margParams_ORI = run(problem_setting, problem_ORI, nRep);
        %--------------------------------------------%
        margPred_ORI = predSFkernel(xORI, iFamily_ORI, margParams_ORI, flag_plot_tuning);
        %--------------------------------------------%
        margR2_ORI = 1-sumsqr(margORI-margPred_ORI)/sumsqr(margORI-mean(margORI));

        xSF = axis_tuning{2};
        xSF_ln = 2.^xSF; % fit in linear SF space
        fxn_tuningLoss_SF = @(param_est) sum((margSF - predSFkernel(xSF_ln, iFamily_SF, param_est, flag_plot_tuning)).^2);
        problem_SF = createOptimProblem('fmincon','objective', fxn_tuningLoss_SF,'x0', (ub_full_all{iFamily_SF}+lb_full_all{iFamily_SF})/2,'lb',lb_full_all{iFamily_SF},'ub',ub_full_all{iFamily_SF},'options',options_fmin);
        margParams_SF = run(problem_setting, problem_SF, nRep);
        %--------------------------------------------%
        margPred_SF = predSFkernel(xSF_ln, iFamily_SF, margParams_SF, flag_plot_tuning);
        %--------------------------------------------%
        margR2_SF = 1-sumsqr(margSF-margPred_SF)/sumsqr(margSF-mean(margSF));

        % Store
        sep_allIter(iIter, iDataset) = sep;

        margORI_allIter(iIter, iDataset, :) = margORI;
        margPred_ORI_allIter(iIter, iDataset, :) = margPred_ORI;
        margParams_ORI_allIter(iIter, iDataset, :) = margParams_ORI;
        margR2_ORI_allIter(iIter, iDataset) = margR2_ORI;

        margSF_allIter(iIter, iDataset, :) = margSF;
        margPred_SF_allIter(iIter, iDataset, :) = margPred_SF;
        margParams_SF_allIter(iIter, iDataset, :) = margParams_SF;
        margR2_SF_allIter(iIter, iDataset) = margR2_SF;
    end % iDataset

    %% 5-8. Metrics, DV, binning, data structs
    d_full  = struct('resp', resp_full_rand,  'iPRS', iPRS_full_rand,  'respC', respC_full_rand,  'cst', cst_full_rand,  'RT', RT_full_rand);
    d_tmpl  = struct('resp', resp_tmpl_rand,  'iPRS', iPRS_tmpl_rand,  'respC', respC_tmpl_rand,  'cst', cst_tmpl_rand,  'RT', RT_tmpl_rand);
    d_train = struct('resp', resp_train_rand, 'iPRS', iPRS_train_rand, 'respC', respC_train_rand, 'cst', cst_train_rand, 'RT', RT_train_rand, 'iPair', iPair_train_rand);
    d_test  = struct('resp', resp_test_rand,  'iPRS', iPRS_test_rand,  'respC', respC_test_rand,  'cst', cst_test_rand,  'RT', RT_test_rand,  'iPair', iPair_test_rand);
    [data_train, data_test, metrics] = fxn_compDV_packData( ...
        d_full, d_tmpl, d_train, d_test, ...
        e3D_train_rand, e3D_test_rand, ...
        flag_whitenDV, lambda_whiten, mu_cov_full, W_white_full, ...
        template_notNormed, convolveType, IVType, ORI_bound, nBins);
    data_train_allIter{iIter} = data_train;
    data_test_allIter{iIter}  = data_test;
    data_metrics_allIter(iIter, :, :) = metrics;

    % Store templates
    template_tmpl_allIter(iIter, :, :) = template_notNormed_tmpl;
    template_full_allIter(iIter, :, :) = template_notNormed_full;
    nBasisORI_tmpl_allIter(iIter, :) = nBasisORI_tmpl;
    nBasisSF_tmpl_allIter(iIter, :) = nBasisSF_tmpl;
    basisFxnORI_tmpl_allIter{iIter} = basisFxnORI_tmpl;
    basisFxnSF_tmpl_allIter{iIter} = basisFxnSF_tmpl;
    basisWidthScaleORI_tmpl_allIter(iIter) = basisWidthScaleORI_tmpl;
    basisWidthScaleSF_tmpl_allIter(iIter) = basisWidthScaleSF_tmpl;
    ridge_tmpl_allIter(iIter) = ridge_tmpl;

    %% Save a progress report in the folder to indicate the finished iteration and time spent
    time_progress = datetime('now');
    time_progress = ceil(minutes(time_progress-time_start)); % round up to minutes
    save(nameFile_progress, 'iIter')

    % rename (to avoid saving one file for each iteration)
    nameFile_progress_new = sprintf('%s_%d_%dmin.mat', nameFile_compIV_A1, iIter, time_progress);
    movefile(nameFile_progress, nameFile_progress_new);
    nameFile_progress = nameFile_progress_new;

end % end for iIter

fprintf('\n\n%s: A1 All iterations done.\n\n', datetime('now'))

%% SAVE  A1
% save(nameFile_compIV_A1, 'time_progress', 'template_ideal', 'criterion_z*', '*_allIter', 'names*', 'flag*', 'ratio_split', 'ORI_bound', '*Type');
save(nameFile_compIV_A1, 'time_progress', 'template_ideal', '*_allIter', 'names*', 'flag*', 'ratio_split', 'ORI_bound', '*Type');

% Delete the progress report
if exist(nameFile_progress_new, 'file')
    delete(nameFile_progress_new);
end

fprintf('\n\n%s: A1 Outputs saved.\n\n', datetime('now'))

%% Plot A1 (optional)
if flag_plot_compIV
    %-------------------%
    NOMplot_compIV;
    %-------------------%
    fprintf('\n\n%s: A1 Plots created.\n\n', datetime('now'))
end
close all;

%% End timing for A1
time_end = datetime('now');
fprintf('%s: A1 Compute DV done.\n\n', time_end)
elapsed = time_end - time_start;
fprintf('A1 Time used: %s\n\n\n\n', char(elapsed));

%% A2
iModelA_fit=2;
nameFile_compIV_A2 = sprintf('%s/n%d_J%d_A%d_compIV', nameFolder_NOM_save, nIter, iJob, iModelA_fit);

%% Process the true template
% Average the derived template (tmpl set) across iterations
template_tmpl_ave_A1 = squeeze(mean(template_tmpl_allIter, 1, 'omitnan'));
peak_tmpl_A1 = max(abs(template_tmpl_ave_A1(:)));
% Rescale the true/ideal template to match the peak of the derived template
% (apply scaling right after template_notNormed is loaded below)
        if ~isnumeric(isubj) % for IO, you may reload 'template_true' from disk
            load(sprintf('%s/truth.mat', nameFolder_OOD_load), 'template_true');
            template_notNormed = template_true;
        else
            template_notNormed = template_ideal;
        end

if isfinite(peak_tmpl_A1) && peak_tmpl_A1 > 0
    peak_true = max(abs(template_notNormed(:)));
    if isfinite(peak_true) && peak_true > 0
        template_notNormed = template_notNormed * (peak_tmpl_A1 / peak_true);
    end
end

%% A2: Main loop over iterations
% fprintf('%s: Creating empty placeholders for running iterations.\n\n', datetime('now'))
data_metrics_allIter = nan(nIter, nDatasets_full, nMetrics); % see fxn_getMetrics for all 11 metrics
data_train_allIter = cell(nIter, 1);
data_test_allIter = data_train_allIter;
template_tmpl_allIter = nan(nIter, nORI, nSF);
template_full_allIter = template_tmpl_allIter;
nameFile_progress = [nameFile_compIV_A2, '.mat'];

fprintf('%s: A2 Started running %d iterations.\n\n', datetime('now'), nIter)

for iIter = 1:nIter
    fprintf('\n%s: %d...', datetime('now'), iIter);

    % Deterministic randomness for THIS iteration (global index across jobs)
    stream.Substream = S_seed.iterIdxList(iIter);

    %% 1. Resample trials into FULL or TEMPLATE/TRAIN / TEST
    %----------------%
    fxn_resampleTrials;
    %----------------%

    %% 2. Select which trials to use to estimate the template (TEMPLATE set)
    [useIdx_full, useIdx_tmpl] = selectTemplateIndices(itype_template, iPRS_full_rand, iPRS_tmpl_rand);

    % Use the FULL selected set to define whitening for DV in this iteration
    e3D_full_rand_sel = e3D_full_rand(useIdx_full, :, :);
    cst_full_rand_sel = cst_full_rand(useIdx_full);
    iPRS_full_rand_sel = iPRS_full_rand(useIdx_full);

    if flag_standEnergy
        e3D_full_rand_norm = normEnergy(e3D_full_rand_sel, cst_full_rand_sel, iPRS_full_rand_sel);
    else
        e3D_full_rand_norm = e3D_full_rand_sel;
    end

    [mu_cov_full, ~, W_white_full] = computeWhiteningParams(e3D_full_rand_norm, lambda_whiten, eps_whiten);

    %% 5-8. Metrics, DV, binning, data structs
    d_full  = struct('resp', resp_full_rand,  'iPRS', iPRS_full_rand,  'respC', respC_full_rand,  'cst', cst_full_rand,  'RT', RT_full_rand);
    d_tmpl  = struct('resp', resp_tmpl_rand,  'iPRS', iPRS_tmpl_rand,  'respC', respC_tmpl_rand,  'cst', cst_tmpl_rand,  'RT', RT_tmpl_rand);
    d_train = struct('resp', resp_train_rand, 'iPRS', iPRS_train_rand, 'respC', respC_train_rand, 'cst', cst_train_rand, 'RT', RT_train_rand, 'iPair', iPair_train_rand);
    d_test  = struct('resp', resp_test_rand,  'iPRS', iPRS_test_rand,  'respC', respC_test_rand,  'cst', cst_test_rand,  'RT', RT_test_rand,  'iPair', iPair_test_rand);
    [data_train, data_test, metrics] = fxn_compDV_packData( ...
        d_full, d_tmpl, d_train, d_test, ...
        e3D_train_rand, e3D_test_rand, ...
        flag_whitenDV, lambda_whiten, mu_cov_full, W_white_full, ...
        template_notNormed, convolveType, IVType, ORI_bound, nBins);
    data_train_allIter{iIter} = data_train;
    data_test_allIter{iIter}  = data_test;
    data_metrics_allIter(iIter, :, :) = metrics;

    % Store templates
    % if iModelA_fit==2
    template_notNormed_tmpl = template_notNormed;
    template_notNormed_full = template_notNormed;
    % end
    template_tmpl_allIter(iIter, :, :) = template_notNormed_tmpl;
    template_full_allIter(iIter, :, :) = template_notNormed_full;

    %% A2: Save a progress report in the folder to indicate the finished iteration and time spent
    time_progress = datetime('now');
    time_progress = ceil(minutes(time_progress-time_start)); % round up to minutes
    save(nameFile_progress, 'iIter')

    % rename (to avoid saving one file for each iteration)
    nameFile_progress_new = sprintf('%s_%d_%dmin.mat', nameFile_compIV_A2, iIter, time_progress);
    movefile(nameFile_progress, nameFile_progress_new);
    nameFile_progress = nameFile_progress_new;

end % end for iIter

fprintf('\n\n%s: A2 All iterations done.\n\n', datetime('now'))

%% SAVE for A2
save(nameFile_compIV_A2, 'time_progress', 'template_ideal', '*_allIter', 'names*', 'flag*', 'ratio_split', 'ORI_bound', '*Type');

% Delete the progress report
if exist(nameFile_progress_new, 'file')
    delete(nameFile_progress_new);
end

fprintf('\n\n%s: A2 Outputs saved.\n\n', datetime('now'))
end

%% HELPER
function [data_train, data_test, metrics] = fxn_compDV_packData( ...
    d_full, d_tmpl, d_train, d_test, ...
    e3D_train_rand, e3D_test_rand, ...
    flag_whitenDV, lambda_whiten, mu_cov_full, W_white_full, ...
    template_notNormed, convolveType, IVType, ORI_bound, nBins)

% Sections 5-8 shared by A1 and A2 loops:
%   behavioral metrics, DV from raw/whitened energy, bin IV, compile data structs.
% d_full/tmpl : struct with fields resp, iPRS, respC, cst, RT
% d_train/test: struct with fields resp, iPRS, respC, cst, RT, iPair

%% 5. Compute behavioral metrics
metrics_full  = fxn_getMetrics(d_full.resp,  d_full.iPRS,  d_full.respC,  d_full.cst,  d_full.RT);
metrics_tmpl  = fxn_getMetrics(d_tmpl.resp,  d_tmpl.iPRS,  d_tmpl.respC,  d_tmpl.cst,  d_tmpl.RT);
metrics_train = fxn_getMetrics(d_train.resp, d_train.iPRS, d_train.respC, d_train.cst, d_train.RT);
metrics_test  = fxn_getMetrics(d_test.resp,  d_test.iPRS,  d_test.respC,  d_test.cst,  d_test.RT);
% [dprime, criterion, pC, pHit, pFA, nanmean(respC), pYES, mean(1./contrast), median(RT)]
metrics = [metrics_full; metrics_tmpl; metrics_train; metrics_test];
criterion_z_train = metrics_train(2);
criterion_z_test  = metrics_test(2);

%% 6. Compute decision variable (DV) from raw or whitened energy
flag_permT = 0;
if flag_whitenDV && ~isnan(lambda_whiten)
    e3D_train_forDV = whiten_e3D(e3D_train_rand, mu_cov_full, W_white_full);
    e3D_test_forDV  = whiten_e3D(e3D_test_rand,  mu_cov_full, W_white_full);
else
    e3D_train_forDV = e3D_train_rand;
    e3D_test_forDV  = e3D_test_rand;
end
DV_train    = fxn_getIV_v3(e3D_train_forDV, template_notNormed, convolveType, IVType, flag_permT, ORI_bound);
DV_test     = fxn_getIV_v3(e3D_test_forDV,  template_notNormed, convolveType, IVType, flag_permT, ORI_bound);
nData_train = length(DV_train);
nData_test  = length(DV_test);

%% 7. Bin IV values
[nTrials_PRS_allBins_train, ~, iTrial4Bin_PRS_train] = histcounts(DV_train(d_train.iPRS == 1), nBins);
[nTrials_ABS_allBins_train, ~, iTrial4Bin_ABS_train] = histcounts(DV_train(d_train.iPRS == 0), nBins);
[nTrials_PRS_allBins_test,  ~, iTrial4Bin_PRS_test]  = histcounts(DV_test(d_test.iPRS  == 1), nBins);
[nTrials_ABS_allBins_test,  ~, iTrial4Bin_ABS_test]  = histcounts(DV_test(d_test.iPRS  == 0), nBins);
[nTrials_allBins_train, ~, iTrial4Bin_train] = histcounts(DV_train, nBins);
[nTrials_allBins_test,  ~, iTrial4Bin_test]  = histcounts(DV_test,  nBins);

%% 8. Compile data structs
data_train = struct();
data_train.criterion_z         = criterion_z_train;
data_train.ndata               = nData_train;
data_train.IV                  = DV_train;
data_train.nTrials_allBins     = nTrials_allBins_train;
data_train.iTrial4Bin          = iTrial4Bin_train;
data_train.nTrials_PRS_allBins = nTrials_PRS_allBins_train;
data_train.iTrial4Bin_PRS      = iTrial4Bin_PRS_train;
data_train.nTrials_ABS_allBins = nTrials_ABS_allBins_train;
data_train.iTrial4Bin_ABS      = iTrial4Bin_ABS_train;
data_train.iPRS                = d_train.iPRS;
data_train.resp                = d_train.resp;
data_train.cst                 = d_train.cst;
data_train.iPair               = d_train.iPair;
data_train.respC               = d_train.respC;
data_train.RT                  = d_train.RT;
data_train.metrics_sim         = metrics_train; % keep field name 'metrics_sim'

data_test = struct();
data_test.criterion_z         = criterion_z_test;
data_test.ndata               = nData_test;
data_test.IV                  = DV_test;
data_test.nTrials_allBins     = nTrials_allBins_test;
data_test.iTrial4Bin          = iTrial4Bin_test;
data_test.nTrials_PRS_allBins = nTrials_PRS_allBins_test;
data_test.iTrial4Bin_PRS      = iTrial4Bin_PRS_test;
data_test.nTrials_ABS_allBins = nTrials_ABS_allBins_test;
data_test.iTrial4Bin_ABS      = iTrial4Bin_ABS_test;
data_test.iPRS                = d_test.iPRS;
data_test.resp                = d_test.resp;
data_test.cst                 = d_test.cst;
data_test.iPair               = d_test.iPair;
data_test.respC               = d_test.respC;
data_test.RT                  = d_test.RT;
data_test.metrics_sim         = metrics_test; % keep field name 'metrics_sim'
end

function e3D_white = whiten_e3D(e3D_in, mu_cov, W_white)
[nT, nOri, nSf] = size(e3D_in);
X = reshape(e3D_in, [nT, nOri * nSf]);
Xw = (X - mu_cov) * W_white;
e3D_white = reshape(Xw, [nT, nOri, nSf]);
end

function [useIdx_full, useIdx_tmpl] = selectTemplateIndices(itype_template, iPRS_full_rand, iPRS_tmpl_rand)
% Shared selector for template/full trial subsets used by A1 and A2.
switch itype_template
    case 1 % PRS trials only
        useIdx_full = (iPRS_full_rand == 1);
        useIdx_tmpl = (iPRS_tmpl_rand == 1);
    case 2 % ABS trials only
        useIdx_full = (iPRS_full_rand == 0);
        useIdx_tmpl = (iPRS_tmpl_rand == 0);
    otherwise % 3 = both PRS and ABS
        useIdx_full = true(size(iPRS_full_rand));
        useIdx_tmpl = true(size(iPRS_tmpl_rand));
end
end

function [mu_cov, Sigma, W_white] = computeWhiteningParams(e3D_norm, lambda_whiten, eps_whiten)
% Compute mean/covariance and optional whitening transform in channel space.
[nTrials_cov, nOri_cov, nSf_cov] = size(e3D_norm);
e3D_vec = reshape(e3D_norm, [nTrials_cov, nOri_cov * nSf_cov]);
mu_cov = mean(e3D_vec, 1);
Sigma = cov(e3D_vec - mu_cov, 1); % population covariance

if isnan(lambda_whiten)
    W_white = [];
    return;
end

Sigma_shrink = (1 - lambda_whiten) * Sigma + lambda_whiten * mean(diag(Sigma)) * eye(size(Sigma, 1));
[V, D] = eig((Sigma_shrink + Sigma_shrink') / 2);
d = diag(D);
d(d < eps_whiten) = eps_whiten;
W_white = V * diag(1 ./ sqrt(d)) * V';
end

function kernel_corr = invcov_kernel(kernel_raw, W_inv, nOri, nSf)
k_raw = kernel_raw(:);
k_corr = W_inv * k_raw;
kernel_corr = reshape(k_corr, [nOri, nSf]);
end

function [e3D_out, resp_out] = applyClassWeightsForRC(e3D_in, resp_in, iPRS_in, wABS, wPRS, weightScale)
% Replicate trials to approximate class weights for regression-based template estimation.
%   ABS trials (iPRS==0) are replicated round(wABS*weightScale) times.
%   PRS trials (iPRS==1) are replicated round(wPRS*weightScale) times.
nABS = round(wABS * weightScale);
nPRS = round(wPRS * weightScale);
idxABS = find(iPRS_in == 0);
idxPRS = find(iPRS_in == 1);
idx_rep = [repmat(idxABS, nABS, 1); repmat(idxPRS, nPRS, 1)];
e3D_out = e3D_in(idx_rep, :, :);
resp_out = resp_in(idx_rep);
end
