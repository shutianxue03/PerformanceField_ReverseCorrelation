function OOD_NOM_Trialwise_compIV_A12(nBasisORI, nBasisSF, isubj, iLocComb, lambda_whiten, flag_regressType, nIter, nJob, iJob)
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

ratio_split = [.6, .3, .1]; % proportion of trials in template set (for RC), training set (for estimating parameters) and test set (for metric predictions)
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
    candidateORI = 3:8;
    candidateSF = 3:8;
    candidateBasisFamilyORI = {'circ_gaussian', 'vonmises'};
    % candidateBasisFamilySF  = {'gaussianLog2', 'asymGaussianLog2', 'LogParabola', 'asymLogParabola'};
    candidateBasisFamilySF  = {'asymGaussianLog2', 'asymLogParabola'};
    candidateRidge  = [0, 10.^(-3:2)];
    % opts.ridge = ridge;
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
    '\n - Energy source: %d (1=TARGET, 2=NOISE)', ...
    '\n - Convolve type: %s', ...
    '\n - IV type: %s\n\n'], ...
    subjName, ...
    iLocComb, namesLocComb{iLocComb}, ...
    nIter, ...
    nblocks, nORI, nSF, ...
    namesType{itype_template}, flag_PatchMode, ...
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

fprintf('%s: ideal template created.\n\n', datetime('now'))

%% A1
% Output file name (before estimation)
iModelA_fit=1;
nameFile_compIV_A1 = sprintf('%s/n%d_J%d_A%d_compIV', nameFolder_NOM_save, nIter, iJob, iModelA_fit);

%% MAIN LOOP over iterations 
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
    switch itype_template
        case 1 % PRS trials only
            useIdx_full = (iPRS_full_rand == 1);
            useIdx_tmpl = (iPRS_tmpl_rand == 1);
        case 2 % ABS trials only
            useIdx_full = (iPRS_full_rand == 0);
            useIdx_tmpl = (iPRS_tmpl_rand == 0);
        otherwise % 3 = both PRS and ABS (here we just keep everything)
            useIdx_full = true(size(iPRS_full_rand));
            useIdx_tmpl = true(size(iPRS_tmpl_rand));
    end

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
    [nTrials_tmpl_cov, nORI, nSF] = size(e3D_tmpl_rand_norm);
    [nTrials_full_cov, ~,   ~  ] = size(e3D_full_rand_norm);

    e3D_tmpl_vec = reshape(e3D_tmpl_rand_norm, [nTrials_tmpl_cov, nORI * nSF]); % trials x channels
    e3D_full_vec = reshape(e3D_full_rand_norm, [nTrials_full_cov, nORI * nSF]); % trials x channels

    % Mean and covariance of channel-energy vectors across trials
    mu_cov_tmpl = mean(e3D_tmpl_vec, 1);
    mu_cov_full = mean(e3D_full_vec, 1);

    Sigma_tmpl = cov(e3D_tmpl_vec - mu_cov_tmpl, 1);   % population covariance
    Sigma_full = cov(e3D_full_vec - mu_cov_full, 1);   % population covariance

    % -------------------------------------------------------------------------
    % 2. Choose covariance-handling strategy
    % -------------------------------------------------------------------------
    if isnan(lambda_whiten)
        % % =====================================================================
        % % Way 2: estimate template in ORIGINAL feature space,
        % %        then apply inverse-covariance correction to the recovered template
        % % =====================================================================
        % % Optional: no shrinkage in inverse-covariance mode
        % Sigma_tmpl_use = Sigma_tmpl;
        % Sigma_full_use = Sigma_full;
        % 
        % % --- TEMPLATE set inverse-covariance correction matrix ---
        % [V_tmpl, D_tmpl] = eig((Sigma_tmpl_use + Sigma_tmpl_use') / 2);
        % d_tmpl = diag(D_tmpl);
        % d_tmpl(d_tmpl < eps_whiten) = eps_whiten;
        % W_inv_tmpl = V_tmpl * diag(1 ./ d_tmpl) * V_tmpl';
        % 
        % % --- FULL set inverse-covariance correction matrix ---
        % [V_full, D_full] = eig((Sigma_full_use + Sigma_full_use') / 2);
        % d_full = diag(D_full);
        % d_full(d_full < eps_whiten) = eps_whiten;
        % W_inv_full = V_full * diag(1 ./ d_full) * V_full';
        % 
        % % Keep original features unchanged
        % e3D_tmpl_use = e3D_tmpl_rand_norm;
        % e3D_full_use = e3D_full_rand_norm;
        % 
        % % Estimate template in original feature space
        % switch flag_regressType
        %     case 1  % Univariate
        %         template_tmpl_raw = SX_sim07_RC(e3D_tmpl_use, resp_tmpl_rand_sel);
        %         template_full_raw = SX_sim07_RC(e3D_full_use, resp_full_rand_sel);
        % 
        %     case 2  % Multivariate + smoothing
        %         out = SX_RC_smoothBasis_circORI_logSF(e3D_tmpl_use, resp_tmpl_rand_sel, axis_tuning{1}, axis_tuning{2}, opts);
        %         template_tmpl_raw = out.template2D;
        % 
        %         out = SX_RC_smoothBasis_circORI_logSF(e3D_full_use, resp_full_rand_sel, axis_tuning{1}, axis_tuning{2}, opts);
        %         template_full_raw = out.template2D;
        % end
        % 
        % % Apply inverse-covariance correction using each dataset's own covariance
        % template_tmpl_raw = invcov_kernel(template_tmpl_raw, W_inv_tmpl, nORI, nSF);
        % template_full_raw = invcov_kernel(template_full_raw, W_inv_full, nORI, nSF);

    else
        % =====================================================================
        % Way 1: whiten the channel-energy predictors before template estimation
        % =====================================================================

        % Shrink covariance toward scaled identity
        Sigma_shrink_tmpl = (1 - lambda_whiten) * Sigma_tmpl + lambda_whiten * mean(diag(Sigma_tmpl)) * eye(size(Sigma_tmpl, 1));
        Sigma_shrink_full = (1 - lambda_whiten) * Sigma_full + lambda_whiten * mean(diag(Sigma_full)) * eye(size(Sigma_full, 1));

        % lambda_whiten = 0 --> use the full empirical covariance
        % lambda_whiten = 1 --> use a scaled identity covariance (ignore channel correlations)

        % --- TEMPLATE set whitening matrix ---
        [V_tmpl, D_tmpl] = eig((Sigma_shrink_tmpl + Sigma_shrink_tmpl') / 2);
        d_tmpl = diag(D_tmpl);
        d_tmpl(d_tmpl < eps_whiten) = eps_whiten;
        W_white_tmpl = V_tmpl * diag(1 ./ sqrt(d_tmpl)) * V_tmpl';

        % --- FULL set whitening matrix ---
        [V_full, D_full] = eig((Sigma_shrink_full + Sigma_shrink_full') / 2);
        d_full = diag(D_full);
        d_full(d_full < eps_whiten) = eps_whiten;
        W_white_full = V_full * diag(1 ./ sqrt(d_full)) * V_full';

        % Apply whitening to each dataset using its own mean and covariance
        e3D_tmpl_use = whiten_e3D(e3D_tmpl_rand_norm, mu_cov_tmpl, W_white_tmpl);
        e3D_full_use = whiten_e3D(e3D_full_rand_norm, mu_cov_full, W_white_full);

        % Estimate template in whitened feature space
        switch flag_regressType
            case 1  % Univariate
                template_tmpl_raw = SX_sim07_RC(e3D_tmpl_use, resp_tmpl_rand_sel);
                template_full_raw = SX_sim07_RC(e3D_full_use, resp_full_rand_sel);

            case 2  % Multivariate + smoothing
                opts.template_ideal = template_ideal;
                out = SX_RC_selectBasis_cv(e3D_tmpl_use, resp_tmpl_rand_sel, axis_tuning{1}, axis_tuning{2}, candidateORI, candidateSF, candidateBasisFamilyORI, candidateBasisFamilySF, candidateRidge, opts);
                template_tmpl_raw = out.template2D;
                nBasisORI_tmpl = out.nBasisORI;
                nBasisSF_tmpl = out.nBasisSF;
                basisFxnORI_tmpl = out.basisFamilyORI;
                basisFxnSF_tmpl = out.basisFamilySF;
                ridge_tmpl = out.ridge;

                out = SX_RC_selectBasis_cv(e3D_full_use, resp_full_rand_sel, axis_tuning{1}, axis_tuning{2}, candidateORI, candidateSF, candidateBasisFamilyORI, candidateBasisFamilySF, candidateRidge, opts);
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

    %% 5. Compute behavioral metrics for all data sets
    %--------------------------------------------%
    metrics_full = fxn_getMetrics(resp_full_rand, iPRS_full_rand, respC_full_rand, cst_full_rand, RT_full_rand);
    metrics_tmpl = fxn_getMetrics(resp_tmpl_rand, iPRS_tmpl_rand, respC_tmpl_rand, cst_tmpl_rand, RT_tmpl_rand);
    metrics_train = fxn_getMetrics(resp_train_rand, iPRS_train_rand, respC_train_rand, cst_train_rand, RT_train_rand);
    metrics_test = fxn_getMetrics(resp_test_rand, iPRS_test_rand, respC_test_rand, cst_test_rand, RT_test_rand);
    %--------------------------------------------%

    % metrics: [dprime, criterion, [pC, pHit, pFA], nanmean(respC), pYES, mean(1./contrast), median(RT)];
    metrics = [metrics_full; metrics_tmpl; metrics_train; metrics_test];

    % Save criterion in z units; later used to convert to criterion in IV
    criterion_z_train = metrics_train(2);
    criterion_z_test = metrics_test(2);

    clear data; % size of IV differs across subjects and  iterations

    %% 6. Compute decision variable (DV) from template and energy (TEST set)
    flag_permT=0;
    %---------------%
    DV_train = fxn_getIV_v3(e3D_train_rand, template_notNormed, convolveType, IVType, flag_permT, ORI_bound);
    DV_test = fxn_getIV_v3(e3D_test_rand, template_notNormed, convolveType, IVType, flag_permT, ORI_bound);
    %---------------%
    nData_train = length(DV_train);
    nData_test = length(DV_test);

    %% 7. Bin IV values
    % Separately for PRS and ABS trials
    [nTrials_PRS_allBins_train, ~, iTrial4Bin_PRS_train] = histcounts(DV_train(iPRS_train_rand == 1), nBins);
    [nTrials_ABS_allBins_train, ~, iTrial4Bin_ABS_train] = histcounts(DV_train(iPRS_train_rand == 0), nBins);
    [nTrials_PRS_allBins_test, ~, iTrial4Bin_PRS_test] = histcounts(DV_test(iPRS_test_rand == 1), nBins);
    [nTrials_ABS_allBins_test, ~, iTrial4Bin_ABS_test] = histcounts(DV_test(iPRS_test_rand == 0), nBins);

    % For all trials combined
    [nTrials_allBins_train, ~, iTrial4Bin_train] = histcounts(DV_train, nBins);
    [nTrials_allBins_test, ~, iTrial4Bin_test] = histcounts(DV_test, nBins);

    %% 8. Compile data struct for this  iteration
    % Store data_train
    data_train = struct();
    data_train.criterion_z = criterion_z_train;
    data_train.ndata = nData_train;
    data_train.IV = DV_train;
    data_train.nTrials_allBins = nTrials_allBins_train;
    data_train.iTrial4Bin = iTrial4Bin_train;
    data_train.nTrials_PRS_allBins = nTrials_PRS_allBins_train;
    data_train.iTrial4Bin_PRS = iTrial4Bin_PRS_train;
    data_train.nTrials_ABS_allBins = nTrials_ABS_allBins_train;
    data_train.iTrial4Bin_ABS = iTrial4Bin_ABS_train;

    data_train.iPRS = iPRS_train_rand;
    data_train.resp = resp_train_rand;
    data_train.cst = cst_train_rand;
    data_train.iPair = iPair_train_rand;
    data_train.respC = respC_train_rand;
    data_train.RT = RT_train_rand;
    data_train.metrics_sim = metrics_train; % keep field name 'metrics_sim'
    data_train_allIter{iIter} = data_train;

    % Store data_test
    data_test = struct();
    data_test.criterion_z = criterion_z_test;
    data_test.ndata = nData_test;
    data_test.IV = DV_test;
    data_test.nTrials_allBins = nTrials_allBins_test;
    data_test.iTrial4Bin = iTrial4Bin_test;
    data_test.nTrials_PRS_allBins = nTrials_PRS_allBins_test;
    data_test.iTrial4Bin_PRS = iTrial4Bin_PRS_test;
    data_test.nTrials_ABS_allBins = nTrials_ABS_allBins_test;
    data_test.iTrial4Bin_ABS = iTrial4Bin_ABS_test;

    data_test.iPRS = iPRS_test_rand;
    data_test.resp = resp_test_rand;
    data_test.cst = cst_test_rand;
    data_test.iPair = iPair_test_rand;
    data_test.respC = respC_test_rand;
    data_test.RT = RT_test_rand;
    data_test.metrics_sim = metrics_test; % keep field name 'metrics_sim'
    % Store
    data_test_allIter{iIter} = data_test;

    % Store metrics of three datasets together (for fast plotting in NOMplot_compIV)
    data_metrics_allIter(iIter, :, :) = metrics;

    % Store templates
    if iModelA_fit==2
        template_notNormed_tmpl = template_notNormed;
        template_notNormed_full = template_notNormed;
    end
    template_tmpl_allIter(iIter, :, :) = template_notNormed_tmpl;
    template_full_allIter(iIter, :, :) = template_notNormed_full;
    nBasisORI_tmpl_allIter(iIter, :) = nBasisORI_tmpl;
    nBasisSF_tmpl_allIter(iIter, :) = nBasisSF_tmpl;
    basisFxnORI_tmpl_allIter{iIter} = basisFxnORI_tmpl;
    basisFxnSF_tmpl_allIter{iIter} = basisFxnSF_tmpl;
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

%% SAVE 
save(nameFile_compIV_A1, 'time_progress', 'template_ideal', 'criterion_z*', '*_allIter', 'names*', 'flag*', 'ratio_split', 'ORI_bound', '*Type');

% Delete the progress report
if exist(nameFile_progress_new, 'file')
    delete(nameFile_progress_new);
end

fprintf('\n\n%s: A1 Outputs saved.\n\n', datetime('now'))

%% Plot (optional) 
if flag_plot_compIV
    %-------------------%
    NOMplot_compIV;
    %-------------------%
    fprintf('\n\n%s: A1 Plots created.\n\n', datetime('now'))
end
close all;

%% End timing 
time_end = datetime('now');
fprintf('%s: A1 Compute DV done.\n\n', time_end)
elapsed = time_end - time_start;
fprintf('A1 Time used: %s\n\n\n\n', char(elapsed));

%% A2
% Output file name (before estimation)
iModelA_fit=2;
nameFile_compIV_A2 = sprintf('%s/n%d_J%d_A%d_compIV', nameFolder_NOM_save, nIter, iJob, iModelA_fit);

%% Process the true template
% Average the derived template (tmpl set) across iterations
template_tmpl_ave_A1 = squeeze(mean(template_tmpl_allIter, 1, 'omitnan'));
peak_tmpl_A1 = max(abs(template_tmpl_ave_A1(:)));
% Rescale the true/ideak template to match the peak of the derived template
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

%% MAIN LOOP over iterations 
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
    switch itype_template
        case 1 % PRS trials only
            useIdx_full = (iPRS_full_rand == 1);
            useIdx_tmpl = (iPRS_tmpl_rand == 1);
        case 2 % ABS trials only
            useIdx_full = (iPRS_full_rand == 0);
            useIdx_tmpl = (iPRS_tmpl_rand == 0);
        otherwise % 3 = both PRS and ABS (here we just keep everything)
            useIdx_full = true(size(iPRS_full_rand));
            useIdx_tmpl = true(size(iPRS_tmpl_rand));
    end

    %% 5. Compute behavioral metrics for all data sets
    %--------------------------------------------%
    metrics_full = fxn_getMetrics(resp_full_rand, iPRS_full_rand, respC_full_rand, cst_full_rand, RT_full_rand);
    metrics_tmpl = fxn_getMetrics(resp_tmpl_rand, iPRS_tmpl_rand, respC_tmpl_rand, cst_tmpl_rand, RT_tmpl_rand);
    metrics_train = fxn_getMetrics(resp_train_rand, iPRS_train_rand, respC_train_rand, cst_train_rand, RT_train_rand);
    metrics_test = fxn_getMetrics(resp_test_rand, iPRS_test_rand, respC_test_rand, cst_test_rand, RT_test_rand);
    %--------------------------------------------%

    % metrics: [dprime, criterion, [pC, pHit, pFA], nanmean(respC), pYES, mean(1./contrast), median(RT)];
    metrics = [metrics_full; metrics_tmpl; metrics_train; metrics_test];

    % Save criterion in z units; later used to convert to criterion in IV
    criterion_z_train = metrics_train(2);
    criterion_z_test = metrics_test(2);

    clear data; % size of IV differs across subjects and  iterations

    %% 6. Compute decision variable (DV) from template and energy (TEST set)
    flag_permT=0;

    %---------------%
    DV_train = fxn_getIV_v3(e3D_train_rand, template_notNormed, convolveType, IVType, flag_permT, ORI_bound);
    DV_test = fxn_getIV_v3(e3D_test_rand, template_notNormed, convolveType, IVType, flag_permT, ORI_bound);
    %---------------%
    nData_train = length(DV_train);
    nData_test = length(DV_test);

    %% 7. Bin IV values
    % Separately for PRS and ABS trials
    [nTrials_PRS_allBins_train, ~, iTrial4Bin_PRS_train] = histcounts(DV_train(iPRS_train_rand == 1), nBins);
    [nTrials_ABS_allBins_train, ~, iTrial4Bin_ABS_train] = histcounts(DV_train(iPRS_train_rand == 0), nBins);
    [nTrials_PRS_allBins_test, ~, iTrial4Bin_PRS_test] = histcounts(DV_test(iPRS_test_rand == 1), nBins);
    [nTrials_ABS_allBins_test, ~, iTrial4Bin_ABS_test] = histcounts(DV_test(iPRS_test_rand == 0), nBins);

    % For all trials combined
    [nTrials_allBins_train, ~, iTrial4Bin_train] = histcounts(DV_train, nBins);
    [nTrials_allBins_test, ~, iTrial4Bin_test] = histcounts(DV_test, nBins);

    %% 8. Compile data struct for this  iteration
    % Store data_train
    data_train = struct();
    data_train.criterion_z = criterion_z_train;
    data_train.ndata = nData_train;
    data_train.IV = DV_train;
    data_train.nTrials_allBins = nTrials_allBins_train;
    data_train.iTrial4Bin = iTrial4Bin_train;
    data_train.nTrials_PRS_allBins = nTrials_PRS_allBins_train;
    data_train.iTrial4Bin_PRS = iTrial4Bin_PRS_train;
    data_train.nTrials_ABS_allBins = nTrials_ABS_allBins_train;
    data_train.iTrial4Bin_ABS = iTrial4Bin_ABS_train;

    data_train.iPRS = iPRS_train_rand;
    data_train.resp = resp_train_rand;
    data_train.cst = cst_train_rand;
    data_train.iPair = iPair_train_rand;
    data_train.respC = respC_train_rand;
    data_train.RT = RT_train_rand;
    data_train.metrics_sim = metrics_train; % keep field name 'metrics_sim'
    data_train_allIter{iIter} = data_train;

    % Store data_test
    data_test = struct();
    data_test.criterion_z = criterion_z_test;
    data_test.ndata = nData_test;
    data_test.IV = DV_test;
    data_test.nTrials_allBins = nTrials_allBins_test;
    data_test.iTrial4Bin = iTrial4Bin_test;
    data_test.nTrials_PRS_allBins = nTrials_PRS_allBins_test;
    data_test.iTrial4Bin_PRS = iTrial4Bin_PRS_test;
    data_test.nTrials_ABS_allBins = nTrials_ABS_allBins_test;
    data_test.iTrial4Bin_ABS = iTrial4Bin_ABS_test;

    data_test.iPRS = iPRS_test_rand;
    data_test.resp = resp_test_rand;
    data_test.cst = cst_test_rand;
    data_test.iPair = iPair_test_rand;
    data_test.respC = respC_test_rand;
    data_test.RT = RT_test_rand;
    data_test.metrics_sim = metrics_test; % keep field name 'metrics_sim'
    % Store
    data_test_allIter{iIter} = data_test;

    % Store metrics of three datasets together (for fast plotting in NOMplot_compIV)
    data_metrics_allIter(iIter, :, :) = metrics;

    % Store templates
    if iModelA_fit==2
        template_notNormed_tmpl = template_notNormed;
        template_notNormed_full = template_notNormed;
    end
    template_tmpl_allIter(iIter, :, :) = template_notNormed_tmpl;
    template_full_allIter(iIter, :, :) = template_notNormed_full;

    %% Save a progress report in the folder to indicate the finished iteration and time spent
    time_progress = datetime('now');
    time_progress = ceil(minutes(time_progress-time_start)); % round up to minutes
    save(nameFile_progress, 'iIter')

    % rename (to avoid saving one file for each iteration)
    nameFile_progress_new = sprintf('%s_%d_%dmin.mat', nameFile_compIV_A2, iIter, time_progress);
    movefile(nameFile_progress, nameFile_progress_new);
    nameFile_progress = nameFile_progress_new;

end % end for iIter

fprintf('\n\n%s: A2 All iterations done.\n\n', datetime('now'))

%% SAVE 
save(nameFile_compIV_A2, 'time_progress', 'template_ideal', 'criterion_z*', '*_allIter', 'names*', 'flag*', 'ratio_split', 'ORI_bound', '*Type');

% Delete the progress report
if exist(nameFile_progress_new, 'file')
    delete(nameFile_progress_new);
end

fprintf('\n\n%s: A2 Outputs saved.\n\n', datetime('now'))
end

%% HELPER
function e3D_white = whiten_e3D(e3D_in, mu_cov, W_white)
[nT, nOri, nSf] = size(e3D_in);
X = reshape(e3D_in, [nT, nOri * nSf]);
Xw = (X - mu_cov) * W_white;
e3D_white = reshape(Xw, [nT, nOri, nSf]);
end

function kernel_corr = invcov_kernel(kernel_raw, W_inv, nOri, nSf)
k_raw = kernel_raw(:);
k_corr = W_inv * k_raw;
kernel_corr = reshape(k_corr, [nOri, nSf]);
end