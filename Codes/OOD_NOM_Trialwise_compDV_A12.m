function OOD_NOM_Trialwise_compDV_A12(isubj, iLocComb, lambda_whiten, nIter, nJob, iJob)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OOD_NOM_Trialwise_compDV.m
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

% Add paths for custom functions (client, using absolute paths)
addpath(genpath('Codes/'));

%% Global settings
% Define directories
nORI=19;

% --------------%
SX_RC1_setting; % defines nameFolder_*, nORI, nSF, namesLocComb, namesModelA, etc.
% --------------%
flag_regressType = flag_regressType;

% Make tuning axes explicit variables (avoid brace-index parsing issues).
nORI=nORI;
nSF=nSF;
axisORI = filtersOri_all - 90;
axisSF = filtersSF_all_log;

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
DVType = 1; % 1=sum of the dot product/convolution; 2=max; 3=normalized
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
ORI_bound = [5, 14]; % orientation window, passed to fxn_getDV_v3

% whiten features
eps_whiten = 1e-3;  % floor for eigenvalues

% Setting for multivariate regression with smoothing
% 1=Univariate; 2=MultiSmooth and Univariate
if flag_regressType == 2
    basisCfg = SX_RC_getBasisSettings();

    opts = struct();
    candidateORI = basisCfg.candidateORI;
    candidateSF = basisCfg.candidateSF;
    candidateBasisFamilyORI = {basisCfg.basisFamilyORI};
    candidateBasisFamilySF = {basisCfg.basisFamilySF};
    candidateRidge = basisCfg.candidateRidge;
    opts.basisWidthScaleORI = basisCfg.basisWidthScaleORI; % set vector to search ORI width scales; typical range: [0.5, 0.9]; higher values = wider basis functions = stronger smoothing
    opts.basisWidthScaleSF = basisCfg.basisWidthScaleSF; % set vector to search SF width scales; typical range: [0.5, 0.9]
    opts.asymSF_rightLeftRatio = basisCfg.asymSF_rightLeftRatio; % typical range: [1.1, 1.5]
    opts.nFolds = basisCfg.nFolds; % number of folds for cross-validation (CV) to select the best model
    opts.link = basisCfg.link;
    opts.sigmaORI_deg = basisCfg.sigmaORI_deg;
    opts.sigmaSF_log2 = basisCfg.sigmaSF_log2;
    opts.oriPeriod_deg = basisCfg.oriPeriod_deg;
    opts.zscorePredictor = basisCfg.zscorePredictor;
    opts.maxIter = basisCfg.maxIter;
    opts.tol = basisCfg.tol;
end

% Settings for fitting tuning functions
nRep = 20;
iFamily_ORI = 10; % input to predKernelSF, 1=scaled gaussian, 8=DoG, 10=von Mises,
iFamily_SF = 14; % 2=log parabola, 14=asymmetric Gaussian in log2 space,
problem_setting = MultiStart('StartPointsToRun', 'bounds','UseParallel', 1, 'Display', 'off');

% Local fmincon options used by runMultistartFmincon (parfor-safe path).
if exist('optimoptions', 'file')
    options_fmin = optimoptions('fmincon', 'Display', 'off');
else
    options_fmin = optimset('Display', 'off');
end

% Cache bounds as plain numeric vectors (avoid brace indexing in loop expressions).
ubORI = ub_full_all{iFamily_ORI};
lbORI = lb_full_all{iFamily_ORI};
ubSF = ub_full_all{iFamily_SF};
lbSF = lb_full_all{iFamily_SF};

flag_plot_tuning = 0;

iSess_start = 1; % first session included
convolveType = 1; % DV from 1=cross-correlation; 2=convolution (fxn_getDV_v3)
flag_standEnergy = 1; % 1=z-score energy before RC
flag_plot_compDV = 1; % plot DV distributions and kernels at the end
if strcmp(str_envir,'HPC'), flag_plot_compDV = 0; end % don't plot when running on HPC

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

if ~strcmp(str_envir, 'HPC') && isempty(dir(nameFolder_Figures_perSubj))
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
    '\n - DV type: %s\n\n'], ...
    subjName, ...
    iLocComb, namesLocComb{iLocComb}, ...
    nIter, ...
    nblocks, nORI, nSF, ...
    namesType{itype_template}, wABS, wPRS, weightScale_RC, flag_PatchMode, ...
    namesConvolveType{convolveType}, namesDVType{DVType});

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

%% [A1] Setup
% Output file name (before estimation)
iModelA_fit = 1;
nameFile_compDV_A1 = sprintf('%s/n%d_J%d_A%d_compDV', nameFolder_NOM_save, nIter, iJob, iModelA_fit);

%% [A1] Main loop over iterations
% fprintf('%s: Creating empty placeholders for running iterations.\n\n', datetime('now'))
data_metrics_allIter = nan(nIter, nDatasets_full, nMetrics); % see fxn_getMetrics for all 11 metrics
data_train_allIter = cell(nIter, 1); % tmpl-split template (mainstream output)
data_test_allIter  = cell(nIter, 1);
template_tmpl_allIter = nan(nIter, nORI, nSF);
template_full_allIter = template_tmpl_allIter;
nBasisORI_tmpl_allIter = nan(nIter, 1);
nBasisSF_tmpl_allIter = nan(nIter, 1);
basisFxnORI_tmpl_allIter = cell(nIter, 1);
basisFxnSF_tmpl_allIter  = cell(nIter, 1);
ridge_tmpl_allIter = nBasisSF_tmpl_allIter;
basisWidthScaleORI_tmpl_allIter = nBasisSF_tmpl_allIter;
basisWidthScaleSF_tmpl_allIter  = nBasisSF_tmpl_allIter;
asymSF_rightLeftRatio_tmpl_allIter = nBasisSF_tmpl_allIter;

sep_allIter = nan(nIter, 2); % 1=template set, 2=full set
margORI_allIter = nan(nIter, 2, nORI);
margPred_ORI_allIter = nan(nIter, 2, nORI);
margParams_ORI_allIter = nan(nIter, 2, length(namesParams_all{iFamily_ORI}));
margR2_ORI_allIter = nan(nIter, 2);
margSF_allIter = nan(nIter, 2, nSF);
margPred_SF_allIter = nan(nIter, 2, nSF);
margParams_SF_allIter = nan(nIter, 2, length(namesParams_all{iFamily_SF}));
margR2_SF_allIter = nan(nIter, 2);

fprintf('%s: A1 Started running %d iterations.\n\n', datetime('now'), nIter)

parfor iIter = 1:nIter
    fprintf('\n%s: %d...', datetime('now'), iIter);

    % Deterministic randomness for THIS iteration (global index across jobs)
    iterStream = RandStream('Threefry', 'Seed', S_seed.grandSeed);
    iterStream.Substream = S_seed.iterIdxList(iIter);
    RandStream.setGlobalStream(iterStream);

    %% [A1] 1. Resample trials into FULL or TEMPLATE/TRAIN / TEST
    R = fxn_resampleTrials(dataMatrix, e3D_allT, iSess_select, iLoc_all, ratio_split);
    e3D_full_rand = R.e3D_full_rand;
    iPRS_full_rand = R.iPRS_full_rand;
    resp_full_rand = R.resp_full_rand;
    cst_full_rand = R.cst_full_rand;
    respC_full_rand = R.respC_full_rand;
    iPair_full_rand = R.iPair_full_rand;
    RT_full_rand = R.RT_full_rand;

    e3D_tmpl_rand = R.e3D_tmpl_rand;
    iPRS_tmpl_rand = R.iPRS_tmpl_rand;
    resp_tmpl_rand = R.resp_tmpl_rand;
    cst_tmpl_rand = R.cst_tmpl_rand;
    respC_tmpl_rand = R.respC_tmpl_rand;
    iPair_tmpl_rand = R.iPair_tmpl_rand;
    RT_tmpl_rand = R.RT_tmpl_rand;

    e3D_train_rand = R.e3D_train_rand;
    iPRS_train_rand = R.iPRS_train_rand;
    resp_train_rand = R.resp_train_rand;
    cst_train_rand = R.cst_train_rand;
    respC_train_rand = R.respC_train_rand;
    iPair_train_rand = R.iPair_train_rand;
    RT_train_rand = R.RT_train_rand;

    e3D_test_rand = R.e3D_test_rand;
    iPRS_test_rand = R.iPRS_test_rand;
    resp_test_rand = R.resp_test_rand;
    cst_test_rand = R.cst_test_rand;
    respC_test_rand = R.respC_test_rand;
    iPair_test_rand = R.iPair_test_rand;
    RT_test_rand = R.RT_test_rand;

    %% [A1] 2. Select which trials to use to estimate the template (TEMPLATE set)
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

    %% [A1] 2b. Build fixed transform from ABS template trials only
    % Fixed transform per iteration:
    % raw energy -> contrast-specific ABS z-score -> optional whitening.
    Tfix = buildFixedTransformFromABS(e3D_tmpl_rand, cst_tmpl_rand, iPRS_tmpl_rand, lambda_whiten, eps_whiten);

    % Apply the SAME transform to all relevant sets.
    e3D_tmpl_forRC = applyFixedTransform(e3D_tmpl_rand_sel, cst_tmpl_rand_sel, Tfix);
    e3D_full_forRC = applyFixedTransform(e3D_full_rand_sel, cst_full_rand_sel, Tfix);
    e3D_train_forDV = applyFixedTransform(e3D_train_rand, cst_train_rand, Tfix);
    e3D_test_forDV = applyFixedTransform(e3D_test_rand, cst_test_rand, Tfix);

    %% [A1] 3. Estimate template from transformed energy

    % Safe defaults for basis function outputs (case 2 only)
    nBasisORI_tmpl = nan;
    nBasisSF_tmpl = nan;
    basisFxnORI_tmpl = '';
    basisFxnSF_tmpl = '';
    basisWidthScaleORI_tmpl = nan;
    basisWidthScaleSF_tmpl = nan;
    asymSF_rightLeftRatio_tmpl = nan;
    ridge_tmpl = nan;

    % Apply class weighting to z-scored energy (if wABS>.5, ABS is emphasized)
    [e3D_tmpl_forRC, resp_tmpl_forRC] = applyClassWeightsForRC(e3D_tmpl_forRC, resp_tmpl_rand_sel, iPRS_tmpl_rand_sel, wABS, wPRS, weightScale_RC);
    [e3D_full_forRC, resp_full_forRC] = applyClassWeightsForRC(e3D_full_forRC, resp_full_rand_sel, iPRS_full_rand_sel, wABS, wPRS, weightScale_RC);

    % Estimate template from raw energy via regression
    switch flag_regressType
        case 1  % Univariate
            template_tmpl_raw = SX_sim07_RC(e3D_tmpl_forRC, resp_tmpl_forRC);
            template_full_raw = SX_sim07_RC(e3D_full_forRC, resp_full_forRC);

        case 2  % Multivariate + smoothing
            opts_local = opts;
            opts_local.template_ideal = template_ideal;
            out = SX_RC_selectBasis_cv(e3D_tmpl_forRC, resp_tmpl_forRC, axisORI, axisSF, candidateORI, candidateSF, candidateBasisFamilyORI, candidateBasisFamilySF, candidateRidge, opts_local);
            template_tmpl_raw = out.template2D;

            nBasisORI_tmpl = out.nBasisORI;
            nBasisSF_tmpl = out.nBasisSF;
            basisFxnORI_tmpl = out.basisFamilyORI;
            basisFxnSF_tmpl = out.basisFamilySF;
            basisWidthScaleORI_tmpl = out.basisWidthScaleORI;
            basisWidthScaleSF_tmpl = out.basisWidthScaleSF;
            asymSF_rightLeftRatio_tmpl = out.asymSF_rightLeftRatio;
            ridge_tmpl = out.ridge;

            out = SX_RC_selectBasis_cv(e3D_full_forRC, resp_full_forRC, axisORI, axisSF, candidateORI, candidateSF, candidateBasisFamilyORI, candidateBasisFamilySF, candidateRidge, opts_local);
            template_full_raw = out.template2D;

    end % switch

    %% [A1] 4. Regularize the derived template
    template_full_raw = fxn_getTemplate(template_full_raw, templateType, flag_plot_template);
    template_tmpl_raw = fxn_getTemplate(template_tmpl_raw, templateType, flag_plot_template);

    %% [A1] Check how channel correlation affects template estimation
    % -----PCA-based template recovery =====
    if flag_plot_compDV
        % NOMplot_PCA
    end

    % ----- Debug: check whether channels are too correlated thus redundant
    % !!!! Generate one figure per iteration !!!!
    if flag_plot_compDV
        % NOMplot_checkChannelCorr
    end

    % ----- Compare the univariate regression and multivariation regression with smoothing
    % One of the two may not be generated; need to pause and generate manually
    if flag_plot_compDV
        % template_smooth = template_full_raw;
        % template_noSmoothing = SX_sim07_RC(e3D_full_rand_norm, resp_full_rand_sel);
        % RCplot_compTempSmoothing
    end

    %% [A1] 5. Marginalization, fit tuning functions, and calculate DV

    % Build behavioral-data structs once (independent of template).
    d_full  = struct('resp', resp_full_rand,  'iPRS', iPRS_full_rand,  'respC', respC_full_rand,  'cst', cst_full_rand,  'RT', RT_full_rand);
    d_tmpl  = struct('resp', resp_tmpl_rand,  'iPRS', iPRS_tmpl_rand,  'respC', respC_tmpl_rand,  'cst', cst_tmpl_rand,  'RT', RT_tmpl_rand);
    d_train = struct('resp', resp_train_rand, 'iPRS', iPRS_train_rand, 'respC', respC_train_rand, 'cst', cst_train_rand, 'RT', RT_train_rand, 'iPair', iPair_train_rand);
    d_test  = struct('resp', resp_test_rand,  'iPRS', iPRS_test_rand,  'respC', respC_test_rand,  'cst', cst_test_rand,  'RT', RT_test_rand,  'iPair', iPair_test_rand);

    for iDataset = 1:2 % 1=template set; 2=full set
        DV_train_override = [];
        DV_test_override = [];

        switch iDataset
            case 1
                template_use_raw = template_tmpl_raw;
            case 2
                template_use_raw = template_full_raw;
        end

        % Calculate separability
        templateType_recon = 2; % 2=reconstructed template
        %--------------------------------------------%
        template_recon = fxn_getTemplate(template_use_raw, templateType_recon, flag_plot_template);
        %--------------------------------------------%
        sep = corr(template_use_raw(:), template_recon(:));

        % Marginalize the template into 1D ORI and SF profiles
        margORI = mean(template_use_raw, 2)';
        margSF = mean(template_use_raw, 1);

        % Fitting
        if any(isnan(margORI)), error('ALERT: NaN in margORI!'), end

        xORI = axisORI;
        fxn_tuningLoss_ORI = @(param_est) sum((margORI - predSFkernel(xORI, iFamily_ORI, param_est, flag_plot_tuning)).^2);
        margParams_ORI = runMultistartFmincon(fxn_tuningLoss_ORI, (ubORI+lbORI)/2, lbORI, ubORI, nRep, options_fmin);
        %--------------------------------------------%
        margPred_ORI = predSFkernel(xORI, iFamily_ORI, margParams_ORI, flag_plot_tuning);
        %--------------------------------------------%
        margR2_ORI = 1-sumsqr(margORI-margPred_ORI)/sumsqr(margORI-mean(margORI));

        xSF = axisSF;
        xSF_ln = 2.^xSF; % fit in linear SF space
        fxn_tuningLoss_SF = @(param_est) sum((margSF - predSFkernel(xSF_ln, iFamily_SF, param_est, flag_plot_tuning)).^2);
        margParams_SF = runMultistartFmincon(fxn_tuningLoss_SF, (ubSF+lbSF)/2, lbSF, ubSF, nRep, options_fmin);
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

        %% [A1] 6. Metrics, DV, binning, data structs (tmpl split only)
        if iDataset == 1
            if flag_regressType == 2
                basisOpts_use = struct();
                basisOpts_use.nBasisORI = nBasisORI_tmpl;
                basisOpts_use.nBasisSF = nBasisSF_tmpl;
                basisOpts_use.basisFamilyORI = basisFxnORI_tmpl;
                basisOpts_use.basisFamilySF = basisFxnSF_tmpl;
                basisOpts_use.basisWidthScaleORI = basisWidthScaleORI_tmpl;
                basisOpts_use.basisWidthScaleSF = basisWidthScaleSF_tmpl;
                basisOpts_use.asymSF_rightLeftRatio = asymSF_rightLeftRatio_tmpl;
                basisOpts_use.oriPeriod_deg = opts.oriPeriod_deg;

                Z_train = SX_RC_basisProject(e3D_train_forDV, axisORI, axisSF, basisOpts_use);
                Z_test = SX_RC_basisProject(e3D_test_forDV, axisORI, axisSF, basisOpts_use);
                beta_basis = SX_RC_basisProject(reshape(template_use_raw, [1, nORI, nSF]), axisORI, axisSF, basisOpts_use);
                DV_train_override = Z_train * beta_basis(:);
                DV_test_override = Z_test * beta_basis(:);
            end

            [data_train, data_test, metrics] = fxn_compDV_packData( ...
                d_full, d_tmpl, d_train, d_test, ...
                e3D_train_forDV, e3D_test_forDV, ...
                template_use_raw, convolveType, DVType, ORI_bound, nBins, ...
                DV_train_override, DV_test_override);
            data_train_allIter{iIter} = data_train;
            data_test_allIter{iIter}  = data_test;
            data_metrics_allIter(iIter, :, :) = metrics;
        end
    end % iDataset

    %% [A1] Store template-related values
    template_tmpl_allIter(iIter, :, :) = template_tmpl_raw;
    template_full_allIter(iIter, :, :) = template_full_raw;
    nBasisORI_tmpl_allIter(iIter, :) = nBasisORI_tmpl;
    nBasisSF_tmpl_allIter(iIter, :) = nBasisSF_tmpl;
    basisFxnORI_tmpl_allIter{iIter} = basisFxnORI_tmpl;
    basisFxnSF_tmpl_allIter{iIter} = basisFxnSF_tmpl;
    basisWidthScaleORI_tmpl_allIter(iIter) = basisWidthScaleORI_tmpl;
    basisWidthScaleSF_tmpl_allIter(iIter) = basisWidthScaleSF_tmpl;
    asymSF_rightLeftRatio_tmpl_allIter(iIter) = asymSF_rightLeftRatio_tmpl;
    ridge_tmpl_allIter(iIter) = ridge_tmpl;

end % end parfor iIter

fprintf('\n\n%s: A1 All iterations done.\n\n', datetime('now'))

time_progress = ceil(minutes(datetime('now')-time_start)); % round up to minutes

%% [A1] Save
save(nameFile_compDV_A1, 'time_progress', 'template_ideal', '*_allIter', 'names*', 'flag*', 'ratio_split', 'ORI_bound', '*Type');

fprintf('\n\n%s: A1 Outputs saved.\n\n', datetime('now'))

%% [A1] Plot
if flag_plot_compDV
    %-------------------%
    NOMplot_compDV;
    %-------------------%
    fprintf('\n\n%s: A1 Plots created.\n\n', datetime('now'))
end
close all;

%% [A1] End timing
time_end = datetime('now');
fprintf('%s: A1 Compute DV done.\n\n', time_end)
elapsed = time_end - time_start;
fprintf('A1 Time used: %s\n\n\n\n', char(elapsed));

%% [A2] Setup
iModelA_fit=2;
nameFile_compDV_A2 = sprintf('%s/n%d_J%d_A%d_compDV', nameFolder_NOM_save, nIter, iJob, iModelA_fit);

%% [A2] Process the true template
% Average the derived template (tmpl set) across iterations
template_tmpl_ave_A1 = squeeze(mean(template_tmpl_allIter, 1, 'omitnan'));
peak_tmpl_A1 = max(abs(template_tmpl_ave_A1(:)));
% Rescale the true/ideal template to match the peak of the derived template
% (apply scaling right after template_notNormed is loaded below)
if ~isnumeric(isubj) % for IO, you may reload 'template_true' from disk
    load(sprintf('%s/truth.mat', nameFolder_OOD_load), 'template_true');
    template_use = template_true;
else
    template_use = template_ideal;
end

if isfinite(peak_tmpl_A1) && peak_tmpl_A1 > 0
    peak_true = max(abs(template_use(:)));
    if isfinite(peak_true) && peak_true > 0
        template_use = template_use * (peak_tmpl_A1 / peak_true);
    end
end

%% [A2] Main loop over iterations
% fprintf('%s: Creating empty placeholders for running iterations.\n\n', datetime('now'))
data_metrics_allIter = nan(nIter, nDatasets_full, nMetrics); % see fxn_getMetrics for all 11 metrics
data_train_allIter = cell(nIter, 1);
data_test_allIter  = cell(nIter, 1);
template_tmpl_allIter = nan(nIter, nORI, nSF);
template_full_allIter = template_tmpl_allIter;
fprintf('%s: A2 Started running %d iterations.\n\n', datetime('now'), nIter)

parfor iIter = 1:nIter
    fprintf('\n%s: %d...', datetime('now'), iIter);

    % Deterministic randomness for THIS iteration (global index across jobs)
    iterStream = RandStream('Threefry', 'Seed', S_seed.grandSeed);
    iterStream.Substream = S_seed.iterIdxList(iIter);
    RandStream.setGlobalStream(iterStream);

    %% [A2] 1. Resample trials into FULL or TEMPLATE/TRAIN / TEST
    R = fxn_resampleTrials(dataMatrix, e3D_allT, iSess_select, iLoc_all, ratio_split);
    e3D_full_rand = R.e3D_full_rand;
    iPRS_full_rand = R.iPRS_full_rand;
    resp_full_rand = R.resp_full_rand;
    cst_full_rand = R.cst_full_rand;
    respC_full_rand = R.respC_full_rand;
    iPair_full_rand = R.iPair_full_rand;
    RT_full_rand = R.RT_full_rand;

    e3D_tmpl_rand = R.e3D_tmpl_rand;
    iPRS_tmpl_rand = R.iPRS_tmpl_rand;
    resp_tmpl_rand = R.resp_tmpl_rand;
    cst_tmpl_rand = R.cst_tmpl_rand;
    respC_tmpl_rand = R.respC_tmpl_rand;
    iPair_tmpl_rand = R.iPair_tmpl_rand;
    RT_tmpl_rand = R.RT_tmpl_rand;

    e3D_train_rand = R.e3D_train_rand;
    iPRS_train_rand = R.iPRS_train_rand;
    resp_train_rand = R.resp_train_rand;
    cst_train_rand = R.cst_train_rand;
    respC_train_rand = R.respC_train_rand;
    iPair_train_rand = R.iPair_train_rand;
    RT_train_rand = R.RT_train_rand;

    e3D_test_rand = R.e3D_test_rand;
    iPRS_test_rand = R.iPRS_test_rand;
    resp_test_rand = R.resp_test_rand;
    cst_test_rand = R.cst_test_rand;
    respC_test_rand = R.respC_test_rand;
    iPair_test_rand = R.iPair_test_rand;
    RT_test_rand = R.RT_test_rand;

    %% [A2] 2. Build fixed transform from ABS template trials only
    [useIdx_full, useIdx_tmpl] = selectTemplateIndices(itype_template, iPRS_full_rand, iPRS_tmpl_rand);

    e3D_tmpl_rand_sel = e3D_tmpl_rand(useIdx_tmpl, :, :);
    cst_tmpl_rand_sel = cst_tmpl_rand(useIdx_tmpl);
    iPRS_tmpl_rand_sel = iPRS_tmpl_rand(useIdx_tmpl);

    e3D_full_rand_sel = e3D_full_rand(useIdx_full, :, :);
    cst_full_rand_sel = cst_full_rand(useIdx_full);
    iPRS_full_rand_sel = iPRS_full_rand(useIdx_full);

    Tfix = buildFixedTransformFromABS(e3D_tmpl_rand, cst_tmpl_rand, iPRS_tmpl_rand, lambda_whiten, eps_whiten);
    e3D_train_forDV = applyFixedTransform(e3D_train_rand, cst_train_rand, Tfix);
    e3D_test_forDV = applyFixedTransform(e3D_test_rand, cst_test_rand, Tfix);

    %% [A2] 5-8. Metrics, DV, binning (using centered raw energy)
    d_full  = struct('resp', resp_full_rand,  'iPRS', iPRS_full_rand,  'respC', respC_full_rand,  'cst', cst_full_rand,  'RT', RT_full_rand);
    d_tmpl  = struct('resp', resp_tmpl_rand,  'iPRS', iPRS_tmpl_rand,  'respC', respC_tmpl_rand,  'cst', cst_tmpl_rand,  'RT', RT_tmpl_rand);
    d_train = struct('resp', resp_train_rand, 'iPRS', iPRS_train_rand, 'respC', respC_train_rand, 'cst', cst_train_rand, 'RT', RT_train_rand, 'iPair', iPair_train_rand);
    d_test  = struct('resp', resp_test_rand,  'iPRS', iPRS_test_rand,  'respC', respC_test_rand,  'cst', cst_test_rand,  'RT', RT_test_rand,  'iPair', iPair_test_rand);

    template_transformed = convertTemplateRawToTransformed(template_use, Tfix);
    DV_train_override = [];
    DV_test_override = [];

    if flag_regressType == 2
        basisOpts_A2 = struct();
        basisOpts_A2.nBasisORI = mode(nBasisORI_tmpl_allIter(isfinite(nBasisORI_tmpl_allIter)));
        basisOpts_A2.nBasisSF = mode(nBasisSF_tmpl_allIter(isfinite(nBasisSF_tmpl_allIter)));
        basisOpts_A2.basisFamilyORI = modeStringCell(basisFxnORI_tmpl_allIter);
        basisOpts_A2.basisFamilySF = modeStringCell(basisFxnSF_tmpl_allIter);
        basisOpts_A2.basisWidthScaleORI = mode(basisWidthScaleORI_tmpl_allIter(isfinite(basisWidthScaleORI_tmpl_allIter)));
        basisOpts_A2.basisWidthScaleSF = mode(basisWidthScaleSF_tmpl_allIter(isfinite(basisWidthScaleSF_tmpl_allIter)));
        basisOpts_A2.asymSF_rightLeftRatio = mode(asymSF_rightLeftRatio_tmpl_allIter(isfinite(asymSF_rightLeftRatio_tmpl_allIter)));
        basisOpts_A2.oriPeriod_deg = opts.oriPeriod_deg;

        Z_train = SX_RC_basisProject(e3D_train_forDV, axisORI, axisSF, basisOpts_A2);
        Z_test = SX_RC_basisProject(e3D_test_forDV, axisORI, axisSF, basisOpts_A2);
        beta_basis = SX_RC_basisProject(reshape(template_transformed, [1, nORI, nSF]), axisORI, axisSF, basisOpts_A2);
        DV_train_override = Z_train * beta_basis(:);
        DV_test_override = Z_test * beta_basis(:);
    end

    [data_train, data_test, metrics] = fxn_compDV_packData( ...
        d_full, d_tmpl, d_train, d_test, ...
        e3D_train_forDV, e3D_test_forDV, ...
        template_transformed, convolveType, DVType, ORI_bound, nBins, ...
        DV_train_override, DV_test_override);
    data_train_allIter{iIter} = data_train;
    data_test_allIter{iIter}  = data_test;
    data_metrics_allIter(iIter, :, :) = metrics;

    % Store templates
    % if iModelA_fit==2
    template_tmpl = template_transformed;
    template_full = template_transformed;
    % end
    template_tmpl_allIter(iIter, :, :) = template_tmpl;
    template_full_allIter(iIter, :, :) = template_full;

end % end parfor iIter

fprintf('\n\n%s: A2 All iterations done.\n\n', datetime('now'))

time_progress = ceil(minutes(datetime('now')-time_start)); % round up to minutes

%% [A2] SAVE
save(nameFile_compDV_A2, 'time_progress', 'template_ideal', '*_allIter', 'names*', 'flag*', 'ratio_split', 'ORI_bound', '*Type');

fprintf('\n\n%s: A2 Outputs saved.\n\n', datetime('now'))
end

%% HELPER
function [data_train, data_test, metrics] = fxn_compDV_packData( ...
    d_full, d_tmpl, d_train, d_test, ...
    e3D_train_base, e3D_test_base, ...
    template_notNormed, convolveType, DVType, ORI_bound, nBins, DV_train_override, DV_test_override)
% e3D_train_base and e3D_test_base must already be centered: E_raw - mu0

% 1. Compute behavioral metrics
metrics_full  = fxn_getMetrics(d_full.resp,  d_full.iPRS,  d_full.respC,  d_full.cst,  d_full.RT);
metrics_tmpl  = fxn_getMetrics(d_tmpl.resp,  d_tmpl.iPRS,  d_tmpl.respC,  d_tmpl.cst,  d_tmpl.RT);
metrics_train = fxn_getMetrics(d_train.resp, d_train.iPRS, d_train.respC, d_train.cst, d_train.RT);
metrics_test  = fxn_getMetrics(d_test.resp,  d_test.iPRS,  d_test.respC,  d_test.cst,  d_test.RT);
% [dprime, criterion, pC, pHit, pFA, nanmean(respC), pYES, mean(1./contrast), median(RT)]
metrics = [metrics_full; metrics_tmpl; metrics_train; metrics_test];
criterion_z_train = metrics_train(2);
criterion_z_test  = metrics_test(2);

% 2. Compute decision variable in transformed coordinate system.

flag_permT = 0;
if nargin >= 13 && ~isempty(DV_train_override) && ~isempty(DV_test_override)
    DV_train = DV_train_override(:);
    DV_test = DV_test_override(:);
else
    DV_train = fxn_getDV_v3(e3D_train_base, template_notNormed, convolveType, DVType, flag_permT, ORI_bound);
    DV_test  = fxn_getDV_v3(e3D_test_base,  template_notNormed, convolveType, DVType, flag_permT, ORI_bound);
end

nData_train = length(DV_train);
nData_test  = length(DV_test);

% 3. Bin DV values
[nTrials_PRS_allBins_train, ~, iTrial4Bin_PRS_train] = histcounts(DV_train(d_train.iPRS == 1), nBins);
[nTrials_ABS_allBins_train, ~, iTrial4Bin_ABS_train] = histcounts(DV_train(d_train.iPRS == 0), nBins);
[nTrials_PRS_allBins_test,  ~, iTrial4Bin_PRS_test]  = histcounts(DV_test(d_test.iPRS  == 1), nBins);
[nTrials_ABS_allBins_test,  ~, iTrial4Bin_ABS_test]  = histcounts(DV_test(d_test.iPRS  == 0), nBins);
[nTrials_allBins_train, ~, iTrial4Bin_train] = histcounts(DV_train, nBins);
[nTrials_allBins_test,  ~, iTrial4Bin_test]  = histcounts(DV_test,  nBins);

% 4 Compile data structs
data_train = struct();
data_train.criterion_z         = criterion_z_train;
data_train.ndata               = nData_train;
data_train.DV                  = DV_train;
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
data_test.DV                  = DV_test;
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
[nTrials, nORI, nSF] = size(e3D_in);
X = reshape(e3D_in, [nTrials, nORI * nSF]);
Xw = (X - mu_cov) * W_white;
e3D_white = reshape(Xw, [nTrials, nORI, nSF]);
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

function Tfix = buildFixedTransformFromABS(e3D_tmpl_rand, cst_tmpl_rand, iPRS_tmpl_rand, lambda_whiten, eps_whiten)
% Estimate transform parameters from ABS template trials only.
idxABS = (iPRS_tmpl_rand == 0);
assert(any(idxABS), 'No ABS trials available for fixed transform estimation.');

[e3D_abs_z, normStats_abs] = normEnergy(e3D_tmpl_rand(idxABS, :, :), cst_tmpl_rand(idxABS), zeros(sum(idxABS), 1));
[mu_cov_abs, ~, W_white_abs] = computeWhiteningParams(e3D_abs_z, lambda_whiten, eps_whiten);

Tfix = struct();
Tfix.normStats = normStats_abs;
Tfix.mu_cov_abs = mu_cov_abs;
Tfix.W_white = W_white_abs;
Tfix.useWhiten = ~isempty(W_white_abs);
end

function e3D_out = applyFixedTransform(e3D_in, cst_in, Tfix)
% Apply contrast-specific z-scoring (ABS-derived) and optional whitening.
[nTrials, nORI, nSF] = size(e3D_in);
e3D_z = nan(size(e3D_in));
eps_sigma = 1e-8;

cst_unique = Tfix.normStats.cst_unique(:)';
for iT = 1:nTrials
    [~, idxC] = min(abs(cst_unique - cst_in(iT)));
    mu_i = squeeze(Tfix.normStats.mu_byGroup(idxC, 1, :, :));
    sd_i = squeeze(Tfix.normStats.sigma_byGroup(idxC, 1, :, :));
    if any(~isfinite(mu_i(:))) || any(~isfinite(sd_i(:)))
        mu_i = squeeze(Tfix.normStats.mu_global_3D);
        sd_i = Tfix.normStats.sigma_global_2D;
    end
    sd_i(~isfinite(sd_i) | sd_i < eps_sigma) = eps_sigma;
    e3D_z(iT, :, :) = (squeeze(e3D_in(iT, :, :)) - mu_i) ./ sd_i;
end

if Tfix.useWhiten
    Xz = reshape(e3D_z, [nTrials, nORI * nSF]);
    Xw = bsxfun(@minus, Xz, Tfix.mu_cov_abs) * Tfix.W_white;
    e3D_out = reshape(Xw, [nTrials, nORI, nSF]);
else
    e3D_out = e3D_z;
end
end

function template_transformed = convertTemplateRawToTransformed(template_raw, Tfix)
% Map raw-centered template into transformed coordinate system.
sigma_vec = Tfix.normStats.sigma_global_2D(:);
sigma_vec(~isfinite(sigma_vec) | sigma_vec < 1e-8) = 1e-8;
t_z = template_raw(:) .* sigma_vec;
if Tfix.useWhiten
    t_trans = Tfix.W_white \ t_z;
else
    t_trans = t_z;
end
template_transformed = reshape(t_trans, size(template_raw));
end

function out = modeStringCell(C)
% Return modal non-empty string from a cell array.
S = string(C(:));
S = S(strlength(S) > 0 & S ~= "<missing>");
if isempty(S)
    out = '';
    return;
end
u = unique(S);
cnt = zeros(size(u));
for i = 1:numel(u)
    cnt(i) = sum(S == u(i));
end
[~, iMax] = max(cnt);
out = char(u(iMax));
end

function setupParforEnvironment(codeRoot)
% Use absolute paths and synchronize worker path/source visibility.

addpath(codeRoot);
addpath(genpath(fullfile(codeRoot, 'fxn_exp')));
addpath(genpath(fullfile(codeRoot, 'fxn_NOM')));
addpath(genpath(fullfile(codeRoot, 'fxn_RCplot')));
addpath(genpath(fullfile(codeRoot, 'fxn_analysis_RC_v2')));
addpath(genpath(fullfile(codeRoot, 'SX_toolbox', 'bads-master')));

pool = gcp('nocreate');
if isempty(pool)
    return;
end

cmds = {
    sprintf('addpath(''%s'');', codeRoot)
    sprintf('addpath(genpath(''%s''));', fullfile(codeRoot, 'fxn_exp'))
    sprintf('addpath(genpath(''%s''));', fullfile(codeRoot, 'fxn_NOM'))
    sprintf('addpath(genpath(''%s''));', fullfile(codeRoot, 'fxn_RCplot'))
    sprintf('addpath(genpath(''%s''));', fullfile(codeRoot, 'fxn_analysis_RC_v2'))
    sprintf('addpath(genpath(''%s''));', fullfile(codeRoot, 'SX_toolbox', 'bads-master'))
    };

for iCmd = 1:numel(cmds)
    pctRunOnAll(cmds{iCmd});
end

attachFiles = {
    fullfile(codeRoot, 'OOD_NOM_Trialwise_compDV_A12.m')
    fullfile(codeRoot, 'fxn_analysis_RC_v2', 'SX_RC_basisProject.m')
    fullfile(codeRoot, 'fxn_analysis_RC_v2', 'normEnergy.m')
    };
attachFiles = attachFiles(cellfun(@(f) exist(f, 'file') == 2, attachFiles));
if ~isempty(attachFiles)
    addAttachedFiles(pool, attachFiles);
    updateAttachedFiles(pool);
end
end
%
% function kernel_corr = invcov_kernel(kernel_raw, W_inv, nOri, nSf)
% k_raw = kernel_raw(:);
% k_corr = W_inv * k_raw;
% kernel_corr = reshape(k_corr, [nOri, nSf]);
% end

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

function bestX = runMultistartFmincon(objFun, x0, lb, ub, nRep, options_fmin)
% Parfor-safe replacement for MultiStart.run.
nRep = max(1, round(nRep));
bestX = x0;
bestF = inf;

for iRep = 1:nRep
    if iRep == 1
        xStart = x0;
    else
        xStart = lb + rand(size(x0)) .* (ub - lb);
    end

    [xCand, fCand] = fmincon(objFun, xStart, [], [], [], [], lb, ub, [], options_fmin);
    if isfinite(fCand) && fCand < bestF
        bestF = fCand;
        bestX = xCand;
    end
end
end

function template_z = whiteTemplate_to_zTemplate(template_white, W_white, nOri, nSf)
% Convert a full template estimated in whitened-energy space back to the
% centered standardized-energy space.

beta_z = W_white * template_white(:);
template_z = reshape(beta_z, [nOri, nSf]);
end
