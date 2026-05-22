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

ratio_split = [.8, .15, .05]; % proportion of trials in template set (for RC), training set (for estimating parameters) and test set (for metric predictions)
assert(abs(sum(ratio_split)-1)<1e-10)
ORI_bound = [5, 14]; % orientation window, passed to fxn_getDV_v3

% whiten features
eps_whiten = 1e-3;  % floor for eigenvalues

% Setting for multivariate regression with smoothing
% 1=Univariate; 2=MultiSmooth and Univariate
if flag_regressType == 2
    basisCfg = struct();
    basisCfg.basisFamilyORI = 'vonmises';
    basisCfg.basisFamilySF = 'asymGaussianLog2';
    basisCfg.asymSF_rightLeftRatio = 1.2;
    basisCfg.oriPeriod_deg = 180;
    basisCfg.nBasisORI = 3; go back to the values reported in the paper
    basisCfg.nBasisSF = 7;
    basisCfg.basisWidthORI = .1;
    basisCfg.basisWidthSF = .1;
    basisCfg.Ridge = 100;
    basisCfg.nFolds = 5;
    basisCfg.link = 'probit';
    basisCfg.sigmaORI_deg = [];
    basisCfg.sigmaSF_log2 = [];
    basisCfg.zscorePredictor = true;
    basisCfg.maxIter = 100;
    basisCfg.tol = 1e-6;
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
% flag_standEnergy = 1; % 1=z-score energy before RC
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
    '\n - Energy source: %d (1=TARGET, 2=NOISE)', ...
    '\n - Convolve type: %s', ...
    '\n - DV type: %s\n\n'], ...
    subjName, ...
    iLocComb, namesLocComb{iLocComb}, ...
    nIter, ...
    nblocks, nORI, nSF, ...
    namesType{itype_template}, flag_PatchMode, ...
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

%% Create Gabor filters used for energy computation
[filter_sin, filter_cos] = SX_sim02_setFilters(stim, noise.filtersSF_all, filtersOri_all, fxn_getSigma_SPdomain, 0);
fprintf('%s: Filter banks (nORI=%d, nSF=%d) created and saved.\n\n', datetime('now'), length(filtersOri_all), length(noise.filtersSF_all))

%% Define the IDEAL template (ORI x SF)
% The ideal template is the energy of the gabor at the averaged CST across signal-prs trials
gaborCST = mean(cst_nonrand(cst_nonrand>0));
templateType_true = 1; % doesn't matter much for the IO Gabor template
stim.gaborCST = gaborCST;

template_gabor_ideal = exp_CreateGabor(stim, stim.gaborCST);
%--------------------------------------------%
template_ideal = SX_RC4_Energy_parfor(stim.mask, {template_gabor_ideal}, filter_sin, filter_cos);
%--------------------------------------------%
template_ideal = squeeze(template_ideal); % remove singleton dim
%--------------------------------------------%
template_ideal = fxn_getTemplate(template_ideal, templateType_true, 0);
%--------------------------------------------%

% Do L2 normalization
% MAKE SURE THIS IS the same as OOD_xx_compDV_A12 when creating the ideal template
template_ideal = template_ideal / norm(template_ideal(:));

fprintf('%s: Ideal template created.\n\n', datetime('now'))

%% Stage overview
% A1: estimate trial-wise templates from data and compute DVs.
% A2: keep the same trial resampling and downstream DV/metric computation,
% but skip RC estimation entirely and use the fixed ideal/true template in
% raw channel-energy space.
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
basisWidthORI_tmpl_allIter = nBasisSF_tmpl_allIter;
basisWidthSF_tmpl_allIter  = nBasisSF_tmpl_allIter;
asymSF_rightLeftRatio_tmpl_allIter = nBasisSF_tmpl_allIter;
% nLL_tmpl_allIter         = nan(nIter, 1); % training NLL of selected model (tmpl set)
nLL_tmpl_allIter = nan(nIter, 1); % mean held-out NLL from CV model selection
% nLL_cv_se_tmpl_allIter   = nan(nIter, 1); % SE of held-out NLL across CV folds
pseudoR2_Tjur_tmpl_allIter = nan(nIter, 1); % Tjur R² on training set

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

for iIter = 1:nIter
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
    resp_tmpl_rand_sel = resp_tmpl_rand(useIdx_tmpl);
    iPRS_tmpl_rand_sel = iPRS_tmpl_rand(useIdx_tmpl);

    e3D_full_rand_sel = e3D_full_rand(useIdx_full, :, :);
    resp_full_rand_sel = resp_full_rand(useIdx_full);
    iPRS_full_rand_sel = iPRS_full_rand(useIdx_full);

    %% [A1] 2b. Build fixed transform from ABS template trials only
    Tfix = fxn_buildFixedTransformFromABS(e3D_tmpl_rand, iPRS_tmpl_rand, lambda_whiten, eps_whiten);

    % Apply the transform ONLY to inputs of RC template estimation.
    e3D_tmpl_forRC = applyFixedTransform(e3D_tmpl_rand_sel, Tfix);
    e3D_full_forRC = applyFixedTransform(e3D_full_rand_sel, Tfix);

    %% [A1] 3. Estimate template from whitened energy
    % RC regression operates in the transformed (whitened) space.
    % The resulting template_*_white lives in that same space.

    % Safe defaults for basis function outputs (case 2 only)
    nBasisORI_tmpl = nan;
    nBasisSF_tmpl = nan;
    basisFxnORI_tmpl = '';
    basisFxnSF_tmpl = '';
    basisWidthORI_tmpl = nan;
    basisWidthSF_tmpl = nan;
    asymSF_rightLeftRatio_tmpl = nan;
    ridge_tmpl = nan;
    nLL_tmpl = nan;
    nLL_cv_mean_tmpl = nan;
    nLL_cv_se_tmpl = nan;
    pseudoR2_Tjur_tmpl = nan;

    % Estimate template in whitened space via regression
    switch flag_regressType
        case 1  % Univariate
            template_tmpl_white = SX_sim07_RC(e3D_tmpl_forRC, resp_tmpl_rand_sel);
            template_full_white = SX_sim07_RC(e3D_full_forRC, resp_full_rand_sel);

        case 2  % Multivariate + smoothing
            basisCfg_thisIter = basisCfg;
            basisCfg_thisIter.template_ideal = template_ideal;

            % For the "template set"
            out = SX_RC_selectBasis_cv(e3D_tmpl_forRC, resp_tmpl_rand_sel, axisORI, axisSF, basisCfg_thisIter);

            % Store the outputs of the best fit
            template_tmpl_white = out.template2D;
            nLL_cv_mean_tmpl = out.bestRes.nLL_mean; % This is the mean held-out nLL across CV folds for the selected model, so primary for model quality
            pseudoR2_Tjur_tmpl = out.pseudoR2_Tjur;
            % nLL_tmpl = out.nLL; % This is training nLL from the best fit, so only secondary

            % Store the selected hyperparameters for reference
            nBasisORI_tmpl = out.nBasisORI;
            nBasisSF_tmpl = out.nBasisSF;
            basisFxnORI_tmpl = out.basisFamilyORI;
            basisFxnSF_tmpl = out.basisFamilySF;
            basisWidthORI_tmpl = out.basisWidthORI;
            basisWidthSF_tmpl = out.basisWidthSF;
            asymSF_rightLeftRatio_tmpl = out.asymSF_rightLeftRatio;
            ridge_tmpl = out.Ridge;

            % Same thing for the full set, but just the best template
            out = SX_RC_selectBasis_cv(e3D_full_forRC, resp_full_rand_sel, axisORI, axisSF, basisCfg_thisIter);
            template_full_white = out.template2D;

    end % switch

    %% [A1] 3b. Convert estimated template back to RAW channel-energy space
    % The DV is computed with raw energy times a raw-space template, so
    % undo the (z-score + whitening) mapping that was applied to RC inputs.
    template_tmpl_raw = fxn_convertTemp2Raw(template_tmpl_white, Tfix);
    template_full_raw = fxn_convertTemp2Raw(template_full_white, Tfix);

    template_tmpl_raw = template_tmpl_white;
    template_full_raw = template_full_white;

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
            % DV is ALWAYS computed from RAW channel energy times the
            % raw-space template. No basis-projection override.
            [data_train, data_test, metrics] = fxn_compDV_packData( ...
                d_full, d_tmpl, d_train, d_test, ...
                e3D_train_rand, e3D_test_rand, ...
                template_use_raw, convolveType, DVType, ORI_bound, nBins, ...
                [], []);
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
    basisWidthORI_tmpl_allIter(iIter) = basisWidthORI_tmpl;
    basisWidthSF_tmpl_allIter(iIter) = basisWidthSF_tmpl;
    asymSF_rightLeftRatio_tmpl_allIter(iIter) = asymSF_rightLeftRatio_tmpl;
    ridge_tmpl_allIter(iIter) = ridge_tmpl;

    % GoF
    nLL_tmpl_allIter(iIter) = nLL_cv_mean_tmpl;
    % nLL_cv_se_tmpl_allIter(iIter) = nLL_cv_se_tmpl;
    pseudoR2_Tjur_tmpl_allIter(iIter) = pseudoR2_Tjur_tmpl;

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
% A2 uses the ideal/true template DIRECTLY in raw channel-energy space.
% No estimation, no peak-rescaling, no whitening, no basis projection.
if ~isnumeric(isubj) % for IO, load 'template_true' from truth.mat
    load(sprintf('%s/truth.mat', nameFolder_OOD_load), 'template_true');
    template_use = template_true;
else
    template_use = template_ideal;
end

% template_true is NOT euqal to template_ideal !!!

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

    %% [A2] 2. (No transform) — A2 uses raw channel energy for DV
    % To keep parity with A1's `selectTemplateIndices` plumbing (some
    % downstream code references *_sel variables), compute them but do NOT
    % build/apply any whitening transform on the DV inputs.
    [useIdx_full, useIdx_tmpl] = selectTemplateIndices(itype_template, iPRS_full_rand, iPRS_tmpl_rand);

    e3D_tmpl_rand_sel = e3D_tmpl_rand(useIdx_tmpl, :, :);
    cst_tmpl_rand_sel = cst_tmpl_rand(useIdx_tmpl);
    iPRS_tmpl_rand_sel = iPRS_tmpl_rand(useIdx_tmpl);

    e3D_full_rand_sel = e3D_full_rand(useIdx_full, :, :);
    cst_full_rand_sel = cst_full_rand(useIdx_full);
    iPRS_full_rand_sel = iPRS_full_rand(useIdx_full);

    %% [A2] 5-8. Metrics, DV, binning (raw energy, raw template)
    d_full  = struct('resp', resp_full_rand,  'iPRS', iPRS_full_rand,  'respC', respC_full_rand,  'cst', cst_full_rand,  'RT', RT_full_rand);
    d_tmpl  = struct('resp', resp_tmpl_rand,  'iPRS', iPRS_tmpl_rand,  'respC', respC_tmpl_rand,  'cst', cst_tmpl_rand,  'RT', RT_tmpl_rand);
    d_train = struct('resp', resp_train_rand, 'iPRS', iPRS_train_rand, 'respC', respC_train_rand, 'cst', cst_train_rand, 'RT', RT_train_rand, 'iPair', iPair_train_rand);
    d_test  = struct('resp', resp_test_rand,  'iPRS', iPRS_test_rand,  'respC', respC_test_rand,  'cst', cst_test_rand,  'RT', RT_test_rand,  'iPair', iPair_test_rand);

    [data_train, data_test, metrics] = fxn_compDV_packData( ...
        d_full, d_tmpl, d_train, d_test, ...
        e3D_train_rand, e3D_test_rand, ...
        template_use, convolveType, DVType, ORI_bound, nBins, ...
        [], []);
    data_train_allIter{iIter} = data_train;
    data_test_allIter{iIter}  = data_test;
    data_metrics_allIter(iIter, :, :) = metrics;

    % Store templates (raw space)
    template_tmpl = template_use;
    template_full = template_use;
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

%%
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

%%
function e3D_out = applyFixedTransform(e3D_in, Tfix)
% Apply per-channel z-scoring and (optionally) the whitening matrix Q.
% Used ONLY as a preconditioner for RC template estimation.
%   x_z     = (x_raw - mu_abs) ./ sigma_abs
%   x_white = (x_z   - mu_cov_abs) * Q             (only if useWhiten)
[nTrials, nORI, nSF] = size(e3D_in);
X = reshape(e3D_in, [nTrials, nORI * nSF]);
X_z = (X - Tfix.mu_abs) ./ Tfix.sigma_abs;

if Tfix.useWhiten
    Xw = bsxfun(@minus, X_z, Tfix.mu_cov_abs) * Tfix.Q;
    e3D_out = reshape(Xw, [nTrials, nORI, nSF]);
else
    e3D_out = reshape(X_z, [nTrials, nORI, nSF]);
end
end

%%
function template_converted = fxn_convertTemp2Raw(template_white, Tfix)
% Map an RC-estimated template from the whitened/z-scored space

% Inverse mapping implemented below.
sigma_vec = Tfix.sigma_abs(:);
sigma_vec(~isfinite(sigma_vec) | sigma_vec < 1e-8) = 1e-8;

temp_white_vec = template_white(:);
if Tfix.useWhiten
    % Q is symmetric SPD-like; Q * t_white inverts the (Q \ t_z) forward step.
    t_z = Tfix.Q * temp_white_vec;
else
    t_z = temp_white_vec;
end
t_raw = t_z ./ sigma_vec;
template_converted = reshape(t_raw, size(template_white));
end

%%
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
