function OOD_NOM_Trialwise_compIV(isubj, iLocComb, iModelA, nIter, iJob, nJob)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OOD_NOM_Trialwise_compIV.m
%
% Trial-wise noisy observer model – PRE-ESTIMATION STAGE
%
% This function:
% 1. Splits trials into template, training and test sets (via fxn_resampleTrials).
% 2. Estimates a template from the template set (using stim energy + resp).
% 3. Computes behavioral metrics in the test set (these are the metrics the model will later try to match).
% 4. Computes, normalizes, and bins empirical internal variables (IVs) for the test set.

% The output is a .mat file containing, for each  iteration:
% - data_allB{ii}: struct with IVs, binning info, responses, etc.
% - data_metrics_allB(ii,:): behavioral metrics per  iteration
% - kernel2D_allB(ii,:,:): 2D templates (ORI x SF) per  iteration
%
% INPUTS
% isubj : index of subject or IO info
% - numeric: index into subjList (human observers)
% - cell: {nameIO, ...} for ideal observer
% iLocComb : location combination index
% 1=fovea; 8=periF (6 deg ecc);
% 6=HM; 7=VM; 5=LVM; 3=UVM
% iModelA : template model
% 1 = RC-derived template
% 2 = IO template
% 3 = permuted template
% IVType : how IV is computed from template × energy
% 1 = sum of dot product / convolution
% 2 = max
% 3 = normalized, etc. (see fxn_getIV_v3)
% templateType : how to extract 1D template from 2D kernel
% 1 = raw
% 2 = reconstructed kernel
% 3 = mirrored template
% flag_PatchMode : 1 = use target patch energy, 2 = use noise patch energy
% itype_template : which trials to use for estimating the template
% 1 = PRS (signal-present) trials only
% 2 = ABS (signal-absent) trials only
% 3 = both PRS and ABS (not implemented here, falls through)
% nIter : number of resampling iterations
%
% OUTPUTS (saved to disk)
% nameFile_compIV = '.../n%d_A%d_compIV.mat', containing:
% - c_zscore : SDT criterion (z units) per  iteration
% - data_allIter : cell array of "data" structs per  iteration
% - data_metrics_allIter : behavioral metrics per  iteration
% - Template_train_allIter : 2D kernels per  iteration (using trials in the training set)
% - Template_full_allIter : 2D kernels per  iteration (using all trials)
% - names*, flag*, ratio_train, ORI_bound, *Type, etc.
%
% Notes:
% - This is the "before estimation" stage; model fits (iModelB) happen in
% OOD_NOM_Trialwise_Est.m
% - IVs are computed in fxn_getIV_v3, typically as template ⊙ energy.
%
% Created by Shutian Xue on August 27, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all;
warning off; % (You may want to remove this once things are stable.)
format compact;
time_start = datetime('now')

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

%% -------------------- Deterministic RNG (grand seed + per-iteration substreams) -------------------- %%
S_seed = GetGrandSeed(nIter, iJob, nJob, nameFolder_Data);

% One RNG stream for the whole job; each iteration uses its own Substream
stream = RandStream('Threefry', 'Seed', S_seed.grandSeed);
RandStream.setGlobalStream(stream);

%% -------------------- General parameters -------------------- %%

IVType = 1; % 1=sum of the dot product/convolution; 2=max; 3=normalized
templateType = 3; % (1) raw (2) reconstructed kernel (3) mirrored template
itype_template = 2; % 1=estimate template from PRS trials, ABS trials, or BOTH trials
flag_PatchMode = 1; % if flag_PatchMode == 1, patchMode = 'T'; else, patchMode = 'N'; end
flag_plot_template = 0;

% Settings for fitting tuning functions
nRep = 20;
iFamily_ORI = 1; % 1=scaled gaussian, 8=DoG
iFamily_SF = 2; % 2=log parabola
problem_setting = MultiStart('StartPointsToRun', 'bounds','UseParallel', 1, 'Display', 'off');
flag_plot_tuning = 0;

iSess_start = 1; % first session included
convolveType = 1; % IV from 1=cross-correlation; 2=convolution (fxn_getIV_v3)
flag_standEnergy = 1; % 1=z-score energy before RC
flag_plot_compIV = 0; % plot IV distributions and kernels at the end
% ratio_train = 3/4; % proportion of trials in training set (for RC)
ratio_split = [.6, .3, .1]; % proportion of trials in template set (for RC), training set (for estimating parameters) and test set (for metric predictions)
ORI_bound = [10, 20]; % orientation window, passed to fxn_getIV_v3

if flag_PatchMode == 1
    namePatchMode = 'T'; % target patch
else
    namePatchMode = 'N'; % noise patch
end
if ~strcmp('HPC', str_envir), flag_plot_compIV = 1; end % don't plot when running on HPC

%% -------------------- Subject / IO and paths -------------------- %%
if isnumeric(isubj) % Human subjects
    subjList = {'YK','SP','SX','LS','RE','MD','AS','HL','FH','HA','CS','DT','DU','RC','SR'};
    nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195, 0];

    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);

    % Folder to load behav + energy
    nameFolder_OOD_load = sprintf('%s/%s%d', nameFolder_Data_OOD, subjName, nblocks);

    % Folder to save NOM trial-wise results
    nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocComb);

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

    % Behavioral measurements
    load(sprintf('%s/behavMeas.mat', nameFolder_OOD_load), 'dataMatrix');

    % Energy
    load(sprintf('%s/energy_%s_%d_%d.mat', nameFolder_OOD_load, namePatchMode, nORI, nSF), ...
        'e3D_target_allT', 'e3D_noise_allT', 'noise', 'filtersOri_all', 'stim', 'nBins');
end

if isempty(dir(nameFolder_NOM_save))
    mkdir(nameFolder_NOM_save);
end

% fprintf('\nBehav and energy LOADED\n\n');

% Output file name (before estimation)
nameFile_compIV = sprintf('%s/n%d_J%d_A%d_compIV', nameFolder_NOM_save, nIter, iJob, iModelA);

%% Print run info
fprintf('\n==================================\nStep 1: Compute IVs and derive templates\n================================== \n\n')
fprintf(['\nSubject/IO name: %s ' ...
    '\n - L%d [%s]', ...
    '\n - A%d [%s]', ...
    '\n - Number of iterations = %d', ...
    '\n - Number of blocks = %d', ...
    '\n - nORI=%d | nSF=%d', ...
    '\n - Template derived from %s trials', ...
    '\n - Energy source: %d (1=TARGET, 2=NOISE)', ...
    '\n - Convolve type: %s', ...
    '\n - IV type: %s\n\n'], ...
    subjName, ...
    iLocComb, namesLocComb{iLocComb}, ...
    iModelA, namesModelA{iModelA}, ...
    nIter, ...
    nblocks, nORI, nSF, ...
    namesType{itype_template}, flag_PatchMode, ...
    namesConvolveType{convolveType}, namesIVType{IVType});

%% -------------------- Choose energy source for NOM -------------------- %%
switch flag_PatchMode
    case 1 % target patch energy
        % e3D_allT_train = e3D_target_allT;
        % e3D_allT_test = e3D_target_allT;
        e3D_allT = e3D_target_allT;
    case 2 % noise patch energy
        % e3D_allT_train = e3D_noise_allT;
        % e3D_allT_test = e3D_noise_allT;
        e3D_allT = e3D_noise_allT;
    otherwise
        error('flag_PatchMode must be 1 (target) or 2 (noise).');
end

%% -------------------- Sanity check: iPRS consistent within pairs -------------------- %%
% Col#1 = itrial_allT; Col#6 = iPRS; Col#8 = iPair;
for ipair = 1:max(dataMatrix(:, 8))
    iIter = find(dataMatrix(:, 8) == ipair);
    iPRS_A = dataMatrix(dataMatrix(:, 1) == iIter(1), 6);
    iPRS_B = dataMatrix(dataMatrix(:, 1) == iIter(2), 6);
    if iPRS_A ~= iPRS_B
        fprintf('WARNING: when ipair=%d, iPRS_A=%d, iPRS_B=%d\n', ipair, iPRS_A, iPRS_B);
    end
end

%% -------------------- Session and location selection -------------------- %%
nSess = length(unique(dataMatrix(:, 2)));

% Assume equal #trials per location; use location==1 as representative.
nTrials_perSingleLoc = sum(dataMatrix(:, 5) == 1);

% Sessions included (can be subsetted here)
iSess_select = iSess_start:nSess;

% Special case: subject 'RE' missing some sessions
if strcmp(subjName, 'RE')
    iSess_select = [1:15, 17:23, 25:nSess];
    nTrials_perSingleLoc = length(iSess_select) * 100;
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

% %% -------------------- Templare/Train/Test split per location -------------------- %%
% nTemp_perSingle = round(nTrials_perSingleLoc * ratio_split(1));
% nTrain_perSingle = round(nTrials_perSingleLoc * ratio_split(2));
% nTest_perSingle = round(nTrials_perSingleLoc * ratio_split(3));
% 
% % Ensure even #trials in TRAIN (so we keep full pairs)
% assert(rem(nTemp_perSingle, 2)==0, 'ALERT: Number of trials is not an even number!!')
% % nTemp_perSingle = nTemp_perSingle + 1;
% % nTest_perSingle = nTest_perSingle - 1;
% % end


%% -------------------- Ideal template -------------------- %%
% For human observers, we still use this "ideal" Gabor-energy template when
% iModelA == 3 (IO template). Otherwise, templates are derived from RC.

% Precompute non-random pools 
e3D_nonrand = [];
resp_nonrand = [];
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

% Define the TRUE template (ORI x SF)
template_gabor_true = exp_CreateGabor(stim, stim.gaborCST);
template_true = SX_RC4_Energy_parfor(stim.mask, {template_gabor_true}, filter_sin, filter_cos);
template_true = squeeze(template_true); % remove singleton dim
template_true = fxn_getTemplate(template_true, templateType_true, 0);

%% -------------------- MAIN LOOP over  iterations -------------------- %%
data_metrics_allIter = nan(nIter, nDatasets_full, nMetrics); % see fxn_getMetrics for all 11 metrics
data_train_allIter = cell(nIter, 1);
data_test_allIter = data_train_allIter;
template_tmpl_allIter = nan(nIter, nORI, nSF);
template_full_allIter = template_tmpl_allIter;

sep_allIter = nan(nIter, 2); % 1=template set, 2=full set
margORI_allIter = nan(nIter, 2, nORI);
margPred_ORI_allIter = nan(nIter, 2, nORI);
margParams_ORI_allIter = nan(nIter, 2, length(namesParams_all{iFamily_ORI}));
margR2_ORI_allIter = nan(nIter, 2);
margSF_allIter = nan(nIter, 2, nSF);
margPred_SF_allIter = nan(nIter, 2, nSF);
margParams_SF_allIter = nan(nIter, 2, length(namesParams_all{iFamily_SF}));
margR2_SF_allIter = nan(nIter, 2);

fprintf('\n\nRunning nIter = %d: ', nIter);

for iIter = 1:nIter
    fprintf('%d... ', iIter);

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
            useIdx_train = (iPRS_train_rand == 1);
            useIdx_test = (iPRS_test_rand == 1);
        case 2 % ABS trials only
            useIdx_full = (iPRS_full_rand == 0);
            useIdx_tmpl = (iPRS_tmpl_rand == 0);
            useIdx_train = (iPRS_train_rand == 0);
            useIdx_test = (iPRS_test_rand == 0);
        otherwise % 3 = both PRS and ABS (here we just keep everything)
            useIdx_full = true(size(iPRS_full_rand));
            useIdx_tmpl = true(size(iPRS_tmpl_rand));
            useIdx_train = true(size(iPRS_train_rand));
            useIdx_test = true(size(iPRS_test_rand));
    end

    % sel for "selected"
    e3D_tmpl_rand_sel = e3D_tmpl_rand(useIdx_tmpl, :, :);
    cst_tmpl_rand_sel = cst_tmpl_rand(useIdx_tmpl);
    resp_tmpl_rand_sel = resp_tmpl_rand(useIdx_tmpl);
    iPRS_tmpl_rand_sel = iPRS_tmpl_rand(useIdx_tmpl);
    RT_tmpl_rand_sel = RT_tmpl_rand(useIdx_tmpl);

    e3D_full_rand_sel = e3D_full_rand(useIdx_full, :, :);
    cst_full_rand_sel = cst_full_rand(useIdx_full);
    resp_full_rand_sel = resp_full_rand(useIdx_full);
    iPRS_full_rand_sel = iPRS_full_rand(useIdx_full);
    RT_full_rand_sel = RT_full_rand(useIdx_full);

    % Standardize energy if requested
    if flag_standEnergy
        e3D_tmpl_rand_norm = normEnergy(e3D_tmpl_rand_sel, cst_tmpl_rand_sel, iPRS_tmpl_rand_sel);
        e3D_full_rand_norm = normEnergy(e3D_full_rand_sel, cst_full_rand_sel, iPRS_full_rand_sel);
    else
        e3D_tmpl_rand_norm = e3D_tmpl_rand_sel;
        e3D_full_rand_norm = e3D_full_rand_sel;
    end

    %% 3. Estimate the template
    % from TEMPLATE set
    template_tmpl = SX_sim07_RC(filtersSF_all, filtersOri_all, e3D_tmpl_rand_norm, resp_tmpl_rand_sel);
    template_tmpl = fxn_getTemplate(template_tmpl, templateType, flag_plot_template);
    
    % from the FULL set
    template_full = SX_sim07_RC(filtersSF_all, filtersOri_all, e3D_full_rand_norm, resp_full_rand_sel);
    template_full = fxn_getTemplate(template_full, templateType, flag_plot_template);
    
    %% 4. Marginalization and fit tuning functions
    for iDataset = 1:2
        switch iDataset
            case 1
                template_InUse = template_tmpl;
            case 2
                template_InUse = template_full;
        end

        % Calculate separability
        templateType_recon = 2; % 2=reconstructed template
        template_recon = fxn_getTemplate(template_InUse, templateType_recon, flag_plot_template);
        sep = corr(template_InUse(:), template_recon(:));

        % Marginalize the template into 1D ORI and SF profiles
        margORI = mean(template_InUse, 2)';
        margSF = mean(template_InUse, 1);

        % Fitting
        xORI = axis_tuning{1};
        fxn_tuningLoss_ORI = @(param_est) sum((margORI - predSFkernel(xORI, iFamily_ORI, param_est, flag_plot_tuning)).^2);
        problem_ORI = createOptimProblem('fmincon','objective', fxn_tuningLoss_ORI,'x0', (ub_full_all{iFamily_ORI}+lb_full_all{iFamily_ORI})/2,'lb',lb_full_all{iFamily_ORI},'ub',ub_full_all{iFamily_ORI},'options',options_fmin);
        margParams_ORI = run(problem_setting, problem_ORI, nRep);
        margPred_ORI = predSFkernel(xORI, iFamily_ORI, margParams_ORI, flag_plot_tuning);
        margR2_ORI = 1-sumsqr(margORI-margPred_ORI)/sumsqr(margORI-mean(margORI));

        xSF = axis_tuning{2};
        xSF_ln = 2.^xSF; % fit in linear SF space
        fxn_tuningLoss_SF = @(param_est) sum((margSF - predSFkernel(xSF_ln, iFamily_SF, param_est, flag_plot_tuning)).^2);
        problem_SF = createOptimProblem('fmincon','objective', fxn_tuningLoss_SF,'x0', (ub_full_all{iFamily_SF}+lb_full_all{iFamily_SF})/2,'lb',lb_full_all{iFamily_SF},'ub',ub_full_all{iFamily_SF},'options',options_fmin);
        margParams_SF = run(problem_setting, problem_SF, nRep);
        margPred_SF = predSFkernel(xSF_ln, iFamily_SF, margParams_SF, flag_plot_tuning);
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
    metrics_full = fxn_getMetrics(resp_full_rand, iPRS_full_rand, respC_full_rand, cst_full_rand, RT_full_rand);
    metrics_tmpl = fxn_getMetrics(resp_tmpl_rand, iPRS_tmpl_rand, respC_tmpl_rand, cst_tmpl_rand, RT_tmpl_rand);
    metrics_train = fxn_getMetrics(resp_train_rand, iPRS_train_rand, respC_train_rand, cst_train_rand, RT_train_rand);
    metrics_test = fxn_getMetrics(resp_test_rand, iPRS_test_rand, respC_test_rand, cst_test_rand, RT_test_rand);

    % metrics: [dprime, criterion, [pC, pHit, pFA], nanmean(respC), pYES, mean(1./contrast), median(RT)];
    metrics = [metrics_full; metrics_tmpl; metrics_train; metrics_test];

    % Save criterion in z units; later used to convert to criterion in IV
    c_zscore_train = metrics_train(2);
    c_zscore_test = metrics_test(2);
    % CONFIRM if c_zscore should be from TRAIN set or TEST set!!

    clear data; % size of IV differs across subjects and  iterations

    %% 6. Compute internal variable (IV) from template and energy (TEST set)
    % 5.1 Derive template
    if iModelA == 1 % Use RC-derived template
    else % Use the ideal template (Gabor energy profile)
        if ~isnumeric(isubj) % for IO, you may reload 'template_true' from disk
            load(sprintf('%s/truth.mat', nameFolder_OOD_load), 'template_true');
        end
        template_tmpl = template_true;
    end

    % 5.2 Compute IV for training and test trials
    %---------------%
    IV_train = fxn_getIV_v3(iModelA, e3D_train_rand, convolveType, IVType, template_tmpl, ORI_bound);
    IV_test = fxn_getIV_v3(iModelA, e3D_test_rand, convolveType, IVType, template_tmpl, ORI_bound);

    %---------------%
    nData_train = length(IV_train);
    nData_test = length(IV_test);

    %% 7. Bin IV values
    % Separately for PRS and ABS trials
    [nTrials_PRS_allBins_train, ~, iTrial4Bin_PRS_train] = histcounts(IV_train(iPRS_train_rand == 1), nBins);
    [nTrials_ABS_allBins_train, ~, iTrial4Bin_ABS_train] = histcounts(IV_train(iPRS_train_rand == 0), nBins);
    [nTrials_PRS_allBins_test, ~, iTrial4Bin_PRS_test] = histcounts(IV_test(iPRS_test_rand == 1), nBins);
    [nTrials_ABS_allBins_test, ~, iTrial4Bin_ABS_test] = histcounts(IV_test(iPRS_test_rand == 0), nBins);

    % For all trials combined
    [nTrials_allBins_train, ~, iTrial4Bin_train] = histcounts(IV_train, nBins);
    [nTrials_allBins_test, ~, iTrial4Bin_test] = histcounts(IV_test, nBins);

    %% 8. Compile data struct for this  iteration
    % Store data_train
    data_train = struct();
    data_train.c_zscore = c_zscore_train;
    data_train.ndata = nData_train;
    data_train.IV = IV_train;
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
    data_test.c_zscore = c_zscore_test;
    data_test.ndata = nData_test;
    data_test.IV = IV_test;
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
    template_tmpl_allIter(iIter, :, :) = template_tmpl;
    template_full_allIter(iIter, :, :) = template_full;

end % end for iIter

% fprintf('\n\n[L%d ModelA%d] ALL  iterations DONE\n', iLocComb, iModelA);

%% -------------------- SAVE -------------------- %%
save(nameFile_compIV, 'template_true', 'c_zscore*', '*_allIter', 'names*', 'flag*', 'ratio_split', 'ORI_bound', '*Type');
fprintf('\n========== Binned IV saved ==========\n\n\n\n\n');

%% -------------------- Plot (optional) -------------------- %%
if flag_plot_compIV
    NOMplot_compIV;
end
close all;

%% -------------------- End timing -------------------- %%
time_end = datetime('now')
elapsed = time_end - time_start;
fprintf('\n\nDONE (time used: %s)\n', char(elapsed));


end
