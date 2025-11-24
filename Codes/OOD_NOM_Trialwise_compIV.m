function OOD_NOM_Trialwise_compIV(isubj, iLocComb, iModelA, nIterations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OOD_NOM_Trialwise_compIV.m
%
% Trial-wise noisy observer model – PRE-ESTIMATION STAGE
%
% This function:
%   1. Splits trials into training and test sets (via fxn_resampleTrials).
%   2. Estimates a template from the training set (using stim energy + resp).
%   3. Computes behavioral metrics in the test set
%      (these are the metrics the model will later try to match).
%   4. Computes, normalizes, and bins empirical internal variables (IVs)
%      for the test set.
%
% The output is a .mat file containing, for each iteration:
%   - data_allB{ii}: struct with IVs, binning info, responses, etc.
%   - data_metrics_allB(ii,:): behavioral metrics per iteration
%   - kernel2D_allB(ii,:,:): 2D templates (ORI x SF) per iteration
%
% INPUTS
%   isubj          : index of subject or IO info
%                    - numeric: index into subjList (human observers)
%                    - cell: {nameIO, ...} for ideal observer
%   iLocComb       : location combination index
%                    1=fovea; 8=periF (6 deg ecc);
%                    6=HM; 7=VM; 5=LVM; 3=UVM
%   iModelA        : template model
%                    1 = RC-derived template
%                    2 = IO template
%                    3 = permuted template
%   IVType         : how IV is computed from template × energy
%                    1 = sum of dot product / convolution
%                    2 = max
%                    3 = normalized, etc. (see fxn_getIV_v3)
%   templateType   : how to extract 1D template from 2D kernel
%                    1 = raw
%                    2 = reconstructed kernel
%                    3 = mirrored template
%   flag_PatchMode : 1 = use target patch energy, 2 = use noise patch energy
%   itype_template : which trials to use for estimating the template
%                    1 = PRS (signal-present) trials only
%                    2 = ABS (signal-absent) trials only
%                    3 = both PRS and ABS (not implemented here, falls through)
%   nIterations    : number of resampling iterations (bootstraps)
%
% OUTPUTS (saved to disk)
%   nameFile_compIV = '.../n%d_A%d_compIV.mat', containing:
%     - c_zscore          : SDT criterion (z units) per iteration
%     - data_allB         : cell array of "data" structs per iteration
%     - data_metrics_allB : behavioral metrics per iteration
%     - kernel2D_allB     : 2D kernels per iteration
%     - names*, flag*, ratio_train, ORI_bound, *Type, etc.
%
% Notes:
%   - This is the "before estimation" stage; model fits (iModelB) happen in
%     OOD_NOM_Trialwise_Est.m
%   - IVs are computed in fxn_getIV_v3, typically as template ⊙ energy.
%
% Created by Shutian Xue on August 27, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all;
warning off;              % (You may want to remove this once things are stable.)
format compact;
time_start = datetime('now')
rng(123);   % define see for reproducibility

addpath(genpath('fxn_exp'));
addpath(genpath('fxn_NOM'));
addpath(genpath('fxn_RCplot'));
addpath(genpath('fxn_analysis_RC_v2'));
addpath(genpath('SX_toolbox/bads-master'));

%% -------------------- General parameters -------------------- %%
% Global settings for RC / NOM analyses
% --------------%
SX_RC1_setting;   % defines nameFolder_*, nORI, nSF, namesLocComb, namesModelA, etc.
% --------------%
IVType = 1;            % 1=sum of the dot product/convolution; 2=max; 3=normalized
templateType = 3; % (1) raw (2) reconstructed kernel (3) mirrored template
itype_template = 2; % 1=estimate template from PRS trials, ABS trials, or BOTH trials
flag_PatchMode = 1; % if flag_PatchMode == 1, patchMode = 'T'; else, patchMode = 'N'; end

iSess_start        = 1;    % first session included
convolveType       = 1;    % IV from 1=cross-correlation; 2=convolution (fxn_getIV_v3)
flag_standEnergy   = 1;    % 1=z-score energy before RC
flag_plot_compIV = 1;   % plot IV distributions and kernels at the end
ratio_train        = 3/4;  % proportion of trials in training set (for RC)
ORI_bound          = [10, 20];  % orientation window, passed to fxn_getIV_v3
if flag_PatchMode == 1
    namePatchMode = 'T';   % target patch
else
    namePatchMode = 'N';   % noise patch
end

%% -------------------- Subject / IO and paths -------------------- %%
if isnumeric(isubj)  % Human subjects
    subjList = {'YK','SP','SX','LS','RE','MD','AS','HL','FH','HA','CS','DT','DU','RC','SR'};
    nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195, 0];

    subjName = subjList{isubj};
    nblocks  = nblocks_allSubj(isubj);

    % Folder to load behav + energy
    nameFolder_OOD_load = sprintf('%s/%s%d', nameFolder_Data_OOD, subjName, nblocks);
    
    % Folder to save NOM trial-wise results
    nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocComb);
    
    % Behavioral measurements
    load(sprintf('%s/%s%d/%s_behavMeas.mat', nameFolder_Data_OOD, subjName, nblocks, subjName), 'dataMatrix');  % ensure dataMatrix is loaded

    % Energy
    load(sprintf('%s/%s_energy_%s_%d_%d.mat', nameFolder_OOD_load, subjName, namePatchMode, nORI, nSF), ...
        'e3D_target_allT', 'e3D_noise_allT', 'noise', 'filtersOri_all', 'stim', 'nBins');

else  % Ideal observer (IO)
    subjName = isubj{1};  % "nameIO" from OOD_sim
    nblocks  = 0;

    % Folder to load behav + energy
    nameFolder_OOD_load = sprintf('%s/%s', nameFolder_Data_OOD, subjName);
    
    % Folder to save NOM trial-wise results
    nameFolder_NOM_save = sprintf('%s/%s',   nameFolder_Data_NOM_Trialwise, subjName);

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
nameFile_compIV = sprintf('%s/n%d_A%d_compIV', nameFolder_NOM_save, nIterations, iModelA);

%% Print run info
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
    nIterations, ...
    nblocks, nORI, nSF, ...
    namesType{itype_template}, flag_PatchMode, ...
    namesConvolveType{convolveType}, namesIVType{IVType});

%% -------------------- Choose energy source for NOM -------------------- %%
switch flag_PatchMode
    case 1  % target patch energy
        e3D_allT_train = e3D_target_allT;
        e3D_allT_test  = e3D_target_allT;
    case 2  % noise patch energy
        e3D_allT_train = e3D_noise_allT;
        e3D_allT_test  = e3D_noise_allT;
    otherwise
        error('flag_PatchMode must be 1 (target) or 2 (noise).');
end

%% -------------------- Sanity check: iPRS consistent within pairs -------------------- %%
% Col#1 = itrial_allT; Col#6 = iPRS; Col#8 = iPair;
for ipair = 1:max(dataMatrix(:, 8))
    ii = find(dataMatrix(:, 8) == ipair);
    iPRS_A = dataMatrix(dataMatrix(:, 1) == ii(1), 6);
    iPRS_B = dataMatrix(dataMatrix(:, 1) == ii(2), 6);
    if iPRS_A ~= iPRS_B
        fprintf('WARNING: when ipair=%d, iPRS_A=%d, iPRS_B=%d\n', ipair, iPRS_A, iPRS_B);
    end
end

%% -------------------- Session and location selection -------------------- %%
nSess = length(unique(dataMatrix(:, 2)));

% Assume equal #trials per location; use location==1 as representative.
ntrials_perSingleLoc = sum(dataMatrix(:, 5) == 1);

% Sessions included (can be subsetted here)
iSess_select = iSess_start:nSess;

% Special case: subject 'RE' missing some sessions
if strcmp(subjName, 'RE')
    iSess_select         = [1:15, 17:23, 25:nSess];
    ntrials_perSingleLoc = length(iSess_select) * 100;
end

% Map iLocComb to actual location indices
if iLocComb < 6
    iLoc_all = iLocComb;
else
    switch iLocComb
        case 6, iLoc_all = [2, 4];   % horizontal meridian
        case 7, iLoc_all = [3, 5];   % vertical meridian
        case 8, iLoc_all = 2:5;      % all perifoveal locations
        otherwise
            error('Unknown iLocComb = %d', iLocComb);
    end
end

%% -------------------- Train/test split per location -------------------- %%
ntrain_perSingle = ntrials_perSingleLoc * ratio_train;
ntest_perSingle  = ntrials_perSingleLoc * (1 - ratio_train);

% Ensure even #trials in TRAIN (so we keep full pairs)
if rem(ntrain_perSingle, 2)
    ntrain_perSingle = ntrain_perSingle + 1;
    ntest_perSingle  = ntest_perSingle - 1;
end

%% -------------------- Precompute non-random pools -------------------- %%
% These pools collect all relevant trials (across selected sessions and
% locations). Each iteration will only resample pairs from these pools.

dataMtx_nonrand   = [];
e3D_nonrand_train = [];
e3D_nonrand_test  = [];
iPRS_nonrand      = [];
iPair_nonrand     = [];
resp_nonrand      = [];
cst_nonrand       = [];

for iSS = 1:length(iSess_select)
    inSess = (dataMatrix(:, 2) == iSess_select(iSS));
    inLoc  = ismember(dataMatrix(:, 5), iLoc_all);  % supports single or multiple locations
    indLoc = inSess & inLoc;

    dataMtx_nonrand   = [dataMtx_nonrand;   dataMatrix(indLoc, :)];
    e3D_nonrand_train = cat(1, e3D_nonrand_train, e3D_allT_train(indLoc, :, :));
    e3D_nonrand_test  = cat(1, e3D_nonrand_test,  e3D_allT_test(indLoc, :, :));
    iPRS_nonrand      = [iPRS_nonrand;      dataMatrix(indLoc, 6)];
    iPair_nonrand     = [iPair_nonrand;     dataMatrix(indLoc, 8)];
    resp_nonrand      = [resp_nonrand;      dataMatrix(indLoc, 9)];
    cst_nonrand       = [cst_nonrand;       dataMatrix(indLoc,11)];
end

cst_unik_nonrand   = unique(cst_nonrand);      %#ok<NASGU> % (available if needed)
iPair_unik_nonrand = unique(iPair_nonrand);    %#ok<NASGU>

%% -------------------- Ideal template -------------------- %%
% For human observers, we still use this "ideal" Gabor-energy template when
% iModelA == 3 (IO template). Otherwise, templates are derived from RC.

gaborCST         = mean(cst_nonrand);
templateType_true = 1;  % doesn't matter much for the IO Gabor template

stim.gaborCST = gaborCST;

% Create a pool of Gabor filters
[filter_sin, filter_cos] = SX_sim02_setFilters(stim, noise.filtersSF_all, ...
    filtersOri_all, fxn_getSigma_SPdomain, 0);

% Define the TRUE template (ORI x SF)
template_gabor_true = exp_CreateGabor(stim, stim.gaborCST);
template_true       = SX_RC4_Energy_parfor(stim.mask, {template_gabor_true}, ...
    filter_sin, filter_cos);
template_true       = squeeze(template_true);               % remove singleton dim
template_true       = fxn_getTemplate(template_true, templateType_true, 0);

%% -------------------- MAIN LOOP over iterations -------------------- %%
nmetrics          = 9;
data_allB         = cell(nIterations, 1);
data_metrics_allB = nan(nIterations, nmetrics);
kernel2D_allB     = nan(nIterations, nORI, nSF);

fprintf('\n\nRunning ni = %d: ', nIterations);

for ii = 1:nIterations

    %% 1. Resample trials into TRAIN / TEST
    % This function uses the precomputed non-random pools above and creates:
    %   e3D_train, iPRS_train, resp_train, cst_train, respC_train
    %   e3D_test,  iPRS_test,  resp_test,  cst_test,  respC_test, iPair_test
    %----------------%
    fxn_resampleTrials;
    %----------------%

    %% 2. Select which trials to use to estimate the template (TRAIN set)
    switch itype_template
        case 1  % PRS trials only
            useIdx    = (iPRS_train == 1);
        case 2  % ABS trials only
            useIdx    = (iPRS_train == 0);
        otherwise
            % 3 = both PRS and ABS (here we just keep everything)
            useIdx    = true(size(iPRS_train));
    end

    e3D_train_sel = e3D_train(useIdx, :, :);
    cst_train_sel = cst_train(useIdx);
    resp_train_sel = resp_train(useIdx);
    iPRS_train_sel = iPRS_train(useIdx);

    % Standardize energy if requested
    if flag_standEnergy
        e3D_train_norm = normEnergy(e3D_train_sel, cst_train_sel, iPRS_train_sel);
    else
        e3D_train_norm = e3D_train_sel;
    end

    %% 3. Estimate 2D kernel (RC) from TRAIN set
    kernels2D = SX_sim07_RC(filtersSF_all, filtersOri_all, e3D_train_norm, resp_train_sel);

    % Normalize the derived template to a fixed amplitude range
    kernels2D = kernels2D / max(kernels2D(:)) * 0.2;  % match OOD_sim scaling

    % Store kernel for this iteration
    kernel2D_allB(ii, :, :) = kernels2D;

    %% 4. Compute test-set behavioral metrics (not binned by IV)
    pYES = mean(resp_test == 1);
    pHit = mean(iPRS_test == 1 & resp_test == 1) / mean(iPRS_test == 1);
    pFA  = mean(iPRS_test == 0 & resp_test == 1) / mean(iPRS_test == 0);
    pC   = mean( (iPRS_test == 0 & resp_test == 0) | (iPRS_test == 1 & resp_test == 1) );

    [dprime, criterion] = SX_sim06_SDT(pHit, pFA);

    % metrics_test:
    %   1 = d'
    %   2 = criterion (c)
    %   3 = pC
    %   4 = pHit
    %   5 = pFA
    %   6 = pA (mean consistency across all pairs)
    %   7–8 = (reserved for pA_PRS, pA_ABS if desired)
    %   9 = pYES
    metrics_test = [dprime, criterion, [pC, pHit, pFA], nanmean(respC_test), pYES];

    % Save criterion in z units; later used to convert to criterion in IV
    c_zscore = criterion;

    clear data;  % size of IV differs across subjects and iterations

    %% 5. Compute internal variable (IV) from template and energy (TEST set)
    % 5.1 Derive template
    if iModelA == 1 % Use RC-derived template
        flag_plot_template = 0;
        template = fxn_getTemplate(kernels2D, templateType, flag_plot_template);
    else % Use the ideal template (Gabor energy profile)
        if ~isnumeric(isubj)  % for IO, you may reload 'template_true' from disk
            load(sprintf('%s/truth.mat', nameFolder_OOD_load), 'template_true');
        end
        template = template_true;
    end

    % 5.2 Compute IV for test trials
    %---------------%
    IV_test = fxn_getIV_v3(iModelA, e3D_test, convolveType, IVType, template, ORI_bound);
    %---------------%
    ndata = length(IV_test);

    %% 6. Bin IV values
    % Separately for PRS and ABS trials
    [nTrials_PRS_allBins, ~, iTrial4Bin_PRS] = histcounts(IV_test(iPRS_test == 1), nBins);
    [nTrials_ABS_allBins, ~, iTrial4Bin_ABS] = histcounts(IV_test(iPRS_test == 0), nBins);

    % For all trials combined
    [nTrials_allBins, ~, iTrial4Bin] = histcounts(IV_test, nBins);

    %% 7. Compile data struct for this iteration
    data.ndata              = ndata;
    data.IV                 = IV_test;

    data.nTrials_allBins    = nTrials_allBins;
    data.iTrial4Bin         = iTrial4Bin;
    data.nTrials_PRS_allBins = nTrials_PRS_allBins;
    data.iTrial4Bin_PRS     = iTrial4Bin_PRS;
    data.nTrials_ABS_allBins = nTrials_ABS_allBins;
    data.iTrial4Bin_ABS     = iTrial4Bin_ABS;

    data.iPRS               = iPRS_test;
    data.resp               = resp_test;
    data.cst                = cst_test;
    data.iPair              = iPair_test;
    data.respC              = respC_test;

    data.metrics_sim        = metrics_test;  % keep field name 'metrics_sim'
    data.template           = template;

    data_allB{ii}           = data;
    data_metrics_allB(ii,:) = metrics_test;

    fprintf('%d ', ii);
end  % end for ii

% fprintf('\n\n[L%d ModelA%d] ALL iterations DONE\n', iLocComb, iModelA);

%% -------------------- SAVE -------------------- %%
save(nameFile_compIV, 'template_true', 'c_zscore', '*_allB', 'names*', 'flag*', 'ratio_train', 'ORI_bound', '*Type');
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
