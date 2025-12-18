function OOD_NOM_Trialwise_compIV_temp(isubj, iLocComb, iModelA, nBoot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OOD_NOM_Trialwise_compIV.m
%
% Trial-wise noisy observer model – PRE-ESTIMATION STAGE
%
% This function:
% 1. Splits trials into training and test sets (via fxn_resampleTrials).
% 2. Estimates a template from the training set (using stim energy + resp).
% 3. Computes behavioral metrics in the test set
% (these are the metrics the model will later try to match).
% 4. Computes, normalizes, and bins empirical internal variables (IVs)
% for the test set.
%
% The output is a .mat file containing, for each  bootstrap:
% - data_allB{ii}: struct with IVs, binning info, responses, etc.
% - data_metrics_allB(ii,:): behavioral metrics per  bootstrap
% - kernel2D_allB(ii,:,:): 2D templates (ORI x SF) per  bootstrap
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
% nBoot : number of resampling bootstraps
%
% OUTPUTS (saved to disk)
% nameFile_compIV = '.../n%d_A%d_compIV.mat', containing:
% - c_zscore : SDT criterion (z units) per  bootstrap
% - data_allBoot : cell array of "data" structs per  bootstrap
% - data_metrics_allBoot : behavioral metrics per  bootstrap
% - Template_train_allBoot : 2D kernels per  bootstrap (using trials in the training set)
% - Template_full_allBoot : 2D kernels per  bootstrap (using all trials)
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
rng(123); % define see for reproducibility

addpath(genpath('fxn_exp'));
addpath(genpath('fxn_NOM'));
addpath(genpath('fxn_RCplot'));
addpath(genpath('fxn_analysis_RC_v2'));
addpath(genpath('SX_toolbox/bads-master'));

%% -------------------- General parameters -------------------- %%
% Global settings for RC / NOM analyses
% --------------%
SX_RC1_setting; % defines nameFolder_*, nORI, nSF, namesLocComb, namesModelA, etc.
% --------------%
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
ratio_train = 3/4; % proportion of trials in training set (for RC)
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
nameFile_compIV = sprintf('%s/n%d_A%d_compIV', nameFolder_NOM_save, nBoot, iModelA);

%% Print run info
fprintf('\n==================================\nStep 1: Compute IVs and derive templates\n================================== \n\n')
fprintf(['\nSubject/IO name: %s ' ...
    '\n - L%d [%s]', ...
    '\n - A%d [%s]', ...
    '\n - Number of bootstraps = %d', ...
    '\n - Number of blocks = %d', ...
    '\n - nORI=%d | nSF=%d', ...
    '\n - Template derived from %s trials', ...
    '\n - Energy source: %d (1=TARGET, 2=NOISE)', ...
    '\n - Convolve type: %s', ...
    '\n - IV type: %s\n\n'], ...
    subjName, ...
    iLocComb, namesLocComb{iLocComb}, ...
    iModelA, namesModelA{iModelA}, ...
    nBoot, ...
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
    iBoot = find(dataMatrix(:, 8) == ipair);
    iPRS_A = dataMatrix(dataMatrix(:, 1) == iBoot(1), 6);
    iPRS_B = dataMatrix(dataMatrix(:, 1) == iBoot(2), 6);
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
    iSess_select = [1:15, 17:23, 25:nSess];
    ntrials_perSingleLoc = length(iSess_select) * 100;
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

%% -------------------- Train/test split per location -------------------- %%
ntrain_perSingle = ntrials_perSingleLoc * ratio_train;
ntest_perSingle = ntrials_perSingleLoc * (1 - ratio_train);

% Ensure even #trials in TRAIN (so we keep full pairs)
if rem(ntrain_perSingle, 2)
    ntrain_perSingle = ntrain_perSingle + 1;
    ntest_perSingle = ntest_perSingle - 1;
end


%% -------------------- MAIN LOOP over  bootstraps -------------------- %%
nMetrics = 11; % see fxn_getMetrics
data_metrics_allBoot = nan(nBoot, 2, nMetrics); %for 3: 1=train, 2=test

fprintf('\n\nRunning nBoot = %d: ', nBoot);

for iBoot = 1:nBoot

    fprintf('%d... ', iBoot);
    
    %% 1. Resample trials into FULL or TRAIN / TEST
    % For combined locations, trials are resampled within each location, and DOWNsampled (to equate number of pairs per loc, like for fovea (1 single loc) vs. perifovea (4 single loc)).
    %----------------%
    fxn_resampleTrials;
    %----------------%

    %% 2. Select which trials to use to estimate the template (TRAIN set)
    switch itype_template
        case 1 % PRS trials only
            useIdx_train = (iPRS_train_rand == 1);
            % useIdx_full = (iPRS_full_rand == 1);
            % useIdx_test = (iPRS_test_rand == 1);
        case 2 % ABS trials only
            useIdx_train = (iPRS_train_rand == 0);
            % useIdx_full = (iPRS_full_rand == 0);
            % useIdx_test = (iPRS_test_rand == 0);
        otherwise % 3 = both PRS and ABS (here we just keep everything)
            useIdx_train = true(size(iPRS_train_rand));
            % useIdx_full = true(size(iPRS_full_rand));
            % useIdx_test = true(size(iPRS_test_rand));
    end

    % sel for "selected"
    % e3D_train_rand_sel = e3D_train_rand(useIdx_train, :, :);
    cst_train_rand_sel = cst_train_rand(useIdx_train);
    resp_train_rand_sel = resp_train_rand(useIdx_train);
    respC_train_rand_sel = respC_train_rand(useIdx_train);
    iPRS_train_rand_sel = iPRS_train_rand(useIdx_train);
    RT_train_rand_sel = RT_train_rand(useIdx_train);

    %% 4. Compute test-set behavioral metrics (not binned by IV)
    metrics_full = fxn_getMetrics(resp_full_rand, iPRS_full_rand, respC_full_rand, cst_full_rand, RT_full_rand);
    metrics_test = fxn_getMetrics(resp_test_rand, iPRS_test_rand, respC_test_rand, cst_test_rand, RT_test_rand);
    % metrics_train = fxn_getMetrics(resp_train_rand_sel, iPRS_train_rand_sel, respC_train_rand_sel, cst_train_rand_sel, RT_train_rand_sel);
    metrics = [metrics_full; metrics_test];

    %% 7. Compile data struct for this  bootstrap
    data_metrics_allBoot(iBoot, :, :) = metrics;

end % end for ii

% fprintf('\n\n[L%d ModelA%d] ALL  bootstraps DONE\n', iLocComb, iModelA);

%% -------------------- SAVE -------------------- %%
save(nameFile_compIV, 'data_metrics_allBoot', '-append');
fprintf('\n========== Binned IV saved ==========\n\n\n\n\n');


%% -------------------- End timing -------------------- %%
time_end = datetime('now')
elapsed = time_end - time_start;
fprintf('\n\nDONE (time used: %s)\n', char(elapsed));


end
