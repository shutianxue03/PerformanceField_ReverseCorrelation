%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OOD_NOM_Trialwise_beforeEst.m
% -------------------------------------------------------------------------
% Created by Shutian Xue on August 27, 2025
% This script operates the following things:
% 1. split trials into training and testing set
% 2. estimate template from the training set (based on simulated stim energy and simulated response)
% 3. calculate behav metrics of the testing set (ie., metrics to be matched by model pred)
% 4. calculate, normalize, and bin empirical IVs from the testing set (in fxn_getIV3, by the settings we have, just a simple dot multiplication)
%
% Inputs:
%     - isubj: subject index or IO info
%     - iLocComb: location combination index
%     - iModelA: model type (1=core, 2=permuted, 3=IO template)
%     - ni: number of iterations
%
% Outputs:
%   - Saves results (data_allB, data_metrics_allB, etc.) to a MAT file
%   - Optionally visualizes IV distributions and metrics if flag_plotIVsDist is set

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OOD_NOM_Trialwise_beforeEst(isubj, iLocComb, iModelA, IVType, ni, flag_logIV)

% The previous version is in OOD_NOM_trialWise
% INPUT
%      isubj: index of subj
%      iLocComb: 1=Fovea, 8=periF(6 deg ecc), 6=HM, 7=VM, 5=LVM, 3=UVM
%      iModelA: 1=core model, 2=permuted template, 3=use IO template
%      ni: number of iterations

close all, warning off, format compact
time_start = datetime('now')

addpath(genpath('fxn_NOM'))
addpath(genpath('fxn_RCplot'))
addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('SX_toolbox/bads-master'))

%% general params
%--------------%
SX_RC1_setting
%--------------%
% flag_logIV = 1;

iSess_start = 1; % from which session data is taken into account
templateType = 1; % (1) raw (2) reconstructed kernel (3) mirrored template
itype_template = 2; % 1=estimate template from PRS trials, ABS trials, or BOTH trials
convolveType = 1; % IV is calcluated by 1=cross correlation of template and stim energy; 2=convolution [fxn_getIV_v3]
% IVType = 1; % 1=sum, 2=the max [fxn_getIV_v3]
flag_PatchMode = 2; %1=use target patches; 2=use noise patches
flag_standEnergy = 1; % 1=standardize (z-score) energy
flag_plot_beforeEst = 1;
ratio_train = 3/4; % proportion of trials in the training set (to derive 2D kernels); others are in the test set (to estimate params)
ORI_bound = [10, 20];

%% Set names of directories
if isnumeric(isubj) % Human subjects
    subjList =             {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC','SR'};
    nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195, 0];
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    nameFolder_OOD_load = sprintf('%s/%s', nameFolder_Data_OOD, subjName); % To Load behav & energy
    nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocComb); % To Save results
else % IO
    subjName = isubj{1};
    criterion_true = isubj{2};
    nblocks=0;
    nameFolder_OOD_load = sprintf('%s/%s', nameFolder_Data_OOD, subjName);% To Load behav & energy
    nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, subjName);% To Save results
end
if isempty(dir(nameFolder_NOM_save)), mkdir(nameFolder_NOM_save), end

% Set names of directories for saving data
nameFile_beforeEst = sprintf('%s/n%d_A%d_beforeEst', nameFolder_NOM_save, ni, iModelA);

% PRINT
fprintf('\nSubject/IO name: %s [nblocks = %d] [nORI=%d | nSF=%d]\n - Number of iterations (ni) = %d\n - Location: %s\n - A%d [%s]\n - Template derived from %s trials\n - Energy source: %d (1=TARGET, 2=NOISE)\n - Convolve type: %s\n - IV type: %s\n', ...
    subjName, nblocks, nORI, nSF, ni, namesLocComb{iLocComb}, iModelA, namesModelA{iModelA}, namesType{itype_template}, flag_PatchMode, namesConvolveType{convolveType}, namesIVType{IVType})

%% Load
if flag_PatchMode == 1, namePatchMode = 'T'; else, namePatchMode = 'N'; end
if isnumeric(isubj)
    % behavioral measurement
    load(sprintf('%s/%s%d/%s_behavMeas.mat', nameFolder_Data_OOD, subjName, nblocks, subjName));
    % energy
    load(sprintf('%s/%s%d/%s_energy_%s_%d_%d.mat', nameFolder_OOD_load, subjName, nblocks, subjName, namePatchMode, nORI, nSF));
else
    % behavioral measurement
    load(sprintf('%s/behavMeas.mat', nameFolder_OOD_load));
    % energy
    load(sprintf('%s/energy_%s_%d_%d.mat', nameFolder_OOD_load, namePatchMode, nORI, nSF));
end

fprintf('\nBehav and energy LOADED\n\n')

% Determine the type of energy to use based on flag_PatchMode
switch flag_PatchMode
    case 1
        e3D_allT_train = e3D_target_allT; % e3D_allT_train is used in fxn_resampleTrials
        e3D_allT_test = e3D_target_allT; % e3D_allT_test is used in fxn_resampleTrials
    case 2
        e3D_allT_train = e3D_noise_allT; % e3D_allT_train is used in fxn_resampleTrials
        e3D_allT_test = e3D_noise_allT; % e3D_all
end

% Check if the iPRS of each pair are the same
for ipair = 1:max(dataMatrix(:, 8))
    ii = find(dataMatrix(:, 8)==ipair); % Col#8 is ipair
    iPRS_A = dataMatrix(dataMatrix(:, 1)==ii(1), 6); % % Col#1 is itrial_allT, Col#6 is iPRS
    iPRS_B = dataMatrix(dataMatrix(:, 1)==ii(2), 6);
    if iPRS_A ~= iPRS_A, fprintf('when ipair=%d, iPrsA=%d iPrsB=%d\n', ii, iPRS_A, iPRS_B), end
end

%%
% Define session and trial parameters for selected location
nSess = length(unique(dataMatrix(:, 2)));
ntrials_perSingleLoc = sum(dataMatrix(:, 5)==1);
iSess_select = iSess_start:nSess; % the index of selected sessions, used in fxn_resampleTrials

% Handle special cases for subject 'RE' due to missing sessions
if strcmp(subjName, 'RE'), iSess_select = [1:15, 17:23, 25:nSess]; ntrials_perSingleLoc = length(iSess_select)*100; end

% Assign location combination based on iLocComb input
if iLocComb < 6, iLoc_all = iLocComb;
else, switch iLocComb , case 6, iLoc_all = [2,4]; case 7, iLoc_all = [3,5]; case 8, iLoc_all = 2:5; end
end

% Set number of trials for training and test groups
ntrain_perSingle = ntrials_perSingleLoc * ratio_train;
ntest_perSingle = ntrials_perSingleLoc * (1-ratio_train);
if rem(ntrain_perSingle, 2), ntrain_perSingle = ntrain_perSingle+1;ntest_perSingle = ntest_perSingle-1; end % Ensure even split between training and test groups

%% MAIN LOOP
% preallocate
nmetrics = 9; 
data_allB = cell(ni, 1);
data_metrics_allB = nan(ni, nmetrics);
kernel2D_all = nan(ni, nORI, nSF);
fprintf('[L%d ModelA%d] Running ni = %d: ', iLocComb, iModelA, ni)

for ii = 1:ni
    fprintf('%d ', ii)
    % Resample trials for current iteration
    %----------------%
    fxn_resampleTrials
    %----------------%

    % standardize the energy for the training set
    switch itype_template
        case 1
            e3D_train = e3D_train(iPRS_train, :, :);
            cst_train = cst_train(iPRS_train);
            resp_train = resp_train(iPRS_train);
            iPRS_train = iPRS_train(iPRS_train);
        case 2
            e3D_train = e3D_train(iPRS_train==0, :, :); % only use signal-ABS trials to derive the template
            cst_train = cst_train(iPRS_train==0);
            resp_train = resp_train(iPRS_train==0);
            iPRS_train = iPRS_train(iPRS_train==0);
    end

    if flag_standEnergy
        e3D_train_norm = normEnergy(e3D_train, cst_train, iPRS_train);
    else
        e3D_train_norm = e3D_train;
    end

    % Calculate 2D kernel from the TRAINING set
    kernel2D = SX_sim07_RC(filtersSF_all, filtersOri_all, e3D_train_norm, resp_train); % 5 seconds

    % Normalize the derived template
    kernel2D = kernel2D / max(kernel2D(:)) * 0.2; % scale the signal energy to roughly match the range of templates derived from subbj data

    % size of kernel2D: nORI x nSF
    kernel2D_all(ii, :, :) = kernel2D; % the ave will be plotted in NOMplot_beforeEst.m

    % Calculate performance metrics for the TEST set
    pYES = mean(resp_test==1);
    pHit = mean(iPRS_test==1 & resp_test==1)*2;
    pFA = mean(iPRS_test==0 & resp_test==1)*2;
    pC = (pHit+1-pFA)/2;
    [d,c] = SX_sim06_SDT(pHit, pFA);
    metrics_test = [d,c, [pC, pHit, pFA], nanmean(respC_test), pYES]; % dprime, criterion, pC, pHit, pFA, pA, pA_PRS, pA_ABS, pYES

    % Later, this c (in unit of z-score) will be converted to criterion in
    % unit of IV on each trial, so that criterion is no longer a parameter
    c_zscore = c; 
    % c_IV = IV+c_zscore*sqrt(SDadd^2+(IV*Nmul).^2);

    clear data % as the size of IV differs across subjects

    %%%%%%%%%%%%%%%%%%%%%%%%
    %  calculate internal variable
    %%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Derive the Template from the TRAINING set, or just the energy profile of the gabor 
    if iModelA== 3 % use IO template (the energy profile of the signal)
        load(sprintf('%s/signalEnergy.mat', nameFolder_Data_OOD), 'template_true'); % nameFolder_Data is created in SX_RC1_setting
        template = template_true;
    else
        flag_plot = 0;
        template = fxn_getTemplate(kernel2D, templateType, flag_plot);
    end

    % 2. Calculate the Internal variable (IV)
    IV = fxn_getIV_v3(iModelA, e3D_test, convolveType, IVType, template, ORI_bound);
    % 
    % % Normalize IV
    % IV = (IV - mean(IV)) / std(IV);

    % Take log form if needed
    if flag_logIV, IV = log10(IV); end
    ndata = length(IV);

    % Bin IV values
    [nTrials_allBins, ~, iTrial4Bin] = histcounts(IV, nBins);

    % Compile data for the current iteration
    data.ndata = ndata;
    data.IV=IV;
    data.nTrials_allBins = nTrials_allBins;
    data.iTrial4Bin = iTrial4Bin;
    data.iPRS = iPRS_test;
    data.resp = resp_test;
    data.cst = cst_test;
    data.iPair = iPair_test;
    data.respC = respC_test;
    data.metrics_sim = metrics_test; % keep the field name 'metrics_sim'!!
    data_allB{ii} = data;
    data_metrics_allB(ii, :) = metrics_test;

end % end of ii

fprintf('\n\n[L%d ModelA%d] ALL iterations DONE\n', iLocComb, iModelA)

% SAVE
save(nameFile_beforeEst, 'c_zscore', '*_allB', 'names*', 'flag*', 'ratio_train', 'ORI_bound', '*Type')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_plot_beforeEst, NOMplot_beforeEst, end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_end = datetime('now')

time_end - time_start

fprintf('\n========== Binned IV saved ==========\n\n\n\n\n')
