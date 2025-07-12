function OOD_MC_current(isubj, iiLoc_all_all, flag_PatchMode, flag_cutMapping, flag_mirrorMapping, ifeature, ifamily, MCmode)

% Purpose: run model comparison to decide which fxn captures the ORI/SF kernel better

% Inputs:
% isubj
% iiLoc_all_all: which two loc to be fitted together
% flag_PatchMode: 1=energy derived from target, 2=from noise patches (default)
% flag_cutMapping: 1=cut extremes; 0=NOT
% flag_mirrorMapping: 1=mirror 2D kernel matrix over ORI=0
% ifeature: 1=ORI, 2=SF
% ifamily: 1=Gaussian, 8=DoG // 2=log parabola, 3=truncated log parabola
% MCmodel: see var 'titlesMCmode' defined below

close all
clc

addpath(genpath('Data_MC'))
addpath(genpath('Data_OOD'))
addpath(genpath('fxn_MC'))
addpath(genpath('fxn_analysis_RC_v2'))

time_start = datetime('now')

% 10-fold CV: ~30 mins
% LOOCV: ~20 mins
% IC: ~1 min

%% general params
%--------------%
SX_RC1_setting
%--------------%

nORI = 29; nSF = 29;
iType = 2;% extract kernels derived from target-abs trials
flagPlot = 0;
flag_block200 = 0; % 1; all observers are forced to have 200 blocks; 0=no

if flag_PatchMode == 1, namePatchMode = 'T'; else, namePatchMode = 'N'; end

iLoc_all_all = [1,8; 6,7; 5,3; 2,4];
iLoc_all = iLoc_all_all(iiLoc_all_all, :); nLoc2 = length(iLoc_all);
switch MCmode, case 1, ncv = 10; case 2, ncv = 2; end % number of rounds of cross validation
nrep = 20; % number of reps to find the global minimum
k = 10; % k-fold cross validation

subjList =             {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC', 'SR'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195];

subjName = subjList{isubj};
nblocks = nblocks_allSubj(isubj);

% calculation of IC comes from Fernandez2022
getBIC = @(nLL, nparamsAll, ndata) ndata * log(nLL) + nparamsAll * log(ndata);
getAIC = @(err, nparamsAll, ndata) ndata * log(err)+ 2*nparamsAll;
getAICc = @(err, nparamsAll, ndata) ndata * log(err)+ 2*nparamsAll + 2*nparamsAll*(nparamsAll+1)/(ndata-nparamsAll-1);

options = optimoptions('fmincon','MaxIterations',5000,'Display','off');

titlesMCmode = {'10-fold cross validation (deviance is nLL)', ...
    'Leave-One-Out CV (deviance is nLL)', ...
    'Information criterion (deviance is nLL)'};
namesYLabel = {'Deviance', 'Deviance', 'IC'};

namesIC = {'AIC', 'AICc', 'BIC'};
nIC = length(namesIC);

fprintf('\n==========================\n')
fprintf('%s\n    >> %s vs. %s (%s)\n    >> cut = %d (1=cut the edge)\n    >> mirror=%d (1=mirror over ORI=0)\n    >> [%s] - Family%d - %s\n    >> %s\n', ...
    subjName, namesLocComb{iLoc_all(1)}, namesLocComb{iLoc_all(2)}, namePatchMode, ...
    flag_cutMapping, flag_mirrorMapping, ...
    namesFeature{ifeature}, ifamily, namesFamily_all{ifamily}, ...
    titlesMCmode{MCmode})
fprintf('============================\n')

%% model setting
% names of params
namesParams = namesParams_all{ifamily};
nparams_full = length(namesParams);

% decide number of models & the param index (1 = specified, 0 = shared)
nCands = 2^nparams_full;
paramInd_all = fxn_getParamInd(nparams_full);

% extract ub and lb
ub_full = ub_full_all{ifamily};
lb_full = lb_full_all{ifamily};

%% extract x (y is extracted for each subj)
if (ifeature == 2) && (~sum(ifamily == [4, 5, 9])) % fitting SF except (4) raised gaussian (5) double exponential (9) gaussian to SF
    axis_tuning_ln = 2.^axis_tuning{ifeature};
else, axis_tuning_ln = axis_tuning{ifeature};
end
nfilters = length(axis_tuning_ln);

%% name of the folder and files
if flag_mirrorMapping, nameM = '_m'; else, nameM = ''; end
if flag_cutMapping, nameC = 'c'; else, nameC = ''; end

nameFolder_MC = sprintf('Data_MC/%d%d%s%s/%s/L%d%d', nORI, nSF, nameM, nameC, subjList{isubj}, iLoc_all);
dirFolder_MC = dir(nameFolder_MC); if isempty(dirFolder_MC), mkdir(nameFolder_MC), end
nameFile_MC = sprintf('%s/Family%d_mode%d.mat', nameFolder_MC, ifamily, MCmode);

if MCmode == 3, nIC_ = nIC; else, nIC_ = 1; end % number of information criterion

%% load kernels
nameKernel = sprintf('Data_OOD/%s%d/%s_kernel_%s_%d_%d.mat', subjName, nblocks, subjName, namePatchMode, nORI, nSF);
fprintf('  Loading kernels...'), tic
load(nameKernel, 'kernels2D_perComb')
dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)

%% extract y (i.e., kernels)
if ifeature==1, yData = nan(nLoc2, nORI); else, yData = nan(nLoc2, nSF); end
for iiLoc = 1:nLoc2
    kk = squeeze(kernels2D_perComb{iLoc_all(iiLoc)}(iType, :, :));
    % cut
    if flag_cutMapping, kk = kk(cut_ORI, cut_SF); end
    % mirror
    if flag_mirrorMapping
        indMir = (nORI-1)/2;
        kk_left = kk(1:indMir, :);
        kk_right = kk(indMir+2:end, :);
        kk_mid = kk(indMir+1,:);
        kk_ave = (kk_left + flip(kk_right))/2;
        kk_mir = [kk_ave; kk_mid; flip(kk_ave)];
        kk = kk_mir;
    end
    if ifeature==1
        yData(iiLoc, :) = mean(kk, 2);
    else
        yData(iiLoc, :) = mean(kk, 1);
    end
end

%% empty folders
irank_allCand = nan(nIC_, nCands); % the sorted model index (best to worst model);
iBest_allCand = nan(nIC_, 1);
dev_allCand = nan(nCands, 3);
params_est_allCand = cell(1, nCands);

%% loop through each candidate model
for iCand = 1:nCands % loop thru model candiadates // cannot and no need to 'parfor', because fitting used parfor
    fprintf('M#%d/%d...\n', iCand, nCands)
    tic
    paramInd = paramInd_all(iCand, :);
    nparamsAll = sum(nLoc2.^paramInd);
    [ub, lb, params0] = MC_getLimits(paramInd, ub_full, lb_full, nLoc2);
    
    switch MCmode
        case 1 % 10-fold Cross validation
            [param_est_allK, dev_test_allK, y_pred_allK] = fxn_runCV(MCmode, k, ncv, nfilters, yData, paramInd, ifamily, options, params0, lb, ub, nrep, axis_tuning_ln);
            dev_allCand(iCand, 1) = nanmedian(dev_test_allK);
            dev_allCand(iCand, 2) = nanstd(dev_test_allK);
            params_est_allCand{iCand} = nanmean(param_est_allK, 1);
            
        case 2 % leave-one-out CV
            [param_est_allK, dev_test_allK, y_pred_allK] = fxn_runCV(MCmode, nfilters, ncv, nfilters, yData, paramInd, ifamily, options, params0, lb, ub, nrep, axis_tuning_ln);
            dev_allCand(iCand, 1) = nanmedian(dev_test_allK);
            dev_allCand(iCand, 2) = nanstd(dev_test_allK);
            params_est_allCand{iCand} = nanmean(param_est_allK, 1);
            
        case 3 % IC (AIC, AICc and BIC, based on nLL)
            fitModel = 2; % 1: dev is RSS; 2: dev is nLL
            fxn_getDev = @(kernelParams) MC_getDev_CV(paramInd, kernelParams, axis_tuning_ln, yData, ifamily, fitModel);
            problem = createOptimProblem('fmincon','objective', fxn_getDev,'x0',params0,'lb',lb,'ub',ub,'options',options);
            ms = MultiStart('StartPointsToRun', 'bounds', 'UseParallel', 1, 'Display', 'off');
            [params_est, dev] = run(ms, problem, nrep);
            params_est_allCand{iCand} = params_est;
            % predict
            y_data_pred = MC_predKernel(paramInd, axis_tuning_ln, 2, nparams_full, params_est, ifamily);
            % compile
            dev_allCand(iCand, :) = [getAIC(dev, nparamsAll, nfilters*nLoc2), getAICc(dev, nparamsAll, nfilters*nLoc2), getBIC(dev, nparamsAll, nfilters*nLoc2)];
    end % switch imodel
    
    % quick plot of data and prediction of each candidate model
    %         if imodel==7, MCplot_IDVDCandidate, end
end % end of imodel

% mode1: 44 sec per model
% mode2: 22 sec per model
% mode 3: .5 sec

%% compile outputs
indCands = 1:nCands;
if flagPlot == 1
    figure('Position', [0 0 400 * nIC_, 300]), hold on
end

for iIC = 1:nIC_
    [~, irank] = sort(squeeze(dev_allCand(:, iIC))); % sort based on averaged dev (for CV) for IC
    indBEST = irank(1);
    irank_allCand(iIC, :) = irank;
    iBest_allCand(iIC) = irank(1);
    
    if flagPlot == 1
        % plot the dev/IC of each model ina ranked way of the current subj
        subplot(1,nIC_, iIC), hold on
        dev_ = squeeze(dev_allCand(irank, iIC));
        dev_delta = dev_-min(dev_);
        bar(indCands, dev_delta)
        xticks(indCands)
        xticklabels(indCands(irank))
        xlabel('Candidate model #')
        if MCmode == 3, title(sprintf('%s-%s\n The best model is #%d [%s]', subjName, namesIC{iIC}, indBEST, num2str(paramInd_all(indBEST, :))))
        else, title(sprintf('%s\n The best model is #%d [%s]', subjName, indBEST, num2str(paramInd_all(indBEST, :))))
        end
    end % if flagPlot
end % iIC

% plot averaged dev and freq of model chosen to be the best
% if flagPlot == 1, MCplot_perMCmode, end

%% save
save(nameFile_MC, 'ifeature', 'ifamily',  'MCmode', 'paramInd_all', '*_allCand')

%%
time_end = datetime('now')

time_end - time_startw