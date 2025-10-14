%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Script name: OOD_NOM_Trialwise_Est.m
% Script type: function
% Author: Shutian Xue
% Last updated: 08/27/2025

% Description:
% This script operate the following things:
% 1. fit NOM to IVs of the testing set (derived in OOD_xx_beforeEst) and estimate parameters
% 2. predict behav metrics using binned empirical IVs and estimated parameters
% Key functions: fxn_getError_v5(), PR_pred_v5()

% Inputs:
%      isubj: index of subject (for human subject) or a cell array with subject name and criterion (for IO)
%      iLocComb: for human subjects only: 1=Fovea, 8=periF(6 deg ecc), 6=HM, 7=VM, 5=LVM, 3=UVM
%      iModelA: 1=core model, 2=permuted template, 3=use IO template
%      iModelB: 1=only constant noise; 2=only induced noise; 3=no internal noise (only lapse rate and criterion); 4=constant noise + criterion; 5=induced noise + criterion
%      ni: number of iterations
% Outputs:
%      indicated by variable nameFile_Est: the estimated parameters and prediction metrics, saved in Data_NOM_Trialwise

% This function is part of the OOD_NOM_Trialwise pipeline.

function OOD_NOM_Trialwise_Est(isubj, iLocComb, iModelA, iModelB, ni, flag_fminconORbads)


close all, warning off, format compact
time_start = datetime('now')

addpath(genpath('fxn_NOM'))
addpath(genpath('fxn_RCplot'))
addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('SX_toolbox/bads-master'))

%% Set parameters
%--------------%
SX_RC1_setting
%--------------%
flag_plot = 1;
flag_plotPerIter = 0;

% Set upper and lower bounds of each parameter
SDadd_lb = 1e-5; SDadd_ub = 10;  % Additive noise
Nmul_lb = 1e-5; Nmul_ub = 1;  % Multiplicative noise
% SDadd_lb = -eps; SDadd_ub = eps;  % Additive noise
% Nmul_lb = -eps; Nmul_ub = eps;  % Multiplicative noise
SDadd0 = mean([SDadd_lb, SDadd_ub]);
Nmul0 = mean([Nmul_lb, Nmul_ub]);

% Set optimization options for BADS (Bayesian Adaptive Direct Search)
options_bads = bads('defaults'); % Load default settings
options_bads.Display = 'none';   % Alternatively, use this to suppress the output
options_bads.Verbosity = 0;

% Set optimization options for fmincon
options_fmin = optimoptions('fmincon', 'MaxIterations', 1e4, 'Display','off');

%% Set up file paths and names
if isnumeric(isubj)
    subjList =             {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC','SR'};
    nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195, 0];
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    nameFolder_NOM = sprintf('%s/ORI%dSF%d/%s/L%d', nameFolder_NOM0, nORI, nSF, subjName, iLocComb);
    % nameFolder_OOD_load = sprintf('%s/%s', nameFolder_Data_OOD, subjName); % To Load behav & energy
    nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocComb); % To Save results
else % IO
    subjName = isubj{1};
    criterion_true = isubj{2};
    nblocks=0;
    % nameFolder_NOM = sprintf('%s/ORI%dSF%d/IO/%s_c%.1f', nameFolder_NOM0, nORI, nSF, subjName, criterion_true);% To Save results
    nameFolder_OOD_load = sprintf('%s/%s', nameFolder_Data_OOD, subjName); % To load behav & energy
    nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, subjName); % To save results for IO
end

if isempty(dir(nameFolder_NOM_save)), mkdir(nameFolder_NOM_save), end

% Define names of files to save
nameFile_beforeEst = sprintf('%s/n%d_A%d_beforeEst', nameFolder_NOM_save, ni, iModelA);
nameFile_Est = sprintf('%s/n%d_A%dB%d', nameFolder_NOM_save, ni, iModelA, iModelB);

% Print information
fprintf('\n%s [nblocks = %d] [ORI%d SF%d]\n - ni = %d\n - Loc: %s\n - MODEL: [A%d] %s & [B%d] %s\n\n', ...
    subjName, nblocks, nORI, nSF, ni, namesLocComb{iLocComb}, iModelA, namesModelA{iModelA}, iModelB, namesModelB{iModelB})

%% Load data
load(nameFile_beforeEst, 'data_allB', 'c_zscore')

%% Set up the data-restricted upper and lower bound
switch iModelB
    case 1, params0 = [lapse0, SDadd0, c0]; params_lb = [lapse_lb, SDadd_lb, c_lb]; params_ub = [lapse_ub, SDadd_ub, c_ub];
    case 2, params0 = [lapse0, Nmul0, c0]; params_lb = [lapse_lb, Nmul_lb, c_lb]; params_ub = [lapse_ub, Nmul_ub, c_ub];
    case 3, params0 = [lapse0, c0]; params_lb=[lapse_lb, c_lb]; params_ub=[lapse_ub, c_ub];
    case 4, params0 = [SDadd0]; params_lb = [SDadd_lb]; params_ub = [SDadd_ub]; % not fitting lapse rate
    case 5, params0 = [Nmul0]; params_lb = [Nmul_lb]; params_ub = [Nmul_ub]; % not fitting lapse rate
end

%% Preallocate variables
nParams = length(namesParamsModel_all{iModelB});
nParams=1;
params_est_allB = nan(ni, nParams);
% params_est_allB_fmincon = params_est_allB;
nLL_allB = nan(ni,1);
pred_metrics_allB = cell(ni, 1);

%% Loop through each iteration
fprintf('  ======== ESTIMATION ======== \n[L%d ModelA%dB%d] Running ni = %d: ', iLocComb, iModelA, iModelB, ni)
for ii = 1:ni
    fprintf('%d ', ii)

    data = data_allB{ii};

    % Fit the model to trial-wise data using BADS optimizer
    % fxn_estParams = @(params) fxn_getError_v5(iModelB, params, data, 0);
    fxn_estParams = @(params) fxn_getError_v6(iModelB, params, data, c_zscore); % criterion is converted from c in zscore unit

    if flag_fminconORbads == 1
        % Use fmincon to estimate (faster)
        [params_est, nLL] = fmincon(fxn_estParams, params0, [], [], [], [], params_lb, params_ub, [], options_fmin);
    else % Use bads to estimate (slower)
        [params_est, nLL] = bads(fxn_estParams, params0, params_lb, params_ub, [], [], [], options_bads);
    end

    params_est_allB(ii, :) = params_est;
     nLL_allB(ii) = nLL;

    % Make predictions on binned IVs based on estimated parameters
    pred = PR_pred_v6(iModelB, nBins, params_est, data, c_zscore, flag_plotPerIter);
    pred_metrics_allB{ii} = pred;

end % end of ii

% pred.criterion_IV is an output of PR_pred_v6()

fprintf('\n\nALL iterations DONE\n')

if flag_plot, NOMplot_Est; end

%% Save results
% copy and rename the xx_beforeEst.mat file
copyfile([nameFile_beforeEst, '.mat'], [nameFile_Est, '.mat']);
save(nameFile_Est, '*_allB', '-append')

time_end = datetime('now')

time_end - time_start

fprintf('\n========== Prediction saved ==========\n\n\n\n\n')
