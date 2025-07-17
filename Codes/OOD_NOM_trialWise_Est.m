
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Script name: OOD_NOM_Trialwise_Est.m
% Script type: function
% Author: Shutian Xue
% Date created: 
% Last updated: 07/15/2025

% Description:
%   This function estimates parameters for the trial-wise Noisy Observer Model (NOM) using behavioral and energy data.
%   It supports both human subjects and ideal observer (IO) simulations, fitting various model variants to trial-wise data.
%   The function uses Bayesian Adaptive Direct Search (BADS) for parameter optimization and saves the estimated parameters and prediction metrics.
%
% Inputs:
%      isubj: index of subject (for human subject) or a cell array with subject name and criterion (for IO)
%      iLocComb: for human subjects only: 1=Fovea, 8=periF(6 deg ecc), 6=HM, 7=VM, 5=LVM, 3=UVM
%      iModelA: 1=core model, 2=permuted template, 3=use IO template
%      iModelB: 1=only constant noise; 2=only induced noise; 3=no internal noise (only lapse rate and criterion); 4=constant noise + criterion; 5=induced noise + criterion
%      ni: number of iterations     
% Outputs:
%      indicated by variable nameFile_Est: the estimated parameters and prediction metrics, saved in Data_NOM_Trialwise

% This function is part of the OOD_NOM_Trialwise pipeline.

function OOD_NOM_Trialwise_Est(isubj, iLocComb, iModelA, iModelB, ni)


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

% Set upper and lower bounds of each parameter
SDadd_lb = 1e-5;  SDadd_ub = 2;  % Additive noise
Nmul_lb = 1e-5; Nmul_ub = 2;  % Multiplicative noise

% Set optimization options for BADS (Bayesian Adaptive Direct Search)
options_bads = bads('defaults'); % Load default settings
options_bads.Display = 'none';   % Alternatively, use this to suppress the output
options_bads.Verbosity = 0;

% Set optimization options for fmincon
options_fmin = optimoptions('fmincon','MaxIterations', 1e4, 'Display','off');

%% Set up file paths and names
if isnumeric(isubj)
    subjList =             {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC','SR', 'IO'};
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
% nameFile_beforeEst = sprintf('%s/n%d_A%d_beforeEst', nameFolder_NOM, ni, iModelA);
% nameFile_Est = sprintf('%s/n%d_A%dB%d', nameFolder_NOM, ni, iModelA, iModelB);
nameFile_beforeEst = sprintf('%s/n%d_A%d_beforeEst', nameFolder_NOM_save, ni, iModelA);
nameFile_Est = sprintf('%s/n%d_A%dB%d', nameFolder_NOM_save, ni, iModelA, iModelB);

% Print information
fprintf('\n%s [nblocks = %d] [ORI%d SF%d]\n - ni = %d\n - Loc: %s\n - MODEL: [A%d] %s & [B%d] %s\n\n', ...
    subjName, nblocks, nORI, nSF, ni, namesLocComb{iLocComb}, iModelA, namesModelA{iModelA}, iModelB, namesModelB{iModelB})

%% MAIN LOOP
% Preallocate variables
nParams = length(namesParamsModel_all{iModelB});
params_est_allB = nan(ni, nParams);
params_est_allB_fmincon = params_est_allB;
nLL_allB = nan(ni,1);
pred_metrics_allB = cell(ni, 1);

fprintf('Running ni = %d: ', ni)

% Load data
load(nameFile_beforeEst, 'data_allB')

% Set up parameters upper and lower bound based on iModelB
c_lb = min(data_allB{1}.IV); c_ub = max(data_allB{1}.IV);
c0 = mean([c_lb, c_ub]);
SDadd0 = mean([SDadd_lb, SDadd_ub]);
Nmul0 = mean([Nmul_lb, Nmul_ub]);

% Compile parameter bounds
switch iModelB
    case 1, params0 = [lapse0, SDadd0, c0]; params_lb = [lapse_lb, SDadd_lb, c_lb]; params_ub = [lapse_ub, SDadd_ub, c_ub];
    case 2, params0 = [lapse0, Nmul0, c0]; params_lb = [lapse_lb, Nmul_lb, c_lb]; params_ub = [lapse_ub, Nmul_ub, c_ub];
    case 3, params0 = [lapse0, c0]; params_lb=[lapse_lb, c_lb]; params_ub=[lapse_ub, c_ub];
    case 4, params0 = [SDadd0, c0]; params_lb = [SDadd_lb, c_lb]; params_ub = [SDadd_ub, c_ub]; % not fitting lapse rate
    case 5, params0 = [Nmul0, c0]; params_lb = [Nmul_lb, c_lb]; params_ub = [Nmul_ub, c_ub]; % not fitting lapse rate
end

% Loop through each iteration
for ii = 1:ni

    data = data_allB{ii};

    % Fit the model to trial-wise data using BADS optimizer
    fxn_estParams = @(params) fxn_getError_v5(iModelB, params, data, 0);
    %     [params_est_fmincon, nLL] = fmincon(fxn_estParams, params0, [], [], [], [], params_lb, params_ub, [], options_fmin);
    [params_est, nLL] = bads(fxn_estParams, params0, params_lb, params_ub, [], [], [], options_bads);
    %     params_est_allB_fmincon(ii, :) = params_est_fmincon;
    params_est_allB(ii, :) = params_est;
    nLL_allB(ii) = nLL;

    % Make predictions on binned IVs based on estimated parameters
    pred = PR_pred_v5(iModelB, nBins, params_est, data, 0);
    pred_metrics_allB{ii} = pred;

    %     fprintf('%d ', ii)
end % end of ii

fprintf('\n\nALL iterations DONE\n')

%% quick plot
% pA_data_allB = nan(ni, nBins);
% pA_pred_med = pA_data_allB;
% IV_allB = pA_data_allB;
% figure, hold on
% for ii=1:ni
%     IV_allB(ii, :) = pred_metrics_allB{ii}.metrics.IV_allBins;
%     pA_data_allB(ii, :) = pred_metrics_allB{ii}.metrics.pA_data_allBins;
%     pA_pred_allB(ii, :) = pred_metrics_allB{ii}.metrics.pA_pred_allBins;
%     pC_data_allB(ii, :) = pred_metrics_allB{ii}.metrics.pC_data_allBins;
%     pC_pred_allB(ii, :) = pred_metrics_allB{ii}.metrics.pC_pred_allBins;
%
% %     plot(IV_allB(ii, :), pC_data_allB(ii, :), 'o')
% %     plot(IV_allB(ii, :), pC_pred_allB(ii, :), '-')
% %     ylim([.5, 1])
% %     pause
% end
% IV_med = getCI(IV_allB, 1, 1);
% pA_data_med = getCI(pC_data_allB, 1, 1);
% pA_pred_med = getCI(pC_pred_allB, 1, 1);
%
% figure, hold on
% plot(IV_med, pA_data_med, 'o')
% plot(IV_med, pA_pred_med, '-')

%% Save results
% copy and rename the xx_beforeEst.mat file
copyfile([nameFile_beforeEst, '.mat'], [nameFile_Est, '.mat']);
save(nameFile_Est, '*_allB', '-append')

time_end = datetime('now')

time_end - time_start

fprintf('\n========== Prediction saved ==========\n\n\n\n\n')
