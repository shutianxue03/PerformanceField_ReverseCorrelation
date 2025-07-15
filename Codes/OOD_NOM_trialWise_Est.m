function OOD_NOM_Trialwise_Est(isubj, iLocComb, iModelA, iModelB, ni, nameFolder_NOM0, nORI)

% To estimate parameters for trial-wise NOM 

% INPUT
%      isubj: index of subj
%      iLocComb: 1=Fovea, 8=periF(6 deg ecc), 6=HM, 7=VM, 5=LVM, 3=UVM
%      ni: number of iterations
%      iModelA: 1=core model, 2=permuted template, 3=use IO template
%      iModelB: 1= only constant noise, 2= only induced noise; 3=no internal noise (only lapse rate and criterion)
%      ni: number of iterations

close all, warning off, format compact
time_start = datetime('now')

addpath(genpath('Data_OOD'))
addpath(genpath('Data_model'))
addpath(genpath('fxn_NOM'))
addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('SX_toolbox/bads-master'))

%% General params
%--------------%
SX_RC1_setting
%--------------%

SDadd_lb = 1e-5; SDadd_ub = 2; Nmul_lb = 1e-5; Nmul_ub = 2;  % ub and lb of criterion see below

% Set optimization options for BADS (Bayesian Adaptive Direct Search)
options_bads = bads('defaults'); % Load default settings
options_bads.Display = 'none';   % Alternatively, use this to suppress the output
options_bads.Verbosity = 0;

options_fmin = optimoptions('fmincon','MaxIterations', 1e4, 'Display','off');

%% Define names of directories
if isnumeric(isubj)
    subjList =             {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC','SR', 'IO'};
    nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195, 0];
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    nameFolder_NOM = sprintf('%s/ORI%dSF%d/%s/L%d', nameFolder_NOM0, nORI, nSF, subjName, iLocComb);
else % IO
    subjName = isubj{1};
    criterion_true = isubj{2};
    nblocks=0;
    nameFolder_NOM = sprintf('%s/ORI%dSF%d/IO/%s_c%.1f', nameFolder_NOM0, nORI, nSF, subjName, criterion_true);% To Save results
end

if isempty(dir(nameFolder_NOM)), mkdir(nameFolder_NOM), end

nameFileModelIDVD_beforeEst = sprintf('%s/n%d_A%d_beforeEst', nameFolder_NOM, ni, iModelA);
nameFileModelIDVD = sprintf('%s/n%d_A%dB%d', nameFolder_NOM, ni, iModelA, iModelB);

% print
fprintf('\n%s [nblocks = %d] [ORI%d SF%d]\n - ni = %d\n - Loc: %s\n - MODEL: [A%d] %s & [B%d] %s\n\n', ...
    subjName, nblocks, nORI, nSF, ni, namesLocComb{iLocComb}, iModelA, namesModelA{iModelA}, iModelB, namesModelB{iModelB})

%% MAIN LOOP
% prelocate 
nParams = length(namesParamsModel_all{iModelB});
params_est_allB = nan(ni, nParams);
params_est_allB_fmincon = params_est_allB;
nLL_allB = nan(ni,1);
pred_metrics_allB = cell(ni, 1);

fprintf('Running ni = %d: ', ni)

% load data
    load(nameFileModelIDVD_beforeEst, 'data_allB')
    
% Set up parameters upper and lower bound based on iModelB
c_lb = min(data_allB{1}.IV); c_ub = max(data_allB{1}.IV); 
c0 = mean([c_lb, c_ub]); 
SDadd0 = mean([SDadd_lb, SDadd_ub]); 
Nmul0 = mean([Nmul_lb, Nmul_ub]); 

switch iModelB
    case 1, params0 = [lapse0, SDadd0, c0]; params_lb = [lapse_lb, SDadd_lb, c_lb]; params_ub = [lapse_ub, SDadd_ub, c_ub];
    case 2, params0 = [lapse0, Nmul0, c0]; params_lb = [lapse_lb, Nmul_lb, c_lb]; params_ub = [lapse_ub, Nmul_ub, c_ub];
    case 3, params0 = [lapse0, c0]; params_lb=[lapse_lb, c_lb]; params_ub=[lapse_ub, c_ub];
    case 4, params0 = [SDadd0, c0]; params_lb = [SDadd_lb, c_lb]; params_ub = [SDadd_ub, c_ub]; % not fitting lapse rate
    case 5, params0 = [Nmul0, c0]; params_lb = [Nmul_lb, c_lb]; params_ub = [Nmul_ub, c_ub]; % not fitting lapse rate 
end

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

%% SAVE
% copy and rename the xx_beforeEst.mat file
copyfile([nameFileModelIDVD_beforeEst, '.mat'], [nameFileModelIDVD, '.mat']);

save(nameFileModelIDVD, '*_allB', '-append')

time_end = datetime('now')

time_end - time_start

fprintf('\n========== Prediction saved ==========\n\n\n\n\n')
