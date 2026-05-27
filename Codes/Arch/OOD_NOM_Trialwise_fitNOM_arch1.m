function OOD_NOM_Trialwise_fitNOM_arch1(isubj, iLocComb, iModelA, iModelB, nBoot)
%==========================================================================%
% OOD_NOM_Trialwise_fitNOM.m
%--------------------------------------------------------------------------
% Author: Shutian Xue
% Last update: 2025-08-27
%
% Part of the OOD_NOM_Trialwise pipeline.
%
% This function:
% 1) Fits a trial-wise noisy observer model (NOM; Model B) to internal
% variables (IVs) for the TEST set, precomputed in
% OOD_NOM_Trialwise_compIV.m.
% 2) Uses the fitted parameters to predict behavioral metrics (pYES, pC,
% pA) based on binned empirical IVs.
%
% Key model functions:
% - fxn_getError_v7(): computes negative log-likelihood (pYES + pA)
% - PR_pred_v7(): predicts metrics from IVs and fitted parameters
%
% INPUTS
% isubj : subject index (for human data), or cell array {name, c_true}
% for IO simulations.
% iLocComb : location combination index
% iModelA : template / IV model index (defined in SX_RC1_setting)
% iModelB : internal noise / correlation model index (see fxn_getError_v7)
% nBoot : number of resampling bootstraps (same as in *_compIV)
%
% OUTPUTS
% Results are saved to:
% nameFile_fitNOM = '<Data_NOM_Trialwise>/<subj>/L<iLocComb>/n<nBoot>_A<iModelA>B<iModelB>.mat'
%
% Saved variables (appended to *_compIV.mat):
% params_est_allB : [nBoot x nParams] fitted parameters
% nLL_allB : [nBoot x 1] negative log-likelihood
% pred_metrics_allB : {nBoot x 1} predictions from PR_pred_v7
%
%==========================================================================%

clc; close all;
warning off;
format compact;
time_start = datetime('now')
rng(123); % define see for reproducibility

addpath(genpath('fxn_exp'));
addpath(genpath('fxn_NOM'));
addpath(genpath('fxn_RCplot'));
addpath(genpath('fxn_analysis_RC_v2'));
addpath(genpath('SX_toolbox/bads-master'));

%% General settings (from master config)
%--------------%
SX_RC1_setting; % defines nORI, nSF, namesLocComb, namesModelA, namesModelB, nBins, etc.
%--------------%
flag_fminconORbads = 1; % 1 = use fmincon (faster, local); 2 = use BADS (slower, more robust)
flag_plot_allIter = 0; % 1 = make summary plots across bootstraps
flag_plot_perIter = 0; % 1 = plot per-bootstrap fits (can be slow)
if ~strcmp('HPC', str_envir), flag_plot_allIter = 1; end % don't plot when running on HPC

%% Set up file paths and names
if isnumeric(isubj)
    % Human subjects
    subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC', 'SR'};
    nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195];

    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);

    nameFolder_OOD_load = sprintf('%s/%s', nameFolder_Data_OOD, subjName); % behav & energy root
    nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocComb);
else
    % Ideal observer or simulated observer
    subjName = isubj{1};
    criterion_true = isubj{2}; %#ok<NASGU> % kept for compatibility, not used here
    nblocks = 0;

    nameFolder_OOD_load = sprintf('%s/%s', nameFolder_Data_OOD, subjName); % behav & energy root
    nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, subjName);
end

if isempty(dir(nameFolder_NOM_save))
    mkdir(nameFolder_NOM_save);
end

% File names:
nameFile_compIV = sprintf('%s/n%d_A%d_compIV', nameFolder_NOM_save, nBoot, iModelA);
nameFile_fitNOM = sprintf('%s/n%d_A%dB%d', nameFolder_NOM_save, nBoot, iModelA, iModelB);

%% Print header
fprintf('\n===========================\nStep 2: Fit NOM and predict metrics \n=========================== \n\n')
fprintf(['\nSubject/IO name: %s ' ...
    '\n - L%d [%s]', ...
    '\n - A%d [%s]', ...
    '\n - B%d [%s]', ...
    '\n - Number of bootstraps = %d\n\n'], ...
    subjName, ...
    iLocComb, namesLocComb{iLocComb}, ...
    iModelA, namesModelA{iModelA}, ...
    iModelB, namesModelB{iModelB}, ...
    nBoot);

%% Load trial-wise data and criterion from xx_compIV.mat file
load(nameFile_compIV, 'data_allBoot', 'c_zscore');
fprintf('\n Loaded xx_compIV.mat for "data_allBoot" and "c_zscore"\n')

%% Parameter vectors per Model B
switch iModelB
    case 1  % FullModel: Induced + constant independent + constant shared
        params0   = [NOMp1_0,   NOMp2_0,   NOMp3_0];
        params_lb = [NOMp1_lb, NOMp2_lb, NOMp3_lb];
        params_ub = [NOMp1_ub, NOMp2_ub, NOMp3_ub];

    case 2  % NoSharedN: Induced + constant independent (no shared noise)
        params0   = [NOMp1_0,   NOMp2_0];
        params_lb = [NOMp1_lb, NOMp2_lb];
        params_ub = [NOMp1_ub, NOMp2_ub];

    case 3  % NoInducedN: constant independent + constant shared (no induced noise)
        params0   = [NOMp2_0,   NOMp3_0];
        params_lb = [NOMp2_lb, NOMp3_lb];
        params_ub = [NOMp2_ub, NOMp3_ub];

    case 4  % NoIdpdtN: Induced + constant shared (no independent constant noise)
        params0   = [NOMp1_0,   NOMp3_0];
        params_lb = [NOMp1_lb, NOMp3_lb];
        params_ub = [NOMp1_ub, NOMp3_ub];

    case 5  % JustIdpdtN: constant independent only
        params0   = [NOMp2_0];
        params_lb = [NOMp2_lb];
        params_ub = [NOMp2_ub];

    case 6  % JustInducedN: induced only
        params0   = [NOMp1_0];
        params_lb = [NOMp1_lb];
        params_ub = [NOMp1_ub];

    case 7  % JustSharedM: constant shared only
        params0   = [NOMp3_0];
        params_lb = [NOMp3_lb];
        params_ub = [NOMp3_ub];

    otherwise
        error('OOD_NOM_Trialwise_fitNOM: Unknown iModelB = %d', iModelB);
end

% switch iModelB
%     case 1 % Full model: (1) Nmul + (2) SDadd + (3) rho [the best model]
%         params0 = [NOMp1_0, NOMp2_0, NOMp3_0];
%         params_lb = [NOMp1_lb, NOMp2_lb, NOMp3_lb];
%         params_ub = [NOMp1_ub, NOMp2_ub, NOMp3_ub];
% 
%     case 2 % Nmul + SDadd (no rho)
%         params0 = [NOMp1_0, NOMp2_0];
%         params_lb = [NOMp1_lb, NOMp2_lb];
%         params_ub = [NOMp1_ub, NOMp2_ub];
% 
%     case 3 % SDadd + rho (no Nmul)
%         params0 = [NOMp2_0, NOMp3_0];
%         params_lb = [NOMp2_lb, NOMp3_lb];
%         params_ub = [NOMp2_ub, NOMp3_ub];
% 
%     case 4 % Nmul + rho (no SDadd)
%         params0 = [NOMp1_0, NOMp3_0];
%         params_lb = [NOMp1_lb, NOMp3_lb];
%         params_ub = [NOMp1_ub, NOMp3_ub];
% 
%     case 5 % SDadd only
%         params0 = [NOMp2_0];
%         params_lb = [NOMp2_lb];
%         params_ub = [NOMp2_ub];
% 
%     case 6 % Nmul only
%         params0 = [NOMp1_0];
%         params_lb = [NOMp1_lb];
%         params_ub = [NOMp1_ub];
% 
%     case 7 % rho only (no internal noise) [the worst model]
%         params0 = [NOMp3_0];
%         params_lb = [NOMp3_lb];
%         params_ub = [NOMp3_ub];
% 
%     otherwise
%         error('OOD_NOM_Trialwise_fitNOM: Unknown iModelB = %d', iModelB);
% end

%% Preallocate outputs
nParams = length(namesModelBparams{iModelB});
params_est_allBoot = nan(nBoot, nParams);
nLL_allBoot = nan(nBoot, 1);
pred_metrics_allBoot = cell(nBoot, 1);

%% Main estimation loop across bootstraps

fprintf('\nRunning nBoot = %d: ', nBoot);

for iBoot = 1:nBoot

    fprintf('%d... ', iBoot);

    % Extract data
    data = data_allBoot{iBoot};

    % if any(iModelB == [1,3,4,7]) % Fit all params separately: for ModelBs with rho
    %     % Objective function for optimizer: nLL from trial-wise pYES + pairwise pA
    %     iStep = 1;
    %     fxn_estParams = @(paramsNOM) fxn_getError_v7(iModelB, paramsNOM, data, c_zscore, iStep);
    %     switch iModelB
    %         case 1, params0_step1 = params0(1:2); params_lb_step1 = params_lb(1:2); params_ub_step1 = params_ub(1:2);
    %         case 3, params0_step1 = params0(1); params_lb_step1 = params_lb(1); params_ub_step1 = params_ub(1);
    %         case 4, params0_step1 = params0(1); params_lb_step1 = params_lb(1); params_ub_step1 = params_ub(1);
    %         case 7, error('ALERT: we are no longer fitting iModelB=7 (no noise)')
    %     end
    %     % Estimate parameters
    %     if flag_fminconORbads == 1 % Faster, local search
    %         [params_est_step1, nLL_step1] = fmincon(fxn_estParams, params0_step1, [], [], [], [], params_lb_step1, params_ub_step1, [], options_fmin);
    %     else % BADS: more robust global + local search
    %         [params_est_step1, nLL_step1] = bads(fxn_estParams, params0_step1, params_lb_step1, params_ub_step1, [], [], [], options_bads);
    %     end
    % 
    %     iStep = 2;
    %     fxn_estParams = @(paramsNOM) fxn_getError_v7(iModelB, paramsNOM, data, c_zscore, iStep, params_est_step1);
    %     switch iModelB
    %         case 1, params0_step2 = params0(3); params_lb_step2 = params_lb(3); params_ub_step2 = params_ub(3);
    %         case 3, params0_step2 = params0(2); params_lb_step2 = params_lb(2); params_ub_step2 = params_ub(2);
    %         case 4, params0_step2 = params0(2); params_lb_step2 = params_lb(2); params_ub_step2 = params_ub(2);
    %         case 7, error('ALERT: we are no longer fitting iModelB=7 (no noise)')
    %     end
    %     % Estimate parameters
    %     if flag_fminconORbads == 1 % Faster, local search
    %         [params_est_step2, nLL_step2] = fmincon(fxn_estParams, params0_step2, [], [], [], [], params_lb_step2, params_ub_step2, [], options_fmin);
    %     else % BADS: more robust global + local search
    %         [params_est_step2, nLL_step2] = bads(fxn_estParams, params0_step2, params_lb_step2, params_ub_step2, [], [], [], options_bads);
    %     end
    %     params_est = [params_est_step1, params_est_step2];
    %     nLL = nLL_step1+nLL_step2;

    % else % Fit all params together
        iStep = 0;
        % fxn_estParams = @(paramsNOM) fxn_getError_v7(iModelB, paramsNOM, data, c_zscore, iStep); % the full model has induced, constant noise and rho (two passes are correlated)
        fxn_estParams = @(paramsNOM) fxn_getError_v8(iModelB, paramsNOM, data, c_zscore, iStep); % the full model has induced, constant noise (shared and independent across passes)

        % Estimate parameters
        if flag_fminconORbads == 1 % Faster, local search
            [params_est, nLL] = fmincon(fxn_estParams, params0, [], [], [], [], params_lb, params_ub, [], options_fmin);
        else % BADS: more robust global + local search
            [params_est, nLL] = bads(fxn_estParams, params0, params_lb, params_ub, [], [], [], options_bads);
        end
    % end

    % Compile
    params_est_allBoot(iBoot, :) = params_est(:).';
    nLL_allBoot(iBoot) = nLL;

    % Predict binned metrics from estimated parameters
    % pred = PR_pred_v7(iModelB, params_est, data, c_zscore, nBins, flag_plot_perIter);
    pred = PR_pred_v8(iModelB, params_est, data, c_zscore, nBins, flag_plot_perIter);
    pred_metrics_allBoot{iBoot} = pred;

end % end of iBoot

% fprintf('\n\nAll bootstraps DONE\n');

%% Save results (append onto *_compIV.mat)
% copyfile([nameFile_compIV, '.mat'], [nameFile_fitNOM, '.mat']); % backup structure from compIV
save(nameFile_fitNOM, '*_allBoot');
fprintf('\n========== Fitting saved ==========\n\n\n\n\n');

%%
if flag_plot_allIter
    % NOMplot_fitNOM;
end
close all;

%% Timing info
time_end = datetime('now')
elapsed = time_end - time_start;
fprintf('\n\nDONE (time used: %s)\n', char(elapsed));

end
