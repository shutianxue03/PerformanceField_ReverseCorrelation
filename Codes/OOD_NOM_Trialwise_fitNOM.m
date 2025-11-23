function OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, iModelA, iModelB, nIterations)
%==========================================================================%
% OOD_NOM_Trialwise_fitNOM.m
%--------------------------------------------------------------------------
% Author:      Shutian Xue
% Last update: 2025-08-27
%
% Part of the OOD_NOM_Trialwise pipeline.
%
% This function:
%   1) Fits a trial-wise noisy observer model (NOM; Model B) to internal
%      variables (IVs) for the TEST set, precomputed in
%      OOD_NOM_Trialwise_compIV.m.
%   2) Uses the fitted parameters to predict behavioral metrics (pYES, pC,
%      pA) based on binned empirical IVs.
%
% Key model functions:
%   - fxn_getError_v7():   computes negative log-likelihood (pYES + pA)
%   - PR_pred_v7():        predicts metrics from IVs and fitted parameters
%
% INPUTS
%   isubj        : subject index (for human data), or cell array {name, c_true}
%                  for IO simulations.
%   iLocComb     : location combination index
%                    1 = Fovea
%                    8 = Perifovea (6° ecc; pooled)
%                    6 = HM (Left + Right horizontal meridian)
%                    7 = VM (Upper + Lower vertical meridian)
%                    5 = LVM
%                    3 = UVM
%   iModelA      : template / IV model index (defined in SX_RC1_setting)
%                    1 = core model
%                    2 = permuted template
%                    3 = IO template
%   iModelB      : internal noise / correlation model index (see fxn_getError_v7)
%                    1 = Nmul + SDadd + rho
%                    2 = Nmul + SDadd (no rho)
%                    3 = rho only (no internal noise)
%                    4 = SDadd + rho (no Nmul)
%                    5 = SDadd only
%                    6 = Nmul + rho (no SDadd)
%                    7 = Nmul only
%   nIterations  : number of resampling iterations (same as in *_compIV)
%
% OUTPUTS
%   Results are saved to:
%       nameFile_fitNOM = '<Data_NOM_Trialwise>/<subj>/L<iLocComb>/n<nIterations>_A<iModelA>B<iModelB>.mat'
%
%   Saved variables (appended to *_compIV.mat):
%       params_est_allB   : [nIterations x nParams] fitted parameters
%       nLL_allB          : [nIterations x 1] negative log-likelihood
%       pred_metrics_allB : {nIterations x 1} predictions from PR_pred_v7
%
%==========================================================================%

close all;
warning off;
format compact;
time_start = datetime('now');
rng(123);   % define see for reproducibility

addpath(genpath('fxn_NOM'));
addpath(genpath('fxn_RCplot'));
addpath(genpath('fxn_analysis_RC_v2'));
addpath(genpath('SX_toolbox/bads-master'));

%% General settings (from master config)
%--------------%
SX_RC1_setting;  % defines nORI, nSF, namesLocComb, namesModelA, namesModelB, nBins, etc.
%--------------%
flag_fminconORbads = 2; % 1 = use fmincon (faster, local); 2 = use BADS    (slower, more robust)
flag_plot_allIter = 1;   % 1 = make summary plots across iterations
flag_plot_perIter = 0;   % 1 = plot per-iteration fits (can be slow)

%% Set up file paths and names
if isnumeric(isubj)
    % Human subjects
    subjList =        {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC', 'SR'};
    nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195, 0];

    subjName = subjList{isubj};
    nblocks  = nblocks_allSubj(isubj);

    nameFolder_OOD_load = sprintf('%s/%s', nameFolder_Data_OOD, subjName);  % behav & energy root
    nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocComb);
else
    % Ideal observer or simulated observer
    subjName       = isubj{1};
    criterion_true = isubj{2}; %#ok<NASGU> % kept for compatibility, not used here
    nblocks        = 0;

    nameFolder_OOD_load = sprintf('%s/%s', nameFolder_Data_OOD, subjName);  % behav & energy root
    nameFolder_NOM_save = sprintf('%s/%s',     nameFolder_Data_NOM_Trialwise, subjName);
end

if isempty(dir(nameFolder_NOM_save))
    mkdir(nameFolder_NOM_save);
end

% File names:
nameFile_compIV = sprintf('%s/n%d_A%d_compIV', nameFolder_NOM_save, nIterations, iModelA);
nameFile_fitNOM       = sprintf('%s/n%d_A%dB%d',        nameFolder_NOM_save, nIterations, iModelA, iModelB);

%% Print header
fprintf(['\nSubject/IO name: %s ' ...
    '\n - L%d [%s]', ...
    '\n - A%d [%s]', ...
    '\n - B%d [%s]', ...
    '\n - Number of iterations = %d'], ...
    subjName, ...
    iLocComb, namesLocComb{iLocComb}, ...
    iModelA, namesModelA{iModelA}, ...
    iModelB, namesModelB{iModelB}, ...
    nIterations);

%% Load trial-wise data and criterion
load(nameFile_compIV, 'data_allB', 'c_zscore');

%% Parameter vectors per Model B
switch iModelB
    case 1  % Nmul + SDadd + rho
        params0   = [Nmul0,   SDadd0,   rho0];
        params_lb = [Nmul_lb, SDadd_lb, rho_lb];
        params_ub = [Nmul_ub, SDadd_ub, rho_ub];

    case 2  % Nmul + SDadd (no rho)
        params0   = [Nmul0,   SDadd0];
        params_lb = [Nmul_lb, SDadd_lb];
        params_ub = [Nmul_ub, SDadd_ub];

    case 3  % rho only (no internal noise)
        params0   = [rho0];
        params_lb = [rho_lb];
        params_ub = [rho_ub];

    case 4  % SDadd + rho (no Nmul)
        params0   = [SDadd0,   rho0];
        params_lb = [SDadd_lb, rho_lb];
        params_ub = [SDadd_ub, rho_ub];

    case 5  % SDadd only
        params0   = [SDadd0];
        params_lb = [SDadd_lb];
        params_ub = [SDadd_ub];

    case 6  % Nmul + rho (no SDadd)
        params0   = [Nmul0,   rho0];
        params_lb = [Nmul_lb, rho_lb];
        params_ub = [Nmul_ub, rho_ub];

    case 7  % Nmul only
        params0   = [Nmul0];
        params_lb = [Nmul_lb];
        params_ub = [Nmul_ub];

    otherwise
        error('OOD_NOM_Trialwise_fitNOM: Unknown iModelB = %d', iModelB);
end

%% Preallocate outputs

nParams            = length(namesParamsModel_all{iModelB});
params_est_allB    = nan(nIterations, nParams);
nLL_allB           = nan(nIterations, 1);
pred_metrics_allB  = cell(nIterations, 1);

%% Main estimation loop across iterations

% fprintf('\n\nRunning ni = %d: ', nIterations);

for ii = 1:nIterations
    fprintf('%d ', ii);

    data = data_allB{ii};

    % Objective function for optimizer: nLL from trial-wise pYES + pairwise pA
    fxn_estParams = @(paramsNOM) fxn_getError_v7(iModelB, paramsNOM, data, c_zscore);

    % Estimate parameters
    if flag_fminconORbads == 1 % Faster, local search
        [params_est, nLL] = fmincon(fxn_estParams, params0, [], [], [], [], params_lb, params_ub, [], options_fmin);
    else % BADS: more robust global + local search
        [params_est, nLL] = bads(fxn_estParams, params0, params_lb, params_ub, [], [], [], options_bads);
    end

    % Compile
    params_est_allB(ii, :) = params_est(:).';
    nLL_allB(ii)           = nLL;

    % Predict binned metrics from estimated parameters
    pred = PR_pred_v7(iModelB, params_est, data, c_zscore, nBins, flag_plot_perIter);
    pred_metrics_allB{ii} = pred;

    fprintf('%d ', ii);
end % end of ii

% fprintf('\n\nAll iterations DONE\n');

%% Plot summary across iterations (optional)

if flag_plot_allIter
    NOMplot_fitNOM;
end
close all;

%% Save results (append onto *_compIV.mat)

copyfile([nameFile_compIV, '.mat'], [nameFile_fitNOM, '.mat']);  % backup structure from compIV
save(nameFile_fitNOM, '*_allB', '-append');

%% Timing info
time_end = datetime('now');
elapsed = time_end - time_start;
fprintf('\n\nDONE (time used: %s)\n', char(elapsed));

end
