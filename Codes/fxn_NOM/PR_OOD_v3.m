% Model Recovery
% estimate internal-external noise ratio by minimizing the nLL calculated
% based on trial-wise responses% the original version: free parameters are alpha_PRS, alpha_ABS and thresh
% v1-5 are diff attenpts to get rid of the tradeoff between parameters by testing variations of the model

% version 3: two params: [alpha, thresh]

addpath(genpath('Data_OOD'))
addpath(genpath('fxn_model'))
addpath(genpath('fxn_analysis_RC'))

clc
close all
clear 

%% simulation SETTING
truth.mu_PRS = 0.64; % mean of the internal variable distribution [PRS]
truth.sigma_PRS = 0.09; % std of the internal variable distribution [PRS], treated as the external noise
truth.mu_ABS = 0.53;
truth.sigma_ABS = 0.08;
params_true = [1.3, 0.6];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ni = 100; % number of bootstrapping
nrep = 20; % number of repetitions to search for the global min
ntrials_sim = 1e5; % [simulation] num of trials simulated to generate the data set
truth.nparams_model = 2;
errorComp = 1:8;
fitMode = 2; % (1) fit 8 metrics (2) fit trial-wise responses
nCountdowns = 10; % number of countdowns

%% model fitting SETTING
namesMetrics_short = {'dprime',  'Criterion', 'pC', 'pHit', 'pFA', 'pA', 'pA [PRS]', 'pA [ABS]'};
nmetrics = length(namesMetrics_short);
options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
options_GS = optimoptions(@fmincon,'Algorithm','sqp');

%% define the range
params0 = params_true;
diff = [2,2];
params_lb = params0 - diff;
params_ub = params0 + diff;
% params_lb = [0, 0, 0];
% params_ub = [2, 2, 2];

%% empty folders for all iterations
% simulation
IV_sim_allB = nan(ni, 2, ntrials_sim);
IV_noisy_sim_allB = IV_sim_allB;
resp_sim_allB = IV_sim_allB;
metric_sim_allB = nan(ni, 8);

% estimated params
params_est_allB = nan(ni, truth.nparams_model);
nLL_allB = nan(ni, 1);
% predictions
IV_pred_allB = nan(ni, 2, ntrials_sim);
resp_pred_allB = IV_sim_allB;
metric_pred_allB = nan(ni, 8);

%%
for ii = 1:ni
    %    simulation  %
    sim = PR_sim_v3(ntrials_sim, truth, params_true, 0);
    IV_sim_allB(ii, :, :) = sim.IV;    % nB x 2 ntrials_sim
    IV_noisy_sim_allB(ii, :, :) = sim.IV_noisy;    % nB x 2 ntrials_sim
    resp_sim_allB(ii, :, :) = sim.resp;    % nB x 2 ntrials_sim
    metric_sim_allB(ii, :) = sim.metrics_sim; % nB x 8 (dprime, criterion, pC, pHit, pFA, pA, pA_PRS, pA_ABS)
    
    %  model fitting  %
    fxn_estParams = @(params) fxn_getError_v3(params, sim);
    
    % since the param search showed that local min = global min, just use fmincon
    [params_est, nLL] = fmincon(fxn_estParams, params0, [], [], [], [], params_lb, params_ub, [], options);
    
    % problem for the MultiStart program
%     problem_MS = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub, 'options',options);
%     ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel',1);
%     [params_est, nLL] = run(ms_ML, problem_MS, nrep);
    %     % problem for the Global Search program
    %     problem_GS = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub, 'options', options_GS);
    %     gs = GlobalSearch('Display', 'final');
    %     [params_est_GS, nLL_GS]  = run(gs, problem_GS);
    
    params_est_allB(ii, :) = params_est;
    nLL_allB(ii) = nLL;
    
    %   Prediction  %
    pred = PR_pred_v3(params_est, sim);
    metric_pred_allB(ii, :) = pred.metrics;
    
    if ~mod(ii, ni/nCountdowns), fprintf('%d ',nCountdowns-round(ii/ni*nCountdowns)+1), end
end % end of ii

%%
sim = PR_sim_v3(ntrials_sim, truth, params_true, 0); % get sim otherwise the sim inside the parfor loop is not accessible
nLL_true = fxn_getError_v3(params_true, sim, truth);

%% plot params estimated
% plot the estimated alpha_ABS = ratio * alpha_PRS
figure('Position', [2000 0 1000 200])
titles = {'alpha [PRS]', 'thresh'};
PRplot_params_est_v3(IV_sim_allB, params_true, truth, params_est_allB,  nLL_allB, nLL_true, titles)

%% plot predictions
PRplot_pred(sim, metric_sim_allB, metric_pred_allB, namesMetrics_short)
set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)

%% plot corr between params estimated
figure('Position', [2000 0 800 200])
PRplot_paramsCorr_v3(params_est_allB, params_true, titles, ni)
% PRplot_paramsCorr(params_est_allB, params_true, titles, ni, params_lb, params_ub)
