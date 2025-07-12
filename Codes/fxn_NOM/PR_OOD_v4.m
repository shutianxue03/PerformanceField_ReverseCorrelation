% Model Recovery
% estimate internal-external noise ratio by minimizing the nLL calculated
% based on trial-wise responses
% version 4: three parameters: (1) N_mul, (2) additive noise (3) thresh

addpath(genpath('Data_OOD'))
addpath(genpath('fxn_model'))
addpath(genpath('fxn_analysis_RC'))

clc
close all
clear all

% when fitting is weirdly bad,  just CLEAR ALL!!!

% %% simulation SETTING
truth.mu_PRS = 2.3;
truth.sigma_PRS = 0.44;
truth.mu_ABS = 1.4;
truth.sigma_ABS = 0.29;

truth.mu_PRS = .24; % mean of the internal variable distribution [PRS]
truth.mu_ABS = .19;
truth.sigma_PRS = 0.046; % std of the internal variable distribution [PRS], treated as the external noise
truth.sigma_ABS = 0.045;
flagIncludePA=1;
flag_pAbuffer = 1;
ntrials_sim = 1e3; % [simulation] num of trials simulated to generate the data set
iModelB = 3; % (1) N_mul, sd_add, thresh (2) sd_add, thresh (3) N_mul, thresh (4) ratio, thresh
params_true = [0, 1.2];
% params_true = [.5, 1.85];
% define the range
params0 = [1, 1.85]; 
params_lb = [1e-5, 1.4];
params_ub = [2, 2.3]; 
flagConstrainThresh = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ni = 1e3; % number of bootstrapping
nrep = 20; % number of repetitions to search for the global min
truth.ndata_PRS = ntrials_sim;

if flagIncludePA
    namesModelB = {'Add. + Multi. + pA buffer + thresh', 'Add. + pA buffer + thresh', 'Multi. + pA buffer + thresh'};
    namesParamsModel_all = {{'Multi noise', 'Additive noise', 'pA buffer', 'Threshold'}, ...
        {'Additive noise', 'pA buffer', 'Threshold'}, {'Multi noise','pA buffer', 'Threshold'}};
    namesParams2 = {'Thresh', 'mu [PRS]', 'sigma [PRS]','mu [ABS]', 'sigma [ABS]', 'pA buffer'};
else
    namesModelB = {'Add. + Multi. + thresh', 'Add. + thresh.', 'Multi. + thresh'};
    namesParamsModel_all = {{'Multi noise', 'Additive noise', 'Threshold'}, ...
        {'Additive noise', 'Threshold'}, {'Multi noise', 'Threshold'}};
    namesParams2 = {'Thresh', 'mu [PRS]', 'sigma [PRS]','mu [ABS]', 'sigma [ABS]'};
    
end

namesParamsModel = namesParamsModel_all{iModelB}; assert(length(namesParamsModel) == length(params_true))
truth.nparams_model = length(namesParamsModel);
fitMode = 2; % (1) fit 8 metrics (2) fit trial-wise responses
nCountdowns = 10; % number of countdowns

%% model fitting SETTING
namesMetrics_short = {'dprime',  'Criterion', 'pC', 'pHit', 'pFA', 'pA', 'pA [PRS]', 'pA [ABS]'};
nmetrics = length(namesMetrics_short);
options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
options_GS = optimoptions(@fmincon,'Algorithm','sqp');

%% empty folders for all iterations
% simulation
metric_sim_allB = nan(ni, 8);
nLL_true_allB = nan(ni,1);
thresh_est_allB = nan(ni, 1);
% estimated params
params_est_allB = nan(ni, truth.nparams_model);
nLL_allB = nan(ni, 1);
% predictions
metric_pred_allB = nan(ni, 8);

%%
fprintf('Running...')
for ii = 1:ni
    %    SIMULATION  %
    sim = PR_sim(iModelB, params_true, truth);
    data.IV{1} = sim.IV(:, 1);
    data.IV{2} = sim.IV(:, 2);
    data.resp{1} = sim.resp(:, 1);
    data.resp{2} = sim.resp(:, 2);
    data.m_pA = sim.m_pA;
    data.n_pA = sim.n_pA;
    data.metrics_sim = sim.metrics_sim;
    metric_sim_allB(ii, :) = sim.metrics_sim; % nB x 8 (dprime, criterion, pC, pHit, pFA, pA, pA_PRS, pA_ABS)
    nLL_true_allB(ii) = fxn_getError_PR(iModelB, params_true, data); % get the true nLL of the current data set
    
    % get IV median & std
    mu_PRS = median(data.IV{1});
    mu_ABS = median(data.IV{2});
    sigma_PRS = std(data.IV{1});
    sigma_ABS = std(data.IV{2});

    % estimate the thresh given the measured criterion, 
    % to constrain the range of thresh when estimating params
    data.noisy_mu_PRS = median(sim.IV_noisy(:, 1));
    data.noisy_mu_ABS = median(sim.IV_noisy(:, 2));
    data.noisy_sigma_PRS = std(sim.IV_noisy(:, 1));
    data.noisy_sigma_ABS = std(sim.IV_noisy(:, 2));
    fxn_estThresh = @(thresh_est) fxn_getErrorThresh(thresh_est, data, sim.metrics_sim(2));
    [thresh_est, ~] = fmincon(fxn_estThresh,(mu_PRS + mu_ABS)/2, [], [], [], [], min([mu_PRS, mu_ABS]), max([mu_PRS,  mu_ABS]), [], options);
    thresh_est_allB(ii) = thresh_est;
    buffer = .1;
    %     buffer = min([sigma_PRS, sigma_ABS]);
    if flagConstrainThresh
    params0(end) = thresh_est;
    params_lb(end) = thresh_est - buffer;
    params_ub(end) = thresh_est + buffer;
    end
    %  MODEL FITTING %
    fxn_estParams = @(params) fxn_getError_PR(iModelB, flagIncludePA, params, data);
    problem_MS = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub, 'options',options);
    ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel',1);
    [params_est, nLL] = run(ms_ML, problem_MS, nrep);
    params_est_allB(ii, :) = params_est;
    nLL_allB(ii) = nLL;
    
    %   PREDICTION  %
    pred = PR_pred(iModelB, flagIncludePA, params_est, data);
    metric_pred_allB(ii, :) = pred.metrics;
    
    if ~mod(ii, ni/nCountdowns), fprintf('%d ', nCountdowns-round(ii/ni*nCountdowns)+1), end
end % end of ii

fprintf('\nDONE')

%%  plot
% Fig 1. plot params estimated
PRplot_params_est
saveas(gcf, sprintf('publishedPDFs/model/PR_M%d_Fig1_B1e%d_n1e%d.jpg', iModelB, log10(ni), log10(ntrials_sim)))

% Fig 2. pred vs. data
PRplot_pred
saveas(gcf, sprintf('publishedPDFs/model/PR_M%d_Fig2_B1e%d_n1e%d.jpg', iModelB, log10(ni), log10(ntrials_sim)))

% Fig 3. corr among params estimated
PRplot_paramsCorr
saveas(gcf, sprintf('publishedPDFs/model/PR_M%d_Fig3_B1e%d_n1e%d.jpg', iModelB, log10(ni), log10(ntrials_sim)))

%% To estimate the thresh given a criterion
function c_pred = predCriterion(thresh, data)
pHit = 1-normcdf(thresh, data.noisy_mu_PRS, data.noisy_sigma_PRS);
pFA = 1-normcdf(thresh, data.noisy_mu_ABS, data.noisy_sigma_ABS);
[~, c_pred] = SX_sim06_SDT(pHit, pFA);
end

function nLL = fxn_getErrorThresh(thresh, data, c_measured)
c_pred = predCriterion(thresh, data);
nLL = -log(normpdf(c_pred, c_measured));
end

