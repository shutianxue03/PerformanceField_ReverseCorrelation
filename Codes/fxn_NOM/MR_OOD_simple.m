% Model Recovery
% estimate internal-external noise ratio by minimizing the calculated pA and estimated pA


addpath(genpath('Data_OOD'))
addpath(genpath('fxn_model'))
addpath(genpath('fxn_analysis_RC'))

clc
close all

%% simulation SETTING
truth.mu_PRS = 0.64; % mean of the internal variable distribution [PRS]
truth.sigma_PRS = 0.09; % std of the internal variable distribution [PRS], treated as 5he external noise
truth.mu_ABS = 0.53;
truth.sigma_ABS = 0.08;
truth.alpha_PRS = 2; % alpha_PRS; alpha_ABS is always prop. to alpha_PRS
truth.alpha_ABS = 1.5; % alpha_PRS; alpha_ABS is always prop. to alpha_PRS
truth.thresh =  0.55; % applied to the NON-standardized IV distirbution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nB = 100; % number of bootstrapping
nrep = 10; % 30
ntrials_sim = 1e6; % [simulation] num of trials simulated to generate the data set
truth.nparams_model = 3;

%% model fitting SETTING
namesMetrics_short = {'dprime',  'Criterion', 'pC', 'pHit', 'pFA', 'pA', 'pA [PRS]', 'pA [ABS]'};
nmetrics = length(namesMetrics_short);
options = optimoptions('fmincon','MaxIterations',5000,'Display','off');

% define the range
switch truth.nparams_model
    case 1
        params0 = .5;
        params_lb = 0;
        params_ub = 2;
    case 2
        params0 = [.5, 0];
        params_lb = [0, -1];
        params_ub = [2, 1];
    case 3
        params0 = [.5, .5, 0];
        params_lb = [0, 0, -1];
        params_ub = [3, 3, 2];
end

assert(truth.nparams_model == length(params_ub));

%% empty folders for all iterations
% simulation
IV_sim_allB = nan(nB, 2, ntrials_sim);
resp_sim_allB = IV_sim_allB;
metric_sim_allB = nan(nB, 8);

% estimated params
params_est_allB = nan(nB, truth.nparams_model);

% predictions
IV_pred_allB = nan(nB, 2, ntrials_sim);
resp_pred_allB = IV_sim_allB;
metric_pred_allB = nan(nB, 8);

%% params_true
switch truth.nparams_model
    case 1, params_true = alpha;
    case 2, params_true = [alpha, truth.thresh];
    case 3, params_true = [truth.alpha_PRS, truth.alpha_ABS, truth.thresh];
end

%%
for iB = 1:nB
    %    simulation  %   
    sim = MR_sim(ntrials_sim, truth, params_true, 0);
    IV_sim_allB(iB, :, :) = sim.IV;    % nB x 2 ntrials_sim
    resp_sim_allB(iB, :, :) = sim.resp;    % nB x 2 ntrials_sim
    metric_sim_allB(iB, :) = sim.metrics_math; % nB x 8 (dprime, criterion, pC, pHit, pFA, pA, pA_PRS, pA_ABS)
    
    %  model fitting  %
    fxn_estParams = @(params) fxn_getError(params, sim);
    
    problem_ML = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub, 'options',options);
    ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel',1);
    params_est = run(ms_ML, problem_ML, nrep);
    params_est_allB(iB, :) = params_est;
    
    %   Prediction  %
    pred = MR_pred(params_est, sim);
    metric_pred_allB(iB, :) = pred;
    
    if ~mod(iB, nB/10), fprintf('%d ',10-round(iB/nB*10)+1), end
end % end of ii

%% plot
figure('Position', [2000 500 1000 200])

%%%%%%%%%%%%%%
% plot data %
%%%%%%%%%%%%%%
subplot(1,4,1), hold on
histogram(squeeze(IV_sim_allB(1,1,:)), 'DisplayStyle', 'stair', 'edgecolor', 'r', 'Normalization', 'Probability')
histogram(squeeze(IV_sim_allB(1,2,:)), 'DisplayStyle', 'stair', 'edgecolor', 'b', 'Normalization', 'Probability')
xline(params_true(3), 'linewidth', 2);
ylabel('Probability')
xlabel('Internal variable')
legend({'IV [PRS]','IV [ABS]','True thresh'},'Location','south')

%%%%%%%%%%%%%%%%%%%%%%%%%
% plot estimated params %
%%%%%%%%%%%%%%%%%%%%%%%%%
titles_paramsEST = {'alpha [PRS]', 'alpha [ABS]', 'thresh'};
% if truth.nparams_model == 3, np= 2; else, np = 1;end
for ip = 1:truth.nparams_model
    subplot(1,4, ip+1), hold on
    
    % true value
    bar(1, params_true(ip), 'FaceColor','w','EdgeColor','k');
    
    % estimation
    params_est = params_est_allB(:, ip);
    params_est_med = median(params_est);
    params_est_CI68 = quantile(params_est, [.16, .84]);
    errorbar(1, params_est_med, params_est_med-params_est_CI68(1), params_est_CI68(2)-params_est_med, 'ok-')
    
    xticks([]),xticklabels([])
    if ip ~= 3, ylim([0, 3]), else, ylim([truth.thresh-.1, truth.thresh+.1]), end
    if ip == 3, legend({'Truth', 'Estimation'}), end
    
    title(sprintf('Estimated %s', titles_paramsEST{ip}))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot predicted metrics %
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',[2000 300 1000 150])
titles_long = {'Sensitivity (d'')',  'Criterion', 'Proportion correct', 'Hit rate', 'False alarm rate', 'Prop. agreement', 'Prop. agreement [PRS]', 'Prop. agreement [ABS]'};
ylims = [0, 2; -1,1; 0,1;0,1;0,1;.5, 1;.5,1;.5,1];
for im = 1:8
    MR_plotPRED(im, namesMetrics_short{im} , metric_sim_allB(:,im), metric_pred_allB(:,im), ylims(im, :))
    xticks([]),xticklabels([])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot corr between estimated params
labels = {'Alpha [PRS]', 'Alpha [ABS]', 'Thresh'};
cc_all = [1,2;2,3;1,3];
figure('Position',[2000 0 800 200])
for ic = 1:3
    ind = cc_all(ic, :);
    x = params_est_allB(:,ind(1));
    y = params_est_allB(:,ind(2));
    subplot(1,3,ic), hold on
    plot(x, y,'o')
    plot(params_true(ind(1)),params_true(ind(2)), 'r*')
    xlabel(labels{ind(1)})
    ylabel(labels{ind(2)})
    [r,p] = corr(x, y);
    title(sprintf('r=%.3f, p=%.3f',r, p))
    
end
sgtitle(sprintf('Correlation between estimated params (nB = %d)', nB))

%% plot correlation between the simulated and predicted metrics
for im = 1:nmetrics
    subplot(1,nmetrics, im)
    plot(metric_sim_allB(:, im), metric_pred_allB(:, 8), 'o')
end

%%
function MR_plotPRED(iplot, title_, data, pred, ylim_)

subplot(1,8,iplot), hold on

% get ave and CI
data_med = median(data);
pred_med = median(pred);
pred_CI68 = quantile(pred, [.16, .84]) ;

% plot
bar(1, data_med)
errorbar(1, pred_med, pred_med - pred_CI68(1), pred_CI68(2)-pred_med, 'ok')
ylim(ylim_)
title(title_)
end

