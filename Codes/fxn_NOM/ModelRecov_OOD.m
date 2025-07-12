% Model Recovery
% estimate internal-external noise ratio by minimizing the calculated pA and estimated pA


addpath(genpath('Data_OOD'))
addpath(genpath('fxn_model'))
addpath(genpath('fxn_analysis_RC'))

%% shell
alpha_all = [.1, .2, .5, .7, .9, 1, 1.2, 1.5]; % alpha_PRS; alpha_ABS is always prop. to alpha_PRS
% alpha_all = [1, 1.2, 1.5]; % alpha_PRS; alpha_ABS is always prop. to alpha_PRS
nalpha = length(alpha_all);

for ialpha = 1:nalpha
    nB = 100; % number of bootstrapping
    nrep = 10; % 30
    truth.ratio_PRSABS = .5; % alpha_ABS/alpha_PRS
    truth.nparams_model = 3;
    % 1 = alpha (assumed the same for PRS and ABS)
    % 2 = alpha (assumed the same for PRS and ABS) and thresh
    % 3 = alpha_PRS, alpha_ABS and thresh
    errorComp = 'dprime, criterion, pA3, pC3';
    errorFxn = '@(pred, data) sum(-log(normpdf(pred, data, 1)), ''all'')';
    
    %% simulation SETTING
    truth.mu_PRS = .65; % mean of the internal variable distribution [PRS]
    truth.sigma_PRS = .09; % std of the internal variable distribution [PRS], treated as 5he external noise
    truth.mu_ABS = .5;
    truth.sigma_ABS = .075;
    truth.thresh = .6; % applied to the NON-standardized IV distirbution
    ntrials_sim = 1e4; % [simulation] num of trials simulated to generate the data set
    
    %% model fitting SETTING
    options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
    normIV = 1;
    nLoc = 1;
    
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
    
    %% ALPHA
    alpha = alpha_all(ialpha);
    truth.alpha_PRS = alpha;
    truth.alpha_ABS = truth.alpha_PRS*truth.ratio_PRSABS;
    fprintf('\n\nAlpha = %.1f starting... (%d/%d)\n', alpha, ialpha, nalpha)
    
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
    
    %%
    fprintf('Simulation and fitting started')
    for iB = 1:nB
        %%    simulation  %
        switch truth.nparams_model
            case 1, params_true = alpha;
            case 2, params_true = [alpha, truth.thresh];
            case 3, params_true = [truth.alpha_PRS, truth.alpha_ABS, truth.thresh];
        end
        
        if iB==1, plotFlag = 1; else, plotFlag = 0; end
        sim = MR_sim(ntrials_sim, truth, params_true, plotFlag);
        IV_sim_allB(iB, :, :) = sim.IV;    % nB x 2 ntrials_sim
        resp_sim_allB(iB, :, :) = sim.resp;    % nB x 2 ntrials_sim
        metric_sim_allB(iB, :) = sim.metrics; % nB x 8 (dprime, criterion, pC, pHit, pFA, pA, pA_PRS, pA_ABS)
        
        %%  model fitting  %
        fxn_estParams = @(params) fxn_getError(params, sim);
        
        problem_ML = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub, 'options',options);
        ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel',1);
        params_est = run(ms_ML, problem_ML, nrep);
        %         params_est = fmincon(fxn_estParams, params0, [], [], [], [], params_lb, params_ub, [], options);
        params_est_allB(iB, :) = params_est;
        
        %%   Prediction  %
        pred = MR_pred(params_est, sim);
        metric_pred_allB(iB, :) = pred;
        
        if ~mod(iB, nB/10), fprintf('='), end
    end % end of ii
    
    fprintf('Done\n')
    
    % save for each alpha
    save(sprintf('Data_OOD/ModelRecov_alpha%d', round(alpha*10)), 'truth', 'IV_pred_allB', 'resp_pred_allB', 'metric_sim_allB', 'metric_pred_allB', 'params_est_allB')
end
