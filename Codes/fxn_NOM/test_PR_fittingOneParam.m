
% fit one of all parameter while fixing the others at one time
% check whether the nLL are the same
% check for tradeoffs between parameters

clc
close all

%% simulation SETTING
truth.mu_PRS = 0.64; % mean of the internal variable distribution [PRS]
truth.sigma_PRS = 0.09; % std of the internal variable distribution [PRS], treated as the external noise
truth.mu_ABS = 0.53;
truth.sigma_ABS = 0.08;
truth.alpha_PRS = 1.2; % alpha_PRS; alpha_ABS is always prop. to alpha_PRS
truth.alpha_ABS = 1.5; % alpha_PRS; alpha_ABS is always prop. to alpha_PRS
truth.thresh =  0.57; % applied to the NON-standardized IV distirbution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ni = 100; % number of bootstrapping
nrep = 10; % number of repetitions to search for the global min
ntrials_sim = 1e4; % [simulation] num of trials simulated to generate the data set
truth.nparams_model = 3;
errorComp = 1:8;
fitMode = 2; % (1) fit 8 metrics (2) fit trial-wise responses
nCountdowns = 10; % number of countdowns

%% model fitting SETTING
titles = {'PRS alpha', 'ABS alpha', 'thresh'};
namesMetrics_short = {'dprime',  'Criterion', 'pC', 'pHit', 'pFA', 'pA', 'pA [PRS]', 'pA [ABS]'};
nmetrics = length(namesMetrics_short);
options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
options_GS = optimoptions(@fmincon,'Algorithm','sqp');

%% 
fileName = 'fxn_model/fittingOneParam.mat';
fileDir = dir(fileName);

if isempty(fileDir)
    params_est_allB = nan(truth.nparams_model, ni);
    nLL_allB = nan(truth.nparams_model, ni);
    for onlyFit = 1
        switch onlyFit
            case 1
                params0 = [.5, .5, truth.thresh];
                params_lb = [.001, .001, truth.thresh];
                params_ub = [5, 5, truth.thresh];
            case 2
                params0 = [truth.alpha_PRS, .5, truth.thresh];
                params_lb = [truth.alpha_PRS, .001, truth.thresh];
                params_ub = [truth.alpha_PRS, 5, truth.thresh];
            case 3
                params0 = [truth.alpha_PRS, truth.alpha_ABS, .5];
                params_lb = [truth.alpha_PRS, truth.alpha_ABS, -1];
                params_ub = [truth.alpha_PRS, truth.alpha_ABS, 2];
        end
        
        for ii = 1:ni
            %    simulation  %
            sim = PR_sim(ntrials_sim, truth, [truth.alpha_PRS, truth.alpha_ABS, truth.thresh], 0);
            
            %  model fitting  %
            fxn_estParams = @(params) fxn_getError(fitMode, params, sim, truth, errorComp);
            
            % since the param search showed that local min = global min, just use fmincon
            %         [params_est, nLL] = fmincon(fxn_estParams, params0, [], [], [], [], params_lb, params_ub, [], options);
            
            % problem for the MultiStart program
            problem_MS = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub, 'options',options);
            ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel',1);
            [params_est_MS, nLL_MS] = run(ms_ML, problem_MS, nrep);
            %     % problem for the Global Search program
            %             problem_GS = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub, 'options', options_GS);
            %             gs = GlobalSearch('Display', 'final');
            %             [params_est_GS, nLL_GS]  = run(gs, problem_GS);
            params_est_allB(onlyFit, ii) = params_est_MS(onlyFit);
            nLL_allB(onlyFit, ii) = nLL_MS;
            
            if ~mod(ii, ni/nCountdowns), fprintf('%d ',nCountdowns-round(ii/ni*nCountdowns)+1), end
        end % end of ii
        fprintf('\n')
    end
    save(fileName, 'params_est_allB', 'nLL_allB')
    fprintf('DONE\n')
else
    load(fileName)
    fprintf('LOADED\n')
end

%%
figure('Position', [2000 200 800 200])
for ip = 1
    subplot(1,4,ip)
    histogram(params_est_allB(ip, :))
    switch ip, case 1, true = truth.alpha_PRS; case 2, true = truth.alpha_ABS; case 3, true = truth.thresh; end
    xline(median(params_est_allB(ip, :)), 'k--', 'linewidth', 2);
    xline(true, 'r-', 'linewidth', 2);
    title(sprintf('sd = %.4f', std(params_est_allB(ip, :))))
    xlabel(titles{ip})
end
legend({'estimation','median', 'true value'}, 'Location', 'south')

subplot(1,4,4), hold on
legends = cell(1,3);
for ip = 1:3
    histogram(nLL_allB(ip, :))
    legends{ip} = sprintf('ONLY fitting %s', titles{ip});
end
legend(legends, 'Location', 'south')
