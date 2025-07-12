% Parameter recovery
% with lapse rate
% For each model variation, simulate data given the true parameters, and fit the model to the data, to
% examine if the true params can be recovered
% also examine whether there is correlation between estimates

addpath(genpath('fxn_NOM'))
addpath(genpath('bads-master'))
clc, close all

% settings
options = bads('defaults');
options.Display = 'off';

%% Data
load('NOM_data')
mu_IV_PRS = median(data.IV{1});
SD_IV_PRS = std(data.IV{1});
mu_IV_ABS = median(data.IV{2});
SD_IV_ABS = std(data.IV{2});

% figure
% hold on
% histogram(data.IV{1}, 'facecolor', 'r')
% histogram(data.IV{2}, 'facecolor', 'b')
% legend({'Signal-present', 'Signal-absent'})
% xline(mu_IV_PRS, 'r', 'linewidth', 2);
% xline(mu_IV_ABS, 'b', 'linewidth', 2);
% xlabel('Internal variable')

%% define true parameter
clc
true_lapse = .01;
true_Nmul = 1;
true_SDadd = .2;
true_crit = 1;

nSamples = 100;
if nSamples>1, flag_plot=0; else, flag_plot=1; end

lb_lapse=eps; ub_lapse=.05;
lb_Nmul=0; ub_Nmul=2;
lb_SDadd=0; ub_SDadd=1;
lb_crit=0; ub_crit=2;

nData_PRS = 5e3;%length(IV_PRS);
nData_ABS = 5e3;%length(IV_ABS);

namesModelsB = {'Full', 'NoNmul', 'NoAdd', 'NoIN'};
iModelB = input('     >>> iModelB: ');
fprintf('          Chosen model: %s\n\n', namesModelsB{iModelB})

switch iModelB
    case 1, params_true = [true_lapse, true_Nmul, true_SDadd, true_crit];
        params_lb = [lb_lapse, lb_Nmul, lb_SDadd, lb_crit];
        params_ub = [ub_lapse, ub_Nmul, ub_SDadd, ub_crit];
        namesParams = {'lapse rate', 'Nmul', 'SDadd', 'Criterion'};
    case 2, params_true = [true_lapse, true_SDadd, true_crit];
        params_lb = [lb_lapse, lb_SDadd, lb_crit];
        params_ub = [ub_lapse, ub_SDadd, ub_crit];
        namesParams = {'lapse rate', 'SDadd', 'Criterion'};
    case 3, params_true = [true_lapse, true_Nmul, true_crit];
        params_lb = [lb_lapse, lb_Nmul, lb_crit];
        params_ub = [ub_lapse, ub_Nmul, ub_crit];
        namesParams = {'lapse rate', 'Nmul', 'Criterion'};
    case 4, params_true = [true_lapse, true_crit];
        params_lb = [lb_lapse, lb_crit];
        params_ub = [ub_lapse, ub_crit];
        namesParams = {'lapse rate', 'Criterion'};
end
params0 = params_true;
nParams = length(params_true);

params_est_allSets = nan(nSamples, nParams);
nLL_allSets = nan(nSamples, 1);

%% simulate one dataset
for iSet = 1:nSamples
    params_est = params_true;
    switch iModelB
        case 1 % multiplicative, additive noise, thresh
            lambda = params_est(1);
            Nmul = params_est(2);
            SDadd = params_est(3);
            crit = params_est(4);
        case 2 % additive noise, thresh
            Nmul = 0;
            lambda = params_est(1);
            SDadd = params_est(2);
            crit = params_est(3);
        case 3 % multiplicative noise, thresh
            SDadd = 0;
            lambda = params_est(1);
            Nmul = params_est(2);
            crit = params_est(3);
        case 4
            Nmul=0;
            SDadd=0;
            lambda = params_est(1);
            crit = params_est(2);
    end
    
    % simulate IVs
    IV_sim_PRS = randn(nData_PRS, 1) * SD_IV_PRS + mu_IV_PRS;
    IV_sim_ABS = randn(nData_ABS, 1) * SD_IV_ABS + mu_IV_ABS;
    % simulate noisy IVs
    IV_noisy_sim_PRS = IV_sim_PRS + randn(nData_PRS, 1).*(IV_sim_PRS*Nmul) + randn(nData_PRS, 1)*SDadd;
    IV_noisy_sim_ABS = IV_sim_ABS + randn(nData_ABS, 1).*(IV_sim_ABS*Nmul) + randn(nData_ABS, 1)*SDadd;
    
    % Simulate base responses (without lapse)
    resp_sim_PRS_base = IV_noisy_sim_PRS > crit;
    resp_sim_ABS_base = IV_noisy_sim_ABS > crit;
    
    % Generate random responses for lapse trials
    random_resp_PRS = rand(nData_PRS, 1) > 0.5;  % random 0 or 1 (50% probability)
    random_resp_ABS = rand(nData_ABS, 1) > 0.5;  % random 0 or 1 (50% probability)
    
    % Apply lapse rate (introduce random responses based on lapse rate)
    is_lapse_PRS = rand(nData_PRS, 1) < lambda;  % Identify trials with lapses
    is_lapse_ABS = rand(nData_ABS, 1) < lambda;  % Identify trials with lapses
    
    % Final responses with lapse incorporated
    resp_sim_PRS = resp_sim_PRS_base;  % Start with the base responses
    resp_sim_PRS(is_lapse_PRS) = random_resp_PRS(is_lapse_PRS);  % Replace with random responses for lapse trials
    resp_sim_ABS = resp_sim_ABS_base;  % Start with the base responses
    resp_sim_ABS(is_lapse_ABS) = random_resp_ABS(is_lapse_ABS);  % Replace with random responses for lapse trials
    
    % organize simulated data
    data_sim.IV_PRS = IV_sim_PRS;
    data_sim.IV_ABS = IV_sim_ABS;
    data_sim.IV_noisy_PRS = IV_noisy_sim_PRS;
    data_sim.IV_noisy_ABS = IV_noisy_sim_ABS;
    data_sim.resp_PRS = resp_sim_PRS;
    data_sim.resp_ABS = resp_sim_ABS;
    
    % plot
    if flag_plot
        figure, hold on
        histogram(IV_sim_PRS, 'FaceColor', 'r', 'Normalization', 'probability')
        histogram(IV_sim_ABS, 'FaceColor', 'b', 'Normalization', 'probability')
        histogram(IV_noisy_sim_PRS, 'FaceColor', 'm', 'Normalization', 'probability')
        histogram(IV_noisy_sim_ABS, 'FaceColor', 'c', 'Normalization', 'probability')
        xlabel('Internal variable space')
        xline(crit, 'k-', 'linewidth', 2);
        legend({'IV [PRS]', 'IV [ABS]',  'Noisy IV [PRS]', 'Noisy IV [ABS]'})
    end
    
    %% fit data
    warning off
    fxn_estParams = @(params) fxn_getError_PR(iModelB, params, data_sim);
    % problem_ML = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub);
    % ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel', 1);
    % [params_est, nLL] = run(ms_ML, problem_ML, nrep);
    %     [params_est, nLL] = fmincon(fxn_estParams, params0, [],[],[],[], params_lb, params_ub);
    [params_est, nLL] = bads(fxn_estParams, params0, params_lb, params_ub, params_lb, params_ub, [], options);
    params_est_allSets(iSet, :) = params_est;
    nLL_allSets(iSet, :) = nLL;
    
    fprintf('\nSet #%d/%d', iSet, nSamples)
    
    %% plot
    if flag_plot
        switch iModelB
            case 1 % multiplicative, additive noise, thresh
                lambda = params_est(1);
                Nmul = params_est(2);
                SDadd = params_est(3);
                crit = params_est(4);
            case 2 % additive noise, thresh
                Nmul = 0;
                lambda = params_est(1);
                SDadd = params_est(2);
                crit = params_est(3);
            case 3 % multiplicative noise, thresh
                SDadd = 0;
                lambda = params_est(1);
                Nmul = params_est(2);
                crit = params_est(3);
            case 4
                Nmul=0;
                SDadd=0;
                lambda = params_est(1);
                crit = params_est(2);
        end
        
        % simulate noisy IVs based on estimates
        IV_noisy_est_PRS = IV_sim_PRS + randn(nData_PRS, 1).*(IV_sim_PRS*Nmul) + randn(nData_PRS, 1)*SDadd;
        IV_noisy_est_ABS = IV_sim_ABS + randn(nData_ABS, 1).*(IV_sim_ABS*Nmul) + randn(nData_ABS, 1)*SDadd;
        % simulate responses based on estimates
        resp_est_PRS = IV_noisy_est_PRS > crit;
        resp_est_ABS = IV_noisy_est_ABS > crit;
        
        % plot
        if flag_plot
            figure, hold on
            % histogram(IV_sim_PRS, 'FaceColor', 'r')
            histogram(IV_noisy_sim_PRS, 'FaceColor', 'm', 'EdgeColor', 'w')
            histogram(IV_noisy_est_PRS, 'FaceColor', 'w', 'EdgeColor', 'm')
            
            % histogram(IV_sim_ABS, 'FaceColor', 'b')
            histogram(IV_noisy_sim_ABS, 'FaceColor', 'c', 'EdgeColor', 'w')
            histogram(IV_noisy_est_ABS, 'FaceColor', 'w', 'EdgeColor', 'c')
            
            xline(crit, 'k--', 'linewidth', 2);
            xline(true_crit, 'k-', 'linewidth', 2);
            legend show
            % legend({'IV [PRS]', 'IV [ABS]',  'Noisy IV [PRS]', 'Noisy IV [ABS]'})
        end
    end
end % iSet

fprintf('\n ALL DONE\n')

%% Figure 1: estimated vs. true parameters
figure('Position', [0 200 1e3 200])

[param_est_med, ~, ~, param_est_SEM]= getCI(params_est_allSets, 1, 1);

for iParam = 1:nParams
    subplot(1, nParams, iParam), hold on
    % idvd estimate
    plot(1+randn(1,nSamples)/10, params_est_allSets(:, iParam), 'k.')
    % true param
    plot(1, params_true(iParam), 'ro')
    % median and error bar of the estimate
    plot(1, param_est_med(iParam), '.', 'HandleVisibility', 'off')
    errorbar(1, param_est_med(iParam), param_est_SEM(iParam), 'r', 'CapSize', 0)
    
    xticks(1), xticklabels(namesParams{iParam}), xlim([0, 2])
    ylim([params_lb(iParam), params_ub(iParam)])
    
    if iParam==nParams
        legend({'Estimate per data set', 'True value', 'Med and 1SEM'})
    end
end

sgtitle(sprintf('True vs. Estimated Parameters\nModel %d: %s (nSamples=%d)', iModelB, namesModelsB{iModelB}, nSamples))
set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
% set(findall(gcf, '-property', 'fontsize'), 'fontsize',15)

saveas(gcf, sprintf('Fig/M%d-%s_TRUEvsEST.jpg', iModelB, namesModelsB{iModelB}))

%% Figure 2: check correlation between estimates
% just check corr between Nmul and SDadd for now, i..e, only for iModelB=1
% later, I should do this for every pair of parameters for all models

indParams_corr = getAllUniquePairs(1:nParams);
nCorr = size(indParams_corr, 1);

figure('Position', [0 200 1e3 400])
for iCorr=1:nCorr
    iParam1 = indParams_corr(iCorr, 1);
    iParam2 = indParams_corr(iCorr, 2);
    
    subplot(2,3,iCorr), hold on, box on
    plot(params_est_allSets(:, iParam1), params_est_allSets(:, iParam2), 'k.')
    plot(params_true(iParam1), params_true(iParam2), 'ro', 'MarkerSize', 10)
    xlim([params_lb(iParam1), params_ub(iParam1)])
    ylim([params_lb(iParam2), params_ub(iParam2)])
    xlabel(namesParams{iParam1})
    ylabel(namesParams{iParam2})
    [r, p] = corr(params_est_allSets(:, iParam1), params_est_allSets(:, iParam2));
    title(sprintf('r=%.2f, p=%.3f', r, p))
%     title(sprintf('%s vs. %s: r=%.2f, p=%.3f', namesParams{iParam1}, namesParams{iParam2}, r, p))

    if iCorr==1
        legend({'Estimate per data set', 'True value'})
    end
    
end
sgtitle(sprintf('Correlation between estimates\nModel %d: %s (nSamples=%d)', iModelB, namesModelsB{iModelB}, nSamples))
set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
set(findall(gcf, '-property', 'fontsize'), 'fontsize',12)
saveas(gcf, sprintf('Fig/M%d-%s_Corr.jpg', iModelB, namesModelsB{iModelB}))

