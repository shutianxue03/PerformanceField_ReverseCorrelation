% Parameter recovery
% with lapse rate
% For each model variation, simulate data given the true parameters, and fit the model to the data, to
% examine if the true params can be recovered
% also examine whether there is correlation between estimates

addpath(genpath('fxn_NOM'))
clc, close all

%% Data
load('NOM_data')
mu_IV_PRS = median(data.IV{1});
SD_IV_PRS = std(data.IV{1});
mu_IV_ABS = median(data.IV{2});
SD_IV_ABS = std(data.IV{2});

figure
hold on
histogram(data.IV{1}, 'facecolor', 'r')
histogram(data.IV{2}, 'facecolor', 'b')
legend({'Signal-present', 'Signal-absent'})
xline(mu_IV_PRS, 'r', 'linewidth', 2);
xline(mu_IV_ABS, 'b', 'linewidth', 2);
xlabel('Internal variable')

%% define true parameter
% true_lapse = 0;
true_Nmul = 1;
true_SDadd = .5;
true_crit = 1;

% lb_lapse=0; ub_lapse=.05;
lb_Nmul=0; ub_Nmul=.2;
lb_SDadd=0; ub_SDadd=1;
lb_crit=0; ub_crit=2;

nData_PRS = 5e3;%length(IV_PRS);
nData_ABS = 5e3;%length(IV_ABS);

iModelB = 1; % TEMPORARY!!

switch iModelB
    case 1, params_true = [true_Nmul, true_SDadd, true_crit];
        params_lb = [lb_Nmul, lb_SDadd, lb_crit];
        params_ub = [ub_Nmul, ub_SDadd, ub_crit];
    case 2, params_true = [true_SDadd, true_crit];
        params_lb = [lb_SDadd, lb_crit];
        params_ub = [ub_SDadd, ub_crit];
    case 3, params_true = [true_Nmul, true_crit];
        params_lb = [lb_Nmul, lb_crit];
        params_ub = [ub_Nmul, ub_crit];
    case 4, params_true = [true_crit];
        params_lb = [lb_crit];
        params_ub = [ub_crit];
end
params0 = params_true;


%% simulate data
% using raw data
params_est = params_true;

switch iModelB
    case 1 % multiplicative, additive noise, thresh
%         lambda = params_est(1);
        Nmul = params_est(1);
        SDadd = params_est(2);
        crit = params_est(3);
    case 2 % additive noise, thresh
        Nmul = 0;
%         lambda = params_est(1);
        SDadd = params_est(1);
        crit = params_est(2);
    case 3 % multiplicative noise, thresh
        SDadd = 0;
%         lambda = params_est(1);
        Nmul = params_est(1);
        crit = params_est(2);
    case 4
        Nmul=0;
        SDadd=0;
%         lambda = params_est(1);
        crit = params_est(1);
end

% simulate IVs
IV_sim_PRS = randn(nData_PRS, 1) * SD_IV_PRS + mu_IV_PRS;
IV_sim_ABS = randn(nData_ABS, 1) * SD_IV_ABS + mu_IV_ABS;
% simulate noisy IVs
IV_noisy_sim_PRS = IV_sim_PRS + randn(nData_PRS, 1).*(IV_sim_PRS*Nmul) + randn(nData_PRS, 1)*SDadd;
IV_noisy_sim_ABS = IV_sim_ABS + randn(nData_ABS, 1).*(IV_sim_ABS*Nmul) + randn(nData_ABS, 1)*SDadd;

% Simulate base responses (without lapse)
resp_sim_PRS = IV_noisy_sim_PRS > crit;
resp_sim_ABS = IV_noisy_sim_ABS > crit;

% organize simulated data
data_sim.IV_PRS = IV_sim_PRS;
data_sim.IV_ABS = IV_sim_ABS;
data_sim.IV_noisy_PRS = IV_noisy_sim_PRS;
data_sim.IV_noisy_ABS = IV_noisy_sim_ABS;
data_sim.resp_PRS = resp_sim_PRS;
data_sim.resp_ABS = resp_sim_ABS;

% plot
figure, hold on
histogram(IV_sim_PRS, 'FaceColor', 'r')
histogram(IV_sim_ABS, 'FaceColor', 'b')
histogram(IV_noisy_sim_PRS, 'FaceColor', 'm')
histogram(IV_noisy_sim_ABS, 'FaceColor', 'c')
xline(crit, 'k-', 'linewidth', 2);
legend({'IV [PRS]', 'IV [ABS]',  'Noisy IV [PRS]', 'Noisy IV [ABS]'})

%% fit data
clc
warning off
fxn_estParams = @(params) fxn_getError_PR_noLapse(iModelB, params, data_sim);
% problem_ML = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub);
% ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel', 1);
% [params_est, nLL] = run(ms_ML, problem_ML, nrep);
[params_est, nLL] = fmincon(fxn_estParams, params0, [],[],[],[], params_lb, params_ub);

params_est
params_true

%% plot
switch iModelB
    case 1 % multiplicative, additive noise, thresh
%         lambda = params_est(1);
        Nmul = params_est(1);
        SDadd = params_est(2);
        crit = params_est(3);
    case 2 % additive noise, thresh
        Nmul = 0;
%         lambda = params_est(1);
        SDadd = params_est(1);
        crit = params_est(2);
    case 3 % multiplicative noise, thresh
        SDadd = 0;
%         lambda = params_est(1);
        Nmul = params_est(1);
        crit = params_est(2);
    case 4
        Nmul=0;
        SDadd=0;
%         lambda = params_est(1);
        crit = params_est(1);
end

% simulate noisy IVs based on estimates
IV_noisy_est_PRS = IV_sim_PRS + randn(nData_PRS, 1).*(IV_sim_PRS*Nmul) + randn(nData_PRS, 1)*SDadd;
IV_noisy_est_ABS = IV_sim_ABS + randn(nData_ABS, 1).*(IV_sim_ABS*Nmul) + randn(nData_ABS, 1)*SDadd;
% simulate responses based on estimates
resp_est_PRS = IV_noisy_est_PRS > crit;
resp_est_ABS = IV_noisy_est_ABS > crit;

% plot
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

