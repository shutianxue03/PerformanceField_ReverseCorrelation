
function pred = PR_pred_v3(iModelB, params_est, data)
% function metrics_pred = MR_pred(mu_PRS, sigma_PRS, mu_ABS, sigma_ABS, params)

switch iModelB
    case 1, % multiplicative, additive noise, thresh
        alpha = params_est(1);
        std_add = params_est(2);
        thresh = params_est(3);
    case 2, % additive noise, thresh
        alpha = 0;
        std_add = params_est(1);
        thresh = params_est(2);
    case 3, % multiplicative noise, thresh
        std_add = 0;
        alpha = params_est(1);
        thresh = params_est(2);
end

%% extract data
% mu_PRS = truth.mu_PRS;
% mu_ABS = truth.mu_ABS;
% sigma_PRS = truth.sigma_PRS;
% sigma_ABS = truth.sigma_ABS;

mu_PRS = median(data.IV{1});
sigma_PRS = std(data.IV{1});
mu_ABS = median(data.IV{2});
sigma_ABS = std(data.IV{2});
nPRS = length(data.IV{1});
nABS = length(data.IV{2});

%% predict metrics by math
% metrics = [dprime_math, criterion_math, pC_math, pHit_math, pFA_math, pA_math, pA_PRS_math, pA_ABS_math];

metrics = PR_predMetricsMath_v3(mu_PRS, mu_ABS, sigma_PRS, sigma_ABS, params_est, iModelB);

%% predicts responses 
resp = [(data.IV{1} + randn(nPRS, 1) * alpha * sigma_PRS + randn(nPRS, 1) * std_add) > thresh; ...
    (data.IV{2}+ randn(nABS, 1) * alpha * sigma_ABS + randn(nABS, 1) * std_add) > thresh];

%% compile
pred.metrics = metrics;
pred.resp = resp;
