
function pred = MR_pred(params_est, data)
% function metrics_pred = MR_pred(mu_PRS, sigma_PRS, mu_ABS, sigma_ABS, params)

nparams_model = length(params_est);
switch nparams_model
    case 1
        alpha_PRS = params_est(1);
        alpha_ABS = alpha_PRS;
        thresh = 0;
    case 2
        alpha_PRS = params_est(1);
        alpha_ABS = alpha_PRS;
        thresh = params_est(2);
    case 3
        alpha_PRS = params_est(1);
        alpha_ABS = params_est(2);
        thresh = params_est(3);
end

%% extract data
% mu_PRS = truth.mu_PRS;
% mu_ABS = truth.mu_ABS;
% sigma_PRS = truth.sigma_PRS;
% sigma_ABS = truth.sigma_ABS;

mu_PRS = median(data.IV(1,:));
sigma_PRS = std(data.IV(1,:));
mu_ABS = median(data.IV(2,:));
sigma_ABS = std(data.IV(2,:));

%% predict metrics by math
pred = MR_predMetricsMath(mu_PRS, mu_ABS, sigma_PRS, sigma_ABS, alpha_PRS, alpha_ABS, thresh);

