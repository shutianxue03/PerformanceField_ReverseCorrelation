
function pred = PR_pred_v4(iModelB, params_est, data)
% function metrics_pred = MR_pred(mu_PRS, sigma_PRS, mu_ABS, sigma_ABS, params)

switch iModelB
    case 1 % multiplicative, additive noise, thresh
        N_mul = params_est(1);
        sd_add = params_est(2);
        thresh = params_est(3);
    case 2 % additive noise, thresh
        N_mul = 0;
        sd_add = params_est(1);
        thresh = params_est(2);
    case 3 % multiplicative noise, thresh
        sd_add = 0;
        N_mul = params_est(1);
        thresh = params_est(2);
    case 4
        N_mul=0;
        sd_add=0;
        thresh = params_est(1);
end

%% extract data
IV_PRS = data.IV{1};% .^ gamma;
IV_ABS = data.IV{2};%.^ gamma;

mu_PRS = median(IV_PRS);
mu_ABS = median(IV_ABS);

sigma_PRS = std(IV_PRS);
sigma_ABS = std(IV_ABS);

nPRS = length(IV_PRS);
nABS = length(IV_ABS);

%% calculate responses
IV_noisy_PRS = IV_PRS * (1+N_mul) + randn(nPRS, 1) * sd_add;
IV_noisy_ABS = IV_ABS * (1+N_mul) + randn(nABS, 1) * sd_add;
resp = [IV_noisy_PRS > thresh; IV_noisy_ABS > thresh];

%% predict metrics by math
mu_PRS_pred = mu_PRS*(1+N_mul);
mu_ABS_pred = mu_ABS*(1+N_mul);
sigma_PRS_pred = sqrt((1+N_mul)^2 * sigma_PRS^2 + sd_add^2); % the same
sigma_ABS_pred = sqrt((1+N_mul)^2 * sigma_ABS^2 + sd_add^2); % the same
% predict pHit/pFA/pC3
pHit_pred = 1-normcdf(thresh, mu_PRS_pred, sigma_PRS_pred);
pFA_pred = 1-normcdf(thresh, mu_ABS_pred, sigma_ABS_pred);
% pC_pred = (pHit_pred+1-pFA_pred)/2;
pC_pred = pHit_pred * nPRS/(nPRS+nABS) + (1-pFA_pred) * nABS/(nPRS+nABS);

%% predict dprime/criterion
[dprime_pred , criterion_pred] = SX_sim06_SDT(pHit_pred, pFA_pred);

% predict pA by math
PDF_PRS = @(x) normpdf(x, mu_PRS_pred, sigma_PRS_pred);
PDF_ABS = @(x) normpdf(x, mu_ABS_pred, sigma_ABS_pred);
CDF_PRS = normcdf(thresh, mu_PRS_pred, sigma_PRS_pred);
CDF_ABS = normcdf(thresh, mu_ABS_pred, sigma_ABS_pred);

fxn_PRS = @(x) PDF_PRS(x) .* (CDF_PRS.^2 + (1-CDF_PRS).^2);
fxn_ABS = @(x) PDF_ABS(x) .* (CDF_ABS.^2 + (1-CDF_ABS).^2);

pA_PRS_pred = integral(fxn_PRS, -inf, inf);
pA_ABS_pred = integral(fxn_ABS, -inf, inf);
pA_pred = pA_PRS_pred * nPRS/(nPRS+nABS) + pA_ABS_pred * nABS/(nPRS+nABS);

% compile
metrics = [dprime_pred, criterion_pred, pC_pred, pHit_pred, pFA_pred, pA_pred, pA_PRS_pred, pA_ABS_pred];

%% compile
pred.metrics = metrics;
pred.resp = resp;

