
function pred = PR_pred_v2(iModelB, params_est, data)
% add gamma parameter

N_mul = params_est(1);
sd_add = params_est(2);
gamma = params_est(3);
thresh = params_est(4);

%% extract data
IV_PRS = data.IV{1} .^ gamma;
IV_ABS = data.IV{2}.^ gamma;

mu_PRS = median(IV_PRS);
mu_ABS = median(IV_ABS);

sigma_PRS = std(IV_PRS);
sigma_ABS = std(IV_ABS);

nPRS = length(IV_PRS);
nABS = length(IV_ABS);

%% predicts responses
if iModelB==4
    IV_PRS_noiseAdd = IV_PRS +  (ratio .* IV_PRS);
    IV_ABS_noiseAdd = IV_ABS + (ratio .* IV_ABS);
else
    IV_PRS_noiseAdd = IV_PRS + randn(nPRS, 1) * N_mul * sigma_PRS + randn(nPRS, 1) * sd_add;
    IV_ABS_noiseAdd = IV_ABS + randn(nABS, 1) * N_mul * sigma_ABS + randn(nABS, 1) * sd_add;
end
% quickPlot_IV_noiseAdded

%%
resp = [ IV_PRS_noiseAdd > thresh; IV_ABS_noiseAdd > thresh];

%% predict metrics by math
if iModelB==4
    sigma_IE_PRS = sqrt(1+ratio^2) *sigma_PRS;
    sigma_IE_ABS = sqrt(1+ratio^2) *sigma_ABS;
else
    %     sigma_IE_PRS = std(IV_PRS_noiseAdd);
    %     sigma_IE_ABS = std(IV_ABS_noiseAdd);
    sigma_IE_PRS = sqrt((1+N_mul^2) * sigma_PRS^2 + sd_add^2);
    sigma_IE_ABS = sqrt((1+N_mul^2) * sigma_ABS^2 + sd_add^2);
end

% predict pHit/pFA/pC3
pHit_math = 1-normcdf(thresh, mu_PRS, sigma_IE_PRS);
pFA_math = 1-normcdf(thresh, mu_ABS, sigma_IE_ABS);
pC_math = (pHit_math+1-pFA_math)/2;

%% predict dprime/criterion
[dprime_math , criterion_math] = SX_sim06_SDT(pHit_math, pFA_math);

% predict pA by math
PDF_PRS = @(x) normpdf(x, mu_PRS, sigma_IE_PRS);
PDF_ABS = @(x) normpdf(x, mu_ABS, sigma_IE_ABS);
CDF_PRS = normcdf(thresh, mu_PRS, sigma_IE_PRS);
CDF_ABS = normcdf(thresh, mu_ABS, sigma_IE_ABS);

fxn_PRS = @(x) PDF_PRS(x) .* (CDF_PRS.^2 + (1-CDF_PRS).^2);
fxn_ABS = @(x) PDF_ABS(x) .* (CDF_ABS.^2 + (1-CDF_ABS).^2);

pA_PRS_math = integral(fxn_PRS, -inf, inf);% + pA_buffer;
pA_ABS_math = integral(fxn_ABS, -inf, inf);%+ pA_buffer;
pA_math = pA_PRS_math * nPRS/(nPRS+nABS) + pA_ABS_math * nABS/(nPRS+nABS);

% compile
metrics = [dprime_math, criterion_math, pC_math, pHit_math, pFA_math, pA_math, pA_PRS_math, pA_ABS_math];

%% compile
pred.metrics = metrics;
pred.resp = resp;


