
% estimate the thresh, mu_PRS/ABS and sigma_PRS/ABS from the metrics directly
% thresh = params(1);
% mu_PRS = params(2);
% sigma_PRS = params(3);
% mu_ABS = params(4);
% sigma_ABS  = params(5);

fxn_estMUSIGMA = @(p) fxn_getError(p, metrics_test);
p_0 = [0, 1, 0, 1, 0];
p_lb = [-10, -10, 0, -10, 0];
p_ub = [10, 10, 10, 10, 10];
options = optimoptions('fmincon','MaxIterations', 1e4, 'Display','off');
problem_ML = createOptimProblem('fmincon','objective', fxn_estMUSIGMA,'x0', p_0, 'lb', p_lb, 'ub', p_ub, 'options', options);
ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel', 1);
[p_est, nLL] = run(ms_ML, problem_ML, nrep);
metrics_pred = fxn_getM(p_est);

%%
function nLL = fxn_getError(params, metrics_data)
metrics_pred = fxn_getM(params);
nLL = sum(-log(normpdf(metrics_pred, metrics_data)));
end

%%
function metrics = fxn_getM(params)
thresh = params(1);
mu_PRS = params(2);
sigma_PRS = params(3);
mu_ABS = params(4);
sigma_ABS  = params(5);

pA_buffer=0;

% predict pHit/pFA/pC3
pHit_math = 1-normcdf(thresh, mu_PRS, sigma_PRS);
pFA_math = 1-normcdf(thresh, mu_ABS, sigma_ABS);
pC_math = (pHit_math+1-pFA_math)/2;

% predict dprime/criterion
[dprime_math , criterion_math] = SX_sim06_SDT(pHit_math, pFA_math);

% predict pA by math
PDF_PRS = @(x) normpdf(x, mu_PRS, sigma_PRS);
PDF_ABS = @(x) normpdf(x, mu_ABS, sigma_ABS);
CDF_PRS = normcdf(thresh, mu_PRS, sigma_PRS);
CDF_ABS = normcdf(thresh, mu_ABS, sigma_ABS);

fxn_PRS = @(x) PDF_PRS(x) .* (CDF_PRS.^2 + (1-CDF_PRS).^2);
fxn_ABS = @(x) PDF_ABS(x) .* (CDF_ABS.^2 + (1-CDF_ABS).^2);
pA_PRS_math = integral(fxn_PRS, -inf, inf) + pA_buffer;
pA_ABS_math = integral(fxn_ABS, -inf, inf) + pA_buffer;
pA_math = mean([pA_PRS_math, pA_ABS_math]);

% compile
metrics = [dprime_math, criterion_math, pC_math, pHit_math, pFA_math, pA_math, pA_PRS_math, pA_ABS_math];
end
