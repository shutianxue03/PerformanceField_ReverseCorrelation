
mu_PRS = median(IV_PRS);
mu_ABS = median(IV_ABS);
sigma_PRS = std(IV_PRS);
sigma_ABS = std(IV_ABS);

fxn_estThresh = @(thresh_est) fxn_getErrorThresh(thresh_est, data, metrics_test(2));
[thresh_est, ~] = fmincon(fxn_estThresh,(mu_PRS + mu_ABS)/2, [], [], [], [], min([mu_PRS, mu_ABS]), max([mu_PRS,  mu_ABS]), [], options);

% 'end': assuming that 'threshold' is the last free model parameter
if flagConstrainThresh
    buffer = min([sigma_PRS, sigma_ABS]);
    params0(end) = thresh_est; % starting point will be the estimated criterion given the distrbution of raw IVs
    params_lb(end) = thresh_est - buffer; % so that the threshold will not be left to the peak of ABS or right to the peak of PRS
    params_ub(end) = thresh_est + buffer;
else
    params0(end) = mean([mu_PRS, mu_ABS]);
    params_lb(end) = min([mu_PRS, mu_ABS]);
    params_ub(end) = max([mu_PRS, mu_ABS]);
end

%% To estimate the thresh given a criterion
function c_pred = predCriterion(thresh, data)
pHit = 1-normcdf(thresh, data.noisy_mu_PRS, data.noisy_sigma_PRS);
pFA = 1-normcdf(thresh, data.noisy_mu_ABS, data.noisy_sigma_ABS);
[~, c_pred] = SX_sim06_SDT(pHit, pFA);
end

%%
function nLL = fxn_getErrorThresh(thresh, data, c_measured)
c_pred = predCriterion(thresh, data);
nLL = -log(normpdf(c_pred, c_measured));
end
