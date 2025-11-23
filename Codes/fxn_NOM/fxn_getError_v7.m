function nLL = fxn_getError_v7(iModelB, params_est, data, c_zscore)
% Compute negative log-likelihood for the trial-wise noisy observer model,

% Extract data
resp_allT = data.resp(:); % 1 = YES, 0 = NO

% Predict pYES and pA based on params
[pYES_pred_allT, pA_pred_allPairs, consistency_allPairs] = fxn_predMetrics(iModelB, params_est, data, c_zscore);

% Calculate nLL based on pYES 
nLL_pYES = -sum( resp_allT .* log(pYES_pred_allT) + (1 - resp_allT) .* log(1 - pYES_pred_allT) );

% Calculate nLL based on pA
nLL_pA = -sum( consistency_allPairs .* log(pA_pred_allPairs) + (1 - consistency_allPairs) .* log(1 - pA_pred_allPairs) );

% Combine
nLL = nLL_pYES + nLL_pA;

