function nLL = fxn_getError_v8(iModelB, params_est, data, c_zscore, iStep, params_fromStep1)
% Compute negative log-likelihood for the trial-wise noisy observer model,

if nargin == 5, params_fromStep1 = []; end
% Extract data
resp_allT = data.resp(:); % 1 = YES, 0 = NO

%% Compute nLL based on iStep (0=combine pYES and pA; 1=pYES only (step 1); 2=pA only (step 2))
switch iStep
    case 0
        % Predict pYES and pA based on params
        [pYES_pred_allT, pA_pred_allPairs, consistency_allPairs] = fxn_predMetrics_v2(iModelB, params_est, data, c_zscore);

        % Calculate nLL based on pYES
        nLL_pYES = -sum( resp_allT .* log(pYES_pred_allT) + (1 - resp_allT) .* log(1 - pYES_pred_allT) );

        % Calculate nLL based on pA
        nLL_pA = -sum( consistency_allPairs .* log(pA_pred_allPairs) + (1 - consistency_allPairs) .* log(1 - pA_pred_allPairs) );

        % Combine
        nLL = nLL_pYES + nLL_pA;

    case 1
        % Step 1: Predict pYES and estimate only the internal noise parameter
        % append rho=0
        if find(iModelB==[1,3,4,7]), params_est = [params_est, 0]; end
        [pYES_pred_allT] = fxn_predMetrics(iModelB, params_est, data, c_zscore);
        % Calculate nLL based on pYES
        nLL = -sum( resp_allT .* log(pYES_pred_allT) + (1 - resp_allT) .* log(1 - pYES_pred_allT) );
        
    case 2
        % Step 2: Predict pA and estimate only rho
        params_est = [params_fromStep1, params_est];
        [~, pA_pred_allPairs, consistency_allPairs] = fxn_predMetrics(iModelB, params_est, data, c_zscore);

        % Calculate nLL based on pA
        nLL = -sum( consistency_allPairs .* log(pA_pred_allPairs) + (1 - consistency_allPairs) .* log(1 - pA_pred_allPairs) );
end

