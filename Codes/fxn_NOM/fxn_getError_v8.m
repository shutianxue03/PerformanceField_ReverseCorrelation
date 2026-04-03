function nLL = fxn_getError_v8(iModelB, params_est, data, Cz_emp, C_contribution, iStep, params_fromStep1)
% Compute negative log-likelihood for the trial-wise noisy observer model,

if nargin == 6, params_fromStep1 = []; end

% Extract data
resp_allT = data.resp(:); % 1 = YES, 0 = NO

%% Compute nLL based on iStep (0=combine pYES and pA; 1=pYES only (step 1); 2=pA only (step 2))
switch iStep
    case 0
        % Predict pYES and pA based on params
        %---------------------------%
        [pYES_pred_allT, pA_pred_allPairs, consistency_allPairs, Cz_fit] = fxn_predMetrics_v3(iModelB, params_est, data);
        %---------------------------%

        % Calculate nLL based on pYES
        nLL_pYES = -sum( resp_allT .* log(pYES_pred_allT) + (1 - resp_allT) .* log(1 - pYES_pred_allT) );

        % Calculate nLL based on pA
        nLL_pA = -sum( consistency_allPairs .* log(pA_pred_allPairs) + (1 - consistency_allPairs) .* log(1 - pA_pred_allPairs) );

        % Combine
        nLL = nLL_pYES + nLL_pA;

        % Add criterion penalty (actually does not make any difference...)
        % penaltyCrit = (Cz_fit - Cz_emp)^2;
        % sigmaC = 1;
        % nLL_C = 0.5*log(2*pi*sigmaC^2) + 0.5*((Cz_emp - Cz_fit)^2)/(sigmaC^2);
        % wCrit = C_contribution*nLL_comb / ((1-C_contribution)*nLL_C); % calculate the weight of criterion loss
        % nLL = nLL_comb + wCrit*nLL_C;

    case 1
        % Step 1: Predict pYES and estimate only the internal noise parameter
        switch iModelB
            case 1, params_est = [params_est([1,2]), 0, params_est(3)]; 
            case 3, params_est = [params_est(1), 0, params_est(2)]; 
            case 4, params_est = [params_est(1), 0, params_est(2)]; 
        end
        %---------------------------%
        [pYES_pred_allT] = fxn_predMetrics_v3(iModelB, params_est, data);
        %---------------------------%
        % Calculate nLL based on pYES
        nLL = -sum( resp_allT .* log(pYES_pred_allT) + (1 - resp_allT) .* log(1 - pYES_pred_allT) );
        
    case 2
        % Step 2: Predict pA and estimate only Nshared
        switch iModelB
            case 1, params_est = [params_fromStep1([1,2]), params_est, params_fromStep1(3)];
            case 3, params_est = [params_fromStep1(1), params_est, params_fromStep1(2)]; 
            case 4, params_est = [params_fromStep1(1), params_est, params_fromStep1(2)]; 
        end
        %---------------------------%
        [~, pA_pred_allPairs, consistency_allPairs] = fxn_predMetrics_v3(iModelB, params_est, data);
        %---------------------------%

        % Calculate nLL based on pA
        nLL = -sum( consistency_allPairs .* log(pA_pred_allPairs) + (1 - consistency_allPairs) .* log(1 - pA_pred_allPairs) );
end
