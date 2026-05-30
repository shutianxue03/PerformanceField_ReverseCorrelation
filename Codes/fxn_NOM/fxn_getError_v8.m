function [nLL, nLL_pYES, nLL_pA] = fxn_getError_v8(iModelB, params_est, data, C_contribution, iStep, params_fromStep1)
% Compute negative log-likelihood for the trial-wise noisy observer model,

if nargin < 6, params_fromStep1 = []; end

% Optional diagnostic terms (returned when requested by caller).
nLL_pYES = NaN;
nLL_pA = NaN;

% Extract data
resp_allT = data.resp(:); % 1 = YES, 0 = NO

% Get Cz_emp
pHit_emp = mean(data.resp(data.iPRS == 1));
pFA_emp  = mean(data.resp(data.iPRS == 0));
[~, Cz_emp] = SX_sim06_SDT(pHit_emp, pFA_emp);
data.Cz_emp = Cz_emp;

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

        % Combine behavioral likelihood terms
        nLL = nLL_pYES + nLL_pA;

        % % Add criterion penalty in z units, if available
        % if isfield(data, 'Cz_emp') && ~isempty(data.Cz_emp) && isfinite(data.Cz_emp)
        %     sigmaC = 0.2; % smaller sigmaC → stricter penalty; larger sigmaC → weaker penalty;e similar to sd
        %     % Calculate nLL, assuming Cz_emp is normally distibuted, centered at Cz_fit and has SD=sigmaC
        %     nLL_C = 0.5 * ((Cz_fit - data.Cz_emp) / sigmaC)^2;
        % 
        %     % convex combination-like weighting
        %     nTrials = numel(resp_allT);
        %     nLL_behav = (nLL_pYES + nLL_pA) / nTrials;
        %     nLL = (1 - C_contribution) * nLL_behav + C_contribution * nLL_C;
        % end

    case 1
        % Step 1: fit private-noise params (Nshared forced to 0)
        switch iModelB
            case 1, params_est = [params_est([1,2]), 0, params_est(3)]; % input=[Nmul,Nadd,criterion]  → full=[Nmul,Nadd,Nshared=0,criterion]
            case 2, params_est = [params_est(1), 0, params_est(2)];    % input=[Nadd,criterion] → model2=[Nadd,Nshared=0,criterion]
            case 3, params_est = [params_est(1), 0, params_est(2)];    % input=[Nmul,criterion] → model3=[Nmul,Nshared=0,criterion]
        end
        %---------------------------%
        [pYES_pred_allT, ~, ~, Cz_fit] = fxn_predMetrics_v3(iModelB, params_est, data);
        %---------------------------%
        % Calculate nLL based on pYES
        nLL_pYES = -sum( resp_allT .* log(pYES_pred_allT) + (1 - resp_allT) .* log(1 - pYES_pred_allT) );
        % Combine
        nLL = nLL_pYES;

        % Add criterion penalty (same as case 0), so criterion is constrained
        % toward the empirical value even during step 1 of two-step fitting
        if isfield(data, 'Cz_emp') && ~isempty(data.Cz_emp) && isfinite(data.Cz_emp)
            sigmaC = 0.5;
            nLL_C = 0.5 * ((Cz_fit - data.Cz_emp) / sigmaC)^2;
            nTrials = numel(resp_allT);
            nLL_behav = nLL_pYES / nTrials;
            nLL = (1 - C_contribution) * nLL_behav + C_contribution * nLL_C;
        end
        
    case 2
        % Step 2: Predict pA and estimate only Nshared
        switch iModelB
            case 1, params_est = [params_fromStep1([1,2]), params_est, params_fromStep1(3)];
            case 2, params_est = [params_fromStep1(1), params_est, params_fromStep1(2)]; 
            case 3, params_est = [params_fromStep1(1), params_est, params_fromStep1(2)]; 
        end
        %---------------------------%
        [~, pA_pred_allPairs, consistency_allPairs] = fxn_predMetrics_v3(iModelB, params_est, data);
        %---------------------------%

        % Calculate nLL based on pA
        nLL_pA = -sum( consistency_allPairs .* log(pA_pred_allPairs) + (1 - consistency_allPairs) .* log(1 - pA_pred_allPairs) );
        nLL = nLL_pA;
end
