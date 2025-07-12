
function error = fxn_getError(fitMode, params, data, errorComp)


% make prediction
pred = PR_pred(params, data);

% calculate error
switch fitMode
%     case 1 % Method 1: fit metrics
%         
%         getError = @(pred, data) sum(-log(normpdf(pred, data, 1)), 'all');
%         error = 0;
%         for imetric = errorComp % (1) dprime, (2) criterion, (3) pC, (4) pHit, (5) pFA, (6) pA, (7) pA_PRS, (8) pA_ABS
%             %     error = error + getError(pred(ierror), sim.metrics_math(ierror));
%             error = error + getError(pred.metrics(imetric), data.metrics_sim(imetric));
%         end
        
    case 2 % method 2: fit trial-wise responses (Fernandez et al 2022)
        % pC is from the predicted respones
%         (1) dprime, (2) criterion, (3) pC, (4) pHit, (5) pFA, (6) pA, (7) pA_PRS, (8) pA_ABS

% way 1: get nLL for PRS and ABS separately
        pC_PRS = pred.metrics(4); % pC comes from the prediction give a param comb
        pC_ABS = 1-pred.metrics(5);
        m = length(data.resp(1,:)); % m and n are the data/simulation
        n_PRS = sum(data.resp(1,:));
        n_ABS = m - sum(data.resp(2,:));
        fxn_nLL = @(pC, m, n) -(n*log(pC) + (m-n) * log(1-pC));
        nLL_PRS = fxn_nLL(pC_PRS, m, n_PRS); if nLL_PRS == Inf, nLL_PRS = 1/eps; end
        nLL_ABS = fxn_nLL(pC_ABS, m, n_ABS); if nLL_ABS == Inf, nLL_ABS = 1/eps; end
        error = nLL_PRS + nLL_ABS;

 % way 2: pC refers to the accuracy of all trials
%         pC = pred.metrics(3);
%         % m and n are the data
%         m_ = length(data.resp(1,:)); 
%         m = m_*2; % number of total trials
%         n_PRS = sum(data.resp(1,:));
%         n_ABS = m_ - sum(data.resp(2,:));
%         n = n_PRS + n_ABS; % number of correct trials
%         fxn_nLL = @(pC, m, n) -(n*log(pC) + (m-n) * log(1-pC));
%         nLL = fxn_nLL(pC, m, n); if nLL == Inf, nLL = 1/eps; end
%         error = nLL;
        
        % if isnan(error), fprintf('ERROR: The calculated error is nan!'), end
end