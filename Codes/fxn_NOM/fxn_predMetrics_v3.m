function [pYES_pred_allT, pA_pred_allPairs, consistency_allPairs, Cz_fit] = fxn_predMetrics_v3(iModelB, params_est, data)
% fxn_predMetrics_v3
% Shared core that computes:
% - sigma_pred_allT (trial-wise DV SD)
% - pYES_pred_allT (trial-wise P(YES))
% - pA_pred_allPairs (pair-wise P(agreement))
% - consistency_allPairs (empirical agreement per pair)
%
% INPUTS
% data.IV : nTrials x 1 internal variable
% data.resp : nTrials x 1 responses (1 = YES, 0 = NO)
% data.iPair : nTrials x 1 pair ID for double-pass
% iModelB : index of noise model variant
% params_est : fitted parameters for the chosen model
% c_zscore : criterion in z units (scaled by sigma_pred and added to median(IV))
%
% New parameterization:
% Nmul : multiplicative (Multi) noise coefficient
% SDidpdt : Additive noise (across passes)
% SDshared : Shared noise (across passes)
% criterion
%
% DV correlation across passes is implied by SDshared, not fit as rho.

%% -------- 1. Parse model parameters (noise, lapse) --------
lambda = 0; % Set lapse rate at 0

% Default values (in case a component is absent in a reduced model)
M = 0;
sigmaAdd = 0;
sigmaShared = 0;
criterion = 0;

switch iModelB
    case 1 % FullModel: Multi + Additive + Shared
        M = params_est(1);
        sigmaAdd = params_est(2);
        sigmaShared = params_est(3);
        criterion = params_est(4);

    case 2 % NoSharedN: Multi + Additive (no shared var)
        M = params_est(1);
        sigmaAdd = params_est(2);
        % sigmaShared = 0;
        criterion = params_est(3);

    case 3 % NoMultiN: Additive + Shared (no Multi var)
        % Nmul = 0;
        sigmaAdd = params_est(1);
        sigmaShared = params_est(2);
        criterion = params_est(3);

    case 4 % NoPrivN: Multi + Shared (no private var)
        M = params_est(1);
        % sigmaAdd = 0;
        sigmaShared = params_est(2);
        criterion = params_est(3);

    otherwise
        error('fxn_predMetrics_v2: Unknown iModelB = %d', iModelB);
end

%% -------- 2. Extract data --------
DV_allT = data.IV(:); % nTrials x 1
resp_allT = data.resp(:); % 1 = YES, 0 = NO
iPair_allT = data.iPair(:); % pair id per trial

nTrials = numel(DV_allT);
if numel(resp_allT) ~= nTrials || numel(iPair_allT) ~= nTrials
    error('fxn_predMetrics_v2: Data fields IV, resp, iPair must have same length.');
end

%% -------- 3. Predict sigma and criterion (per trial) --------

% Constant variance (shared + independent)
varConst = sigmaShared^2 + sigmaAdd^2;

% Trial-wise total noise SD (constant + Multi)
sigma_pred_allT = sqrt(varConst + (M .* DV_allT).^2 ); % nTrials x 1

% Avoid exactly zero variance (mvncdf & normcdf can be unhappy)
sigma_pred_allT = max(sigma_pred_allT, 1e-6);

%% -------- 4. Predicted pYES per trial --------

pYES_pred_allT = lambda/2 + (1 - lambda) .* (1 - normcdf(criterion, DV_allT, sigma_pred_allT));

% Safety: clamp probabilities away from 0 and 1
pYES_pred_allT = min(max(pYES_pred_allT, eps), 1 - eps);

%% Compute fitted criterion in z units from predicted pHit and pFA
pHit_fit = mean(pYES_pred_allT(data.iPRS == 1));
pFA_fit  = mean(pYES_pred_allT(data.iPRS == 0));

[~, Cz_fit] = SX_sim06_SDT(pHit_fit, pFA_fit);

%% -------- 5. Build pair-wise structure for agreement (pA) --------

iPairUnik = unique(iPair_allT);
nPairs = numel(iPairUnik);

consistency_allPairs = nan(nPairs, 1); % 1 if two responses in the pair agree, 0 otherwise
pairIdx = nan(nPairs, 2); % indices of the two trials in each pair

for iUnik = 1:nPairs
    idx = find(iPair_allT == iPairUnik(iUnik));
    if numel(idx) ~= 2
        error('fxn_predMetrics_v2: Each iPair should have exactly 2 trials (found %d).', numel(idx));
    end

    indPassA = idx(1);
    indPassB = idx(2);

    pairIdx(iUnik,:) = [indPassA, indPassB];

    respA = resp_allT(indPassA);
    respB = resp_allT(indPassB);
    consistency_allPairs(iUnik) = (respA == respB);
end

%% -------- 6. Predicted pA from shared / independent decision variables --------

pA_pred_allPairs = nan(nPairs, 1);

% Regions of integration for the bivariate Gaussian:
% YES = DV - criterion > 0, NO = DV - criterion < 0
YY_lb = [0, 0]; YY_ub = [Inf, Inf]; % both YES
NN_lb = [-Inf, -Inf]; NN_ub = [0, 0]; % both NO

for iUnik = 1:nPairs
    indPassA = pairIdx(iUnik,1);
    indPassB = pairIdx(iUnik,2);

    % Center DVs relative to criterion:
    mu_passA = DV_allT(indPassA) - criterion;
    mu_passB = DV_allT(indPassB) - criterion;

    % Total variance per pass (constant (shared + private) + multiplicative)
    var_total_passA = varConst + (M * DV_allT(indPassA))^2;
    var_total_passB = varConst + (M * DV_allT(indPassB))^2;

    % Enforce strictly positive variances
    var_total_passA = max(var_total_passA, 1e-10);
    var_total_passB = max(var_total_passB, 1e-10);

    % Covariance comes only from the shared constant noise
    CovAB = sigmaShared^2;
    
    % --- CAP: enforce positive definiteness ---
    % Need CovAB^2 < varA * varB; sufficient to enforce CovAB < min(varA,varB)."
    CovAB = min(CovAB, 0.999 * sqrt(var_total_passA* var_total_passB));

    % Create the covariance matrix
    CovMtx = [ var_total_passA, CovAB; ...
        CovAB, var_total_passB ];

    % P(both YES)
    P_YY = mvncdf(YY_lb, YY_ub, [mu_passA, mu_passB], CovMtx);

    % P(both NO)
    P_NN = mvncdf(NN_lb, NN_ub, [mu_passA, mu_passB], CovMtx);

    % Total agreement probability for the pair
    pA_pred_allPairs(iUnik) = P_YY + P_NN;
end

% Clamp pA_pred away from 0 and 1 for log
pA_pred_allPairs = min(max(pA_pred_allPairs, eps), 1 - eps);