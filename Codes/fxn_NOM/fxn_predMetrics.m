function [pYES_pred_allT, pA_pred_allPairs, consistency_allPairs] = fxn_predMetrics(iModelB, params_est, data, c_zscore)
% compute_pMetrics_sharedDV
% Shared core that computes:
%   - sigma_pred (trial-wise DV SD)
%   - criterion_IV (trial-wise criterion in IV units)
%   - pYES_pred (trial-wise P(YES))
%   - pA_pair (pair-wise P(agreement))
%   - pA_trial (trial-wise P(agreement), same value for both passes in a pair)
%
% INPUTS
%   IV        : nTrials x 1 internal variable
%   iPair     : nTrials x 1 pair ID for double-pass
%   Nmul      : multiplicative noise coefficient
%   SDadd     : additive noise (constant)
%   rhoDV     : decision-variable correlation across passes
%   lambda    : lapse rate (can be 0 or empty)
%   c_zscore  : criterion in z units (scaled by sigma_pred and added to median(IV))
%
% OUTPUTS
%   pYES_pred   : nTrials x 1 predicted P(YES)
%   pA_trial    : nTrials x 1 predicted P(agreement) per trial
%   pA_pair     : nPairs x 1 predicted P(agreement) per pair
%   sigma_pred  : nTrials x 1 predicted DV SD
%   criterion_IV: nTrials x 1 criterion in IV units

%% -------- 1. Parse model parameters (noise, lapse, rhoDV) --------
lambda = 0; % Set lapse rate at 0

switch iModelB
    case 1 % Full model: induced and constant noise
        Nmul = params_est(1);
        SDadd = params_est(2);
        rhoDV = params_est(3);

    case 2 % No Rho
        Nmul = params_est(1);
        SDadd = params_est(2);
        rhoDV = 0;

    case 3 % No induced noise
        Nmul = 0;
        SDadd = params_est(1);
        rhoDV = params_est(2);

    case 4 % No induced noise or rho
        Nmul = 0;
        SDadd = params_est(1);
        rhoDV = 0;

    case 5 % No constant noise
        Nmul = params_est(1);
        SDadd = 0;
        rhoDV = params_est(2);

    case 6 % No constant noise or rho
        Nmul = params_est(1);
        SDadd = 0;
        rhoDV = 0;

    case 7 % No Noise [the worst model]
        Nmul = 0;
        SDadd = 0;
        rhoDV = params_est(1);

    otherwise
        error('fxn_getError_v7: Unknown iModelB = %d', iModelB);
end

%% -------- 2. Extract data --------
IV_allT = data.IV(:); % nTrials x 1
resp_allT = data.resp(:); % 1 = YES, 0 = NO
iPair_allT = data.iPair(:); % pair id per trial

nTrials = numel(IV_allT);
if numel(resp_allT) ~= nTrials || numel(iPair_allT) ~= nTrials
    error('fxn_getError_v7: Data fields IV, resp, iPair must have same length.');
end

%% -------- 3. Predict sigma and criterion (per trial) --------

% Trial-wise total noise SD
sigma_pred_allT = sqrt((IV_allT .* Nmul).^2 + SDadd.^2); % nTrials x 1

% Avoid exactly zero variance (mvncdf & normcdf can be unhappy)
sigma_pred_allT = max(sigma_pred_allT, 1e-6);

% Criterion in IV units (could be changed to a scalar if desired)
criterion_allT = median(IV_allT) + c_zscore .* sigma_pred_allT;
% Alternatively: criterion_IV = c_zscore; % if you want fixed criterion

%% -------- 4. Predicted pYES per trial --------

pYES_pred_allT = lambda/2 + (1 - lambda) .* (1 - normcdf(criterion_allT, IV_allT, sigma_pred_allT));

% Safety: clamp probabilities away from 0 and 1 for log
pYES_pred_allT = min(max(pYES_pred_allT, eps), 1 - eps);

%% -------- 5. Build pair-wise structure for agreement (pA) --------

iPairUnik = unique(iPair_allT);
nPairs = numel(iPairUnik);

consistency_allPairs = nan(nPairs, 1); % 1 if two responses in the pair agree, 0 otherwise

% Create placeholder for indices of the two trials in each pair for pA_pred
pairIdx = nan(nPairs, 2);

for iUnik = 1:nPairs
    idx = find(iPair_allT == iPairUnik(iUnik));
    if numel(idx) ~= 2
        error('fxn_getError_v7: Each iPair should have exactly 2 trials (found %d).', numel(idx));
    end

    PassA = idx(1);
    PassB = idx(2);

    pairIdx(iUnik,:) = [PassA, PassB];

    respA = resp_allT(PassA);
    respB = resp_allT(PassB);
    consistency_allPairs(iUnik) = (respA == respB);
end

%% -------- 6. Predicted pA from correlated decision variables --------

pA_pred_allPairs = nan(nPairs, 1);

% Regions of integration for the bivariate Gaussian:
% YY (both responses = YES) corresponds to X_A > 0 (after shifting according to criterion) and X_B > 0
YY_lb = [0, 0];
YY_ub = [Inf, Inf];
% NN (both responses = No) corresponds to X_A < 0 and X_B < 0
NN_lb = [-Inf, -Inf];
NN_ub = [0, 0];

for iUnik = 1:nPairs
    PassA = pairIdx(iUnik,1);
    PassB = pairIdx(iUnik,2);

    % Center the decision variables relative to the criterion:
    % X = DV - criterion, so X > 0 means YES, X < 0 means NO.
    muA = IV_allT(PassA) - criterion_allT(PassA);
    muB = IV_allT(PassB) - criterion_allT(PassB);

    sigmaA_allPairs = sigma_pred_allT(PassA);
    sigmaB_allPairs = sigma_pred_allT(PassB);

    % Create the covariance matrix for the joint distribution of (X_A, X_B):
    % Variances: sigmaA^2, sigmaB^2
    % Covariance: rhoDV * sigmaA * sigmaB
    CovMtx = [ sigmaA_allPairs^2, rhoDV * sigmaA_allPairs * sigmaB_allPairs; ...
        rhoDV * sigmaA_allPairs*sigmaB_allPairs, sigmaB_allPairs^2 ];

    % Calculate the probability both passes yield YES (X_A > 0, X_B > 0)
    P_YY = mvncdf(YY_lb, YY_ub, [muA, muB], CovMtx);

    % Calculate the probability both passes yield No (X_A < 0, X_B < 0)
    P_NN = mvncdf(NN_lb, NN_ub, [muA, muB], CovMtx);

    % Combine
    pA_pred_allPairs(iUnik) = P_YY + P_NN;
end

% Clamp pA_pred away from 0 and 1 for log
pA_pred_allPairs = min(max(pA_pred_allPairs, eps), 1 - eps);
