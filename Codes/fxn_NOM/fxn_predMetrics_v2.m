function [pYES_pred_allT, pA_pred_allPairs, consistency_allPairs] = fxn_predMetrics_v2(iModelB, params_est, data, c_zscore)
% fxn_predMetrics_v2
% Shared core that computes:
% - sigma_pred_allT (trial-wise DV SD)
% - criterion_allT (trial-wise criterion in IV units)
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
% Nmul : multiplicative (induced) noise coefficient
% SDidpdt : constant independent noise (across passes)
% SDshared : constant shared noise (across passes)
%
% DV correlation across passes is implied by SDshared, not fit as rho.

%% -------- 1. Parse model parameters (noise, lapse) --------
lambda = 0; % Set lapse rate at 0

% Default values (in case a component is absent in a reduced model)
Nmul = 1e-10;
sigma_idpdt = 1e-10;
sigma_shared = 1e-10;

switch iModelB
    case 1 % FullModel: Induced + constant independent + constant shared
        Nmul = params_est(1);
        sigma_idpdt = params_est(2);
        sigma_shared = params_est(3);

    case 2 % NoSharedN: Induced + constant independent (no shared noise)
        Nmul = params_est(1);
        sigma_idpdt = params_est(2);
        % sigma_shared = 0;

    case 3 % NoInducedN: constant independent + constant shared (no induced noise)
        % Nmul = 0;
        sigma_idpdt = params_est(1);
        sigma_shared = params_est(2);

    case 4 % NoIdpdtN: Induced + constant shared (no independent constant noise)
        Nmul = params_est(1);
        % sigma_idpdt = 0;
        sigma_shared = params_est(2);

    case 5 % JustIdpdtN: constant independent only
        Nmul = 0;
        sigma_idpdt = params_est(1);
        % sigma_shared = 0;

    case 6 % JustInducedN: induced only
        Nmul = params_est(1);
        % sigma_idpdt = 0;
        % sigma_shared = 0;

    case 7 % JustSharedM: constant shared only
        % Nmul = 0;
        % sigma_idpdt = 0;
        sigma_shared = params_est(1);

    otherwise
        error('fxn_predMetrics_v2: Unknown iModelB = %d', iModelB);
end

%% -------- 2. Extract data --------
IV_allT = data.IV(:); % nTrials x 1
resp_allT = data.resp(:); % 1 = YES, 0 = NO
iPair_allT = data.iPair(:); % pair id per trial

nTrials = numel(IV_allT);
if numel(resp_allT) ~= nTrials || numel(iPair_allT) ~= nTrials
    error('fxn_predMetrics_v2: Data fields IV, resp, iPair must have same length.');
end

%% -------- 3. Predict sigma and criterion (per trial) --------

% Constant variance (shared + independent)
sigma_const = sigma_shared^2 + sigma_idpdt^2;

% Trial-wise total noise SD (constant + induced)
sigma_pred_allT = sqrt( sigma_const + (Nmul .* IV_allT).^2 ); % nTrials x 1

% Avoid exactly zero variance (mvncdf & normcdf can be unhappy)
sigma_pred_allT = max(sigma_pred_allT, 1e-6);

% Criterion in IV units (location-dependent via sigma_pred_allT)
criterion_allT = median(IV_allT) + c_zscore .* sigma_pred_allT;
% Alternatively: criterion_allT = c_zscore * ones(size(IV_allT));

%% -------- 4. Predicted pYES per trial --------

pYES_pred_allT = lambda/2 + (1 - lambda) .* ...
    (1 - normcdf(criterion_allT, IV_allT, sigma_pred_allT));

% Safety: clamp probabilities away from 0 and 1
pYES_pred_allT = min(max(pYES_pred_allT, eps), 1 - eps);

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

    PassA = idx(1);
    PassB = idx(2);

    pairIdx(iUnik,:) = [PassA, PassB];

    respA = resp_allT(PassA);
    respB = resp_allT(PassB);
    consistency_allPairs(iUnik) = (respA == respB);
end

%% -------- 6. Predicted pA from shared / independent decision variables --------

pA_pred_allPairs = nan(nPairs, 1);

% Regions of integration for the bivariate Gaussian:
% YES = DV - criterion > 0, NO = DV - criterion < 0
YY_lb = [0, 0]; YY_ub = [Inf, Inf]; % both YES
NN_lb = [-Inf, -Inf]; NN_ub = [0, 0]; % both NO

for iUnik = 1:nPairs
    PassA = pairIdx(iUnik,1);
    PassB = pairIdx(iUnik,2);

    % Center DVs relative to criterion:
    mu_passA = IV_allT(PassA) - criterion_allT(PassA);
    mu_passB = IV_allT(PassB) - criterion_allT(PassB);

    % Total variance per pass: constant (shared + independent) + induced
    IV_passA = IV_allT(PassA);
    IV_passB = IV_allT(PassB);

    sigma_passA2 = sigma_const + (Nmul * IV_passA)^2;
    sigma_passB2 = sigma_const + (Nmul * IV_passB)^2;

    % Enforce strictly positive variances
    sigma_passA2 = max(sigma_passA2, 1e-10);
    sigma_passB2 = max(sigma_passB2, 1e-10);

    % Covariance comes only from the shared constant noise
    % (independent constant + induced noise are independent across passes)
    CovAB = sigma_shared^2;

    % --- CAP: enforce positive definiteness ---
    % The shared variance (CovAB^2 ) has to be smaller than the total
    % variance (sigmaA2^2 * sigmaB2^2)
    CovAB_max = 0.999 * sqrt(sigma_passA2 * sigma_passB2); % 0.999 to stay strictly inside
    if abs(CovAB) > CovAB_max
        CovAB = sign(CovAB) * CovAB_max;
    end

    % Create the covariance matrix
    CovMtx = [ sigma_passA2, CovAB; ...
        CovAB, sigma_passB2 ];

    % P(both YES)
    P_YY = mvncdf(YY_lb, YY_ub, [mu_passA, mu_passB], CovMtx);

    % P(both NO)
    P_NN = mvncdf(NN_lb, NN_ub, [mu_passA, mu_passB], CovMtx);

    % Total agreement probability for the pair
    pA_pred_allPairs(iUnik) = P_YY + P_NN;
end

% Clamp pA_pred away from 0 and 1 for log
pA_pred_allPairs = min(max(pA_pred_allPairs, eps), 1 - eps);