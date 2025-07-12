
function nLL = fxn_getError_v5(iModelB, params_est, data, flag_plot)
% To calculate negative log-likelihood (nLL) for different models.

% INPUT:
%   iModelB - Model type for internal noise structure and decision criterion
%   params_est - Estimated parameters for model
%   data - Data structure with fields for internal variables (IV) and responses
%   flag_plot - Flag to enable or disable plotting (unused in current code)

%% Model setup: Define parameters based on model type (iModelB)
switch iModelB
    case 1 % lapse rate, additive noise, criterion
        N_mul = 0;
        lambda = params_est(1);
        SD_add = params_est(2);
        criterion = params_est(3);
    case 2 % lapse rate, multiplicative noise, criterion
        SD_add = 0;
        lambda = params_est(1);
        N_mul = params_est(2);
        criterion = params_est(3);
    case 3 % lapse rate, criterion (no internal noise)
        N_mul=0;
        SD_add=0;
        lambda = params_est(1);
        criterion = params_est(2);
    case 4 % additive noise, criterion
        N_mul=0;
        SD_add=params_est(1);
        criterion = params_est(2);
    case 5 % multiplicative noise, criterion
        N_mul=params_est(1);
        SD_add=0;
        criterion = params_est(2);
end

%% Extract data variables
IV = data.IV;
resp_data = data.resp;
iPair = data.iPair;

%% Predict metrics
% Compute predicted standard deviation (sigma) for each trial
sigma_pred = sqrt((IV*N_mul).^2 + SD_add^2);

% Calculate probability of responding "Present" (pYES)
if iModelB<=3
    pYES_pred = lambda/2+(1-lambda)*(1-normcdf(criterion, IV, sigma_pred));
else
    pYES_pred = 1-normcdf(criterion, IV, sigma_pred);
end

% Calculate probablity of responding consistently (pA)
iPairUnik = unique(iPair); % Unique pair identifiers
nUnik = length(iPairUnik); % Number of unique pairs
pYES_pred_pair = nan(nUnik, 1); % preallocate predicted probabilities for pairs
consistency = pYES_pred_pair; % preallocate consistency for pairs (1=consistent response, 0=inconsistent)
for iUnik = 1:nUnik
    iTrialAB = find(iPair == iPairUnik(iUnik)); % Find trials in each pair
    respA = resp_data(iTrialAB(1)); % Response of the first trial in this pair
    respB = resp_data(iTrialAB(2)); % Response of the second trial in this pair
    consistency(iUnik) = respA == respB; % Check if responses in this pair are consistent
    pYES_pred_A = pYES_pred(iTrialAB(1)); pYES_pred_B = pYES_pred(iTrialAB(2));assert(pYES_pred_A == pYES_pred_B, 'Predicted probabilities for pair trials should match');
    pYES_pred_pair(iUnik) = pYES_pred_A; % Store predicted probability for pair
end

% Predict the probability of paired responses being consistent
pA_pred = pYES_pred_pair.^2 + (1-pYES_pred_pair).^2;

% Calculate nLL by fitting predicted pYES to binary responses
if iModelB<=3
    % Use lapse rate for models with iModelB <= 3
    nLL_YES = -sum(resp_data.* log(pYES_pred) + (1 - resp_data) .* log(1 - pYES_pred));
else
    % Avoid zero probabilities by adding a small value (eps)
    nLL_YES = -sum(resp_data.* log(pYES_pred+eps) + (1 - resp_data) .* log(1 - pYES_pred+eps));
end

% Set final nLL output
nLL = nLL_YES;

%% BELOW ARE NOT IN USE
% Calculate nLL by fitting predicted pA to paired consistency responses
nLL_consistent = -sum(consistency.* log(pA_pred) + (1 - consistency) .* log(1 - pA_pred));

% nLL = nLL_consistent;

% Factorial log function for nLL normalization (currently unused in nLL)
fxn_factorialLog = @(m, n) sum(log(1:m)) - sum(log(1:n)) - sum(log(1:(m - n)));
% Example of full nLL calculation with factorial terms
% nLL = nLL_YES + fxn_factorialLog(length(resp_data), sum(resp_data)) + nLL_consistent + fxn_factorialLog(length(consistency), sum(consistency));
