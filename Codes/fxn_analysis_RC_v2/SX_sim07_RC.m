function [kernel2D, intercept2D, R2, R2_Tjur, pValues, pCat, yfit] = SX_sim07_RC(filtersSF_all, filtersOri_all, e3D_allT, resp_allT, fxnLink)
% SX_sim07_RC: Computes regression model and metrics for each orientation and spatial frequency channel
%
% This function applies a generalized linear model (GLM) to calculate kernel coefficients,
% intercepts, R-squared values, and additional statistical metrics for each combination of
% orientation and spatial frequency (SF) filters across trials. By default, it uses a probit
% link function for binary response data.
%
% INPUTS:
%   filtersSF_all    - Array of spatial frequency filters
%   filtersOri_all   - Array of orientation filters
%   e3D_allT              - 3D matrix of energy values (nTrials x ORI x SF)
%   resp_allT             - Vector of binary responses for each trial
%   fxnLink          - Link function for the GLM ('probit' or 'logit') [optional, default='probit']

% OUTPUTS:
%   kernel2D         - 2D matrix of GLM coefficients (ORI x SF) for predictors
%   intercept2D      - 2D matrix of intercept values (ORI x SF) from the GLM
%   R2               - 2D matrix of standard R-squared values for model fit (ORI x SF)
%   R2_Tjur          - 2D matrix of Tjur R-squared values, representing model quality (ORI x SF)
%   pValues          - 3D matrix of p-values for intercept and slope (ORI x SF x 2)
%   pCat             - 2D matrix of categorical prediction accuracy (ORI x SF)
%   yfit             - 3D matrix of fitted responses for each trial (trials x ORI x SF)

% Set default link function to 'probit' if not specified
if nargin < 5, fxnLink = 'probit'; end

% Extract the number of trials, orientations, and spatial frequencies
nTrials = length(resp_allT);
nSF = length(filtersSF_all);
nORI = length(filtersOri_all);

opts = statset('glmfit');
opts.MaxIter = 1000;  % Increase from default (100)

%% multivariate GLM
% % x: n x p; y: nx1; beta: (p+1) x 1
x = nan(nTrials, nORI * nSF); ii=1; indORI = nan(nORI * nSF, 1); indSF = nan(nORI * nSF, 1);
% reorganize
for iORI = 1:nORI, for iSF = 1:nSF, x(:, ii) = e3D_allT(:, iORI, iSF); indORI(ii) = iORI; indSF(ii) = iSF; ii=ii+1; end, end
%------------------------------------------------------
beta = glmfit(x, resp_allT, 'binomial', 'Link', fxnLink);
beta_ = glmfit(zscore(x), resp_allT,  'binomial');
%------------------------------------------------------
beta = beta(2:end); % the 1st element is the intercept
kernel2D = nan(nORI, nSF); for ii=1:nORI*nSF, kernel2D(indORI(ii), indSF(ii)) = beta(ii); end % reshape

%% Initialize empty containers for output variables

kernel2D = nan(nORI, nSF);               % Coefficients for each ORI and SF
intercept2D = kernel2D;                  % Intercepts for each ORI and SF
R2 = kernel2D;                           % Standard R-squared values
R2_Tjur = kernel2D;                      % Tjur R-squared values
pValues = nan(nORI, nSF, 2);             % p-values for intercept and slope
pCat = kernel2D;                         % Prediction accuracy for categories
yfit = nan(nTrials, nORI, nSF);          % Fitted responses

% Loop through each orientation and spatial frequency
warning('off', 'stats:glmfit:IterationLimit');

for iORI = 1:nORI
    parfor iSF = 1:nSF

        % Extract predictor variable for the current orientation and SF
        x = e3D_allT(:, iORI, iSF);
        try
            % Fit a GLM to predict response using specified link function
            [beta, ~, stats] = glmfit(x, resp_allT, 'binomial', 'Link', fxnLink, 'Options', opts);
            % for super easy trials, The GLM can draw a hyperplane that
            % separates all 1s and 0s perfectly. As a result, the estimated
            % slope (beta(2)) goes to ±∞ to make the predicted
            % probabilities exactly 0 or 1. 'LikelihoodPenalty' adds a
            % Bayesian prior (Jeffreys prior) that pulls coefficients back
            % from ±∞, giving finite and well-behaved estimates.

            % Store the intercept and coefficient (kernel) values
            intercept2D(iORI, iSF) = beta(1);
            kernel2D(iORI, iSF) = beta(2);

            % Calculate standard R-squared based on residuals
            R2(iORI, iSF) = 1 - var(stats.resid) / var(resp_allT);

            % Store p-values for intercept and slope
            pValues(iORI, iSF, :) = stats.p;

            % Generate fitted values for the model
            yfit_ = glmval(beta, x, fxnLink);
            yfit(:, iORI, iSF) = yfit_;

            % Calculate additional metrics (Tjur R-squared and categorical prediction accuracy)
            SS_res = sumsqr(stats.resid);
            SS_mod = sumsqr(mean(resp_allT) - yfit(:, iORI, iSF));
            SS_total = sumsqr(mean(resp_allT) - resp_allT);
            R2_mod = SS_mod / SS_total;
            R2_res = 1 - SS_res / SS_total;
            R2_Tjur_ = (R2_mod + R2_res) / 2;
            pCat_ = (sum(resp_allT(yfit_ >= 0.5) == 1) + sum(resp_allT(yfit_ <= 0.5) == 0)) / length(resp_allT);

            % Store Tjur R-squared and categorical accuracy
            R2_Tjur(iORI, iSF) = R2_Tjur_;
            pCat(iORI, iSF) = pCat_;

            % Uncomment to enable a quick plot of regression results
            % quickPlot_regression
        catch ME
            warning('GLM failed for ORI=%d SF=%d: %s', iORI, iSF, ME.message)
            % Leave as NaN
        end
    end % End iSF loop
end % End iORI loop

