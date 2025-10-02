function [kernel2D] = SX_RC(filtersSF_all, filtersOri_all, e3D_allT, resp_allT, flag_GLM)
% SX_RC: Computes regression model and metrics for each orientation and spatial frequency channel
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
%   flag_GLM - multi/univariate GLM

% OUTPUTS:
%   kernel2D         - 2D matrix of GLM coefficients (ORI x SF) for predictors


% Set default link function to 'probit' if not specified
fxnLink = 'probit';

% Extract the number of trials, orientations, and spatial frequencies
nTrials = length(resp_allT);
nSF = length(filtersSF_all);
nORI = length(filtersOri_all);

opts = statset('glmfit');
opts.MaxIter = 1000;  % Increase from default (100)

switch flag_GLM
    case 'multivariate'
        %% multivariate GLM
        beta_mtx = glmfit(reshape(e3D_allT, nTrials, nORI*nSF), resp_allT, 'binomial', 'Link', fxnLink);
        kernel2D = reshape(beta_mtx(2:end), nORI, nSF);

    case 'univariate'
        kernel2D = nan(nORI, nSF);               % Coefficients for each ORI and SF
        for iORI = 1:nORI
            parfor iSF = 1:nSF
                % Fit a GLM to predict response using specified link function
                beta = glmfit(e3D_allT(:, iORI, iSF), resp_allT, 'binomial', 'Link', fxnLink, 'Options', opts);
                % Store the intercept and coefficient (kernel) values
                kernel2D(iORI, iSF) = beta(2);
            end % End iSF loop
        end % End iORI loop
end
