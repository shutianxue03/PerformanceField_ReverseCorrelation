function out = SX_RC_smoothBasis_circORI_logSF(e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, opts)

% Smooth-basis RC / GLM for ORI x SF energy templates
%
% Main idea:
% Instead of estimating one independent weight for every ORI x SF channel
% (which can become noisy and unstable), we force the template to live in a
% smaller, smoother space spanned by basis functions.
%
% Here:
%   - ORI is modeled with circular Gaussian bumps
%   - SF is modeled with Gaussian bumps in log2(SF)
%
% So the final template is NOT a free grid of weights.
% It is reconstructed from a small set of smooth basis coefficients.
%
% This is inspired by the logic of GAM / smooth-basis estimation:
%   1) represent the template in a smooth low-dimensional basis
%   2) fit binary responses with a GLM
%   3) reconstruct the smooth template surface
%
% Important:
% This code does NOT implement a full GAM roughness penalty
% (e.g., second-derivative penalty).
% The smoothness comes mainly from:
%   - restricting the template to a smooth basis
%   - using only a small number of basis functions
% The explicit penalty in this code is a ridge penalty on the fitted basis
% coefficients (see ridgeMat below).

if nargin < 5
    opts = struct();
end

%% ---------------- defaults ----------------
% link:
%   'probit' or 'logit' for the binomial GLM
% nBasisORI / nBasisSF:
%   number of basis functions along each dimension
% sigmaORI_deg / sigmaSF_log2:
%   width of the Gaussian basis functions
% oriPeriod_deg:
%   180 for orientation, 360 if direction is modeled instead
% zscorePredictor:
%   whether to z-score the projected basis predictors before fitting
% ridge:
%   ridge penalty strength on fitted coefficients
% maxIter / tol:
%   settings for IRLS optimization

axis_ori_deg = axis_ori_deg(:);
axis_sf_log2 = axis_sf_log2(:);

[nTrials, nORI, nSF] = size(e3D_allT);

if length(axis_ori_deg) ~= nORI
    error('axis_ori_deg length does not match nORI');
end
if length(axis_sf_log2) ~= nSF
    error('axis_sf_log2 length does not match nSF');
end
if length(resp_allT) ~= nTrials
    error('resp_allT length does not match nTrials');
end

resp_allT = resp_allT(:);

%% ---------------- ORI basis: circular Gaussian bumps ----------------
% Why circular Gaussian bumps?
% ORI is periodic: 0 deg and 180 deg should be neighbors, not endpoints.
% A circular basis respects that geometry and avoids boundary artifacts.
%
% centersORI_deg:
%   centers of the basis bumps around the circular ORI axis
centersORI_deg = linspace(0, opts.oriPeriod_deg, opts.nBasisORI + 1);
centersORI_deg(end) = []; % remove duplicate endpoint

% If user does not specify the width, choose it based on basis spacing.
% Wider bumps = smoother / lower-dimensional model.
if isempty(opts.sigmaORI_deg)
    if opts.nBasisORI > 1
        d = opts.oriPeriod_deg / opts.nBasisORI;
        sigmaORI_deg = d * 0.8;
    else
        sigmaORI_deg = opts.oriPeriod_deg / 4;
    end
else
    sigmaORI_deg = opts.sigmaORI_deg;
end

% Bori is [nORI x Kori]
% Each column is one circular Gaussian basis function over ORI.
Bori = circular_gaussian_basis(axis_ori_deg, centersORI_deg, sigmaORI_deg, opts.oriPeriod_deg);

% Normalize columns so basis functions have comparable scale.
Bori = normalize_columns(Bori);

%% ---------------- SF basis: Gaussian bumps in log2(SF) ----------------
% Why Gaussian bumps in log2(SF)?
% SF tuning is usually smoother and more interpretable on a logarithmic
% scale (octaves) than in linear cpd units.
%
% centersSF_log2:
%   basis centers along the log-SF axis
centersSF_log2 = linspace(min(axis_sf_log2), max(axis_sf_log2), opts.nBasisSF);

if isempty(opts.sigmaSF_log2)
    if opts.nBasisSF > 1
        d = mean(diff(centersSF_log2));
        sigmaSF_log2 = d * 0.8;
    else
        sigmaSF_log2 = range(axis_sf_log2) / 4 + eps;
    end
else
    sigmaSF_log2 = opts.sigmaSF_log2;
end

% Bsf is [nSF x Ksf]
% Each column is one Gaussian basis function over log-SF.
Bsf = gaussian_basis(axis_sf_log2, centersSF_log2, sigmaSF_log2);
Bsf = normalize_columns(Bsf);

% Just for readability/output: convert basis centers back to cpd
centersSF_cpd = 2 .^ centersSF_log2;

Kori = size(Bori, 2);
Ksf = size(Bsf, 2);
Ktot = Kori * Ksf;

%% ---------------- project each trial into tensor-product basis ----------------
% This is the key "subspace" step.
%
% Original trial data:
%   Ei = ORI x SF energy map for one trial
%
% Instead of fitting Ei directly channel by channel, we project Ei onto the
% tensor-product basis defined by Bori and Bsf:
%
%   Zi = Bori' * Ei * Bsf
%
% So Zi is a low-dimensional summary of Ei in the smooth basis space.
%
% If Ei is [nORI x nSF], then Zi is [Kori x Ksf].
% We then vectorize Zi into one predictor vector per trial.
Z = zeros(nTrials, Ktot);

for iTrial = 1:nTrials
    Ei = squeeze(e3D_allT(iTrial, :, :)); % [nORI x nSF]
    Zi = Bori' * Ei * Bsf;                % [Kori x Ksf]
    Z(iTrial, :) = Zi(:)';
end

%% ---------------- standardize predictors ----------------
% Standardizing helps optimization and makes ridge regularization more
% comparable across coefficients.
muZ = mean(Z, 1);
sdZ = std(Z, 0, 1);
sdZ(sdZ == 0) = 1;

if opts.zscorePredictor
    Zfit = (Z - muZ) ./ sdZ;
else
    Zfit = Z;
    muZ = zeros(1, Ktot);
    sdZ = ones(1, Ktot);
end

%% ---------------- fit ridge-regularized binomial GLM via IRLS ----------------
% Here we fit:
%
%   response ~ binomial(link^{-1}(beta0 + Zfit * beta))
%
% where beta are the coefficients on the low-dimensional smooth basis
% predictors.
%
% IMPORTANT:
% The explicit penalty here is RIDGE, not a true GAM roughness penalty.
%
% ridgeMat penalizes large coefficients:
%   penalty ~ ridge * sum(beta_j^2)
%
% This stabilizes the fit, especially when predictors are correlated.
%
% But it does NOT directly penalize wiggliness.
% Smoothness mainly comes from the basis restriction above.
X = [ones(nTrials,1), Zfit];
p = size(X,2);

beta = zeros(p,1);

% Do not penalize the intercept
ridgeMat = diag([0; ones(p-1,1)]) * opts.ridge;

dev_prev = Inf;

for iter = 1:opts.maxIter
    eta = X * beta;

    % Convert linear predictor to response probability
    switch lower(opts.link)
        case 'logit'
            mu = 1 ./ (1 + exp(-eta));
            gprime = mu .* (1 - mu);
        case 'probit'
            mu = normcdf(eta);
            gprime = normpdf(eta);
        otherwise
            error('Unknown link: use ''probit'' or ''logit''');
    end

    mu = min(max(mu, 1e-8), 1 - 1e-8);

    % IRLS working weights and pseudo-response
    W = (gprime.^2) ./ (mu .* (1 - mu) + eps);
    z = eta + (resp_allT - mu) ./ (gprime + eps);

    WX = bsxfun(@times, X, sqrt(W));
    wz = sqrt(W) .* z;

    % Ridge-regularized weighted least squares update
    beta_new = (WX' * WX + ridgeMat) \ (WX' * wz);

    % Binomial deviance
    dev = -2 * sum(resp_allT .* log(mu) + (1 - resp_allT) .* log(1 - mu));

    % if opts.verbose && (iter == 1 || mod(iter,10)==0)
    %     fprintf('Iter %d, deviance = %.6f\n', iter, dev);
    % end

    if abs(dev_prev - dev) < opts.tol
        beta = beta_new;
        break;
    end

    beta = beta_new;
    dev_prev = dev;
end

%% ---------------- final predictions ----------------
eta = X * beta;
switch lower(opts.link)
    case 'logit'
        yhat = 1 ./ (1 + exp(-eta));
    case 'probit'
        yhat = normcdf(eta);
end
yhat = min(max(yhat, 1e-8), 1 - 1e-8);

deviance = -2 * sum(resp_allT .* log(yhat) + (1 - resp_allT) .* log(1 - yhat));

%% ---------------- reconstruct template ----------------
% beta(2:end) are the fitted coefficients in basis space.
% These correspond to the tensor-product basis dimensions, not directly to
% ORI x SF channels.
beta0 = beta(1);
theta_vec = beta(2:end);

% Undo predictor scaling so the reconstructed template lives on the original
% projected-feature scale rather than the z-scored one.
theta_vec_unscaled = theta_vec ./ sdZ(:);

% Reshape back into [Kori x Ksf] coefficient matrix
Theta = reshape(theta_vec_unscaled, [Kori, Ksf]);

% Reconstruct smooth template in the original ORI x SF grid:
%
%   template2D = Bori * Theta * Bsf'
%
% So the final template is a weighted sum of smooth basis functions.
template2D = Bori * Theta * Bsf';

%% ---------------- fit summary ----------------
idx1 = resp_allT == 1;
idx0 = resp_allT == 0;
pseudoR2_Tjur = mean(yhat(idx1)) - mean(yhat(idx0));

%% ---------------- output ----------------
out = struct();
out.template2D = template2D;      % reconstructed smooth ORI x SF template
out.Theta = Theta;                % coefficients in basis space
out.beta0 = beta0;                % intercept
out.Bori = Bori;                  % ORI basis functions
out.Bsf = Bsf;                    % SF basis functions
out.centersORI_deg = centersORI_deg;
out.centersSF_cpd = centersSF_cpd;
out.centersSF_log2 = centersSF_log2;
out.yhat = yhat;                  % fitted response probabilities
out.deviance = deviance;
out.pseudoR2_Tjur = pseudoR2_Tjur;
out.Z = Z;                        % projected trial predictors before z-scoring
out.muZ = muZ;
out.sdZ = sdZ;
out.axis_sf_log2 = axis_sf_log2;
end

%% HELPERS

function B = gaussian_basis(x, centers, sigma)
% Standard Gaussian bump basis on a line
% Input:
%   x       : positions where basis is evaluated
%   centers : centers of basis bumps
%   sigma   : width of each bump
x = x(:);
centers = centers(:)';
B = exp(-0.5 * ((x - centers) ./ sigma).^2);
end

function B = circular_gaussian_basis(x_deg, centers_deg, sigma_deg, period_deg)
% Circular Gaussian bump basis
% Distance is computed on a circle, so 0 deg and 180 deg are neighbors
% when period_deg = 180.
x_deg = x_deg(:);
centers_deg = centers_deg(:)';
B = zeros(length(x_deg), length(centers_deg));

for k = 1:length(centers_deg)
    d = mod(x_deg - centers_deg(k) + period_deg/2, period_deg) - period_deg/2;
    B(:,k) = exp(-0.5 * (d ./ sigma_deg).^2);
end
end

function B = normalize_columns(B)
% Normalize each basis function to unit L2 norm
% This keeps one basis column from dominating just because it is larger in
% raw magnitude.
nrm = sqrt(sum(B.^2, 1));
nrm(nrm == 0) = 1;
B = bsxfun(@rdivide, B, nrm);
end