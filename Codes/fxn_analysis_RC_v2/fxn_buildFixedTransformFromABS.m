function Tfix = fxn_buildFixedTransformFromABS(e3D_allT, iPRS_allT, lambda_whiten, eps_whiten)
% FXN_BUILDFIXEDTRANSFORMFROMABS
%   Build a fixed (z-score + optional whitening) transform from ABS trials.

idxABS = (iPRS_allT == 0);
assert(any(idxABS), 'No ABS trials available for fixed transform estimation.');

nTrials_all = size(e3D_allT, 1);
assert(numel(iPRS_allT) == nTrials_all, 'iPRS_allT length must match size(e3D_allT,1).');

[nTrials_abs, nORI, nSF] = size(e3D_allT(idxABS, :, :));
e3D_allAbs = reshape(e3D_allT(idxABS, :, :), [nTrials_abs, nORI * nSF]);

%% Per-channel mean and SD from ABS trials only
mu_abs = mean(e3D_allAbs, 1);
sigma_abs = std(e3D_allAbs, [], 1);
sigma_abs(~isfinite(sigma_abs) | sigma_abs < 1e-8) = 1e-8;

% Z-score ABS trials (used for whitening matrix estimation only)
e3D_allAbs_z = (e3D_allAbs - mu_abs) ./ sigma_abs;

%% Optional whitening matrix from z-scored ABS trials
disableWhiten = isnan(lambda_whiten) || isequal(lambda_whiten, 99);
if disableWhiten
    mu_cov_abs = [];
    Q = [];
else
    [mu_cov_abs, Q] = fxn_getWhiteningMatrix(e3D_allAbs_z, lambda_whiten, eps_whiten);
end

%% Pack
Tfix = struct();
Tfix.mu_abs     = mu_abs;        % [1, nORI*nSF]
Tfix.sigma_abs  = sigma_abs;     % [1, nORI*nSF]
Tfix.mu_cov_abs = mu_cov_abs;    % [1, nORI*nSF]
Tfix.Q          = Q;             % [nChan, nChan] or []
Tfix.useWhiten  = ~disableWhiten;

end % fxn_buildFixedTransformFromABS

%% -----------------------------------------------------------------------
function [mu_cov, Q] = fxn_getWhiteningMatrix(X_z, lambda_whiten, eps_whiten)
% Build a symmetric whitening matrix from z-scored data via eigendecomp,
% with diagonal shrinkage and an eigenvalue floor.

mu_cov = mean(X_z, 1);
Sigma  = cov(X_z, 1);

if isnan(lambda_whiten)
    Q = [];
    return;
end

Sigma_shrink = (1 - lambda_whiten) * Sigma + lambda_whiten * mean(diag(Sigma)) * eye(size(Sigma, 1));
[V, D] = eig((Sigma_shrink + Sigma_shrink') / 2);
d = diag(D);
d(d < eps_whiten) = eps_whiten;
Q = V * diag(1 ./ sqrt(d)) * V';
end
