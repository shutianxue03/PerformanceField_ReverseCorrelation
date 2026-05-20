function Tfix = fxn_buildFixedTransformFromABS(e3D, iPRS, lambda_whiten, eps_whiten)

% Estimate fixed transform parameters from ABS trials only.
% Z-scores each channel by a single mean and SD (pooled across all ABS trials,
% ignoring contrast level), then optionally computes a whitening matrix Q.

idxABS = (iPRS == 0);
assert(any(idxABS), 'No ABS trials available for fixed transform estimation.');

[nTrials, nORI, nSF] = size(e3D(idxABS, :, :));
e3D_abs = reshape(e3D(idxABS, :, :), [nTrials, nORI * nSF]);

% One mean and SD per channel across all ABS trials
mu_abs = mean(e3D_abs, 1);
sigma_abs = std(e3D_abs, [], 1);
sigma_abs(~isfinite(sigma_abs) | sigma_abs < 1e-8) = 1e-8;

e3D_abs_z = (e3D_abs - mu_abs) ./ sigma_abs;

Q = fxn_getWhiteningMatrix(e3D_abs_z, lambda_whiten, eps_whiten);

Tfix = struct();
Tfix.mu_abs = mu_abs; % [1, nORI*nSF]
Tfix.sigma_abs = sigma_abs; % [1, nORI*nSF]
Tfix.Q = Q;

end

%%
function Q = fxn_getWhiteningMatrix(e3D_zscored, lambda_whiten, eps_whiten)

% Reshape e3D_norm to 2D: trials x (ORI*SF)
[nTrials_cov, nORI_cov, nSF_cov] = size(e3D_zscored);
e3D_norm_vec = reshape(e3D_zscored, [nTrials_cov, nORI_cov * nSF_cov]);

% Obtain covariance Sigma
Sigma = cov(e3D_norm_vec, 1); % Sigma is MN x MN (M=nORI_cov, N=nSF_cov)

% Calculate the whitening matrix Q via eigen-decomposition on Sigma
if isnan(lambda_whiten)
    Q = [];
else
    Sigma_shrink = (1 - lambda_whiten) * Sigma + lambda_whiten * mean(diag(Sigma)) * eye(size(Sigma, 1));
    [V, D] = eig((Sigma_shrink + Sigma_shrink') / 2);
    d = diag(D);
    d(d < eps_whiten) = eps_whiten;
    Q = V * diag(1 ./ sqrt(d)) * V';
end

end
