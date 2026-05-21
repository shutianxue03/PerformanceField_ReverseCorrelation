function Tfix = fxn_buildFixedTransformFromABS(e3D_allT, iPRS_allT, CST_allT, lambda_whiten, eps_whiten, flag_meanMode, flag_plot)

% Estimate fixed transform parameters from ABS trials only.
% Z-scores each channel, then optionally computes a whitening matrix Q.
% Mean/SD mode is controlled by flag_meanMode:
%   'abs_global' = global from ABS trials only (default; current behavior)
%   'all_global' = global from all trials
%   'abs_per_contrast' = ABS-only, separately per contrast (requires CST_allT)

if ~(ischar(flag_meanMode) || (isstring(flag_meanMode) && isscalar(flag_meanMode)))
    error('flag_meanMode must be a string: ''abs_global'', ''all_global'', or ''abs_per_contrast''.');
end

idxABS = (iPRS_allT == 0);
assert(any(idxABS), 'No ABS trials available for fixed transform estimation.');

nTrials_all = size(e3D_allT, 1);
assert(numel(iPRS_allT) == nTrials_all, 'iPRS_allT length must match size(e3D_allT,1).');

[nTrials_abs, nORI, nSF] = size(e3D_allT(idxABS, :, :));
e3D_allAbs = reshape(e3D_allT(idxABS, :, :), [nTrials_abs, nORI * nSF]);
e3D_allT = reshape(e3D_allT, [nTrials_all, nORI * nSF]);

%% Global ABS stats are always saved as a fallback/reference.
mu_abs_global = mean(e3D_allAbs, 1);
sigma_abs_global = std(e3D_allAbs, [], 1);
sigma_abs_global(~isfinite(sigma_abs_global) | sigma_abs_global < 1e-8) = 1e-8;

mu_abs = mu_abs_global;
sigma_abs = sigma_abs_global;
mu_abs_perCST = [];
sigma_abs_perCST = [];
cst_levels = [];

switch flag_meanMode
    case 'abs_global'
        % Global from ABS trials only.
        e3D_allAbs_z = (e3D_allAbs - mu_abs) ./ sigma_abs;

    case 'all_global'
        % Global from all trials.
        mu_abs = mean(e3D_allT, 1);
        sigma_abs = std(e3D_allT, [], 1);
        sigma_abs(~isfinite(sigma_abs) | sigma_abs < 1e-8) = 1e-8;
        e3D_allAbs_z = (e3D_allAbs - mu_abs) ./ sigma_abs;

    case 'abs_per_contrast'
        % Per-contrast from ABS trials only.
        assert(~isempty(CST_allT), 'flag_meanMode=''abs_per_contrast'' requires CST_allT.');
        assert(numel(CST_allT) == nTrials_all, 'CST_allT length must match size(e3D_allT,1).');

        cst_abs = CST_allT(idxABS);
        cst_levels = unique(cst_abs(:));
        cst_levels = cst_levels(isfinite(cst_levels));
        assert(~isempty(cst_levels), 'No finite contrast levels found in ABS trials.');

        nCst = numel(cst_levels);
        nChan = size(e3D_allAbs, 2);
        mu_abs_perCST = nan(nCst, nChan);
        sigma_abs_perCST = nan(nCst, nChan);
        e3D_allAbs_z = nan(size(e3D_allAbs));

        for iCst = 1:nCst
            idxThis = (cst_abs == cst_levels(iCst));
            Xi = e3D_allAbs(idxThis, :);
            mu_i = mean(Xi, 1);
            sigma_i = std(Xi, [], 1);
            sigma_i(~isfinite(sigma_i) | sigma_i < 1e-8) = 1e-8;

            mu_abs_perCST(iCst, :) = mu_i;
            sigma_abs_perCST(iCst, :) = sigma_i;
            e3D_allAbs_z(idxThis, :) = (Xi - mu_i) ./ sigma_i;
        end

    otherwise
        error('Unknown flag_meanMode=%s. Use ''abs_global'', ''all_global'', or ''abs_per_contrast''.', flag_meanMode);
end

%% Obtain whitening matrix (Q)
[mu_cov_abs, Q] = fxn_getWhiteningMatrix(e3D_allAbs_z, lambda_whiten, eps_whiten);

% Store all fixed transform parameters in a struct for output.
Tfix = struct();
Tfix.mu_abs = mu_abs; % [1, nORI*nSF]
Tfix.sigma_abs = sigma_abs; % [1, nORI*nSF]
Tfix.mu_abs_global = mu_abs_global; % [1, nORI*nSF]
Tfix.sigma_abs_global = sigma_abs_global; % [1, nORI*nSF]
Tfix.flag_meanMode = flag_meanMode;
Tfix.mu_abs_perCST = mu_abs_perCST; % [nCst, nORI*nSF] when flag_meanMode=3
Tfix.sigma_abs_perCST = sigma_abs_perCST; % [nCst, nORI*nSF] when flag_meanMode=3
Tfix.cst_levels = cst_levels; % [nCst, 1] when flag_meanMode=3
Tfix.mu_cov_abs = mu_cov_abs; % [1, nORI*nSF]
Tfix.Q = Q;
Tfix.useWhiten = ~isempty(Q);

%% Visualize the fixed transform parameters for sanity check.
if flag_plot
    set(0, 'DefaultFigureVisible', 'on')
    figure('Position', [120 120 1300 650]);

    % Histograms of mu_abs
    subplot(2,3,1);
    histogram(Tfix.mu_abs(:), 40, 'FaceColor', [0.2 0.5 0.8]);
    xlabel('mu\_abs'); ylabel('Count');
    title(sprintf('Tfix.mu\_abs | expected: channel baseline means\nmode=%s', Tfix.flag_meanMode), 'Interpreter', 'none');

    % Histograms of sigma_abs
    subplot(2,3,2);
    histogram(Tfix.sigma_abs(:), 40, 'FaceColor', [0.3 0.7 0.4]);
    xlabel('sigma\_abs'); ylabel('Count');
    title('Tfix.sigma\_abs | expected: > 0 (floor at 1e-8)', 'Interpreter', 'none');

    % Histograms of mu_cov_abs
    subplot(2,3,3);
    histogram(Tfix.mu_cov_abs(:), 40, 'FaceColor', [0.7 0.45 0.2]);
    xline(0, 'r--', 'LineWidth', 1.5);
    xlabel('mu\_cov\_abs'); ylabel('Count');
    title('Tfix.mu\_cov\_abs | expected: concentrated near 0', 'Interpreter', 'none');

    % Visualization of Q
    subplot(2,3,4);
    if Tfix.useWhiten
        imagesc(Tfix.Q); colorbar; axis square;
        title('Tfix.Q (whitening matrix) | expected: symmetric-like, finite', 'Interpreter', 'none');
    else
        axis off;
        text(0.5, 0.5, 'Tfix.Q is empty\n(expected when whitening disabled)', ...
            'HorizontalAlignment', 'center', 'FontSize', 11);
        title('Tfix.Q | expected: empty if lambda\_whiten is NaN', 'Interpreter', 'none');
    end

    % Histogram of eigenvalues of Q
    subplot(2,3,5);
    if Tfix.useWhiten
        eigQ = eig((Tfix.Q + Tfix.Q')/2);
        histogram(eigQ, 40, 'FaceColor', [0.5 0.4 0.75]);
        xline(0, 'r--', 'LineWidth', 1.5);
        xlabel('eig(Q)'); ylabel('Count');
        title('eig(Tfix.Q) | expected: mostly positive (SPD whitening)', 'Interpreter', 'none');
    else
        axis off;
        text(0.5, 0.5, 'No eigenvalues to show\n(Q is empty)', ...
            'HorizontalAlignment', 'center', 'FontSize', 11);
        title('eig(Tfix.Q) | expected: not applicable when Q is empty', 'Interpreter', 'none');
    end

    % Bar plot of trial counts per contrast level (only relevant for per-contrast mode)
    subplot(2,3,6);
    if strcmp(Tfix.flag_meanMode, 'abs_per_contrast') && ~isempty(Tfix.cst_levels)
        nByCST = zeros(numel(Tfix.cst_levels), 1);
        cst_abs_tmp = CST_allT(iPRS_allT==0);
        for iC = 1:numel(Tfix.cst_levels)
            nByCST(iC) = sum(cst_abs_tmp == Tfix.cst_levels(iC));
        end
        bar(Tfix.cst_levels, nByCST, 'FaceColor', [0.2 0.6 0.6]);
        xlabel('Contrast level (ABS only)'); ylabel('nTrials');
        title('Per-contrast mode | expected: each used level has enough trials', 'Interpreter', 'none');
    else
        axis off;
        text(0.5, 0.5, 'Per-contrast stats not used\n(expected unless mode=abs_per_contrast)', ...
            'HorizontalAlignment', 'center', 'FontSize', 11);
        title('Per-contrast summary | expected: inactive in global modes', 'Interpreter', 'none');
    end

    sgtitle(sprintf('Tfix diagnostics | mode=%s | useWhiten=%d | lambda=%.3g', Tfix.flag_meanMode, Tfix.useWhiten, lambda_whiten), 'Interpreter', 'none');
    set(0, 'DefaultFigureVisible', 'off')
end

end % end of this function

%%
function [mu_cov, Q] = fxn_getWhiteningMatrix(e3D_allT_zscored, lambda_whiten, eps_whiten)

% Reshape e3D_allT_norm to 2D: trials x (ORI*SF)
[nTrials_cov, nORI_cov, nSF_cov] = size(e3D_allT_zscored);
e3D_allT_norm_vec = reshape(e3D_allT_zscored, [nTrials_cov, nORI_cov * nSF_cov]);

% Obtain covariance Sigma
mu_cov = mean(e3D_allT_norm_vec, 1);
Sigma = cov(e3D_allT_norm_vec, 1); % Sigma is MN x MN (M=nORI_cov, N=nSF_cov)

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
