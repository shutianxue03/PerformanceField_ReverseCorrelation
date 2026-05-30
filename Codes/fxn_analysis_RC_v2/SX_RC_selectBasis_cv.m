function output = SX_RC_selectBasis_cv(e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, opts)
% ============================================================
% Select ORI / SF basis counts, basis families, and ridge by
% cross-validated held-out negative log likelihood.
% Also conduct basis reconstruction (reversing the z-scoring of the basis predictors)
% Cached-Z + precomputed fold-standardization version.
% ============================================================

if nargin < 5
    error('Need e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, and opts');
end

%% Extract response vector
resp_allT = resp_allT(:);
nTrials = size(e3D_allT, 1);

if length(resp_allT) ~= nTrials
    error('resp_allT length does not match number of trials');
end

%% Extract hyperparameters
candidateORI = opts.nBasisORI(:)';
candidateSF = opts.nBasisSF(:)';
candidateBasisFamilyORI = opts.basisFamilyORI;
candidateBasisFamilySF = opts.basisFamilySF;
candidateBasisFamilyORI = cellstr(candidateBasisFamilyORI); % turn into cellstr if it's a single string for easier handling in comb loop
candidateBasisFamilySF = cellstr(candidateBasisFamilySF);
candidateBasisWidthORI = opts.basisWidthORI(:)';
candidateBasisWidthSF = opts.basisWidthSF(:)';
candidateAsymSF_rightLeftRatio = opts.asymSF_rightLeftRatio(:)';
candidateRidge = opts.Ridge(:)';

%% Create CV folds
nFolds = opts.nFolds;
nTrials_fold = numel(resp_allT);
foldID = nan(nTrials_fold,1);

if isfield(opts, 'rngSeed') && ~isempty(opts.rngSeed)
    rng(opts.rngSeed);
end

if opts.stratifyByResp % To ensure balanced response labels in CV folds, we stratify by response class before assigning fold IDs
    idx0 = find(resp_allT == 0);
    idx1 = find(resp_allT == 1);

    idx0 = idx0(randperm(numel(idx0)));
    idx1 = idx1(randperm(numel(idx1)));

    labels0 = repmat(1:nFolds, 1, ceil(numel(idx0)/nFolds));
    labels0 = labels0(1:numel(idx0));
    foldID(idx0) = labels0;

    labels1 = repmat(1:nFolds, 1, ceil(numel(idx1)/nFolds));
    labels1 = labels1(1:numel(idx1));
    foldID(idx1) = labels1;
else % No stratification, just random assignment to folds
    idx = randperm(nTrials_fold);
    labels = repmat(1:nFolds, 1, ceil(nTrials_fold/nFolds));
    labels = labels(1:nTrials_fold);
    foldID(idx) = labels;
end

foldID = foldID(:);

% precompute fold masks once
idxTrain_all = cell(nFolds,1);
idxTest_all  = cell(nFolds,1);
for iFold = 1:nFolds
    idxTest_all{iFold} = (foldID == iFold);
    idxTrain_all{iFold} = ~idxTest_all{iFold};
end

%% Build all non-ridge combinations
comb = struct([]);
iComb = 0;

for iFamORI = 1:numel(candidateBasisFamilyORI)
    for iFamSF = 1:numel(candidateBasisFamilySF)
        for iORIcand = 1:numel(candidateORI)
            for iSFcand = 1:numel(candidateSF)
                for iBWori = 1:numel(candidateBasisWidthORI)
                    for iBWsf = 1:numel(candidateBasisWidthSF)
                        for iAsym = 1:numel(candidateAsymSF_rightLeftRatio)
                            iComb = iComb + 1;
                            comb(iComb).nBasisORI = candidateORI(iORIcand);
                            comb(iComb).nBasisSF  = candidateSF(iSFcand);
                            comb(iComb).basisFamilyORI = char(candidateBasisFamilyORI{iFamORI});
                            comb(iComb).basisFamilySF  = char(candidateBasisFamilySF{iFamSF});
                            comb(iComb).basisWidthORI = candidateBasisWidthORI(iBWori);
                            comb(iComb).basisWidthSF  = candidateBasisWidthSF(iBWsf);
                            comb(iComb).asymSF_rightLeftRatio = candidateAsymSF_rightLeftRatio(iAsym);
                        end
                    end
                end
            end
        end
    end
end

nCombBase = numel(comb);
nRidge = numel(candidateRidge);
nRes = nCombBase * nRidge;

% preallocate result struct
res(nRes) = struct( ...
    'nBasisORI', [], ...
    'nBasisSF', [], ...
    'basisFamilyORI', '', ...
    'basisFamilySF', '', ...
    'basisWidthORI', [], ...
    'basisWidthSF', [], ...
    'asymSF_rightLeftRatio', [], ...
    'Ridge', [], ...
    'nLL_mean', [], ...
    'nLL_se', [] );

iRes = 0;

%% Main search
for iComb = 1:nCombBase

    % Extract the basis parameters for this combination and prepare the basis matrices
    opts_combo = opts;
    opts_combo.nBasisORI = comb(iComb).nBasisORI;
    opts_combo.nBasisSF  = comb(iComb).nBasisSF;
    opts_combo.basisFamilyORI = comb(iComb).basisFamilyORI;
    opts_combo.basisFamilySF  = comb(iComb).basisFamilySF;
    opts_combo.basisWidthORI = comb(iComb).basisWidthORI;
    opts_combo.basisWidthSF  = comb(iComb).basisWidthSF;
    opts_combo.asymSF_rightLeftRatio = comb(iComb).asymSF_rightLeftRatio;

    % Project all trials once using the shared basis-projector
    %--------------------------%
    Z_all = SX_RC_basisProject(e3D_allT, axis_ori_deg, axis_sf_log2, opts_combo);
    %--------------------------%

    % Precompute fold-specific standardized design matrices once
    Xtrain_all = cell(nFolds,1);
    Xtest_all  = cell(nFolds,1);
    rtrain_all = cell(nFolds,1);
    rtest_all  = cell(nFolds,1);

    for iFold = 1:nFolds
        idxTrain = idxTrain_all{iFold};
        idxTest  = idxTest_all{iFold};

        Z_train = Z_all(idxTrain, :);
        Z_test  = Z_all(idxTest, :);

        % fold-specific standardization and design matrix construction
        %--------------------------%
        [Xtrain_all{iFold}, muZ_fold, sdZ_fold] = build_design_matrix(Z_train, opts_combo.zscorePredictor);
        Xtest_all{iFold} = build_design_matrix(Z_test, opts_combo.zscorePredictor, muZ_fold, sdZ_fold);
        %--------------------------%

        rtrain_all{iFold} = resp_allT(idxTrain);
        rtest_all{iFold}  = resp_allT(idxTest);
    end

    % Now sweep ridge cheaply
    for iRidge = 1:nRidge
        Ridge = candidateRidge(iRidge);
        nLL_test_allFolds = nan(nFolds,1);

        % optional warm start across folds is usually not worth it;
        % warm start across ridge values within a fold would require
        % storing one beta per fold. Keep simple first.
        for iFold = 1:nFolds

            % Fit on train, predict on test using precomputed beta from precomputed X
            %--------------------------%
            beta = fit_ridge_glm_from_X(Xtrain_all{iFold}, rtrain_all{iFold}, Ridge, opts.link, opts.maxIter, opts.tol);
            %--------------------------%

            % Predict on test set
            %--------------------------%
            pred = predict_ridge_glm_from_X(Xtest_all{iFold}, rtest_all{iFold}, beta, opts.link);
            %--------------------------%

            % Store the test nLL for this fold
            nLL_test_allFolds(iFold) = pred.nLL;
        end

        % Store the results for this combination and ridge value
        iRes = iRes + 1;
        res(iRes).nBasisORI = comb(iComb).nBasisORI;
        res(iRes).nBasisSF  = comb(iComb).nBasisSF;
        res(iRes).basisFamilyORI = comb(iComb).basisFamilyORI;
        res(iRes).basisFamilySF  = comb(iComb).basisFamilySF;
        res(iRes).basisWidthORI = comb(iComb).basisWidthORI;
        res(iRes).basisWidthSF  = comb(iComb).basisWidthSF;
        res(iRes).asymSF_rightLeftRatio = comb(iComb).asymSF_rightLeftRatio;
        res(iRes).Ridge = Ridge;
        res(iRes).nLL_mean = mean(nLL_test_allFolds, 'omitnan');
        res(iRes).nLL_se   = std(nLL_test_allFolds, 'omitnan') / sqrt(sum(isfinite(nLL_test_allFolds)));
    end % iRidge
end % iComb

%% Choose the best by lowest mean nLL
[~, idxBest] = min([res.nLL_mean]);
bestRes = res(idxBest);

%% Store and refit the best model on the full dataset (no CV) to get the final template and training nLL for reference
opts_best = opts;
opts_best.nBasisORI = bestRes.nBasisORI;
opts_best.nBasisSF  = bestRes.nBasisSF;
opts_best.basisFamilyORI = bestRes.basisFamilyORI;
opts_best.basisFamilySF  = bestRes.basisFamilySF;
opts_best.basisWidthORI = bestRes.basisWidthORI;
opts_best.basisWidthSF  = bestRes.basisWidthSF;
opts_best.asymSF_rightLeftRatio = bestRes.asymSF_rightLeftRatio;
opts_best.Ridge = bestRes.Ridge;

%--------------------------%
out_best = SX_RC_fit_smoothBasis(e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, opts_best);
%--------------------------%

%% Store the outputs of the best fit
output = [];

output.template2D = out_best.template2D;
output.nLL = out_best.nLL; % This is training nLL from the best fit, so only secondary
output.deviance = out_best.deviance;
output.pseudoR2_Tjur = out_best.pseudoR2_Tjur;
output.bestRes = bestRes; % GoF from CV, more informative for model selection than the training nLL

% Store the selected hyperparameters for reference
output.nBasisORI = opts_best.nBasisORI;
output.nBasisSF = opts_best.nBasisSF;
output.basisFamilyORI = opts_best.basisFamilyORI;
output.basisFamilySF = opts_best.basisFamilySF;
output.basisWidthORI = opts_best.basisWidthORI;
output.basisWidthSF = opts_best.basisWidthSF;
output.asymSF_rightLeftRatio = opts_best.asymSF_rightLeftRatio;
output.Ridge = opts_best.Ridge;

%% Visualize GoF
% set(0, 'DefaultFigureVisible', 'on')
% bestPerFam = plot_basis_cv_heatmap(res);
%
% % Recovered template
% % Ideal template is projected into selected basis functions
% plot_best_template_marginalized(out_best.template2D, axis_ori_deg, axis_sf_log2, opts.template_ideal, bestRes);
%
% % close all

end

%% Helper: standardize predictors and build design matrix
function [X, muZ, sdZ] = build_design_matrix(Z, zscorePredictor, muZ, sdZ)

if nargin < 3
    muZ = mean(Z, 1);
    sdZ = std(Z, 0, 1);
    sdZ(sdZ < 1e-8) = 1;
end

if zscorePredictor
    Zfit = (Z - muZ) ./ sdZ;
else
    Zfit = Z;
    if nargin < 3
        muZ = zeros(1, size(Z,2));
        sdZ = ones(1, size(Z,2));
    end
end

X = [ones(size(Zfit,1),1), Zfit];
end

%% Helper: fit GLM from precomputed X
function beta = fit_ridge_glm_from_X(X, resp_allT, Ridge, linkName, maxIter, tol)

resp_allT = resp_allT(:);
p = size(X,2);

beta = zeros(p,1);
ridgeMat = diag([0; ones(p-1,1)]) * Ridge;
dev_prev = Inf;

for iIter = 1:maxIter
    eta = X * beta;
    eta = max(min(eta, 8), -8);

    switch lower(linkName)
        case 'logit'
            mu = 1 ./ (1 + exp(-eta));
            gprime = mu .* (1 - mu);
        case 'probit'
            mu = normcdf(eta);
            gprime = normpdf(eta);
        otherwise
            error('Unknown link: use ''probit'' or ''logit''');
    end

    % Define bound mu and gprime to avoid numerical issues
    mu = min(max(mu, 1e-8), 1 - 1e-8);
    gprime = max(gprime, 1e-6);

    % Compute IRLS weights and adjusted response
    W = (gprime.^2) ./ (mu .* (1 - mu) + eps);
    W = max(W, 1e-8);

    % Adjusted response for IRLS
    z = eta + (resp_allT - mu) ./ (gprime + eps);

    % Update beta using the weighted least squares solution with ridge penalty
    WX = X .* sqrt(W);
    wz = sqrt(W) .* z;

    % Add ridge penalty to the normal equations
    beta_new = (WX' * WX + ridgeMat) \ (WX' * wz);

    % Compute deviance for convergence check
    dev = -2 * sum(resp_allT .* log(mu) + (1 - resp_allT) .* log(1 - mu));

    if abs(dev_prev - dev) < tol
        beta = beta_new;
        break;
    end

    % Update beta and deviance for next iteration
    beta = beta_new;
    dev_prev = dev;
end
end

%% Helper: predict from precomputed X
function pred = predict_ridge_glm_from_X(X, resp, beta, linkName)

resp = resp(:);

% Compute linear predictor
eta = X * beta;
eta = max(min(eta, 8), -8);

% Compute predicted probabilities based on the link function
switch lower(linkName)
    case 'logit'
        yhat = 1 ./ (1 + exp(-eta));
    case 'probit'
        yhat = normcdf(eta);
    otherwise
        error('Unknown link: use ''probit'' or ''logit''');
end

yhat = min(max(yhat, 1e-8), 1 - 1e-8);

pred = struct();
pred.yhat = yhat;
pred.nLL = -sum(resp .* log(yhat) + (1 - resp) .* log(1 - yhat));
pred.deviance = -2 * sum(resp .* log(yhat) + (1 - resp) .* log(1 - yhat));
pred.eta = eta;
end

%%
function out = SX_RC_fit_smoothBasis(e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, opts)

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

% build basis matrices for the selected hyperparameters (important to do this before any CV-based standardization to avoid data leakage)
[Bori, Bsf] = build_basis_matrices_from_opts(axis_ori_deg, axis_sf_log2, opts);

% Optional: keep these only if you still want them in output
centersORI_deg = [];
centersSF_log2 = [];
centersSF_cpd  = [];

Kori = size(Bori, 2);
Ksf = size(Bsf, 2);

% Project trials into basis space
Z = SX_RC_basisProject(e3D_allT, axis_ori_deg, axis_sf_log2, opts);

% Fit ridge-regularized binomial GLM via IRLS (iteratively reweighted least squares)
%--------------------------%
[X, muZ, sdZ] = build_design_matrix(Z, opts.zscorePredictor);
%--------------------------%

%--------------------------%
beta = fit_ridge_glm_from_X(X, resp_allT, opts.Ridge, opts.link, opts.maxIter, opts.tol);
%--------------------------%

% Final predictions on training data
predTrain = predict_ridge_glm_from_X(X, resp_allT, beta, opts.link);
%--------------------------%

% Reconstruct the recovered template
beta0 = beta(1);
theta_vec = beta(2:end);
theta_vec_unscaled = theta_vec ./ sdZ(:); % unscale by predictor SD to get back to original basis space
Theta = reshape(theta_vec_unscaled, [Kori, Ksf]);
template2D = Bori * Theta * Bsf';

% Fit summary
idx1 = resp_allT == 1;
idx0 = resp_allT == 0;
pseudoR2_Tjur = mean(predTrain.yhat(idx1)) - mean(predTrain.yhat(idx0));

% Store output
out = struct();
out.template2D = template2D;
out.Theta = Theta;
out.beta0 = beta0;
out.beta = beta;               % important for held-out prediction
out.Bori = Bori;
out.Bsf = Bsf;
out.centersORI_deg = centersORI_deg;
out.centersSF_cpd = centersSF_cpd;
out.centersSF_log2 = centersSF_log2;
out.yhat = predTrain.yhat;
out.nLL = predTrain.nLL;
out.deviance = predTrain.deviance;
out.pseudoR2_Tjur = pseudoR2_Tjur;
out.Z = Z;
out.muZ = muZ;
out.sdZ = sdZ;
out.axis_ori_deg = axis_ori_deg;
out.axis_sf_log2 = axis_sf_log2;
out.link = opts.link;
out.zscorePredictor = opts.zscorePredictor;
out.opts = opts;
end

%%
function Z = SX_RC_basisProject(e3D_allT, axis_ori_deg, axis_sf_log2, basisOpts)
% SX_RC_basisProject
% Project ORI x SF energy maps into the tensor-product basis used by SX_RC_selectBasis_cv.
% Basis preparation happens here; predSFkernel only evaluates per-kernel equations.
% Called by SX_RC_selectBasis_cv, which calls predSFkernel for basis preparation and evaluation.

[Bori, Bsf] = build_basis_matrices_from_opts(axis_ori_deg, axis_sf_log2, basisOpts);

nTrials = size(e3D_allT, 1);
Kori = size(Bori, 2);
Ksf = size(Bsf, 2);
Z = zeros(nTrials, Kori * Ksf);

for iTrial = 1:nTrials
    Ei = squeeze(e3D_allT(iTrial, :, :));
    Zi = Bori' * Ei * Bsf;
    Z(iTrial, :) = Zi(:)';
end
end

%%
function [Bori, Bsf] = build_basis_matrices_from_opts(axis_ori_deg, axis_sf_log2, opts)
%----------- ORI -----------
% Define centers for ORI basis functions
centersORI_deg = linspace(0, opts.oriPeriod_deg, opts.nBasisORI + 1);
centersORI_deg(end) = []; % remove the last one to avoid duplication at the period boundary

% Calculate spacing between centers (in degrees)
dist_ori = opts.oriPeriod_deg / opts.nBasisORI; %

% Determine the  width parameter
str_family = opts.basisFamilyORI;
switch lower(char(opts.basisFamilyORI))
    case 'circ_gaussian_basis'
        ori_param2 = dist_ori * opts.basisWidthORI;
    case 'von_mises_basis'
        sigma_rad = deg2rad(dist_ori * opts.basisWidthORI);
        ori_param2 = 1 / max(sigma_rad.^2, 1e-6);
    otherwise
        error('Unknown ORI basis family: %s', opts.basisFamilyORI);
end

% Compute ORI basis matrix
Bori = zeros(numel(axis_ori_deg), numel(centersORI_deg));
for iCenter = 1:numel(centersORI_deg)
    Bori(:, iCenter) = predSFkernel(axis_ori_deg, str_family, [centersORI_deg(iCenter), ori_param2, opts.oriPeriod_deg], 0);
end

%----------- SF -----------
% Define centers for SF basis functions in both log2 and linear (cpd) space
centersSF_log2 = linspace(min(axis_sf_log2), max(axis_sf_log2), opts.nBasisSF);
axis_sf_cpd = 2 .^ axis_sf_log2;
centersSF_cpd = logspace(log10(min(axis_sf_cpd)), log10(max(axis_sf_cpd)), opts.nBasisSF);

% Calculate spacing between centers in log2 and log10 space
dist_log2 = mean(diff(centersSF_log2));
dist_log10cpd = mean(diff(log10(centersSF_cpd)));

% Determine the parameters for SF basis functions
sf_family = opts.basisFamilySF;
switch lower(char(opts.basisFamilySF))
    case 'log2_gaussian'
        sf_centers = centersSF_log2;
        sf_param2 = dist_log2 * opts.basisWidthSF;
        sf_param3 = [];
    case 'asym_log2_gaussian'
        sf_centers = centersSF_log2;
        sf_param2 = dist_log2 * opts.basisWidthSF;
        sf_param3 = sf_param2 * opts.asymSF_rightLeftRatio;
    case 'log_parabola_basis'
        sf_centers = centersSF_cpd;
        sf_param2 = dist_log10cpd * 1.2;
        sf_param3 = [];
    case 'asym_log_parabola_basis'
        sf_centers = centersSF_cpd;
        sf_param2 = dist_log10cpd;
        sf_param3 = sf_param2 * opts.asymSF_rightLeftRatio;
    otherwise
        error('Unknown SF basis family: %s', opts.basisFamilySF);
end

% Compute SF basis matrix
Bsf = zeros(numel(axis_sf_log2), numel(sf_centers));
if isempty(sf_param3) % symmetric case
    for iCenter = 1:numel(sf_centers)
        Bsf(:, iCenter) = predSFkernel(axis_sf_log2, sf_family, [sf_centers(iCenter), sf_param2], 0);
    end
else % asymmetric case
    for iCenter = 1:numel(sf_centers)
        Bsf(:, iCenter) = predSFkernel(axis_sf_log2, sf_family, [sf_centers(iCenter), sf_param2, sf_param3], 0);
    end
end

%----------- Normalization -----------
Bori = normalize_columns_local(Bori);
Bsf = normalize_columns_local(Bsf);
end

%%
function B = normalize_columns_local(B)
nrm = sqrt(sum(B.^2, 1));
nrm(nrm == 0) = 1;
B = bsxfun(@rdivide, B, nrm);
end
