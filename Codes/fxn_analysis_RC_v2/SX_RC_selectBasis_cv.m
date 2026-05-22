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

opts = fill_default_opts(opts);

resp_allT = resp_allT(:);
nTrials = size(e3D_allT, 1);

if length(resp_allT) ~= nTrials
    error('resp_allT length does not match number of trials');
end

%% Extract (range of) hyperparameters
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

%% ---------------- create CV folds ----------------
nFolds = opts.nFolds;
foldID = make_cv_folds(resp_allT, nFolds, opts);

% precompute fold masks once
idxTrain_all = cell(nFolds,1);
idxTest_all  = cell(nFolds,1);
for iFold = 1:nFolds
    idxTest_all{iFold} = (foldID == iFold);
    idxTrain_all{iFold} = ~idxTest_all{iFold};
end

%% ---------------- build all non-ridge combinations ----------------
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

%% ---------------- main search ----------------
for iComb = 1:nCombBase

    opts_basis = opts;
    opts_basis.nBasisORI = comb(iComb).nBasisORI;
    opts_basis.nBasisSF  = comb(iComb).nBasisSF;
    opts_basis.basisFamilyORI = comb(iComb).basisFamilyORI;
    opts_basis.basisFamilySF  = comb(iComb).basisFamilySF;
    opts_basis.basisWidthORI = comb(iComb).basisWidthORI;
    opts_basis.basisWidthSF  = comb(iComb).basisWidthSF;
    opts_basis.asymSF_rightLeftRatio = comb(iComb).asymSF_rightLeftRatio;

    % project all trials once using the shared basis-projector
    Z_all = SX_RC_basisProject(e3D_allT, axis_ori_deg, axis_sf_log2, opts_basis);

    % precompute fold-specific standardized design matrices once
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
        [Xtrain_all{iFold}, Xtest_all{iFold}] = prepare_fold_design_matrices(Z_train, Z_test, opts.zscorePredictor);

        rtrain_all{iFold} = resp_allT(idxTrain);
        rtest_all{iFold}  = resp_allT(idxTest);
    end

    % now sweep ridge cheaply
    for iRidge = 1:nRidge
        ridge = candidateRidge(iRidge);
        nLL_test_allFolds = nan(nFolds,1);

        % optional warm start across folds is usually not worth it;
        % warm start across ridge values within a fold would require
        % storing one beta per fold. Keep simple first.
        for iFold = 1:nFolds
            % fit on train, predict on test using precomputed X
            beta = fit_ridge_glm_from_X(Xtrain_all{iFold}, rtrain_all{iFold}, ridge, opts.link, opts.maxIter, opts.tol);
            % Prediction on test set using the beta fitted on the train set of this fold
            pred = predict_ridge_glm_from_X(Xtest_all{iFold}, rtest_all{iFold}, beta, opts.link);
            % Store the test nLL for this fold
            nLL_test_allFolds(iFold) = pred.nLL;
        end

        iRes = iRes + 1;
        res(iRes).nBasisORI = comb(iComb).nBasisORI;
        res(iRes).nBasisSF  = comb(iComb).nBasisSF;
        res(iRes).basisFamilyORI = comb(iComb).basisFamilyORI;
        res(iRes).basisFamilySF  = comb(iComb).basisFamilySF;
        res(iRes).basisWidthORI = comb(iComb).basisWidthORI;
        res(iRes).basisWidthSF  = comb(iComb).basisWidthSF;
        res(iRes).asymSF_rightLeftRatio = comb(iComb).asymSF_rightLeftRatio;
        res(iRes).Ridge = ridge;
        res(iRes).nLL_mean = mean(nLL_test_allFolds, 'omitnan');
        res(iRes).nLL_se   = std(nLL_test_allFolds, 'omitnan') / sqrt(sum(isfinite(nLL_test_allFolds)));
    end % iRidge
end % iComb

%% ---------------- choose best by lowest mean nLL ----------------
[~, idxBest] = min([res.nLL_mean]);
bestRes = res(idxBest);

%% ---------------- refit on all data using full model ----------------
opts_best = opts;
opts_best.nBasisORI = bestRes.nBasisORI;
opts_best.nBasisSF  = bestRes.nBasisSF;
opts_best.basisFamilyORI = bestRes.basisFamilyORI;
opts_best.basisFamilySF  = bestRes.basisFamilySF;
opts_best.basisWidthORI = bestRes.basisWidthORI;
opts_best.basisWidthSF  = bestRes.basisWidthSF;
opts_best.asymSF_rightLeftRatio = bestRes.asymSF_rightLeftRatio;
opts_best.Ridge = bestRes.Ridge;

out_best = SX_RC_fit_smoothBasis(e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, opts_best);

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
%%
end


%% Helper: precompute standardized train/test X once per fold
function [Xtrain, Xtest] = prepare_fold_design_matrices(Ztrain, Ztest, zscorePredictor)

muZ = mean(Ztrain, 1);
sdZ = std(Ztrain, 0, 1);
sdZ(sdZ < 1e-8) = 1;

if zscorePredictor
    Ztrain_fit = (Ztrain - muZ) ./ sdZ;
    Ztest_fit  = (Ztest  - muZ) ./ sdZ;
else
    Ztrain_fit = Ztrain;
    Ztest_fit  = Ztest;
end

Xtrain = [ones(size(Ztrain_fit,1),1), Ztrain_fit];
Xtest  = [ones(size(Ztest_fit,1),1),  Ztest_fit];
end

%% Helper: fit GLM from precomputed X
function beta = fit_ridge_glm_from_X(X, resp_allT, ridge, linkName, maxIter, tol)

resp_allT = resp_allT(:);
p = size(X,2);

beta = zeros(p,1);
ridgeMat = diag([0; ones(p-1,1)]) * ridge;
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

    mu = min(max(mu, 1e-8), 1 - 1e-8);
    gprime = max(gprime, 1e-6);

    W = (gprime.^2) ./ (mu .* (1 - mu) + eps);
    W = max(W, 1e-8);

    z = eta + (resp_allT - mu) ./ (gprime + eps);

    WX = X .* sqrt(W);
    wz = sqrt(W) .* z;

    beta_new = (WX' * WX + ridgeMat) \ (WX' * wz);

    dev = -2 * sum(resp_allT .* log(mu) + (1 - resp_allT) .* log(1 - mu));

    if abs(dev_prev - dev) < tol
        beta = beta_new;
        break;
    end

    beta = beta_new;
    dev_prev = dev;
end
end

%% Helper: predict from precomputed X
function pred = predict_ridge_glm_from_X(X, resp, beta, linkName)

resp = resp(:);

eta = X * beta;
eta = max(min(eta, 8), -8);

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
end

%%
function out = SX_RC_fit_smoothBasis(e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, opts)

opts = fill_default_opts(opts);

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

%--------------------%
Bori = predSFkernel('make_basis_ori', axis_ori_deg, opts, 0);
Bsf  = predSFkernel('make_basis_sf', axis_sf_log2, opts, 0);
%--------------------%

% Optional: keep these only if you still want them in output
centersORI_deg = [];
centersSF_log2 = [];
centersSF_cpd  = [];

Kori = size(Bori, 2);
Ksf = size(Bsf, 2);
Ktot = Kori * Ksf;

% ---------------- project trials into basis space ----------------
Z = SX_RC_basisProject(e3D_allT, axis_ori_deg, axis_sf_log2, opts);

% ---------------- standardize predictors ----------------
muZ = mean(Z, 1);
sdZ = std(Z, 0, 1);
sdZ(sdZ < 1e-8) = 1;

if opts.zscorePredictor
    Zfit = (Z - muZ) ./ sdZ;
else
    Zfit = Z;
    muZ = zeros(1, Ktot);
    sdZ = ones(1, Ktot);
end

% ---------------- fit ridge-regularized binomial GLM via IRLS ----------------
X = [ones(nTrials,1), Zfit];
p = size(X,2);

beta = zeros(p,1);
ridgeMat = diag([0; ones(p-1,1)]) * opts.ridge;
dev_prev = Inf;

for iIter = 1:opts.maxIter

    eta = X * beta;
    eta = max(min(eta, 8), -8);

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
    gprime = max(gprime, 1e-6);

    W = (gprime.^2) ./ (mu .* (1 - mu) + eps);
    W = max(W, 1e-8);

    z = eta + (resp_allT - mu) ./ (gprime + eps);

    WX = bsxfun(@times, X, sqrt(W));
    wz = sqrt(W) .* z;

    beta_new = (WX' * WX + ridgeMat) \ (WX' * wz);

    dev = -2 * sum(resp_allT .* log(mu) + (1 - resp_allT) .* log(1 - mu));

    if abs(dev_prev - dev) < opts.tol
        beta = beta_new;
        break;
    end

    beta = beta_new;
    dev_prev = dev;
end % iIter

% ---------------- final predictions on training data ----------------
predTrain = compute_predictions_from_beta(Z, resp_allT, beta, muZ, sdZ, opts.link, opts.zscorePredictor);

% ---------------- reconstruct template ----------------
beta0 = beta(1);
theta_vec = beta(2:end);
theta_vec_unscaled = theta_vec ./ sdZ(:); % unscale by predictor SD to get back to original basis space
Theta = reshape(theta_vec_unscaled, [Kori, Ksf]);
template2D = Bori * Theta * Bsf';

% ---------------- fit summary ----------------
idx1 = resp_allT == 1;
idx0 = resp_allT == 0;
pseudoR2_Tjur = mean(predTrain.yhat(idx1)) - mean(predTrain.yhat(idx0));

% ---------------- output ----------------
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

%% Helper: INTERNAL PREDICTION HELPER
function pred = compute_predictions_from_beta(Z, resp, beta, muZ, sdZ, linkName, zscorePredictor)

if zscorePredictor
    Zfit = (Z - muZ) ./ sdZ;
else
    Zfit = Z;
end

X = [ones(size(Zfit,1),1), Zfit];
eta = X * beta;
eta = max(min(eta, 8), -8);

switch lower(linkName)
    case 'logit'
        yhat = 1 ./ (1 + exp(-eta));
    case 'probit'
        yhat = normcdf(eta);
    otherwise
        error('Unknown link: use ''probit'' or ''logit''');
end

yhat = min(max(yhat, 1e-8), 1 - 1e-8);

nLL = -sum(resp .* log(yhat) + (1 - resp) .* log(1 - yhat));
deviance = -2 * sum(resp .* log(yhat) + (1 - resp) .* log(1 - yhat));

pred = struct();
pred.yhat = yhat;
pred.nLL = nLL;
pred.deviance = deviance;
pred.eta = eta;
end

%% Helper: MAKE CV FOLDS
function foldID = make_cv_folds(resp_allT, nFolds, opts)

nTrials = numel(resp_allT);
foldID = nan(nTrials,1);

if isfield(opts, 'rngSeed') && ~isempty(opts.rngSeed)
    rng(opts.rngSeed);
end

if opts.stratifyByResp
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
else
    idx = randperm(nTrials);
    labels = repmat(1:nFolds, 1, ceil(nTrials/nFolds));
    labels = labels(1:nTrials);
    foldID(idx) = labels;
end

foldID = foldID(:);
end

%% Helper: DEFAULT OPTS
function opts = fill_default_opts(opts)

if nargin < 1 || isempty(opts)
    opts = struct();
end

if ~isfield(opts, 'link') || isempty(opts.link), opts.link = 'probit'; end
if ~isfield(opts, 'nBasisORI') || isempty(opts.nBasisORI), opts.nBasisORI = 6; end
if ~isfield(opts, 'nBasisSF') || isempty(opts.nBasisSF), opts.nBasisSF = 5; end
if ~isfield(opts, 'sigmaORI_deg'), opts.sigmaORI_deg = []; end
if ~isfield(opts, 'sigmaSF_log2'), opts.sigmaSF_log2 = []; end
if ~isfield(opts, 'oriPeriod_deg') || isempty(opts.oriPeriod_deg), opts.oriPeriod_deg = 180; end
if ~isfield(opts, 'zscorePredictor') || isempty(opts.zscorePredictor), opts.zscorePredictor = true; end
if ~isfield(opts, 'ridge') || isempty(opts.ridge), opts.ridge = 10; end
if ~isfield(opts, 'maxIter') || isempty(opts.maxIter), opts.maxIter = 100; end
if ~isfield(opts, 'tol') || isempty(opts.tol), opts.tol = 1e-6; end
if ~isfield(opts, 'nFolds') || isempty(opts.nFolds), opts.nFolds = 5; end
if ~isfield(opts, 'foldID'), opts.foldID = []; end
if ~isfield(opts, 'rngSeed'), opts.rngSeed = []; end
if ~isfield(opts, 'stratifyByResp') || isempty(opts.stratifyByResp), opts.stratifyByResp = true; end
if ~isfield(opts, 'basisFamilyORI') || isempty(opts.basisFamilyORI), opts.basisFamilyORI = 'circ_gaussian'; end
if ~isfield(opts, 'basisFamilySF') || isempty(opts.basisFamilySF), opts.basisFamilySF = 'gaussianLog2'; end
if ~isfield(opts, 'kappaORI'), opts.kappaORI = []; end
if ~isfield(opts, 'basisWidthORI') || isempty(opts.basisWidthORI), opts.basisWidthORI = opts.basisWidth; end
if ~isfield(opts, 'basisWidthSF') || isempty(opts.basisWidthSF), opts.basisWidthSF = opts.basisWidth; end
if ~isfield(opts, 'asymSF_rightLeftRatio') || isempty(opts.asymSF_rightLeftRatio), opts.asymSF_rightLeftRatio = 1.5; end
if ~isfield(opts, 'widthSF_logParabola'), opts.widthSF_logParabola = []; end
if ~isfield(opts, 'widthSF_logParabola_left'),  opts.widthSF_logParabola_left  = []; end
if ~isfield(opts, 'widthSF_logParabola_right'), opts.widthSF_logParabola_right = []; end
end
