function output = SX_RC_selectBasis_cv(e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, ...
    candidateORI, candidateSF, candidateBasisFamilyORI, candidateBasisFamilySF, candidateRidge, opts)
% ============================================================
% Select ORI / SF basis counts, basis families, and ridge by
% cross-validated held-out negative log likelihood.
% Also conduct basis reconstruction (reversing the z-scoring of the basis predictors)
% Cached-Z + precomputed fold-standardization version.
% ============================================================

if nargin < 10
    error(['Need e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, ', ...
        'candidateORI, candidateSF, candidateBasisFamilyORI, candidateBasisFamilySF, candidateRidge, opts']);
end

opts = fill_default_opts(opts);

resp_allT = resp_allT(:);
nTrials = size(e3D_allT, 1);

if length(resp_allT) ~= nTrials
    error('resp_allT length does not match number of trials');
end

if ischar(candidateBasisFamilyORI) || isstring(candidateBasisFamilyORI)
    candidateBasisFamilyORI = cellstr(candidateBasisFamilyORI);
end
if ischar(candidateBasisFamilySF) || isstring(candidateBasisFamilySF)
    candidateBasisFamilySF = cellstr(candidateBasisFamilySF);
end

% candidate basis-width scales (allow scalar or vector); keep backward compatibility
candidateBasisWidthScaleORI = get_candidate_basis_width(opts, 'basisWidthScaleORI');
candidateBasisWidthScaleSF  = get_candidate_basis_width(opts, 'basisWidthScaleSF');
candidateAsymSF_rightLeftRatio = get_candidate_asym_sf_ratio(opts);

% ---------------- create CV folds ----------------
if isfield(opts, 'foldID') && ~isempty(opts.foldID)
    foldID = opts.foldID(:);
    if length(foldID) ~= nTrials
        error('opts.foldID length does not match number of trials');
    end
    nFolds = max(foldID);
else
    nFolds = opts.nFolds;
    foldID = make_cv_folds(resp_allT, nFolds, opts);
end

% precompute fold masks once
idxTrain_all = cell(nFolds,1);
idxTest_all  = cell(nFolds,1);
for iFold = 1:nFolds
    idxTest_all{iFold} = (foldID == iFold);
    idxTrain_all{iFold} = ~idxTest_all{iFold};
end

% ---------------- build all non-ridge combinations ----------------
comb = struct([]);
iComb = 0;

for iFamORI = 1:numel(candidateBasisFamilyORI)
    for iFamSF = 1:numel(candidateBasisFamilySF)
        for iORIcand = 1:numel(candidateORI)
            for iSFcand = 1:numel(candidateSF)
                for iBWori = 1:numel(candidateBasisWidthScaleORI)
                    for iBWsf = 1:numel(candidateBasisWidthScaleSF)
                        for iAsym = 1:numel(candidateAsymSF_rightLeftRatio)
                            iComb = iComb + 1;
                            comb(iComb).nBasisORI = candidateORI(iORIcand);
                            comb(iComb).nBasisSF  = candidateSF(iSFcand);
                            comb(iComb).basisFamilyORI = char(candidateBasisFamilyORI{iFamORI});
                            comb(iComb).basisFamilySF  = char(candidateBasisFamilySF{iFamSF});
                            comb(iComb).basisWidthScaleORI = candidateBasisWidthScaleORI(iBWori);
                            comb(iComb).basisWidthScaleSF  = candidateBasisWidthScaleSF(iBWsf);
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
    'basisWidthScaleORI', [], ...
    'basisWidthScaleSF', [], ...
    'asymSF_rightLeftRatio', [], ...
    'ridge', [], ...
    'nLL_mean', [], ...
    'nLL_se', [] );

iRes = 0;

% ---------------- main search ----------------
for iComb = 1:nCombBase

    opts_basis = opts;
    opts_basis.nBasisORI = comb(iComb).nBasisORI;
    opts_basis.nBasisSF  = comb(iComb).nBasisSF;
    opts_basis.basisFamilyORI = comb(iComb).basisFamilyORI;
    opts_basis.basisFamilySF  = comb(iComb).basisFamilySF;
    opts_basis.basisWidthScaleORI = comb(iComb).basisWidthScaleORI;
    opts_basis.basisWidthScaleSF  = comb(iComb).basisWidthScaleSF;
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
            beta = fit_ridge_glm_from_X(Xtrain_all{iFold}, rtrain_all{iFold}, ridge, opts.link, opts.maxIter, opts.tol);

            pred = predict_ridge_glm_from_X(Xtest_all{iFold}, rtest_all{iFold}, beta, opts.link);

            nLL_test_allFolds(iFold) = pred.nLL;
        end

        iRes = iRes + 1;
        res(iRes).nBasisORI = comb(iComb).nBasisORI;
        res(iRes).nBasisSF  = comb(iComb).nBasisSF;
        res(iRes).basisFamilyORI = comb(iComb).basisFamilyORI;
        res(iRes).basisFamilySF  = comb(iComb).basisFamilySF;
        res(iRes).basisWidthScaleORI = comb(iComb).basisWidthScaleORI;
        res(iRes).basisWidthScaleSF  = comb(iComb).basisWidthScaleSF;
        res(iRes).asymSF_rightLeftRatio = comb(iComb).asymSF_rightLeftRatio;
        res(iRes).ridge = ridge;
        res(iRes).nLL_mean = mean(nLL_test_allFolds, 'omitnan');
        res(iRes).nLL_se   = std(nLL_test_allFolds, 'omitnan') / sqrt(sum(isfinite(nLL_test_allFolds)));
    end % iRidge
end % iComb

% ---------------- choose best by lowest mean nLL ----------------
[~, idxBest] = min([res.nLL_mean]);
bestRes = res(idxBest);

% ---------------- refit on all data using full model ----------------
opts_best = opts;
opts_best.nBasisORI = bestRes.nBasisORI;
opts_best.nBasisSF  = bestRes.nBasisSF;
opts_best.basisFamilyORI = bestRes.basisFamilyORI;
opts_best.basisFamilySF  = bestRes.basisFamilySF;
opts_best.basisWidthScaleORI = bestRes.basisWidthScaleORI;
opts_best.basisWidthScaleSF  = bestRes.basisWidthScaleSF;
opts_best.asymSF_rightLeftRatio = bestRes.asymSF_rightLeftRatio;
opts_best.ridge = bestRes.ridge;

out_best = SX_RC_fit_smoothBasis(e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, opts_best);

output = [];
output.template2D = out_best.template2D;
output.nBasisORI = opts_best.nBasisORI;
output.nBasisSF = opts_best.nBasisSF;
output.basisFamilyORI = opts_best.basisFamilyORI;
output.basisFamilySF = opts_best.basisFamilySF;
output.basisWidthScaleORI = opts_best.basisWidthScaleORI;
output.basisWidthScaleSF = opts_best.basisWidthScaleSF;
output.asymSF_rightLeftRatio = opts_best.asymSF_rightLeftRatio;
output.ridge = opts_best.ridge;
output.bestRes = bestRes;
% output.res = res; % enable only if needed for diagnostics

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


%% Helper: FIT FUNCTION
% function fitOut = fit_ridge_glm_from_Z(Z, resp_allT, ridge, linkName, zscorePredictor, maxIter, tol)
%
% resp_allT = resp_allT(:);
% nTrials = size(Z, 1);
%
% % standardize predictors
% muZ = mean(Z, 1);
% sdZ = std(Z, 0, 1);
% sdZ(sdZ < 1e-8) = 1;
%
% if zscorePredictor
%     Zfit = (Z - muZ) ./ sdZ;
% else
%     Zfit = Z;
%     muZ = zeros(1, size(Z,2));
%     sdZ = ones(1, size(Z,2));
% end
%
% X = [ones(nTrials,1), Zfit];
% p = size(X,2);
%
% beta = zeros(p,1);
% ridgeMat = diag([0; ones(p-1,1)]) * ridge;
% dev_prev = Inf;
%
% for iIter = 1:maxIter
%     eta = X * beta;
%     eta = max(min(eta, 8), -8);
%
%     switch lower(linkName)
%         case 'logit'
%             mu = 1 ./ (1 + exp(-eta));
%             gprime = mu .* (1 - mu);
%         case 'probit'
%             mu = normcdf(eta);
%             gprime = normpdf(eta);
%         otherwise
%             error('Unknown link: use ''probit'' or ''logit''');
%     end
%
%     mu = min(max(mu, 1e-8), 1 - 1e-8);
%     gprime = max(gprime, 1e-6);
%
%     W = (gprime.^2) ./ (mu .* (1 - mu) + eps);
%     W = max(W, 1e-8);
%
%     z = eta + (resp_allT - mu) ./ (gprime + eps);
%
%     WX = X .* sqrt(W);
%     wz = sqrt(W) .* z;
%
%     beta_new = (WX' * WX + ridgeMat) \ (WX' * wz);
%
%     dev = -2 * sum(resp_allT .* log(mu) + (1 - resp_allT) .* log(1 - mu));
%
%     if abs(dev_prev - dev) < tol
%         beta = beta_new;
%         break;
%     end
%
%     beta = beta_new;
%     dev_prev = dev;
% end
%
% fitOut = struct();
% fitOut.beta = beta;
% fitOut.muZ = muZ;
% fitOut.sdZ = sdZ;
% end

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
if ~isfield(opts, 'basisWidthScale') || isempty(opts.basisWidthScale), opts.basisWidthScale = 0.8; end
if ~isfield(opts, 'basisWidthScaleORI') || isempty(opts.basisWidthScaleORI), opts.basisWidthScaleORI = opts.basisWidthScale; end
if ~isfield(opts, 'basisWidthScaleSF') || isempty(opts.basisWidthScaleSF), opts.basisWidthScaleSF = opts.basisWidthScale; end
if ~isfield(opts, 'asymSF_rightLeftRatio') || isempty(opts.asymSF_rightLeftRatio), opts.asymSF_rightLeftRatio = 1.5; end
if ~isfield(opts, 'widthSF_logParabola'), opts.widthSF_logParabola = []; end
if ~isfield(opts, 'widthSF_logParabola_left'),  opts.widthSF_logParabola_left  = []; end
if ~isfield(opts, 'widthSF_logParabola_right'), opts.widthSF_logParabola_right = []; end
end


%% Helper: normalize candidate basis-width option
function vals = get_candidate_basis_width(opts, fieldName)

if isfield(opts, fieldName) && ~isempty(opts.(fieldName))
    vals = opts.(fieldName);
elseif isfield(opts, 'basisWidthScale') && ~isempty(opts.basisWidthScale)
    vals = opts.basisWidthScale;
else
    vals = 0.8;
end

vals = vals(:)';
end

%% Helper: normalize candidate asymmetry ratio option
function vals = get_candidate_asym_sf_ratio(opts)

if isfield(opts, 'asymSF_rightLeftRatio') && ~isempty(opts.asymSF_rightLeftRatio)
    vals = opts.asymSF_rightLeftRatio;
else
    vals = 1.5;
end

vals = vals(:)';
end

%% Check helper: plot nLL
% function bestPerFam = plot_basis_cv_heatmap(res)
%
% candidateORI = unique([res.nBasisORI]);
% candidateSF  = unique([res.nBasisSF]);
%
% famORI = unique(string({res.basisFamilyORI}));
% famSF  = unique(string({res.basisFamilySF}));
%
% nRows = numel(famORI);
% nCols = numel(famSF);
%
% figure('Position', [100 100 800*nCols 600*nRows]);
% tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
%
% allVals = [res.nLL_mean];
% cMin = min(allVals);
% cMax = max(allVals);
% bestPerFam = nan(nRows, nCols, 2);
%
% for iFamORI = 1:nRows
%     for iFamSF = 1:nCols
%
%         nexttile; hold on;
%
%         thisFamORI = famORI(iFamORI);
%         thisFamSF  = famSF(iFamSF);
%
%         M = nan(numel(candidateORI), numel(candidateSF));
%
%         idxPanel = strcmp(string({res.basisFamilyORI}), thisFamORI) & ...
%             strcmp(string({res.basisFamilySF}),  thisFamSF);
%
%         resPanel = res(idxPanel);
%
%         for i = 1:numel(resPanel)
%             r = resPanel(i);
%
%             iRow = find(candidateORI == r.nBasisORI, 1, 'first');
%             iCol = find(candidateSF  == r.nBasisSF,  1, 'first');
%
%             M(iRow, iCol) = r.nLL_mean;
%         end
%
%         % find minimum nLL in this panel
%         [minVal, idxMin] = min(M(:), [], 'omitnan');
%         [iRowMin, iColMin] = ind2sub(size(M), idxMin);
%         best_nORI = candidateORI(iRowMin);
%         best_nSF  = candidateSF(iColMin);
%         bestPerFam(iFamORI, iFamSF, 1:2) = [best_nORI, best_nSF];
%
%         % Plot
%         imagesc(1:numel(candidateSF), 1:numel(candidateORI), M);
%         set(gca, 'YDir', 'normal', ...
%             'XTick', 1:numel(candidateSF), ...
%             'XTickLabel', string(candidateSF), ...
%             'YTick', 1:numel(candidateORI), ...
%             'YTickLabel', string(candidateORI), ...
%             'FontSize', 10);
%         axis square;
%         clim([cMin cMax]);
%         xlabel('nBasisSF');
%         ylabel('nBasisORI');
%         title(sprintf(['ORI: %s | SF: %s\n' 'min nLL = %.2f | nORI = %d, nSF = %d'], thisFamORI, thisFamSF, minVal, best_nORI, best_nSF), 'Interpreter', 'none');
%
%         % annotate cell values
%         for iRow = 1:size(M,1)
%             for iCol = 1:size(M,2)
%                 if isfinite(M(iRow, iCol))
%                     text(iCol, iRow, sprintf('%.2f', M(iRow, iCol)), ...
%                         'HorizontalAlignment', 'center', ...
%                         'VerticalAlignment', 'middle', ...
%                         'FontSize', 8, 'Color', 'k');
%                 end
%             end % iCol
%         end % iRow
%
%     end % iSF
% end % iORI
%
% cb = colorbar;
% cb.Layout.Tile = 'east';
% cb.Label.String = 'Mean held-out nLL';
%
% sgtitle('Cross-validated basis selection landscape');
% end

%%
% function plot_best_template_marginalized(template_est_2D, axis_ori_deg, axis_sf_log2, template_true, bestRes)
%
% nORI = numel(axis_ori_deg);
% nSF  = numel(axis_sf_log2);
%
% hasTruth = ~isempty(template_true);
%
% % ---------------- estimated template ----------------
% margORI_est = mean(template_est_2D, 2)';
% margSF_est  = mean(template_est_2D, 1);
%
% if max(abs(margORI_est)) > 0
%     margORI_est = margORI_est ./ max(abs(margORI_est));
% end
% if max(abs(margSF_est)) > 0
%     margSF_est = margSF_est ./ max(abs(margSF_est));
% end
%
% % ---------------- true template + true projected into selected basis ----------------
% if hasTruth
%     template_true_2D = reshape(template_true, [nORI, nSF]);
%
%     margORI_true = mean(template_true_2D, 2)';
%     margSF_true  = mean(template_true_2D, 1);
%
%     if max(abs(margORI_true)) > 0
%         margORI_true = margORI_true ./ max(abs(margORI_true));
%     end
%     if max(abs(margSF_true)) > 0
%         margSF_true = margSF_true ./ max(abs(margSF_true));
%     end
%
%     % Build opts from bestRes for basis-only projection
%     opts_proj = struct();
%     opts_proj.nBasisORI = bestRes.nBasisORI;
%     opts_proj.nBasisSF  = bestRes.nBasisSF;
%     opts_proj.basisFamilyORI = bestRes.basisFamilyORI;
%     opts_proj.basisFamilySF  = bestRes.basisFamilySF;
%     opts_proj.oriPeriod_deg = 180;
%     opts_proj.sigmaORI_deg = [];
%     opts_proj.sigmaSF_log2 = [];
%     opts_proj.kappaORI = [];
%     opts_proj.widthSF_logParabola = [];
%
%     out_proj = project_trueTemplate_to_basis(template_true_2D, axis_ori_deg, axis_sf_log2, opts_proj);
%     template_proj_2D = out_proj.template_proj_2D;
%
%     margORI_proj = mean(template_proj_2D, 2)';
%     margSF_proj  = mean(template_proj_2D, 1);
%
%     if max(abs(margORI_proj)) > 0
%         margORI_proj = margORI_proj ./ max(abs(margORI_proj));
%     end
%     if max(abs(margSF_proj)) > 0
%         margSF_proj = margSF_proj ./ max(abs(margSF_proj));
%     end
% end
%
% % ---------------- plot ----------------
% figure('Position', [100 100 1000 400]);
% tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
%
% nexttile; hold on;
% if hasTruth
%     % True template in original space
%     plot(axis_ori_deg, margORI_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
%     % True template projected into selected basis function
%     plot(axis_ori_deg, margORI_proj, 'b--', 'LineWidth', 1.8, 'DisplayName', 'True \rightarrow basis');
% end
% % Estimated template, in basis function
% plot(axis_ori_deg, margORI_est, 'r--', 'LineWidth', 1.8, 'DisplayName', 'Est.');
%
% legend('show', 'Location', 'best');
% box on;
% xlabel('Orientation (deg)');
% ylabel('Marginalized weight');
% ylim([-0.2 1.05]);
% title(sprintf('ORI | %s | %s | nORI=%d, nSF=%d', ...
%     bestRes.basisFamilyORI, bestRes.basisFamilySF, bestRes.nBasisORI, bestRes.nBasisSF), ...
%     'Interpreter', 'none');
%
% nexttile; hold on;
% if hasTruth
%     plot(2.^axis_sf_log2, margSF_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
%     plot(2.^axis_sf_log2, margSF_proj, 'b--', 'LineWidth', 1.8, 'DisplayName', 'True \rightarrow basis');
% end
% plot(2.^axis_sf_log2, margSF_est, 'r--', 'LineWidth', 1.8, 'DisplayName', 'Est.');
% legend('show', 'Location', 'best');
% box on;
% set(gca, 'XScale', 'log');
% xlabel('Spatial frequency (cpd)');
% ylabel('Marginalized weight');
% ylim([-0.2 1.05]);
% title(sprintf('SF | %s | %s | nORI=%d, nSF=%d', ...
%     bestRes.basisFamilyORI, bestRes.basisFamilySF, bestRes.nBasisORI, bestRes.nBasisSF), ...
%     'Interpreter', 'none');
%
% sgtitle('Best recovered template (marginalized)');
% end

%%
% function out = project_trueTemplate_to_basis(template_true_2D, axis_ori_deg, axis_sf_log2, opts)

% % Build basis
% Bori = predSFkernel('make_basis_ori', axis_ori_deg, opts, 0);
% Bsf  = predSFkernel('make_basis_sf', axis_sf_log2, opts, 0);

% Kori = size(Bori, 2);
% Ksf  = size(Bsf, 2);

% % Solve least-squares projection:
% % vec(T) ≈ kron(Bsf, Bori) * vec(Theta)
% A = kron(Bsf, Bori);                  % [nORI*nSF x Kori*Ksf]
% t = template_true_2D(:);              % [nORI*nSF x 1]

% theta_vec = A \ t;                    % least-squares solution
% Theta = reshape(theta_vec, [Kori, Ksf]);

% template_proj_2D = Bori * Theta * Bsf';

% out = struct();
% out.template_proj_2D = template_proj_2D;
% out.Theta = Theta;
% out.Bori = Bori;
% out.Bsf = Bsf;
% end
