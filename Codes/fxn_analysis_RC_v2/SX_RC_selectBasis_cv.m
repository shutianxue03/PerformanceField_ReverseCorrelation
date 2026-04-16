function output = SX_RC_selectBasis_cv(e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, candidateORI, candidateSF, opts)
% ============================================================
% Select number of ORI / SF basis functions by cross-validated
% held-out negative log likelihood for step 1 only.
% Created by Shutian Xue on 04/16/2026
% ============================================================

if nargin < 6
    error('Need e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, candidateORI, candidateSF, opts');
end

opts = fill_default_opts(opts);

resp_allT = resp_allT(:);
nTrials = size(e3D_allT, 1);

if length(resp_allT) ~= nTrials
    error('resp_allT length does not match number of trials');
end

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

% ---------------- outer grid search ----------------
nORIcand = numel(candidateORI);
nSFcand  = numel(candidateSF);

res = struct([]);
iRes = 0;

% fprintf('\n============================================================\n');
% fprintf('%s\n', string(datetime('now')));
% fprintf('Selecting basis counts by %d-fold CV\n', nFolds);
% fprintf('Candidate ORI: %s\n', mat2str(candidateORI));
% fprintf('Candidate SF : %s\n', mat2str(candidateSF));
% fprintf('============================================================\n');

for iORIcand = 1:nORIcand
    for iSFcand = 1:nSFcand

        nBasisORI = candidateORI(iORIcand);
        nBasisSF  = candidateSF(iSFcand);

        % fprintf('\nTesting nBasisORI = %d, nBasisSF = %d\n', nBasisORI, nBasisSF);

        opts_this = opts;
        opts_this.nBasisORI = nBasisORI;
        opts_this.nBasisSF  = nBasisSF;

        nLL_test_allFolds = nan(nFolds,1);
        dev_test_allFolds = nan(nFolds,1);
        acc_test_allFolds = nan(nFolds,1);

        for iFold = 1:nFolds
            idxTrain = foldID ~= iFold;
            idxTest  = foldID == iFold;

            e3D_train = e3D_allT(idxTrain, :, :);
            r_train   = resp_allT(idxTrain);

            e3D_test = e3D_allT(idxTest, :, :);
            r_test   = resp_allT(idxTest);

            %-----------------------%
            % Fit to training set with the selected basis functions combo
            fitOut = SX_RC_fit_smoothBasis_circORI_logSF(e3D_train, r_train, axis_ori_deg, axis_sf_log2, opts_this);
            %-----------------------%

            %-----------------------%
            % Predict for test set
            predOut = SX_RC_predict_smoothBasis_circORI_logSF(e3D_test, r_test, fitOut);
            %-----------------------%

            nLL_test_allFolds(iFold) = predOut.nLL;
            dev_test_allFolds(iFold) = predOut.deviance;
            acc_test_allFolds(iFold) = mean((predOut.yhat >= 0.5) == r_test);

            % fprintf('  Fold %d/%d: nLL = %.3f, dev = %.3f, acc = %.3f\n', iFold, nFolds, predOut.nLL, predOut.deviance, acc_test_allFolds(iFold));
        end % iFold

        iRes = iRes + 1;
        res(iRes).nBasisORI = nBasisORI;
        res(iRes).nBasisSF  = nBasisSF;
        res(iRes).nLL_test_allFolds = nLL_test_allFolds;
        res(iRes).dev_test_allFolds = dev_test_allFolds;
        res(iRes).acc_test_allFolds = acc_test_allFolds;

        res(iRes).nLL_mean = mean(nLL_test_allFolds, 'omitnan');
        res(iRes).nLL_se   = std(nLL_test_allFolds, 'omitnan') / sqrt(sum(isfinite(nLL_test_allFolds)));

        res(iRes).dev_mean = mean(dev_test_allFolds, 'omitnan');
        res(iRes).dev_se   = std(dev_test_allFolds, 'omitnan') / sqrt(sum(isfinite(dev_test_allFolds)));

        res(iRes).acc_mean = mean(acc_test_allFolds, 'omitnan');
        res(iRes).acc_se   = std(acc_test_allFolds, 'omitnan') / sqrt(sum(isfinite(acc_test_allFolds)));
    end % iSFcand
end % iORIcand

% ---------------- choose best by lowest mean nLL ----------------
nLL_all = [res.nLL_mean];
[~, idxBest] = min(nLL_all);

results = struct();
results.all = res;
results.best = res(idxBest);

results.table = struct2table(rmfield(res, {'nLL_test_allFolds','dev_test_allFolds','acc_test_allFolds'}));
results.foldID = foldID;
results.opts = opts;

% fprintf('\n============================================================\n');
% fprintf('Best basis pair:\n');
% fprintf('  nBasisORI = %d\n', results.best.nBasisORI);
% fprintf('  nBasisSF  = %d\n', results.best.nBasisSF);
% fprintf('  mean held-out nLL = %.3f\n', results.best.nLL_mean);
% fprintf('============================================================\n');

%% Refit on all data using the best basis counts
opts_best = opts;
opts_best.nBasisORI = results.best.nBasisORI;
opts_best.nBasisSF  = results.best.nBasisSF;

out_best= SX_RC_fit_smoothBasis_circORI_logSF(e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, opts_best);
output = [];
output.template2D = out_best.template2D;
output.nBasisORI = opts_best.nBasisORI;
output.nBasisSF = opts_best.nBasisSF;

end

%% Local helper functions
% ============================================================
% FIT FUNCTION
% ============================================================
function out = SX_RC_fit_smoothBasis_circORI_logSF(e3D_allT, resp_allT, axis_ori_deg, axis_sf_log2, opts)

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

%% ---------------- ORI basis ----------------
centersORI_deg = linspace(0, opts.oriPeriod_deg, opts.nBasisORI + 1);
centersORI_deg(end) = [];

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

Bori = circular_gaussian_basis(axis_ori_deg, centersORI_deg, sigmaORI_deg, opts.oriPeriod_deg);
Bori = normalize_columns(Bori);

%% ---------------- SF basis ----------------
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

Bsf = gaussian_basis(axis_sf_log2, centersSF_log2, sigmaSF_log2);
Bsf = normalize_columns(Bsf);

centersSF_cpd = 2 .^ centersSF_log2;

Kori = size(Bori, 2);
Ksf = size(Bsf, 2);
Ktot = Kori * Ksf;

%% ---------------- project trials into basis space ----------------
Z = project_trials_to_basis(e3D_allT, Bori, Bsf);

%% ---------------- standardize predictors ----------------
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

%% ---------------- fit ridge-regularized binomial GLM via IRLS ----------------
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
end

%% ---------------- final predictions on training data ----------------
predTrain = compute_predictions_from_beta(Z, resp_allT, beta, muZ, sdZ, opts.link, opts.zscorePredictor);

%% ---------------- reconstruct template ----------------
beta0 = beta(1);
theta_vec = beta(2:end);
theta_vec_unscaled = theta_vec ./ sdZ(:);
Theta = reshape(theta_vec_unscaled, [Kori, Ksf]);
template2D = Bori * Theta * Bsf';

%% ---------------- fit summary ----------------
idx1 = resp_allT == 1;
idx0 = resp_allT == 0;
pseudoR2_Tjur = mean(predTrain.yhat(idx1)) - mean(predTrain.yhat(idx0));

%% ---------------- output ----------------
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

% ============================================================
% PREDICT FUNCTION
% ============================================================
function pred = SX_RC_predict_smoothBasis_circORI_logSF(e3D_test, resp_test, fitOut)

resp_test = resp_test(:);

Ztest = project_trials_to_basis(e3D_test, fitOut.Bori, fitOut.Bsf);

pred = compute_predictions_from_beta( ...
    Ztest, resp_test, fitOut.beta, fitOut.muZ, fitOut.sdZ, ...
    fitOut.link, fitOut.zscorePredictor);

pred.Ztest = Ztest;
end

% ============================================================
% INTERNAL PREDICTION HELPER
% ============================================================
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

% ============================================================
% PROJECT ENERGY MAPS INTO TENSOR-PRODUCT BASIS
% ============================================================
function Z = project_trials_to_basis(e3D_allT, Bori, Bsf)

nTrials = size(e3D_allT, 1);
Kori = size(Bori, 2);
Ksf  = size(Bsf, 2);
Ktot = Kori * Ksf;

Z = zeros(nTrials, Ktot);

for iTrial = 1:nTrials
    Ei = squeeze(e3D_allT(iTrial, :, :));
    Zi = Bori' * Ei * Bsf;
    Z(iTrial, :) = Zi(:)';
end
end

% ============================================================
% MAKE CV FOLDS
% ============================================================
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

% ============================================================
% DEFAULT OPTS
% ============================================================
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
if ~isfield(opts, 'ridge') || isempty(opts.ridge), opts.ridge = 1; end
if ~isfield(opts, 'maxIter') || isempty(opts.maxIter), opts.maxIter = 100; end
if ~isfield(opts, 'tol') || isempty(opts.tol), opts.tol = 1e-6; end
if ~isfield(opts, 'nFolds') || isempty(opts.nFolds), opts.nFolds = 5; end
if ~isfield(opts, 'foldID'), opts.foldID = []; end
if ~isfield(opts, 'rngSeed'), opts.rngSeed = []; end
if ~isfield(opts, 'stratifyByResp') || isempty(opts.stratifyByResp), opts.stratifyByResp = true; end
end

% ============================================================
% HELPERS
% ============================================================
function B = gaussian_basis(x, centers, sigma)
x = x(:);
centers = centers(:)';
B = exp(-0.5 * ((x - centers) ./ sigma).^2);
end

function B = circular_gaussian_basis(x_deg, centers_deg, sigma_deg, period_deg)
x_deg = x_deg(:);
centers_deg = centers_deg(:)';
B = zeros(length(x_deg), length(centers_deg));
for k = 1:length(centers_deg)
    d = mod(x_deg - centers_deg(k) + period_deg/2, period_deg) - period_deg/2;
    B(:,k) = exp(-0.5 * (d ./ sigma_deg).^2);
end
end

function B = normalize_columns(B)
nrm = sqrt(sum(B.^2, 1));
nrm(nrm == 0) = 1;
B = bsxfun(@rdivide, B, nrm);
end