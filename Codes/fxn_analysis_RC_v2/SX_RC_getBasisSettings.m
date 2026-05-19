function cfg = SX_RC_getBasisSettings(overrides)
% SX_RC_getBasisSettings
% Centralized defaults for basis projection + smooth-basis CV fitting.
%
% Optional input:
%   overrides (struct) - fields to override any defaults below.

if nargin < 1 || isempty(overrides)
    overrides = struct();
end

cfg = struct();

% Shared basis family geometry
cfg.basisFamilyORI = 'vonmises';
cfg.basisFamilySF = 'asymGaussianLog2';
cfg.asymSF_rightLeftRatio = 1.2;
cfg.oriPeriod_deg = 180;

% Default basis counts/widths
cfg.nBasisORI = 4;
cfg.nBasisSF = 6;
cfg.basisWidthScaleORI = 0.9;
cfg.basisWidthScaleSF = 0.6;

% CV / ridge defaults for multivariate smooth regression
cfg.candidateORI = cfg.nBasisORI;
cfg.candidateSF = cfg.nBasisSF;
cfg.candidateRidge = 100;
cfg.nFolds = 5;
cfg.link = 'probit';
cfg.sigmaORI_deg = [];
cfg.sigmaSF_log2 = [];
cfg.zscorePredictor = true;
cfg.maxIter = 100;
cfg.tol = 1e-6;

% Apply caller overrides
fns = fieldnames(overrides);
for iF = 1:numel(fns)
    cfg.(fns{iF}) = overrides.(fns{iF});
end

% Keep candidate defaults in sync with basis counts unless explicitly overridden.
if ~isfield(overrides, 'candidateORI') && isfield(overrides, 'nBasisORI')
    cfg.candidateORI = cfg.nBasisORI;
end
if ~isfield(overrides, 'candidateSF') && isfield(overrides, 'nBasisSF')
    cfg.candidateSF = cfg.nBasisSF;
end
end