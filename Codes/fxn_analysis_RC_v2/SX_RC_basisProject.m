function Z = SX_RC_basisProject(e3D_allT, axis_ori_deg, axis_sf_log2, basisOpts)
% SX_RC_basisProject
% Project ORI x SF energy maps into the tensor-product basis used by SX_RC_selectBasis_cv.

basisOpts = fill_default_opts_local(basisOpts);

Bori = make_basis_ori_local(axis_ori_deg, basisOpts);
Bsf = make_basis_sf_local(axis_sf_log2, basisOpts);

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

function opts = fill_default_opts_local(opts)
if nargin < 1 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'nBasisORI') || isempty(opts.nBasisORI), opts.nBasisORI = 6; end
if ~isfield(opts, 'nBasisSF') || isempty(opts.nBasisSF), opts.nBasisSF = 5; end
if ~isfield(opts, 'basisFamilyORI') || isempty(opts.basisFamilyORI), opts.basisFamilyORI = 'circ_gaussian'; end
if ~isfield(opts, 'basisFamilySF') || isempty(opts.basisFamilySF), opts.basisFamilySF = 'gaussianLog2'; end
if ~isfield(opts, 'basisWidthScaleORI') || isempty(opts.basisWidthScaleORI), opts.basisWidthScaleORI = 0.8; end
if ~isfield(opts, 'basisWidthScaleSF') || isempty(opts.basisWidthScaleSF), opts.basisWidthScaleSF = 0.8; end
if ~isfield(opts, 'asymSF_rightLeftRatio') || isempty(opts.asymSF_rightLeftRatio), opts.asymSF_rightLeftRatio = 1.5; end
if ~isfield(opts, 'sigmaORI_deg'), opts.sigmaORI_deg = []; end
if ~isfield(opts, 'sigmaSF_log2'), opts.sigmaSF_log2 = []; end
if ~isfield(opts, 'kappaORI'), opts.kappaORI = []; end
if ~isfield(opts, 'sigmaSF_log2_left'), opts.sigmaSF_log2_left = []; end
if ~isfield(opts, 'sigmaSF_log2_right'), opts.sigmaSF_log2_right = []; end
if ~isfield(opts, 'oriPeriod_deg') || isempty(opts.oriPeriod_deg), opts.oriPeriod_deg = 180; end
end

function Bori = make_basis_ori_local(axis_ori_deg, opts)
axis_ori_deg = axis_ori_deg(:);
centersORI_deg = linspace(0, opts.oriPeriod_deg, opts.nBasisORI + 1);
centersORI_deg(end) = [];

switch lower(opts.basisFamilyORI)
    case 'circ_gaussian'
        if isempty(opts.sigmaORI_deg)
            if opts.nBasisORI > 1
                d = opts.oriPeriod_deg / opts.nBasisORI;
                sigmaORI_deg = d * opts.basisWidthScaleORI;
            else
                sigmaORI_deg = opts.oriPeriod_deg / 4;
            end
        else
            sigmaORI_deg = opts.sigmaORI_deg;
        end
        Bori = circular_gaussian_basis_local(axis_ori_deg, centersORI_deg, sigmaORI_deg, opts.oriPeriod_deg);

    case 'vonmises'
        if isempty(opts.kappaORI)
            if opts.nBasisORI > 1
                d_deg = opts.oriPeriod_deg / opts.nBasisORI;
                sigma_rad = deg2rad(d_deg * opts.basisWidthScaleORI);
                kappaORI = 1 / max(sigma_rad.^2, 1e-6);
            else
                kappaORI = 2;
            end
        else
            kappaORI = opts.kappaORI;
        end
        Bori = vonmises_basis_deg_local(axis_ori_deg, centersORI_deg, kappaORI, opts.oriPeriod_deg);

    otherwise
        error('Unknown ORI basis family: %s', opts.basisFamilyORI);
end

Bori = normalize_columns_local(Bori);
end

function B = circular_gaussian_basis_local(x_deg, centers_deg, sigma_deg, period_deg)
x_deg = x_deg(:);
centers_deg = centers_deg(:)';
B = zeros(length(x_deg), length(centers_deg));
for k = 1:length(centers_deg)
    d = mod(x_deg - centers_deg(k) + period_deg/2, period_deg) - period_deg/2;
    B(:, k) = exp(-0.5 * (d ./ sigma_deg).^2);
end
end

function B = vonmises_basis_deg_local(x_deg, centers_deg, kappa, period_deg)
x_deg = x_deg(:);
centers_deg = centers_deg(:)';
B = zeros(length(x_deg), length(centers_deg));
for k = 1:length(centers_deg)
    d = mod(x_deg - centers_deg(k) + period_deg/2, period_deg) - period_deg/2;
    theta = d / period_deg * 2*pi;
    B(:,k) = exp(kappa * cos(theta));
end
end

function Bsf = make_basis_sf_local(axis_sf_log2, opts)
axis_sf_log2 = axis_sf_log2(:);
axis_sf_cpd = 2 .^ axis_sf_log2;

switch lower(opts.basisFamilySF)
    case 'gaussianlog2'
        centersSF_log2 = linspace(min(axis_sf_log2), max(axis_sf_log2), opts.nBasisSF);
        if isempty(opts.sigmaSF_log2)
            if opts.nBasisSF > 1
                d = mean(diff(centersSF_log2));
                sigmaSF_log2 = d * opts.basisWidthScaleSF;
            else
                sigmaSF_log2 = range(axis_sf_log2) / 4 + eps;
            end
        else
            sigmaSF_log2 = opts.sigmaSF_log2;
        end
        Bsf = gaussian_basis_local(axis_sf_log2, centersSF_log2, sigmaSF_log2);

    case 'asymgaussianlog2'
        centersSF_log2 = linspace(min(axis_sf_log2), max(axis_sf_log2), opts.nBasisSF);
        if isempty(opts.sigmaSF_log2_left)
            if opts.nBasisSF > 1
                d = mean(diff(centersSF_log2));
                sigma_left = d * opts.basisWidthScaleSF;
            else
                sigma_left = range(axis_sf_log2) / 4 + eps;
            end
        else
            sigma_left = opts.sigmaSF_log2_left;
        end
        if isempty(opts.sigmaSF_log2_right)
            sigma_right = sigma_left * opts.asymSF_rightLeftRatio;
        else
            sigma_right = opts.sigmaSF_log2_right;
        end
        Bsf = zeros(length(axis_sf_log2), numel(centersSF_log2));
        for k = 1:numel(centersSF_log2)
            Bsf(:,k) = asym_gaussian_log2_basis_local(axis_sf_log2, centersSF_log2(k), sigma_left, sigma_right);
        end

    case 'logparabola'
        centersSF_cpd = logspace(log10(min(axis_sf_cpd)), log10(max(axis_sf_cpd)), opts.nBasisSF);
        if opts.nBasisSF > 1
            d = mean(diff(log10(centersSF_cpd)));
            widthSF = max(d * 1.2, 1e-3);
        else
            widthSF = 0.3;
        end
        Bsf = zeros(length(axis_sf_cpd), numel(centersSF_cpd));
        for k = 1:numel(centersSF_cpd)
            Bsf(:,k) = logparabola_basis_local(axis_sf_cpd, centersSF_cpd(k), widthSF);
        end

    case 'asymlogparabola'
        centersSF_cpd = logspace(log10(min(axis_sf_cpd)), log10(max(axis_sf_cpd)), opts.nBasisSF);
        if opts.nBasisSF > 1
            d = mean(diff(log10(centersSF_cpd)));
            width_left = max(d * 1.0, 1e-3);
        else
            width_left = 0.3;
        end
        width_right = width_left * opts.asymSF_rightLeftRatio;
        Bsf = zeros(length(axis_sf_cpd), numel(centersSF_cpd));
        for k = 1:numel(centersSF_cpd)
            Bsf(:,k) = asym_logparabola_basis_local(axis_sf_cpd, centersSF_cpd(k), width_left, width_right);
        end

    otherwise
        error('Unknown SF basis family: %s', opts.basisFamilySF);
end

Bsf = normalize_columns_local(Bsf);
end

function B = gaussian_basis_local(x, centers, sigma)
x = x(:);
centers = centers(:)';
B = exp(-0.5 * ((x - centers) ./ sigma).^2);
end

function y = asym_gaussian_log2_basis_local(x_log2, center_log2, sigma_left, sigma_right)
x_log2 = x_log2(:);
y = zeros(size(x_log2));
idxLeft = x_log2 <= center_log2;
idxRight = x_log2 > center_log2;
y(idxLeft) = exp(-0.5 * ((x_log2(idxLeft) - center_log2) ./ sigma_left).^2);
y(idxRight) = exp(-0.5 * ((x_log2(idxRight) - center_log2) ./ sigma_right).^2);
end

function y = logparabola_basis_local(x_cpd, peakSF, width)
x_cpd = max(x_cpd(:), eps);
y = 10 .^ (-(log10(x_cpd ./ peakSF) ./ width).^2);
end

function y = asym_logparabola_basis_local(x_cpd, peakSF, width_left, width_right)
x_cpd = max(x_cpd(:), eps);
y = zeros(size(x_cpd));
idxLeft = x_cpd <= peakSF;
idxRight = x_cpd > peakSF;
y(idxLeft) = 10 .^ (-(log10(x_cpd(idxLeft) ./ peakSF) ./ width_left).^2);
y(idxRight) = 10 .^ (-(log10(x_cpd(idxRight) ./ peakSF) ./ width_right).^2);
end

function B = normalize_columns_local(B)
nrm = sqrt(sum(B.^2, 1));
nrm(nrm == 0) = 1;
B = bsxfun(@rdivide, B, nrm);
end
