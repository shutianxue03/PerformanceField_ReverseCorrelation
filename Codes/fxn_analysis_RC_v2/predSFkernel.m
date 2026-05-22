
function y = predSFkernel(x, iFamily, params, flag_plot)

% Utility mode: shared basis construction for RC projection code.
if ischar(x) || isstring(x)
    cmd = lower(string(x));
    switch cmd
        case "make_basis_ori"
            y = make_basis_ori_shared(iFamily, params);
        case "make_basis_sf"
            y = make_basis_sf_shared(iFamily, params);
        otherwise
            error('Unknown utility command for predSFkernel: %s', cmd);
    end
    return;
end

SX_normPDF = @(x,mu,sigma) exp(-((x-mu)/sigma).^2);

switch iFamily
    case 1 % gaussian
        gain = params(1);
        width_deg = params(2);
        base = params(3);
        y = gain * SX_normPDF(x, 0, width_deg) + base;

        % periodDeg = 180;
        % d = mod(x + periodDeg/2, periodDeg) - periodDeg/2;
        % y = gain * exp(-0.5 * (d ./ width_deg).^2) + base;

    case 2 % log parabola
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        base = params(4);
        y = gain * 10.^(-(log10(x / peakSF) / width).^2) + base;

    case 3 % truncated log parabola
        % x is linear!! peakSF and width are also linear!!
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        base = params(4);
        lambda  = params(5);
        y = gain * 10.^(-(log10(x / peakSF) / width).^2) + base;
        y((x < peakSF) & (y < lambda)) = lambda;

    case 4 % raised gaussian
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        base = params(4);
        power  = params(5);
        y = gain * exp(-1/2*(x-peakSF).^2/width^2) .^power + base;

    case 5 % double exponential (JigoCarrasco 2020 for neutral CSF, eq 2)
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        y = gain * x .^ (peakSF / width) .* exp(-x/width);

    case 6 % skewed gaussian
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        y = 2*SX_normPDF(x, peakSF, width).*SX_normPDF(gain*x, peakSF, width);

    case 7 % truncated raised gaussian
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        base = params(4);
        power  = params(5);
        lambda = params(6);
        y = gain * SX_normPDF(x,peakSF,width).^ power + base;
        y(x < peakSF & y < lambda) = lambda;

    case 8 % difference of gaussians
        % from jigo2023
        gain1 = params(1);
        gain2 = params(2);
        sigma1 = params(3);
        sigma_r = params(4);
        baseline = params(5);
        y = gain1*(SX_normPDF(x, 0, sigma1) - gain2*SX_normPDF(x, 0, sigma1*sigma_r)) + baseline;

    case 9 % gaussian to SF
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        base = params(4);
        y = gain * SX_normPDF(x, peakSF, width) + base;

    case 10 % von Mises for ORI (180-deg periodic)
        gain = params(1);
        k = params(2); % concentration parameter
        base = params(3);

        theta = x / 180 * 2 * pi; % match period=180 mapping used in basis construction
        y = gain * v_vonmisespdf(theta, 0, k) + base;

    case 11 % gaussian with mean free to vary
        peak = params(1);
        gain = params(2);
        width = params(3);
        base = params(4);
        y = gain * SX_normPDF(x, peak, width) + base;

    case 12 % double peak
        peakSF1 = params(1);
        gain1 = params(2);
        width1 = params(3);
        base1 = params(4);
        peakSF2 = params(5);
        gain2 = params(6);
        width2 = params(7);
        base2 = params(8);
        y1 = gain1 * 10.^(-(log10(x / peakSF1) / width1).^2) + base1;
        y2 = gain2 * 10.^(-(log10(x / peakSF2) / width2).^2) + base2;
        y = y1 + y2;

    case 13 % raised DoG
        gain1 = params(1);
        gain2 = params(2);
        sigma1 = params(3);
        sigma_r = params(4);
        power = params(5);
        baseline = params(6);
        y = gain1*(SX_normPDF(x, 0, sigma1).^power - gain2*SX_normPDF(x, 0, sigma1*sigma_r))+ baseline;

    case 14 % asymmetric Gaussian for SF
        peakSF = params(1);
        gain = params(2);
        width_left = params(3);
        width_right = params(4);
        base = params(5);
        y = zeros(size(x));
        idxLeft = x <= peakSF;
        idxRight = x > peakSF;
        y(idxLeft) = gain * SX_normPDF(x(idxLeft), peakSF, width_left) + base;
        y(idxRight) = gain * SX_normPDF(x(idxRight), peakSF, width_right) + base;
end

%% Plot
if flag_plot
    figure, hold on
    plot(x, y, 'o-')
    % ylimit = ylim;
    switch iFamily
        case {1,8,10,11,13}
            xline(0, 'k--');
            xlabel('Orientation (º)')
            xlim([-90, 90])
            xticks(-90:45:90)
        otherwise
            xline(1, 'k--');
            xlabel('SF (cpd)')
            xlim([1,4])
            xticks([1,2,4]), xticklabels([1,2,4])
    end
    ylabel('sensitivity kernel (a.u.)')

    text(3.5, 1, sprintf('f0 = %.2f\nsigma = %.2f\nalpha = %.2f\nb = %.2f\nlambda = %.2f', params))
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
end

end

function opts = fill_basis_opts_shared(opts)
if nargin < 1 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'nBasisORI') || isempty(opts.nBasisORI), opts.nBasisORI = 6; end
if ~isfield(opts, 'nBasisSF') || isempty(opts.nBasisSF), opts.nBasisSF = 5; end
if ~isfield(opts, 'basisFamilyORI') || isempty(opts.basisFamilyORI), opts.basisFamilyORI = 'circ_gaussian'; end
if ~isfield(opts, 'basisFamilySF') || isempty(opts.basisFamilySF), opts.basisFamilySF = 'gaussianLog2'; end
if ~isfield(opts, 'basisWidthScaleORI') || isempty(opts.basisWidthScaleORI), opts.basisWidthScaleORI = 0.8; end
if ~isfield(opts, 'basisWidthScaleSF') || isempty(opts.basisWidthScaleSF), opts.basisWidthScaleSF = 0.8 ; end
if ~isfield(opts, 'asymSF_rightLeftRatio') || isempty(opts.asymSF_rightLeftRatio), opts.asymSF_rightLeftRatio = 1.5; end
if ~isfield(opts, 'sigmaORI_deg'), opts.sigmaORI_deg = []; end
if ~isfield(opts, 'sigmaSF_log2'), opts.sigmaSF_log2 = []; end
if ~isfield(opts, 'kappaORI'), opts.kappaORI = []; end
if ~isfield(opts, 'sigmaSF_log2_left'), opts.sigmaSF_log2_left = []; end
if ~isfield(opts, 'sigmaSF_log2_right'), opts.sigmaSF_log2_right = []; end
if ~isfield(opts, 'oriPeriod_deg') || isempty(opts.oriPeriod_deg), opts.oriPeriod_deg = 180; end
if ~isfield(opts, 'widthSF_logParabola'), opts.widthSF_logParabola = []; end
if ~isfield(opts, 'widthSF_logParabola_left'), opts.widthSF_logParabola_left = []; end
if ~isfield(opts, 'widthSF_logParabola_right'), opts.widthSF_logParabola_right = []; end
end

function Bori = make_basis_ori_shared(axis_ori_deg, opts)
opts = fill_basis_opts_shared(opts);
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
        Bori = circular_gaussian_basis_shared(axis_ori_deg, centersORI_deg, sigmaORI_deg, opts.oriPeriod_deg);

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
        Bori = vonmises_basis_deg_shared(axis_ori_deg, centersORI_deg, kappaORI, opts.oriPeriod_deg);

    otherwise
        error('Unknown ORI basis family: %s', opts.basisFamilyORI);
end

Bori = normalize_columns_shared(Bori);
end

function Bsf = make_basis_sf_shared(axis_sf_log2, opts)
opts = fill_basis_opts_shared(opts);
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
        Bsf = gaussian_basis_shared(axis_sf_log2, centersSF_log2, sigmaSF_log2);

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
            Bsf(:,k) = asym_gaussian_log2_basis_shared(axis_sf_log2, centersSF_log2(k), sigma_left, sigma_right);
        end

    case 'logparabola'
        centersSF_cpd = logspace(log10(min(axis_sf_cpd)), log10(max(axis_sf_cpd)), opts.nBasisSF);
        if ~isfield(opts, 'widthSF_logParabola') || isempty(opts.widthSF_logParabola)
            if opts.nBasisSF > 1
                d = mean(diff(log10(centersSF_cpd)));
                widthSF = max(d * 1.2, 1e-3);
            else
                widthSF = 0.3;
            end
        else
            widthSF = opts.widthSF_logParabola;
        end
        Bsf = zeros(length(axis_sf_cpd), numel(centersSF_cpd));
        for k = 1:numel(centersSF_cpd)
            Bsf(:,k) = logparabola_basis_shared(axis_sf_cpd, centersSF_cpd(k), widthSF);
        end

    case 'asymlogparabola'
        centersSF_cpd = logspace(log10(min(axis_sf_cpd)), log10(max(axis_sf_cpd)), opts.nBasisSF);
        if ~isfield(opts, 'widthSF_logParabola_left') || isempty(opts.widthSF_logParabola_left)
            if opts.nBasisSF > 1
                d = mean(diff(log10(centersSF_cpd)));
                width_left = max(d * 1.0, 1e-3);
            else
                width_left = 0.3;
            end
        else
            width_left = opts.widthSF_logParabola_left;
        end
        if ~isfield(opts, 'widthSF_logParabola_right') || isempty(opts.widthSF_logParabola_right)
            width_right = width_left * opts.asymSF_rightLeftRatio;
        else
            width_right = opts.widthSF_logParabola_right;
        end
        Bsf = zeros(length(axis_sf_cpd), numel(centersSF_cpd));
        for k = 1:numel(centersSF_cpd)
            Bsf(:,k) = asym_logparabola_basis_shared(axis_sf_cpd, centersSF_cpd(k), width_left, width_right);
        end

    otherwise
        error('Unknown SF basis family: %s', opts.basisFamilySF);
end

Bsf = normalize_columns_shared(Bsf);
end

function B = circular_gaussian_basis_shared(x_deg, centers_deg, sigma_deg, period_deg)
x_deg = x_deg(:);
centers_deg = centers_deg(:)';
B = zeros(length(x_deg), length(centers_deg));
for k = 1:length(centers_deg)
    d = mod(x_deg - centers_deg(k) + period_deg/2, period_deg) - period_deg/2;
    B(:, k) = exp(-0.5 * (d ./ sigma_deg).^2);
end
end

function B = vonmises_basis_deg_shared(x_deg, centers_deg, kappa, period_deg)
x_deg = x_deg(:);
centers_deg = centers_deg(:)';
B = zeros(length(x_deg), length(centers_deg));
for k = 1:length(centers_deg)
    d = mod(x_deg - centers_deg(k) + period_deg/2, period_deg) - period_deg/2;
    theta = d / period_deg * 2*pi;
    B(:,k) = exp(kappa * cos(theta));
end
end

function B = gaussian_basis_shared(x, centers, sigma)
x = x(:);
centers = centers(:)';
B = exp(-0.5 * ((x - centers) ./ sigma).^2);
end

function y = asym_gaussian_log2_basis_shared(x_log2, center_log2, sigma_left, sigma_right)
x_log2 = x_log2(:);
y = zeros(size(x_log2));
idxLeft = x_log2 <= center_log2;
idxRight = x_log2 > center_log2;
y(idxLeft) = exp(-0.5 * ((x_log2(idxLeft) - center_log2) ./ sigma_left).^2);
y(idxRight) = exp(-0.5 * ((x_log2(idxRight) - center_log2) ./ sigma_right).^2);
end

function y = logparabola_basis_shared(x_cpd, peakSF, width)
x_cpd = max(x_cpd(:), eps);
y = 10 .^ (-(log10(x_cpd ./ peakSF) ./ width).^2);
end

function y = asym_logparabola_basis_shared(x_cpd, peakSF, width_left, width_right)
x_cpd = max(x_cpd(:), eps);
y = zeros(size(x_cpd));
idxLeft = x_cpd <= peakSF;
idxRight = x_cpd > peakSF;
y(idxLeft) = 10 .^ (-(log10(x_cpd(idxLeft) ./ peakSF) ./ width_left).^2);
y(idxRight) = 10 .^ (-(log10(x_cpd(idxRight) ./ peakSF) ./ width_right).^2);
end

function B = normalize_columns_shared(B)
nrm = sqrt(sum(B.^2, 1));
nrm(nrm == 0) = 1;
B = bsxfun(@rdivide, B, nrm);
end
