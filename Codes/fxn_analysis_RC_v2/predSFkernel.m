function y = predSFkernel(x, str_family, params, flag_plot)

SX_normPDF = @(xv,mu,sigma) exp(-((xv-mu)/sigma).^2);

switch str_family
    case 'gaussian_zero_mean' % ORI
        gain = params(1);
        width_deg = params(2);
        base = params(3);
        y = gain * SX_normPDF(x, 0, width_deg) + base;

    case 'log_parabola' % SF
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        base = params(4);
        y = gain * 10.^(-(log10(x / peakSF) / width).^2) + base;

    case 'log_parabola_truncated' % SF
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        base = params(4);
        lambda  = params(5);
        y = gain * 10.^(-(log10(x / peakSF) / width).^2) + base;
        y((x < peakSF) & (y < lambda)) = lambda;

    case 'raised_gaussian' % ORI
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        base = params(4);
        power  = params(5);
        y = gain * exp(-1/2*(x-peakSF).^2/width^2) .^power + base;

    case 'double_exponential' % ORI
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        y = gain * x .^ (peakSF / width) .* exp(-x/width);

    case 'skewed_gaussian' % ORI
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        y = 2*SX_normPDF(x, peakSF, width).*SX_normPDF(gain*x, peakSF, width);

    case 'raised_gaussian_truncated' % SF
        peakSF = params(1);
        gain = params(2);
        width = params(3);
        base = params(4);
        power  = params(5);
        lambda = params(6);
        y = gain * SX_normPDF(x,peakSF,width).^ power + base;
        y(x < peakSF & y < lambda) = lambda;

    case 'DoG' % Difference of Gaussians
        gain1 = params(1);
        gain2 = params(2);
        sigma1 = params(3);
        sigma_r = params(4);
        baseline = params(5);
        y = gain1*(SX_normPDF(x, 0, sigma1) - gain2*SX_normPDF(x, 0, sigma1*sigma_r)) + baseline;

    case 'von_mises' % ORI
        gain = params(1);
        k = params(2);
        base = params(3);
        theta = x / 180 * 2 * pi;
        y = gain * v_vonmisespdf(theta, 0, k) + base;

    case 'double_peak' % SF
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

    case 'raised_DoG' % SF
        gain1 = params(1);
        gain2 = params(2);
        sigma1 = params(3);
        sigma_r = params(4);
        power = params(5);
        baseline = params(6);
        y = gain1*(SX_normPDF(x, 0, sigma1).^power - gain2*SX_normPDF(x, 0, sigma1*sigma_r))+ baseline;

    case 'asym_gaussian' % SF
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

    case 'circ_gaussian_basis' % ORI
        center_deg = params(1);
        sigma_deg = params(2);
        period_deg = params(3);
        distance = mod(x - center_deg + period_deg/2, period_deg) - period_deg/2;
        y = exp(-0.5 * (distance ./ sigma_deg).^2);
        y = y(:);

    case 'von_mises_basis' % ORI
        center_deg = params(1);
        kappa = params(2);
        period_deg = params(3);
        distance = mod(x - center_deg + period_deg/2, period_deg) - period_deg/2;
        theta = distance / period_deg * 2*pi;
        y = exp(kappa * cos(theta));
        y = y(:);

    case 'log2_gaussian' % SF
        center_log2 = params(1);
        sigma_log2 = params(2);
        y = exp(-0.5 * ((x - center_log2) ./ sigma_log2).^2);
        y = y(:);

    case 'asym_log2_gaussian' % SF
        center_log2 = params(1);
        sigma_left = params(2);
        sigma_right = params(3);
        x_log2 = x(:);
        y = zeros(size(x_log2));
        idxLeft = x_log2 <= center_log2;
        idxRight = x_log2 > center_log2;
        y(idxLeft) = exp(-0.5 * ((x_log2(idxLeft) - center_log2) ./ sigma_left).^2);
        y(idxRight) = exp(-0.5 * ((x_log2(idxRight) - center_log2) ./ sigma_right).^2);

    case 'log_parabola_basis' % SF
        peakSF_cpd = params(1);
        width = params(2);
        x_cpd = max((2 .^ x(:)), eps);
        y = 10 .^ (-(log10(x_cpd ./ peakSF_cpd) ./ width).^2);

    case 'asym_log_parabola_basis' % SF
        peakSF_cpd = params(1);
        width_left = params(2);
        width_right = params(3);
        x_cpd = max((2 .^ x(:)), eps);
        y = zeros(size(x_cpd));
        idxLeft = x_cpd <= peakSF_cpd;
        idxRight = x_cpd > peakSF_cpd;
        y(idxLeft) = 10 .^ (-(log10(x_cpd(idxLeft) ./ peakSF_cpd) ./ width_left).^2);
        y(idxRight) = 10 .^ (-(log10(x_cpd(idxRight) ./ peakSF_cpd) ./ width_right).^2);

    otherwise
        error('Unknown family: %s', string(str_family));
end

if flag_plot
    figure, hold on
    plot(x, y, 'o-')
    switch str_family
        case {'gaussian_zero_mean','raised_gaussian','double_exponential','skewed_gaussian','DoG','von_mises','circ_gaussian_kernel','von_mises_basis'}
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
