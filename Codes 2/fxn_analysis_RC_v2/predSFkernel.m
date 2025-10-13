
function y = predSFkernel(x, ifamily, params, plotFlag)

SX_normPDF = @(x,mu,sigma) exp(-((x-mu)/sigma).^2);

switch ifamily
    case 1 % gaussian
        gain = params(1);
        width = params(2);
        base = params(3);
        y = gain * SX_normPDF(x, 0, width) + base;
        
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

    case 10 % von mises for OR gain = params(1);
        gain = params(1);
        k = params(2); % concentration parameter
        base = params(3);
        
        y = gain * v_vonmisespdf(x/180*pi, 0, k) + base;
    case 11 % gaussian with mean free to vary
        peak = params(1);
        gain = params(2);
        width = params(3);
        base = params(4);
        %         kernel = gain * exp(-1/2*x.^2/width^2) + base;
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
        %         y = gain1 * 10.^(-(log10(x / peakSF1) / width1).^2) + gain2 * 10.^(-(log10(x / peakSF2) / width2).^2) + base;
        y1= gain1 * 10.^(-(log10(x / peakSF1) / width1).^2)+base1;
        y2= gain2 * 10.^(-(log10(x / peakSF2) / width2).^2)+base2;
        y = y1+y2;
%         figure, hold on, plot(x, y1, 'b--'), plot(x, y2, 'r--'), plot(x, y, 'k-')
    case 13 % raised DoG
        gain1 = params(1);
        gain2 = params(2);
        sigma1 = params(3);
        sigma_r = params(4);
        power = params(5);
        baseline = params(6);
        y = gain1*(SX_normPDF(x, 0, sigma1).^power - gain2*SX_normPDF(x, 0, sigma1*sigma_r))+ baseline;
end

if plotFlag
    figure, hold on
    plot(x, y, 'o-')
    ylimit = ylim;
    plot([2,2], ylimit, 'Color', [.5,.5,.5])
    xlabel('SF channels (cpd)')
    ylabel('sensitivity kernel (a.u.)')
    xlim([1,4])
    xticks([1,2,4]), xticklabels([1,2,4])
    text(3.5, 1, sprintf('f0 = %.2f\nsigma = %.2f\nalpha = %.2f\nb = %.2f\nlambda = %.2f', params))
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
end
