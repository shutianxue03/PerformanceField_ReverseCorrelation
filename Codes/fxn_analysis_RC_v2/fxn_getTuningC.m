
function tuningC = fxn_getTuningC(x, iFeature, iFamily, pred, params)

% function tuningC = fxn_getTuningC(x, ifeature, pred, param)
% extract the characteristics of the tuning function
% ORI:
%   1. peak: the height of the top (assumed to center at the signal ORI, i.e., 0 deg)
%   2. half-width: the range of ORI at half height (i.e., (max-min)/2)
%   3. bottom: the left/rightmost point of the tuning curve; the same as the estimated baseline (or b)
%   4. half-height
% SF
% if ifamily=3 (with truncation):1-10
% if ifamily=2 (no truncation):1-9 ()
%   1. peakSF: at which SF the peak is found [ln, to be consistent with the estimated params]
% . 2. peak: the height of the top
%   3. full width (in octave): the range of ORI at half height on log scale
%   4. baseline: the rightmost point of the tuning curve; NOT sthe same as the estimated baseline (or b)
%   5. truncation: the leftmost point of the tuning curve
%   6/5. left width
%   7/6. right width
%   8/7. bottom
%   9/8. half height
%   10/9. full width (in cpd): on linear scale

nfilters = length(x);

if iFeature == 1
    switch iFamily
        case 1 % scaled Gaussian
            % ORI-peak and bottom
            peakAmp_ORI = pred(ceil(nfilters/2));
            bottom_ORI = pred(1);
            width_ORI = params(2)*sqrt(2*log(2));
            % compile
            tuningC = [peakAmp_ORI, width_ORI, bottom_ORI];

        case 8 % DoG
            baseline = pred(1);

            % SX_normPDF = @(x,mu,sigma) exp(-((x-mu)/sigma).^2);
            % fxn = @(x, params) params(1)*(exp(-(x/params(3)).^2) - params(2)*exp(-(x/params(4)).^2)) + params(5);
            % peak amp
            peakAmp_ORI = pred(ceil(nfilters/2));
            % trough depth
            trough_depth = min(pred);
            %% trough ori
            trough_ORI = inverseFxn(trough_depth, iFamily, params, linspace(0, 90, 1e4));
            % ORI - width
            pred_max = max(pred); pred_min = min(pred);
            pred_half_max_min = (pred_max - pred_min)/2 + pred_min; % the y coordinate of the width
            %             pred_half_max_base = (pred_max - baseline)/2 + baseline; % the y coordinate of the width
            width_ORI_min = 2*inverseFxn(pred_half_max_min, iFamily, params, linspace(0, 90, 1e4));
            %             width_ORI_base = 2*inverseFxn(pred_half_max_base, ifamily, params, linspace(0, 90, 1e4));
            % compile
            tuningC = [peakAmp_ORI, trough_ORI, trough_depth, width_ORI_min, baseline];
            %%
            % figure, hold on,
            % plot(x, pred)
            % xline(peak_ORI);
            % xline(trough_ORI);
            % yline(params(5));
            % yline(pred_half_max_min);
            % yline(pred_half_max_base);
    end


else % iF=2, SF
    x_ln = linspace(1,4, 1e3);

    switch iFamily
        case 12
            pred_itp = predSFkernel(x_ln, iFamily, params, 0);
            localMax_ = islocalmax(pred_itp);
            localMax = find(localMax_==1);
            if isempty(localMax), disp('NO local max'), [~, peakSF1_ln] = max(pred_itp); peakSF1_ln = x_ln(peakSF1_ln);
            else, peakSF1_ln = x_ln(localMax(1));
            end
            if length(localMax)<2
                if abs(peakSF1_ln-params(5)) <= abs(peakSF1_ln-params(1))
                    peakSF2_ln = params(1);
                else
                    peakSF2_ln = params(5);
                end
                if peakSF1_ln>peakSF2_ln
                    peakSF1_ln_ = peakSF2_ln;
                    peakSF2_ln = peakSF1_ln;
                    peakSF1_ln = peakSF1_ln_;
                end
            else, peakSF2_ln = x_ln(localMax(2));
            end
            peakAmp1 = predSFkernel(peakSF1_ln, 12, params, 0);
            peakAmp2 = predSFkernel(peakSF2_ln, 12, params, 0);

            if params(1) <= params(5)
                % pred1 = predSFkernel(x_ln, 2, params(1:4), 0);
                % pred2 = predSFkernel(x_ln, 2, params(5:8), 0);
                %===============%
                [bw_log1, bw_ln1] = getSFbandwidth(params(1:4));
                [bw_log2, bw_ln2] = getSFbandwidth(params(5:8));
                %===============%
            else
                % pred1 = predSFkernel(x_ln, 2, params(5:8), 0);
                % pred2 = predSFkernel(x_ln, 2, params(1:4), 0);
                %===============%
                [bw_log1, bw_ln1] = getSFbandwidth(params(5:8));
                [bw_log2, bw_ln2] = getSFbandwidth(params(1:4));
                %===============%
            end
            assert(peakSF1_ln <= peakSF2_ln)

            tuningC = [peakSF1_ln, peakAmp1, bw_log1, peakSF2_ln, peakAmp2, bw_log2];

            %% plot for F12
            figure, hold on
            pred = predSFkernel(x_ln, iFamily, params, 0);
            plot(x_ln, pred, 'k-', 'linewidth', 4), yline(0, 'k-');
            pred1 = predSFkernel(x_ln, 2, params(1:4), 0);
            pred2 = predSFkernel(x_ln, 2, params(5:8), 0);
            plot(x_ln, pred1, 'r-', 'linewidth', 2)
            plot(x_ln, pred2, 'b-', 'linewidth', 2)

            localMax = islocalmax(pred);
            plot(x_ln(localMax), pred(localMax), 'r*', 'MarkerSize', 30)
            plot(([peakSF1_ln, peakSF1_ln]), [0, peakAmp1], 'k-', 'linewidth', 2)
            plot(([peakSF2_ln, peakSF2_ln]), [0, peakAmp2], 'k-', 'linewidth', 2)
            errorbar(peakSF1_ln, peakAmp1, bw_ln1/2, 'r', 'horizontal')
            errorbar(peakSF2_ln, peakAmp2, bw_ln2/2, 'b', 'horizontal')

        case 2
            peakSF = params(1);
            peakAmp = max(pred);
            baseline = params(4);
            %===============%
            [width_full_log] = getSFbandwidth(params);
            %===============%
            tuningC = [peakSF, peakAmp, width_full_log, baseline];
            %%
        case 3
            peakAmp = max(pred);
            peakSF = params(1);
            % bottom_SF = pred(end);
            baseline = params(4);
            % peakSF and peak (from the true fxn without the trunc)
            % pred_true = predSFkernel(x_ln, 2, params(1:4), 0); % as in some subj, the trunc is even higher that the peak of the true kernel
            %===============%
            [width_full_log] = getSFbandwidth(params);
            %===============%
            trunc_SF = pred(1);
            tuningC = [peakSF, peakAmp, width_full_log, baseline, trunc_SF];%, ...
            %                 width_L_SF, width_R_SF, bottom_SF, pred_half, width_full_ln];
    end
end % if ifeature==1

end % fxn

%%
function x = inverseFxn(y, ifamily, params, x_all)
dev = nan(length(x_all), 1);
for ii = 1:length(x_all)
    x_p = x_all(ii);
    %     y_sim = fxn(x_p, params);
    y_sim = predSFkernel(x_p, ifamily, params, 0);
    dev(ii) = abs(y-y_sim);
end
[~, ix] = min(dev);
x = x_all(ix);
end


function [width_full_oct, width_full_cpd, fL, fR] = getSFbandwidth(params)
% getSFbandwidth  Bandwidth (FWHH) for SF tuning in OCTAVES and CPD
%
% Model in code:
%   Ksf(f) = base + gain * 10.^(-(log10(f/peakSF)/width).^2)
%
% params = [peakSF, gain, width, base]
%
% Outputs:
%   width_full_oct : full width at half height in octaves
%   width_full_cpd : full width at half height in cycles/deg (linear scale)
%   fL, fR         : half-height frequencies (left/right), in cpd

peakSF = params(1);
w      = params(3);

% Half-height constant for base-10 log-parabola:
% 10^(-(log10(f/f0)/w)^2) = 1/2  ->  (log10(f/f0)/w)^2 = log10(2)
c = sqrt(log10(2));

fL = peakSF * 10^(-w * c);
fR = peakSF * 10^(+w * c);

width_full_oct = log2(fR / fL);   % octaves
width_full_cpd = fR - fL;         % linear bandwidth (cpd)

end


% function [width_full_log, width_full_ln] = getSFbandwidth(params, pred_true, x_)
% [pred_max, imax] = max(pred_true);
% % SF-width
% pred_min = min(pred_true);
% pred_half = (pred_max - pred_min)/2 + pred_min; % the y coordinate of the width
% x_half_L_SF = params(1) * 10^(-(params(3)*sqrt(-log10((pred_half-params(4))/params(2))))); % linear scale
% x_half_R_SF = params(1) * 10^(params(3)*sqrt(-log10((pred_half-params(4))/params(2)))); % linear scale
% % get bandiwdth in octave unit (on log scale)
% width_full_log = log2(x_half_R_SF/x_half_L_SF);
% % get bandiwdth in cpd unit (on linear scale)
% width_L_SF = x_(imax)-x_half_L_SF; % left, linear scale
% width_R_SF = x_half_R_SF-x_(imax); % right, linear scale
% width_full_ln = width_L_SF + width_R_SF; % linear scale
% end

