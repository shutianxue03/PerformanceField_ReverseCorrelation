function [filter_sin, filter_cos] = SX_sim02_setFilters(stim, filtersSF_all, filtersOri, fxn_getSigma_SPdomain, flag_plotFilterPixels)

% This script is located in fxn_analysis_RC_v2

cstOfFilter = 1; % the contrast of the filter, 100%

nfiltersSF = length(filtersSF_all);
nfiltersOri = length(filtersOri);

filter_sin = cell(nfiltersOri, nfiltersSF);
filter_cos = filter_sin;

if flag_plotFilterPixels
    figure('Position', [0 0 2e3 800])
    iplot = 1;
    oo_all = 1;
    ff_all = 1:7:29;
    nSF = length(ff_all);
else
    oo_all = 1:nfiltersOri;
    ff_all = 1:nfiltersSF;
end


for oo = oo_all
    stim.targetOri = filtersOri(oo);
    for ff = ff_all
        stim.gaborSF = filtersSF_all(ff); % linear SF in unit of ppd

        % adjust the bandwidth to ensure constant bandwidth in octave
        %         stim.gaborSD = stim.gaborSF * bw_scaling;
        stim.gaborSD = fxn_getSigma_SPdomain(stim.gaborSF);

        for n = 1:2 % sin and cos
            stim.phase = pi/2 * (n-1); % sin: 0, cos: pi/2
            % use SX's code
            %------------------------------%
            patch = exp_CreateGabor(stim, cstOfFilter, 0);% .* stim.mask;
            %------------------------------%
            if n==1, filter_sin{oo, ff} = patch; else, filter_cos{oo, ff} = patch; end
        end

        if flag_plotFilterPixels
            % subplot(nfiltersSF, nfiltersOri, iplot)
            % Spatial domain
            subplot(2, nSF, iplot)
            imagesc(patch)
            axis square
            title(sprintf('ORI = %dº, SF = %.2f cpd', stim.targetOri, stim.gaborSF))
            xlabel('x pixel')
            ylabel('y pixel')

            % Fourier domain
            subplot(2, nSF, iplot+nSF)

            % Compute 2D FFT and magnitude spectrum
            fft_patch = fftshift(fft2(patch));
            mag_fft = log(abs(fft_patch) + 1);  % add 1 to avoid log(0)

            % Frequency axis in cycles per degree (cpd)
            x_freqs_cpd = linspace(-stim.ppd/2, stim.ppd/2, stim.gabor_psz);
            imagesc(x_freqs_cpd, x_freqs_cpd, mag_fft)
            axis square
            title('FFT magnitude (log-scaled)')

            iplot = iplot + 1;
            
        end
    end
end

if flag_plotFilterPixels,
    sgtitle('Selected Gabor Filters (gabor SD is normalized by center SF)')
    set(findall(gcf, '-property', 'fontsize'), 'fontsize', 20)
    saveas(gcf, 'Fig/Sim_NOM_TrialWise/example_gabor_filters_SDnormalized.jpg')    
end