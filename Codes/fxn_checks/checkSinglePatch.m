

% inspect each stim
close all
noise.noiseCST = .1;
stim.gaborCST = 1;
sz = stim.aper_sz;

% create a patch
[patch_both, gaborPhase, noise_both, noise_mask_both] = SX_sim03_setStim(1, stim, noise);
figure('Position', [0 0 1500 350])
subplot(2,7,1), imshow(patch_both{1}), title('PRS')
subplot(2,7,8), imshow(patch_both{2}), title('ABS')


for n=1:2
    pp = patch_both{n};

    % get energy
    lumiBG = pp(1,1);
    energy = nan(nfiltersOri, nfiltersSF);
    for ifilterOri = 1:nfiltersOri
        for ifilterSF = 1:nfiltersSF
            [a,b] = SX_sim04_computeEnergy(pp, lumiBG, filter_sin{ifilterOri, ifilterSF}, filter_cos{ifilterOri, ifilterSF}, 0, 1, filtersSF_all(ifilterSF), filtersOri_all(ifilterOri));
            energy(ifilterOri, ifilterSF) = a;
        end
    end
    
    % plot marginalized energy of ORI
    if n==1, subplot(2,7,2), else, subplot(2,7,9), end, hold on
    plot(filtersOri_all-90, mean(energy, 2)), xlabel('ORI')
    xline(0, 'r');
    xlim([-80, 80])
    xticks(-80:40:80)
    title('marginalized energy')
    
    % plot marginalized energy of SF
    if n==1, subplot(2,7,3), else, subplot(2,7,10), end, hold on
    plot(filtersSF_all_log, mean(energy, 1)), xlabel('SF')
    xline(1, 'r');
    xticks([0,1,2]), xticklabels([1,2,4])
    
    % 2D FFT
    fft2D = abs(fft2(pp));
    if n==1, subplot(2,7,4), else, subplot(2,7,11), end, hold on
    imagesc(fftshift(fft2D(2:end,2:end))), axis square
    title('FFT 2D')
    
    % 1D FFT
    fft1D = mean(abs(fft2(pp)),2);
    if n==1, subplot(2,7,5), else, subplot(2,7,12), end, hold on
    stem(log2((1:9)/sz), fft1D(2:10), 'k')
    xline(1, 'r');
    xticks([0,1,2]), xticklabels([1,2,4])
    xlim([0,2])
    xlabel('SF'), ylabel('power')
    title('mean(abs(fft2D)')
    
    % 1D FFT
    fft1D = abs(mean(fft2(pp),2));
    if n==1, subplot(2,7,6), else, subplot(2,7,13), end, hold on
    stem(log2((1:9)/sz), fft1D(2:10), 'k')
    xline(1, 'r');
    xticks([0,1,2]), xticklabels([1,2,4])
    xlim([0,2])
    xlabel('SF'), ylabel('power')
    title('abs(mean(fft2D))')
    
    % radial average
    [x_rad, ave_rad] = radial_profile(fftshift(fft2D(2:end, 2:end)), 1);
    if n==1, subplot(2,7,7), else, subplot(2,7,14), end, hold on
    stem(x_rad(1:10), ave_rad(1:10), 'k')
%     xline(2, 'r');
    xticks([1,5,10])
%     xlim([1,4])
    xlabel('unknown'), ylabel('power')
    title('radial average')
    
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',12)
sgtitle(sprintf('noise cst = %.1f, Gabor cst = %.1f', noise.noiseCST, stim.gaborCST), 'FontSize',20)

function [x_rad, ave_rad]=radial_profile(data, radial_step)
%main axii specified:
x=(1:size(data,2))-size(data,2)/2;
y=(1:size(data,1))-size(data,1)/2;
% coordinate grid:
[X,Y]=meshgrid(x,y);
% creating circular layers
Z_integer=round(abs(X+1i*Y)/radial_step)+1;
% % illustrating the principle:
% figure;imagesc(Z_integer.*data)
% very fast MatLab calculations:
x_rad = accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
ave_rad = accumarray(Z_integer(:),data(:),[],@mean);
end