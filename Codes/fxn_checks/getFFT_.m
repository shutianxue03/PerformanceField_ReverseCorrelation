function [fft2D, fft1D] = getFFT_(pp, margDim, iplot)
fft2D = abs(fft2(pp));
fft1D = squeeze(mean(fft2D, margDim));

if iplot
    sz = size(pp, 1);
    x = (1:20)/2;
    yind = 2:21;
    subplot(5,3,iplot), imagesc(fftshift(fft2D(2:sz, 2:sz))), axis square
    subplot(5,3,iplot+1), stem(x, fft1D(yind)), xticks([1:5, 10]), xline(2, 'r'); xline(4, 'g');
end
