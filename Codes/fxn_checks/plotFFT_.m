
function plotFFT_(pp_fft2_2D, margDim)
pp_fft2_2D_ave = squeeze(mean(pp_fft2_2D, 1));
pp_fft2_1D_ave = squeeze(mean((pp_fft2_2D_ave), margDim));

figure('position', [0 800 600 300])

subplot(1,2,1)
imagesc(fftshift(pp_fft2_2D_ave(2:end, 2:end)))
% imagesc(fftshift(pp_fft2_2D_ave))
axis square
colorbar, colormap gray
title('FFT of patches')

subplot(1,2,2)
% stem((1:20)/2, pp_fft2_1D_ave(2:21))
stem((1:20)/2, pp_fft2_1D_ave(2:21))
xline(2, 'r');
xline(4, 'g');
xticks([1:4, 5, 10])
% ylim([0, 30])
set(findall(gcf, '-property', 'FontSize'), 'FontSize',20)
end
