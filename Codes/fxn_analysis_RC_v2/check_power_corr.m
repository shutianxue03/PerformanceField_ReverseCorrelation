clc
close all

pp_all = patch0_AF_allT;
% pp_all = patch_all;
% get fft 1D
nfft = 10;
fft_1D_all = nan(ntrials, nfft);

figure, fft_2D = abs(fft2(pp_all{1})); imagesc((fft_2D(2:end, 2:end)))
for itrial = 1:ntrials
    fft_2D = abs(fft2(pp_all{itrial}));
%     fft_1D_all(itrial, :) = fft_2D(1, 2:nfft+1);
    fft_1D_all(itrial, :) = fft_2D(2:nfft+1, 1);
end

% get correlation
r = nan(nfft, nfft);
for i1 = 1:nfft
    fft_fix = fft_1D_all(:, i1);
    for i2 = 1:nfft
        fft_comp = fft_1D_all(:, i2);
        r(i1,i2) = corr(fft_fix, fft_comp);
    end
end

%%%%% 
% plot %
%%%%% 
x = 1:nfft;
figure('Position', [2000 0 800 300])

% power of each trial
subplot(1,2,1), hold on
plot(x, fft_1D_all', 'color', [.5, .5, .5])
plot(x, mean(fft_1D_all, 1), 'k', 'linewidth', 2)
xticks(x(1:3:end)), xticklabels(round(x(1:3:end)/sz, 1))
xlabel('SF'), ylabel('power'), title('Power spectrum of each trial')
xlim([0,nfft])
xline(sz, 'r'); xline(4*sz, 'r');
xline(sz*sf, 'c');

% correlation
subplot(1,2,2)
imagesc(r), colorbar
xticks(x(1:3:end)), xticklabels(round(x(1:3:end)/sz, 1))
yticks(x(1:3:end)), yticklabels(round(x(1:3:end)/sz, 1))
title('Correlation')
set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
sgtitle('Invalid - PRSENT trials')
% sgtitle('simulated ABSENT trials (no smoothing)')

%%
figure, hold on
plot(mean(fft_1D_all, 1), 'k', 'linewidth', 2)
for it = 1:ntrials
    plot(1:nfft, fft_1D_all(it, :), 'color', [.5, .5, .5])
    %     ylim([0, 200])
    title(['itrial = ', num2str(it)])
    xticks(1:nfft), xticklabels((1:nfft)/sz)
    xlabel('SF'), ylabel('power')
    waitforbuttonpress
end

