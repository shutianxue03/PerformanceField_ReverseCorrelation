

figure('Position', [0 0 1500 600])

sz = stim.gabor_sz;
for iLoc = 1:nLoc8

    % 1. CI
    for n=1
        switch n, case 1, CI_ = CI_mean{iLoc}; case 2, CI_ = CI_var{iLoc}; end    
        subplot(3, nLoc8, iLoc)
        imshow(CI_ + 0.5)
        if iLoc== 1, ylabel(namesCI{n}), end
        title(namesLocComb{iLoc})
    end
    
    % 2. 2D-FFT of CI
    CI_fft2_2D = abs(fft2(CI_));
    subplot(3,nLoc8, iLoc+nLoc8), imagesc(fftshift(CI_fft2_2D)), axis square
    
    % 3. 1D FFT of CI
    CI_fft2_1D = CI_fft2_2D(:, 1);
    subplot(3, nLoc8, iLoc+nLoc8*2)
    stem((1:18)/sz, CI_fft2_1D(1:18))
    xline(2, 'r');
    xline(4, 'r');
    xticks([1,2,4])
%     ylim([0, 5])
    
    
    %     for n=1:4
    %         switch n, case 1, pp = CI_mean{iLoc}; case 2, pp = CI_mean_norm{iLoc}; case 3, pp = CI_var{iLoc}; case 4, pp = CI_var_norm{iLoc}; end
    %         subplot(nLoc, 4, n+4*(iLoc-1))
    %         imshow(pp + 0.25)
    %         if iLoc== 1, title(titles_CI{n}), end
    %     end
end

clear pp
sgtitle([subjName, ' - classification image'])
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)

