% iframe = 1;
nframes = size(noiseImg_, 3);
figure('Position', [3e3 0 1000 500]),
for iframe = 1:nframes
    subplot(2,4,1), imshow(gabor), title('Gabor ONLY')
    subplot(2,4,2), imshow(squeeze(noiseImg(:, :, iframe, td.targetLoc))), title('Noise ONLY')
    subplot(2,4,3), imshow(squeeze(noiseImg_(:, :, iframe, td.targetLoc))), title('Gabor+Noise')
    subplot(2,4,4), imshow(squeeze(stimImg(:, :, iframe, td.targetLoc))), title('Gabor+Noise (rescaled)')
    
    subplot(2,4,5), imagesc(gabor), axis square, colorbar
    subplot(2,4,6), imagesc(squeeze(noiseImg(:, :, iframe, td.targetLoc))), axis square, colorbar
    subplot(2,4,7), imagesc(squeeze(noiseImg_(:, :, iframe, td.targetLoc))), axis square, colorbar
    subplot(2,4,8), imagesc(squeeze(stimImg(:, :, iframe, td.targetLoc))), axis square, colorbar
    
    sgtitle(sprintf('Noise level = %d%% [#%d]\nFrame #%d', ...
        round(params.extNoiseLvl(td.extNoiseLvl)*100), td.extNoiseLvl, iframe))
    pause
end