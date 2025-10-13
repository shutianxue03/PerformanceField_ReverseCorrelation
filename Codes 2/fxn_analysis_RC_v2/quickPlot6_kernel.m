% close all
figure('Position', [0, 0, 1800, 600])
for iType = 1:nTypes
    % 2D
    subplot(3, nTypes, iType), hold on
    kk = squeeze(kernels2D(iType, :, :));
    imagesc(kk)
    axis square, colorbar
    title(sprintf('%s-%s\nSep = %d%%', namesType{iType},namesLocComb{iLocComb}, round(sep(iType)*100)))
    
    % ORI marginalized kernels
    subplot(3, nTypes, iType+3), hold on
    plot(axis_tuning{1}, mean(kk, 2))
    xlabel('ORI')
    
    % SF marginalized kernels
    subplot(3, nTypes, iType+6), hold on
    plot(axis_tuning{2}, mean(kk, 1))
    xlabel('SF (log)')
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
