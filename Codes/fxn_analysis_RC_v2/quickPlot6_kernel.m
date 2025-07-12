% close all
figure('Position', [0, 0, 1800, 600])
for itype = 1:ntypes
    % 2D
    subplot(3, ntypes, itype), hold on
    kk = squeeze(kernels2D(itype, :, :));
    imagesc(kk)
    axis square, colorbar
    title(sprintf('%s-%s\nSep = %d%%', namesType{itype},namesLocComb{iLocComb}, round(sep(itype)*100)))
    
    % ORI marginalized kernels
    subplot(3, ntypes, itype+3), hold on
    plot(axis_tuning{1}, mean(kk, 2))
    xlabel('ORI')
    
    % SF marginalized kernels
    subplot(3, ntypes, itype+6), hold on
    plot(axis_tuning{2}, mean(kk, 1))
    xlabel('SF (log)')
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
