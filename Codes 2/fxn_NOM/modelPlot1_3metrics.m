sz_marker = 8;

ticks_allM = {0:.6:3, .5:.1:1, .5:.1:1};
lim_allM = [-.12, 3.12; .48, 1.02; .48, 1.02];

mm = 1:8;
for iim = 1:length(imetric_plot)
    im = mm(imetric_plot(iim));
    subplot(1,3,iim), hold on, box on
    for iiLoc = 1:nLocComb
        errorbar(data_med_allSubj(:, im, iiLoc), pred_med_allSubj(:, im, iiLoc), ...
            pred_neg_allSubj(:, im, iiLoc), pred_pos_allSubj(:, im, iiLoc), ...
            data_neg_allSubj(:, im, iiLoc), data_pos_allSubj(:, im, iiLoc), '.', 'color', colors2{iiLoc}, 'CapSize', 0)
        %         for isubj = 1:nsubj, plot(data_med_allSubj(isubj, im, iiLoc), pred_med_allSubj(isubj, im, iiLoc), [colors2{iiLoc}, markers_allSubj{isubj}]), end
        plot(data_med_allSubj(:, im, iiLoc), pred_med_allSubj(:, im, iiLoc), ...
            'o','MarkerSize', sz_marker, 'MarkerEdgeColor', 'w', 'MarkerFaceColor',  colors2{iiLoc})
        
    end
    plot([lim_allM(iim, 1), lim_allM(iim, 2)], [lim_allM(iim, 1), lim_allM(iim, 2)], 'k-')
    xlim(lim_allM(iim, :))
    ylim(lim_allM(iim, :))
    xticks(ticks_allM{iim})
    yticks(ticks_allM{iim})
    axis square
    title(namesMetrics{im})
end % iim

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
sgtitle(sprintf('%s model\nModel structure: %s\n%s // %s', namesModelA{iModelA}, namesModelB{iModelB}, namesConvolveType{convolveType}, namesIVType{IVType}))
saveas(gcf, sprintf('talk2022/fig/model/n12_Fig1_%dmetrics_A%dB%d_conv%d_IV%d.jpg', length(imetric_plot), iModelA, iModelB, convolveType, IVType))
