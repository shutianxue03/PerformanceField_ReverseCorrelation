

input('WARNING: Is pA_allSubj defined correctly? ');
figure('Position', [0 500 2000 200])
for iLoc = 1:nLoc8
    y = pA_allSubj(:, iLoc, 1);
    x = pA_allSubj(:, iLoc, 2);
    
    subplot(1,8,iLoc), hold on
    plot(x,y,'o')
    errorbar(mean(x), mean(y), std(x)/sqrt(nsubj), 'horizontal', 'CapSize', 0, 'linewidth', 2, 'color', colors_comb(iLoc, :))
    errorbar(mean(x), mean(y), std(y)/sqrt(nsubj), 'vertical', 'CapSize', 0, 'linewidth', 2, 'color', colors_comb(iLoc, :))
    
    plot([.5, 1], [.5, 1], '-k')
    ylim([.5, 1])
    xlim([.5, 1])
    xticks([.5, .75, 1]), yticks([.5, .75, 1])
    if iLoc == 1
        ylabel('pA [PRS]'), xlabel('pA [ABS]')
    end
    axis square, box on
    title(namesLocComb{iLoc})
end