

% plot data of two loci (dots) and predictions (number of lines depends on the type of model)

figure('Position', [3e3, 200, 600, 500]), hold on

for iline = 1:2
    % data
    plot(axis_tuning_ln, yData(iline, :), 'o', 'color', colors2{iline})
    
    % prediction on all observations
    if MCmode <= 2
        plot(axis_tuning_ln,  squeeze(median(y_pred_allK(:, iline, :),1)) , '-', 'color', colors2{iline})
    else
        plot(axis_tuning_ln,  y_data_pred(iline, :) , '-', 'color', colors2{iline}, 'MarkerSize',  5)
        % prediction on all observations
    end
    
    %     ylim([-.03,.1])
    %     yline(0,'k-');
    %     xlim([.9, 4.1])
    %     xticks([1,2,4])
    %     xticklabels([1,2,4])
end % end of iline

%% title
if MCmode <= 2  % 10 fold CV
    title(sprintf('%s - model #%d [%s] \ndev of the test set: %.5f (%.5f)', ...
        subjName, imodel, num2str(paramInd), nanmedian(dev_test_allK), nanstd(dev_test_allK)), 'FontSize',20)
else % IC
    %     title(sprintf('%s %s \ndev: %.5f\ninfo criterion: %.5f', ...
    %         namesMCmode{MCmode}, num2str(paramInd), dev, dev_ave(imodel)), 'FontSize',20)
    title(sprintf('%s - model #%d [%s]\nAIC (%.2f) AICc (%.2f) BIC (%.2f)', ...
        subjName, imodel, num2str(paramInd), squeeze(dev_allSubj(isubj, imodel, :))), 'FontSize',20)
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize',14)


%% inset for CV only
if MCmode <= 2 % 10 fold CV
    % define the coordinate of the inset figure (values are prop. in the figure)
    xstart=.7; xend=.9; ystart=.5; yend=.6;
    axes('position',[xstart ystart xend-xstart yend-ystart ])
    
    histogram(dev_test_allK, 10)
    xlabel('deviance')
    xline(nanmean(dev_test_allK), 'r-');
    xline(nanmedian(dev_test_allK), 'r--');
end
