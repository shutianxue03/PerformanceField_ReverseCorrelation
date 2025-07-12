figure('Position', [0 0 2000 600])
    xitv = 0;%1/nplots/10*2;
    yitv = 0;%1/50*2;
    xsz = 1/nSess - xitv;
    ysz = 1/5- yitv;
    ymax = max(kernel_allLoc_allB(:));
    ymin = min(kernel_allLoc_allB(:));
    
    for iLoc = 1:5
        kernel = squeeze(kernel_allLoc_allB(:, iLoc, :));
        nSess = 50 - sum(isnan(kernel(iLoc,:)));
        for isp = 1:nSess
            ax = subplot(5, nSess, isp + (iLoc-1)*nSess); hold on
            %         plot(filterSF_all_log, kernel, 'color', [.5,.5,.5])  % plot all kernels in gray
            plot(filterSF_all_log, kernel(:,isp), '-r', 'linewidth', 2)              % plot the current kernel in red
            plot(filterSF_all_log, kernel_allLoc(iLoc, :), '-k', 'linewidth', 1.5)  % plot the kernel based on all data
            plot(log([2,2]), [ymin, ymax], '-', 'color',[.5,.5,.5]) % the tgt SF
            plot(log([1,4]), [0 0], '-', 'color',[.5,.5,.5]) % kernel = 0
            
            if isp == 1, ylabel(locNames{iLoc}, 'FontSize',15) , end
            if iLoc == 1, title(['S', num2str(isp)], 'FontSize', 15), end
            text(log(1.2), -.3, sprintf('%.2f\n%.2f\n%.4f', dprime_allLoc_allB(isp, iLoc), criterion_allLoc_allB(isp, iLoc), exp(x_cstAll_log(isp, iLoc))) )
            
            %         ax.Position = [xsz * (isp-1), 1-(ysz)*iLoc, xsz, ysz];
            xticklabels([]), yticklabels([])
            xlim(log([1,4])), ylim([ymin, ymax])
        end
        
    end
    % set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
    sgtitle(sprintf('kernel of each plot - %s', subjName), 'FontSize',20)
    
    saveas(gcf, ['html/kernel_byB_', subjName], 'jpg')