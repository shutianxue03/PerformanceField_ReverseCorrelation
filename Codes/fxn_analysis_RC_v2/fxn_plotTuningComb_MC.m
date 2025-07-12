% 
% 
% 
% 
% % ax = subplot(nrow+1, ncol, imodel);
% % ax_all(imodel) = ax;
% % figure('Position', [0 400 400 400])
% 
% if ifeature == 1
%     xticks([-80:40:80])
% else
%     xline(1,'k');
%     xticks([0,1,2]), xticklabels([1,2,4])
% end
% 
% plot(xaxis_itp, yData(1, :), 'ro', 'MarkerSize', 3)
% plot(xaxis_itp, y_pred(1, :), 'r-')
% plot(xaxis_itp, yData(2, :), 'bo', 'MarkerSize', 3)
% plot(xaxis_itp, y_pred(2, :), 'b-')
% yline(0, 'k');
% 
% ylim([ymin_, ymax_])
% % title(sprintf('M%d [%s] %.0f', imodel, num2str(modelInd), AICc(imodel)))
% title(sprintf('%s %d', subjList{isubj}, round(AICc_allSubj(isubj))))
% set(findall(gcf, '-property', 'FontSize'), 'FontSize',12)
% 
% % params_fit = fxn_organizeParams(printFlag, namesParams, nparams_full, nLocTrue, modelInd, kernelParams_est);
% % 
% % for iparam = 1:nparams_full
% %     params_fit_all{iparam} = [params_fit_all{iparam}, params_fit{iparam}];
% % end
% 
% 
