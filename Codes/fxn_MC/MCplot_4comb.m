%         
% printFlag = 0; % 0=do NOT print the estimated params
% imodel_opt = nan(ncomb , ntypes);
% imodel_opt_IC = nan(ncomb , ntypes);
% % imodel_opt_IC_delta = nan(ncomb , ntypes);
% imodel_opt_paramInd = cell(ncomb , ntypes);
% 
% for icomb = 1:ncomb
%     for itype = 1:ntypes
%         
%         IC_delta = IC_delta_all{icomb, itype};
%         [IC_delta_sort, iIC_delta_sort] = sort(IC_delta);
%         IC = IC_all{icomb, itype};
%         IC_sort = IC(iIC_delta_sort);
%         nmodels_s = min([5, nmodels]); % only show the first 5 best-fitting models
%         
%         % record the optimal model 
%         imodel_opt(icomb, itype) = iIC_delta_sort(1);
%         imodel_opt_IC(icomb, itype) = IC_sort(1);
% %         imodel_opt_IC_delta(icomb, itype) = IC_delta_sort(1);
%         imodel_opt_paramInd{icomb, itype} = paramInd_all(iIC_delta_sort(1), :);
% 
%         figure('Position', [2500 0 nsubj*250 nmodels_s*200])
%         for imodel = 1:nmodels_s
%             imodel_true = iIC_delta_sort(imodel);
%             yData_allSubj = yData_allSubj_all{icomb, itype, imodel_true};
%             yPred_allSubj = yPred_allSubj_all{icomb, itype, imodel_true};
%             IC_allSubj = IC_allSubj_all{icomb, itype, imodel_true};
%             
%             ymin_ = min([yData_allSubj(:); yPred_allSubj(:)]);
%             ymax_= max([yData_allSubj(:); yPred_allSubj(:)]);
%             
%             for isubj = 1:nsubj
%                 subplot(nmodels_s+1, nsubj, nsubj*(imodel-1)+isubj), hold on
%                 
%                 yData = squeeze(yData_allSubj(isubj, :, :));
%                 y_pred = squeeze(yPred_allSubj(isubj, :, :)); % could change to y_pred_nLL if needed
%                 
%                 buffer = .05;
%                 if ifeature == 1
%                     xline(0,'k');
%                     xticks([-80,0,80])
%                     xaxis = xaxis_itp;
%                 else
%                     xline(1,'k');
%                     xticks([0,1,2]), xticklabels([1,2,4])
%                     xaxis = log2(xaxis_itp);
%                 end
%                 
%                 for n=1:2
%                     plot(xaxis, yData(n, :), 'o', 'MarkerSize', 3, 'color', colors_comb(combInd(icomb, n), :))
%                     plot(xaxis, y_pred(n, :), '-', 'color', colors_comb(combInd(icomb, n), :))
%                     % R2
%                     R2 = getR2(yData(n, :), y_pred(n, :));
%                     text(xaxis(2), ymax_ - buffer*n, sprintf('%d%%', round(R2*100)), 'color', colors_comb(combInd(icomb, n), :))
%                 end
%                 yline(0, 'k');
%                 
%                 
%                 % highlight the subj with the best fitting
%                 if IC_allSubj(isubj) == min(IC_allSubj)
%                     ax = gca;
%                     set(ax,'color',[1 .8 .75])
%                     set(ax,'box','on')
%                 end
%                 
%                 ylim([ymin_, ymax_])
%                 
%                 % xticks and yticks
%                 if (isubj > 1) || (imodel > 1)
%                     xticks([])
%                     yticks([])
%                 end
%                 
%                 % y label: M# and [xx]
%                 if isubj==1
%                     ylabel(sprintf('M%d [%s]', imodel_true, num2str(paramInd_all(imodel_true, :))))
%                 end
%                 
%                 % title: nameSubj, AICc
%                 title(sprintf('%s %d', subjList{isubj}, round(IC_allSubj(isubj))))
%                 
%             end % end of isubj
%         end % end of imodel
%         
%         %%%%%%%%%
%         %    plot AICc    %
%         %%%%%%%%%
%         subplot(nmodels_s+1, nsubj, [nsubj*nmodels_s+1:(nmodels_s+1)*nsubj]), hold on
%         bar(IC_delta_sort, 'BarWidth', .5, 'EdgeColor', colorsType{itype}, 'FaceColor', 'w', 'linewidth', 2)
%         for im = 1:nmodels
%         text(im, IC_delta_sort(im)+std(IC_delta_sort)/2, num2str(round(IC_sort(im))), 'HorizontalAlignment', 'center', 'fontsize', 15)
%         end
%         
%         % x tick labels
%         for imodel = 1:nmodels
%             xticklabels_ = [];
%             for ip = 1:nparams_full
%                 xticklabels_ = sprintf('%s\\newline%d', xticklabels_ , paramInd_all(iIC_delta_sort(imodel), ip));
%             end
%             xticklabels_all{imodel} = xticklabels_;
%         end
%         
%         xticks(1:nmodels)
%         xticklabels(xticklabels_all)
%         
%         ax = gca;
%         ax.YGrid = 'ON';
%         ax.LineWidth=2;
%         ax.TickLength=[0 0];
%         
%         % plotIDVD = 0;
%         % flagLongXTicks = 1;% show the param index in x-axis(1 for free and 0 for fixed)
%         % ==================
%         % RCplot_AICc
%         % ==================
%         ylabel(sprintf('Delta %s',IC_type))
%         
%         set(findall(gcf, '-property', 'FontSize'), 'FontSize',14)
%         sgtitle(sprintf('[%s-M%d-%s] %s vs. %s', namesFeature{ifeature}, imodelKernel, namesType{itype}, namesLocComb{combInd(icomb, [1,2])}), 'fontsize', 25)
%         
%     end % end of itype
% end % end of imodel
% 
% fprintf('Plotting done')
% 
% 
% %%   plot all the fitted params
% %             plotFittedParams_MC
% 
