% 
% 
% 
% buffer = .05; % the distance between two R2 text
% nmodels_s2 = 8; % the number of models to be saved fo their IC values
% 
% for icomb = 1%:ncomb
%     for itype = 1:ntypes
%         
%         %%%%%%%%%
%         % extract data %
%         %%%%%%%%%
%         IC_allSubj = IC_allSubj_all{icomb, itype};
%         yData_allSubj = yData_allSubj_all{icomb, itype};
%         yPred_allSubj = yPred_allSubj_all{icomb, itype};
%         
%         %%%%%%%%%%%%%%%%%%%%%%%
%         % select the 8 best models based on ave  %
%         %%%%%%%%%%%%%%%%%%%%%%%
%         % get ave and std
%         IC_AVE = mean(IC_allSubj,1);
%         IC_SEM = std(IC_allSubj,[],1)/sqrt(nsubj);
%         % sort
%         [IC_AVE_sort, iIC_AVE_sort] = sort(IC_AVE);
%         IC_SEM_sort = IC_SEM(iIC_AVE_sort);
%         IC_allSubj_sort = IC_allSubj(:, iIC_AVE_sort);
%         % save (for cross-family comparison)
%         IC_AVE_sort_all{imodel_, iIC, icomb, itype} = IC_AVE_sort;
%         iIC_AVE_sort_all{imodel_, iIC, icomb, itype} = iIC_AVE_sort;
%         IC_SEM_sort_all{imodel_, iIC, icomb, itype} = IC_SEM_sort;
%         IC_allSubj_sort_all{imodel_, iIC, icomb, itype} = IC_allSubj_sort;
%         paramInd_all_sort_all{imodel_, iIC, icomb} = paramInd_all(iIC_AVE_sort, :);
%         
%         if plotFlag
%             data_toPlot = [yData_allSubj(:); yPred_allSubj(:)];
%             ymin_ = min(data_toPlot);
%             ymax_ = max(data_toPlot);
%             
%             % do NOT print the figure in the monitor, otherwise the fig
%             % would be too narrow
%             figure('Position', [0 0 nsubj*250 nmodels_s*200])
%             
%             for isubj = 1:nsubj
%                 %%%%%%%%%%
%                 %  extract data  %
%                 %%%%%%%%%%
%                 y_perSubj = squeeze(yData_allSubj(isubj, :, :));
%                 yPred_allModels = squeeze(yPred_allSubj(isubj, :, :, :));
%                 IC_allModels = IC_allSubj(isubj, :);
%                 [IC_sort, iIC_sort] = sort(IC_allModels);
%                 
%                 for imodel_s = 1:nmodels_s
%                     subplot(nmodels_s+1, nsubj, (imodel_s-1)*nsubj+isubj), hold on
%                     if ifeature == 1
%                         xline(0,'k');
%                         xticks([-80,0,80])
%                         xaxis = xaxis_itp;
%                     else
%                         xline(1,'k');
%                         xticks([0,1,2]), xticklabels([1,2,4])
%                         xaxis = log2(xaxis_itp);
%                     end
%                     
%                     for n=1:2
%                         yData_toPlot = y_perSubj(n, :);
%                         yPred_toPlot = squeeze(yPred_allModels(iIC_sort(imodel_s), n, :));
%                         plot(xaxis, yData_toPlot, 'o', 'MarkerSize', 3, 'color', colors_comb(combInd(icomb, n), :))
%                         plot(xaxis, yPred_toPlot, '-', 'color', colors_comb(combInd(icomb, n), :))
%                         % R2
%                         R2 = getR2(yData_toPlot, yPred_toPlot');
%                         text(xaxis(2), ymax_ - buffer*n, sprintf('%d%%', round(R2*100)), 'color', colors_comb(combInd(icomb, n), :))
%                     end
%                     % ylim
%                     ylim([ymin_, ymax_])
%                     % extra line
%                     yline(0, 'k--');
%                     % ylabel
%                     if imodel_s == 1, ylabel(subjList{isubj}), end
%                     % title
%                     title(sprintf('%s (%d)', num2str(paramInd_all(iIC_sort(imodel_s), :)), round(IC_sort(imodel_s))))
%                     % xticks and yticks
%                     if imodel_s > 1
%                         xticks([])
%                         yticks([])
%                     end
%                 end % end of imodel_s
%                 
%             end % end of isubj
%             
%             %%%%%%%
%             % . plot IC . %
%             %%%%%%%
%             % get the min
%             IC_AVE_min = min(IC_AVE_sort);
%             
%             % plot average
%             subplot(nmodels_s+1, nsubj, nmodels_s*nsubj+3:(nmodels_s+1)*nsubj-3), hold on
%             bar(1:nmodels_s2, IC_AVE_sort(1:nmodels_s2) - IC_AVE_min, 'barwidth', .5, 'edgecolor', colorsType{itype}, 'facecolor', 'w')
%             errorbar(1:nmodels_s2, IC_AVE_sort(1:nmodels_s2) - IC_AVE_min, IC_SEM_sort(1:nmodels_s2), '.k', 'CapSize', 0)
%             
%             % plot IDVD data
%             for isubj = 1:nsubj
%                 plot((1:nmodels_s2)-.2, IC_allSubj_sort(isubj,1:nmodels_s2) - IC_AVE_min, marks_allSubj{isubj}, 'color', ones(1,3)/2)
%             end
%             
%             % x tick labels and text
%             for imodel = 1:nmodels_s2
%                 text(imodel, IC_AVE_sort(imodel) - IC_AVE_min + 10, sprintf('%.2f', IC_AVE_sort(imodel) - IC_AVE_min), 'fontsize', 10)
%                 xticklabels_ = [];
%                 for ip = 1:nparams_full
%                     xticklabels_ = sprintf('%s\\newline%d', xticklabels_ , paramInd_all(iIC_AVE_sort(imodel), ip));
%                 end
%                 xticklabels_all{imodel} = xticklabels_;
%             end
%             
%             xticks(1:nmodels)
%             xticklabels(xticklabels_all)
%             
%             ylabel(['delta ' IC_type])
%             ax = gca;
%             ax.YGrid = 'ON';
%             ax.LineWidth=2;
%             ax.TickLength=[0 0];
%             
%             %%%
%             set(findall(gcf, '-property', 'FontSize'), 'FontSize',14)
%             sgtitle(sprintf('%s-M%d\n[%s] %s vs. %s', namesFeature{ifeature}, imodelKernel, namesType{itype}, namesLocComb{combInd(icomb, [1,2])}), 'fontsize', 25)
%         end
%     end
% end
% 
% 
