% 
% % plot the 2D kernel (kernels2D) and the variance explained from GLM (var_exp)
% 
% %% set parans
% clear xticks_
% xticks_{1} = filtersSF_all_log_cut([1,round(nfiltersSF_cut/2),nfiltersSF_cut]);
% xticks_{2} = [-80, -40, 0, 40, 80];
% xticklabels_{1} = round(filtersSF_all_cut([1,round(nfiltersSF_cut/2),nfiltersSF_cut]),2);
% xticklabels_{2} = [-80, -40, 0, 40, 80];
% 
% xlabels = {'SF (cpd)', 'orientation (deg)'};
% 
% sz_ticks = 15;
% sz_label = 20;
% sz_title = 20;
% sz_text = 14;
% 
% %% WEIGHT the kernel by number of trials
% kernels2D_ave_mirror = 0;
% if nsubj>1
%     for isubj = 1:nsubj
%         kernels2D_ave_mirror = kernels2D_ave_mirror + squeeze(kernels2D_allSubj_mirror(isubj, :, :, :, :)) * ntrialsProp(isubj);
%     end
% else, kernels2D_ave_mirror = squeeze(kernels2D_allSubj_mirror(isubj, :, :, :, :));
% end
% 
% %%
% for itype = 1:ntypes
%     for n=1:2 % 1=fovea+peri, 2=peri only
%         if n == 1, iLocPlot = 1:nLoc; else, iLocPlot = 2:5;  end
%         
%         % decide the axis range of
%         kk = squeeze(kernels2D_ave_mirror(itype, iLocPlot, :, :));
%         caxisMin = min(kk(:));
%         caxisMax = max(kk(:));
%         
%         figure('Position', [0 0 1000 800])
%         for iLoc = iLocPlot
%             subplot(3,3, subplot_locs(iLoc)), hold on
%             
%             imagesc(filtersOri_all-90, flip(filtersSF_all_log_cut), flip(squeeze(kernels2D_ave_mirror(itype, iLoc , :, :))'))
%             xline(0, 'r-', 'linewidth', 2);
%             yline(1, 'r-', 'linewidth', 2);
%             cc = colorbar; if iLoc == 3, cc.Label.String = 'kernel (a.u.)'; end
%             caxis([caxisMin, caxisMax])
%             axis square
%             
%             if iLoc == 2, xlabel(xlabels{2}), ylabel(xlabels{1}), end
%             xticks(xticks_{2})
%             xticklabels(xticklabels_{2})
%             xlim(xticks_{2}([1,end]))
%             
%             yticks(xticks_{1}) % yticks vector must be increasing
%             yticklabels(xticklabels_{1})
%             ylim(xticks_{1}([1,end]))
%             
%             title(locNames{iLoc})
%             
%         end
%         set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
%         if ~publishFlag
%             if nsubj==1, saveas(gcf, sprintf('publishedPDFs/fig/%s_2D_%s.jpg', subjName, typeNames{itype}))
%             else, saveas(gcf, sprintf('publishedPDFs/fig/n%d_2D_%s.jpg', nsubj, typeNames{itype}))
%             end
%         else
%             if nsubj == 1, sgtitle(sprintf('[%s] %s: %d trials per loc', typeNames{itype}, subjName, nAllTrials), 'FontSize',25)
%             else, sgtitle(sprintf('[%s] n = %d (%d +- %d)',  typeNames{itype}, nsubj, round(mean(nAllTrials)), round(std(nAllTrials)/sqrt(nsubj))), 'FontSize',25)
%             end
%         end
%         saveas(gcf, sprintf('publishedPDFs/fig/fig_%s%d.jpg', typeNames{itype}, n))
%     end
% end
% 
% %% plot R2_Tjur
% if nsubj==1
%     for itype = 1:ntypes
%         R2_Tjur__ = R2_Tjur(itype, :, :, :);
%         caxisMax = max(R2_Tjur__(:));
%         caxisMin = min(R2_Tjur__(:));
%         figure('Position', [0 0 1000 800])
%         for iLoc = 1:nLoc
%             
%             R2_Tjur_ = squeeze(R2_Tjur(itype, iLoc , :, :));
%             subplot(3,3, subplot_locs(iLoc)), hold on
%             
%             imagesc(filtersOri_all-90, flip(filtersSF_all_log_cut), flip(R2_Tjur_'))
%             xline(0, 'r-', 'linewidth', 2);
%             yline(1, 'r-', 'linewidth', 2);
%             colorbar
%             caxis([caxisMin, caxisMax])
%             axis square
%             
%             if iLoc == 2, xlabel(xlabels{2}), ylabel(xlabels{1}), end
%             xticks(xticks_{2})
%             xticklabels(xticklabels_{2})
%             xlim(xticks_{2}([1,end]))
%             
%             yticks(xticks_{1}) % yticks vector must be increasing
%             yticklabels(xticklabels_{1})
%             ylim(xticks_{1}([1,end]))
%             
%             title(locNames{iLoc})
%             
%         end
%         
%         set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
%         sgtitle(sprintf('[%s] R2 Tjur\n%s: %d trials per loc', typeNames{itype}, subjName, nAllTrials), 'FontSize',25)
%     end
% end
% 
% %% plot p values (for beta(2), i.e., slope)
% beta_ind =2; % 1=intercept, 2=slope
% pCriterion = 0.001; % the p value below which is considered significant
% if nsubj==1
%     for itype = 1:ntypes
%         figure('Position', [0 0 1000 800])
%         for iLoc = 1:nLoc
%             
%             pValues_ = squeeze(pValues(itype, iLoc , :, :, beta_ind));
%             subplot(3,3, subplot_locs(iLoc)), hold on
%             
%             % plot the p values
% %             imagesc(filtersOri_all-90, flip(filtersSF_all_log_cut), flip(pValues_'))
%             % plot significant results
%             imagesc(filtersOri_all-90, flip(filtersSF_all_log_cut), flip(pValues_') <= pCriterion)
%             xline(0, 'r-', 'linewidth', 2);
%             yline(1, 'r-', 'linewidth', 2);
%             colorbar
%             caxis([0, .002])
%             axis square
%             
%             if iLoc == 2, xlabel(xlabels{2}), ylabel(xlabels{1}), end
%             xticks(xticks_{2})
%             xticklabels(xticklabels_{2})
%             xlim(xticks_{2}([1,end]))
%             
%             yticks(xticks_{1}) % yticks vector must be increasing
%             yticklabels(xticklabels_{1})
%             ylim(xticks_{1}([1,end]))
%             
%             title(locNames{iLoc})
%             
%         end
%         
%         set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
%         sgtitle(sprintf('[%s] p value of slope estimation\n%s: %d trials per loc', typeNames{itype}, subjName, nAllTrials), 'FontSize',25)
%     end
% end
% 
% %% combine locations (each cell: nsubj x ntypes x nfiltersOri x nfiltersSF)
% kernels2D_allSubj_comp = cell(length(locNames_comp), 1);
% kernelsOri_allSubj_comp = cell(length(locNames_comp), 1);
% kernelsSF_allSubj_comp = cell(length(locNames_comp), 1);
% 
% % each single loc (1-5)
% for iLoc = 1:nLoc
%     kernels2D_allSubj_comp{iLoc} = squeeze(kernels2D_allSubj_mirror(:, :, iLoc, :, :));
%     kernelsOri_allSubj_comp{iLoc} = squeeze(mean(kernels2D_allSubj_mirror(:, :, iLoc, :, :), 5));
%     kernelsSF_allSubj_comp{iLoc} = squeeze(mean(kernels2D_allSubj_mirror(:, :, iLoc, :, :), 4));
% end
% 
% % HM (6)
% kernels2D_allSubj_comp{6} = squeeze(mean(kernels2D_allSubj_mirror(:, :, [2,4], :, :), 3));
% kernelsOri_allSubj_comp{6} = squeeze(mean(mean(kernels2D_allSubj_mirror(:, :, [2,4], :, :), 3), 5));
% kernelsSF_allSubj_comp{6} = squeeze(mean(mean(kernels2D_allSubj_mirror(:, :, [2,4], :, :), 3), 4));
% 
% % VM (7)
% kernels2D_allSubj_comp{7} = squeeze(mean(kernels2D_allSubj_mirror(:, :, [3,5], :, :), 3));
% kernelsOri_allSubj_comp{7} = squeeze(mean(mean(kernels2D_allSubj_mirror(:, :, [3,5], :, :), 3), 5));
% kernelsSF_allSubj_comp{7} = squeeze(mean(mean(kernels2D_allSubj_mirror(:, :, [3,5], :, :), 3), 4));
% 
% % periphery (8)
% kernels2D_allSubj_comp{8} = squeeze(mean(kernels2D_allSubj_mirror(:, :, 2:5, :, :), 3));
% kernelsOri_allSubj_comp{8} = squeeze(mean(mean(kernels2D_allSubj_mirror(:, :, 2:5, :, :), 3), 5));
% kernelsSF_allSubj_comp{8} = squeeze(mean(mean(kernels2D_allSubj_mirror(:, :, 2:5, :, :), 3), 4));
% 
