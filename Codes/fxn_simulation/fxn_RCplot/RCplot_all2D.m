% 
% clc
% pCriterion = .001;
% LocFlag = 2; % 1= plot only [fovea+peri], 2= plot that and [peri only]
% 
% %% Kernels
% data2D_ave = getCI(kernels2D, errType, ntrialsProp);
% for itype = 1:ntypes
%     data2D_ave_ = squeeze(data2D_ave(itype, :, :, :));
%     data2D = squeeze(kernels2D(:, itype, :, :, :));
%     
%     for n=1:LocFlag
%         if n == 1, iLocPlot = 1:nLoc5; else, iLocPlot = 2:5;  end
%         % decide the axis range of colorbar
%         data2D_toPlot_ave = data2D_ave_(iLocPlot, :, :);
%         
%         caxisLim = [min(data2D_toPlot_ave(:)), max(data2D_toPlot_ave(:))];
%         figure('Position', [0 0 1000 800])
%         for iLoc = iLocPlot
%             subplot(3,3, subplot_locs(iLoc)), hold on
%             % get data
%             data2D_ave_perLoc = squeeze(data2D_ave_(iLoc , :, :));
%             data2D_perLoc = squeeze(data2D(:, iLoc, :, :));
%             % get outline
%             if nB>1, outline = getOutline(nfiltersOri, nfiltersSF, data2D_perLoc, nB, pCriterion);
%             else, outline = [];
%             end
%             % plot
%             RCplot_2Dkernel(data2D_ave_perLoc', xtlabels_tuning, xticks_contour_tuning, caxisLim, outline')
%             
%             if iLoc == 2, xlabel('ORI (deg)'), ylabel('SF (cpd)'), end
%             title(namesLoc2D{iLoc})
%         end
%         set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
%         make_sgtitle(sprintf('[%s-%s]', title_, namesType{itype}), subjName, nB, nsubj, nAllTrials_allSubj)
%     end
% end
% 
% %%  intercept
% clc
% data2D_ave = getCI(intercept2D, errType, ntrialsProp);
% iLocPlot = 1:nLoc5;
% for itype = 1:ntypes
%     data2D_ave_ = squeeze(data2D_ave(itype, :, :, :));
%     caxisLim = [min(data2D_ave_(:)), max(data2D_ave_(:))];
%     figure('Position', [0 0 1000 800])
%     for iLoc = iLocPlot
%         subplot(3,3, subplot_locs(iLoc)), hold on
%         % get data
%         data2D_ave_perLoc = squeeze(data2D_ave_(iLoc , :, :));
%         data2D_perLoc = squeeze(data2D(:, iLoc , :, :));
%         % get outline
%         if nB>1, outline = getOutline(nfiltersOri, nfiltersSF, data2D_perLoc, nB, pCriterion);
%         else, outline = [];
%         end
%         % plot
%         RCplot_2Dkernel(data2D_ave_perLoc', xtlabels_tuning, xticks_contour_tuning, caxisLim, outline')
% 
%         if iLoc == 2, xlabel('ORI (deg)'), ylabel('SF (cpd)'), end
%         title(namesLoc2D{iLoc})
%     end
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
%     make_sgtitle(sprintf('[%s] intercept', namesType{itype}), subjName, nB, nsubj, nAllTrials_allSubj)
% end
% 
% %%  Tjur's D
% clc
% R2_Tjur_ave = getCI(R2_Tjur, errType, ntrialsProp);
% iLocPlot = 1:nLoc5;
% for itype = 1:ntypes
%     data2D_ave_ = squeeze(R2_Tjur_ave(itype, :, :, :));
%     caxisLim = [min(data2D_ave_(:)), max(data2D_ave_(:))];
%     figure('Position', [0 0 1000 800])
%     for iLoc = iLocPlot
%         subplot(3,3, subplot_locs(iLoc)), hold on
%         % get data
%         data2D_ave_ = squeeze(data2D_ave_(iLoc , :, :));
%         data2D_ = squeeze(data2D(:, iLoc , :, :));
%         % get outline
%         if nB>1, outline = getOutline(nfiltersOri, nfiltersSF, data2D_, nB, pCriterion);
%         else, outline = [];
%         end
%         % plot
%         RCplot_2Dkernel(data2D_ave_', xtlabels_tuning, xticks_contour_tuning, caxisLim, outline')
% 
%         if iLoc == 2, xlabel('ORI (deg)'), ylabel('SF (cpd)'), end
%         title(namesLoc2D{iLoc})
%     end
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
%     make_sgtitle(sprintf('[%s] R2 Tjur [probit regression]', namesType{itype}), subjName, nB, nsubj, nAllTrials_allSubj)
% end
% 
% %% p values (for beta(2), i.e., slope)
% clc
% plotPBorder = 1; % 1=plot the delineation of a certain value on top of imagesc
% beta_ind = 2; % 1=intercept, 2=slope
% pCriterion = 0.001; % the p value below which is considered significant
% pValues_ave = getCI(pValues, errType, ntrialsProp);
% for itype = 1:ntypes
%     data2D_ave_ = squeeze(pValues_ave(itype, : , :, :, beta_ind));
%     caxisLim = [min(data2D_ave_(:)), max(data2D_ave_(:))];
%     figure('Position', [0 0 1000 800])
%     for iLoc = iLocPlot
%         subplot(3,3, subplot_locs(iLoc)), hold on
%         % get data
%         data2D_ave_ = squeeze(data2D_ave_(iLoc , :, :));
%         % get outline
%         if nB>1, outline = data2D_ave_ < pCriterion;
%         else, outline=[];
%         end
%         % plot
%         RCplot_2Dkernel(data2D_ave_', xtlabels_tuning, xticks_contour_tuning, caxisLim, outline')
%         if iLoc == 2, xlabel('ORI (deg)'), ylabel('SF (cpd)'), end
%         title(namesLoc2D{iLoc})
%     end
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
%     make_sgtitle(sprintf('[%s] p value [probit regression]', namesType{itype}),subjName, nB, nsubj, nAllTrials_allSubj)
% end
% 
% %% proportion of correct Categorization in the probit regression
% pCat_ave = getCI(pCat, errType, ntrialsProp);
% iLocPlot = 1:nLoc5;
% for itype = 1:ntypes
%     data2D_ave_ = squeeze(pCat_ave(itype, :, :, :));
% %     caxisLim = [min(data2D_ave(:)), max(data2D_ave(:))];
%     caxisLim = [0,1];
%     figure('Position', [0 0 1000 800])
%     for iLoc = iLocPlot
%         subplot(3,3, subplot_locs(iLoc)), hold on
%         % get data
%         data2D_ave_ = squeeze(data2D_ave_(iLoc , :, :));
%         data2D_ = squeeze(data2D(:, iLoc , :, :));
%         % get outline
%         if nB>1, outline = getOutline(nfiltersOri, nfiltersSF, data2D_, nB, pCriterion);
%         else, outline = [];
%         end
%         % plot
%         RCplot_2Dkernel(data2D_ave_', xtlabels_tuning, xticks_contour_tuning, caxisLim, outline')
% 
%         if iLoc == 2, xlabel('ORI (deg)'), ylabel('SF (cpd)'), end
%         title(namesLoc2D{iLoc})
%     end
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
%     make_sgtitle(sprintf('[%s] Prop. correct categorization [probit regression]', namesType{itype}), subjName, nB, nsubj, nAllTrials_allSubj)
% end
% 
% %% plot 2D kernels and stats of combined location
% clc
% pCriterion = .001;
% for itype = 1:ntypes
%     figure('Position', [0 0 1600 1200])
%     for icomb = 1:ncomb4
%         data2D1_ave = squeeze(data2D_ave(itype, combInd(icomb, 1), :, :));
%         data2D2_ave = squeeze(data2D_ave(itype, combInd(icomb, 2), :, :));
%         
%         data2D1 = squeeze(data2D(:, itype, combInd(icomb, 1), :, :));
%         data2D2 = squeeze(data2D(:, itype, combInd(icomb, 2), :, :));
%         
%         % 2D kernels of Fov/HM/LVM/Left
%         if nB>1, outline = getOutline(nfiltersOri, nfiltersSF, data2D1, nB, pCriterion);
%         else, outline = [];
%         end
%         caxisLim = [min(data2D1_ave(:)), max(data2D1_ave(:))];
%         subplot(3,ncomb4, icomb), hold on
%         RCplot_2Dkernel(data2D1_ave', xtlabels_tuning, xticks_contour_tuning, caxisLim, outline')
%         title(namesLocComb{combInd(icomb, 1)})
%         
%         % 2D kernels of Peri/VM/UVM/Right
%         caxisLim = [min(data2D2_ave(:)), max(data2D2_ave(:))];
%         if nB>1, outline = getOutline(nfiltersOri, nfiltersSF, data2D2, nB, pCriterion);
%         else, outline = [];
%         end
%         subplot(3,ncomb4, icomb+ncomb4), hold on
%         RCplot_2Dkernel(data2D2_ave', xtlabels_tuning, xticks_contour_tuning, caxisLim, outline')
%         title(namesLocComb{combInd(icomb, 2)})
%         
%         % diff
%         data2D1_allB = squeeze(data2D(:, itype, combInd(icomb, 1), :, :));
%         data2D2_allB = squeeze(data2D(:, itype, combInd(icomb, 2), :, :));
%         data2D_diff_allB = data2D1_allB - data2D2_allB;
%         if nB>1, data2D_diff_ave = squeeze(mean(data2D_diff_allB, 1));
%         else, data2D_diff_ave = data2D_diff_allB;
%         end
%         caxisLim = [min(data2D_diff_ave(:)), max(data2D_diff_ave(:))];
%         % get outline
%         if nB>1, outline = getOutline(nfiltersOri, nfiltersSF, data2D_diff_allB, nB, pCriterion);
%         else, outline = [];
%         end
%         subplot(3,ncomb4, icomb+ncomb4*2), hold on
%         RCplot_2Dkernel(data2D_diff_ave', xtlabels_tuning, xticks_contour_tuning, caxisLim, outline')
%         title(sprintf('%s minus %s', namesLocComb{combInd(icomb, 1)}, namesLocComb{combInd(icomb, 2)}))
%         
%     end
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize',14)
%     sgtitle(namesType{itype}, 'FontSize',20)
% end
% 
% %% Tjur R2 of the combined loci
% clc
% for itype = 1:ntypes
%     figure('Position', [0 0 1600 800])
%     for icomb = 1:ncomb4
%         data2D1_ave = squeeze(R2_Tjur_ave(itype, combInd(icomb, 1) , :, :));
%         data2D2_ave = squeeze(R2_Tjur_ave(itype, combInd(icomb, 2) , :, :));
%         
%         data2D1 = squeeze(R2_Tjur(:, itype, combInd(icomb, 1) , :, :));
%         data2D2 = squeeze(R2_Tjur(:, itype, combInd(icomb, 2) , :, :));
%         
%         % 2D kernels of Fov/HM/LVM/Left
%         if nB>1, outline = getOutline(nfiltersOri, nfiltersSF, data2D1, nB, pCriterion);
%         else, outline = [];
%         end
%         caxisLim = round([min(data2D1_ave(:)), max(data2D1_ave(:))],2);
%         subplot(2,ncomb4, icomb), hold on
%         RCplot_2Dkernel(data2D1_ave', xtlabels_tuning, xticks_contour_tuning, caxisLim, outline')
%         title(namesLocComb{combInd(icomb, 1)})
%         
%         % 2D kernels of Peri/VM/UVM/Right
%         caxisLim = round([min(data2D2_ave(:)), max(data2D2_ave(:))],2);
%         if nB>1, outline = getOutline(nfiltersOri, nfiltersSF, data2D2, nB, pCriterion);
%         else, outline = [];
%         end
%         subplot(2, ncomb4, icomb+ncomb4), hold on
%         RCplot_2Dkernel(data2D2_ave', xtlabels_tuning, xticks_contour_tuning, caxisLim, outline')
%         title(namesLocComb{combInd(icomb, 2)})
%         
%     end
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize',14)
%     sgtitle(sprintf('[%s] Tjur R2 of combined location', namesType{itype}), 'FontSize',20)
% end
% 
% %% p values of the combined loci
% clc
% for itype = 1:ntypes
%     figure('Position', [0 0 1600 800])
%     for icomb = 1:ncomb4
%         data2D1_ave = squeeze(pValues_ave(itype, combInd(icomb, 1) , :, :, beta_ind));
%         data2D2_ave = squeeze(pValues_ave(itype, combInd(icomb, 2) , :, :, beta_ind));
%         
%         % 2D kernels of Fov/HM/LVM/Left
%         if nB>1, outline = data2D1_ave < pCriterion;
%         else, outline = [];
%         end
%         caxisLim = [min(data2D1_ave(:)), max(data2D1_ave(:))];
%         subplot(2,ncomb4, icomb), hold on
%         RCplot_2Dkernel(data2D1_ave', xtlabels_tuning, xticks_contour_tuning, caxisLim, outline')
%         title(namesLocComb{combInd(icomb, 1)})
%         
%         % 2D kernels of Peri/VM/UVM/Right
%         caxisLim = [min(data2D2_ave(:)), max(data2D2_ave(:))];
%         if nB>1, outline = data2D2_ave < pCriterion;
%         else, outline = [];
%         end
%         subplot(2, ncomb4, icomb+ncomb4), hold on
%         RCplot_2Dkernel(data2D2_ave', xtlabels_tuning, xticks_contour_tuning, caxisLim, outline')
%         title(namesLocComb{combInd(icomb, 2)})
%         
%     end
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize',14)
%     sgtitle(sprintf('[%s] p values of combined location', namesType{itype}), 'FontSize',20)
% end
