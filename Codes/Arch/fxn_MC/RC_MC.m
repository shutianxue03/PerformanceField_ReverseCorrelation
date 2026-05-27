% clc
% warning('off')
% 
% plotFlag = 0;
% publishFlag = 0;
% flagML = 0; % 0=do NOT do MLE, thus no nLL, AIC_ML or BIC calculated
% filesGROUPName = sprintf('Data_OOD/n%d', nsubj);
% if ~exist('margK_allSubj', 'var'), load(sprintf('Data_OOD/n%d_raw.mat', nsubj), 'margK_allSubj', 'ntrials_allSubj'), end
% 
% if publishFlag
%     publishOptionsMC = publishOptions;
%     publishOptionsMC.outputDir = 'publishedPDFs/MC/';
% end
% 
% nmodels_s = 5; % the number of best models to be displayed for each subj
% IC_type_all = {'AIC_LS' 'AIC_ML', 'AICc', 'BIC'};
% nIC = length(IC_type_all); % AIC and AICc
% 
% ifeature = 2;%input('Which feature to fit (1 = ORI, 2 = SF): \n');
% if ifeature == 1, ifamily_all = [1,8]; else, ifamily_all = 3; end
% nfamily_all = length(ifamily_all);
% 
% % empty containers (sorted according to model selection)
% IC_AVE_sort_all = cell(nfamily_all, nIC, ncomb, ntypes);
% iIC_AVE_sort_all = IC_AVE_sort_all;
% IC_SEM_sort_all = IC_AVE_sort_all;
% IC_allSubj_sort_all = IC_AVE_sort_all;
% paramInd_all_sort_all = cell(nfamily_all, nIC, ncomb);
% 
% for iff = 1:nfamily_all
%     ifamily = ifamily_all(iff);
%     
%     %%%%%%%%%
%     %     Fitting       %
%     %%%%%%%%%
%     MCFileName = sprintf('Data/MC_n%d_IDVD_%s_M%d.mat', nsubj , namesFeature{ifeature}, ifamily);
%     MCFileDir = dir(MCFileName);
%     
%     if isempty(MCFileDir)
%         %         MC_4comb %  compare 2 loci within each of 4 pairs, fit all observers together
%         MC_4comb_IDVD % compare 2 loci within each of 4 pairs, fit all observers individually
%         %         MC_4comb_CV % cross validation
%         %         MC_5Loc % compare 5 loc (under construction)
%         
%     else
%         load(MCFileName)
%         fprintf('MC file loaded: %s - [M%d] "%s"\n', namesFeature{ifeature}, ifamily, namesFamily_all{ifamily_all(iff)})
%     end
%     
%     %%%%%%%%%
%     %  PLOTTING  %
%     %%%%%%%%%
%     IC_allSubj_all_allIC = {AIC_LS_allSubj_all, AIC_ML_allSubj_all, AICc_allSubj_all, BIC_allSubj_all};
%     
%     for iIC = 1:nIC % AIC,AICc, BIC
%         if iIC == 3
%             yPred_allSubj_all = yPred_ML_allSubj_all;
%         else
%             yPred_allSubj_all = yPred_LS_allSubj_all;
%         end
%         
%         IC_allSubj_all = IC_allSubj_all_allIC{iIC};
%         IC_type = IC_type_all{iIC};
%         
%         if publishFlag
%             publish('MCplot_4comb_IDVD', publishOptionsMC);
%             movefile([publishOptionsMC.outputDir,'MCplot_4comb_IDVD.pdf'], sprintf('%sMC_n%d_%s_M%d_%s.pdf',  publishOptionsMC.outputDir, nsubj, namesFeature{ifeature}, ifamily, IC_type));
%             close all
%             fprintf('Publish DONE - %s\n', IC_type)
%         else
%             MCplot_4comb_IDVD
%         end
%     end % end of iIC
% end % end of imodel_
% 
% save(MCFileName, '*sort_all', '-append')
% 
% %% compare diff families
% nmodels_s3 = 3; % number of models to be shown eventually
% % must be ODD numbers!!
% barwid = .3;
% 
% for iIC = 1:nIC % AIC,AICc, BIC
%     figure('Position', [0 0 2000 1000])
%     for icomb = 1:ncomb
%         for itype = 1:ntypes
%             
%             % empty containers
%             y_all = nan(nfamily_all, nmodels_s3);
%             yerr_all = y_all;
%             xbar_all = y_all;
%             
%             subplot(ntypes, ncomb, (itype-1)*ncomb+icomb), hold on
%             
%             for iff = 1:nfamily_all
%                 ifamily = ifamily_all(iff);
%                 y_all(iff, :) = IC_AVE_sort_all{iff, iIC, icomb, itype}(1:nmodels_s3);
%                 
%                 %                 iIC_AVE_sort_all %= IC_AVE_sort_all;
%                 %                 IC_allSubj_sort_all %= IC_AVE_sort_all;
%             end
%             
%             y_all_norm = y_all-min(y_all(:));
%             
%             xticklabels_all = cell(nmodels_s3, nfamily_all);
%             for  iff = 1:nfamily_all
%                 yerr = IC_SEM_sort_all{iff, iIC, icomb, itype}(1:nmodels_s3);
%                 xbar = (iff - barwid*floor(nmodels_s3/2)) : barwid : (iff + barwid*floor(nmodels_s3/2));
%                 xbar_all(iff, :) = xbar;
%                 paramInd_all_sort = paramInd_all_sort_all{iff, iIC, icomb};
%                 nparams_full = size(paramInd_all_sort,2);
%                 
%                 for imodels_s3 = 1:nmodels_s3
%                     xticklabels_ = [];
%                     for ip = 1:nparams_full
%                         xticklabels_ = sprintf('%s\\newline%d', xticklabels_ , paramInd_all_sort(imodels_s3, ip));
%                     end
%                     xticklabels_all{imodels_s3, iff} = xticklabels_;
%                 end
%                 bar(xbar, y_all_norm(iff, :), 'barwidth', barwid*2.5)
%                 errorbar(xbar, y_all_norm(iff, :), yerr, '.k', 'Capsize', 0, 'handlevisibility', 'off')
%                 
%                 % highlight the optimal model with red thick error bars
%                 for ii = 1:nmodels_s3
%                     if abs(y_all_norm(iff, ii)) < eps, errorbar(xbar(ii), y_all_norm(iff, ii), yerr(ii), '.r', 'linewidth', 2, 'Capsize', 0, 'handlevisibility', 'off'), end
%                 end
%             end
%             
%             % xticks
%             xbar_all_ = xbar_all';
%             xticks(xbar_all_(:))
%             xticklabels(xticklabels_all)
%             
%             if itype+icomb == 2
%                 legend(namesModel_all(ifamily_all), 'Location', 'best')
%             end
%             title(sprintf('[%s] %s vs. %s', namesType{itype}, namesLocComb{combInd(icomb, [1,2])}))
%             
%         end % end of itype
%     end % end of icomb
%     
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize',14)
%     sgtitle(sprintf('%s', IC_type_all{iIC}), 'FontSize',20)
%     saveas(gcf, sprintf('publishedPDFs/MC/MC_n%d_%s_%s.jpg', nsubj, namesFeature{ifeature}, IC_type_all{iIC}))
%     
% end %end of iIC
% 

