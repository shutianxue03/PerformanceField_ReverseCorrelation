% 
% % close all
% clc
% close all
% 
% % analysize the CV output of all subj
% ifeature = 2;
% ifamily = 4;
% plotIDVDFlag = 0;
% 
% subjList = {'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'YK', 'HL', 'FH', 'HA'};
% nsubj = length(subjList);
% 
% SX_RC1_setting
% sz_label = 30;
% sz_ticks = 25;
% 
% if plotIDVDFlag
%     load(sprintf('Data_OOD/n%d_B1.mat', nsubj)) % to load margK_allSubj
% end
% 
% ylim_ = [0, .015];
% yticks_ = 0:.005:.015;
% yticklabels_ = {'0', '0.005', '0.010', '0.015'};
% 
% 
% %%
% nparams_full = length(namesParams_all{ifamily});
% if (ifamily == 3) || (ifamily == 7) % the truncation term is always left free
%     nmodels = 2^(nparams_full-1);
%     paramInd_all_ = fxn_getParamInd(nparams_full-1);
%     paramInd_all = [paramInd_all_, ones(nmodels,1)];
% else
%     nmodels = 2^nparams_full;
%     paramInd_all = fxn_getParamInd(nparams_full);
% end
% 
% % empty containers
% RSS_allSubj = nan(nsubj, nmodels);
% ibest_allSubj = nan(1, nmodels);
% ibest_mean_record = zeros(1,nmodels);
% ibest_med_record = zeros(1,nmodels);
% freq_allSubj = nan(nsubj, nmodels);
% 
% for isubj = 1:nsubj
%     subjName = subjList{isubj};
%     load(sprintf('Data_MC/CV_%s_%s_M%d.mat', subjName, namesFeature{ifeature}, ifamily))
%     RSS_ave_all = RSS_ML_ave_all;
%     RSS_std_all = RSS_ML_std_all;
%     assert(size(RSS_ave_all, 2) == nmodels)
%     
%     [ncv, nmodels] = size(RSS_ML_ave_all);
%     
%     %% get freq of ibest
%     [min_, ii] = min(RSS_ML_ave_all, [], 2);
%     freq = nan(1,nmodels);
%     for imodel = 1:nmodels
%         freq(imodel) = sum(ii == imodel);
%     end
%     freq = freq/ncv;
%     freq_allSubj(isubj, :) = freq;
%     %% plot the RSS ave and std of each model
%     RSS_med = median(RSS_ave_all);
%     RSS_mean = mean(RSS_ave_all);
%     RSS_std = median(RSS_std_all);
%     
%     RSS_allSubj(isubj, :) = RSS_med; % use the median, which is more stable
%     [~, ibest] = min(RSS_med);
%     ibest_allSubj(isubj) = ibest;
%     [RSS_ave_s ,iRSS_ave_s]= sort(RSS_med);
%     
%     % get the best model given median and mean
%     [~, ibest_med] = min(RSS_med);
%     [~, ibest_mean] = min(RSS_mean);
%     ibest_mean_record(ibest_mean) = ibest_mean_record(ibest_mean)+1;
%     ibest_med_record(ibest_med) = ibest_med_record(ibest_med)+1;
%     fprintf('Best model by mean is %d by median is %d\n', ibest_mean, ibest_med)
%     
%     % plot
%     if plotIDVDFlag == 1
%         figure
%         subplot(1,2,2), hold on
%         bar(1:nmodels, RSS_ave_s, 'FaceColor', 'w', 'EdgeColor', 'k')
%         errorbar(1:nmodels, RSS_ave_s , RSS_std(iRSS_ave_s), 'k.', 'CapSize', 0)
%         xlabel('Model #')
%         ylabel('Mean RSS')
%         xticks(1:nmodels)
%         title(subjName)
%         
%         % x tick labels and text
%         xtl_sort = cell(1,nmodels); % xTickLabels
%         for imodel = 1:nmodels
%             xtl_ = [];
%             for ip = 1:nparams_full
%                 xtl_ = sprintf('%s\\newline%d', xtl_ , paramInd_all(iRSS_ave_s(imodel), ip));
%             end
%             xtl_sort{imodel} = xtl_;
%         end
%         xticks(1:nmodels)
%         xticklabels(xtl_sort)
%         
%         subplot(1,2,1), hold on
%         plot(squeeze(margK_allSubj{ifeature}(isubj, 3, 1, 1, :)))
%         plot(squeeze(margK_allSubj{ifeature}(isubj, 3, 1, 2, :)))
%         
%         set(findall(gcf, '-property', 'FontSize'), 'FontSize',14)
%         sgtitle(subjName)
%     end
% end % end of isubj
% 
% %% freq of each model being selected
% freq_allSubj_ave = mean(freq_allSubj);
% freq_allSubj_sem = std(freq_allSubj)/sqrt(nsubj);
% % ibest_mean_freq = ibest_mean_record/nsubj;
% % ibest_med_freq = ibest_med_record/nsubj;
% 
% %% plot group average (for poster)
% RSS_allSubj_ave = mean(RSS_allSubj);
% RSS_allSubj_std = std(RSS_allSubj)/sqrt(nsubj);
% [RSS_allSubj_ave_s, iRSS_allSubj]= sort(RSS_allSubj_ave);
% 
% figure('Position', [0 0 1000 800])
% % for suppFlag = [0,1] % 0 = for poster; 1=for supplememtary
%     suppFlag = 1;
%     subplot(2,1,1), hold on
%     bar(1:nmodels, RSS_allSubj_ave_s - min(RSS_allSubj_ave_s), 'FaceColor', 'w', 'EdgeColor', 'k', 'barwidth', .5)
%     errorbar(1:nmodels, RSS_allSubj_ave_s - min(RSS_allSubj_ave_s) , RSS_allSubj_std(iRSS_allSubj), 'k.', 'CapSize', 0)
% %     xlabel('Candidate model index', 'fontsize', sz_label)
%     ylabel('10-fold cross-validated $\Delta$MSE', 'fontsize', sz_label, 'interpreter', 'latex')
%     
%     % x tick labels and text
% %     xtl_sort = cell(1,nmodels);
% %     for imodel = 1:nmodels
% %         xtl_ = [];
% %         for ip = 1:nparams_full
% %             xtl_ = sprintf('%s\\newline%d', xtl_ , paramInd_all(iRSS_allSubj(imodel), ip));
% %         end
% %         xtl_sort{imodel} = xtl_;
% %     end
% %     xticks(1:nmodels)
%     
% %     set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',2)
%     
%     ax = gca;
%     ax.XAxis.FontSize = sz_ticks;
%     ax.YAxis.FontSize = sz_ticks;
%     ax.LineWidth = 3.5;
%     
% %     if suppFlag
%         xticklabels([])
%         title(sprintf('%s - %s\nMSE of the best model = %.4f', namesFeature{ifeature}, namesFamily_all{ifamily}, RSS_ave_s(1)), 'fontsize', 20)
%         
%         saveas(gcf,sprintf( 'VSS2022/supp/MC_%s_M%d.jpg', namesFeature{ifeature}, ifamily ))
% %     else
% %         xticklabels(1:nmodels)
% %         ylim(ylim_)
% %         yticks(yticks_)
% %         yticklabels(yticklabels_)
% %         saveas(gcf,sprintf( 'VSS2022/fig/MC_%s_M%d.jpg', namesFeature{ifeature}, ifamily ))
% %     end
% % end
% 
% %% plot the frequency of each model being selected (for VSS supplement)
% % figure('Position', [0 0 1000 500]), hold on
% subplot(2,1,2), hold on
% bar(1:nmodels, freq_allSubj_ave(iRSS_allSubj), 'FaceColor', 'w', 'EdgeColor', 'k', 'barwidth', .5)
% errorbar(1:nmodels, freq_allSubj_ave(iRSS_allSubj), freq_allSubj_sem(iRSS_allSubj), '.k', 'CapSize', 0)
% xtl_sort = cell(1,nmodels);
% for imodel = 1:nmodels
%     xtl_ = [];
%     for ip = 1:nparams_full
%         xtl_ = sprintf('%s\\newline%d', xtl_ , paramInd_all(iRSS_allSubj(imodel), ip));
%     end
%     xtl_sort{imodel} = xtl_;
% end
% xticks(1:nmodels)
% xticklabels(xtl_sort)
% 
% title('Frequency of being selected as the best model')
% 
% set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',2)
% set(findall(gcf, '-property', 'fontsize'), 'fontsize',20)
% 
% saveas(gcf,sprintf( 'VSS2022/supp/MC_%s_M%d.jpg', namesFeature{ifeature}, ifamily ))
% 
