% 
% % ModelPlot_local
% % plot model predictions for one particular variation for all/one subjects
% 
% %% get median and CI
% if isubj>1
%     for isubj = 1:nsubj
%         [data_F_med, ~, ~, data_F_neg, data_F_pos] = getCI(data_metrics_allSubj_F{isubj}); data_F_med_allSubj(isubj, :) = data_F_med; data_F_neg_allSubj(isubj,:) = data_F_neg; data_F_pos_allSubj(isubj,:) = data_F_pos;
%         [pred_F_med, ~, ~, pred_F_neg, pred_F_pos] = getCI(pred_allB_allSubj_F{isubj}); pred_F_med_allSubj(isubj, :) = pred_F_med; pred_F_neg_allSubj(isubj,:) = pred_F_neg; pred_F_pos_allSubj(isubj,:) = pred_F_pos;
%         [params_F_med, ~, ~, params_F_neg, params_F_pos] = getCI(params_est_allSubj_F{isubj}); params_F_med_allSubj(isubj, :) = params_F_med; params_F_neg_allSubj(isubj,:) = params_F_neg; params_F_pos_allSubj(isubj,:) = params_F_pos;
%         if plotLoc2
%             [data_P_med, ~, ~, data_P_neg, data_P_pos] = getCI(data_metrics_allSubj_P{isubj}); data_P_med_allSubj(isubj, :) = data_P_med; data_P_neg_allSubj(isubj,:) = data_P_neg; data_P_pos_allSubj(isubj,:) = data_P_pos;
%             [pred_P_med, ~, ~, pred_P_neg, pred_P_pos] = getCI(pred_allB_allSubj_P{isubj}); pred_P_med_allSubj(isubj, :) = pred_P_med; pred_P_neg_allSubj(isubj,:) = pred_P_neg; pred_P_pos_allSubj(isubj,:) = pred_P_pos;
%             [params_P_med, ~, ~, params_P_neg, params_P_pos] = getCI(params_est_allSubj_P{isubj}); params_P_med_allSubj(isubj, :) = params_P_med; params_P_neg_allSubj(isubj,:) = params_P_neg; params_P_pos_allSubj(isubj,:) = params_P_pos;
%         end
%     end
% else
%     data_F_med_allSubj = data_metrics_allSubj_F{isubj}; data_F_neg_allSubj = zeros(size(data_metrics_allSubj_F{isubj})); data_F_pos_allSubj = zeros(size(data_metrics_allSubj_F{isubj}));
%     pred_F_med_allSubj = pred_allB_allSubj_F{1}; pred_F_neg_allSubj = zeros(size(data_metrics_allSubj_F{isubj})); pred_F_pos_allSubj = zeros(size(data_metrics_allSubj_F{isubj}));
%     params_F_med_allSubj = params_est_allSubj_F{1}; params_F_neg_allSubj = zeros(size(data_metrics_allSubj_F{isubj})); params_F_pos_allSubj = zeros(size(data_metrics_allSubj_F{isubj}));
%     if plotLoc2
%         data_P_med_allSubj = data_metrics_allSubj_P{isubj}; data_P_neg_allSubj = zeros(size(data_metrics_allSubj_P{isubj})); data_P_pos_allSubj = zeros(size(data_metrics_allSubj_P{isubj}));
%         pred_P_med_allSubj = pred_allB_allSubj_P{1}; pred_P_neg_allSubj = zeros(size(data_metrics_allSubj_P{isubj})); pred_P_pos_allSubj = zeros(size(data_metrics_allSubj_P{isubj}));
%         params_P_med_allSubj = params_est_allSubj_P{1}; params_P_neg_allSubj = zeros(size(data_metrics_allSubj_P{isubj})); params_P_pos_allSubj = zeros(size(data_metrics_allSubj_P{isubj}));
%     end
% end
% 
% %% Fig 1: metrics data vs. pred
% xlim_allM = [.7, 1.6; -1, 1; .5, 1; .5, 1; .1, .5; .5, .9; .5, .9;.5, .9];
% figure('Position', [0 0 1000 600])
% 
% for im = 1:nmetrics
%     subplot(2,4,im), hold on
%     errorbar(data_F_med_allSubj(:, im), pred_F_med_allSubj(:, im), pred_F_neg_allSubj(:, im), pred_F_pos_allSubj(:, im), data_F_neg_allSubj(:, im), data_F_pos_allSubj(:, im), 'r.', 'CapSize', 0)
%     if plotLoc2
%         errorbar(data_P_med_allSubj(:, im), pred_P_med_allSubj(:, im), pred_P_neg_allSubj(:, im), pred_P_pos_allSubj(:, im), data_P_neg_allSubj(:, im), data_P_pos_allSubj(:, im), 'r.', 'CapSize', 0)
%     end
%     plot([xlim_allM(im, 1), xlim_allM(im, 2)], [xlim_allM(im, 1), xlim_allM(im, 2)], 'k-')
%     xlim(xlim_allM(im, :))
%     ylim(xlim_allM(im, :))
%     axis square
%     title(namesMetrics{im})
% end
% set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
% 
% saveas(gcf, 'publishedPDFs/model/n12_Fig1_data_vs_pred.jpg')
% 
% %% Fig 2: pA vs. params
% figure('Position', [0 0 1000 500])
% for ip = 1:nparams_model
%     subplot(1,nparams_model,ip), hold on
%     %     errorbar(data_F_med_allSubj(:, 6), params_F_med_allSubj(:, ip), data_F_neg_allSubj(:, 6), data_F_pos_allSubj(:, 6), params_F_neg_allSubj(:, ip), params_F_pos_allSubj(:, ip), 'r.', 'CapSize', 0)
%     %     errorbar(data_P_med_allSubj(:, 6), params_P_med_allSubj(:, ip), data_P_neg_allSubj(:, 6), data_P_pos_allSubj(:, 6), params_P_neg_allSubj(:, ip), params_P_pos_allSubj(:, ip), 'b.', 'CapSize', 0)
%     if isubj>1
%         for isubj = 1:nsubj
%             plot(data_F_med_allSubj(isubj, 6), params_F_med_allSubj(isubj, ip), ['r', markers_allSubj{isubj}])
%             if plotLoc2
%                 plot(data_P_med_allSubj(isubj, 6), params_P_med_allSubj(isubj, ip), ['b', markers_allSubj{isubj}])
%             end
%         end
%     else
%         plot(data_F_med_allSubj(:, 6), params_F_med_allSubj(:, ip), ['r', markers_allSubj{isubj}])
%         if plotLoc2
%             plot(data_P_med_allSubj(:, 6), params_P_med_allSubj(:, ip), ['b', markers_allSubj{isubj}])
%         end
%     end
%     %     [r, p] = corr([data_F_med_allSubj(:, 6); data_P_med_allSubj(:, 6)], [params_F_med_allSubj(:, ip);params_F_med_allSubj(:, ip)]);
%     [r_F, p_F] = corr(data_F_med_allSubj(:, 6), params_F_med_allSubj(:, ip));
%     if plotLoc2
%         [r_P, p_P] = corr(data_P_med_allSubj(:, 6), params_P_med_allSubj(:, ip));
%     else
%         r_P =nan; p_P=nan;
%     end
%     axis square
%     xlabel('pA')
%     ylabel(namesParamsModel{ip})
%     xlim([.55, .85])
%     title(sprintf('[Fovea] r=%.3f, p=%.3f\n[Peri] r=%.3f, p=%.3f', r_F, p_F, r_P,p_P))
% end
% set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
% saveas(gcf, 'publishedPDFs/model/n12_Fig2_pA_vs_params.jpg')
% 
% %% Fig 3. params, F vs.P
% if plotLoc2
%     figure('Position', [3e3 0 1000 500])
%     for ip = 1:nparams_model
%         xmin = min([params_F_med_allSubj(:, ip); params_P_med_allSubj(:, ip)]);
%         xmax = max([params_F_med_allSubj(:, ip); params_P_med_allSubj(:, ip)]);
%         subplot(1,nparams_model,ip), hold on
%         
%         if isubj>1
%             for isubj = 1:nsubj
%                 errorbar(params_F_med_allSubj(isubj, ip), params_P_med_allSubj(isubj, ip), ...
%                     params_F_neg_allSubj(isubj, ip), params_F_pos_allSubj(isubj, ip), ...
%                     params_P_neg_allSubj(isubj, ip), params_P_pos_allSubj(isubj, ip), '.k', 'Capsize', 0, 'HandleVisibility','off')
%                 plot(params_F_med_allSubj(isubj, ip), params_P_med_allSubj(isubj, ip), ['k', markers_allSubj{isubj}], 'MarkerFaceColor', 'w')
%             end
%             errorbar(mean(params_F_med_allSubj(:, ip)), mean(params_P_med_allSubj(:, ip)), std(params_F_med_allSubj(:, ip))/sqrt(nsubj), 'horizontal', 'm.', 'Capsize', 0, 'linewidth', 2)
%             errorbar(mean(params_F_med_allSubj(:, ip)), mean(params_P_med_allSubj(:, ip)), std(params_P_med_allSubj(:, ip))/sqrt(nsubj), 'vertical', 'm.', 'Capsize', 0, 'linewidth', 2)
%             
%         else
%             errorbar(params_F_med_allSubj(:, ip), params_P_med_allSubj(:, ip), params_F_neg_allSubj(:, ip), params_F_pos_allSubj(:, ip), params_P_neg_allSubj(:, ip), params_P_pos_allSubj(:, ip), '.k', 'Capsize', 0, 'HandleVisibility','off')
%             plot(params_F_med_allSubj(:, ip), params_P_med_allSubj(:, ip), '.k')
%             
%             errorbar(mean(params_F_med_allSubj(:, ip)), mean(params_P_med_allSubj(:, ip)), std(params_F_med_allSubj(:, ip))/sqrt(nsubj), 'horizontal', 'm.', 'Capsize', 0, 'linewidth', 2)
%             errorbar(mean(params_F_med_allSubj(:, ip)), mean(params_P_med_allSubj(:, ip)), std(params_P_med_allSubj(:, ip))/sqrt(nsubj), 'vertical', 'm.', 'Capsize', 0, 'linewidth', 2)
%             
%         end
%         %     if ip == 1, legend(subjList), end
%         plot([xmin, xmax], [xmin, xmax], 'k-')
%         xlim([xmin, xmax])
%         ylim([xmin, xmax])
%         xlabel('Fovea')
%         ylabel('Periphery')
%         axis square
%         
%         [~, p, CI, stats] = ttest(params_F_med_allSubj(:, ip), params_P_med_allSubj(:, ip));
%         title(sprintf('%s\nt(%d) = %.3f, p = %.3f', namesParamsModel{ip}, stats.df, stats.tstat, p))
%     end
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
%     saveas(gcf, 'publishedPDFs/model/n12_Fig3_params_F_vs_P.jpg')
% end