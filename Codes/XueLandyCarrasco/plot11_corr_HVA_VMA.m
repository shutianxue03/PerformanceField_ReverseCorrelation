% 
% sz_marker = 20;
% 
% load(sprintf('n11_B%d_CS_L67', nB), 'cs_allSubj_')
% HVA_cs = (cs_allSubj_(:, :, 1) - cs_allSubj_(:, :, 2))./(cs_allSubj_(:, :, 1)+cs_allSubj_(:, :, 2));
% 
% load(sprintf('n11_B%d_CS_L53', nB), 'cs_allSubj_')
% VMA_cs = (cs_allSubj_(:, :, 1) - cs_allSubj_(:, :, 2))./(cs_allSubj_(:, :, 1)+cs_allSubj_(:, :, 2));
% 
% load(sprintf('n11_B%d_L67_%d_%d_m_ORI%dSF%d', nB, nORI, nSF, ifamily_perF), 'margTuningC*')
% HVA_ORI_allTunC = (margTuningC_ORI_allSubj(:, :, 1, :, :) - margTuningC_ORI_allSubj(:, :, 2, :, :))./(margTuningC_ORI_allSubj(:, :, 1, :, :) + margTuningC_ORI_allSubj(:, :, 2, :, :));
% HVA_SF_allTunC = (margTuningC_SF_allSubj(:, :, 1, :, :) - margTuningC_SF_allSubj(:, :, 2, :, :))./(margTuningC_SF_allSubj(:, :, 1, :, :) + margTuningC_SF_allSubj(:, :, 2, :, :));
% nTunC_ORI = size(HVA_ORI_allTunC, 5);
% nTunC_SF = size(HVA_SF_allTunC, 5);
% 
% load(sprintf('n11_B%d_L53_%d_%d_m_ORI%dSF%d', nB, nORI, nSF, ifamily_perF), 'margTuningC*')
% VMA_ORI_allTunC = (margTuningC_ORI_allSubj(:, :, 1, :, :) - margTuningC_ORI_allSubj(:, :, 2, :, :))./(margTuningC_ORI_allSubj(:, :, 1, :, :) + margTuningC_ORI_allSubj(:, :, 2, :, :));
% VMA_SF_allTunC = (margTuningC_SF_allSubj(:, :, 1, :, :) - margTuningC_SF_allSubj(:, :, 2, :, :))./(margTuningC_SF_allSubj(:, :, 1, :, :) + margTuningC_SF_allSubj(:, :, 2, :, :));
% 
% load('model_n11_n1000_L67_N_conv1_IV1', '*params_med*')
% HVA_NOM = squeeze((NOM_params_med_allSubj(:, 1, :) - NOM_params_med_allSubj(:, 2, :))./(NOM_params_med_allSubj(:, 1, :) + NOM_params_med_allSubj(:, 2, :)));
% 
% load('model_n11_n1000_L53_N_conv1_IV1', '*params_med*')
% VMA_NOM = squeeze((NOM_params_med_allSubj(:, 1, :) - NOM_params_med_allSubj(:, 2, :))./(NOM_params_med_allSubj(:, 1, :) + NOM_params_med_allSubj(:, 2, :)));
% 
% %% select subj and get med
% itype = 2;
% HVA_cs_med = getCI(HVA_cs(indSubj, :), 1, 2);
% VMA_cs_med = getCI(VMA_cs(indSubj, :), 1, 2);
% 
% HVA_ORI_med = getCI(HVA_ORI_allTunC(indSubj, :, :, itype, :), 1, 2);
% HVA_SF_med = getCI(HVA_SF_allTunC(indSubj, :, :, itype, :), 1, 2);
% 
% VMA_ORI_med = getCI(VMA_ORI_allTunC(indSubj, :, :, itype, :), 1, 2);
% VMA_SF_med = getCI(VMA_SF_allTunC(indSubj, :, :, itype, :), 1, 2);
% 
% folderName = sprintf('%s/%s/corr_HVA_VMA/%s/%s/', nameFigFolder, name_numFilters_Fitting, namesType{itype});
% folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
% 
% %% compare HVA and VMA
% sz_marker_idvd = 20;
% wd_bar = .5; % the width of the bar (not the bar edge!!)
% wd_barEdge = 5;
% wd_idvd = 2;
% wd_ref = 3;
% sz_title = 10;
% sz_ticks = 50;
% 
% figure('Position', [0 0 6e2 7e2]), hold on
% [HVA_cs_ave, ~, ~, HVA_cs_sem] = getCI(HVA_cs_med, 2, 1);
% [VMA_cs_ave, ~, ~, VMA_cs_sem] = getCI(VMA_cs_med, 2, 1);
% bar(1:2, [HVA_cs_ave, VMA_cs_ave], 'BarWidth', wd_bar, 'linewidth', wd_barEdge, 'FaceColor', 'w')
% errorbar(1:2, [HVA_cs_ave, VMA_cs_ave], [HVA_cs_sem, VMA_cs_sem], '.k','CapSize', 0, 'linewidth', wd_barEdge)
% for isubj=1:nsubj
%     plot([1.35, 1.65], [HVA_cs_med(isubj), VMA_cs_med(isubj)],  '-', 'color', ones(1,3)*.7, 'markerfacecolor', 'w', 'linewidth', wd_idvd, 'markersize', sz_marker_idvd)
% end
% xticks([1,2]), xticklabels({'HVA', 'VMA'})
% [~, p, ~, stats] = ttest(HVA_cs_med, VMA_cs_med);
% title(sprintf('t=%.2f, p=%.3f (%d/%d)', p, stats.tstat, sum(HVA_cs_med>=VMA_cs_med), nsubj))
% ax = gca;
% ax.XAxis.FontSize = sz_ticks;
% ax.YAxis.FontSize = sz_ticks;
% ax.LineWidth = wd_barEdge;
% yticks([0, .1, .2])
% get_ydiff = @(lb, ub) ub-(ub-lb)/6;
% diff_y = get_ydiff(0, .2);
% mm_diff = HVA_cs_med-VMA_cs_med;
% [mm_diff_ave, ~, ~, mm_diff_SEM] = getCI(mm_diff, 2, 1);
% errorbar(1.5, diff_y, mm_diff_SEM, 'k.', 'CapSize', 0, 'linewidth', wd_barEdge)
% plot([1,2], [diff_y, diff_y], 'k-', 'linewidth', wd_barEdge)
% 
% saveas(gcf, 'XueCarrasco/fig/RC/behav/HVA_vs_VMA.jpg')
% 
% %% plot CS
% figure('Position', [ 0 0 500 500]), hold on
% for isubj = 1:nsubj
%     plot(HVA_cs_med(isubj), VMA_cs_med(isubj), [markers_allSubj{isubj}, 'k'], 'MarkerSize', sz_marker)
% end
% % linear regression
% x = HVA_cs_med;
% lm = polyfit(x, VMA_cs_med, 1);
% x_lm2 = linspace(min(x) - std(x), max(x) + std(x), 2);
% yfit = polyval(lm, x_lm2);
% plot(x_lm2, yfit, 'k', 'linewidth', 2, 'handlevisibility', 'off');
% 
% % corr
% text_title=[];
% for ieff = 1:3
%     [coeff, p] = corr(HVA_cs_med, VMA_cs_med,'Type', corr_type_all{ieff}, 'Tail',tail_type);
%     text_title = [text_title, sprintf('%s=%.2f, p=%.3f\n', namesCoeff_all{ieff}, coeff, p)];
% end
% if p<.05, text_sig = '_sig'; else, text_sig=''; end
% %
% xline(0, 'color', ones(1,3)/2);
% yline(0, 'color', ones(1,3)/2);
% axis square
% title(sprintf('CS\n%s', text_title))
% xlabel('HVA')
% ylabel('VMA')
% set(findall(gcf, '-property', 'fontsize'), 'fontsize',20)
% set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 2)
% saveas(gcf, sprintf('%sn%d_CS%s.jpg', folderName, nsubj, text_sig))
% 
% %% plot tunC
% for ifeature=1:2
%     switch ifeature
%         case 1, nTunC = nTunC_ORI; HVA_med = HVA_ORI_med; VMA_med = VMA_ORI_med;
%         case 2, nTunC = nTunC_SF; HVA_med = HVA_SF_med; VMA_med = VMA_SF_med;
%     end
%     
%     for iT = 1:nTunC
%         figure('Position', [ 0 0 500 500]), hold on
%         for isubj = 1:nsubj
%             plot(HVA_med(isubj, iT), VMA_med(isubj, iT), [markers_allSubj{isubj}, 'k'], 'MarkerSize', sz_marker)
%             
%         end
%         % linear regression
%         x = HVA_med(:, iT);
%         lm = polyfit(x, VMA_med(:, iT), 1);
%         x_lm2 = linspace(min(x) - std(x), max(x) + std(x), 2);
%         yfit = polyval(lm, x_lm2);
%         plot(x_lm2, yfit, 'k', 'linewidth', 2, 'handlevisibility', 'off');
%         
%         % corr
%         text_title=[];
%         for ieff = 1:3
%             [coeff, p] = corr(HVA_med(:, iT), VMA_med(:, iT),'Type', corr_type_all{ieff}, 'Tail', tail_type);
%             text_title = [text_title, sprintf('%s=%.2f, p=%.3f\n', namesCoeff_all{ieff}, coeff, p)];
%         end
%         if p<.05, text_sig = '_sig'; else, text_sig=''; end
%         %
%         xline(0, 'color', ones(1,3)/2);
%         yline(0, 'color', ones(1,3)/2);
%         axis square
%         title(sprintf('%s %s\n%s', namesFeature{ifeature}, namesTunC_unit_perF{ifamily_perF(ifeature), 2}{iT}, text_title))
%         xlabel('HVA')
%         ylabel('VMA')
%         set(findall(gcf, '-property', 'fontsize'), 'fontsize',20)
%         set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 2)
%         saveas(gcf, sprintf('%sn%d_%s%d%s.jpg', folderName, nsubj, namesFeature{ifeature}, iT, text_sig))
%     end % iT
%     
%     
% end % ifeature
% 
% %% plot NOM params
% % for ip = 1:2
% %     figure('Position', [ 0 0 500 500]), hold on
% %     for isubj = 1:nsubj
% %         plot(HVA_NOM(isubj, ip), VMA_NOM(isubj, ip), [markers_allSubj{isubj}, 'k'], 'MarkerSize', sz_marker)
% %
% %     end
% %     % linear regression
% %     x = HVA_NOM(:, ip);
% %     lm = polyfit(x, VMA_NOM(:, ip), 1);
% %     x_lm2 = linspace(min(x) - std(x), max(x) + std(x), 2);
% %     yfit = polyval(lm, x_lm2);
% %     plot(x_lm2, yfit, 'k', 'linewidth', 2, 'handlevisibility', 'off');
% %
% %     % corr
% %     [coeff, p] = corr(HVA_NOM(:, ip), VMA_NOM(:, ip));
% %     if p<.05, text_sig = '_sig'; else, text_sig=''; end
% %     %
% %     xline(0, 'color', ones(1,3)/2);
% %     yline(0, 'color', ones(1,3)/2);
% %     axis square
% %     title(sprintf('%s \nr=%.2f, p=%.3f', namesParamsModel_all{1}{ip}, coeff, p))
% %     xlabel('HVA')
% %     ylabel('VMA')
% %     set(findall(gcf, '-property', 'fontsize'), 'fontsize',30)
% %     set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 2)
% %     saveas(gcf, sprintf('%sn%d_NOM%d%s.jpg', folderName, nsubj, ip, text_sig))
% % end % ip
