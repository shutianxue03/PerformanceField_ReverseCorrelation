% 
% %% plot settings
% sz_ticks = 15;
% sz_title = 20;
% 
% %%
% 
% 
% 
% 
% subplot(nparams, ncomp, iplotInd(iparam, icomp))
% hold on, box on, grid on
% 
% 
% 
% % SF peak, need to be converted to log scale
% if (ifeature == 2) && (iparam == 1), params_all(:, iparam) = log2(params_all(:, iparam)); end
% 
% % plot params
% bar(1, params_all(1, iparam), 'FaceColor', color1, 'edgecolor', color1, 'faceAlpha', .5, 'barwidth', .3)
% bar(2, params_all(2, iparam), 'FaceColor', color2, 'edgecolor', color2, 'faceAlpha', .5, 'barwidth', .3)
% 
% xticks([])
% ax = gca; ax.FontSize = sz_ticks;
% xlim([.5, 2.5])
% 
% if (ifeature == 2) && (iparam == 1)
%     yticks(log2(yticks_p{ifeature}(iparam, :)))
%     ylim(log2(yticks_p{ifeature}(iparam, [1,end])))
%     yticklabels(round(yticks_p{ifeature}(iparam, :), 1))
%     yline(1);
% else
%     yticks(yticks_p{ifeature}(iparam, :))
%     ylim(yticks_p{ifeature}(iparam, [1,end]))
%     yticklabels(round(yticks_p{ifeature}(iparam, :), 2))
% end
% 
% set(gca, 'YGrid', 'on', 'XGrid', 'off')
% 
% title(paramsNames{iparam}, 'fontsize', sz_title)
% % end
% 
% set(findall(gcf, '-property', 'linewidth'), 'linewidth',2)
% saveas(gcf, sprintf('publishedPDFs/fig/SP_params_%s_comb%d.jpg', kernelNames{ifeature}, icomp))
% 
