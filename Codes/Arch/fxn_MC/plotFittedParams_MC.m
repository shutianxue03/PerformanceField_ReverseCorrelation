% 
% 
% % model_problematic = input('Enter the problematic model index in a vector: ');
% imodel_problematic = [];
% figure('Position', [0 0 1000 800])
% 
% for iparam = 1:nparams_full
%     subplot(nparams_full, 1,iparam), hold on
%     
%     nx = length(params_fit_all{iparam});
%     x = 1:nx;
%     
%     % bar-show all params
%     bar(params_fit_all{iparam});
%     yticks('auto')
%     ax = gca;
%     ax.YGrid = 'ON'; ax.YMinorGrid = 'ON';
%     ax.LineWidth=2;
%     ax.TickLength=[0 0];% 
%     ylabel(namesParams{iparam})
%     
%     % list the model# of ALL models
%     xticklabels_all = [];
%     for im = 1:imodel
%         xticklabels_ = ones(1,2^(paramInd_all(im, iparam))) * im;
%         xticklabels_all = [xticklabels_all, xticklabels_];
%     end
%     xticks(1:nx)
%     xticklabels(xticklabels_all)
%     
%     % highlight the bars of BAD models
%     for ip = imodel_problematic
%         istart = sum(2.^paramInd_all(1:(ip-1), iparam))+1;
%         iend = sum(2.^paramInd_all(1:ip, iparam));
%         for ibar = istart : iend
%             bar(x(ibar), params_fit_all{iparam}(ibar), 'r')
%         end
%     end
%     
% end
% set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
% 
