
% bm stands for behavioral measurements

close all

sz_title = 55;
sz_label = 40;
sz_ticks = 50;
wd_border = 4;
sz_marker = 800;

% switch im
%     case 1
%         title_full = 'd''';
%         title_ = 'd''';
%         yline_ = dprime_theo;
%         yticks_ = [1, 1.2, 1.4];
%         ydiffErr = 1.3;
%         ylim_ = [.96, 1.44];
%         ylim_all = [.72, 1.68]; % axis limit in the scatter plot
%         yticks_all = .8:.4:1.6; % axis ticks in the scatter plot
%         yticklabels_ = yticks_all;
%     case 2
%         title_full = 'Criterion';
%         title_ = 'Criterion';
% %         data_ = criterion_allT_allSubj;
%         yline_ = 0;
%         ylim_ = [-.6, .6];
%         yticks_ = [-.5, 0, .5];
%         ydiffErr = .2;
%         ylim_all = [-.6, .6]; % axis limit in the scatter plot
%         yticks_all = -.5:.5:.5; % axis ticks in the scatter plot
%         yticklabels_ = yticks_all;
%     case 3
%         title_full = 'RT (ms)';
%         title_ = 'RT';
% %         data_ = squeeze(100*exp(data_aveB_allSubj(:, 5, :))+500);
%         ylim_ = [506.8, 521.2];
%         yticks_ = [508, 514, 520];
%         ydiffErr = 515;
%         ylim_all = [498, 522]; % axis limit in the scatter plot
%         yticks_all = 500:10:520; % axis ticks in the scatter plot
%         yticklabels_ = yticks_all;
%     case 4
%         title_full = 'Titrated gabor contrast (%)';
%         title_ = 'cst';
% %         data_ = cst_ave_allSubj;
%         ylim_ = log10([25, 55]/100);
%         yticks_ = log10([25, 40, 55]/100);
%         yticklabels_ = [25, 40, 55];
%         ydiffErr = -.35;
%         ylim_all = [log10(.23), log10(.57)]; % axis limit in the scatter plot
%         yticks_all = log10(.25:.15:.55); % axis ticks in the scatter plot
%     case 5
%         title_full = 'Resp consistency (%)';
%         title_ = 'pA';
% %         data_ = squeeze(pA_ave_allSubj(:, 3, :))*100;
%         ylim_ = [60, 90];
%         yticks_ = 60:15:90;
%         ydiffErr = 80;
%         ylim_all = [58, 82]; % axis limit in the scatter plot
%         yticks_all = 60:10:80; % axis ticks in the scatter plot
%         yticklabels_ = yticks_all;
%
% end

bm_ave_perS = squeeze(data_aveB_allSubj(:, im, :));
bm_sem_perS= squeeze(data_semB_allSubj(:, im, :));

if nsubj==1
    bm_ave_perS = bm_ave_perS';
    bm_sem_perS = bm_sem_perS';
end
bm_aveS = mean(bm_ave_perS, 1);
bm_semS = std(bm_ave_perS, [], 1)/sqrt(nsubj); % not SEM

% ydiff_buffer = (yticks_(end) - yticks_(1))/10;

for icomb = 1:ncomb4
    subplot(1,ncomb4, icomb), hold on
    % ind of locations
    ind1 = combInd(icomb, 1);
    ind2 = combInd(icomb, 2);
    
    % get the difference between two loc
    %     diff_sem = std(bm_ave_perS(:,ind1) - bm_ave_perS(:, ind2))/sqrt(nsubj);
    %
    %     % paired t-test
    %     if nsubj>1
    %         if imeasure == 3
    %             [~, p,~, ~] = ttest(exp(bm_ave_perS(:,ind1)), exp(bm_ave_perS(:,ind2)));
    %         else
    %             [~, p,~, ~] = ttest(bm_ave_perS(:,ind1), bm_ave_perS(:,ind2));
    %         end
    %     end
    
    %% scatter plot
    subplot(1,ncomb, icomb)
    hold on, box on
    ymin_all = min(bm_ave_perS(:));
    ymax_all = max(bm_ave_perS(:));
    
    % IDVD data (black dot and errorbar for block sem)
    errorbar(bm_ave_perS(:, ind1), bm_ave_perS(:, ind2), bm_sem_perS(:, ind1), 'horizontal', '.k', 'CapSize', 0)
    errorbar(bm_ave_perS(:, ind1), bm_ave_perS(:, ind2), bm_sem_perS(:, ind2), 'vertical', '.k', 'CapSize', 0)
    scatter(bm_ave_perS(:, ind1), bm_ave_perS(:, ind2), sz_marker, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k')
    
    % group average (purple cross
    if nsubj>1
        errorbar(bm_aveS(:,ind1), bm_aveS(:,ind2), bm_aveS(:,ind1), 'horizontal', '.','color', [.75,0,.75], 'CapSize', 0, 'linewidth', 8)
        errorbar(bm_aveS(:,ind1), bm_aveS(:,ind2), bm_aveS(:,ind2), 'vertical', '.','color', [.75,0,.75], 'CapSize', 0, 'linewidth', 8)
    end
    
    % diagonal line
    %     plot(ylim_all, ylim_all, 'k-', 'linewidth', wd_border)
    
    % extra line
    %     if imeasure == 2 % criterion
    %         yline(0, '-', 'color', ones(1,3)*.5, 'linewidth', wd_border/2);
    %         xline(0, '-', 'color', ones(1,3)*.5, 'linewidth', wd_border/2);
    %     end
    
    % limit
    %     xlim(ylim_all)
    %     ylim(ylim_all)
    %
    %     % ticks and ticklabels
    %     xticks(yticks_all)
    %     yticks(yticks_all)
    %     switch imeasure
    %         case 4 % CST
    %             xticklabels(round(10.^(yticks_all)*100))
    %             yticklabels(round(10.^(yticks_all)*100))
    %     end
    %
    %     % labels
    %     if sum(im == [1,4])
    %         xlabel(names2Loc{1}, 'FontSize', sz_ticks)
    %         ylabel(names2Loc{2}, 'FontSize', sz_ticks)
    %     end
    %     axis square
    %     ax = gca;
    %     ax.XAxis.FontSize = sz_ticks;
    %     ax.YAxis.FontSize = sz_ticks;
    %     ax.LineWidth = wd_border;
    %
    %     title(title_full, 'FontSize', sz_title)
    
end % end of icomb


