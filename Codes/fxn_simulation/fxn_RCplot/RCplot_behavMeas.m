
sz_marker =200;
lw = 2;

%% 1. performance change by session
% if nsubj==1
%     figure('Position', [0 500 1500 300]), hold on, grid on
%     
%     for iLoc = 1:nLoc5
%         if iLoc== 4, plot(bm_perSess(:,iLoc), '--', 'color',colors_comb(iLoc, :), 'LineWidth', .5),
%         else, plot(bm_perSess(:,iLoc), 'color',colors_comb(iLoc, :), 'LineWidth', .5), end
%     end
%     
%     ylim([min(bm_perSess(:)), max(bm_perSess(:))])
%     if im == 4 % CST
%         s_ = quantile(bm_perSess(:), linspace(0,1,7));
%         yticks(s_)
%         yticklabels(round(10.^s_*100))
%         ylim(s_([1,end]))
%     end
%     
%     xlabel('Session #')
%     ylabel(title_)
%     xtickInd = round(quantile(1:nSess, linspace(0,1,5)));
%     if nSess > 1, xlim([0, nSess + 2])
%         if nSess>2, xticks(xtickInd), xticklabels(xtickInd), else, xticks([1,2]), end
%     else, xlim([.5,2]), xticks(1)
%     end
%     
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
%     make_sgtitle(title_, subjName, nB, nsubj, nAllTrials)
% end


%% scatter plot
figure('Position', [0 0 1500 400])
for icomb = 1:ncomb4
    subplot(1,ncomb4, icomb)
    hold on, box on
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
    
    % IDVD data (black dot and errorbar for block sem)
    errorbar(bm_ave_perSubj(:, ind1), bm_ave_perSubj(:, ind2), bm_std_perSubj(:, ind1), 'horizontal', '.k', 'CapSize', 0)
    errorbar(bm_ave_perSubj(:, ind1), bm_ave_perSubj(:, ind2), bm_std_perSubj(:, ind2), 'vertical', '.k', 'CapSize', 0)
    scatter(bm_ave_perSubj(:, ind1), bm_ave_perSubj(:, ind2), sz_marker, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k')
    
    % group average (purple cross
    if nsubj>1
        errorbar(bm_aveSubj_(:,ind1), bm_aveSubj_(:,ind2), bm_semSubj_(:,ind1), 'horizontal', '.','color', [.75,0,.75], 'CapSize', 0, 'linewidth', lw)
        errorbar(bm_aveSubj_(:,ind1), bm_aveSubj_(:,ind2), bm_semSubj_(:,ind2), 'vertical', '.','color', [.75,0,.75], 'CapSize', 0, 'linewidth', lw)
    end
    
    % diagonal line
    plot(lim_scatter, lim_scatter, 'k-', 'linewidth', lw)
    
    % extra line
    if ~isnan(yline_)
    yline(yline_, '-', 'color', ones(1,3)*.5);
    xline(yline_, '-', 'color', ones(1,3)*.5);
    end
    
    % limit
    xlim(lim_scatter)
    ylim(lim_scatter)
    %
    % ticks and ticklabels
    xticks(ticks_scatter)
    yticks(ticks_scatter)
    switch im
        case 3 % RT
            xticklabels(round(exp(ticks_scatter)*100+500))
            yticklabels(round(exp(ticks_scatter)*100+500))
        case 4 % CST
            xticklabels(round(10.^(ticks_scatter)*100))
            yticklabels(round(10.^(ticks_scatter)*100))
    end
    
    xlabel(namesLocComb{ind1})
    ylabel(namesLocComb{ind2})
    
    axis square
end % end of icomb

set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
make_sgtitle(title_, subjName, nB, nsubj, nAllTrials)

