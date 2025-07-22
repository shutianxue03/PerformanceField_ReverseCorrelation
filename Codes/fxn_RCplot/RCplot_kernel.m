

figure('Position', [0 800 nLoc*300 300])
ylimit = [-.1, .4];

itype = 2; % target-absent

for iLoc = 1:nLoc
    if isubjPlot, x = squeeze(kernel_allLoc(itype, iLoc,:)); x_pred = squeeze(kernel_pred_allLoc(itype, iLoc,:));
    else, x = squeeze(kernel_allSubj(:, itype, iLoc,:)); x_pred = squeeze(kernel_pred_allSubj(:,itype, iLoc,:));
    end
    
    subplot(1, nLoc, iLoc), hold on
    
    if isubjPlot % each subj
        plot(filterSF_all_log, x, 'ko', 'handlevisibility', 'off')
        plot(filterSF_all_log, x_pred, 'k-', 'handlevisibility', 'off')
    else % all subjs
        weights = ntrials_allSubj/sum(ntrials_allSubj);
        % plot individual data
        for isubj_p = 1:nsubj
            mkSize = weights(isubj_p) * 15;
            plot(filterSF_all_log, x(isubj_p, :), 'o', 'MarkerEdgeColor', colors(:, isubj_p), 'MarkerSize', mkSize, 'handlevisibility', 'off')
%             plot(filterSF_all_log, x_pred(isubj_p, :), '-', 'color', colors(:, isubj_p), 'handlevisibility', 'off')
        end
        
        % plot ave data
        if sum(find(iLoc == [1,3,5])), kernel_ave = weights' * x; %nanmean(x, 1);
        else, kernel_ave = nanmean(x, 1);
        end
        kernel_sem = nanstd(x, [], 1)/sqrt(nsubj);
        [kernel_ave_pred, params_ave_est, nLL_ave] = SX_sim08_fit(filterSF_all, kernel_ave, model, fitMode);
%         errorbar(filterSF_all_log, kernel_ave, kernel_sem, 'ok', 'MarkerFaceColor', 'k', 'handlevisibility', 'off')
%         plot(filterSF_all_log, kernel_ave, 'ok', 'MarkerFaceColor', 'k', 'handlevisibility', 'off')
        patch([filterSF_all_log, flip(filterSF_all_log)], [kernel_ave - kernel_sem, flip(kernel_ave + kernel_sem)], [.5,.5,.5], 'EdgeColor','none','FaceAlpha',.5,  'handlevisibility', 'off')
        plot(filterSF_all_log, kernel_ave_pred, '-k', 'MarkerFaceColor', 'k','linewidth', 2, 'handlevisibility', 'off')
    end
    
    if isubjPlot, plot(log2([params_est_allLoc(itype, iLoc, 1), params_est_allLoc(itype, iLoc, 1)]), ylimit, 'k-') % estimated peak
    else, plot(log2([params_ave_est(1), params_ave_est(1)]), ylimit, 'k-'), end % estimated peak
    plot(log2([2,2]), ylimit, 'k', 'color',[.5,.5,.5])
    plot(log2([.8,4.8]), [0 0], 'k', 'color',[.5,.5,.5], 'handlevisibility', 'off')
    
    if iLoc== 1
        xlabel('SF channel (cpd)'), ylabel('kernel (a.u.)')
%         legend('est peak SF', 'target SF (2 cpd)', 'location', 'best')
    end
    
    xticks(log2([1,2,4])), xticklabels([1,2,4])
    xlim(log2([.8,4.8]))
    ylim(ylimit)
    title(locNames{iLoc})
    
    %     text(log(1.5), ylimit(1), sprintf('center=%.2f\nsigma=%.2f\nalpha=%.2f\nb=%.2f\n', params_ave_est))
end

%
if isubjPlot, sgtitle(sprintf('fitted by %s\n%s %d trials per loc', modelName, subjName, nAllTrials))
else, sgtitle(sprintf('fitted by %s\nn = %d (%d \\pm %d trials per loc)', modelName, nsubj, round(mean(nAllTrials)), round(std(nAllTrials))))
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
