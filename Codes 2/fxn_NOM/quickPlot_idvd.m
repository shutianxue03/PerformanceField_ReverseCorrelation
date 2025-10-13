
metrics_sim = nan(ni, nmetrics);
for ii = 1:ni
    metrics_sim(ii, :) = data_allB{ii}.metrics_sim;
end

%%
% dprime, criterion, pC, pHit, pFA, pA, pA_PRS, pA_ABS
ticksLim = {[1,2], [-1, 1], [.5, .9], [.7, 1], [.1, .6], [.5, .9], [.5, .9], [.5, .9]};

%% plot predictions
figure('Position', [3e3 0 1500 200])
for im = 1:nmetrics
    pred_im = squeeze(pred_allB(:, im));
    subplot(2,nmetrics, im), hold on
    histogram(pred_im)
%     xline(data.metrics_sim(im), 'k-', 'linewidth', 2);
    title(namesMetrics{im})
    
    subplot(2,nmetrics, im+nmetrics), hold on
    scatter(metrics_sim(:, im), pred_im)
    xlim(ticksLim{im})
    ylim(ticksLim{im})
    plot(ticksLim{im}, ticksLim{im}, 'k-')
    axis square 
    xlabel('Data (test set)')
    ylabel('prediction')
%     bar(median(metrics_sim(:, im)))
%     errorbar(median(pred_im), std(pred_im), 'o')
%     if find(im == [3,4,6:8]), ylim([.5, 1])
%     elseif im == 1,  ylim([0, 2])
%     elseif im == 2, ylim([-1, 1])
%     elseif im == 5,  ylim([0, 1])
%     end
end
sgtitle(sprintf('%s %s',namesLocComb{iLocComb}, subjName))

%% plot estimated params
figure('Position', [3e3 0 600 200])
for ip = 1:nparams_model
    p = squeeze(params_est_allB(:, ip));
    subplot(2,nparams_model, ip), hold on
    histogram(p)
    
    subplot(2,nparams_model, ip+nparams_model), hold on
    errorbar(median(p), std(p), 'o')
end
sgtitle(sprintf('%s %s',namesLocComb{iLocComb}, subjName))

%% estimated params vs. pA
figure('Position', [3e3 0 600 200])
for ip = 1:nparams_model
    subplot(1,nparams_model, ip), hold on
    plot(metrics_sim(:, 6), params_est_allB(:, ip), 'o')
    [r,p] = corr(metrics_sim(:, 6), params_est_allB(:, ip));
    title(sprintf('%s\nr = %.2f, p = %.2f', namesParamsModel{ip}, r, p))
end


