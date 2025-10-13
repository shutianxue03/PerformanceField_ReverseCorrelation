% ModelPlot_local
% plot model predictions for one particular variation for all/one subjects

%% get median and CI
[data_med, ~, ~, data_neg, data_pos] = getCI(data_metrics_allB);
[pred_med, ~, ~, pred_neg, pred_pos] = getCI(pred_metrics_allB);
[params_est_med, ~, ~, params_est_neg, params_est_pos] = getCI(params_est_allB);
[IC_med, ~, ~, IC_neg, IC_pos] = getCI(IC_allB);

%% Fig 1: metrics data vs. pred
xlim_allM = [.7, 2; -1, 1; .5, 1; .5, 1; .1, .6; .5, .9; .5, .9;.5, .9];
figure('Position', [0 0 1000 600])

for im = 1:nmetrics
    subplot(2,4,im), hold on
    plot(data_metrics_allB(:, im), pred_metrics_allB(:, im), '.', 'color', colors_comb(iLocComb, :))
    errorbar(data_med(:, im), pred_med(:, im), ...
        pred_neg(:, im), pred_pos(:, im), ...
        data_neg(:, im), data_pos(:, im), 'k.', 'CapSize', 0)
    plot([xlim_allM(im, 1), xlim_allM(im, 2)], [xlim_allM(im, 1), xlim_allM(im, 2)], 'k-')
    xlim(xlim_allM(im, :))
    ylim(xlim_allM(im, :))
    axis square
    title(namesMetrics{im})
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)

%% Fig 2: pA vs. params
nparams_est_model = size(params_est_allB, 2);
figure('Position', [0 0 1000 500])
for ip = 1:nparams_est_model
    subplot(1,nparams_est_model,ip), hold on
    plot(data_metrics_allB(:, 6), params_est_allB(:, ip), '.', 'color', colors_comb(iLocComb, :))
    errorbar(data_med(:, 6), params_est_med(:, ip), ...
        params_est_neg(:, ip), params_est_pos(:, ip), data_neg(:, 6), data_pos(:, 6), 'k.', 'CapSize', 0)
    
    [r, p] = corr(data_metrics_allB(:, 6), params_est_allB(:, ip));
    axis square
    xlabel('pA')
    ylabel(namesParams_NOM{ip})
    xlim([.55, .85])
    title(sprintf('r=%.3f, p=%.3f', r, p))
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)

%% Fig 3. IC
figure, hold on
errorbar(1:nIC, IC_med, IC_neg, IC_pos, 'o')
xticks(1:nIC)
xticklabels(namesIC)
xlim([0, nIC+1])
[p,tab] = anova1(IC_allB, [], 'off');
title(sprintf('Main effect of type of IC: p=%.3f', p))

% %% Figure 4: template
% figure
% subplot(1,2,1), RCplot_2Dkernel(squeeze(mean(kernel2D_allB, 1))'), title('Mean of all templates')
% subplot(1,2,2), RCplot_2Dkernel(squeeze(std(kernel2D_allB, [], 1))'), title('SD of all templates')
