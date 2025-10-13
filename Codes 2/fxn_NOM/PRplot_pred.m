% function PRplot_pred(sim, metric_sim_allB, metric_pred_allB, titles)

figure('Position',[0 300 1800 900])
% titles_long = {'Sensitivity (d'')',  'Criterion', 'Proportion correct', 'Hit rate', 'False alarm rate', 'Prop. agreement', 'Prop. agreement [PRS]', 'Prop. agreement [ABS]'};
ylims = [0, 3; -1,1; .5,1; .5,1; 0,.5; .5, 1; .5,1; .5,1];
for im = 1:8
    data = metric_sim_allB(:,im);
    pred = metric_pred_allB(:,im);
    subplot(2,4,im), hold on
    xmax = max([data; pred]) + std([data; pred])*10;
    xmin = min([data; pred]) - std([data; pred])*10;
    
    plot(data, pred, '.')
    [data_med, ~, ~, data_neg, data_pos] = getCI(data);
    [pred_med, ~, ~, pred_neg, pred_pos] = getCI(pred);
    errorbar(data_med, pred_med, pred_neg, pred_pos,data_neg, data_pos, '.r', 'CapSize', 0, 'Linewidth', 2)
%     plot([xmin, xmax], [xmin, xmax], 'k-')
    plot(ylims(im, :), ylims(im, :), 'k-')
%     xline(sim.metrics_math(im), 'r-');
%     yline(sim.metrics_math(im), 'r-');
    
    xlim(ylims(im, :))
    ylim(ylims(im, :))
    if im == 1
        xlabel('Simulated metric')
        ylabel('Predicted metric')
    end
    axis square, box on
    title(namesMetrics_short{im})
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',20)
