function plot_allALPHA(iplot, imetric, title_, data_allALPHA, pred_allALPHA, alpha_all, ylim_)

nalpha = length(alpha_all);
subplot(4,3,iplot), hold on

for ialpha = 1:nalpha
    data = data_allALPHA{ialpha}(:, imetric);
    pred = pred_allALPHA{ialpha}(:, imetric);
    
    % get ave and CI
    data_med = median(data);
    bar(ialpha, data_med)
    pred_med = median(pred);
    pred_CI68 = quantile(pred, [.16, .84]) ;
    errorbar(ialpha, pred_med, pred_med - pred_CI68(1), pred_CI68(2)-pred_med, 'ok')
    
    xticks(1:nalpha)
    xticklabels(alpha_all)
    ylim(ylim_)
end
title(title_)