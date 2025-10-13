% function PRplot_params_est(params_true, truth, params_est_allB, nLL_allB, nLL_true_allB, titles)

figure('Position', [0 0 1000 200])

nplots = truth.nparams_model+2;

%% plot estimated params %
for ip = 1:truth.nparams_model
    subplot(1, nplots, ip), hold on
    
     % estimation
    histogram(params_est_allB(:, ip), 20);
   
    % true value
    xline(params_true(ip), 'r-', 'linewidth', 2);
    
    title(sprintf('Estimated %s', namesParamsModel{ip}))
end

%% plot nLL & median of the true nLL given each data set
[nLL_true_med, nLL_true_lb, nLL_true_ub] = getCI(nLL_true_allB);
subplot(1, nplots, nplots-1), hold on
histogram(nLL_allB, 20)
xline(nLL_true_med, 'r', 'linewidth', 2);
xline(nLL_true_lb, 'r', 'linewidth', 1);
xline(nLL_true_ub, 'r', 'linewidth', 1);
title('nLL')

%% plot estimted thresh
subplot(1, nplots, nplots), hold on
histogram(thresh_est_allB, 20)
xline(params_true(end), 'r', 'linewidth', 2);
title('Estimated thresh given the measured criterion')
