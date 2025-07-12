function PRplot_params_est_v3(IV_sim_allB, params_true, truth, params_est_allB, nLL_allB, nLL_true, titles)

nplots = truth.nparams_model+2;

%% plot data 
subplot(1,nplots,1), hold on
histogram(squeeze(IV_sim_allB(1,1,:)), 'DisplayStyle', 'stair', 'edgecolor', 'r', 'Normalization', 'Probability')
histogram(squeeze(IV_sim_allB(1,2,:)), 'DisplayStyle', 'stair', 'edgecolor', 'b', 'Normalization', 'Probability')
xline(params_true(2), 'linewidth', 2);
ylabel('Probability')
xlabel('Internal variable')
legend({'IV [PRS]','IV [ABS]','True thresh'},'Location','south')

%% plot estimated params %

% if truth.nparams_model == 3, np= 2; else, np = 1;end
for ip = 1:truth.nparams_model
    subplot(1, nplots,  ip+1), hold on
    
     % estimation
    histogram(params_est_allB(:, ip), 20);
   
    % true value
    xline(params_true(ip), 'r-', 'linewidth', 2);
    
    title(sprintf('Estimated %s', titles{ip}))
end

%% plot nLL
subplot(1, nplots, nplots), hold on
histogram(nLL_allB, 20)
xline(nLL_true, 'r', 'linewidth', 2);
title('nLL')

