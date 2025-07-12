% collect the predictions from OOD

clc
close all
alpha_all = [.1, .2, .5, .7, .9, 1, 1.2, 1.5]; % alpha_PRS; alpha_ABS is always prop. to alpha_PRS
% alpha_all = [1, 1.2, 1.5]; % alpha_PRS; alpha_ABS is always prop. to alpha_PRS
nalpha = length(alpha_all);

data_allALPHA = cell(1, nalpha);
pred_allALPHA = data_allALPHA;
params_est_allALPHA = data_allALPHA;

%% load data
for ialpha = 1:nalpha
    % load data
    load(sprintf('Data_OOD/ModelRecov_alpha%d', round(alpha_all(ialpha)*10)))
    
    % save for each alpha
    data_allALPHA{ialpha} = metric_sim_allB;
    pred_allALPHA{ialpha} = metric_pred_allB;
    params_est_allALPHA{ialpha}  = params_est_allB;
end

%% plot predicted metrics
figure('Position', [2000 0 1000 1000])
plot_allALPHA(1, 1, 'Sensitivity (d'')', data_allALPHA, pred_allALPHA, alpha_all, [0, 2])
plot_allALPHA(2, 2, 'Criterion', data_allALPHA, pred_allALPHA, alpha_all, [-1, 1])
plot_allALPHA(4, 3, 'Proportion correct', data_allALPHA, pred_allALPHA, alpha_all, [0, 1])
plot_allALPHA(5, 4, 'Hit rate', data_allALPHA, pred_allALPHA, alpha_all, [0, 1])
plot_allALPHA(6, 5, 'False alarm rate', data_allALPHA, pred_allALPHA, alpha_all, [0, 1])
plot_allALPHA(7, 6, 'Prop. agreement', data_allALPHA, pred_allALPHA, alpha_all, [.5, 1])
plot_allALPHA(8, 7, 'Prop. agreement [Gabor-present trials]', data_allALPHA, pred_allALPHA, alpha_all, [.5, 1])
plot_allALPHA(9, 8, 'Prop. agreement [Gabor-absent trials]', data_allALPHA, pred_allALPHA, alpha_all, [.5, 1])

%% plot estimated params
titles_paramsEST = {'alpha [PRS]', 'alpha [ABS]', 'thresh'};
% if truth.nparams_model == 3, np= 2; else, np = 1;end
for ip = 1:truth.nparams_model
    subplot(4, 3, 9 + ip), hold on
    for ialpha = 1:nalpha
        
        % true value
        switch ip
            case 1, bar(ialpha, alpha_all(ialpha));
            case 2, bar(ialpha, alpha_all(ialpha) * truth.ratio_PRSABS); 
            case 3, yline(truth.thresh, 'k--');
        end
        
        % prediction
        params_est = params_est_allALPHA{ialpha}(:, ip);
        params_est_med = median(params_est);
        params_est_CI68 = quantile(params_est, [.16, .84]) ;
        errorbar(ialpha, params_est_med, params_est_med-params_est_CI68(1), params_est_CI68(2)-params_est_med, 'ok-')
        
    end % end of ialpha
    xticks(1:nalpha)
    xticklabels(alpha_all)
    xlim([0, nalpha+1])
    
    if ip ~= 3, ylim([0, 3]), else, ylim([truth.thresh-.1, truth.thresh+.1]), end
    if ip == 3, legend({'True threshold', 'Estimated threshold'}), end
    
    title(sprintf('Estimated %s', titles_paramsEST{ip}))
end
