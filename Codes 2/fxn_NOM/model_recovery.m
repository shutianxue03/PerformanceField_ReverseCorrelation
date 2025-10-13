% collect the predictions from OOD

clc
close all
alpha_all = [.1, .2, .5, .7, .9, 1, 1.2, 1.5]; % alpha_PRS; alpha_ABS is always prop. tp alpha_PRS
% alpha_all = [.1]; % alpha_PRS; alpha_ABS is always prop. tp alpha_PRS
nalpha = length(alpha_all);

data_allALPHA = cell(1, nalpha);
pred_allALPHA = data_allALPHA;
params_est_allALPHA = data_allALPHA;

%% load data
for ialpha = 1:nalpha 
    % load data
    load(sprintf('Data_OOD/ModelRecov_alpha%d', round(alpha_all(ialpha)*10)))

    % save for each alpha
    data_allALPHA{ialpha} = data_all;
    pred_allALPHA{ialpha} = pred_all;
    params_est_allALPHA{ialpha}  = params_est_all;
end

%% plot predicted metrics
figure('Position', [2000 0 1000 1000])
plot_allALPHA(1, 'dprime', 'Sensitivity (d'')', data_allALPHA, pred_allALPHA, alpha_all, [0, 2])
plot_allALPHA(2, 'criterion', 'Criterion', data_allALPHA, pred_allALPHA, alpha_all, [-1, 1])
plot_allALPHA(4, 'pC', 'Proportion correct', data_allALPHA, pred_allALPHA, alpha_all, [0, 1])
plot_allALPHA(5, 'pHit', 'Hit rate', data_allALPHA, pred_allALPHA, alpha_all, [0, 1])
plot_allALPHA(6, 'pFA', 'False alarm rate', data_allALPHA, pred_allALPHA, alpha_all, [0, 1])
plot_allALPHA(7, 'pA', 'Prop. agreement', data_allALPHA, pred_allALPHA, alpha_all, [.5, 1])
plot_allALPHA(8, 'pA_PRS', 'Prop. agreement [Gabor-present trials]', data_allALPHA, pred_allALPHA, alpha_all, [.5, 1])
plot_allALPHA(9, 'pA_ABS', 'Prop. agreement [Gabor-absent trials]', data_allALPHA, pred_allALPHA, alpha_all, [.5, 1])

%% plot estimated params
titles_paramsEST = {'alpha [PRS]', 'alpha [ABS]', 'thresh'};
% if nparams_model == 3, np= 2; else, np = 1;end
for ip = 1:nparams_model
    subplot(4, 3, 9 + ip), hold on
    for ialpha = 1:nalpha
        % real value
        switch ip
            case 1, bar(ialpha, alpha_all(ialpha));
            case 2, bar(ialpha, alpha_all(ialpha) * ratio_PRSABS); % manually set
            case 3, yline(thresh, 'k--');
        end
        
        % prediction
        params_est = params_est_allALPHA{ialpha};
        pred_med = median(params_est);
        pred_CI68 = quantile(params_est, [.16, .84]) ;
        errorbar(ialpha, pred_med(:, ip), pred_med(:, ip)-pred_CI68(1, ip), pred_CI68(2, ip)-pred_med(:, ip), 'ok-')
        
        xticks(1:nalpha)
        xticklabels(alpha_all)
        xlim([0, nalpha+1])
        if ip == 3, ylim([-.5, .5]), else, ylim([0, 3]), end
        if ip == 3, legend({'True threshold', 'Estimated threshold'}), end
        
    end
    title(sprintf('Estimated %s', titles_paramsEST{ip}))
end

%% helper fxn - plot_allALPHA
function plot_allALPHA(iplot, fieldName, title_, data_allALPHA, pred_allALPHA, alpha_all, ylim_)

nalpha = length(alpha_all);
subplot(4,3,iplot), hold on

for ialpha = 1:nalpha
    data = getfield(data_allALPHA{ialpha}, fieldName);
    pred = getfield(pred_allALPHA{ialpha}, fieldName);
    
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
end