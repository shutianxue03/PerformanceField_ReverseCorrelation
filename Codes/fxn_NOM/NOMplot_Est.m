
% This script generates plots for the estimated parameters and prediction metrics of the Noisy Observer Model (NOM).


pA_data_allB = nan(ni, nBins);
pA_pred_med = pA_data_allB;
IV_allB = pA_data_allB;

figure, hold on
for ii=1:ni
    IV_allB(ii, :) = pred_metrics_allB{ii}.metrics.IV_allBins;
    pA_data_allB(ii, :) = pred_metrics_allB{ii}.metrics.pA_data_allBins;
    pA_pred_allB(ii, :) = pred_metrics_allB{ii}.metrics.pA_pred_allBins;
    pC_data_allB(ii, :) = pred_metrics_allB{ii}.metrics.pC_data_allBins;
    pC_pred_allB(ii, :) = pred_metrics_allB{ii}.metrics.pC_pred_allBins;

    plot(IV_allB(ii, :), pC_data_allB(ii, :), 'o')
    plot(IV_allB(ii, :), pC_pred_allB(ii, :), '-')

    plot(ii, pC_data_allB(ii, :), 'o')
    plot(ii, pC_pred_allB(ii, :), '-')
    ylim([.5, 1])
    xlabel('IV')
    xlabel('IV')
    % pause
end

%% 
figure
for ii=1:ni

    subplot(1,2,1), hold on
    plot(ii, pred_metrics_allB{ii}.metrics.pA_data_allBins, 'o')
    plot(ii, pred_metrics_allB{ii}.metrics.pA_pred_allBins, '-')
    ylim([.5, 1])
    xlabel('Iteration #')
    ylabel('pA')
title

    subplot(1,2,2), hold on
    plot(ii, pred_metrics_allB{ii}.metrics.pA_data_allBins, 'o')
    plot(ii, pred_metrics_allB{ii}.metrics.pA_pred_allBins, '-')
    ylim([.5, 1])
    xlabel('Iteration #')
    ylabel('pA')
    
    % pause
end


IV_med = getCI(IV_allB, 1, 1);
pA_data_med = getCI(pC_data_allB, 1, 1);
pA_pred_med = getCI(pC_pred_allB, 1, 1);

% Plot 
figure, hold on
plot(IV_med, pA_data_med, 'o')
plot(IV_med, pA_pred_med, '-')