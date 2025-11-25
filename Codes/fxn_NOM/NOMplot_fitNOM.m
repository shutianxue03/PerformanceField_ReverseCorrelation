
% This script generates plots for the estimated parameters and prediction metrics of the Noisy Observer Model (NOM).

%% Compile data and pred for all iterations and bins
pYES_data_allBins = nan(nBoot, nBins);
pYES_pred_allBins = pYES_data_allBins;
pC_data_allBins = pYES_data_allBins;
pC_pred_allBins = pYES_data_allBins;
pA_data_allBins = pYES_data_allBins;
pA_pred_allBins = pYES_data_allBins;
nData_allBins = pYES_data_allBins;
IV_allBins_all = pYES_data_allBins;

for iBoot = 1:nBoot
    % IVbin_allB(ii, :) = data_allB{ii}.IV_allBins;
    IV_allBins_all(iBoot, :) = pred.metrics.IV_allBins;

    % dprime_data_allB(ii, :) = pred_metrics_allBoot{ii}.metrics.dprime_data_allBins;
    % dprime_pred_allB(ii, :) = pred_metrics_allBoot{ii}.metrics.dprime_pred_allBins;

    pYES_data_allBins(iBoot, :) = pred_metrics_allBoot{iBoot}.metrics.pYES_data_allBins;
    pYES_pred_allBins(iBoot, :) = pred_metrics_allBoot{iBoot}.metrics.pYES_pred_allBins;

    pC_data_allBins(iBoot, :) = pred_metrics_allBoot{iBoot}.metrics.pC_data_allBins;
    pC_pred_allBins(iBoot, :) = pred_metrics_allBoot{iBoot}.metrics.pC_pred_allBins;

    pA_data_allBins(iBoot, :) = pred_metrics_allBoot{iBoot}.metrics.pA_data_allBins;
    pA_pred_allBins(iBoot, :) = pred_metrics_allBoot{iBoot}.metrics.pA_pred_allBins;

    nData_allBins(iBoot, :) = data_allBoot{iBoot}.nTrials_allBins;
end

%% pYES/pC/pA as a fxn of binned IV
[IV_allBins_med] = getCI(IV_allBins_all, 1, 1);
[nData_allB_med] = getCI(nData_allBins, 1, 1);
% Obtain median and CI of data across iteractions (dot+errorbars)
[pYES_data_med, ~, ~, pYES_data_SEM_neg, pYES_data_SEM_pos] = getCI(pYES_data_allBins, 1, 1);
[pC_data_med, ~, ~, pC_data_SEM_neg, pC_data_SEM_pos] = getCI(pC_data_allBins, 1, 1);
[pA_data_med, ~, ~, pA_data_SEM_neg, pA_data_SEM_pos] = getCI(pA_data_allBins, 1, 1);
% Obtain median and CI of pred across iteractions (line+shaded errorbars)
[pYES_pred_med, pYES_pred_lb, pYES_pred_ub] = getCI(pYES_pred_allBins, 1, 1);
[pC_pred_med, pC_pred_lb, pC_pred_ub] = getCI(pC_pred_allBins, 1, 1);
[pA_pred_med, pA_pred_lb, pA_pred_ub] = getCI(pA_pred_allBins, 1, 1);

sz_scale = 80;

figure('Position', [0 200 300 800])
subplot(3,1,1), hold on
% plot(IV_allBins, pYES_pred_allBins, 'k-')
for iBin=1:nBins
    % Data (dot+errorbars)
    plot(IV_allBins_med(iBin), pYES_data_med(iBin), 'ko', 'MarkerSize', nData_allB_med(iBin)/sz_scale+5)
    errorbar(IV_allBins_med(iBin), pYES_data_med(iBin), pYES_data_SEM_neg(iBin), pYES_data_SEM_pos(iBin), 'k', 'CapSize', 0)
end
% Pred (shaded errorbars)
plot(IV_allBins_med, pYES_pred_med, 'k-')
patch([IV_allBins_med, fliplr(IV_allBins_med)], [pYES_pred_lb, fliplr(pYES_pred_ub)], 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');
% xline(pred.criterion_IV, 'r-');
yline(.5, 'k--');
ylim([0,1])
xlabel('Binned IV')
ylabel('pYES')
% add r and R-squared
r_pYES = corr(IV_allBins_med', pYES_data_med');
R2_pYES = 1 - sum((pYES_data_med - pYES_pred_med).^2) / sum((pYES_data_med - mean(pYES_data_med)).^2);
title(sprintf('pYES (r=%.2f, R²=%.2f)', r_pYES, R2_pYES));

subplot(3,1,2), hold on
for iBin = 1:nBins
    % Data (dot+errorbars)
    plot(IV_allBins_med(iBin), pC_data_med(iBin), 'ko', 'MarkerSize', nData_allB_med(iBin)/sz_scale+5)
    errorbar(IV_allBins_med(iBin), pC_data_med(iBin), pC_data_SEM_neg(iBin), pC_data_SEM_pos(iBin), 'k', 'CapSize', 0)
end
% Pred (shaded errorbars)
plot(IV_allBins_med, pC_pred_med, 'k-')
patch([IV_allBins_med, fliplr(IV_allBins_med)], [pC_pred_lb, fliplr(pC_pred_ub)], 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');
% xline(pred.criterion_IV, 'r-');
yline(.5, 'k--');
yline(mean(pC_data_med), 'k-');
ylim([0,1])
xlabel('Binned IV')
ylabel('pC')
% add r and R-squared
r_pC = corr(pC_pred_med', pC_data_med');
R2_pC = 1 - sum((pC_data_med - pC_pred_med).^2) / sum((pC_data_med - mean(pC_data_med)).^2);
title(sprintf('pC (r=%.2f, R²=%.2f)', r_pC, R2_pC));

subplot(3,1,3), hold on
for iBin = 1:nBins
    % Data (dot+errorbars)
    plot(IV_allBins_med(iBin), pA_data_med(iBin), 'ko', 'MarkerSize', nData_allB_med(iBin)/sz_scale+5)
    errorbar(IV_allBins_med(iBin), pA_data_med(iBin), pA_data_SEM_neg(iBin), pA_data_SEM_pos(iBin), 'k', 'CapSize', 0)
    % Plot predictied pA from measured pYES
    % plot(IV_allBins_med(iBin), pYES_data_med(iBin)^2+(1-pYES_data_med(iBin))^2, 'co', 'MarkerSize', nData_allB_med(iBin)/sz_scale+5)
end
% Pred (shaded errorbars)
plot(IV_allBins_med, pA_pred_med, 'k-');
patch([IV_allBins_med, fliplr(IV_allBins_med)], [pA_pred_lb, fliplr(pA_pred_ub)], 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');
% xline(pred.criterion_IV, 'r-');
ylim([.5,1])
xlabel('Binned IV')
ylabel('pA')
% add r and R-squared
r_pA = corr(pA_pred_med', pA_data_med');
R2_pA = 1 - sum((pA_data_med - pA_pred_med).^2) / sum((pA_data_med - mean(pA_data_med)).^2);
title(sprintf('pA (r=%.2f, R²=%.2f)', r_pA, R2_pA));

sgtitle(sprintf('Figure 1. Metrics as a fxn of binned IV\n%s (ModelA%dB%d %s, niterations=%d)', ...
    subjName, iModelA, iModelB, namesModelB{iModelB}, nBoot))

saveas(gcf, sprintf('%s/3Metrics_%s_L%d_A%dB%d.jpg', nameFolder_Figures_NOM, subjName, iLocComb, iModelA, iModelB))

%% Plot estimated parameters across iterations

figure('Position', [100,100,nParams*400,400]);
for iParam = 1:nParams

    subplot(1, nParams, iParam);
    plot(params_est_allBoot(:, iParam), '-o');
    yline(mean(params_est_allBoot(:, iParam)), 'k-')

    % Load and plot the true param
    % load(sprintf('%s/truth.mat', nameFolder_OOD_load), 'noiseP_true')
    % yline(noiseP_true, 'r-'); % should match noiseP defined in OOD_sim

    ylim([params_lb(iParam), params_ub(iParam)])
    xlabel('Bootstrap');
    ylabel(namesModelBparams{iModelB}{iParam});
    title(['Parameter: ', namesModelBparams{iModelB}{iParam}]);

end

sgtitle(sprintf('Figure 2. Estimated Parameters Across Bootstraps\n%s (ModelA%dB%d %s, niterations=%d)', ...
    subjName, iModelA, iModelB, namesModelB{iModelB}, nBoot))

saveas(gcf, sprintf('%s/4Params_%s_L%d_A%dB%d.jpg', nameFolder_Figures_NOM, subjName, iLocComb, iModelA, iModelB))
