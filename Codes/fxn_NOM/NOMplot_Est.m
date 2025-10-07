
% This script generates plots for the estimated parameters and prediction metrics of the Noisy Observer Model (NOM).

%% Compile data and pred for all iteractions and bins
pYES_data_allB = nan(ni, nBins);
pYES_pred_allB = pYES_data_allB;
pC_data_allB =pYES_data_allB;
pC_pred_allB = pYES_data_allB;
pA_data_allB = pYES_data_allB;
pA_pred_allB = pYES_data_allB;
nData_allB = pYES_data_allB;
IV_allBins_all = pYES_data_allB;

% see OOD_xx_beforeEst Line179: % dprime, criterion, pC, pHit, pFA, pA,
% pA_PRS, pA_ABS, pYES
indPyes = 9;
indPc = 3;
indPa = 6;

for ii=1:ni
    % IVbin_allB(ii, :) = data_allB{ii}.IV_allBins;
    IV_allBins_all(ii, :) = pred.metrics.IV_allBins;
    
    pYES_data_allB(ii, :) = data_allB{ii}.metrics_sim(indPyes);
    pYES_pred_allB(ii, :) = pred_metrics_allB{ii}.metrics.pYES_pred_allBins;

    pC_data_allB(ii, :) = data_allB{ii}.metrics_sim(indPc);
    pC_pred_allB(ii, :) = pred_metrics_allB{ii}.metrics.pC_pred_allBins;

    pA_data_allB(ii, :) = data_allB{ii}.metrics_sim(indPa);
    pA_pred_allB(ii, :) = pred_metrics_allB{ii}.metrics.pA_pred_allBins;

    nData_allB(ii, :) = data_allB{ii}.nTrials_allBins;
end

%% pYES/pC/pA as a fxn of binned IV
[IV_allBins_med] = getCI(IV_allBins_all, 1, 1);
[nData_allB_med] = getCI(nData_allB, 1, 1);
% Obtain median and CI of data across iteractions (dot+errorbars)
[pYES_data_med, ~, ~, pYES_data_SEM_neg, pYES_data_SEM_pos] = getCI(pYES_data_allB, 1, 1);
[pC_data_med, ~, ~, pC_data_SEM_neg, pC_data_SEM_pos] = getCI(pC_data_allB, 1, 1);
[pA_data_med, ~, ~, pA_data_SEM_neg, pA_data_SEM_pos] = getCI(pA_data_allB, 1, 1);
% Obtain median and CI of pred across iteractions (line+shaded errorbars)
[pYES_pred_med, pYES_pred_lb, pYES_pred_ub] = getCI(pYES_pred_allB, 1, 1);
[pC_pred_med, pC_pred_lb, pC_pred_ub] = getCI(pC_pred_allB, 1, 1);
[pA_pred_med, pA_pred_lb, pA_pred_ub] = getCI(pA_pred_allB, 1, 1);

sz_scale = 80;

figure('Position', [0 200 800 300 ])
subplot(1,3,1), hold on
% plot(IV_allBins, pYES_pred_allBins, 'k-')
for iBin=1:nBins
    % Data (dot+errorbars)
    plot(IV_allBins_med(iBin), pYES_data_med(iBin), 'ko', 'MarkerSize', nData_allB_med(iBin)/sz_scale+5)
    errorbar(IV_allBins_med(iBin), pYES_data_med(iBin), pYES_data_SEM_neg(iBin), pYES_data_SEM_pos(iBin), 'k', 'CapSize', 0)
end
% Pred (shaded errorbars)
plot(IV_allBins_med, pYES_pred_med, 'k-')
patch([IV_allBins_med, fliplr(IV_allBins_med)], [pYES_pred_lb, fliplr(pYES_pred_ub)], 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');
xlabel('Binned IV')
ylabel('pYES')
title('pYES')

subplot(1,3,2), hold on
for iBin=1:nBins
    % Data (dot+errorbars)
    plot(IV_allBins_med(iBin), pC_data_med(iBin), 'ko', 'MarkerSize', nData_allB_med(iBin)/sz_scale+5)
    errorbar(IV_allBins_med(iBin), pC_data_med(iBin), pC_data_SEM_neg(iBin), pC_data_SEM_pos(iBin), 'k', 'CapSize', 0)
end
% Pred (shaded errorbars)
plot(IV_allBins_med, pC_pred_med, 'k-')
patch([IV_allBins_med, fliplr(IV_allBins_med)], [pC_pred_lb, fliplr(pC_pred_ub)], 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');
xlabel('Binned IV')
ylabel('pC')
title('pC')

subplot(1,3,3), hold on
for iBin=1:nBins
    % Data (dot+errorbars)
    plot(IV_allBins_med(iBin), pA_data_med(iBin), 'ko', 'MarkerSize', nData_allB_med(iBin)/sz_scale+5)
    errorbar(IV_allBins_med(iBin), pA_data_med(iBin), pA_data_SEM_neg(iBin), pA_data_SEM_pos(iBin), 'k', 'CapSize', 0)
end
% Pred (shaded errorbars)
plot(IV_allBins_med, pA_pred_med, 'k-');
patch([IV_allBins_med, fliplr(IV_allBins_med)], [pA_pred_lb, fliplr(pA_pred_ub)], 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');
xlabel('Binned IV')
ylabel('pA')
title('pA')


%% pYES: data vs. pred (the mean across bins will be 50%)
figure('Position', [100 100 2e3 1e3])
hold on
for iBin=1:nBins
    subplot(2, nBins/2, iBin), hold on
    plot(pYES_data_allB(:, iBin), pC_pred_allB(:, iBin), 'o')

    % plot ave and sd on top of the scatter plot
    errorbar(nanmean(pYES_data_allB(:, iBin)), nanmean(pYES_pred_allB(:, iBin)), nanstd(pYES_data_allB(:, iBin)), 'k', 'horizontal', 'CapSize', 0, 'LineWidth', 2)
    errorbar(nanmean(pYES_data_allB(:, iBin)), nanmean(pYES_pred_allB(:, iBin)), nanstd(pYES_pred_allB(:, iBin)), 'k', 'vertical', 'CapSize', 0, 'LineWidth', 2)

    xlabel('Emp pYES')
    ylabel('Pred pYES')
    plot([0, 1], [0, 1], 'k--')
    xline(.5, 'r--');
    yline(.5, 'r--');
    axis square
    xlim([0,1])
    ylim([0,1])
    % pause
    title(sprintf('Bin #%d (%.1f trials)', iBin, mean(nData_allB(:, iBin))))
end

sgtitle(sprintf('[pYES] NOM ModelA%dB%d [nblocks = %d] [ORI%d SF%d] [niter = %d]', ...
    iModelA, iModelB, nblocks, nORI, nSF, ni))

%% pC: data vs. pred
figure('Position', [100 100 2e3 1e3])
hold on
for iBin=1:nBins
    subplot(2, nBins/2, iBin), hold on
    plot(pC_data_allB(:, iBin), pC_pred_allB(:, iBin), 'o')

    % plot ave and sd on top of the scatter plot
    errorbar(nanmean(pC_data_allB(:, iBin)), nanmean(pC_pred_allB(:, iBin)), nanstd(pC_data_allB(:, iBin)), 'k', 'horizontal', 'CapSize', 0, 'LineWidth', 2)
    errorbar(nanmean(pC_data_allB(:, iBin)), nanmean(pC_pred_allB(:, iBin)), nanstd(pC_pred_allB(:, iBin)), 'k', 'vertical', 'CapSize', 0, 'LineWidth', 2)

    xlabel('Emp pC')
    ylabel('Pred pC')
    plot([0, 1], [0, 1], 'k--')
    xline(.5, 'r--');
    yline(.5, 'r--');
    axis square
    xlim([0,1])
    ylim([0,1])
    % pause
    title(sprintf('Bin #%d (%.1f trials)', iBin, mean(nData_allB(:, iBin))))
end

sgtitle(sprintf('[pC] NOM ModelA%dB%d [nblocks = %d] [ORI%d SF%d] [niter = %d]', ...
    iModelA, iModelB, nblocks, nORI, nSF, ni))

%% pA: data vs. pred
figure('Position', [100 100 2e3 1e3])
hold on
for iBin=1:nBins
    subplot(2, nBins/2, iBin), hold on
    plot(pA_data_allB(:, iBin), pA_pred_allB(:, iBin), 'o')

    % plot ave and sd on top of the scatter plot
    errorbar(nanmean(pA_data_allB(:, iBin)), nanmean(pA_pred_allB(:, iBin)), nanstd(pA_data_allB(:, iBin)), 'k', 'horizontal', 'CapSize', 0, 'LineWidth', 2)
    errorbar(nanmean(pA_data_allB(:, iBin)), nanmean(pA_pred_allB(:, iBin)), nanstd(pA_pred_allB(:, iBin)), 'k', 'vertical', 'CapSize', 0, 'LineWidth', 2)

    xlabel('Emp pA')
    ylabel('Pred pA')
    plot([0, 1], [0, 1], 'k--')
    axis square
    xline(.5, 'r--');
    yline(.5, 'r--');
    xlim([0,1])
    ylim([0,1])
    % pause
    title(sprintf('Bin #%d (%.1f trials)', iBin, mean(nData_allB(:, iBin))))
end

sgtitle(sprintf('[pA] NOM ModelA%dB%d [nblocks = %d] [ORI%d SF%d] [niter = %d]', ...
    iModelA, iModelB, nblocks, nORI, nSF, ni))

%% Plot estimated parameters across iterations
if ismember(iModelB, [4,5]) % so that 1st param is noise P, 2nd param is criterion
    figure('Name','Estimated Parameters','Position',[100,100,1200,400]);
    % for iParam = 1:size(params_est_allB,2)

    iParam = 1;
    subplot(1, 2, iParam);
    plot(params_est_allB(:,iParam), '-o');
    yline(.1, 'r-'); % match noiseP defined in OOD_sim
    xlabel('Iteration');
    ylabel(namesParamsModel_all{iModelB}{iParam});
    title(['Parameter: ', namesParamsModel_all{iModelB}{iParam}]);

    % iParam = 2;
    % subplot(1, 2, iParam);
    % plot(params_est_allB(:,iParam), '-o');
    % yline(criterion_true, 'r-'); % criterion_true
    % xlabel('Iteration');
    % ylabel(namesParamsModel_all{iModelB}{iParam});
    % title(['Parameter: ', namesParamsModel_all{iModelB}{iParam}]);
    % end
    sgtitle('Estimated Parameters Across Iterations');
end