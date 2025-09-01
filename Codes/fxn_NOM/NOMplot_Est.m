
% This script generates plots for the estimated parameters and prediction metrics of the Noisy Observer Model (NOM).

%% Compile data and pred for all iteractions and bins
pC_data_allB = nan(ni, nBins);
pC_pred_allB = pC_data_allB;
pA_data_allB = pC_data_allB;
pA_pred_allB = pC_data_allB;
nData_allB = pC_data_allB;

indPc = 3;
indPa = 6;

for ii=1:ni
    pC_data_allB(ii, :) = data_allB{ii}.metrics_sim(indPc);
    pC_pred_allB(ii, :) = pred_metrics_allB{ii}.metrics.pC_pred_allBins;

    pA_data_allB(ii, :) = data_allB{ii}.metrics_sim(indPa);
    pA_pred_allB(ii, :) = pred_metrics_allB{ii}.metrics.pA_pred_allBins;

    nData_allB(ii, :) = data_allB{ii}.nTrials_allBins;
end

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
    xlim([0,1])
    ylim([0,1])
    % pause
    title(sprintf('Bin #%d', iBin))
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
    yline(.5, 'k-')
    xlabel('Iteration');
    ylabel(namesParamsModel_all{iModelB}{iParam});
    title(['Parameter: ', namesParamsModel_all{iModelB}{iParam}]);

    iParam = 2;
    subplot(1, 2, iParam);
    plot(params_est_allB(:,iParam), '-o');
    yline(-.1, 'k-')
    xlabel('Iteration');
    ylabel(namesParamsModel_all{iModelB}{iParam});
    title(['Parameter: ', namesParamsModel_all{iModelB}{iParam}]);
    % end
    sgtitle('Estimated Parameters Across Iterations');
end