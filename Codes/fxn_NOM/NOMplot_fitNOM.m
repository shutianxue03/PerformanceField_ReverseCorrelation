
% This script generates plots for the estimated parameters and prediction metrics of the Noisy Observer Model (NOM).


%% Compile data and pred for all iterations and bins
set(0, 'DefaultFigureVisible', 'off') % avoid printing figures on the desktop

pred = pred_test;

pYES_data_allBins = nan(nIter, nBins);
pYES_pred_allBins = pYES_data_allBins;
pC_data_allBins = pYES_data_allBins;
pC_pred_allBins = pYES_data_allBins;
pA_data_allBins = pYES_data_allBins;
pA_pred_allBins = pYES_data_allBins;
nData_allBins = pYES_data_allBins;
DV_allBins_all = pYES_data_allBins;

for iIter = 1:nIter

    DV_allBins_all(iIter, :) = pred.metrics.IV_allBins;

    % dprime_data_allB(ii, :) = pred_metrics_allIter{ii}.metrics.dprime_data_allBins;
    % dprime_pred_allB(ii, :) = pred_metrics_allIter{ii}.metrics.dprime_pred_allBins;

    pYES_data_allBins(iIter, :) = pred_metrics_allIter{iIter}.metrics.pYES_data_allBins;
    pYES_pred_allBins(iIter, :) = pred_metrics_allIter{iIter}.metrics.pYES_pred_allBins;

    pC_data_allBins(iIter, :) = pred_metrics_allIter{iIter}.metrics.pC_data_allBins;
    pC_pred_allBins(iIter, :) = pred_metrics_allIter{iIter}.metrics.pC_pred_allBins;

    pA_data_allBins(iIter, :) = pred_metrics_allIter{iIter}.metrics.pA_data_allBins;
    pA_pred_allBins(iIter, :) = pred_metrics_allIter{iIter}.metrics.pA_pred_allBins;

    nData_allBins(iIter, :) = data_test_allIter{iIter}.nTrials_allBins;
end

% fprintf('\n *** Compiling DONE, READY to plot ***\n ')

%% 1. pYES/pC/pA as a fxn of binned DV
[DV_allBins_med] = getCI(DV_allBins_all, 1, 1);
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

for iBin = 1:nBins
    % Data (dot+errorbars)
    plot(DV_allBins_med(iBin), pYES_data_med(iBin), 'ko', 'MarkerSize', nData_allB_med(iBin)/sz_scale+5)
    errorbar(DV_allBins_med(iBin), pYES_data_med(iBin), pYES_data_SEM_neg(iBin), pYES_data_SEM_pos(iBin), 'k', 'CapSize', 0)
end
% Pred (shaded errorbars)
plot(DV_allBins_med, pYES_pred_med, 'k-')
patch([DV_allBins_med, fliplr(DV_allBins_med)], [pYES_pred_lb, fliplr(pYES_pred_ub)], 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');

yline(.5, 'k--');
ylim([0,1])
xlim([0, 30])
xlabel('Binned DV')
ylabel('pYES')
% add r and R-squared
r_pYES = corr(DV_allBins_med', pYES_data_med');
R2_pYES = 1 - sum((pYES_data_med - pYES_pred_med).^2) / sum((pYES_data_med - mean(pYES_data_med)).^2);
title(sprintf('pYES (r=%.2f, R²=%.2f)', r_pYES, R2_pYES));

subplot(3,1,2), hold on
for iBin = 1:nBins
    % Data (dot+errorbars)
    plot(DV_allBins_med(iBin), pC_data_med(iBin), 'ko', 'MarkerSize', nData_allB_med(iBin)/sz_scale+5)
    errorbar(DV_allBins_med(iBin), pC_data_med(iBin), pC_data_SEM_neg(iBin), pC_data_SEM_pos(iBin), 'k', 'CapSize', 0)
end
% Pred (shaded errorbars)
plot(DV_allBins_med, pC_pred_med, 'k-')
patch([DV_allBins_med, fliplr(DV_allBins_med)], [pC_pred_lb, fliplr(pC_pred_ub)], 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');

yline(.5, 'k--');
yline(mean(pC_data_med), 'k-');
ylim([0,1])
xlim([0, 30])
xlabel('Binned DV')
ylabel('pC')
% add r and R-squared
r_pC = corr(pC_pred_med', pC_data_med');
R2_pC = 1 - sum((pC_data_med - pC_pred_med).^2) / sum((pC_data_med - mean(pC_data_med)).^2);
title(sprintf('pC (r=%.2f, R²=%.2f)', r_pC, R2_pC));

subplot(3,1,3), hold on
for iBin = 1:nBins
    % Data (dot+errorbars)
    plot(DV_allBins_med(iBin), pA_data_med(iBin), 'ko', 'MarkerSize', nData_allB_med(iBin)/sz_scale+5)
    errorbar(DV_allBins_med(iBin), pA_data_med(iBin), pA_data_SEM_neg(iBin), pA_data_SEM_pos(iBin), 'k', 'CapSize', 0)
    % Plot predictied pA from measured pYES
    % plot(DV_allBins_med(iBin), pYES_data_med(iBin)^2+(1-pYES_data_med(iBin))^2, 'co', 'MarkerSize', nData_allB_med(iBin)/sz_scale+5)
end
% Pred (shaded errorbars)
plot(DV_allBins_med, pA_pred_med, 'k-');
patch([DV_allBins_med, fliplr(DV_allBins_med)], [pA_pred_lb, fliplr(pA_pred_ub)], 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');

ylim([.5,1])
xlim([0, 30])
xlabel('Binned DV')
ylabel('pA')
% add r and R-squared
r_pA = corr(pA_pred_med', pA_data_med');
R2_pA = 1 - sum((pA_data_med - pA_pred_med).^2) / sum((pA_data_med - mean(pA_data_med)).^2);
title(sprintf('pA (r=%.2f, R²=%.2f)', r_pA, R2_pA));

sgtitle(sprintf('Figure 1. Metrics vs. binned DV\n%s\nModelA%dB%d %s, L%d, nIter=%d', ...
    subjName, iModelA_fit, iModelB_fit, namesModelB{iModelB_fit}, iLocComb, nIter), ...
    'fontsize', 10)

saveas(gcf, sprintf('%s/4Metrics_L%d_A%dB%d.jpg', nameFolder_Figures_perSubj, iLocComb, iModelA_fit, iModelB_fit))

%% 2. Plot estimated parameters across iterations
nParams_full = 4; % Nmul, Nadd, Nshared, criterion_DV

figure('Position', [100,100,nParams_full*400,400]);

for iParam = 1:nParams
    % Load and plot the true param
   if numel(subjName) > 10
        load(sprintf('%s/truth.mat', nameFolder_OOD_load), '*_true')
        switch iModelB_fit
            case 1 % full
                subplot(1, nParams_full, iParam); hold on
                switch iParam
                    case 1, val_yline = Nmul_true;
                    case 2, val_yline = Nadd_true;
                    case 3, val_yline = Nshared_true;
                    case 4, val_yline = criterion_DV_true;
                end
            case 2 % No Nshared
                switch iParam
                    case 1, val_yline = Nmul_true; subplot(1, nParams_full, iParam); hold on
                    case 2, val_yline = Nadd_true; subplot(1, nParams_full, iParam); hold on
                    case 3, val_yline = criterion_DV_true; subplot(1, nParams_full, iParam+1); hold on
                end 
            case 3 % No Nmul
                subplot(1, nParams_full, iParam+1); hold on
                switch iParam
                    case 1, val_yline = Nadd_true;
                    case 2, val_yline = Nshared_true;
                    case 3, val_yline = criterion_DV_true;
                end
            case 4 % No Nadd
                switch iParam
                    case 1, val_yline = Nmul_true; subplot(1, nParams_full, 1); hold on
                    case 2, val_yline = Nshared_true; subplot(1, nParams_full, iParam+1); hold on
                    case 3, val_yline = criterion_DV_true; subplot(1, nParams_full, iParam+1); hold on
                end
        end
        yline(val_yline, 'r-');
    end

    % Plot estimate per iteration
    plot(params_est_allIter(:, iParam), '-o');

    % Plot the average (should be changed to median later)
    yline(nanmean(params_est_allIter(:, iParam)), 'k-')

    ylim([params_lb(iParam), params_ub(iParam)])
    xlabel('Iteration');
    ylabel(namesModelBparams{iModelB_fit}{iParam});
    title(sprintf('Parameter: %s\n %.2f', namesModelBparams{iModelB_fit}{iParam}, nanmean(params_est_allIter(:, iParam))));

    % Set ylimit
    switch iModelB_fit
        case 1 % full model: Nmul, Nadd, Nshared, criterion in DV
            switch iParam
                case 1, ylim([0, 1])
                case 2, ylim([0, 50])
                case 3, ylim([0, 50])
                case 4, ylim([0, 50])
            end
        case 2 % No Nshared
            switch iParam
                case 1, ylim([0, 1])
                case 2, ylim([0, 50])
                case 3, ylim([0, 50])
            end
        case 3 % No Nmul
            switch iParam
                case 1, ylim([0, 50])
                case 2, ylim([0, 50])
                case 3, ylim([0, 50])
            end
        case 4 % No Nadd
            switch iParam
                case 1, ylim([0, 1])
                case 2, ylim([0, 50])
                case 3, ylim([0, 50])
            end
    end
end

sgtitle(sprintf('Figure 2. Estimated parameters across iterations\n%s (ModelA%dB%d %s, L%d, nIter=%d)', ...
    subjName, iModelA_fit, iModelB_fit, namesModelB{iModelB_fit}, iLocComb, nIter))

saveas(gcf, sprintf('%s/5Params_L%d_A%dB%d.jpg', nameFolder_Figures_perSubj, iLocComb, iModelA_fit, iModelB_fit))

%% 3. Freeze other params and vary one param to see its corr with pA
% figure('Position', [0 0 1e3 800])
%
% for iParamFreeze = 1:nParams
%     NOM_grid = linspace(params_lb(iParamFreeze), params_ub(iParamFreeze), 20);
%     pYES_mean = nan(size(NOM_grid));
%     pA_mean = pYES_mean;
%
%     for i = 1:numel(NOM_grid)
%         NOM_test = NOM_grid(i);
%
%         switch nParams
%             case 4
%                 switch iParamFreeze
%                     case 1, params_test  = [NOM_test, params_est(2), params_est(3), params_est(4)]; % iModelB_fit = 1
%                     case 2, params_test  = [params_est(1), NOM_test, params_est(3), params_est(4)]; % iModelB_fit = 1
%                     case 3, params_test  = [params_est(1), params_est(2), NOM_test, params_est(4)]; % iModelB_fit = 1
%                     case 4, params_test  = [params_est(1), params_est(2), params_est(3), NOM_test]; % iModelB_fit = 1
%                 end
%             case 3
%                 switch iParamFreeze
%                     case 1, params_test  = [NOM_test, params_est(2), params_est(3)]; % iModelB_fit = 1
%                     case 2, params_test  = [params_est(1), NOM_test, params_est(3)]; % iModelB_fit = 1
%                     case 3, params_test  = [params_est(1), params_est(2), NOM_test]; % iModelB_fit = 1
%                 end
%             case 2
%                 switch iParamFreeze
%                     case 1, params_test  = [NOM_test, params_est(2)]; % iModelB_fit = 2
%                     case 2, params_test  = [params_est(1), NOM_test]; % iModelB_fit = 2
%                 end
%         end
%
%         % Make prediction
%         switch flag_incluCrit
%             case 0
%                 %---------------------%
%                 [pYES_pred_allT, pA_pred_allPairs, ~] = fxn_predMetrics_v2(iModelB_fit, params_test, data_test, criterion_z_train);
%                 %---------------------%
%             case 1
%                 %---------------------%
%                 [pYES_pred_allT, pA_pred_allPairs, ~, Cz_fit] = fxn_predMetrics_v3(iModelB_fit, params_test, data_test);
%                 %---------------------%
%         end
%         pYES_mean(i) = mean(pYES_pred_allT);
%         pA_mean(i) = mean(pA_pred_allPairs);
%
%         subplot(2,nParams, iParamFreeze)
%         plot(NOM_grid, pYES_mean, '-o');
%         ylim([0, 1])
%         xlabel(sprintf('Only Vary %s', namesModelBparams{iModelB_fit}{iParamFreeze}));
%         ylabel('Mean pYES_{pred}');
%
%         subplot(2,nParams, iParamFreeze+nParams)
%         plot(NOM_grid, pA_mean, '-o');
%         ylim([.5, 1])
%         xlabel(sprintf('Only Vary %s', namesModelBparams{iModelB_fit}{iParamFreeze}));
%         ylabel('Mean pA_{pred}');
%     end
% end
%
% sgtitle(sprintf('Figure 3. Vary one param (and freeze others) to see its corr with pA\n%s (ModelA%dB%d %s, L%d, nIter=%d)', ...
%     subjName, iModelA_fit, iModelB_fit, namesModelB{iModelB_fit}, iLocComb, nIter))
%
% saveas(gcf, sprintf('%s/6FreezeCorr_L%d_A%dB%d.jpg', nameFolder_Figures_perSubj, iLocComb, iModelA_fit, iModelB_fit))
