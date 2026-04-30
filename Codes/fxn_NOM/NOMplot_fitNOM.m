
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
Cz_emp_allIter = nan(nIter, 1);

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

    Cz_emp_allIter(iIter) = pred_metrics_allIter{iIter}.Cz_emp;
end

% fprintf('\n *** Compiling DONE, READY to plot ***\n ')

%% Figure 1: pYES/pC/pA as a fxn of binned DV
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

saveas(gcf, sprintf('%s/21Metrics_L%d_A%dB%d.jpg', nameFolder_Figures_perSubj, iLocComb, iModelA_fit, iModelB_fit))

%% Figure 2: Plot estimated parameters across iterations
nParams_full = 5; % Nmul, Nadd, Nshared, criterion_DV, Empirical Criterion SDT
params_est_allIter = [params_est_allIter, Cz_emp_allIter];

% Load truth once if this is a simulation dataset
hasTruth = numel(subjName) > 10;
if hasTruth
    S_truth = load(fullfile(nameFolder_OOD_load, 'truth.mat'), '*_true');
end

% Define fixed 5-column layout
% Columns: 1=Nmul, 2=Nadd, 3=Nshared, 4=criterion_DV, 5=Empirical Criterion SDT
NOMc_ub = 20;
switch iModelB_fit
    case 1 % Full: Nmul, Nadd, Nshared, criterion_DV
        col_idx      = [1 2 3 4 5];
        param_labels = {'Nmul', 'Nadd', 'Nshared', 'Empirical criterion DV', 'Empirical criterion SDT'};
        ylim_list    = {[NOMp1_lb NOMp1_ub], [NOMp2_lb NOMp2_ub], [NOMp3_lb NOMp3_ub], [0 NOMc_ub], [-1 1]};
        true_vals    = nan(1, 5);
        if hasTruth
            if isfield(S_truth, 'Nmul_true'), true_vals(1) = S_truth.Nmul_true; end
            if isfield(S_truth, 'Nadd_true'), true_vals(2) = S_truth.Nadd_true; end
            if isfield(S_truth, 'Nshared_true'), true_vals(3) = S_truth.Nshared_true; end
            if isfield(S_truth, 'criterion_DV_true'), true_vals(4) = S_truth.criterion_DV_true; end
            if isfield(S_truth, 'cSDT_true'), true_vals(5) = S_truth.cSDT_true; end
        end

    case 2 % No Nmul: Nadd, Nshared, criterion_DV
        col_idx      = [2 3 4 5];
        param_labels = {'Nadd', 'Nshared', 'Empirical Criterion DV', 'Empirical Criterion SDT'};
        ylim_list    = {[NOMp2_lb NOMp2_ub], [NOMp3_lb NOMp3_ub], [0 NOMc_ub], [-1 1]};
        true_vals    = nan(1, 4);
        if hasTruth
            if isfield(S_truth, 'Nadd_true'), true_vals(1) = S_truth.Nadd_true; end
            if isfield(S_truth, 'Nshared_true'), true_vals(2) = S_truth.Nshared_true; end
            if isfield(S_truth, 'criterion_DV_true'), true_vals(3) = S_truth.criterion_DV_true; end
            if isfield(S_truth, 'cSDT_true'), true_vals(4) = S_truth.cSDT_true; end
        end

    case 3 % No Nadd: Nmul, Nshared, criterion_DV
        col_idx      = [1 3 4 5];
        param_labels = {'Nmul', 'Nshared', 'Empirical Criterion DV', 'Empirical Criterion SDT'};
        ylim_list    = {[NOMp1_lb NOMp1_ub], [NOMp3_lb NOMp3_ub], [0 NOMc_ub], [-1 1]};
        true_vals    = nan(1, 4);
        if hasTruth
            if isfield(S_truth, 'Nmul_true'), true_vals(1) = S_truth.Nmul_true; end
            if isfield(S_truth, 'Nshared_true'), true_vals(2) = S_truth.Nshared_true; end
            if isfield(S_truth, 'criterion_DV_true'), true_vals(3) = S_truth.criterion_DV_true; end
            if isfield(S_truth, 'cSDT_true'), true_vals(4) = S_truth.cSDT_true; end
        end

    case 4 % No Nshared: Nmul, Nadd, criterion_DV
        col_idx      = [1 2 4 5];
        param_labels = {'Nmul', 'Nadd', 'Empirical Criterion DV', 'Empirical Criterion SDT'};
        ylim_list    = {[NOMp1_lb NOMp1_ub], [NOMp2_lb NOMp2_ub], [0 NOMc_ub], [-1 1]};
        true_vals    = nan(1, 4);
        if hasTruth
            if isfield(S_truth, 'Nmul_true'), true_vals(1) = S_truth.Nmul_true; end
            if isfield(S_truth, 'Nadd_true'), true_vals(2) = S_truth.Nadd_true; end
            if isfield(S_truth, 'criterion_DV_true'), true_vals(3) = S_truth.criterion_DV_true; end
            if isfield(S_truth, 'cSDT_true'), true_vals(4) = S_truth.cSDT_true; end
        end

    case 5 % Nmul only: Nmul, criterion_DV
        col_idx      = [1 4 5];
        param_labels = {'Nmul', 'Empirical Criterion DV', 'Empirical Criterion SDT'};
        ylim_list    = {[NOMp1_lb NOMp1_ub], [0 NOMc_ub], [-1 1]};
        true_vals    = nan(1, 3);
        if hasTruth
            if isfield(S_truth, 'Nmul_true'), true_vals(1) = S_truth.Nmul_true; end
            if isfield(S_truth, 'criterion_DV_true'), true_vals(2) = S_truth.criterion_DV_true; end
            if isfield(S_truth, 'cSDT_true'), true_vals(3) = S_truth.cSDT_true; end
        end

    case 6 % Nadd only: Nadd, criterion_DV
        col_idx      = [2 4 5];
        param_labels = {'Nadd', 'Empirical Criterion DV', 'Empirical Criterion SDT'};
        ylim_list    = {[NOMp2_lb NOMp2_ub], [0 NOMc_ub], [-1 1]};
        true_vals    = nan(1, 3);
        if hasTruth
            if isfield(S_truth, 'Nadd_true'), true_vals(1) = S_truth.Nadd_true; end
            if isfield(S_truth, 'criterion_DV_true'), true_vals(2) = S_truth.criterion_DV_true; end
            if isfield(S_truth, 'cSDT_true'), true_vals(3) = S_truth.cSDT_true; end
        end

    case 7 % Nshared only: Nshared, criterion_DV
        col_idx      = [3 4 5];
        param_labels = {'Nshared', 'Empirical Criterion DV', 'Empirical Criterion SDT'};
        ylim_list    = {[NOMp3_lb NOMp3_ub], [0 NOMc_ub], [-1 1]};
        true_vals    = nan(1, 3);
        if hasTruth
            if isfield(S_truth, 'Nshared_true'), true_vals(1) = S_truth.Nshared_true; end
            if isfield(S_truth, 'criterion_DV_true'), true_vals(2) = S_truth.criterion_DV_true; end
            if isfield(S_truth, 'cSDT_true'), true_vals(3) = S_truth.cSDT_true; end
        end

    case 8 % Criterion only: criterion_DV
        col_idx      = [4 5];
        param_labels = {'Empirical Criterion DV', 'Empirical Criterion SDT'};
        ylim_list    = {[0 NOMc_ub], [-1 1]};
        true_vals    = nan(1, 2);
        if hasTruth
            if isfield(S_truth, 'criterion_DV_true'), true_vals(1) = S_truth.criterion_DV_true; end
            if isfield(S_truth, 'cSDT_true'), true_vals(2) = S_truth.cSDT_true; end
        end

    otherwise
        error('Unknown iModelB_fit = %d', iModelB_fit);
end

figure('Position', [100, 100, nParams_full * 400, 400]);

for iParam = 1:size(params_est_allIter, 2)

    thisCol   = col_idx(iParam);
    thisLabel = param_labels{iParam};
    thisYLim  = ylim_list{iParam};
    thisVals  = params_est_allIter(:, iParam);

    subplot(1, nParams_full, thisCol); hold on;

    % Plot true value if available
    if hasTruth && ~isnan(true_vals(iParam))
        yline(true_vals(iParam), 'r--', 'LineWidth', 2, 'DisplayName', 'True value');
    end

    % Plot estimate per iteration
    plot(thisVals, 'ko', 'LineWidth', 1, 'DisplayName', 'Estimates');

    % Plot median across iterations
    medVal = nanmedian(thisVals);
    SDVal = std(thisVals);
    yline(medVal, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Med. est.');

    ylim(thisYLim);
    xlabel('Iteration');
    ylabel(thisLabel);
    if hasTruth && ~isnan(true_vals(iParam))
        title(sprintf('Parameter: %s \nTrue: %.2f\nEst.: %.2f (SD=%.2f)', thisLabel, true_vals(iParam), medVal, SDVal));
    else
        title(sprintf('Parameter: %s\nEst.: %.2f (SD=%.2f)', thisLabel, medVal, SDVal));
    end
    box on;

    if iParam==1
        legend('show', 'Location', 'best')
    end
end

sgtitle(sprintf('Figure 2. Estimated parameters across iterations\n%s (ModelA%dB%d %s, L%d, nIter=%d)', ...
    subjName, iModelA_fit, iModelB_fit, namesModelB{iModelB_fit}, iLocComb, nIter));

saveas(gcf, sprintf('%s/22Params_L%d_A%dB%d.jpg', ...
    nameFolder_Figures_perSubj, iLocComb, iModelA_fit, iModelB_fit));

%% Figure 3: Correlation between paired parameters across iterations
params_plot = params_est_allIter(:, 1:numel(param_labels));
nParams_plot = size(params_plot, 2);

% 68% CI summary for each parameter across iterations
param_med = nan(1, nParams_plot);
param_lb68 = nan(1, nParams_plot);
param_ub68 = nan(1, nParams_plot);
for iParam = 1:nParams_plot
    [param_med(iParam), param_lb68(iParam), param_ub68(iParam)] = getCI(params_plot(:, iParam), 1, 1);
end

figure('Position', [100, 100, 320*nParams_plot, 320*nParams_plot]);

for iRow = 1:nParams_plot
    for iCol = 1:nParams_plot
        if iCol >= iRow
        subplot(nParams_plot, nParams_plot, (iRow-1)*nParams_plot + iCol); hold on;

        x = params_plot(:, iCol);
        y = params_plot(:, iRow);
        valid = ~isnan(x) & ~isnan(y);
        x = x(valid);
        y = y(valid);

        if iRow == iCol
            % diagonal: show 1D distribution across iterations
            histogram(x, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
            xline(param_med(iCol), 'k-', 'LineWidth', 1.2);
            xline(param_lb68(iCol), 'k--', 'LineWidth', 1);
            xline(param_ub68(iCol), 'k--', 'LineWidth', 1);
            if hasTruth && ~isnan(true_vals(iCol))
                xline(true_vals(iCol), 'r-', 'LineWidth', 1.2);
            end
            xlabel(param_labels{iCol});
            ylabel('Count');
            title(param_labels{iCol}, 'Interpreter', 'none');
        else
            % scatter across iterations
            plot(x, y, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', [0.65 0.65 0.65], ...
                'MarkerEdgeColor', [0.65 0.65 0.65]);

            % median point with 68% CI error bars
            errorbar(param_med(iCol), param_med(iRow), ...
                param_med(iRow)-param_lb68(iRow), param_ub68(iRow)-param_med(iRow), ...
                param_med(iCol)-param_lb68(iCol), param_ub68(iCol)-param_med(iCol), ...
                'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 0);
            plot(param_med(iCol), param_med(iRow), 'kx', 'MarkerSize', 12, 'LineWidth', 2);

            % true value cross for simulations
            if hasTruth && ~isnan(true_vals(iCol)) && ~isnan(true_vals(iRow))
                plot(true_vals(iCol), true_vals(iRow), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
            end

            % correlation across iterations
            if numel(x) >= 3 && std(x) > 0 && std(y) > 0
                r_xy = corr(x, y, 'rows', 'complete');
                title(sprintf('r = %.2f', r_xy));
            end

            xlabel(param_labels{iCol}, 'Interpreter', 'none');
            ylabel(param_labels{iRow}, 'Interpreter', 'none');

            % if true_vals(iRow)<=0
            %     ylim([-1,1])
            % elseif true_vals(iRow)<=1
            %     ylim([0,1])
            % else
            %     ylim([0,20])
            % end
            % 
            % if true_vals(iCol)<=0
            %     xlim([-1,1])
            % elseif true_vals(iCol)<=1
            %     xlim([0,1])
            % else
            %     xlim([0,20])
            % end

        end

        box on; axis square
        end
    end
end

sgtitle(sprintf(['Figure 3. Pairwise correlation between parameter estimates across iterations\n' ...
    '%s (ModelA%dB%d %s, L%d, nIter=%d)'], ...
    subjName, iModelA_fit, iModelB_fit, namesModelB{iModelB_fit}, iLocComb, nIter));

saveas(gcf, sprintf('%s/23ParamCorr_L%d_A%dB%d.jpg', ...
    nameFolder_Figures_perSubj, iLocComb, iModelA_fit, iModelB_fit));

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
% saveas(gcf, sprintf('%s/24FreezeCorr_L%d_A%dB%d.jpg', nameFolder_Figures_perSubj, iLocComb, iModelA_fit, iModelB_fit))
