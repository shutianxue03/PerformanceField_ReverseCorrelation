%% ============================================================
% simPlot_part2_TemplateRecovery.m
%
% Part 2. Ideal observer with a data-derived template
% Efficient version: summarize each file immediately
%
% Created by Shutian Xue
%% ============================================================

clear; clc; close all;
set(0, 'DefaultFigureVisible', 'off');

%--------------%
SX_RC1_setting;
%--------------%

%% ---------------- settings ----------------
basis_comb = 4:9;
basis_comb = allcomb(basis_comb, basis_comb);
basis_comb = [6,5;6,6;7,5;7,6;8,5;8,6;9,5;9,6];
basis_comb = [6,5];
nComb = size(basis_comb, 1);

marker_list = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h', 'x', '+'};

for iComb=1:nComb
    str_part = sprintf('Part3_nBasis%d%d_noMirror', basis_comb(iComb, :));
    iModelA_fit = 2; % derived template
    iModelB_fit_all = 1:4; % inspect all fitted models

    nameFolder_Data = sprintf('%s/Data_%s', nameFolder_server, str_part);
    nameFolder_Data_OOD = sprintf('%s/Data_OOD_%d%d', nameFolder_Data, nORI, nSF);
    nameFolder_Data_NOM_Trialwise = sprintf('%s/Data_NOM_Trialwise_%d%d', nameFolder_Data, nORI, nSF);

    nameFolder_Figures_part = fullfile(nameFolder_Figures, sprintf('IO_%s', str_part));
    if ~exist(nameFolder_Figures_part, 'dir')
        mkdir(nameFolder_Figures_part);
    end

    metricNames = {'pYES','pC','pA'};

    %% ---------------- build file table once ----------------
    S = build_file_table(nameFolder_Data_OOD, nameFolder_Data_NOM_Trialwise, iModelA_fit);

    if isempty(S)
        error('No valid simulation folders found in %s', nameFolder_Data_OOD);
    end

    signal_unik = unique([S.signalCST]);
    Cz_unik = unique([S.Cz]);
    lambda_unik = unique([S.lambda]);
    cmap_Cz = parula(numel(Cz_unik));

    %% loop over fitted models

    for iModelB_fit = iModelB_fit_all

        fprintf('\n============================================\n');
        fprintf('Part 2: nBasis (ORI=%d, SF=%d) | A=%d, B=%d (%s)\n', basis_comb(iComb, :), iModelA_fit, iModelB_fit, namesModelB{iModelB_fit});
        fprintf('============================================\n');

        nFiles = numel(S);

        % determine parameter names for this fitted model
        if iModelB_fit <= numel(namesModelBparams)
            param_names_all = namesModelBparams{iModelB_fit};
        else
            param_names_all = {};
        end

        % infer number of IN params from model name list:
        % last one should be criterion
        nParams = max(numel(param_names_all) - 1, 1);
        if isempty(param_names_all)
            param_names_all = [arrayfun(@(x) sprintf('Param%d', x), 1:nParams, 'UniformOutput', false), {'criterion'}];
        end
        IN_param_names = param_names_all(1:nParams);

        %% ---------------- summary storage ----------------
        signal_plot = nan(nFiles,1);
        Cz_plot = nan(nFiles,1);
        lambda_plot = nan(nFiles,1);
        keepFile = false(nFiles,1);
        
        % template
        template_rmse_mean = nan(nFiles,1); template_rmse_lb = nan(nFiles,1); template_rmse_ub = nan(nFiles,1);
        template_R2_mean = nan(nFiles,1); template_R2_lb = nan(nFiles,1); template_R2_ub = nan(nFiles,1);
        template_true_plot = nan(nFiles,1);
        template_est_plot = nan(nFiles,1); template_est_lb_plot = nan(nFiles,1); template_est_ub_plot = nan(nFiles,1);

        % Tuning fuctions
        margORI_true_plot = nan(nFiles, nORI);
        margORI_est_plot  = nan(nFiles, nORI);
        margORI_lb_plot   = nan(nFiles, nORI);
        margORI_ub_plot   = nan(nFiles, nORI);

        margSF_true_plot = nan(nFiles, nSF);
        margSF_est_plot  = nan(nFiles, nSF);
        margSF_lb_plot   = nan(nFiles, nSF);
        margSF_ub_plot   = nan(nFiles, nSF);

        % metrics
        for iM = 1:numel(metricNames)
            m = metricNames{iM};
            M.(m).rmse_mean = nan(nFiles,1);
            M.(m).rmse_lb = nan(nFiles,1);
            M.(m).rmse_ub = nan(nFiles,1);
            M.(m).R2_mean = nan(nFiles,1);
            M.(m).R2_lb = nan(nFiles,1);
            M.(m).R2_ub = nan(nFiles,1);

            M.(m).data_plot = nan(nFiles,1);
            M.(m).pred_plot = nan(nFiles,1);
            M.(m).pred_lb_plot = nan(nFiles,1);
            M.(m).pred_ub_plot = nan(nFiles,1);
        end

        % criterion
        criterion_rmse_mean = nan(nFiles,1); criterion_rmse_lb = nan(nFiles,1); criterion_rmse_ub = nan(nFiles,1);
        criterion_true_plot = nan(nFiles,1);
        criterion_est_plot = nan(nFiles,1); criterion_est_lb_plot = nan(nFiles,1); criterion_est_ub_plot = nan(nFiles,1);

        % IN params
        IN_rmse_mean = nan(nFiles,nParams); IN_rmse_lb = nan(nFiles,nParams); IN_rmse_ub = nan(nFiles,nParams);
        IN_true_plot = zeros(nFiles,nParams);
        IN_est_plot = nan(nFiles,nParams); IN_est_lb_plot = nan(nFiles,nParams); IN_est_ub_plot = nan(nFiles,nParams);

        % Full iteration-level estimated parameters for each file
        P_allFiles = cell(nFiles, 1);

        %% ---------------- compile per file ----------------
        fprintf('\n%s: Processing %d files.\n', datetime('now'), nFiles)
        for iFile = 1:nFiles
            fprintf('%d ', iFile)
            % fprintf('%s: %d/%d %s | B=%d\n', datetime('now'), iFile, nFiles, S(iFile).nameIO, iModelB_fit);

            fitPat = sprintf('n*_A%dB%d.mat', iModelA_fit, iModelB_fit);
            fitDir = dir(fullfile(S(iFile).nomFolder, fitPat));
            fitDir = fitDir(~contains({fitDir.name}, 'min'));
            if isempty(fitDir)
                continue;
            end
            fitFile = fullfile(fitDir(1).folder, fitDir(1).name);

            truth = load(S(iFile).truthFile);
            data_compIV = load(S(iFile).compIVFile);
            data_fitNOM = load(fitFile);

            nIter = truth.nIter;

            signal_plot(iFile) = S(iFile).signalCST;
            Cz_plot(iFile) = S(iFile).Cz;
            lambda_plot(iFile) = S(iFile).lambda;

            % ---------- template ----------
            % True template
            template_true = truth.template_true(:)';
            % Normalize true template
            template_true = template_true/max(template_true);
            template_true_plot(iFile) = mean(template_true);

            % Estimated template
            template_est = data_compIV.template_tmpl_allIter;

            if ndims(template_est) == 3
                template_est = reshape(template_est, size(template_est,1), []);
            end

            % Normalize
            template_est = template_est./max(template_est, [], 2);

            template_est_mean_iter = mean(template_est, 2);
            [template_est_plot(iFile), template_est_lb_plot(iFile), template_est_ub_plot(iFile)] = ...
                getCI(template_est_mean_iter, 1, 1);

            % Get RMSE
            template_rmse_iter = sqrt(mean((template_est - template_true).^2, 2));
            [template_rmse_mean(iFile), template_rmse_lb(iFile), template_rmse_ub(iFile)] = getCI(template_rmse_iter, 1, 1);

            % Get R2
            template_R2_iter = nan(nIter,1);
            yTrue = template_true(:);
            for iIter = 1:nIter
                yEst = template_est(iIter,:)';
                sse = sum((yTrue - yEst).^2);
                sst = sum((yTrue - mean(yTrue)).^2);
                template_R2_iter(iIter) = 1 - sse / sst;
            end
            [template_R2_mean(iFile), template_R2_lb(iFile), template_R2_ub(iFile)] = ...
                getCI(template_R2_iter, 1, 1);

            % ---------- marginalized templates: ORI and SF ----------
            template_true_2D = reshape(template_true, [nORI, nSF]);

            % true marginals
            margORI_true = mean(template_true_2D, 2)';   % 1 x nORI
            margSF_true  = mean(template_true_2D, 1);    % 1 x nSF

            margORI_true_plot(iFile, :) = margORI_true/max(margORI_true);
            margSF_true_plot(iFile, :)  = margSF_true/max(margSF_true);

            % estimated marginals per iteration
            margORI_est_iter = nan(nIter, nORI);
            margSF_est_iter  = nan(nIter, nSF);

            for iIter = 1:nIter
                template_est_2D = reshape(template_est(iIter,:), [nORI, nSF]);

                a = mean(template_est_2D, 2)';
                margORI_est_iter(iIter, :) = a/max(a);
                b = mean(template_est_2D, 1);
                margSF_est_iter(iIter, :) = b/max(b);
            end

            % summarize across iterations
            for iCh = 1:nORI
                [margORI_est_plot(iFile, iCh), ...
                    margORI_lb_plot(iFile, iCh), ...
                    margORI_ub_plot(iFile, iCh)] = ...
                    getCI(margORI_est_iter(:, iCh), 1, 1);
            end

            for iCh = 1:nSF
                [margSF_est_plot(iFile, iCh), ...
                    margSF_lb_plot(iFile, iCh), ...
                    margSF_ub_plot(iFile, iCh)] = ...
                    getCI(margSF_est_iter(:, iCh), 1, 1);
            end

            % ---------- metrics ----------
            pred_metrics_allIter = data_fitNOM.pred_metrics_allIter;

            for iM = 1:numel(metricNames)
                m = metricNames{iM};

                rmse_iter = nan(nIter,1);
                R2_iter = nan(nIter,1);
                dataMean_iter = nan(nIter,1);
                predMean_iter = nan(nIter,1);

                for iIter = 1:nIter
                    yData = pred_metrics_allIter{iIter}.metrics.([m '_data_allBins'])(:);
                    yPred = pred_metrics_allIter{iIter}.metrics.([m '_pred_allBins'])(:);

                    good = isfinite(yData) & isfinite(yPred);
                    yData = yData(good);
                    yPred = yPred(good);

                    rmse_iter(iIter) = sqrt(mean((yData - yPred).^2));
                    sse = sum((yData - yPred).^2);
                    sst = sum((yData - mean(yData)).^2);
                    R2_iter(iIter) = 1 - sse / sst;

                    dataMean_iter(iIter) = mean(yData);
                    predMean_iter(iIter) = mean(yPred);
                end

                [M.(m).rmse_mean(iFile), M.(m).rmse_lb(iFile), M.(m).rmse_ub(iFile)] = ...
                    getCI(rmse_iter, 1, 1);
                [M.(m).R2_mean(iFile), M.(m).R2_lb(iFile), M.(m).R2_ub(iFile)] = ...
                    getCI(R2_iter, 1, 1);

                M.(m).data_plot(iFile) = mean(dataMean_iter);
                [M.(m).pred_plot(iFile), M.(m).pred_lb_plot(iFile), M.(m).pred_ub_plot(iFile)] = ...
                    getCI(predMean_iter, 1, 1);
            end % iM

            % ---------- parameters ----------
            if iscell(data_fitNOM.params_est_allIter)
                P = cell2mat(cellfun(@(x) x(:)', data_fitNOM.params_est_allIter, 'UniformOutput', false));
            else
                P = data_fitNOM.params_est_allIter;
            end

            % Save full iteration-level parameter estimates for this file
            P_allFiles{iFile} = P;

            IN_est_iter = P(:,1:nParams);
            IN_rmse_iter = abs(IN_est_iter); % true IN = 0
            for iParam = 1:nParams
                [IN_rmse_mean(iFile,iParam), IN_rmse_lb(iFile,iParam), IN_rmse_ub(iFile,iParam)] = ...
                    getCI(IN_rmse_iter(:,iParam), 1, 1);
                [IN_est_plot(iFile,iParam), IN_est_lb_plot(iFile,iParam), IN_est_ub_plot(iFile,iParam)] = ...
                    getCI(IN_est_iter(:,iParam), 1, 1);
            end

            % ---------- Criterion----------
            criterion_est_iter = P(:,end);

            criterion_true_plot(iFile) = truth.criterion_DV_true;
            criterion_rmse_iter = abs(criterion_est_iter - criterion_true_plot(iFile));

            [criterion_rmse_mean(iFile), criterion_rmse_lb(iFile), criterion_rmse_ub(iFile)] = ...
                getCI(criterion_rmse_iter, 1, 1);
            [criterion_est_plot(iFile), criterion_est_lb_plot(iFile), criterion_est_ub_plot(iFile)] = ...
                getCI(criterion_est_iter, 1, 1);

            keepFile(iFile) = true;

        end % iFile

        % ---------------- keep valid rows only ----------------
        idx = keepFile;

        signal_plot = signal_plot(idx);
        Cz_plot = Cz_plot(idx);
        lambda_plot = lambda_plot(idx);

        % Template
        template_rmse_mean = template_rmse_mean(idx); template_rmse_lb = template_rmse_lb(idx); template_rmse_ub = template_rmse_ub(idx);
        template_R2_mean = template_R2_mean(idx); template_R2_lb = template_R2_lb(idx); template_R2_ub = template_R2_ub(idx);
        template_true_plot = template_true_plot(idx);
        template_est_plot = template_est_plot(idx); template_est_lb_plot = template_est_lb_plot(idx); template_est_ub_plot = template_est_ub_plot(idx);

        % Tuning functions
        margORI_true_plot = margORI_true_plot(idx, :);
        margORI_est_plot  = margORI_est_plot(idx, :);
        margORI_lb_plot   = margORI_lb_plot(idx, :);
        margORI_ub_plot   = margORI_ub_plot(idx, :);

        margSF_true_plot = margSF_true_plot(idx, :);
        margSF_est_plot  = margSF_est_plot(idx, :);
        margSF_lb_plot   = margSF_lb_plot(idx, :);
        margSF_ub_plot   = margSF_ub_plot(idx, :);

        % Metrics
        for iM = 1:numel(metricNames)
            m = metricNames{iM};
            M.(m).rmse_mean = M.(m).rmse_mean(idx);
            M.(m).rmse_lb = M.(m).rmse_lb(idx);
            M.(m).rmse_ub = M.(m).rmse_ub(idx);
            M.(m).R2_mean = M.(m).R2_mean(idx);
            M.(m).R2_lb = M.(m).R2_lb(idx);
            M.(m).R2_ub = M.(m).R2_ub(idx);

            M.(m).data_plot = M.(m).data_plot(idx);
            M.(m).pred_plot = M.(m).pred_plot(idx);
            M.(m).pred_lb_plot = M.(m).pred_lb_plot(idx);
            M.(m).pred_ub_plot = M.(m).pred_ub_plot(idx);
        end

        criterion_rmse_mean = criterion_rmse_mean(idx); criterion_rmse_lb = criterion_rmse_lb(idx); criterion_rmse_ub = criterion_rmse_ub(idx);
        criterion_true_plot = criterion_true_plot(idx);
        criterion_est_plot = criterion_est_plot(idx); criterion_est_lb_plot = criterion_est_lb_plot(idx); criterion_est_ub_plot = criterion_est_ub_plot(idx);

        IN_rmse_mean = IN_rmse_mean(idx,:); IN_rmse_lb = IN_rmse_lb(idx,:); IN_rmse_ub = IN_rmse_ub(idx,:);
        IN_true_plot = IN_true_plot(idx,:); IN_est_plot = IN_est_plot(idx,:); IN_est_lb_plot = IN_est_lb_plot(idx,:); IN_est_ub_plot = IN_est_ub_plot(idx,:);

        %% ---------------- PLOT 1: GoF summary figure ---------------
        fprintf('\n%s: 1. Plot the summary of GoF.\n', datetime('now'))

        h = figure('Position', [100 100 2e3 2e3]);
        hold on

        % Template RMSE
        subplot(4,3,1), hold on; box on;
        fxn_plotPerLambda(0, lambda_plot, Cz_plot, signal_plot, template_rmse_mean, template_rmse_lb, template_rmse_ub, Cz_unik, signal_unik, cmap_Cz);
        title('[RMSE] Template recovery'); ylabel('RMSE');

        % Template R2
        subplot(4,3,2), hold on;
        fxn_plotPerLambda(1, lambda_plot, Cz_plot, signal_plot, template_R2_mean, template_R2_lb, template_R2_ub, Cz_unik, signal_unik, cmap_Cz);
        title('[R^2] Template recovery'); ylabel('R^2');

        % pYES RMSE
        subplot(4,3,4), hold on; box on;
        fxn_plotPerLambda(0, lambda_plot, Cz_plot, signal_plot, M.pYES.rmse_mean, M.pYES.rmse_lb, M.pYES.rmse_ub, Cz_unik, signal_unik, cmap_Cz);
        title('[RMSE] pYES recovery'); ylabel('RMSE');

        % pYES R2
        subplot(4,3,5), hold on;
        fxn_plotPerLambda(1, lambda_plot, Cz_plot, signal_plot, M.pYES.R2_mean, M.pYES.R2_lb, M.pYES.R2_ub, Cz_unik, signal_unik, cmap_Cz);
        title('[R^2] pYES recovery'); ylabel('R^2');

        % Add legend
        fxn_addLegends(Cz_unik, cmap_Cz, signal_unik, lambda_unik, marker_list)

        % pC RMSE
        subplot(4,3,7), hold on; box on;
        fxn_plotPerLambda(0, lambda_plot, Cz_plot, signal_plot, M.pC.rmse_mean, M.pC.rmse_lb, M.pC.rmse_ub, Cz_unik, signal_unik, cmap_Cz);
        title('[RMSE] pC recovery'); ylabel('RMSE');

        % pC R2
        subplot(4,3,8), hold on;
        fxn_plotPerLambda(1, lambda_plot, Cz_plot, signal_plot, M.pC.R2_mean, M.pC.R2_lb, M.pC.R2_ub, Cz_unik, signal_unik, cmap_Cz);
        title('[R^2] pC recovery'); ylabel('R^2');

        % % pA RMSE
        subplot(4,3,10), hold on; box on;
        fxn_plotPerLambda(0, lambda_plot, Cz_plot, signal_plot, M.pA.rmse_mean, M.pA.rmse_lb, M.pA.rmse_ub, Cz_unik, signal_unik, cmap_Cz);
        title('[RMSE] pA recovery'); ylabel('RMSE');

        % pA R2
        subplot(4,3,11), hold on;
        fxn_plotPerLambda(1, lambda_plot, Cz_plot, signal_plot, M.pA.R2_mean, M.pA.R2_lb, M.pA.R2_ub, Cz_unik, signal_unik, cmap_Cz);
        title('[R^2] pA recovery'); ylabel('R^2');

        % Params RMSE
        switch iModelB_fit
            case 1, iPlot_params = [3,6,9 ];
            case 2, iPlot_params = [3,6   ]; % No Nshared
            case 3, iPlot_params = [   6,9]; % No Nmul
            case 4, iPlot_params = [3,   9]; % No Nadd
        end
        for iParam = 1:nParams
            subplot(4,3,iPlot_params(iParam)), hold on; box on;
            fxn_plotPerLambda(0, lambda_plot, Cz_plot, signal_plot, IN_rmse_mean(:,iParam), IN_rmse_lb(:,iParam), IN_rmse_ub(:,iParam), Cz_unik, signal_unik, cmap_Cz);
            title(sprintf('[RMSE] %s recovery', IN_param_names{iParam})); ylabel('RMSE'); xlabel('\lambda');
            ylim([0,1])
        end

        % Criterion RMSE
        subplot(4,3,12), hold on; box on;
        fxn_plotPerLambda(0, lambda_plot, Cz_plot, signal_plot, criterion_rmse_mean, criterion_rmse_lb, criterion_rmse_ub, Cz_unik, signal_unik, cmap_Cz);
        title('[RMSE] Criterion recovery'); ylabel('RMSE'); xlabel('\lambda');


        set(findall(gcf, '-property', 'fontsize'), 'fontsize', 18)
        sgtitle(sprintf('B%d %s (nIter=%d)', iModelB_fit, namesModelB{iModelB_fit}, nIter))
        saveas(h, fullfile(nameFolder_Figures_part, sprintf('Fig1_GoF_B%d.png', iModelB_fit)));
        close(h)

        %% ---------------- PLOT 2: True vs. estimated tuning functions ----------------
        fprintf('\n%s: 2. Plot tuning functions.\n', datetime('now'))
        % ORI
        fxn_tuningFxn( ...
            axis_tuning{1}, ...
            margORI_true_plot, ...
            margORI_est_plot, ...
            margORI_lb_plot, ...
            margORI_ub_plot, ...
            lambda_plot, Cz_plot, signal_plot, ...
            Cz_unik, signal_unik, cmap_Cz, ...
            'Orientation channel', 'Marginalized template weight', ...
            sprintf('Fig2 Marginalized ORI template | %s', namesModelB{iModelB_fit}), ...
            fullfile(nameFolder_Figures_part, sprintf('Fig2_margORI_B%d.png', iModelB_fit)));

        % SF
        fxn_tuningFxn( ...
            axis_tuning{2}, ...
            margSF_true_plot, ...
            margSF_est_plot, ...
            margSF_lb_plot, ...
            margSF_ub_plot, ...
            lambda_plot, Cz_plot, signal_plot, ...
            Cz_unik, signal_unik, cmap_Cz, ...
            'Spatial frequency channel', 'Marginalized template weight', ...
            sprintf('Fig2 Marginalized SF template | %s', namesModelB{iModelB_fit}), ...
            fullfile(nameFolder_Figures_part, sprintf('Fig2_margSF_B%d.png', iModelB_fit)));
        close all

        %% ---------------- PLOT 3: scatter summary figure ----------------
        fprintf('\n%s: 3. Plot the scatter summary of true vs. estimates.\n', datetime('now'))

        h = figure('Position', [100 100 2e3 2e3]);

        % Template
        subplot(2, 4, 1); hold on;
        fxn_plotScatter(template_true_plot, template_est_plot, template_est_lb_plot, template_est_ub_plot, ...
            lambda_plot, Cz_plot, signal_plot, Cz_unik, signal_unik, cmap_Cz);
        xlabel('True template mean');
        ylabel('Estimated template mean');
        title('Template recovery');

        % add legends
        % fxn_addLegends(Cz_unik, cmap_Cz, signal_unik, lambda_unik, marker_list)

        % pYES / pC / pA
        for iM = 1:numel(metricNames)
            m = metricNames{iM};
            subplot(2, 4, iM+1); hold on;
            fxn_plotScatter(M.(m).data_plot, M.(m).pred_plot, M.(m).pred_lb_plot, M.(m).pred_ub_plot, ...
                lambda_plot, Cz_plot, signal_plot, Cz_unik, signal_unik, cmap_Cz);
            xlabel(sprintf('Simulated %s', m));
            ylabel(sprintf('Predicted %s', m));
            title(sprintf('%s recovery', m));
        end

        % IN parameters
        switch iModelB_fit
            case 1, iPlot_params = [5,6,7];
            case 2, iPlot_params = [5,6  ]; % No Nshared
            case 3, iPlot_params = [   6,7]; % No Nmul
            case 4, iPlot_params = [5,  7]; % No Nadd
        end
        for iParam = 1:nParams
            subplot(2, 4, iPlot_params(iParam)); hold on;
            fxn_plotScatter(IN_true_plot(:,iParam), IN_est_plot(:,iParam), IN_est_lb_plot(:,iParam), IN_est_ub_plot(:,iParam), ...
                lambda_plot, Cz_plot, signal_plot, Cz_unik, signal_unik, cmap_Cz);
            xlabel(sprintf('True %s', IN_param_names{iParam}));
            ylabel(sprintf('Estimated %s', IN_param_names{iParam}));
            title(sprintf('%s recovery', IN_param_names{iParam}));
        end

        % criterion
        subplot(2, 4, 8); hold on;
        fxn_plotScatter(criterion_true_plot, criterion_est_plot, criterion_est_lb_plot, criterion_est_ub_plot, ...
            lambda_plot, Cz_plot, signal_plot, Cz_unik, signal_unik, cmap_Cz);
        xlabel('True criterion');
        ylabel('Estimated criterion');
        title('Criterion recovery');

        set(findall(gcf, '-property', 'fontsize'), 'fontsize', 16)
        sgtitle(sprintf('B%d %s (nIter=%d)', iModelB_fit, namesModelB{iModelB_fit}, nIter))

        saveas(h, fullfile(nameFolder_Figures_part, sprintf('Fig3_Scatter_B%d.png', iModelB_fit)));
        close(h)

        %% ---------------- PLOT 4: correlation among estimated parameters across iterations ----------------
        fprintf('\n%s: 4. Plot corr between estimated params across iterations.\n', datetime('now'))

        % P_allFiles{iFile}: [nIter x (nParams+1)]
        % columns = [estimated IN params, estimated criterion]

        param_names_all = [IN_param_names, {'criterion'}];
        nParams_all = numel(param_names_all);

        % -------- compile pooled iteration-level table --------
        params_est_iter_all = [];
        lambda_iter_all = [];
        signal_iter_all = [];
        Cz_iter_all = [];

        for iFile = 1:numel(P_allFiles)

            P = P_allFiles{iFile};
            if isempty(P)
                continue
            end

            nIter_this = size(P, 1);

            params_est_iter_all = [params_est_iter_all; P];
            lambda_iter_all = [lambda_iter_all; repmat(lambda_plot(iFile), nIter_this, 1)];
            signal_iter_all = [signal_iter_all; repmat(signal_plot(iFile), nIter_this, 1)];
            Cz_iter_all = [Cz_iter_all; repmat(Cz_plot(iFile), nIter_this, 1)];
        end % iFile

        % remove rows with NaN in either params or covariates
        goodRow = all(isfinite(params_est_iter_all), 2) & ...
            isfinite(lambda_iter_all) & isfinite(signal_iter_all) & isfinite(Cz_iter_all);

        params_est_iter_all = params_est_iter_all(goodRow, :);
        lambda_iter_all = lambda_iter_all(goodRow);
        signal_iter_all = signal_iter_all(goodRow);
        Cz_iter_all = Cz_iter_all(goodRow);

        % -------- panel layout --------
        nPairs = nchoosek(nParams_all, 2);
        nCols = 2;
        nRows = ceil(nPairs / nCols);

        h = figure('Position', [100 100 1200 350*nRows]);
        tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

        iPanel = 0;
        for iParam1 = 1:nParams_all-1
            for iParam2 = iParam1+1:nParams_all
                iPanel = iPanel + 1;
                nexttile; hold on;

                x_iter = params_est_iter_all(:, iParam1);
                y_iter = params_est_iter_all(:, iParam2);

                % partial correlation controlling for lambda, signal, Cz
                [r_partial, p_partial] = partialcorr( ...
                    x_iter, y_iter, [lambda_iter_all, signal_iter_all, Cz_iter_all]);

                fxn_plotCorrAcrossIter( ...
                    x_iter, y_iter, ...
                    lambda_iter_all, Cz_iter_all, signal_iter_all, ...
                    Cz_unik, signal_unik, cmap_Cz, marker_list, ...
                    param_names_all{iParam1}, param_names_all{iParam2}, ...
                    sprintf('%s vs %s\npartial r = %.2f, p = %.3g', ...
                    param_names_all{iParam1}, param_names_all{iParam2}, ...
                    r_partial, p_partial));
            end % iParam2
        end % iParam1

        fxn_addLegends(Cz_unik, cmap_Cz, signal_unik, unique(lambda_iter_all), marker_list)

        saveas(h, fullfile(nameFolder_Figures_part, sprintf('Fig4_paramCorr_B%d.png', iModelB_fit)));
        close all

        % clear variables to save space
        clear *plot *mean* *lb* *ub*
    end % iModelB_fit

    fprintf('\n\n%s: ALL DONE \n\n', datetime('now'))
end
%% local functions

function S = build_file_table(nameFolder_Data_OOD, nameFolder_Data_NOM_Trialwise, iModelA_fit)
d = dir(nameFolder_Data_OOD);
d = d([d.isdir]);
d = d(~ismember({d.name},{'.','..'}));

S = struct([]);
k = 0;

for i = 1:numel(d)
    nameIO = d(i).name;

    tok = regexp(nameIO, ...
        'IO_cN([-\d\.]+)_cG([-\d\.]+).*_Cz([-\d\.]+)_cont([-\d\.]+)_whiten([-\d\.]+)', ...
        'tokens', 'once');
    if isempty(tok), continue; end

    truthFile = fullfile(nameFolder_Data_OOD, nameIO, 'truth.mat');
    if ~exist(truthFile,'file'), continue; end

    nomFolder = fullfile(nameFolder_Data_NOM_Trialwise, nameIO);
    compPat = sprintf('n*_A%d_compIV.mat', iModelA_fit);
    compDir = dir(fullfile(nomFolder, compPat));
    compDir = compDir(~contains({compDir.name}, 'min'));
    if isempty(compDir), continue; end

    k = k + 1;
    S(k).nameIO = nameIO;
    S(k).signalCST = str2double(tok{2}) / 100;
    S(k).Cz = str2double(tok{3});
    S(k).lambda = str2double(tok{5});
    S(k).truthFile = truthFile;
    S(k).compIVFile = fullfile(compDir(1).folder, compDir(1).name);
    S(k).nomFolder = nomFolder;
end
end
%% 1
function fxn_plotPerLambda(flag_plotR2, lambda_all, Cz_all, signal_all, y_mean, y_lb, y_ub, Cz_unik, signal_unik, cmap_Cz)
signal_min = min(signal_unik);
signal_max = max(signal_unik);

pairMat = [Cz_all(:), signal_all(:)];
[~, ~, pairIdx] = unique(pairMat, 'rows');

for iPair = 1:max(pairIdx)
    idx_pair = pairIdx == iPair;
    [x_sort, ord] = sort(lambda_all(idx_pair));
    y_sort = y_mean(idx_pair);
    y_sort = y_sort(ord);
    plot(x_sort, y_sort, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
end

for iPt = 1:numel(y_mean)
    iCz = find(Cz_unik == Cz_all(iPt), 1, 'first');
    color_face = cmap_Cz(iCz, :);

    if signal_max == signal_min
        gray_edge = 0.4;
    else
        frac = (signal_all(iPt) - signal_min) / (signal_max - signal_min);
        gray_edge = 0.85 - 0.65 * frac;
    end
    color_edge = gray_edge * [1 1 1];

    line([lambda_all(iPt), lambda_all(iPt)], [y_lb(iPt), y_ub(iPt)], 'Color', color_edge, 'LineWidth', 1.3);
    scatter(lambda_all(iPt), y_mean(iPt), 85, 'MarkerFaceColor', color_face, 'MarkerEdgeColor', color_edge, 'LineWidth', 1.4);
end

if flag_plotR2
    ylim([.5, 1])
else
    ylim([0, .1])
    if y_ub>.1
        ylim([0,1])
    end
end
box on; set(gca, 'FontSize', 11); xlabel('\lambda');
end

%% 2
function fxn_tuningFxn(xAxis, ...
    marg_true, marg_est, marg_lb, marg_ub, ...
    lambda_all, Cz_all, signal_all, ...
    Cz_unik, signal_unik, cmap_Cz, ...
    xlab, ylab, ttl, saveName)

marker_list = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h', 'x', '+'};

lambda_unik = unique(lambda_all);
nLambda = numel(lambda_unik);

% choose panel layout
nCols = ceil(sqrt(nLambda));
nRows = ceil(nLambda / nCols);

figure('Position', [100 100 420*nCols 320*nRows]);

tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

signal_min = min(signal_unik);
signal_max = max(signal_unik);
nCh = numel(xAxis);

for iLam = 1:nLambda
    nexttile; hold on;

    if iLam==2
        % add legends once, after tiledlayout is created
        fxn_addLegends(Cz_unik, cmap_Cz, signal_unik, lambda_unik, marker_list);

    end
    idx_lam = lambda_all == lambda_unik(iLam);
    idx_cond = find(idx_lam);

    for iiCond = 1:numel(idx_cond)
        iCond = idx_cond(iiCond);

        iCz = find(Cz_unik == Cz_all(iCond), 1, 'first');
        color_face = cmap_Cz(iCz, :);

        if signal_max == signal_min
            gray_edge = 0.4;
        else
            frac = (signal_all(iCond) - signal_min) / (signal_max - signal_min);
            gray_edge = 0.85 - 0.65 * frac;
        end
        color_edge = gray_edge * [1 1 1];

        % true marginalized template
        plot(xAxis, marg_true(iCond,1:nCh), '-r', 'LineWidth', 2.5);

        % estimated marginalized template
        plot(xAxis, marg_est(iCond,1:nCh), '--', 'Color', color_face, 'LineWidth', 0.8);

        scatter(xAxis, marg_est(iCond,:), 35, ...
            'MarkerFaceColor', color_face, ...
            'MarkerEdgeColor', color_edge, ...
            'LineWidth', 0.9);
    end % ii
    ylim([-.4, 1.2])
    box on;
    set(gca, 'FontSize', 11);
    xlabel(xlab);
    ylabel(ylab);
    title(sprintf('\\lambda = %.1f', lambda_unik(iLam)));
end % iLam

sgtitle(ttl);
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 12)
saveas(gcf, saveName);
end

%% 3
function fxn_plotScatter(x_true, y_est, y_lb, y_ub, lambda_all, Cz_all, signal_all, Cz_unik, signal_unik, cmap_Cz)

signal_min = min(signal_unik);
signal_max = max(signal_unik);

lambda_unik = unique(lambda_all);

% connect points with same lambda
for iLam = 1:numel(lambda_unik)
    idx_lam = lambda_all == lambda_unik(iLam);
    [x_sort, ord] = sort(x_true(idx_lam));
    y_sort = y_est(idx_lam);
    y_sort = y_sort(ord);
    plot(x_sort, y_sort, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
end

for iPt = 1:numel(y_est)
    iCz = find(Cz_unik == Cz_all(iPt), 1, 'first');
    color_face = cmap_Cz(iCz, :);

    if signal_max == signal_min
        gray_edge = 0.4;
    else
        frac = (signal_all(iPt) - signal_min) / (signal_max - signal_min);
        gray_edge = 0.85 - 0.65 * frac;
    end
    color_edge = gray_edge * [1 1 1];

    % Errorbar
    line([x_true(iPt), x_true(iPt)], [y_lb(iPt), y_ub(iPt)], 'Color', color_edge, 'LineWidth', 1.1);

    % Data point
    scatter(x_true(iPt), y_est(iPt), 85, ...
        'MarkerFaceColor', color_face, ...
        'MarkerEdgeColor', color_edge, ...
        'LineWidth', 1.4);
end

xy_all = [x_true(:); y_est(:)];
xy_all = xy_all(isfinite(xy_all));
xy_min = min(xy_all);
xy_max = max(xy_all);

if xy_min == xy_max
    xy_min = xy_min - 0.1;
    xy_max = xy_max + 0.1;
end

% Unity line
plot([xy_min xy_max], [xy_min xy_max], 'k--', 'LineWidth', 1.2);

axis square;
xlim([xy_min xy_max]);
ylim([xy_min xy_max]);
end

%% 4
function fxn_plotCorrAcrossIter(x_all, y_all, ...
    lambda_all, Cz_all, signal_all, ...
    Cz_unik, signal_unik, cmap_Cz, marker_list, ...
    xlab, ylab, ttl)

signal_min = min(signal_unik);
signal_max = max(signal_unik);

lambda_unik = unique(lambda_all);
if numel(lambda_unik) > numel(marker_list)
    error('Not enough marker styles for the number of unique lambda values.');
end

[~, iLam_all] = ismember(lambda_all, lambda_unik);

for iPt = 1:numel(x_all)
    iCz = find(Cz_unik == Cz_all(iPt), 1, 'first');
    color_face = cmap_Cz(iCz, :);

    if signal_max == signal_min
        gray_edge = 0.4;
    else
        frac = (signal_all(iPt) - signal_min) / (signal_max - signal_min);
        gray_edge = 0.85 - 0.65 * frac;
    end
    color_edge = gray_edge * [1 1 1];

    this_marker = marker_list{iLam_all(iPt)};

    scatter(x_all(iPt), y_all(iPt), 18, ...
        'Marker', this_marker, ...
        'MarkerFaceColor', color_face, ...
        'MarkerEdgeColor', color_edge, ...
        'LineWidth', 0.6, ...
        'MarkerFaceAlpha', 0.35, ...
        'MarkerEdgeAlpha', 0.35);
end

xlim([min(x_all(:)), max(x_all(:))]);
ylim([min(y_all(:)), max(y_all(:))]);

box on;
set(gca, 'FontSize', 12);
xlabel(sprintf('Estimated %s', xlab));
ylabel(sprintf('Estimated %s', ylab));
title(ttl);
end

%% Legends
function fxn_addLegends(Cz_unik, cmap_Cz, signal_unik, lambda_unik, marker_list)
% ---------------- Cz legend: face color ----------------
hCz = gobjects(numel(Cz_unik), 1);
for iCz = 1:numel(Cz_unik)
    hCz(iCz) = scatter(nan, nan, 60, ...
        'MarkerFaceColor', cmap_Cz(iCz,:), ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1);
end

leg1 = legend(hCz, ...
    arrayfun(@(x) sprintf('Cz = %.1f', x), Cz_unik, 'UniformOutput', false), ...
    'Location', 'southwest');
title(leg1, 'True criterion');

% ---------------- signal legend: edge grayness ----------------
signal_min = min(signal_unik);
signal_max = max(signal_unik);

ax1 = gca;
ax2 = axes('Position', ax1.Position, 'Color', 'none', 'Visible', 'off');
hold(ax2, 'on');

hSig = gobjects(numel(signal_unik), 1);
for iSig = 1:numel(signal_unik)
    if signal_max == signal_min
        gray_edge = 0.4;
    else
        frac = (signal_unik(iSig) - signal_min) / (signal_max - signal_min);
        gray_edge = 0.85 - 0.65 * frac;
    end
    color_edge = gray_edge * [1 1 1];

    hSig(iSig) = scatter(ax2, nan, nan, 60, ...
        'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', color_edge, ...
        'LineWidth', 1.5);
end

leg2 = legend(ax2, hSig, ...
    arrayfun(@(x) sprintf('signal = %.2f', x), signal_unik, 'UniformOutput', false), ...
    'Location', 'southeast');
title(leg2, 'Signal contrast');

% ---------------- lambda legend: marker shape ----------------
if nargin >= 5 && ~isempty(lambda_unik)
    if numel(lambda_unik) > numel(marker_list)
        error('Not enough marker styles for the number of unique lambda values.');
    end

    ax3 = axes('Position', ax1.Position, 'Color', 'none', 'Visible', 'off');
    hold(ax3, 'on');

    hLam = gobjects(numel(lambda_unik), 1);
    for iLam = 1:numel(lambda_unik)
        hLam(iLam) = scatter(ax3, nan, nan, 60, ...
            'Marker', marker_list{iLam}, ...
            'MarkerFaceColor', [1 1 1], ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.2);
    end

    leg3 = legend(ax3, hLam, ...
        arrayfun(@(x) sprintf('\\lambda = %.1f', x), lambda_unik, 'UniformOutput', false), ...
        'Location', 'northeast');
    title(leg3, 'Whitening');
end

end