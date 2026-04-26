% simPlot3_INrecovery.m
%
% Global heatmaps for recovery metrics across simulation conditions
%
% Organizes conditions as:
% x-axis : Nshared_true
% y-axis : Nmul_true
% columns : Nadd_true
% rows : gaborCST
%
% Values are averaged over Cz_true within each cell.
%
% This script follows the same loading / compile logic as simPlot2,
% but summarizes all conditions into heatmaps.
%
% Created by Shutian Xue

clear; clc; close all;
set(0, 'DefaultFigureVisible', 'off');

%--------------%
SX_RC1_setting;
%--------------%

%% settings
str_part = 'Lambda1'; % <-- change if needed
iModelA_fit = 1; %1=use data-derived template; 2=use true template
mask_pC = [.6, .8]; % only analyze simulated datasets wtih pC falling within this range
nBfit = 4; % number of fitted models

nameFolder_Data = sprintf('%s/Data_%s', nameFolder_server, str_part);
nameFolder_Data_OOD = sprintf('%s/Data_OOD_%d%d', nameFolder_Data, nORI, nSF);
nameFolder_Data_NOM_Trialwise = sprintf('%s/Data_NOM_Trialwise_%d%d', nameFolder_Data, nORI, nSF);

nameFolder_Figures_part3 = fullfile(nameFolder_Figures, sprintf('IO_%s_A%d/Part3', str_part, iModelA_fit));
if ~exist(nameFolder_Figures_part3, 'dir')
    mkdir(nameFolder_Figures_part3);
end

namesMetrics_behav = {'pYES','pC','pA'};

%% find IO folders
nameDir = dir(nameFolder_Data_NOM_Trialwise);
nameDir = nameDir([nameDir.isdir]);
nameDir = nameDir(~ismember({nameDir.name}, {'.', '..'}));
nFiles = numel(nameDir);

if nFiles == 0
    error('No IO folders found in %s', nameFolder_Data_NOM_Trialwise);
end

fprintf('\nFound %d IO folders.\n', nFiles);

%% Parse all conditions once
info_all = cell(nFiles, 1);
keepParse = false(nFiles, 1);

for iFile = 1:nFiles
    nameIO = nameDir(iFile).name;
    info = fxn_parse_nameIO(nameIO);

    if isempty(fieldnames(info))
        fprintf('Could not parse: %s\n', nameIO);
        continue
    end

    info.nameIO = nameIO;
    info.folder_OOD = fullfile(nameFolder_Data_OOD, nameIO);
    info.folder_NOM = fullfile(nameFolder_Data_NOM_Trialwise, nameIO);

    info_all{iFile} = info;
    keepParse(iFile) = true;
end

info_all = info_all(keepParse);
nFiles = numel(info_all);

fprintf('Parsed %d valid IO folders.\n', nFiles);

%% Compile all fitted models

R = struct([]);

fprintf('\n============================================\n');
fprintf('%s\n', datetime('now'))
fprintf('Compiling all fitted models\n');
fprintf('============================================\n');

nFilesKeep = 0; % To count the number of file with pC within the range
for iFile = 1:nFiles
    fprintf('\n%d/%d', iFile, nFiles)
    info = info_all{iFile};

    % =========================================================
    % 1. Load data once
    % =========================================================
    nameFile_truth = fullfile(info.folder_OOD, 'truth.mat');
    if ~exist(nameFile_truth, 'file')
        fprintf('Missing truth.mat: %s\n', nameFile_truth);
        continue
    end
    truth = load(nameFile_truth);

    str_loadCompIV = sprintf('n*_A%d_compIV.mat', iModelA_fit);
    nameDir_compIV = dir(fullfile(info.folder_NOM, str_loadCompIV));
    nameDir_compIV = nameDir_compIV(~contains({nameDir_compIV.name}, 'min'));
    if isempty(nameDir_compIV)
        fprintf('No compIV file found in %s\n', info.folder_NOM);
        continue
    end
    data_compIV = load(fullfile(nameDir_compIV(1).folder, nameDir_compIV(1).name));

    % =========================================================
    % 2. metadata and iteration count
    % =========================================================
    fileBase = info;

    if isfield(truth, 'nIter')
        nIter = truth.nIter;
    else
        template_est_tmp = data_compIV.template_tmpl_allIter;
        % if ndims(template_est_tmp) == 3
        nIter = size(template_est_tmp, 1);
        % else
        % nIter = size(template_est_tmp, 1);
        % end
    end

    % initialize shared fields
    % true values -> template-related values -> tuning-related values ->
    % metric-related values -> parameter-related values -> criterion-related values

    % ---------- true values ----------
    fileBase.lambda_whiten = nan;
    fileBase.C_contribution = nan;

    % ---------- template-related values ----------
    fileBase.template_rmse = nan;
    fileBase.template_R2 = nan;
    fileBase.template_rmse_iter = [];

    fileBase.nBasisORI_mode = nan;
    fileBase.nBasisORI_mode_pct = nan;
    fileBase.nBasisSF_mode = nan;
    fileBase.nBasisSF_mode_pct = nan;
    fileBase.basisFxnORI_mode = "";
    fileBase.basisFxnORI_mode_pct = nan;
    fileBase.basisFxnSF_mode = "";
    fileBase.basisFxnSF_mode_pct = nan;
    fileBase.ridge_mode = "";
    fileBase.ridge_mode_pct = nan;

    % ---------- tuning-related values ----------
    fileBase.margORI_true = nan(1, nORI);
    fileBase.margSF_true = nan(1, nSF);
    fileBase.margORI_est_med = nan(1, nORI);
    fileBase.margSF_est_med = nan(1, nSF);
    fileBase.margORI_est_iter_norm = [];
    fileBase.margSF_est_iter_norm = [];

    % ---------- metric-related values ----------
    fileBase.pYES_data_med = nan; fileBase.pYES_data_lb = nan; fileBase.pYES_data_ub = nan;
    fileBase.pC_data_med = nan; fileBase.pC_data_lb = nan; fileBase.pC_data_ub = nan;
    fileBase.pA_data_med = nan; fileBase.pA_data_lb = nan; fileBase.pA_data_ub = nan;

    fileBase.DV_allBins_med = [];
    fileBase.nTrials_allBins_med = [];
    fileBase.pYES_data_curve_med = [];
    fileBase.pC_data_curve_med = [];
    fileBase.pA_data_curve_med = [];

    % ---------- criterion-related values ----------
    fileBase.criterion_DV_true = nan;
    fileBase.keep_pC = false;

    fileBase.lambda_whiten = truth.lambda_whiten;
    fileBase.C_contribution = truth.C_contribution;

    % =========================================================
    % 1. template recovery
    % =========================================================
    [template_rmse_iter, template_R2_iter] = fxn_templateRecovery_simPlot2Style_iter( ...
        truth, data_compIV, nIter);

    [fileBase.template_rmse, ~, ~] = getCI(template_rmse_iter, 1, 1);
    [fileBase.template_R2, ~, ~] = getCI(template_R2_iter, 1, 1);
    fileBase.template_rmse_iter = template_rmse_iter(:)';

    % ---------- marginalized templates ----------
    template_true_vec = truth.template_true(:)';
    template_true_2D = reshape(template_true_vec, [nORI, nSF]);

    margORI_true = mean(template_true_2D, 2)';
    margSF_true = mean(template_true_2D, 1);

    % Make peak = 1
    margORI_true = margORI_true ./ max(margORI_true);
    margSF_true = margSF_true ./ max(margSF_true);

    fileBase.margORI_true = margORI_true;
    fileBase.margSF_true = margSF_true;

    template_est = data_compIV.template_tmpl_allIter;
    if ndims(template_est) == 3
        template_est = reshape(template_est, size(template_est,1), []);
    end

    nIter_est = size(template_est, 1);
    margORI_est_iter = nan(nIter_est, nORI);
    margSF_est_iter = nan(nIter_est, nSF);

    for iIter = 1:nIter_est
        tmp2D = reshape(template_est(iIter,:), [nORI, nSF]);
        margORI_est_iter(iIter, :) = mean(tmp2D, 2)';
        margSF_est_iter(iIter, :) = mean(tmp2D, 1);
    end

    % Normalize each iteration separately so iteration CI reflects tuning-shape variability.
    margORI_est_iter_norm = margORI_est_iter;
    margSF_est_iter_norm = margSF_est_iter;
    for iIter = 1:nIter_est
        mxORI = max(margORI_est_iter_norm(iIter, :), [], 2);
        mxSF = max(margSF_est_iter_norm(iIter, :), [], 2);
        if isfinite(mxORI) && mxORI ~= 0
            margORI_est_iter_norm(iIter, :) = margORI_est_iter_norm(iIter, :) ./ mxORI;
        end
        if isfinite(mxSF) && mxSF ~= 0
            margSF_est_iter_norm(iIter, :) = margSF_est_iter_norm(iIter, :) ./ mxSF;
        end
    end

    % Take median across iterations
    margORI_est_med = getCI(margORI_est_iter_norm, 1, 1);
    margSF_est_med = getCI(margSF_est_iter_norm, 1, 1);

    % Store
    fileBase.margORI_est_med = margORI_est_med;
    fileBase.margSF_est_med = margSF_est_med;
    fileBase.margORI_est_iter_norm = margORI_est_iter_norm;
    fileBase.margSF_est_iter_norm = margSF_est_iter_norm;

    % =========================================================
    % 2. Basis-settings-related values
    % =========================================================
    nBasisORI_allIter = data_compIV.nBasisORI_tmpl_allIter;
    nBasisSF_allIter = data_compIV.nBasisSF_tmpl_allIter;

    fileBase.nBasisORI_mode = mode(nBasisORI_allIter);
    fileBase.nBasisORI_mode_pct = mean(nBasisORI_allIter == mode(nBasisORI_allIter));

    fileBase.nBasisSF_mode = mode(nBasisSF_allIter);
    fileBase.nBasisSF_mode_pct = mean(nBasisSF_allIter == mode(nBasisSF_allIter));

    S = string(data_compIV.basisFxnORI_tmpl_allIter(:));
    [u, ~, ic] = unique(S);
    counts = accumarray(ic, 1);
    [~, idxMax] = max(counts);
    fileBase.basisFxnORI_mode = u(idxMax);
    fileBase.basisFxnORI_mode_pct = counts(idxMax) / numel(S);

    S = string(data_compIV.basisFxnSF_tmpl_allIter(:));
    [u, ~, ic] = unique(S);
    counts = accumarray(ic, 1);
    [~, idxMax] = max(counts);
    fileBase.basisFxnSF_mode = u(idxMax);
    fileBase.basisFxnSF_mode_pct = counts(idxMax) / numel(S);

    % L2 ridge summary
    S = string(data_compIV.ridge_tmpl_allIter(:));
    [u, ~, ic] = unique(S);
    counts = accumarray(ic, 1);
    [~, idxMax] = max(counts);
    fileBase.ridge_mode = u(idxMax);
    fileBase.ridge_mode_pct = counts(idxMax) / numel(S);

    % =========================================================
    % 3. Metric-related value
    % Use first available Bfit file only because simulated data do not depend on fitted model
    % =========================================================
    found_fit_for_data = false;

    for iModelB_fit_probe = 1:nBfit
        str_loadfitNOM = sprintf('n*_A%dB%d.mat', iModelA_fit, iModelB_fit_probe);
        nameDir_fitNOM_probe = dir(fullfile(info.folder_NOM, str_loadfitNOM));
        nameDir_fitNOM_probe = nameDir_fitNOM_probe(~contains({nameDir_fitNOM_probe.name}, 'min'));

        if isempty(nameDir_fitNOM_probe)
            continue
        end

        data_fitNOM_probe = load(fullfile(nameDir_fitNOM_probe(1).folder, nameDir_fitNOM_probe(1).name));
        pred_metrics_allIter_probe = data_fitNOM_probe.pred_metrics_allIter;
        nIter_probe = numel(pred_metrics_allIter_probe);

        pYES_data_iter = nan(nIter_probe,1);
        pC_data_iter = nan(nIter_probe,1);
        pA_data_iter = nan(nIter_probe,1);

        nBins_curve = numel(pred_metrics_allIter_probe{1}.metrics.IV_allBins);
        DV_allBins_iter = nan(nIter_probe, nBins_curve);
        nTrials_allBins_iter = nan(nIter_probe, nBins_curve);
        pYES_data_curve_iter = nan(nIter_probe, nBins_curve);
        pC_data_curve_iter = nan(nIter_probe, nBins_curve);
        pA_data_curve_iter = nan(nIter_probe, nBins_curve);

        for iIter = 1:nIter_probe
            pYES_data_iter(iIter) = mean(pred_metrics_allIter_probe{iIter}.metrics.pYES_data_allBins, 'omitnan');
            pC_data_iter(iIter) = mean(pred_metrics_allIter_probe{iIter}.metrics.pC_data_allBins, 'omitnan');
            pA_data_iter(iIter) = mean(pred_metrics_allIter_probe{iIter}.metrics.pA_data_allBins, 'omitnan');

            DV_allBins_iter(iIter, :) = pred_metrics_allIter_probe{iIter}.metrics.IV_allBins;
            nTrials_allBins_iter(iIter, :) = pred_metrics_allIter_probe{iIter}.metrics.nTrials_allBins;
            pYES_data_curve_iter(iIter, :) = pred_metrics_allIter_probe{iIter}.metrics.pYES_data_allBins;
            pC_data_curve_iter(iIter, :) = pred_metrics_allIter_probe{iIter}.metrics.pC_data_allBins;
            pA_data_curve_iter(iIter, :) = pred_metrics_allIter_probe{iIter}.metrics.pA_data_allBins;
        end

        [fileBase.pYES_data_med, fileBase.pYES_data_lb, fileBase.pYES_data_ub] = getCI(pYES_data_iter, 1, 1);
        [fileBase.pC_data_med, fileBase.pC_data_lb, fileBase.pC_data_ub] = getCI(pC_data_iter, 1, 1);
        [fileBase.pA_data_med, fileBase.pA_data_lb, fileBase.pA_data_ub] = getCI(pA_data_iter, 1, 1);

        fileBase.DV_allBins_med = nan(1, nBins_curve);
        fileBase.nTrials_allBins_med = nan(1, nBins_curve);
        fileBase.pYES_data_curve_med = nan(1, nBins_curve);
        fileBase.pC_data_curve_med = nan(1, nBins_curve);
        fileBase.pA_data_curve_med = nan(1, nBins_curve);

        for iBin = 1:nBins_curve
            [fileBase.DV_allBins_med(iBin), ~, ~] = getCI(DV_allBins_iter(:, iBin), 1, 1);
            [fileBase.nTrials_allBins_med(iBin), ~, ~] = getCI(nTrials_allBins_iter(:, iBin), 1, 1);
            [fileBase.pYES_data_curve_med(iBin), ~, ~] = getCI(pYES_data_curve_iter(:, iBin), 1, 1);
            [fileBase.pC_data_curve_med(iBin), ~, ~] = getCI(pC_data_curve_iter(:, iBin), 1, 1);
            [fileBase.pA_data_curve_med(iBin), ~, ~] = getCI(pA_data_curve_iter(:, iBin), 1, 1);
        end

        found_fit_for_data = true;
        break
    end

    if ~found_fit_for_data
        fprintf('No fitNOM file found in %s for any Bfit\n', info.folder_NOM);
        continue
    end

    % ---------- criterion true ----------
    if isfield(truth, 'criterion_DV_true')
        fileBase.criterion_DV_true = truth.criterion_DV_true;
    end

    % ---------- filtering rule: simulated pC median only ----------
    fileBase.keep_pC = (fileBase.pC_data_med >= mask_pC(1)) && (fileBase.pC_data_med <= mask_pC(2));
    if fileBase.keep_pC, nFilesKeep = nFilesKeep+1; end
    fprintf(': med sim pC=%.2f, keep? %d', fileBase.pC_data_med, fileBase.keep_pC)

    % =========================================================
    % 4. Fitted-model-related values
    % =========================================================
    for iModelB_fit = 1:nBfit

        str_loadfitNOM = sprintf('n*_A%dB%d.mat', iModelA_fit, iModelB_fit);
        nameDir_fitNOM = dir(fullfile(info.folder_NOM, str_loadfitNOM));
        nameDir_fitNOM = nameDir_fitNOM(~contains({nameDir_fitNOM.name}, 'min'));
        if isempty(nameDir_fitNOM)
            fprintf('No fitNOM file found in %s for Bfit=%d\n', info.folder_NOM, iModelB_fit);
            continue
        end
        data_fitNOM = load(fullfile(nameDir_fitNOM(1).folder, nameDir_fitNOM(1).name));

        Ri = fileBase;
        Ri.iModelB_fit = iModelB_fit;

        % initialize fit-level fields
        % true values -> template-related values -> tuning-related values ->
        % metric-related values -> parameter-related values -> criterion-related values

        % ---------- true values ----------

        % ---------- template-related values ----------
        % inherited from fileBase

        % ---------- tuning-related values ----------
        % inherited from fileBase

        % ---------- metric-related values ----------
        Ri.pYES_rmse = nan; Ri.pYES_R2 = nan;
        Ri.pC_rmse = nan; Ri.pC_R2 = nan;
        Ri.pA_rmse = nan; Ri.pA_R2 = nan;

        Ri.pYES_pred_med = nan; Ri.pYES_pred_lb = nan; Ri.pYES_pred_ub = nan;
        Ri.pC_pred_med = nan; Ri.pC_pred_lb = nan; Ri.pC_pred_ub = nan;
        Ri.pA_pred_med = nan; Ri.pA_pred_lb = nan; Ri.pA_pred_ub = nan;

        Ri.pYES_pred_curve_med = [];
        Ri.pC_pred_curve_med = [];
        Ri.pA_pred_curve_med = [];

        Ri.nLL_med = nan; Ri.nLL_lb = nan; Ri.nLL_ub = nan;

    % ---------- parameter-related values ----------
    Ri.Nmul_rmse = nan; Ri.Nmul_est_med = nan; Ri.Nmul_est_lb = nan; Ri.Nmul_est_ub = nan;
    Ri.Nadd_rmse = nan; Ri.Nadd_est_med = nan; Ri.Nadd_est_lb = nan; Ri.Nadd_est_ub = nan;
    Ri.Nshared_rmse = nan; Ri.Nshared_est_med = nan; Ri.Nshared_est_lb = nan; Ri.Nshared_est_ub = nan;

    % ---------- criterion-related values ----------
        Ri.criterion_DV_rmse = nan;
        Ri.criterion_DV_est_med = nan;
        Ri.criterion_DV_est_lb = nan;
        Ri.criterion_DV_est_ub = nan;

        % ---------- fit score ----------
        if isfield(data_fitNOM, 'nLL_test_allIter')
            [Ri.nLL_med, Ri.nLL_lb, Ri.nLL_ub] = getCI(data_fitNOM.nLL_test_allIter(:), 1, 1);
        end

        pred_metrics_allIter = data_fitNOM.pred_metrics_allIter;
        nIter_fit = numel(pred_metrics_allIter);

        % ---------- metric recovery ----------
        for iMetric = 1:numel(namesMetrics_behav)
            m = namesMetrics_behav{iMetric};
            [rmse_iter, R2_iter] = fxn_metricRecovery_iter(pred_metrics_allIter, m);

            [Ri.(sprintf('%s_rmse', m)), ~, ~] = getCI(rmse_iter, 1, 1);
            [Ri.(sprintf('%s_R2', m)), ~, ~] = getCI(R2_iter, 1, 1);
        end

        % ---------- predicted scalar metrics ----------
        pYES_pred_iter = nan(nIter_fit,1);
        pC_pred_iter = nan(nIter_fit,1);
        pA_pred_iter = nan(nIter_fit,1);

        nBins_curve = numel(pred_metrics_allIter{1}.metrics.IV_allBins);
        pYES_pred_curve_iter = nan(nIter_fit, nBins_curve);
        pC_pred_curve_iter = nan(nIter_fit, nBins_curve);
        pA_pred_curve_iter = nan(nIter_fit, nBins_curve);

        for iIter = 1:nIter_fit
            pYES_pred_iter(iIter) = mean(pred_metrics_allIter{iIter}.metrics.pYES_pred_allBins, 'omitnan');
            pC_pred_iter(iIter) = mean(pred_metrics_allIter{iIter}.metrics.pC_pred_allBins, 'omitnan');
            pA_pred_iter(iIter) = mean(pred_metrics_allIter{iIter}.metrics.pA_pred_allBins, 'omitnan');

            pYES_pred_curve_iter(iIter, :) = pred_metrics_allIter{iIter}.metrics.pYES_pred_allBins;
            pC_pred_curve_iter(iIter, :) = pred_metrics_allIter{iIter}.metrics.pC_pred_allBins;
            pA_pred_curve_iter(iIter, :) = pred_metrics_allIter{iIter}.metrics.pA_pred_allBins;
        end

        [Ri.pYES_pred_med, Ri.pYES_pred_lb, Ri.pYES_pred_ub] = getCI(pYES_pred_iter, 1, 1);
        [Ri.pC_pred_med, Ri.pC_pred_lb, Ri.pC_pred_ub] = getCI(pC_pred_iter, 1, 1);
        [Ri.pA_pred_med, Ri.pA_pred_lb, Ri.pA_pred_ub] = getCI(pA_pred_iter, 1, 1);

        Ri.pYES_pred_curve_med = nan(1, nBins_curve);
        Ri.pC_pred_curve_med = nan(1, nBins_curve);
        Ri.pA_pred_curve_med = nan(1, nBins_curve);

        for iBin = 1:nBins_curve
            [Ri.pYES_pred_curve_med(iBin), ~, ~] = getCI(pYES_pred_curve_iter(:, iBin), 1, 1);
            [Ri.pC_pred_curve_med(iBin), ~, ~] = getCI(pC_pred_curve_iter(:, iBin), 1, 1);
            [Ri.pA_pred_curve_med(iBin), ~, ~] = getCI(pA_pred_curve_iter(:, iBin), 1, 1);
        end

        % ---------- parameter estimates ----------
        if iscell(data_fitNOM.params_est_allIter)
            P = cell2mat(cellfun(@(x) x(:)', data_fitNOM.params_est_allIter, 'UniformOutput', false));
        else
            P = data_fitNOM.params_est_allIter;
        end

        criterion_est_iter = P(:, end);
        criterion_rmse_iter = abs(criterion_est_iter - fileBase.criterion_DV_true);

        [Ri.criterion_DV_rmse, ~, ~] = getCI(criterion_rmse_iter, 1, 1);
        [Ri.criterion_DV_est_med, Ri.criterion_DV_est_lb, Ri.criterion_DV_est_ub] = ...
            getCI(criterion_est_iter, 1, 1);

        param_names_all = namesModelBparams_short{iModelB_fit};
        param_names_IN = param_names_all(1:end-1);

        for iParamIN = 1:numel(param_names_IN)
            pName = param_names_IN{iParamIN};
            trueVal = info.(sprintf('%s_true', pName));

            est_iter = P(:, iParamIN);
            err_iter = abs(est_iter - trueVal);

            [Ri.(sprintf('%s_rmse', pName)), ~, ~] = getCI(err_iter, 1, 1);
            [Ri.(sprintf('%s_est_med', pName)), ...
                Ri.(sprintf('%s_est_lb', pName)), ...
                Ri.(sprintf('%s_est_ub', pName))] = getCI(est_iter, 1, 1);
        end

        R = [R; Ri];
    end % iModelB_fit

end % iFile

fprintf('\n\n%s: All files compiled. \n\n', datetime('now'))

% Post-processing
% Apply filter on simulated pC
R_unfiltered = R;
R = R_unfiltered([R_unfiltered.keep_pC]);
% assert(nFilesKeep == numel(R)/nBfit)
fprintf('\nFiltering by simulated pC in [%.2f, %.2f]: kept %d / %d files.\n', mask_pC(1), mask_pC(2), nFilesKeep, nFiles);

% Define unique variable levels
gaborCST_unik = unique([R.gaborCST]);
Cz_unik = unique([R.Cz_true]);
Nmul_unik = unique([R.Nmul_true]);
Nadd_unik = unique([R.Nadd_true]);
Nshared_unik = unique([R.Nshared_true]);
Bsim_unik = unique([R.iModelB_sim]);
Bfit_unik = unique([R.iModelB_fit]);
lambda_whiten_unik = unique([R.lambda_whiten]);
C_contribution_unik = unique([R.C_contribution]);

% Save the compiled record to /Output
nameFile_R = sprintf('%s/Outputs/R_A%d.mat', nameFolder_server, iModelA_fit);
save(nameFile_R, 'R', 'R_unfiltered', 'lambda_whiten_unik', 'C_contribution_unik')

%% Metric list to plot
error
metricNames_plot = {'template_rmse', 'template_R2', ...
    'pYES_rmse', 'pYES_R2', ...
    'pC_rmse', 'pC_R2', ...
    'pA_rmse', 'pA_R2', ...
    'criterion_DV_rmse'};

param_names_all = namesModelBparams_short{1};
param_names_IN = param_names_all(1:end-1);

for iParam = 1:numel(param_names_IN)
    metricNames_plot{end+1} = sprintf('%s_rmse', param_names_IN{iParam}); %#ok<SAGROW>
end

fprintf('\n\n%s: Convert metrics list to plot. \n\n', datetime('now'))

%% Figure 1: heatmaps of RMSE and R2
fprintf('\n\n%s: Fig 1: Plot heatmaps of RMSE and R2. \n\n', datetime('now'))

for iModelB_fit = 1:numel(Bfit_unik)

    fitModel = Bfit_unik(iModelB_fit);

    idx_fit = [R.iModelB_fit] == fitModel;
    R_fit = R(idx_fit);

    if isempty(R_fit)
        continue
    end

    for iModelB_sim = 1:numel(Bsim_unik)
        fprintf('%s: Bsim=%d, Bfit=%d. \n', datetime('now'), iModelB_sim, iModelB_fit)
        simModel = Bsim_unik(iModelB_sim);

        idx_gen = [R_fit.iModelB_sim] == simModel;
        R_fit_gen = R_fit(idx_gen);

        if isempty(R_fit_gen)
            continue
        end

        for iMetric = 1:numel(metricNames_plot)
            metricName = metricNames_plot{iMetric};

            nRows = numel(gaborCST_unik);
            nCols = numel(Nadd_unik);

            % ---------- first pass: build all heatmaps and get global color range ----------
            M_allPanels = cell(nRows, nCols);
            allVals = [];

            for iSignalCST = 1:numel(gaborCST_unik)
                for iNadd = 1:numel(Nadd_unik)

                    thisG = gaborCST_unik(iSignalCST);
                    thisA = Nadd_unik(iNadd);

                    idx_sub = [R_fit_gen.gaborCST] == thisG & [R_fit_gen.Nadd_true] == thisA;
                    R_sub = R_fit_gen(idx_sub);

                    M = nan(numel(Nmul_unik), numel(Nshared_unik));

                    for iNmul = 1:numel(Nmul_unik)
                        for iNshared = 1:numel(Nshared_unik)
                            idx_cell = [R_sub.Nmul_true] == Nmul_unik(iNmul) & ...
                                [R_sub.Nshared_true] == Nshared_unik(iNshared);

                            if any(idx_cell)
                                vals = [R_sub(idx_cell).(metricName)];
                                M(iNmul, iNshared) = mean(vals, 'omitnan'); % average across Cz_true
                            end
                        end
                    end

                    M_allPanels{iSignalCST, iNadd} = M;
                    allVals = [allVals; M(isfinite(M))];
                end
            end

            if isempty(allVals)
                fprintf('No finite values found for %s, Bsim=%d, Bfit=%d\n', metricName, simModel, fitModel);
                continue
            end

            % ---------- shared color axis ----------
            if contains(metricName, 'R2', 'IgnoreCase', true)
                cMin = 0;
                cMax = 1;
            else
                cMin = min(allVals);
                cMax = max(allVals);

                if cMin == cMax
                    cMin = cMin - eps;
                    cMax = cMax + eps;
                end
            end

            % ---------- second pass: plot using shared color axis ----------
            h = figure('Position', [100 100 350*nCols 300*nRows]);
            tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

            for iSignalCST = 1:numel(gaborCST_unik)
                for iNadd = 1:numel(Nadd_unik)

                    nexttile; hold on;

                    thisG = gaborCST_unik(iSignalCST);
                    thisA = Nadd_unik(iNadd);
                    M = M_allPanels{iSignalCST, iNadd};

                    imagesc(1:numel(Nshared_unik), 1:numel(Nmul_unik), M);
                    set(gca, 'YDir', 'normal', ...
                        'XTick', 1:numel(Nshared_unik), ...
                        'XTickLabel', string(Nshared_unik), ...
                        'YTick', 1:numel(Nmul_unik), ...
                        'YTickLabel', string(Nmul_unik), ...
                        'FontSize', 10);

                    axis square;
                    clim([cMin cMax]);

                    % annotate values
                    for iNmul = 1:numel(Nmul_unik)
                        for iNshared = 1:numel(Nshared_unik)
                            if isfinite(M(iNmul, iNshared))

                                frac = (M(iNmul, iNshared) - cMin) / (cMax - cMin);
                                frac = max(0, min(1, frac));

                                if frac > 0.6
                                    txtColor = 'w';
                                else
                                    txtColor = 'k';
                                end

                                text(iNshared, iNmul, sprintf('%.2f', M(iNmul, iNshared)), ...
                                    'HorizontalAlignment', 'center', ...
                                    'VerticalAlignment', 'middle', ...
                                    'FontSize', 8, ...
                                    'Color', txtColor);
                            end
                        end % iNshared
                    end % iNmul

                    xlabel('Nshared true');
                    ylabel('Nmul true');
                    title(sprintf('gabor = %.1f, Nadd = %g', thisG, thisA));
                end % iNadd
            end % iSignalCST

            sgtitle(sprintf('[%s] | Bsim=%d | Bfit=%d (%s) | mean across Cz', ...
                metricName, simModel, fitModel, namesModelB{fitModel}), ...
                'FontSize', 14, 'FontWeight', 'bold');

            cb = colorbar;
            cb.Layout.Tile = 'east';
            cb.Label.String = metricName;

            saveas(h, fullfile(nameFolder_Figures_part3, sprintf('Fig1_%s_B%d%d.png', metricName, simModel, fitModel)));
            close(h);
        end % iMetric
    end % iModelB_sim
end % iModelB_fit

%% Figure 2: Scatter plots: true vs. estimated values
fprintf('\n\n%s: Fig 2: Scatter plots of true vs. estimated values. \n\n', datetime('now'))

metricNames_scatter = {'pYES', 'pC', 'pA'};
marker_list = {'o','s','^','d','v','>','<','p','h'};

for iModelB_fit = 1:numel(Bfit_unik)

    fitModel = Bfit_unik(iModelB_fit);

    idx_fit = [R.iModelB_fit] == fitModel;
    R_fit = R(idx_fit);

    param_names_all = namesModelBparams_short{fitModel}; % e.g., {'Nmul','Nadd','Nshared','criterion_DV'}
    param_names_IN = param_names_all(1:end-1);

    % ---------- 1. behavioral metrics ----------
    for iMetric = 1:numel(metricNames_scatter)
        mName = metricNames_scatter{iMetric};

        for iModelB_sim = 1:numel(Bsim_unik)
            fprintf('%s: Bsim=%d, Bfit=%d. \n', datetime('now'), iModelB_sim, iModelB_fit)

            simModel = Bsim_unik(iModelB_sim);

            nRows = numel(gaborCST_unik);
            nCols = 1;

            h = figure('Position', [100 100 420*nCols 300*nRows]);
            tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

            for iSignalCST = 1:numel(gaborCST_unik)

                nexttile; hold on;

                idx = [R_fit.gaborCST] == gaborCST_unik(iSignalCST) & [R_fit.iModelB_sim] == simModel;
                Rsub = R_fit(idx);

                x_true = [Rsub.(sprintf('%s_data_med', mName))]';
                y_est = [Rsub.(sprintf('%s_pred_med', mName))]';

                Nshared_vals = [Rsub.Nshared_true]';
                Nmul_vals = [Rsub.Nmul_true]';
                Nadd_vals = [Rsub.Nadd_true]';

                fxn_plotScatterPanel(x_true, y_est, Nshared_vals, Nmul_vals, Nadd_vals, marker_list);

                if strcmp(mName, 'pC'), xlim([.6, .8]), ylim([.6, .8]), end

                title(sprintf('gabor = %.1f | Bsim = %d', gaborCST_unik(iSignalCST), simModel));
                xlabel(sprintf('Simulated %s', mName));
                ylabel(sprintf('Predicted %s', mName));
            end % iSignalCST

            sgtitle(sprintf('%s recovery | Bsim = %d | Bfit = %d (%s)', ...
                mName, simModel, fitModel, namesModelB{fitModel}), ...
                'FontSize', 14, 'FontWeight', 'bold');

            fxn_addScatterLegends_oneFigure(Nshared_unik, Nmul_unik, Nadd_unik, marker_list)

            saveas(h, fullfile(nameFolder_Figures_part3, ...
                sprintf('Fig2_%s_B%d%d.png', mName, simModel, fitModel)));
            close(h);
        end % iModelB_sim
    end % iMetric

    % ---------- 2. criterion ----------
    for iModelB_sim = 1:numel(Bsim_unik)
        simModel = Bsim_unik(iModelB_sim);

        nRows = numel(gaborCST_unik);
        nCols = 1;

        h = figure('Position', [100 100 420*nCols 300*nRows]);
        tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

        for iSignalCST = 1:numel(gaborCST_unik)

            nexttile; hold on;

            idx = [R_fit.gaborCST] == gaborCST_unik(iSignalCST) & [R_fit.iModelB_sim] == simModel;
            Rsub = R_fit(idx);

            x_true = [Rsub.criterion_DV_true]';
            y_est = [Rsub.criterion_DV_est_med]';

            Nshared_vals = [Rsub.Nshared_true]';
            Nmul_vals = [Rsub.Nmul_true]';
            Nadd_vals = [Rsub.Nadd_true]';

            fxn_plotScatterPanel(x_true, y_est, Nshared_vals, Nmul_vals, Nadd_vals, marker_list);

            title(sprintf('gabor = %.1f | Bsim = %d', gaborCST_unik(iSignalCST), simModel));
            xlabel('True criterion');
            ylabel('Estimated criterion');
        end % iSignalCST

        sgtitle(sprintf('Criterion recovery | Bsim = %d | Bfit = %d (%s)', ...
            simModel, fitModel, namesModelB{fitModel}), ...
            'FontSize', 14, 'FontWeight', 'bold');

        fxn_addScatterLegends_oneFigure(Nshared_unik, Nmul_unik, Nadd_unik, marker_list)

        saveas(h, fullfile(nameFolder_Figures_part3, ...
            sprintf('Fig2_criterionDV_B%d%d.png', simModel, fitModel)));
        close(h);
    end % iModelB_sim

    % ---------- 3. IN parameters ----------
    for iParam = 1:numel(param_names_IN)
        pName = param_names_IN{iParam}; % 'Nmul', 'Nadd', or 'Nshared'

        for iModelB_sim = 1:numel(Bsim_unik)
            simModel = Bsim_unik(iModelB_sim);

            nRows = numel(gaborCST_unik);
            nCols = 1;

            h = figure('Position', [100 100 420*nCols 300*nRows]);
            tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

            for iSignalCST = 1:numel(gaborCST_unik)

                nexttile; hold on;

                idx = [R_fit.gaborCST] == gaborCST_unik(iSignalCST) & [R_fit.iModelB_sim] == simModel;
                Rsub = R_fit(idx);

                x_true = [Rsub.(sprintf('%s_true', pName))]';
                y_est = [Rsub.(sprintf('%s_est_med', pName))]';

                Nshared_vals = [Rsub.Nshared_true]';
                Nmul_vals = [Rsub.Nmul_true]';
                Nadd_vals = [Rsub.Nadd_true]';

                fxn_plotScatterPanel(x_true, y_est, Nshared_vals, Nmul_vals, Nadd_vals, marker_list);

                title(sprintf('gabor = %.1f | Bsim = %d', gaborCST_unik(iSignalCST), simModel));
                xlabel(sprintf('True %s', pName));
                ylabel(sprintf('Estimated %s', pName));
            end % iSignalCST

            sgtitle(sprintf('%s recovery | Bsim = %d | Bfit = %d (%s)', pName, simModel, fitModel, namesModelB{fitModel}), ...
                'FontSize', 14, 'FontWeight', 'bold');

            fxn_addScatterLegends_oneFigure(Nshared_unik, Nmul_unik, Nadd_unik, marker_list)

            saveas(h, fullfile(nameFolder_Figures_part3, sprintf('Fig2_%s_B%d%d.png', pName, simModel, fitModel)));
            close(h);
        end % iModelB_sim
    end % iParam
end % iModelB_fit

%% Figure 3: Simulated & estimated tuning functions
% Note that ModelB_fit is not looped
fprintf('\n\n%s: Fig 3: Marginalized tuning functions. \n\n', datetime('now'))

% style for within-panel lines
cmap_Cz = [1,0,0; 0,0,0; 0,0,1];
lineStyles_signal = {'-', '--', ':'};
if numel(gaborCST_unik) > numel(lineStyles_signal)
    error('Not enough line styles for signal contrast levels.');
end

for iModelB_sim = 1:numel(Bsim_unik)
    fprintf('%s: Bsim=%d. \n', datetime('now'), iModelB_sim)
    simModel = Bsim_unik(iModelB_sim);

    for iFeature = 1:2

        nameFeature = namesFeature{iFeature};
        xAxis = axis_tuning{iFeature};

        for iNshared = 1:numel(Nshared_unik)
            thisNshared = Nshared_unik(iNshared);

            nRows = numel(Nmul_unik);
            nCols = numel(Nadd_unik);

            h = figure('Position', [100 100 350*nCols 260*nRows]);
            tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

            for iNmul = 1:numel(Nmul_unik)
                for iNadd = 1:numel(Nadd_unik)

                    nexttile; hold on;

                    thisNmul = Nmul_unik(iNmul);
                    thisNadd = Nadd_unik(iNadd);

                    idx = [R.iModelB_sim] == simModel & ...
                        [R.Nshared_true] == thisNshared & ...
                        [R.Nmul_true] == thisNmul & ...
                        [R.Nadd_true] == thisNadd;

                    Rsub = R(idx);

                    if isempty(Rsub)
                        title(sprintf('Nmul = %g | Nadd = %g\n(no data)', thisNmul, thisNadd));
                        box on;
                        continue
                    end

                    % plot all 9 estimated curves (3 signal × 3 Cz)
                    for iSignalCST = 1:numel(gaborCST_unik)
                        for iCz = 1:numel(Cz_unik)

                            idx_line = [Rsub.gaborCST] == gaborCST_unik(iSignalCST) & ...
                                [Rsub.Cz_true] == Cz_unik(iCz);

                            if ~any(idx_line)
                                continue
                            end

                            Ri = Rsub(find(idx_line, 1, 'first'));

                            if iFeature == 1
                                y_est = Ri.margORI_est_med;
                            else
                                y_est = Ri.margSF_est_med;
                            end

                            plot(xAxis, y_est, ...
                                'Color', cmap_Cz(iCz,:), ...
                                'LineStyle', lineStyles_signal{iSignalCST}, ...
                                'LineWidth', 1.2);
                        end
                    end

                    % overlay one true curve (same for all 9 conditions here)
                    Ri0 = Rsub(1);
                    if iFeature == 1
                        y_true = Ri0.margORI_true;
                    else
                        y_true = Ri0.margSF_true;
                    end

                    plot(xAxis, y_true, 'k-', 'LineWidth', 2.5);

                    box on;
                    set(gca, 'FontSize', 10);

                    if iFeature == 1
                        xlabel('Orientation channel');
                    else
                        xlabel('Spatial frequency channel');
                    end
                    ylabel('Normalized tuning');

                    ylim([-.5, 1])
                    title(sprintf('Nmul = %g | Nadd = %g', thisNmul, thisNadd));
                end
            end % iNmul

            sgtitle(sprintf('Fig 3: %s tuning | Bsim = %d | Nshared = %g', nameFeature, simModel, thisNshared), ...
                'FontSize', 14, 'FontWeight', 'bold');

            % add figure-level legends
            fxn_addTuningLegends_oneFigure(Cz_unik, cmap_Cz, gaborCST_unik, lineStyles_signal)

            saveas(h, fullfile(nameFolder_Figures_part3, sprintf('Fig3_%s_Nshared%g_Bsim%d.png', nameFeature, thisNshared, simModel)));
            close(h);
        end % iNshared
    end % iFeature
end % iModelB_sim

%% Figure 4: simulated vs predicted metrics
fprintf('\n\n%s: Fig 4: Simulated vs predicted metrics.\n\n', datetime('now'))

metricNames_curve = {'pYES', 'pC', 'pA'};

% within-panel style
cmap_Cz = [1,0,0; 0,0,0; 0,0,1]; % 3 Cz levels
marker_list_signal = {'o', 's', '^'}; % 3 signal contrasts
lineStyles_fit = {'-', '--', ':', '-.'}; % 4 fitted models

if numel(gaborCST_unik) > numel(marker_list_signal)
    error('Not enough marker styles for signal contrast levels.');
end
if numel(Bfit_unik) > numel(lineStyles_fit)
    error('Not enough line styles for fitted models.');
end

for iModelB_sim = 1:numel(Bsim_unik)
    fprintf('%s: Bsim=%d. \n', datetime('now'), iModelB_sim)
    simModel = Bsim_unik(iModelB_sim);

    for iNshared = 1:numel(Nshared_unik)
        thisNshared = Nshared_unik(iNshared);

        nRows = numel(Nmul_unik);
        nCols = 3 * numel(Nadd_unik); % 3 metric blocks

        h = figure('Position', [100 100 260*nCols 220*nRows]);
        tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

        for iNmul = 1:numel(Nmul_unik)
            thisNmul = Nmul_unik(iNmul);

            for iMetric = 1:numel(metricNames_curve)
                mName = metricNames_curve{iMetric};

                for iNadd = 1:numel(Nadd_unik)
                    thisNadd = Nadd_unik(iNadd);

                    nexttile; hold on;

                    idx_panel = [R.iModelB_sim] == simModel & ...
                        [R.Nshared_true] == thisNshared & ...
                        [R.Nmul_true] == thisNmul & ...
                        [R.Nadd_true] == thisNadd;

                    Rsub = R(idx_panel);

                    if isempty(Rsub)
                        title(sprintf('%s | Nadd=%g\n(no data)', mName, thisNadd));
                        box on;
                        continue
                    end

                    yDataField = sprintf('%s_data_curve_med', mName);
                    yPredField = sprintf('%s_pred_curve_med', mName);

                    % plot 3 criteria x 3 signal contrasts
                    for iCz = 1:numel(Cz_unik)
                        for iSignal = 1:numel(gaborCST_unik)

                            idx_cond = [Rsub.Cz_true] == Cz_unik(iCz) & ...
                                [Rsub.gaborCST] == gaborCST_unik(iSignal);

                            if ~any(idx_cond)
                                continue
                            end

                            % -------- simulated dots --------
                            % choose first matching row because data curve is the same across Bfit
                            Ri_data = Rsub(find(idx_cond, 1, 'first'));

                            if ~isfield(Ri_data, 'DV_allBins_med') || ~isfield(Ri_data, yDataField)
                                warning('Missing DV/curve fields in R.');
                                continue
                            end

                            xDV = Ri_data.DV_allBins_med;
                            yData = Ri_data.(yDataField);

                            thisColor = cmap_Cz(iCz, :);
                            thisMarker = marker_list_signal{iSignal};

                            plot(xDV, yData, ...
                                'LineStyle', 'none', ...
                                'Marker', thisMarker, ...
                                'MarkerSize', 4, ...
                                'MarkerFaceColor', thisColor, ...
                                'MarkerEdgeColor', thisColor);

                            % -------- predicted lines: one per Bfit --------
                            for iModelB_fit = 1:numel(Bfit_unik)
                                fitModel = Bfit_unik(iModelB_fit);

                                idx_pred = idx_cond & [Rsub.iModelB_fit] == fitModel;
                                if ~any(idx_pred)
                                    continue
                                end

                                Ri_pred = Rsub(find(idx_pred, 1, 'first'));

                                if ~isfield(Ri_pred, yPredField)
                                    continue
                                end

                                yPred = Ri_pred.(yPredField);

                                plot(xDV, yPred, ...
                                    'Color', thisColor, ...
                                    'LineStyle', lineStyles_fit{iModelB_fit}, ...
                                    'LineWidth', 1.0);
                            end % iFit
                        end % iSignal
                    end % iCz

                    box on;
                    set(gca, 'FontSize', 9);

                    % only label outer axes to reduce clutter
                    if iNmul == nRows
                        xlabel('Binned DV');
                    end
                    if iMetric == 1 && iNadd == 1
                        ylabel(sprintf('Nmul = %g', thisNmul));
                    end

                    title(sprintf('%s | Nadd = %g', mName, thisNadd));
                    switch mName
                        case 'pYES'
                            ylim([0,1])
                        case {'pC','pA'}
                            ylim([0.5,1])
                    end
                end % iNadd
            end % iMetric
        end % iNmul

        sgtitle(sprintf('Fig 4: simulated vs predicted metrics | Bsim = %d | Nshared = %g', ...
            simModel, thisNshared), ...
            'FontSize', 14, 'FontWeight', 'bold');

        fxn_addMetricCurveLegends_oneFigure(Cz_unik, cmap_Cz, gaborCST_unik, marker_list_signal, Bfit_unik, lineStyles_fit)

        saveas(h, fullfile(nameFolder_Figures_part3, ...
            sprintf('Fig4_metrics_B%d_Nshared%g.png', simModel, thisNshared)));
        close(h);
    end % iNshared
end % iModelB_sim

%% Figure 5: model recovery
fprintf('\n\n%s: Fig 5: Model recovery (panels are GaborCST x True Cz).\n\n', datetime('now'))

scoreField = 'nLL_med'; % <-- change if needed

if ~isfield(R, scoreField)
    error('R does not contain the field "%s". Add it during compile first.', scoreField);
end

%% ---------- 5A. Win-rate matrix ----------
fprintf('%s: Win-rate matrix. \n', datetime('now'))
h = figure('Position', [100 100 320*numel(Cz_unik) 280*numel(gaborCST_unik)]);
tiledlayout(numel(gaborCST_unik), numel(Cz_unik), 'TileSpacing', 'compact', 'Padding', 'compact');

for iSignal = 1:numel(gaborCST_unik)
    for iCz = 1:numel(Cz_unik)

        nexttile; hold on;

        thisG = gaborCST_unik(iSignal);
        thisCz = Cz_unik(iCz);

        idx_panel = [R.gaborCST] == thisG & [R.Cz_true] == thisCz;
        Rsub = R(idx_panel);

        W = nan(numel(Bfit_unik), numel(Bsim_unik));

        for iModelB_sim = 1:numel(Bsim_unik)
            simModel = Bsim_unik(iModelB_sim);

            idx_sim = [Rsub.iModelB_sim] == simModel;
            Rsim = Rsub(idx_sim);

            if isempty(Rsim)
                continue
            end

            % unique condition groups, collapsing across fitted model
            % panel already fixed gaborCST and Cz_true, so key uses only Nmul/Nadd/Nshared
            condKeys = arrayfun(@(x) sprintf('M%.3f_A%.3f_S%.3f', ...
                x.Nmul_true, x.Nadd_true, x.Nshared_true), ...
                Rsim, 'UniformOutput', false);
            unikKeys = unique(condKeys);

            winCount = zeros(numel(Bfit_unik), 1);
            nCond = 0;

            for iKey = 1:numel(unikKeys)
                idx_key = strcmp(condKeys, unikKeys{iKey});
                Rkey = Rsim(idx_key);

                if numel(Rkey) < 2
                    continue
                end

                scores = nan(numel(Bfit_unik), 1);
                for iModelB_fit = 1:numel(Bfit_unik)
                    fitModel = Bfit_unik(iModelB_fit);
                    idx_fit = [Rkey.iModelB_fit] == fitModel;
                    if any(idx_fit)
                        scores(iModelB_fit) = Rkey(find(idx_fit, 1, 'first')).(scoreField);
                    end
                end

                if all(isnan(scores))
                    continue
                end

                [~, iBest] = min(scores);
                winCount(iBest) = winCount(iBest) + 1;
                nCond = nCond + 1;
            end

            if nCond > 0
                W(:, iModelB_sim) = winCount / nCond;
            end
        end

        imagesc(1:numel(Bsim_unik), 1:numel(Bfit_unik), W);
        set(gca, 'YDir', 'normal', ...
            'XTick', 1:numel(Bsim_unik), ...
            'XTickLabel', string(Bsim_unik), ...
            'YTick', 1:numel(Bfit_unik), ...
            'YTickLabel', string(Bfit_unik), ...
            'FontSize', 10);
        axis square;
        clim([0 1]);

        for iModelB_fit = 1:numel(Bfit_unik)
            for iModelB_sim = 1:numel(Bsim_unik)
                if isfinite(W(iModelB_fit, iModelB_sim))
                    txtColor = 'k';
                    if W(iModelB_fit, iModelB_sim) > 0.6
                        txtColor = 'w';
                    end
                    text(iModelB_sim, iModelB_fit, sprintf('%.2f', W(iModelB_fit, iModelB_sim)), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'FontSize', 8, ...
                        'Color', txtColor);
                end
            end
        end

        xlabel('Bsim');
        ylabel('Bfit');
        title(sprintf('gabor = %.1f | Cz = %g', thisG, thisCz));
    end
end

sgtitle('Fig 5A: Model recovery | win rate', ...
    'FontSize', 14, 'FontWeight', 'bold');

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Win rate';

saveas(h, fullfile(nameFolder_Figures_part3, sprintf('Fig5A_ModelRecovery_WinRate_%s.png', scoreField)));
close(h);

%% ---------- 5B. Median delta-score matrix ----------
fprintf('%s: dnLL matrix. \n', datetime('now'))
h = figure('Position', [100 100 320*numel(Cz_unik) 280*numel(gaborCST_unik)]);
tiledlayout(numel(gaborCST_unik), numel(Cz_unik), 'TileSpacing', 'compact', 'Padding', 'compact');

dnLL_all = cell(numel(gaborCST_unik), numel(Cz_unik));
allVals = [];

for iSignal = 1:numel(gaborCST_unik)
    for iCz = 1:numel(Cz_unik)

        thisG = gaborCST_unik(iSignal);
        thisCz = Cz_unik(iCz);

        idx_panel = [R.gaborCST] == thisG & [R.Cz_true] == thisCz;
        Rsub = R(idx_panel);

        D = nan(numel(Bfit_unik), numel(Bsim_unik));

        for iModelB_sim = 1:numel(Bsim_unik)
            simModel = Bsim_unik(iModelB_sim);

            idx_sim = [Rsub.iModelB_sim] == simModel;
            Rsim = Rsub(idx_sim);

            if isempty(Rsim)
                continue
            end

            condKeys = arrayfun(@(x) sprintf('M%.3f_A%.3f_S%.3f', ...
                x.Nmul_true, x.Nadd_true, x.Nshared_true), ...
                Rsim, 'UniformOutput', false);
            unikKeys = unique(condKeys);

            deltaVals = cell(numel(Bfit_unik), 1);

            for iKey = 1:numel(unikKeys)
                idx_key = strcmp(condKeys, unikKeys{iKey});
                Rkey = Rsim(idx_key);

                scores = nan(numel(Bfit_unik), 1);
                for iModelB_fit = 1:numel(Bfit_unik)
                    fitModel = Bfit_unik(iModelB_fit);
                    idx_fit = [Rkey.iModelB_fit] == fitModel;
                    if any(idx_fit)
                        scores(iModelB_fit) = Rkey(find(idx_fit, 1, 'first')).(scoreField);
                    end
                end

                if all(isnan(scores))
                    continue
                end

                bestScore = min(scores);

                for iModelB_fit = 1:numel(Bfit_unik)
                    if isfinite(scores(iModelB_fit))
                        deltaVals{iModelB_fit}(end+1,1) = scores(iModelB_fit) - bestScore; %#ok<AGROW>
                    end
                end
            end

            for iModelB_fit = 1:numel(Bfit_unik)
                if ~isempty(deltaVals{iModelB_fit})
                    D(iModelB_fit, iModelB_sim) = median(deltaVals{iModelB_fit}, 'omitnan');
                end
            end
        end

        dnLL_all{iSignal, iCz} = D;
        allVals = [allVals; D(isfinite(D))];
    end % iCz
end % iSignal


% Determine ub and lb of color axis
cMin = min(allVals);
cMin = max(0, cMin);
cMax = 150; % cap display range

if cMin >= cMax
    cMin = 0;
    cMax = 150;
end

for iSignal = 1:numel(gaborCST_unik)
    for iCz = 1:numel(Cz_unik)

        nexttile; hold on;

        D = dnLL_all{iSignal, iCz};

        imagesc(1:numel(Bsim_unik), 1:numel(Bfit_unik), D);
        set(gca, 'YDir', 'normal', ...
            'XTick', 1:numel(Bsim_unik), ...
            'XTickLabel', string(Bsim_unik), ...
            'YTick', 1:numel(Bfit_unik), ...
            'YTickLabel', string(Bfit_unik), ...
            'FontSize', 10);
        axis square;
        clim([cMin cMax]);

        for iModelB_fit = 1:numel(Bfit_unik)
            for iModelB_sim = 1:numel(Bsim_unik)
                if isfinite(D(iModelB_fit, iModelB_sim))
                    frac = (D(iModelB_fit, iModelB_sim) - cMin) / (cMax - cMin);
                    frac = max(0, min(1, frac));
                    txtColor = 'k';
                    if frac > 0.6
                        txtColor = 'w';
                    end
                    text(iModelB_sim, iModelB_fit, sprintf('%.2f', D(iModelB_fit, iModelB_sim)), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'FontSize', 8, ...
                        'Color', txtColor);
                end
            end % iModelB_sim
        end % iModelB_fit

        xlabel('Bsim');
        ylabel('Bfit');
        title(sprintf('gabor = %.1f | Cz = %g', gaborCST_unik(iSignal), Cz_unik(iCz)));
    end % iCz
end

sgtitle('Fig 5B: Model recovery | median ΔnLL', ...
    'FontSize', 14, 'FontWeight', 'bold');

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = sprintf('Δ%s from best fit', scoreField);

saveas(h, fullfile(nameFolder_Figures_part3, sprintf('Fig5B_ModelRecovery_Delta%s.png', scoreField)));
close(h);


%% Local functions

function info = fxn_parse_nameIO(nameIO)

info = struct();

tok = regexp(nameIO, ...
    'IO_cN([-\d\.]+)_cG([-\d\.]+)_nT([A-Za-z0-9\.\-]+)_Nm([A-Za-z0-9\.\-]+)_Na([A-Za-z0-9\.\-]+)_Ns([A-Za-z0-9\.\-]+)_Cz([-\d\.]+)_cont([-\d\.]+)_whiten([-\d\.]+)_R(\d+)_\d+\d+_B(\d+)', ...
    'tokens', 'once');

if isempty(tok)
    return
end

info.noiseCST = str2double(tok{1}) / 100;
info.gaborCST = str2double(tok{2}) / 100;
info.nTrials = fxn_parse_num2exp(tok{3});
info.Nmul_true = fxn_parse_num2exp(tok{4});
info.Nadd_true = fxn_parse_num2exp(tok{5});
info.Nshared_true = fxn_parse_num2exp(tok{6});
info.Cz_true = str2double(tok{7});
info.C_contribution = str2double(tok{8});
info.lambda_whiten = str2double(tok{9});
info.flag_regressType = str2double(tok{10});
info.iModelB_sim = str2double(tok{11});
end

%%
function x = fxn_parse_num2exp(str_in)
% parse strings like '8e3', '0', '10', etc.
x = str2double(strrep(str_in, 'p', '.'));
end

%%
function [template_rmse_iter, template_R2_iter] = fxn_templateRecovery_simPlot2Style_iter(truth, data_compIV, nIter)

if isfield(truth, 'template_true_raw')
    template_true = truth.template_true_raw(:)';
elseif isfield(truth, 'template_true')
    template_true = truth.template_true(:)';
else
    error('No template_true found in truth.mat');
end

template_true = template_true / max(template_true);
template_est = data_compIV.template_tmpl_allIter;
if ndims(template_est) == 3
    template_est = reshape(template_est, size(template_est,1), []);
end

rowMax = max(template_est, [], 2);
rowMax(rowMax == 0) = 1;
template_est = template_est ./ rowMax;

template_rmse_iter = sqrt(mean((template_est - template_true).^2, 2));

template_R2_iter = nan(nIter,1);
yTrue = template_true(:);
yTrue = yTrue/max(yTrue);

for iIter = 1:nIter
    yEst = template_est(iIter,:)';
    yEst = yEst/max(yEst);
    sse = sum((yTrue - yEst).^2);
    sst = sum((yTrue - mean(yTrue)).^2);
    template_R2_iter(iIter) = 1 - sse / sst;
    % template_R2_iter(iIter) = corr(yTrue, yEst);
end
end

%%
function [rmse_iter, R2_iter] = fxn_metricRecovery_iter(pred_metrics_allIter, metricName)

nIter = numel(pred_metrics_allIter);
rmse_iter = nan(nIter,1);
R2_iter = nan(nIter,1);

for iIter = 1:nIter
    yData = pred_metrics_allIter{iIter}.metrics.([metricName '_data_allBins'])(:);
    yPred = pred_metrics_allIter{iIter}.metrics.([metricName '_pred_allBins'])(:);
    nTrials = pred_metrics_allIter{iIter}.metrics.nTrials_allBins(:);

    good = isfinite(yData) & isfinite(yPred) & isfinite(nTrials) & (nTrials > 0);
    yData = yData(good);
    yPred = yPred(good);
    w = nTrials(good);

    if isempty(yData)
        continue
    end

    rmse_iter(iIter) = sqrt(sum(w .* (yData - yPred).^2) / sum(w));

    sse = sum((yData - yPred).^2);
    sst = sum((yData - mean(yData)).^2);
    R2_iter(iIter) = 1 - sse / sst;
end
end

%%
function fxn_plotScatterPanel(x_true, y_est, Nshared_vals, Nmul_vals, Nadd_vals, marker_list)

good = isfinite(x_true) & isfinite(y_est) & isfinite(Nshared_vals) & isfinite(Nmul_vals) & isfinite(Nadd_vals);
x_true = x_true(good);
y_est = y_est(good);
Nshared_vals = Nshared_vals(good);
Nmul_vals = Nmul_vals(good);
Nadd_vals = Nadd_vals(good);

if isempty(x_true)
    box on;
    return
end

% map style
Nshared_unik = unique(Nshared_vals);
Nmul_unik = unique(Nmul_vals);
Nadd_unik = unique(Nadd_vals);

cmap = parula(numel(Nshared_unik));

if numel(Nadd_unik) > numel(marker_list)
    error('Not enough marker styles for Nadd levels.');
end

for iPt = 1:numel(x_true)
    iCol = find(Nshared_unik == Nshared_vals(iPt), 1, 'first');
    faceColor = cmap(iCol, :);

    if numel(Nmul_unik) == 1
        gray_edge = 0.4;
    else
        frac = (Nmul_vals(iPt) - min(Nmul_unik)) / (max(Nmul_unik) - min(Nmul_unik));
        gray_edge = 0.85 - 0.65 * frac;
    end
    edgeColor = gray_edge * [1 1 1];

    iMark = find(Nadd_unik == Nadd_vals(iPt), 1, 'first');
    this_marker = marker_list{iMark};

    scatter(x_true(iPt), y_est(iPt), 55, ...
        'Marker', this_marker, ...
        'MarkerFaceColor', faceColor, ...
        'MarkerEdgeColor', edgeColor, ...
        'LineWidth', 1.0);
end

xy_all = [x_true(:); y_est(:)];
xy_all = xy_all(isfinite(xy_all));
xy_min = min(xy_all);
xy_max = max(xy_all);

if xy_min == xy_max
    xy_min = xy_min - 0.1;
    xy_max = xy_max + 0.1;
end

plot([xy_min xy_max], [xy_min xy_max], 'k--', 'LineWidth', 1);

axis square;
xlim([xy_min xy_max]);
ylim([xy_min xy_max]);
box on;
set(gca, 'FontSize', 11);
end

%%
function fxn_addScatterLegends_oneFigure(Nshared_unik, Nmul_unik, Nadd_unik, marker_list)

axMain = gca;

% face color: Nshared
cmap_shared = parula(numel(Nshared_unik));
makeLegend(axMain, ...
    arrayfun(@(i) scatter(nan, nan, 60, ...
    'Marker', 'o', ...
    'MarkerFaceColor', cmap_shared(i,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1), ...
    1:numel(Nshared_unik), 'UniformOutput', false), ...
    arrayfun(@(x) sprintf('%g', x), Nshared_unik, 'UniformOutput', false), ...
    'northoutside', 'horizontal', 'Nshared');

% edge grayness: Nmul
mul_min = min(Nmul_unik);
mul_max = max(Nmul_unik);
makeLegend(newInvisibleAxes(axMain), ...
    arrayfun(@(x) scatter(nan, nan, 60, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', getGrayEdge(x, mul_min, mul_max), ...
    'LineWidth', 1.5), ...
    Nmul_unik, 'UniformOutput', false), ...
    arrayfun(@(x) sprintf('%g', x), Nmul_unik, 'UniformOutput', false), ...
    'southoutside', 'horizontal', 'Nmul');

% marker shape: Nadd
if numel(Nadd_unik) > numel(marker_list)
    error('Not enough marker styles for Nadd levels.');
end

makeLegend(newInvisibleAxes(axMain), ...
    arrayfun(@(i) scatter(nan, nan, 60, ...
    'Marker', marker_list{i}, ...
    'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2), ...
    1:numel(Nadd_unik), 'UniformOutput', false), ...
    arrayfun(@(x) sprintf('%g', x), Nadd_unik, 'UniformOutput', false), ...
    'eastoutside', 'vertical', 'Nadd');

end

%%
function ax = newInvisibleAxes(axRef)
ax = axes('Position', axRef.Position, 'Color', 'none', 'Visible', 'off');
hold(ax, 'on');
end

function c = getGrayEdge(x, xmin, xmax)
if xmax == xmin
    gray_edge = 0.4;
else
    frac = (x - xmin) / (xmax - xmin);
    gray_edge = 0.85 - 0.65 * frac;
end
c = gray_edge * [1 1 1];
end

function makeLegend(ax, hCell, labels, loc, orient, ttl)
h = [hCell{:}];
leg = legend(ax, h, labels, 'Location', loc, 'Orientation', orient);
title(leg, ttl);
end

function fxn_addTuningLegends_oneFigure(Cz_unik, cmap_Cz, gabor_unik, lineStyles_signal)

axMain = gca;

% Legend 1: color = Cz
hCz = gobjects(numel(Cz_unik), 1);
for i = 1:numel(Cz_unik)
    hCz(i) = plot(nan, nan, ...
        'Color', cmap_Cz(i,:), ...
        'LineStyle', '-', ...
        'LineWidth', 2);
end

leg1 = legend(axMain, hCz, ...
    arrayfun(@(x) sprintf('%g', x), Cz_unik, 'UniformOutput', false), ...
    'Location', 'northoutside', ...
    'Orientation', 'horizontal');
title(leg1, 'Cz');

% Legend 2: line style = signal contrast
ax2 = axes('Position', axMain.Position, 'Color', 'none', 'Visible', 'off');
hold(ax2, 'on');

hSig = gobjects(numel(gabor_unik), 1);
for i = 1:numel(gabor_unik)
    hSig(i) = plot(ax2, nan, nan, ...
        'Color', 'k', ...
        'LineStyle', lineStyles_signal{i}, ...
        'LineWidth', 2);
end

leg2 = legend(ax2, hSig, ...
    arrayfun(@(x) sprintf('%g', x), gabor_unik, 'UniformOutput', false), ...
    'Location', 'eastoutside');
title(leg2, 'Signal contrast');

end

function fxn_addMetricCurveLegends_oneFigure(Cz_unik, cmap_Cz, gabor_unik, marker_list_signal, ModelB_fit_unik, lineStyles_fit)

axMain = gca;

% -------- Legend 1: color = Cz --------
hCz = gobjects(numel(Cz_unik), 1);
for i = 1:numel(Cz_unik)
    hCz(i) = plot(nan, nan, ...
        'Color', cmap_Cz(i,:), ...
        'LineStyle', '-', ...
        'LineWidth', 2);
end
leg1 = legend(axMain, hCz, ...
    arrayfun(@(x) sprintf('%g', x), Cz_unik, 'UniformOutput', false), ...
    'Location', 'northoutside', 'Orientation', 'horizontal');
title(leg1, 'Cz');

% -------- Legend 2: marker = signal contrast (dots) --------
ax2 = axes('Position', axMain.Position, 'Color', 'none', 'Visible', 'off');
hold(ax2, 'on');

hSig = gobjects(numel(gabor_unik), 1);
for i = 1:numel(gabor_unik)
    hSig(i) = plot(ax2, nan, nan, ...
        'LineStyle', 'none', ...
        'Marker', marker_list_signal{i}, ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k');
end
leg2 = legend(ax2, hSig, ...
    arrayfun(@(x) sprintf('%g', x), gabor_unik, 'UniformOutput', false), ...
    'Location', 'southoutside', 'Orientation', 'horizontal');
title(leg2, 'Signal CST (dots)');

% -------- Legend 3: line style = Bfit --------
ax3 = axes('Position', axMain.Position, 'Color', 'none', 'Visible', 'off');
hold(ax3, 'on');

hFit = gobjects(numel(ModelB_fit_unik), 1);
for i = 1:numel(ModelB_fit_unik)
    hFit(i) = plot(ax3, nan, nan, ...
        'Color', 'k', ...
        'LineStyle', lineStyles_fit{i}, ...
        'LineWidth', 1.5);
end
leg3 = legend(ax3, hFit, ...
    arrayfun(@(x) sprintf('Bfit = %d', x), ModelB_fit_unik, 'UniformOutput', false), ...
    'Location', 'eastoutside');
title(leg3, 'Predicted line');
end
