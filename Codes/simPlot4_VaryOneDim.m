% simPlot4_VaryOneDim.m
%
% Compiles simulation results and plots parameter / template / metric recovery
% as a function of individual simulation factors (one factor varied per figure panel).
%
% Pipeline:
%   1. Parse IO condition folders and load truth + fitted-model files.
%   2. Build struct array R: one row per (simulation condition x fitted B-model).
%   3. Save R to /Outputs/R_A<n>.mat.
%   4. Generate figures:
%        Fig 1 : basis-function selection rates
%        Fig 2 : template RMSE + tuning-curve recovery
%        Fig 3 : behavioral-metric recovery curves
%        Fig 4 : internal-noise parameter recovery scatter
%        Fig 5 : model-recovery confusion matrices
%
% Created by Shutian Xue

clear; clc; close all;
set(0, 'DefaultFigureVisible', 'off');
flag_whitenDV=0;
%--------------%
SX_RC1_setting;
%--------------%

%% settings
str_part = 'OOD'; % <-- change if needed
iModelA_fit = 1; %1=use data-derived template; 2=use true template
nBfit = 7;          % number of B-model variants to evaluate (each uses a different internal-noise structure)
nBins_Part4 = 3; % define bins for collapsing parameter recovery points; use 3 for main text, 5 for Supp

nameFolder_Data = sprintf('%s/Data_%s', nameFolder_server, str_part);
nameFile_R = sprintf('%s/Outputs/R_A%d.mat', nameFolder_server, iModelA_fit);
nameFolder_Data_OOD = sprintf('%s/Data_OOD_%d%d', nameFolder_Data, nORI, nSF);
nameFolder_Data_NOM_Trialwise = sprintf('%s/Data_NOM_Trialwise_%d%d', nameFolder_Data, nORI, nSF);

namesMetrics_behav = {'pYES','pC','pA'};

%% find IO folders
if ~exist(nameFile_R, 'file')

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
% R: growing struct array; one row per (simulation condition × fitted B-model) pair.
% Fields are appended inside the loop below and the array grows via R = [R; Ri].

R = struct([]);

fprintf('\n============================================\n');
fprintf('%s\n', datetime('now'))
fprintf('Compiling all fitted models\n');
fprintf('============================================\n');

for iFile = 1:nFiles
    fprintf('\n%d/%d: ', iFile, nFiles)
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
    % 2. Condition-level base record
    %    fileBase inherits parsed condition metadata from info, then is
    %    augmented with compiled statistics (template recovery, behavioral
    %    metrics, basis settings).  It is copied into Ri (one row per fitted
    %    B-model) and each Ri is appended to R.
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

    % =========================================================
    % 0. Create empty placeholders for all fields
    % =========================================================
    % true values -> template-related values -> tuning-related values ->
    % metric-related values -> parameter-related values -> criterion-related values

    % ---------- true values ----------
    fileBase.lambda_whiten = nan;
    fileBase.C_contribution = nan;

    % ---------- template-related values ----------
    fileBase.template_rmse_iter = [];
    fileBase.template_R2_iter = [];

    % ---------- Basis-setting-related values ----------
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
    fileBase.basisWidthScaleORI_mode = "";
    fileBase.basisWidthScaleORI_mode_pct = nan;
    fileBase.basisWidthScaleSF_mode = "";
    fileBase.basisWidthScaleSF_mode_pct = nan;
    fileBase.asymSF_rightLeftRatio_mode = "";
    fileBase.asymSF_rightLeftRatio_mode_pct = nan;

    % ---------- tuning-related values ----------
    fileBase.margORI_true = nan(1, nORI);
    fileBase.margSF_true = nan(1, nSF);
    fileBase.margORI_est_iter_norm = [];
    fileBase.margSF_est_iter_norm = [];

    % ---------- metric-related values ----------
    fileBase.pYES_data_iter = [];
    fileBase.pC_data_iter = [];
    fileBase.pA_data_iter = [];

    fileBase.DV_allBins_iter = [];
    fileBase.nTrials_allBins_iter = [];
    fileBase.pYES_data_curve_iter = [];
    fileBase.pC_data_curve_iter = [];
    fileBase.pA_data_curve_iter = [];

    % ---------- criterion-related values ----------
    fileBase.criterion_DV_true = nan;

    fileBase.lambda_whiten = truth.lambda_whiten;
    fileBase.C_contribution = truth.C_contribution;

    % =========================================================
    % 1. template recovery
    % =========================================================
    [template_rmse_iter, template_R2_iter] = fxn_templateRecovery_simPlot2Style_iter( ...
        truth, data_compIV, nIter);

    fileBase.template_rmse_iter = template_rmse_iter(:)';
    fileBase.template_R2_iter = template_R2_iter(:)';

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

    % Store full per-iteration normalized marginals; summary stats are computed later.
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

    % Find the most frequently selected basis family across iterations
    % For ORI
    strVec = string(data_compIV.basisFxnORI_tmpl_allIter(:));
    [uniqueVals, ~, groupIdx] = unique(strVec);
    counts = accumarray(groupIdx, 1);
    [~, idxMax] = max(counts);
    fileBase.basisFxnORI_mode     = uniqueVals(idxMax);
    fileBase.basisFxnORI_mode_pct = counts(idxMax) / numel(strVec);

    % For SF
    strVec = string(data_compIV.basisFxnSF_tmpl_allIter(:));
    [uniqueVals, ~, groupIdx] = unique(strVec);
    counts = accumarray(groupIdx, 1);
    [~, idxMax] = max(counts);
    fileBase.basisFxnSF_mode     = uniqueVals(idxMax);
    fileBase.basisFxnSF_mode_pct = counts(idxMax) / numel(strVec);

    % Modal L2 ridge regularization setting across iterations.
    strVec = string(data_compIV.ridge_tmpl_allIter(:));
    [uniqueVals, ~, groupIdx] = unique(strVec);
    counts = accumarray(groupIdx, 1);
    [~, idxMax] = max(counts);
    fileBase.ridge_mode     = uniqueVals(idxMax);
    fileBase.ridge_mode_pct = counts(idxMax) / numel(strVec);

    % Modal ORI width-scale setting across iterations.
    if isfield(data_compIV, 'basisWidthScaleORI_tmpl_allIter') && ~isempty(data_compIV.basisWidthScaleORI_tmpl_allIter)
        vals = data_compIV.basisWidthScaleORI_tmpl_allIter(:);
        if isnumeric(vals)
            vals = vals(isfinite(vals));
        end
        if ~isempty(vals)
            strVec = string(vals);
            [uniqueVals, ~, groupIdx] = unique(strVec);
            counts = accumarray(groupIdx, 1);
            [~, idxMax] = max(counts);
            fileBase.basisWidthScaleORI_mode     = uniqueVals(idxMax);
            fileBase.basisWidthScaleORI_mode_pct = counts(idxMax) / numel(strVec);
        end
    end

    % Modal SF width-scale setting across iterations.
    if isfield(data_compIV, 'basisWidthScaleSF_tmpl_allIter') && ~isempty(data_compIV.basisWidthScaleSF_tmpl_allIter)
        vals = data_compIV.basisWidthScaleSF_tmpl_allIter(:);
        if isnumeric(vals)
            vals = vals(isfinite(vals));
        end
        if ~isempty(vals)
            strVec = string(vals);
            [uniqueVals, ~, groupIdx] = unique(strVec);
            counts = accumarray(groupIdx, 1);
            [~, idxMax] = max(counts);
            fileBase.basisWidthScaleSF_mode     = uniqueVals(idxMax);
            fileBase.basisWidthScaleSF_mode_pct = counts(idxMax) / numel(strVec);
        end
    end

    % Modal SF asymmetry setting across iterations.
    if isfield(data_compIV, 'asymSF_rightLeftRatio_tmpl_allIter') && ~isempty(data_compIV.asymSF_rightLeftRatio_tmpl_allIter)
        vals = data_compIV.asymSF_rightLeftRatio_tmpl_allIter(:);
        if isnumeric(vals)
            vals = vals(isfinite(vals));
        end
        if ~isempty(vals)
            strVec = string(vals);
            [uniqueVals, ~, groupIdx] = unique(strVec);
            counts = accumarray(groupIdx, 1);
            [~, idxMax] = max(counts);
            fileBase.asymSF_rightLeftRatio_mode     = uniqueVals(idxMax);
            fileBase.asymSF_rightLeftRatio_mode_pct = counts(idxMax) / numel(strVec);
        end
    end

    % =========================================================
    % 3. Metrics-related values
    %    Simulated data metrics are identical regardless of which B-model is fitted,
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
        end %iIter

        fileBase.pYES_data_iter = pYES_data_iter(:)';
        fileBase.pC_data_iter = pC_data_iter(:)';
        fileBase.pA_data_iter = pA_data_iter(:)';
        fileBase.DV_allBins_iter = DV_allBins_iter;
        fileBase.nTrials_allBins_iter = nTrials_allBins_iter;
        fileBase.pYES_data_curve_iter = pYES_data_curve_iter;
        fileBase.pC_data_curve_iter = pC_data_curve_iter;
        fileBase.pA_data_curve_iter = pA_data_curve_iter;

        found_fit_for_data = true;
        break
    end % iModelB_fit_probe

    if ~found_fit_for_data
        fprintf('No fitNOM file found in %s for any Bfit\n', info.folder_NOM);
        continue
    end

    % ---------- criterion true ----------
    if isfield(truth, 'criterion_DV_true')
        fileBase.criterion_DV_true = truth.criterion_DV_true;
    end

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

        % Ri: one record destined for R, for this (simulation condition × B-fit model) pair.
        % It inherits all condition-level fields from fileBase, then adds fit-specific fields.
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
        Ri.pYES_rmse_allIter = []; Ri.pYES_R2_allIter = [];
        Ri.pC_rmse_allIter = []; Ri.pC_R2_allIter = [];
        Ri.pA_rmse_allIter = []; Ri.pA_R2_allIter = [];

        Ri.pYES_pred_iter = [];
        Ri.pC_pred_iter = [];
        Ri.pA_pred_iter = [];

        Ri.pYES_pred_curve_iter = [];
        Ri.pC_pred_curve_iter = [];
        Ri.pA_pred_curve_iter = [];

        Ri.nLL_test_allIter = []; % full per-iteration test nLL vector (for iteration-level win rate)

        % ---------- parameter-related values ----------
        Ri.Nmul_rmse_allIter = []; Ri.Nmul_est_allIter = [];
        Ri.Nadd_rmse_allIter = []; Ri.Nadd_est_allIter = [];
        Ri.Nshared_rmse_allIter = []; Ri.Nshared_est_allIter = [];

        % ---------- criterion-related values ----------
        Ri.criterion_DV_rmse_allIter = [];
        Ri.criterion_DV_est_allIter = [];

        % ---------- NLL ----------
        if isfield(data_fitNOM, 'nLL_test_allIter')
            Ri.nLL_test_allIter = data_fitNOM.nLL_test_allIter(:)'; % store full iter vector
        end

        pred_metrics_allIter = data_fitNOM.pred_metrics_allIter;
        nIter_fit = numel(pred_metrics_allIter);

        % ---------- Metric recovery ----------
        for iMetric = 1:numel(namesMetrics_behav)
            m = namesMetrics_behav{iMetric};
            [rmse_iter, R2_iter] = fxn_metricRecovery_iter(pred_metrics_allIter, m);

            Ri.(sprintf('%s_rmse_allIter', m)) = rmse_iter(:)';
            Ri.(sprintf('%s_R2_allIter', m)) = R2_iter(:)';
        end

        % ---------- Predicted metrics ----------
        pYES_pred_iter = nan(nIter_fit,1);
        pC_pred_iter = nan(nIter_fit,1);
        pA_pred_iter = nan(nIter_fit,1);

        nBins_curve = numel(pred_metrics_allIter{1}.metrics.IV_allBins);
        pYES_pred_curve_iter = nan(nIter_fit, nBins_curve);
        pC_pred_curve_iter = nan(nIter_fit, nBins_curve);
        pA_pred_curve_iter = nan(nIter_fit, nBins_curve);

        for iIter = 1:nIter_fit
            pYES_pred_iter(iIter) = getCI(pred_metrics_allIter{iIter}.metrics.pYES_pred_allBins, 2, 1);
            pC_pred_iter(iIter) = getCI(pred_metrics_allIter{iIter}.metrics.pC_pred_allBins, 2, 1);
            pA_pred_iter(iIter) = getCI(pred_metrics_allIter{iIter}.metrics.pA_pred_allBins, 2, 1);

            pYES_pred_curve_iter(iIter, :) = pred_metrics_allIter{iIter}.metrics.pYES_pred_allBins;
            pC_pred_curve_iter(iIter, :) = pred_metrics_allIter{iIter}.metrics.pC_pred_allBins;
            pA_pred_curve_iter(iIter, :) = pred_metrics_allIter{iIter}.metrics.pA_pred_allBins;
        end % iIter

        Ri.pYES_pred_iter = pYES_pred_iter(:)';
        Ri.pC_pred_iter = pC_pred_iter(:)';
        Ri.pA_pred_iter = pA_pred_iter(:)';

        Ri.pYES_pred_curve_iter = pYES_pred_curve_iter;
        Ri.pC_pred_curve_iter = pC_pred_curve_iter;
        Ri.pA_pred_curve_iter = pA_pred_curve_iter;

        % ---------- Parameter estimates ----------
        % params_mat: [nIter × nParams] matrix of fitted parameter values.
        % Columns 1:(end-1) are the internal-noise parameters (Nmul, Nadd, Nshared);
        % the last column is always the decision criterion.
        if iscell(data_fitNOM.params_est_allIter)
            params_mat = cell2mat(cellfun(@(x) x(:)', data_fitNOM.params_est_allIter, 'UniformOutput', false));
        else
            params_mat = data_fitNOM.params_est_allIter;
        end

        criterion_est_iter  = params_mat(:, end);   % last column = decision criterion
        criterion_rmse_iter = abs(criterion_est_iter - fileBase.criterion_DV_true);

        Ri.criterion_DV_rmse_allIter = criterion_rmse_iter(:)';
        Ri.criterion_DV_est_allIter = criterion_est_iter(:)';

        % param_names_all: all parameter names for this B-model (IN noise params + criterion).
        % param_names_IN:  all except the last (criterion), i.e. only the internal-noise parameters.
        param_names_all = namesModelBparams_short{iModelB_fit};
        param_names_IN  = param_names_all(1:end-1);

        for iParamIN = 1:numel(param_names_IN)
            pName   = param_names_IN{iParamIN};
            trueVal = info.(sprintf('%s_true', pName)); % ground-truth value for this parameter

            est_iter = params_mat(:, iParamIN);
            err_iter = abs(est_iter - trueVal);

            Ri.(sprintf('%s_rmse_allIter', pName)) = err_iter(:)';
            Ri.(sprintf('%s_est_allIter', pName)) = est_iter(:)';
        end

        R = [R; Ri];
    end % iModelB_fit

end % iFile

fprintf('\n\n%s: All files compiled. \n\n', datetime('now'))

% Post-processing
% Summarize all per-iteration fields once before plotting/saving.
R = fxn_attach_medians_from_iter(R, namesMetrics_behav);

% Save the compiled record to /Output
save(nameFile_R, 'R', 'lambda_whiten_unik', 'C_contribution_unik')

else
    % fprintf('Output file already exists, skipping compilation: %s\n', nameFile_R);
    % Load output file
    load(nameFile_R, 'R')
    fprintf('\n\n%s: Loaded compiled record R with %d entries.\n\n', datetime('now'), numel(R));
end % if ~exist(nameFile_R)

%% ---------- unique levels ----------
gaborCST_unik = unique([R.gaborCST]);
cSDT_unik = unique([R.cSDT_true]);
Nmul_unik = unique([R.Nmul_true]);
Nadd_unik = unique([R.Nadd_true]);
Nshared_unik = unique([R.Nshared_true]);
Bsim_unik = unique([R.iModelB_sim]);
Bfit_unik = unique([R.iModelB_fit]);
lambda_whiten_unik = unique([R.lambda_whiten]);
C_contribution_unik = unique([R.C_contribution]);

%% Shared plotting formats
setting = struct();
% ---------- output ----------
nameFolder_Figures_part4 = fullfile(nameFolder_Figures, sprintf('IO_%s_A%d', str_part, iModelA_fit));
if ~exist(nameFolder_Figures_part4, 'dir')
    mkdir(nameFolder_Figures_part4);
end

% ---------- figure sizes ----------
setting.figPos_wide   = [100 100 1500 650];
setting.figPos_medium = [100 100 1300 700];
setting.figPos_tall   = [100 100 1100 900];

% ---------- colors / styles ----------
setting.cmap_Cz = [0.85, 0.20, 0.15; 0.35, 0.35, 0.35; 0.18, 0.40, 0.70];   % low/mid/high Cz
setting.lineStyles_signal = {'-', '--', ':'};
setting.marker_signal = {'o', 's', '^'};
setting.lineStyles_fit = {'-', '--', ':', '-.'};

% ---------- appearance ----------
setting.bandAlpha = 0.16;
setting.lineWidth = 1.6;
setting.trueWidth = 2.8;
setting.markerSize = 6;
setting.fontSize = 10;
setting.errLineWidth = 1.1;
setting.axisLineWidth = 1.0;
setting.gridAlpha = 0.25;
setting.tickDir   = 'out';
setting.fig4_axisBufferProp = 0.5;
setting.fig4_markerSizeMin = 7;
setting.fig4_markerSizeMax = 20;
setting.fig4_markerSizeExp = 0.5;
setting.fig2_top_axisBufferProp = 0.12;

% Golden-angle palette: 24 hues spaced by 1/φ, fixed sat=0.60 val=0.72.
% Scales to any number of groups; index with palette(1:nGrp,:).
nPal = 64;
phi  = (1 + sqrt(5)) / 2;
hues = mod((0:nPal-1)' / phi, 1);
setting.palette = hsv2rgb([hues, repmat(0.60, nPal, 1), repmat(0.72, nPal, 1)]);

% ---------- collapse rules ----------
setting.scalarCollapseFcn = @mean;   % mean across rows
setting.curveCollapseFcn  = @mean;   % mean across rows after collapsing
setting.curveNGrid = 9;              % redefine bins after collapsing

% ---------- figure-specific fixed choices ----------
setting.fig12_Bsim = 1;
setting.fig12_Bfit = 1;

% ---------- field names ----------
setting.fig1_mode_fields = {'nBasisORI_mode', 'nBasisSF_mode', 'basisFxnORI_mode', 'basisFxnSF_mode', ...
    'basisWidthScaleORI_mode', 'basisWidthScaleSF_mode', 'ridge_mode', 'asymSF_rightLeftRatio_mode'};
setting.fig1_mode_labels = {'ORI basis mode', 'SF basis mode', 'ORI basis family mode', 'SF basis family mode', ...
    'ORI width-scale mode', 'SF width-scale mode', 'Ridge mode', 'SF asymmetry mode'};

setting.fig1_pct_fields = {'nBasisORI_mode_pct', 'nBasisSF_mode_pct', 'basisFxnORI_mode_pct', 'basisFxnSF_mode_pct', ...
    'basisWidthScaleORI_mode_pct', 'basisWidthScaleSF_mode_pct', 'ridge_mode_pct', 'asymSF_rightLeftRatio_mode_pct'};
setting.fig1_pct_labels = {'ORI mode selection rate', 'SF mode selection rate', 'ORI basis family mode selection rate', 'SF basis family mode selection rate', ...
    'ORI width-scale mode selection rate', 'SF width-scale mode selection rate', 'Ridge mode selection rate', 'SF asymmetry mode selection rate'};

setting.varFields = {'Nmul_true', 'Nadd_true', 'Nshared_true', 'gaborCST', 'cSDT_true'};
setting.varNames  = {'Nmul', 'Nadd', 'Nshared', 'signalCST', 'Cz'};

setting.fig1_varFields = setting.varFields;
setting.fig1_varNames  = setting.varNames;

setting.fig2_varFields = setting.varFields;
setting.fig2_varNames  = setting.varNames;

setting.metricFields = {'pYES', 'pC', 'pA'};
setting.metricCurveDataFields = {'pYES_data_curve_med', 'pC_data_curve_med', 'pA_data_curve_med'};
setting.metricCurvePredFields = {'pYES_pred_curve_med', 'pC_pred_curve_med', 'pA_pred_curve_med'};

setting.fig4_paramNames = {'Nmul', 'Nadd', 'Nshared', 'criterion_DV'};
setting.fig4_paramTitles = {'Nmul', 'Nadd', 'Nshared', 'criterion DV'};

% Expand plotting styles to match however many levels are present.
base_cmap_Cz = setting.cmap_Cz;
if numel(cSDT_unik) <= size(base_cmap_Cz, 1)
    setting.cmap_Cz = base_cmap_Cz(1:numel(cSDT_unik), :);
else
    nExtraCz = numel(cSDT_unik) - size(base_cmap_Cz, 1);
    setting.cmap_Cz = [base_cmap_Cz; setting.palette(1:nExtraCz, :)];
end

base_marker_signal = setting.marker_signal;
if numel(gaborCST_unik) > numel(base_marker_signal)
    nRep = ceil(numel(gaborCST_unik) / numel(base_marker_signal));
    setting.marker_signal = repmat(base_marker_signal, 1, nRep);
end
setting.marker_signal = setting.marker_signal(1:numel(gaborCST_unik));

base_lineStyles_fit = setting.lineStyles_fit;
if numel(Bfit_unik) > numel(base_lineStyles_fit)
    nRep = ceil(numel(Bfit_unik) / numel(base_lineStyles_fit));
    setting.lineStyles_fit = repmat(base_lineStyles_fit, 1, nRep);
end
setting.lineStyles_fit = setting.lineStyles_fit(1:numel(Bfit_unik));

%% Prepare summary fields for plotting
% Recompute all *_med/*_lb/*_ub fields from per-iteration values so each
% figure section can be run independently after this point.
R = fxn_attach_medians_from_iter(R, namesMetrics_behav);
fprintf('\n%s: Median/CI fields extracted for all metrics.\n', string(datetime('now')))

%% Figure 1: Basis selection
fprintf('\n%s: Fig 1: Basis selection rates\n', string(datetime('now')))

R_fig1 = R([R.iModelB_sim] == setting.fig12_Bsim & [R.iModelB_fit] == setting.fig12_Bfit);

% group labels: signalCST × Cz, one per entry in R_fig1
grpLabels_fig1 = arrayfun(@(r) sprintf('sig=%.3g x cSDT=%.3g', r.gaborCST, r.cSDT_true), ...
    R_fig1, 'UniformOutput', false);

grpLabels_fig1 = arrayfun(@(r) sprintf('sig=%.3g', r.gaborCST), ...
    R_fig1, 'UniformOutput', false);

% grpLabels_fig1 = arrayfun(@(r) sprintf('Nmul=%g x Nadd=%g x Nshared=%g', r.Nmul_true, r.Nadd_true, r.Nshared_true), ...
%     R_fig1, 'UniformOutput', false);

uGrp_fig1    = unique(grpLabels_fig1);
grpCmap_fig1 = setting.palette(1:numel(uGrp_fig1), :);

h = figure('Position', [100 100 1800 820]);

panel_order = [1, 5, 2, 6, 3, 7, 4, 8];
panel_xlabel = {'# ORI basis functions', '# SF basis functions', ...
    'ORI basis family', 'SF basis family', ...
    'ORI width scale', 'SF width scale', ...
    'Ridge penalty', 'asymSF_rightLeftRatio'};
panel_title = {'nBasisORI', 'nBasisSF', ...
    'basisFxnORI', 'basisFxnSF', ...
    'basisWidthScaleORI', 'basisWidthScaleSF', ...
    'ridge', 'asymSF_rightLeftRatio'};

for iPanel = 1:numel(setting.fig1_mode_fields)
    subplot(2,4,panel_order(iPanel)); hold on;
    plot_ranked_categorical({R_fig1.(setting.fig1_mode_fields{iPanel})}, panel_title{iPanel}, grpLabels_fig1, grpCmap_fig1, iPanel == 1);
    ylabel('% iterations selected');
    xlabel(panel_xlabel{iPanel});
    title(panel_title{iPanel});
    fxn_style_ax(gca, setting);
end

sgtitle('Fig 1: rank of selected basis settings across iterations', 'FontWeight', 'bold');

saveas(h, fullfile(nameFolder_Figures_part4, 'FigS1_selectedBasisRank.png'));
close(h);

%% Figure 2: Template recovery
%    Use one representative Bsim/Bfit because these are file-level

fprintf('\n%s: Fig 2: Template recovery\n', string(datetime('now')))

R_fig2 = R([R.iModelB_sim] == setting.fig12_Bsim & [R.iModelB_fit] == setting.fig12_Bfit);

nVars = numel(setting.fig2_varFields);
setting_fig2 = setting;

% Row-1 grouping follows lambda_whiten x C_contribution combinations.
grpLabels_lwcc = arrayfun(@(r) sprintf('lw=%.3g x cont=%.3g', r.lambda_whiten, r.C_contribution), ...
    R_fig2, 'UniformOutput', false);
uGrp_lwcc = unique(grpLabels_lwcc);
grpCmap_lwcc = setting.palette(1:numel(uGrp_lwcc), :);
grpIdx_lwcc = cell(numel(uGrp_lwcc), 1);
for iGrp = 1:numel(uGrp_lwcc)
    grpIdx_lwcc{iGrp} = strcmp(grpLabels_lwcc, uGrp_lwcc{iGrp});
end

fig2Cache = struct([]);
for iCol = 1:nVars
    vField  = setting.fig2_varFields{iCol};
    vName   = setting.fig2_varNames{iCol};
    vLevels = unique([R_fig2.(vField)]);

    fig2Cache(iCol).vField = vField; %#ok<AGROW>
    fig2Cache(iCol).vName = vName;
    fig2Cache(iCol).vLevels = vLevels;
    fig2Cache(iCol).rmse_all = fxn_collapse_scalar_by_x_iterCI( ...
        R_fig2, vField, 'template_rmse_iter', vLevels, setting.scalarCollapseFcn);
    fig2Cache(iCol).rmse_lwcc = cell(numel(uGrp_lwcc), 1);

    for iGrp = 1:numel(uGrp_lwcc)
        Rsub = R_fig2(grpIdx_lwcc{iGrp});
        fig2Cache(iCol).rmse_lwcc{iGrp} = fxn_collapse_scalar_by_x_iterCI( ...
            Rsub, vField, 'template_rmse_iter', vLevels, setting.scalarCollapseFcn);
    end

    [xAxisORI, yTrueORI, curvesORI] = fxn_collapse_tuning_panel(R_fig2, vField, vLevels, 'ORI');
    [xAxisSF,  yTrueSF,  curvesSF]  = fxn_collapse_tuning_panel(R_fig2, vField, vLevels, 'SF');
    fig2Cache(iCol).tuning(1).xAxis = xAxisORI;
    fig2Cache(iCol).tuning(1).yTrue = yTrueORI;
    fig2Cache(iCol).tuning(1).curves = curvesORI;
    fig2Cache(iCol).tuning(2).xAxis = xAxisSF;
    fig2Cache(iCol).tuning(2).yTrue = yTrueSF;
    fig2Cache(iCol).tuning(2).curves = curvesSF;
end

h = figure('Position', [100 100 2400 1050]);
tiledlayout(3, nVars, 'TileSpacing', 'loose', 'Padding', 'loose');
ax_rmse_last = [];

for iCol = 1:nVars
    vName = fig2Cache(iCol).vName;

    % ------ Row 1: RMSE ------
    ax_rmse = nexttile(iCol); hold on;
    ax_rmse_last = ax_rmse;

    % For Nmul/Nadd/Nshared/signalCST/Cz columns, plot lambda x contribution combos.
    if iCol <= 5
        for iGrp = 1:numel(uGrp_lwcc)
            Ssub = fig2Cache(iCol).rmse_lwcc{iGrp};
            if isempty(Ssub.x), continue; end
            fxn_plot_connected_dots_with_err(Ssub.x, Ssub.y, Ssub.lb, Ssub.ub, grpCmap_lwcc(iGrp,:), setting_fig2);
        end
    end

    % Always overlay grand mean.
    S = fig2Cache(iCol).rmse_all;
    if ~isempty(S.x)
        Pthick = setting_fig2; Pthick.lineWidth = setting.trueWidth;
        fxn_plot_connected_dots_with_err(S.x, S.y, S.lb, S.ub, [0 0 0], Pthick);
    end

    fxn_style_ax(gca, setting_fig2);
    xticks(fig2Cache(iCol).vLevels);
    xVals_al = unique(fig2Cache(iCol).vLevels(isfinite(fig2Cache(iCol).vLevels)));
    if isempty(xVals_al)
        xlim([0 1]);
    else
        xMin_al = min(xVals_al); xMax_al = max(xVals_al);
        if numel(xVals_al) >= 2
            buf_al = setting_fig2.fig2_top_axisBufferProp * (xMax_al - xMin_al);
        else
            buf_al = setting_fig2.fig2_top_axisBufferProp * max(abs(xVals_al(1)), 1);
        end
        xlim([xMin_al - buf_al, xMax_al + buf_al]);
    end
    xlabel(vName); ylabel('Template RMSE'); ylim([0 0.2]); title(vName);

    % ------ Rows 2 & 3: ORI then SF tuning ------
    for iRow = 1:2
        nexttile(iRow*nVars + iCol); hold on;

        yTrue = fig2Cache(iCol).tuning(iRow).yTrue;
        curves = fig2Cache(iCol).tuning(iRow).curves;
        if ~isempty(yTrue)
            yTrueFinite = yTrue(isfinite(yTrue));
            if ~isempty(yTrueFinite)
                yTrueScale = max(yTrueFinite);
            else
                yTrueScale = nan;
            end
            if isfinite(yTrueScale) && yTrueScale ~= 0
                yTrueNorm = yTrue / yTrueScale;
            else
                yTrueNorm = yTrue;
            end
            hTrue = plot(axis_tuning{iRow}, yTrueNorm, 'r-', 'LineWidth', setting.trueWidth);
            xIdx = 1:numel(yTrueNorm);

            nLev = numel(curves);
            grayMin = 0.55;
            grayMax = 0.00;
            if nLev >= 2
                grayVals = linspace(grayMin, grayMax, nLev)';
            else
                grayVals = grayMin;
            end

            legH = [hTrue; gobjects(nLev, 1)];
            legStr = [{'True'}; cell(nLev, 1)];
            nLeg = 1;
            for iLev = 1:nLev
                if isempty(curves(iLev).y), continue; end
                nLeg = nLeg + 1;
                cLev = repmat(grayVals(iLev), 1, 3);
                yLev = curves(iLev).y;
                yLevFinite = yLev(isfinite(yLev));
                if ~isempty(yLevFinite)
                    yLevScale = max(yLevFinite);
                else
                    yLevScale = nan;
                end
                if isfinite(yLevScale) && yLevScale ~= 0
                    yLevNorm = yLev / yLevScale;
                else
                    yLevNorm = yLev;
                end
                legH(nLeg) = plot(axis_tuning{iRow}, yLevNorm, 'LineWidth', setting.lineWidth, 'Color', cLev);
                rmse = fxn_curve_rmse(xIdx, yTrueNorm, xIdx, yLevNorm);
                legStr{nLeg} = sprintf('%g  RMSE=%.3f', curves(iLev).level, rmse);
            end
            legend(legH(1:nLeg), legStr(1:nLeg), 'Location', 'south', 'FontSize', setting_fig2.fontSize, 'Box', 'off');
        end

        fxn_style_ax(gca, setting_fig2);
        xticks(axisTicks_tuning{iRow}); xticklabels(axisTL_tuning{iRow});
        ylim([-.2, 1]);
        title(vName); ylabel('Marg. weight (divided by max)'); xlabel(namesFeature{iRow});
    end
end

% Legend for lambda x contribution combo lines (only row 1 condition lines).
if isgraphics(ax_rmse_last)
    h_gl = gobjects(numel(uGrp_lwcc), 1);
    for i_gl = 1:numel(uGrp_lwcc)
        h_gl(i_gl) = plot(ax_rmse_last, nan, nan, '-', 'Color', grpCmap_lwcc(i_gl,:), 'LineWidth', 1.5);
    end
    legend(ax_rmse_last, h_gl, cellstr(uGrp_lwcc), 'Location', 'eastoutside', 'Box', 'off', 'FontSize', 7);
end

sgtitle('Fig 2: Template RMSE + tuning recovery', 'FontWeight', 'bold', 'FontSize', setting_fig2.fontSize);
saveas(h, fullfile(nameFolder_Figures_part4, 'FigS2_TemplateRecovery.png'));
close(h);

%% Figure 3: Metric recovery
%    Simulated curve from Bfit = 1 only; predictions from all Bfit

fprintf('\n%s: Fig 3: Metric recovery\n', string(datetime('now')))

h = figure('Position', [100 100 2e3 1e3]);
tiledlayout(numel(setting.metricFields), numel(Bsim_unik), 'TileSpacing', 'loose', 'Padding', 'loose');

for iRow_metric = 1:numel(setting.metricFields)
    mName = setting.metricFields{iRow_metric};
    dataField = setting.metricCurveDataFields{iRow_metric};
    predField = setting.metricCurvePredFields{iRow_metric};

    for iCol_Bsim = 1:numel(Bsim_unik)
        simModel = Bsim_unik(iCol_Bsim);

    % simulated data are duplicated across Bfit; use Bfit = 1
    R_data = R([R.iModelB_sim] == simModel & [R.iModelB_fit] == 1);

        nexttile; hold on;
        isFirstPanel = (iCol_Bsim == 1) && (iRow_metric == 1);

        % collapsed simulated curve
        [Cdata, nDataBins] = fxn_collapse_curve_with_counts( ...
            R_data, 'DV_allBins_med', dataField, 'nTrials_allBins_med', setting.curveNGrid, setting.curveCollapseFcn);
        dot_nMin = nan;
        dot_nMax = nan;
        dot_sMin = max(4, setting.markerSize * 0.5);
        dot_sMax = max(dot_sMin + 2, setting.markerSize * 1.5);
        if ~isempty(Cdata.x)
            mSizes = setting.markerSize * ones(size(Cdata.x));
            if ~isempty(nDataBins) && numel(nDataBins) == numel(Cdata.x)
                nDot = nDataBins(:)';
                goodN = isfinite(nDot) & nDot > 0;
                if any(goodN)
                    nMin = min(nDot(goodN));
                    nMax = max(nDot(goodN));
                    dot_nMin = nMin;
                    dot_nMax = nMax;
                    if nMax > nMin
                        frac = (nDot - nMin) / (nMax - nMin);
                    else
                        frac = zeros(size(nDot));
                    end
                    frac(~goodN) = 0;
                    mMin = dot_sMin;
                    mMax = dot_sMax;
                    mSizes = mMin + (mMax - mMin) * (frac .^ setting.fig4_markerSizeExp);
                end
            end
            % Plot binned data
            for iDot = 1:numel(Cdata.x)
                plot(Cdata.x(iDot), Cdata.y(iDot), 'ko', ...
                    'MarkerSize', mSizes(iDot), ...
                    'MarkerFaceColor', 'w', 'LineWidth', setting.lineWidth);
            end
        end

        % predicted curves across Bfit — accumulate handles for per-panel legend
        nFit  = numel(Bfit_unik);
        legH   = gobjects(nFit, 1);
        legStr = cell(nFit, 1);
        nLeg   = 0;
        for iFit = 1:nFit
            fitModel = Bfit_unik(iFit);
            R_pred = R([R.iModelB_sim] == simModel & [R.iModelB_fit] == fitModel);

            %=====================%
            Cpred = fxn_collapse_curve(R_pred, 'DV_allBins_med', predField, setting.curveNGrid, setting.curveCollapseFcn);
            %=====================%

            if isempty(Cpred.x), continue; end
            color_pred = [.5, .5, .5];
            if fitModel == simModel
                color_pred = [1, 0, 0];
            end
            %=====================%
            hPred = plot(Cpred.x, Cpred.y, ...
                'Color', color_pred, ...
                'LineStyle', setting.lineStyles_fit{iFit}, ...
                'LineWidth', setting.lineWidth);
            %=====================%

            nLeg = nLeg + 1;
            legH(nLeg) = hPred;
            fitLabel = sprintf('B%d', fitModel);

            if exist('namesModelB', 'var') && ~isempty(namesModelB) && fitModel >= 1 && fitModel <= numel(namesModelB)
                fitLabel = char(string(namesModelB{fitModel}));
            end

            if ~isempty(Cdata.x)
                [r2, ~] = fxn_curve_fit_metrics(Cdata.x, Cdata.y, Cpred.x, Cpred.y, nDataBins);
                r2 = max(r2, 0);
                rmse = fxn_curve_rmse(Cdata.x, Cdata.y, Cpred.x, Cpred.y, nDataBins);

                if isFirstPanel
                    legStr{nLeg} = sprintf('%s | %.3f | %.2f', fitLabel, rmse, r2);
                else
                    legStr{nLeg} = sprintf('%.3f | %.2f', rmse, r2);
                end
            else
                if isFirstPanel
                    legStr{nLeg} = sprintf('%s | nan | nan', fitLabel);
                else
                    legStr{nLeg} = 'nan | nan';
                end
            end
        end
        legH   = legH(1:nLeg);
        legStr = legStr(1:nLeg);

        %=====================%
        fxn_style_ax(gca, setting);
        set(gca, 'YGrid', 'off');
        %=====================%
        % Print x label
        if iRow_metric == numel(setting.metricFields)
            xlabel('Collapsed DV bins');
        end
        % Print y label
        switch mName
            case 'pYES';  metricTitle = 'Detection rate';
            case 'pC';    metricTitle = 'Accuracy';
            case 'pA';    metricTitle = 'Consistency rate';
            otherwise;    metricTitle = mName;
        end
        ylabel(metricTitle);
        % Print title
        title('');

        switch mName
            case 'pYES'
                ylim([0 1]);
            case {'pC','pA'}
                ylim([0.5 1]);
        end

        if ~isempty(legH)
            %=====================%
            legend(legH, legStr, 'Location', 'best', 'FontSize', setting.fontSize-2, 'Box', 'off');
            %=====================%
        end

        if isFirstPanel && ~isempty(Cdata.x)
            axMain = get(h, 'CurrentAxes');
            if isempty(axMain) || ~isgraphics(axMain, 'axes')
                axMain = gca;
            end
            axPos = get(axMain, 'Position');
            axDot = axes('Position', axPos, 'Color', 'none', 'Visible', 'off', ...
                'HitTest', 'off', 'HandleVisibility', 'off');
            hold(axDot, 'on');

            refN = [10, 100, 1000];
            if isfinite(dot_nMin) && isfinite(dot_nMax) && dot_nMax > dot_nMin
                fracRef = (refN - dot_nMin) / (dot_nMax - dot_nMin);
                fracRef = min(max(fracRef, 0), 1);
                dotSizesRef = dot_sMin + (dot_sMax - dot_sMin) * (fracRef .^ setting.fig4_markerSizeExp);
            else
                dotSizesRef = [dot_sMin, (dot_sMin + dot_sMax) / 2, dot_sMax];
            end

            hDot1 = plot(axDot, nan, nan, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', dotSizesRef(1), 'LineWidth', setting.lineWidth);
            hDot2 = plot(axDot, nan, nan, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', dotSizesRef(2), 'LineWidth', setting.lineWidth);
            hDot3 = plot(axDot, nan, nan, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', dotSizesRef(3), 'LineWidth', setting.lineWidth);
            legDot = legend(axDot, [hDot1 hDot2 hDot3], {'10', '100', '1000'}, ...
                'Location', 'northwest', 'FontSize', setting.fontSize-2, 'Box', 'off');
            title(legDot, 'Dot size');
        end
    end % iCol_Bsim
end % iRow_metric

sgtitle('Fig 3: metric recovery', 'FontWeight', 'bold');
saveas(h, fullfile(nameFolder_Figures_part4, 'FigS3_metricRecovery.png'));
close(h);

%% Figure 4: parameter recovery
%    Only matched pairs Bsim = Bfit

fprintf('\n%s: Fig 4: Parameter recovery\n', string(datetime('now')))

setting_fig4 = setting;
setting_fig4.fontSize = setting.fontSize;

% Condition coding in Fig 4:
%   marker shape                    = cSDT_true
%   marker face shade (light→black) = signalCST low→high

% Marker shapes for cSDT levels
markers_cz_fig4 = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};
nCz_fig4  = numel(cSDT_unik);
if nCz_fig4 > numel(markers_cz_fig4)
    markers_cz_fig4 = repmat(markers_cz_fig4, 1, ceil(nCz_fig4 / numel(markers_cz_fig4)));
end
markers_cz_fig4 = markers_cz_fig4(1:nCz_fig4);

% Grayscale ramp for signalCST: index 1 = lightest (0.78), last = black (0)
nSig_fig4 = numel(gaborCST_unik);
faceShadeVals_fig4 = linspace(0.78, 0.0, max(nSig_fig4, 2));
if nSig_fig4 == 1, faceShadeVals_fig4 = 0.0; end
% [sorted low→high so index matches sorted gaborCST_unik]

h = figure('Position', [100 100 1300 900]);
tiledlayout(numel(setting.fig4_paramNames), numel(Bfit_unik), 'TileSpacing', 'compact', 'Padding', 'compact');
ax_fig4_last = gobjects(1);

for iRow_param = 1:numel(setting.fig4_paramNames)
    pName = setting.fig4_paramNames{iRow_param};
    pTitle = setting.fig4_paramTitles{iRow_param};

    for iCol_Bsimfit = 1:numel(Bfit_unik)
        Bsimfit = Bfit_unik(iCol_Bsimfit);

        nexttile; hold on;

        % Select datasets with matched simulating and fitted model
        Rsub = R([R.iModelB_sim] == Bsimfit & [R.iModelB_fit] == Bsimfit);

        % skip parameters absent in reduced models
        if strcmp(pName, 'Nshared') && ~any(strcmp(namesModelBparams_short{Bsimfit}, 'Nshared'))
            axis off; continue
        end
        if strcmp(pName, 'Nmul') && ~any(strcmp(namesModelBparams_short{Bsimfit}, 'Nmul'))
            axis off; continue
        end
        if strcmp(pName, 'Nadd') && ~any(strcmp(namesModelBparams_short{Bsimfit}, 'Nadd'))
            axis off; continue
        end

        % one series per (signalCST x cSDT_true), pooling across lambda_whiten
        for iSig = 1:numel(gaborCST_unik)
            for iCz = 1:numel(cSDT_unik)
                S = fxn_collapse_param_byLevel( ...
                    Rsub, pName, ...
                    setting.scalarCollapseFcn, nBins_Part4, gaborCST_unik(iSig), cSDT_unik(iCz));

                if isempty(S.x), continue; end

                % marker shape = cSDT, face shade = signalCST
                mkrCz   = markers_cz_fig4{iCz};
                faceClr = repmat(faceShadeVals_fig4(iSig), 1, 3);

                good = isfinite(S.x) & isfinite(S.y);
                if any(good)
                    plot(S.x(good), S.y(good), '-', ...
                        'Color', faceClr, ...
                        'LineWidth', 1);
                end

                for iPt = 1:numel(S.x)
                    if ~isfinite(S.x(iPt)) || ~isfinite(S.y(iPt))
                        continue
                    end

                    Pplot = setting;
                    Pplot.fontSize = setting_fig4.fontSize;
                    if strcmp(pName, 'criterion_DV')
                        nCt = max(S.n(iPt), 1);
                        sc  = nCt ^ setting.fig4_markerSizeExp;
                        szMin = setting.fig4_markerSizeMin * 0.5;
                        szMax = setting.fig4_markerSizeMax * 0.5;
                        Pplot.markerSize = szMin + (szMax - szMin) * (sc - 1) / max(sc, 1);
                        Pplot.markerSize = min(max(Pplot.markerSize, szMin), szMax);
                    end

                    %=====================%
                    line([S.x(iPt) S.x(iPt)], [S.lb(iPt) S.ub(iPt)], ...
                        'Color', faceClr, 'LineWidth', Pplot.errLineWidth);
                    scatter(S.x(iPt), S.y(iPt), Pplot.markerSize^2, faceClr, mkrCz, ...
                        'MarkerEdgeColor', 'k', ...
                        'MarkerFaceAlpha', 0.5, ...
                        'LineWidth', Pplot.lineWidth);
                    %=====================%
                end
            end
        end

        %=====================%
        if strcmp(pName, 'criterion_DV')
            S_axi = fxn_bin_criterion_values(Rsub, nBins_Part4, @mean);
            tickVals = S_axi.x(:)';
        else
            tickVals = unique([Rsub.(sprintf('%s_true', pName))]);
        end
        tickVals = unique(tickVals(isfinite(tickVals)));
        if isempty(tickVals)
            mn = 0; mx = 1;
        else
            mn = min(tickVals); mx = max(tickVals);
            if numel(tickVals) >= 2
                buf_axi = setting.fig4_axisBufferProp * (mx - mn);
            else
                buf_axi = setting.fig4_axisBufferProp * max(abs(tickVals(1)), 1);
            end
            mn = mn - buf_axi; mx = mx + buf_axi;
        end
        %=====================%
        plot([mn mx], [mn mx], 'k--', 'LineWidth', 1);
        xlim([mn mx]); ylim([mn mx]);
        if ~isempty(tickVals)
            xticks(tickVals);
            yticks(tickVals);
        end

        %=====================%
        axis square; fxn_style_ax(gca, setting_fig4);
        %=====================%
        set(gca, 'XGrid', 'off', 'YGrid', 'off');
        if iRow_param == numel(setting.fig4_paramNames)
            xlabel('True');
        end
        if iCol_Bsimfit == 1
            ylabel(pTitle);
        end
        Ppanel = fxn_collect_fig4_panel_points( ...
            Rsub, pName, gaborCST_unik, cSDT_unik, setting.scalarCollapseFcn, nBins_Part4);
        [rho_panel, ~, ~, rmse_panel] = fxn_param_recovery_stats(Ppanel.x, Ppanel.y);
        title(sprintf('Bfit=Bsim=%d \nr=%.2f | RMSE=%.3f', Bsimfit, rho_panel, rmse_panel));
        ax_fig4_last = gca;
    end % iCol_Bsimfit
end % iRow_param

sgtitle('Fig 4: Parameter recovery | Matched model pairs', 'FontWeight', 'bold');

if isgraphics(ax_fig4_last)
    pos_fig4leg = ax_fig4_last.Position;
    legFontSz = max(ax_fig4_last.FontSize - 1, 7);

    % Legend 1: marker shape = cSDT_true
    axShape_fig4 = axes('Position', pos_fig4leg, 'Color', 'none', 'Visible', 'off', ...
        'HitTest', 'off', 'HandleVisibility', 'off');
    hold(axShape_fig4, 'on');
    hShape_fig4 = gobjects(nCz_fig4, 1);
    for i_leg = 1:nCz_fig4
        hShape_fig4(i_leg) = plot(axShape_fig4, nan, nan, ...
            'LineStyle', 'none', 'Marker', markers_cz_fig4{i_leg}, ...
            'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k', ...
            'MarkerSize', 7, 'LineWidth', 1.5);
    end
    leg1_fig4 = legend(axShape_fig4, hShape_fig4, ...
        arrayfun(@(x) sprintf('%g', x), cSDT_unik, 'UniformOutput', false), ...
        'Location', 'southeast', 'Orientation', 'vertical', 'Box', 'on', ...
        'FontSize', legFontSz);
    title(leg1_fig4, 'cSDT (shape)', 'Interpreter', 'none');

    % Legend 2: face shade = signalCST (light grey → black)
    axFace_fig4 = axes('Position', pos_fig4leg, 'Color', 'none', 'Visible', 'off', ...
        'HitTest', 'off', 'HandleVisibility', 'off');
    hold(axFace_fig4, 'on');
    hFace_fig4 = gobjects(nSig_fig4, 1);
    for i_leg = 1:nSig_fig4
        fc = repmat(faceShadeVals_fig4(i_leg), 1, 3);
        hFace_fig4(i_leg) = plot(axFace_fig4, nan, nan, ...
            'LineStyle', 'none', 'Marker', 'o', ...
            'MarkerFaceColor', fc, 'MarkerEdgeColor', 'k', ...
            'MarkerSize', 7, 'LineWidth', 1.5);
    end
    leg2_fig4 = legend(axFace_fig4, hFace_fig4, ...
        arrayfun(@(x) sprintf('%g', x), gaborCST_unik, 'UniformOutput', false), ...
        'Location', 'northwest', 'Orientation', 'vertical', 'Box', 'on', ...
        'FontSize', legFontSz);
    title(leg2_fig4, 'signalCST (face shade)', 'Interpreter', 'none');
end
saveas(h, fullfile(nameFolder_Figures_part4, 'FigS4_paramRecovery.png'));
close(h);

%% Figure 5: model recovery
fprintf('\n%s: Fig 5: Model recovery\n', string(datetime('now')))

% REPLACE the block from "W_bestNLL = nan(..." down to the closing "end" of that if-block

% Win rate: for each condition (nameIO), compute the fraction of iterations on which
% each fitted model has the lowest test nLL.  Then average those fractions across
% conditions (equal weight per condition).
W_bestNLL = nan(numel(Bfit_unik), numel(Bsim_unik));
if ~isempty(R) && isfield(R, 'nLL_test_allIter')
    if isfield(R, 'nameIO')
        datasetIDs_wr = string({R.nameIO});
    else
        datasetIDs_wr = string(arrayfun(@(r) sprintf('Nm%g_Na%g_Ns%g_G%g_cSDT%g', ...
            r.Nmul_true, r.Nadd_true, r.Nshared_true, r.gaborCST, r.cSDT_true), R, 'UniformOutput', false));
    end
    for iSim_wr = 1:numel(Bsim_unik)
        simModel_wr = Bsim_unik(iSim_wr);
        idxSim_wr   = [R.iModelB_sim] == simModel_wr;
        idsSim_wr   = unique(datasetIDs_wr(idxSim_wr));

        % winRatePerCond_wr: [nBfit x nConds] — per-condition per-model win rates
        winRatePerCond_wr = nan(numel(Bfit_unik), numel(idsSim_wr));

        for iID_wr = 1:numel(idsSim_wr)
            idxID_wr = idxSim_wr & datasetIDs_wr == idsSim_wr(iID_wr);
            Rsub_wr  = R(idxID_wr);
            if isempty(Rsub_wr), continue; end

            % Build [nBfit x nIter] matrix of test nLL values.
            nIter_wr = max(cellfun(@numel, {Rsub_wr.nLL_test_allIter}));
            if nIter_wr == 0, continue; end
            nllMat_wr = nan(numel(Bfit_unik), nIter_wr);
            for iFit_wr = 1:numel(Bfit_unik)
                idxFit_wr = find([Rsub_wr.iModelB_fit] == Bfit_unik(iFit_wr), 1);
                if isempty(idxFit_wr), continue; end
                v = Rsub_wr(idxFit_wr).nLL_test_allIter(:)';
                nllMat_wr(iFit_wr, 1:numel(v)) = v;
            end

            % Per-iteration winner vote.
            wins_iter_wr  = zeros(numel(Bfit_unik), 1);
            nValidIter_wr = 0;
            for iIt = 1:nIter_wr
                col = nllMat_wr(:, iIt);
                if ~any(isfinite(col)), continue; end
                [~, iBest] = min(col);
                wins_iter_wr(iBest) = wins_iter_wr(iBest) + 1;
                nValidIter_wr = nValidIter_wr + 1;
            end
            if nValidIter_wr > 0
                winRatePerCond_wr(:, iID_wr) = wins_iter_wr / nValidIter_wr;
            end
        end % iID_wr

        % Average win rate equally across conditions.
        validCols_wr = any(isfinite(winRatePerCond_wr), 1);
        if any(validCols_wr)
            W_bestNLL(:, iSim_wr) = mean(winRatePerCond_wr(:, validCols_wr), 2, 'omitnan');
        end
    end % iSim_wr
end

M_curveRMSE = nan(numel(Bfit_unik), numel(Bsim_unik));
M_curveR2   = nan(numel(Bfit_unik), numel(Bsim_unik));
M_paramRMSE = nan(numel(Bfit_unik), numel(Bsim_unik));
M_paramR2   = nan(numel(Bfit_unik), numel(Bsim_unik));

for iSim = 1:numel(Bsim_unik)
    simModel = Bsim_unik(iSim);
    for iFit = 1:numel(Bfit_unik)
        fitModel = Bfit_unik(iFit);
        Rsub = R([R.iModelB_sim] == simModel & [R.iModelB_fit] == fitModel);
        if isempty(Rsub), continue; end

        % metric RMSE
        rmseMetricVals = [[Rsub.pYES_rmse]'; [Rsub.pC_rmse]'; [Rsub.pA_rmse]'];
        rmseMetricVals = rmseMetricVals(isfinite(rmseMetricVals));
        if ~isempty(rmseMetricVals)
            M_curveRMSE(iFit, iSim) = mean(rmseMetricVals, 'omitnan');
        end

        % metric R2
        r2MetricVals = [[Rsub.pYES_R2]'; [Rsub.pC_R2]'; [Rsub.pA_R2]'];
        r2MetricVals = r2MetricVals(isfinite(r2MetricVals));
        if ~isempty(r2MetricVals)
            M_curveR2(iFit, iSim) = mean(r2MetricVals, 'omitnan');
        end

        % parameter RMSE
        rmseParamFields = strcat(namesModelBparams_short{fitModel}, '_rmse');
        rmseParamVals = [];
        for iFld = 1:numel(rmseParamFields)
            if isfield(Rsub, rmseParamFields{iFld})
                vals = [Rsub.(rmseParamFields{iFld})]';
                rmseParamVals = [rmseParamVals; vals(isfinite(vals))]; %#ok<AGROW>
            end
        end
        if ~isempty(rmseParamVals)
            M_paramRMSE(iFit, iSim) = mean(rmseParamVals, 'omitnan');
        end

        % parameter R2: R² of true vs estimated values across conditions
        pNames_fit = namesModelBparams_short{fitModel};
        r2ParamVals = [];
        for iP = 1:numel(pNames_fit)
            pN = pNames_fit{iP};
            if strcmp(pN, 'criterion_DV')
                trueField = 'criterion_DV_true';
                estField  = 'criterion_DV_est_med';
            else
                trueField = sprintf('%s_true', pN);
                estField  = sprintf('%s_est_med', pN);
            end
            if ~isfield(Rsub, trueField) || ~isfield(Rsub, estField), continue; end
            yTrue = [Rsub.(trueField)]';
            yEst  = [Rsub.(estField)]';
            good  = isfinite(yTrue) & isfinite(yEst);
            if sum(good) < 2, continue; end
            yTrue = yTrue(good); yEst = yEst(good);
            sse = sum((yTrue - yEst).^2);
            sst = sum((yTrue - mean(yTrue)).^2);
            if sst > 0
                r2ParamVals(end+1) = 1 - sse/sst; %#ok<AGROW>
            end
        end
        if ~isempty(r2ParamVals)
            M_paramR2(iFit, iSim) = mean(r2ParamVals, 'omitnan');
        end
    end
end

% ---------- delta NLL matrix (median ΔnLL from best fit, aggregated across all conditions) ----------
M_deltaNLL = nan(numel(Bfit_unik), numel(Bsim_unik));

if ~isempty(R) && isfield(R, 'nLL_med')
    if isfield(R, 'nameIO')
        datasetIDs_dn = string({R.nameIO});
    else
        datasetIDs_dn = string(arrayfun(@(r) sprintf('Nm%g_Na%g_Ns%g_G%g_cSDT%g', ...
            r.Nmul_true, r.Nadd_true, r.Nshared_true, r.gaborCST, r.cSDT_true), R, 'UniformOutput', false));
    end

    for iSim_dn = 1:numel(Bsim_unik)
        simModel_dn = Bsim_unik(iSim_dn);
        idxSim_dn   = [R.iModelB_sim] == simModel_dn;
        idsSim_dn   = unique(datasetIDs_dn(idxSim_dn));

        deltaVals_dn = cell(numel(Bfit_unik), 1);

        for iID_dn = 1:numel(idsSim_dn)
            idxID_dn = idxSim_dn & datasetIDs_dn == idsSim_dn(iID_dn);
            Rkey_dn  = R(idxID_dn);
            if isempty(Rkey_dn), continue; end

            scores_dn = nan(numel(Bfit_unik), 1);
            for iFit_dn = 1:numel(Bfit_unik)
                idxFit_dn = [Rkey_dn.iModelB_fit] == Bfit_unik(iFit_dn);
                vals_dn   = [Rkey_dn(idxFit_dn).nLL_med]';
                vals_dn   = vals_dn(isfinite(vals_dn));
                if ~isempty(vals_dn)
                    scores_dn(iFit_dn) = mean(vals_dn, 'omitnan');
                end
            end

            if ~any(isfinite(scores_dn)), continue; end
            bestScore_dn = min(scores_dn);

            for iFit_dn = 1:numel(Bfit_unik)
                if isfinite(scores_dn(iFit_dn))
                    deltaVals_dn{iFit_dn}(end+1, 1) = scores_dn(iFit_dn) - bestScore_dn; %#ok<AGROW>
                end
            end
        end

        for iFit_dn = 1:numel(Bfit_unik)
            if ~isempty(deltaVals_dn{iFit_dn})
                M_deltaNLL(iFit_dn, iSim_dn) = median(deltaVals_dn{iFit_dn}, 'omitnan');
            end
        end
    end
end

% ---------- 2x3 figure ----------
% Layout:
%   col 1: NLL metrics      | row 1: win rate,     row 2: median delta NLL
%   col 2: metric recovery  | row 1: RMSE,          row 2: R2
%   col 2: metric/param RMSE | row 1: Metric RMSE,   row 2: Param RMSE

% Layout: 2x2
%   (1,1) Win rate      (1,2) Metric RMSE
%   (2,1) Median ΔnLL   (2,2) Param RMSE
panels_fig5 = { ...
    W_bestNLL,   'Win rate (% iters, avg over conds)',  [0 1]; ...
    M_deltaNLL,  'ΔNLL from the best model', [ 0 160]; ...
    M_curveRMSE, 'Metrics RMSE',            [0.05, 0.07]; ...
    M_paramRMSE, 'Parameters RMSE',             [1 2.8] };

% tile order in a 2x2 tiledlayout (row-major): (1,1),(1,2),(2,1),(2,2)
% desired: WinRate, MetricRMSE, DeltaNLL, ParamRMSE
tileOrder = [1, 3, 2, 4]; % index into panels_fig5

h = figure('Position', [100 100 900 700]);
tlo = tiledlayout(2, 2, 'TileSpacing', 'loose', 'Padding', 'loose');

% Use model names for axis ticks when available.
xTickLabels_fig5 = arrayfun(@(k) sprintf('B%d', k), Bsim_unik, 'UniformOutput', false);
yTickLabels_fig5 = arrayfun(@(k) sprintf('B%d', k), Bfit_unik, 'UniformOutput', false);
if exist('namesModelB', 'var') && ~isempty(namesModelB)
    for iLab = 1:numel(Bsim_unik)
        k = Bsim_unik(iLab);
        if k >= 1 && k <= numel(namesModelB)
            xTickLabels_fig5{iLab} = char(string(namesModelB{k}));
        end
    end
    for iLab = 1:numel(Bfit_unik)
        k = Bfit_unik(iLab);
        if k >= 1 && k <= numel(namesModelB)
            yTickLabels_fig5{iLab} = char(string(namesModelB{k}));
        end
    end
end

for iTile = 1:4
    iPan = tileOrder(iTile);
    M_plot   = panels_fig5{iPan, 1};
    ttl_plot = panels_fig5{iPan, 2};
    clim_fix = panels_fig5{iPan, 3};

    nexttile; hold on;
    imagesc(1:numel(Bsim_unik), 1:numel(Bfit_unik), M_plot);
    set(gca, 'YDir', 'normal', ...
        'XTick', 1:numel(Bsim_unik), 'XTickLabel', xTickLabels_fig5, ...
        'YTick', 1:numel(Bfit_unik), 'YTickLabel', yTickLabels_fig5, ...
        'XLim', [0.5, numel(Bsim_unik)+0.5], ...
        'YLim', [0.5, numel(Bfit_unik)+0.5]);
    xlabel('Simulating model'); ylabel('Fitting model');
    title(ttl_plot);
    fxn_style_ax(gca, setting);
    axis square;

    % Keep tile sizing under TiledChartLayout control to avoid Position warnings.

    % color axis — use each panel's own data range for full color spectrum
    finVals = M_plot(isfinite(M_plot));
    if ~isempty(finVals)
        cLo = min(finVals); cHi = max(finVals);
        if cLo == cHi, cLo = cLo - eps; cHi = cHi + eps; end
    else
        cLo = 0; cHi = 1;
    end
    clim([cLo cHi]);
    % show colorbar
    colorbar('Location', 'eastoutside');

    % annotate values with adaptive font color (low → white, high → black)
    for iFit_ann = 1:numel(Bfit_unik)
        for iSim_ann = 1:numel(Bsim_unik)
            v = M_plot(iFit_ann, iSim_ann);
            if ~isfinite(v), continue; end
            frac = (v - cLo) / (cHi - cLo);
            frac = max(0, min(1, frac));
            if frac < 0.5
                txtClr = [1 1 1]; % white for low values
            else
                txtClr = [0 0 0]; % black for high values
            end
            text(iSim_ann, iFit_ann, sprintf('%.3f', v), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', 8, 'Color', txtClr);
        end
    end

    % highlight diagonal cells (matched sim/fit model) with thick black borders
    nDiag = min(numel(Bfit_unik), numel(Bsim_unik));
    for iDiag = 1:nDiag
        rectangle('Position', [iDiag - 0.5, iDiag - 0.5, 1, 1], ...
            'EdgeColor', 'k', 'LineWidth', 2.5, 'LineStyle', '-');
    end

end % for iTile

sgtitle('Fig 5: Model recovery', 'FontWeight', 'bold');
saveas(h, fullfile(nameFolder_Figures_part4, 'FigS5_modelRecovery.png'));
close(h);

%% Figure 6: parameter-pair correlations
fprintf('\n%s: Fig 6: Parameter-pair correlations\n', string(datetime('now')))

for iFit = 1:numel(Bfit_unik)
    fitModel = Bfit_unik(iFit);

    % Use matched pairs to keep one coherent parameterization per model.
    Rsub = R([R.iModelB_sim] == fitModel & [R.iModelB_fit] == fitModel);
    if isempty(Rsub), continue; end

    pNames_all = namesModelBparams_short{fitModel};
    % Use per-iteration estimate arrays as source for dots
    iterFields = cell(size(pNames_all));
    for iP = 1:numel(pNames_all)
        if strcmp(pNames_all{iP}, 'criterion_DV')
            iterFields{iP} = 'criterion_DV_est_allIter';
        else
            iterFields{iP} = sprintf('%s_est_allIter', pNames_all{iP});
        end
    end

    keepParam = false(size(pNames_all));
    for iP = 1:numel(pNames_all)
        if ~isfield(Rsub, iterFields{iP}), continue; end
        vals = vertcat(Rsub.(iterFields{iP}));
        keepParam(iP) = any(isfinite(vals(:)));
    end

    pNames = pNames_all(keepParam);
    iterFields = iterFields(keepParam);
    nParam = numel(pNames);
    if nParam < 2, continue; end

    sigVals_all = [Rsub.gaborCST]';
    czVals_all  = [Rsub.cSDT_true]';
    [isSig, sigIdx_all] = ismember(sigVals_all, gaborCST_unik);
    [isCz,  czIdx_all]  = ismember(czVals_all,  cSDT_unik);
    validGroupBase = isfinite(sigVals_all) & isfinite(czVals_all) & isSig & isCz;

    groupRows = cell(numel(gaborCST_unik), numel(cSDT_unik));
    for iSig = 1:numel(gaborCST_unik)
        for iCz = 1:numel(cSDT_unik)
            groupRows{iSig, iCz} = find(validGroupBase & sigIdx_all == iSig & czIdx_all == iCz);
        end
    end

    % Grayscale color per condition based on sumsqr(Nmul_true, Nadd_true, Nshared_true).
    % Black = highest sumsqr, light gray = lowest.
    inFields_gray = {'Nmul_true', 'Nadd_true', 'Nshared_true'};
    sumsqr_cond = zeros(numel(Rsub), 1);
    for iGF = 1:numel(inFields_gray)
        if isfield(Rsub, inFields_gray{iGF})
            v = [Rsub.(inFields_gray{iGF})]';
            v(~isfinite(v)) = 0;
            sumsqr_cond = sumsqr_cond + v .^ 2;
        end
    end
    sq_min = min(sumsqr_cond);
    sq_max = max(sumsqr_cond);
    if sq_max > sq_min
        sq_norm = (sumsqr_cond - sq_min) / (sq_max - sq_min);
    else
        sq_norm = zeros(size(sumsqr_cond));
    end
    lightGrayVal = 0.78;
    % normVal=1 → black (0), normVal=0 → light gray (lightGrayVal)
    condGrayColor = (1 - sq_norm) * lightGrayVal .* ones(numel(Rsub), 3);

    pairIdx = nchoosek(1:nParam, 2);
    nPairs  = size(pairIdx, 1);
    nCols   = min(3, nPairs);
    nRows   = ceil(nPairs / nCols);

    h = figure('Position', [100 100 420*nCols 320*nRows]);
    tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for iPair = 1:nPairs
        iP1 = pairIdx(iPair, 1);
        iP2 = pairIdx(iPair, 2);

        nexttile; hold on;

        % Collect all iter values pooled across conditions for overall r
        xPool = [];
        yPool = [];

        for iSig = 1:numel(gaborCST_unik)
            for iCz = 1:numel(cSDT_unik)
                idxRows = groupRows{iSig, iCz};
                if isempty(idxRows), continue; end

                mkr = setting.marker_signal{iSig};

                for iRow = 1:numel(idxRows)
                    rIdx = idxRows(iRow);
                    col  = condGrayColor(rIdx, :);
                    xIter = Rsub(rIdx).(iterFields{iP1})(:);
                    yIter = Rsub(rIdx).(iterFields{iP2})(:);
                    goodIt = isfinite(xIter) & isfinite(yIter);
                    if sum(goodIt) < 1, continue; end

                    xi = xIter(goodIt);
                    yi = yIter(goodIt);
                    xPool = [xPool; xi]; %#ok<AGROW>
                    yPool = [yPool; yi]; %#ok<AGROW>

                    % dots: one per iteration
                    plot(xi, yi, mkr, ...
                        'LineStyle', 'none', ...
                        'MarkerSize', max(setting.markerSize - 2, 4), ...
                        'MarkerFaceColor', col, ...
                        'MarkerEdgeColor', 'w', ...
                        'LineWidth', 0.5);

                    % regression line per condition
                    if sum(goodIt) >= 2 && numel(unique(xi)) >= 2
                        pfit_cond = polyfit(xi, yi, 1);
                        xfit_cond = linspace(min(xi), max(xi), 60);
                        plot(xfit_cond, polyval(pfit_cond, xfit_cond), '-', ...
                            'Color', col, 'LineWidth', 1.4);
                    end
                end
            end
        end

        if numel(xPool) >= 2
            [r_val, p_val] = corr(xPool, yPool);
            title(sprintf('%s vs %s | r=%.2f, p=%.3g', pNames{iP1}, pNames{iP2}, r_val, p_val));
        else
            title(sprintf('%s vs %s', pNames{iP1}, pNames{iP2}));
        end

        axis square;
        xlabel(sprintf('Est. %s', pNames{iP1}));
        ylabel(sprintf('Est. %s', pNames{iP2}));
        fxn_style_ax(gca, setting);
    end

    % Figure-level legends: grayscale = sumsqr(IN params), marker = signalCST
    axBase = gca;
    pos = axBase.Position;

    % Build a few representative grayscale legend entries
    axColor = axes('Position', pos, 'Color', 'none', 'Visible', 'off', ...
        'HitTest', 'off', 'HandleVisibility', 'off');
    hold(axColor, 'on');
    sq_unik_sorted = unique(sumsqr_cond);
    nGrayLeg = min(numel(sq_unik_sorted), 5);
    legIdxGray = round(linspace(1, numel(sq_unik_sorted), nGrayLeg));
    hGray = gobjects(nGrayLeg, 1);
    grayLegLabels = cell(nGrayLeg, 1);
    for iGL = 1:nGrayLeg
        sqVal = sq_unik_sorted(legIdxGray(iGL));
        normGL = (sq_max > sq_min) * (sqVal - sq_min) / max(sq_max - sq_min, eps);
        grayGL = (1 - normGL) * lightGrayVal;
        hGray(iGL) = plot(axColor, nan, nan, 'o', ...
            'LineStyle', 'none', ...
            'MarkerFaceColor', grayGL * [1 1 1], ...
            'MarkerEdgeColor', grayGL * [1 1 1], ...
            'MarkerSize', 7);
        grayLegLabels{iGL} = sprintf('%.3g', sqVal);
    end
    leg1 = legend(axColor, hGray, grayLegLabels, ...
        'Location', 'northeast', 'Box', 'on', 'FontSize', max(setting.fontSize - 1, 8));
    title(leg1, 'Total internal var.');

    axSig = axes('Position', pos, 'Color', 'none', 'Visible', 'off', ...
        'HitTest', 'off', 'HandleVisibility', 'off');
    hold(axSig, 'on');
    hSig = gobjects(numel(gaborCST_unik), 1);
    for iSig = 1:numel(gaborCST_unik)
        hSig(iSig) = plot(axSig, nan, nan, ...
            'LineStyle', 'none', ...
            'Marker', setting.marker_signal{iSig}, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 7);
    end
    leg2 = legend(axSig, hSig, ...
        arrayfun(@(x) sprintf('%g', x), gaborCST_unik, 'UniformOutput', false), ...
        'Location', 'southeast', 'Box', 'on', 'FontSize', max(setting.fontSize - 1, 8));
    title(leg2, 'signalCST');

    sgtitle(sprintf('Fig 6: parameter-pair correlations | Bfit=Bsim=%d', fitModel), 'FontWeight', 'bold');
    saveas(h, fullfile(nameFolder_Figures_part4, sprintf('FigS6_paramCorr_B%d.png', fitModel)));
    close(h);
end % iModelB_fit

%% Figure 6B: parameter-pair correlation summary by signalCST x cSDT
fprintf('\n%s: Fig 6B: Parameter-pair correlation summary\n', string(datetime('now')))

% Fixed 6 parameter-pair columns from the canonical set used in Fig 6.
pairParams_fig6b = {'Nmul', 'Nadd', 'Nshared', 'criterion_DV'};
pairLabels_fig6b = {'N_{mul}', 'N_{add}', 'N_{shared}', 'Criterion DV'};
pairIdx_fig6b = nchoosek(1:numel(pairParams_fig6b), 2);
pairIdx_fig6b = pairIdx_fig6b([1, 2, 4, 3, 5, 6], :); % rearrange so all xx-cDV pairs are last
nCols_fig6b = size(pairIdx_fig6b, 1); % expected = 6

% One row per matched model (Bfit = Bsim), expected 7 rows.
Bdiag_fig6b = intersect(Bfit_unik(:)', Bsim_unik(:)');
nRows_fig6b = numel(Bdiag_fig6b);

if nRows_fig6b == 0
    fprintf('Skipping Fig 6B: no matched Bfit=Bsim models found.\n');
else
    nSig_fig6b = numel(gaborCST_unik);
    nCz_fig6b = numel(cSDT_unik);
    nBars_fig6b = nSig_fig6b * nCz_fig6b;
    xBar_fig6b = 1:nBars_fig6b;

    barLabels_fig6b = cell(1, nBars_fig6b);
    barSigIdx_fig6b = nan(1, nBars_fig6b);
    barCzIdx_fig6b = nan(1, nBars_fig6b);
    iBar_fig6b = 0;
    for iCz = 1:nCz_fig6b
        for iSig = 1:nSig_fig6b
            iBar_fig6b = iBar_fig6b + 1;
            barSigIdx_fig6b(iBar_fig6b) = iSig;
            barCzIdx_fig6b(iBar_fig6b) = iCz;
            barLabels_fig6b{iBar_fig6b} = sprintf('G%g|C%g', gaborCST_unik(iSig), cSDT_unik(iCz));
        end
    end

    h = figure('Position', [80 80 330*nCols_fig6b 210*nRows_fig6b]);
    tiledlayout(nRows_fig6b, nCols_fig6b, 'TileSpacing', 'compact', 'Padding', 'compact');

    for iRow = 1:nRows_fig6b
        fitModel = Bdiag_fig6b(iRow);
        Rsub = R([R.iModelB_sim] == fitModel & [R.iModelB_fit] == fitModel);

        for iCol = 1:nCols_fig6b
            nexttile((iRow-1)*nCols_fig6b + iCol); hold on;

            iP1 = pairIdx_fig6b(iCol, 1);
            iP2 = pairIdx_fig6b(iCol, 2);
            p1 = pairParams_fig6b{iP1};
            p2 = pairParams_fig6b{iP2};

            if strcmp(p1, 'criterion_DV')
                f1 = 'criterion_DV_est_allIter';
            else
                f1 = sprintf('%s_est_allIter', p1);
            end
            if strcmp(p2, 'criterion_DV')
                f2 = 'criterion_DV_est_allIter';
            else
                f2 = sprintf('%s_est_allIter', p2);
            end

            yMean = nan(nBars_fig6b, 1);
            ySem  = nan(nBars_fig6b, 1);
            nRep  = zeros(nBars_fig6b, 1);

            for iB = 1:nBars_fig6b
                iSig = barSigIdx_fig6b(iB);
                iCz = barCzIdx_fig6b(iB);

                idxRows = find([Rsub.gaborCST] == gaborCST_unik(iSig) & [Rsub.cSDT_true] == cSDT_unik(iCz));
                if isempty(idxRows)
                    continue;
                end

                rVals = nan(numel(idxRows), 1);
                for iRR = 1:numel(idxRows)
                    rIdx = idxRows(iRR);

                    xTmp = Rsub(rIdx).(f1);
                    yTmp = Rsub(rIdx).(f2);
                    xIter = xTmp(:);
                    yIter = yTmp(:);

                    nItCommon = min(numel(xIter), numel(yIter));
                    if nItCommon < 2
                        continue;
                    end
                    xIter = xIter(1:nItCommon);
                    yIter = yIter(1:nItCommon);

                    good = isfinite(xIter) & isfinite(yIter);
                    if sum(good) < 2
                        continue;
                    end

                    xi = xIter(good);
                    yi = yIter(good);
                    if numel(unique(xi)) < 2 || numel(unique(yi)) < 2
                        continue;
                    end

                    rVals(iRR) = corr(xi, yi);
                end

                goodR = isfinite(rVals);
                if any(goodR)
                    yMean(iB) = mean(rVals(goodR));
                    nRep(iB) = sum(goodR);
                    if nRep(iB) >= 2
                        ySem(iB) = std(rVals(goodR), 0) / sqrt(nRep(iB));
                    else
                        ySem(iB) = 0;
                    end
                end
            end

            if all(~isfinite(yMean))
                text(0.5, 0.5, 'N/A', 'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                axis off;
                continue;
            end

            b = bar(xBar_fig6b, yMean, 0.85, 'FaceColor', 'flat', 'EdgeColor', 'none');

            % cSDT -> hue: more negative = bluer, more positive = redder.
            czMin = min(cSDT_unik);
            czMax = max(cSDT_unik);

            % signalCST -> contrast/vividness: low = lighter/paler, high = more vivid.
            sigMin = min(gaborCST_unik);
            sigMax = max(gaborCST_unik);

            for iB = 1:nBars_fig6b
                if ~isfinite(yMean(iB))
                    continue;
                end

                czVal = cSDT_unik(barCzIdx_fig6b(iB));
                sigVal = gaborCST_unik(barSigIdx_fig6b(iB));

                if czMax > czMin
                    czNorm = (czVal - czMin) / (czMax - czMin);
                else
                    czNorm = 0.5;
                end

                if sigMax > sigMin
                    sigNorm = (sigVal - sigMin) / (sigMax - sigMin);
                else
                    sigNorm = 0.5;
                end

                % Hue: blue (2/3) -> red (0) as cSDT goes negative -> positive.
                hCol = (1 - czNorm) * (2 / 3);

                % Saturation/Value: low signal = pale/light, high signal = vivid.
                sCol = 0.25 + 0.70 * sigNorm;   % 0.25 .. 0.95
                vCol = 0.98 - 0.13 * sigNorm;   % 0.98 .. 0.85

                b.CData(iB, :) = hsv2rgb([hCol, sCol, vCol]);
            end

            goodErr = isfinite(yMean) & isfinite(ySem);
            if any(goodErr)
                errorbar(xBar_fig6b(goodErr), yMean(goodErr), ySem(goodErr), 'k.', ...
                    'LineWidth', 1.0, 'CapSize', 0);
            end

            yline(0, 'k-', 'LineWidth', 0.8);
            xlim([0.4, nBars_fig6b + 0.6]);
            ylim([-1, 1]);

            % Visual separators between cSDT groups (each group has nSig_fig6b bars).
            for iSep = 1:max(nCz_fig6b - 1, 0)
                xline(iSep * nSig_fig6b + 0.5, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.9);
            end

            xticks(xBar_fig6b);
            if iRow == nRows_fig6b
                xticklabels(barLabels_fig6b);
                xtickangle(45);
            else
                xticklabels({});
            end

            if iCol == 1
                ylabel({sprintf('B%d', fitModel), 'r mean +/- 1sem'});
            end
            if iRow == 1
                title(sprintf('%s vs %s', pairLabels_fig6b{iP1}, pairLabels_fig6b{iP2}));
            end

            fxn_style_ax(gca, setting);
        end
    end

    sgtitle('Fig 6B: r by signalCST x cSDT | rows: Bfit=Bsim, cols: param pair', 'FontWeight', 'bold');
    saveas(h, fullfile(nameFolder_Figures_part4, 'FigS6B_paramCorr_summaryBars.png'));
    close(h);
end

%% Figure 7: Full vs. reduced models—Iteration-level ΔnLL distributions (Generating model is the full model B1)
%    For each reduced model Br, shows the per-iteration distribution of
%    ΔnLL(Br − B1) aggregated across conditions.
%    A positive ΔnLL means B1 fits the held-out data better on that iteration.
%    If distributions straddle 0, B1 barely wins and identifiability is weak.

fprintf('\n%s: Fig 7: Iteration-level delta-nLL distributions\n', string(datetime('now')))

Br_fig7 = setdiff(Bfit_unik, 1); % reduced models to compare against B1
nReduced_f7 = numel(Br_fig7);

if ~isfield(R, 'nLL_test_allIter') || ~any([R.iModelB_sim] == 1) || nReduced_f7 == 0
    fprintf('Skipping Fig 7: nLL_test_allIter not available or no B1-simulated data.\n');
else
    if isfield(R, 'nameIO')
        datasetIDs_f7 = string({R.nameIO});
    else
        datasetIDs_f7 = string(arrayfun(@(r) sprintf('Nm%g_Na%g_Ns%g_G%g_cSDT%g', ...
            r.Nmul_true, r.Nadd_true, r.Nshared_true, r.gaborCST, r.cSDT_true), R, 'UniformOutput', false));
    end
    idxSim1_f7 = [R.iModelB_sim] == 1;
    idsSim1_f7 = unique(datasetIDs_f7(idxSim1_f7));
    nConds_f7  = numel(idsSim1_f7);

    dnLL_iter_allConds = cell(nReduced_f7, nConds_f7);
    for iCond = 1:nConds_f7
        idxCond = idxSim1_f7 & datasetIDs_f7 == idsSim1_f7(iCond);
        Rcond   = R(idxCond);
        if isempty(Rcond), continue; end

        idxB1 = find([Rcond.iModelB_fit] == 1, 1);
        if isempty(idxB1) || isempty(Rcond(idxB1).nLL_test_allIter), continue; end
        nll_B1 = Rcond(idxB1).nLL_test_allIter(:);

        for iR = 1:nReduced_f7
            Br      = Br_fig7(iR);
            idxBr   = find([Rcond.iModelB_fit] == Br, 1);
            if isempty(idxBr) || isempty(Rcond(idxBr).nLL_test_allIter), continue; end
            nll_Br  = Rcond(idxBr).nLL_test_allIter(:);
            nIt     = min(numel(nll_B1), numel(nll_Br));
            dnLL_iter_allConds{iR, iCond} = nll_Br(1:nIt) - nll_B1(1:nIt);
        end
    end

    h = figure('Position', [100 100 280*nReduced_f7 480]);
    tiledlayout(1, nReduced_f7, 'TileSpacing', 'loose', 'Padding', 'loose');

    for iR = 1:nReduced_f7
        Br = Br_fig7(iR);
        nexttile; hold on;

        condMeds = nan(nConds_f7, 1);
        for iCond = 1:nConds_f7
            v = dnLL_iter_allConds{iR, iCond};
            if ~isempty(v), condMeds(iCond) = median(v, 'omitnan'); end
        end

        good = isfinite(condMeds);
        if any(good)
            xJit = 1 + 0.12 * (rand(sum(good), 1) - 0.5);
            scatter(xJit, condMeds(good), 28, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5);
            [gMed, gLb, gUb] = getCI(condMeds(good), 1, 1);
            errorbar(1, gMed, gMed - gLb, gUb - gMed, 'ko', ...
                'MarkerFaceColor', 'k', 'MarkerSize', 7, 'LineWidth', 1.8, 'CapSize', 0);
            winRate_f7 = mean(condMeds(good) > 0);
            text(1.35, gMed, sprintf('win=%.0f%%', winRate_f7*100), ...
                'FontSize', setting.fontSize, 'VerticalAlignment', 'middle');
        end

        yline(0, 'k--', 'LineWidth', 1.2);
        xlim([0.5 1.8]); xticks([]);
        ylabel('\Delta nLL (B_r - B_1) per condition median');
        title(sprintf('B%d - B1', Br));
        fxn_style_ax(gca, setting);
    end

    sgtitle('Fig 7: Full vs. reduced models\Delta nLL distributions | Bsim=Full model', 'FontWeight', 'bold');
    saveas(h, fullfile(nameFolder_Figures_part4, 'FigS7_dnLL_distributions.png'));
    close(h);
end

%% Figure 8: ΔnLL vs true omitted parameter value (Generating model is the full model B1)
%    For each reduced model Br, plots per-condition median ΔnLL(Br − B1) on the
%    y-axis against the true value of the parameter Br omits on the x-axis.
%    Expected: ΔnLL should increase as the omitted parameter grows away from 0.
%    Weak or flat relationship = identifiability problem even with large true parameter.

fprintf('\n%s: Fig 8: delta-nLL vs true omitted parameter\n', string(datetime('now')))

omitB_fig8      = [2,              3,            4,          5,          5,              6,          6,              7,          7        ];
omitParam_fig8  = {'Nshared',     'Nmul',       'Nadd',     'Nadd',     'Nshared',      'Nmul',     'Nshared',      'Nmul',     'Nadd'   };
trueField_fig8  = {'Nshared_true','Nmul_true',  'Nadd_true','Nadd_true','Nshared_true', 'Nmul_true','Nshared_true', 'Nmul_true','Nadd_true'};

keepPair = ismember(omitB_fig8, Bfit_unik);
omitB_fig8     = omitB_fig8(keepPair);
omitParam_fig8 = omitParam_fig8(keepPair);
trueField_fig8 = trueField_fig8(keepPair);
nPairs_f8 = numel(omitB_fig8);

if ~isfield(R, 'nLL_test_allIter') || ~any([R.iModelB_sim] == 1) || nPairs_f8 == 0
    fprintf('Skipping Fig 8: nLL_test_allIter not available or no valid model pairs.\n');
else
    if isfield(R, 'nameIO')
        datasetIDs_f8 = string({R.nameIO});
    else
        datasetIDs_f8 = string(arrayfun(@(r) sprintf('Nm%g_Na%g_Ns%g_G%g_cSDT%g', ...
            r.Nmul_true, r.Nadd_true, r.Nshared_true, r.gaborCST, r.cSDT_true), R, 'UniformOutput', false));
    end
    idxSim1_f8 = [R.iModelB_sim] == 1;
    idsSim1_f8 = unique(datasetIDs_f8(idxSim1_f8));

    paletteSrc = lines(max(nPairs_f8, 3));

    h = figure('Position', [100 100 380*nPairs_f8 420]);
    tiledlayout(1, nPairs_f8, 'TileSpacing', 'loose', 'Padding', 'loose');

    for iP = 1:nPairs_f8
        Br          = omitB_fig8(iP);
        tField      = trueField_fig8{iP};
        pLabel      = omitParam_fig8{iP};

        xvals = nan(numel(idsSim1_f8), 1);
        yvals = nan(numel(idsSim1_f8), 1);

        for iCond = 1:numel(idsSim1_f8)
            idxCond = idxSim1_f8 & datasetIDs_f8 == idsSim1_f8(iCond);
            Rcond   = R(idxCond);
            if isempty(Rcond), continue; end

            if isfield(Rcond, tField), xvals(iCond) = Rcond(1).(tField); end

            idxB1 = find([Rcond.iModelB_fit] == 1, 1);
            idxBr = find([Rcond.iModelB_fit] == Br, 1);
            if isempty(idxB1) || isempty(idxBr), continue; end
            nll_B1 = Rcond(idxB1).nLL_test_allIter(:);
            nll_Br = Rcond(idxBr).nLL_test_allIter(:);
            if isempty(nll_B1) || isempty(nll_Br), continue; end
            nIt = min(numel(nll_B1), numel(nll_Br));
            yvals(iCond) = median(nll_Br(1:nIt) - nll_B1(1:nIt), 'omitnan');
        end

        good = isfinite(xvals) & isfinite(yvals);
        nexttile; hold on;

        scatter(xvals(good), yvals(good), 40, paletteSrc(iP, :), ...
            'filled', 'MarkerFaceAlpha', 0.75);

        if sum(good) >= 3 && range(xvals(good)) > eps(max(abs(xvals(good)))) * numel(xvals(good))
            pCoef = polyfit(xvals(good) - mean(xvals(good)), yvals(good), 1);
            xFit  = linspace(min(xvals(good)), max(xvals(good)), 60);
            plot(xFit, polyval(pCoef, xFit - mean(xvals(good))), '-', 'Color', paletteSrc(iP,:), 'LineWidth', 1.8);
            [rho_f8, ~] = corr(xvals(good), yvals(good));
            text(0.05, 0.92, sprintf('r = %.2f', rho_f8), 'Units', 'normalized', ...
                'FontSize', setting.fontSize, 'Color', paletteSrc(iP,:));
        end

        yline(0, 'k--', 'LineWidth', 1.2);
        xlabel(sprintf('True %s', pLabel));
        ylabel(sprintf('\\Delta nLL median (B%d - B1)', Br));
        title(sprintf('B%d vs B1 | omitted: %s', Br, pLabel));
        fxn_style_ax(gca, setting);
    end

    sgtitle('Fig 8: \Delta nLL vs true omitted parameter | Bsim=B1', 'FontWeight', 'bold');
    saveas(h, fullfile(nameFolder_Figures_part4, 'FigS8_dnLL_vs_trueParam.png'));
    close(h);
end

%% Figure 9: Correlate template RMSE with each parameter's RMSE across simulating conditions

fprintf('\n%s: Fig 9: Template RMSE vs parameter RMSE correlations\n', string(datetime('now')))

% Use matched pairs (Bsim == Bfit) only, so each row is a consistent condition.
% For each panel: scatter template_rmse (x) vs that parameter's RMSE (y) and annotate r.
% Panels 1-4: RMSE of fitted noise params + DV criterion (never excluded).
% Panel 5:    True SDT criterion |cSDT_true| — shows if criterion bias predicts template error.

% --- RMSE-based panels (one per fitted parameter) ---
rmseParamNames_fig9  = {'Nmul',       'Nadd',       'Nshared',       'criterion_DV'};
rmseParamLabels_fig9 = {'N_{mul}',    'N_{add}',    'N_{shared}',    'Criterion DV'};

% Keep only those with at least one finite value; criterion_DV is always kept.
Rmatch_fig9 = R([R.iModelB_sim] == [R.iModelB_fit]);
keepRmse_fig9 = false(1, numel(rmseParamNames_fig9));
for iP9 = 1:numel(rmseParamNames_fig9)
    fld = sprintf('%s_rmse', rmseParamNames_fig9{iP9});
    if strcmp(rmseParamNames_fig9{iP9}, 'criterion_DV')
        keepRmse_fig9(iP9) = true;   % never exclude criterion DV
    elseif isfield(Rmatch_fig9, fld) && any(isfinite([Rmatch_fig9.(fld)]))
        keepRmse_fig9(iP9) = true;
    end
end
rmseParamNames_fig9  = rmseParamNames_fig9(keepRmse_fig9);
rmseParamLabels_fig9 = rmseParamLabels_fig9(keepRmse_fig9);

nPanels_fig9 = numel(rmseParamNames_fig9) + 1;   % +1 for SDT criterion panel
nCols_fig9   = min(nPanels_fig9, 3);
nRows_fig9   = ceil(nPanels_fig9 / nCols_fig9);
cmap_fig9    = lines(numel(Bfit_unik));

h = figure('Position', [100 100 380*nCols_fig9 320*nRows_fig9]);
tiledlayout(nRows_fig9, nCols_fig9, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panels: template RMSE vs each parameter RMSE ---
for iP9 = 1:numel(rmseParamNames_fig9)
    pName9    = rmseParamNames_fig9{iP9};
    pLabel9   = rmseParamLabels_fig9{iP9};
    rmseField = sprintf('%s_rmse', pName9);

    nexttile; hold on;

    hLines_fig9 = gobjects(numel(Bfit_unik), 1);
    for iBfit9 = 1:numel(Bfit_unik)
        fitModel9 = Bfit_unik(iBfit9);
        Rsub9 = R([R.iModelB_sim] == fitModel9 & [R.iModelB_fit] == fitModel9);
        if isempty(Rsub9) || ~isfield(Rsub9, rmseField), continue; end

        x9 = [Rsub9.template_rmse]';
        y9 = [Rsub9.(rmseField)]';
        good9 = isfinite(x9) & isfinite(y9);
        if sum(good9) < 2, continue; end

        col9 = cmap_fig9(iBfit9, :);
        hLines_fig9(iBfit9) = scatter(x9(good9), y9(good9), 28, col9, ...
            'filled', 'MarkerFaceAlpha', 0.65);

        if numel(unique(x9(good9))) >= 2
            pCoef9 = polyfit(x9(good9), y9(good9), 1);
            xFit9  = linspace(min(x9(good9)), max(x9(good9)), 80);
            plot(xFit9, polyval(pCoef9, xFit9), '-', 'Color', col9, ...
                'LineWidth', 1.6, 'HandleVisibility', 'off');
        end

        [r9, ~] = corr(x9(good9), y9(good9));
        hLines_fig9(iBfit9).DisplayName = sprintf('B%d: r = %.2f', fitModel9, r9);
    end % iBfit9

    xlabel('Template RMSE');
    ylabel(sprintf('%s RMSE', pLabel9));
    title(sprintf('Template RMSE vs %s RMSE', pLabel9));
    fxn_style_ax(gca, setting);
    set(gca, 'YGrid', 'off');
    axis square;
    validH9 = hLines_fig9(isgraphics(hLines_fig9));
    if ~isempty(validH9)
        legend(validH9, 'Location', 'best', 'Box', 'off', 'FontSize', max(setting.fontSize-1,7));
    end
end % iP9

% --- Final panel: template RMSE vs true SDT criterion |cSDT_true| ---
nexttile; hold on;
hLines_sdt = gobjects(numel(Bfit_unik), 1);
for iBfit9 = 1:numel(Bfit_unik)
    fitModel9 = Bfit_unik(iBfit9);
    Rsub9 = R([R.iModelB_sim] == fitModel9 & [R.iModelB_fit] == fitModel9);
    if isempty(Rsub9) || ~isfield(Rsub9, 'cSDT_true'), continue; end

    x9   = [Rsub9.template_rmse]';
    y9   = abs([Rsub9.cSDT_true]');
    good9 = isfinite(x9) & isfinite(y9);
    if sum(good9) < 2, continue; end

    col9 = cmap_fig9(iBfit9, :);
    hLines_sdt(iBfit9) = scatter(x9(good9), y9(good9), 28, col9, ...
        'filled', 'MarkerFaceAlpha', 0.65);

    if numel(unique(x9(good9))) >= 2
        pCoef9 = polyfit(x9(good9), y9(good9), 1);
        xFit9  = linspace(min(x9(good9)), max(x9(good9)), 80);
        plot(xFit9, polyval(pCoef9, xFit9), '-', 'Color', col9, ...
            'LineWidth', 1.6, 'HandleVisibility', 'off');
    end

    [r9, ~] = corr(x9(good9), y9(good9));
    hLines_sdt(iBfit9).DisplayName = sprintf('B%d r=%.2f', fitModel9, r9);
end
xlabel('Template RMSE');
ylabel('|SDT criterion| (|Cz_{true}|)');
title('Template RMSE vs SDT criterion');
fxn_style_ax(gca, setting);
set(gca, 'YGrid', 'off');
axis square;
validH_sdt = hLines_sdt(isgraphics(hLines_sdt));
if ~isempty(validH_sdt)
    legend(validH_sdt, 'Location', 'best', 'Box', 'off', 'FontSize', max(setting.fontSize-1,7));
end

sgtitle('Fig 9: Template RMSE vs Parameter RMSE | Matched Bsim=Bfit', 'FontWeight', 'bold');
saveas(h, fullfile(nameFolder_Figures_part4, 'FigS9_Corr_tempRMSE_paramRMSE.png'));
close(h);

%% Local helpers  (only functions used 2+ times)

function info = fxn_parse_nameIO(nameIO)

info = struct();

tok = regexp(nameIO, ...
    'IO_cN([-\d\.]+)_cG([-\d\.]+)_nT([A-Za-z0-9\.\-]+)_Nm([A-Za-z0-9\.\-]+)_Na([A-Za-z0-9\.\-]+)_Ns([A-Za-z0-9\.\-]+)_cSDT([-\d\.]+)_cont([-\d\.]+)_whiten([-\d\.]+)_R(\d+)_\d+\d+_B(\d+)', ...
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
info.cSDT_true = str2double(tok{7});
info.C_contribution = str2double(tok{8});
info.lambda_whiten = str2double(tok{9});
info.flag_regressType = str2double(tok{10});
info.iModelB_sim = str2double(tok{11});
end

function x = fxn_parse_num2exp(str_in)
% parse strings like '8e3', '0', '10', etc.
x = str2double(strrep(str_in, 'p', '.'));
end

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
end
end

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

function fxn_style_ax(ax, P)
box(ax, 'on');
set(ax, 'FontSize', P.fontSize, 'LineWidth', P.axisLineWidth, ...
    'YGrid', 'on', 'GridAlpha', P.gridAlpha, 'TickDir', P.tickDir);
end

function S = fxn_collapse_scalar_by_x_iterCI(Rin, xField, yIterField, xLevels, collapseFcn)
S = struct('x', [], 'y', [], 'lb', [], 'ub', [], 'n', []);

x = [];
y = [];
lb = [];
ub = [];
n = [];

for i = 1:numel(xLevels)
    idx = [Rin.(xField)] == xLevels(i);
    Rsub = Rin(idx);
    if isempty(Rsub), continue; end

    yRow = nan(numel(Rsub), 1);
    lbRow = nan(numel(Rsub), 1);
    ubRow = nan(numel(Rsub), 1);

    for iRow = 1:numel(Rsub)
        vals = Rsub(iRow).(yIterField);
        vals = vals(:);
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end
        [yRow(iRow), lbRow(iRow), ubRow(iRow)] = getCI(vals, 1, 1);
    end

    good = isfinite(yRow);
    if ~any(good), continue; end

    x(end+1,1) = xLevels(i); %#ok<AGROW>
    y(end+1,1) = collapseFcn(yRow(good), 'omitnan'); %#ok<AGROW>
    lb(end+1,1) = collapseFcn(lbRow(good), 'omitnan'); %#ok<AGROW>
    ub(end+1,1) = collapseFcn(ubRow(good), 'omitnan'); %#ok<AGROW>
    n(end+1,1) = sum(good); %#ok<AGROW>
end

S.x = x;
S.y = y;
S.lb = lb;
S.ub = ub;
S.n = n;
end


function [xAxis, yTrue, curves] = fxn_collapse_tuning_panel(Rin, vField, vLevels, domainName)
xAxis = [];
yTrue = [];
curves = struct('level', {}, 'y', {}, 'lb', {}, 'ub', {}, 'n', {});

if isempty(Rin), return; end

switch domainName
    case 'ORI'
        fieldTrue = 'margORI_true';
        fieldEstIter = 'margORI_est_iter_norm';
        nCh = numel(Rin(1).margORI_true);
    case 'SF'
        fieldTrue = 'margSF_true';
        fieldEstIter = 'margSF_est_iter_norm';
        nCh = numel(Rin(1).margSF_true);
    otherwise
        error('Unknown domain %s', domainName);
end

xAxis = 1:nCh;

T = vertcat(Rin.(fieldTrue));
yTrue = mean(T, 1, 'omitnan');

for i = 1:numel(vLevels)
    idx = [Rin.(vField)] == vLevels(i);
    Rsub = Rin(idx);
    if isempty(Rsub), continue; end

    yRow = nan(numel(Rsub), nCh);
    lbRow = nan(numel(Rsub), nCh);
    ubRow = nan(numel(Rsub), nCh);

    for iRowSub = 1:numel(Rsub)
        Yiter = Rsub(iRowSub).(fieldEstIter);
        if isempty(Yiter), continue; end

        for iCh = 1:nCh
            vals = Yiter(:, iCh);
            vals = vals(isfinite(vals));
            if isempty(vals), continue; end
            [yRow(iRowSub, iCh), lbRow(iRowSub, iCh), ubRow(iRowSub, iCh)] = getCI(vals, 1, 1);
        end
    end

    y = mean(yRow, 1, 'omitnan');
    lb = mean(lbRow, 1, 'omitnan');
    ub = mean(ubRow, 1, 'omitnan');

    curves(end+1).level = vLevels(i); %#ok<AGROW>
    curves(end).y = y;
    curves(end).lb = lb;
    curves(end).ub = ub;
    curves(end).n = sum(any(isfinite(yRow), 2));
end
end

function C = fxn_collapse_curve(Rin, xField, yField, nGrid, collapseFcn)
C = struct('x', [], 'y', [], 'lb', [], 'ub', [], 'n', 0);
if isempty(Rin), return; end

xMin = inf; xMax = -inf;
keep = false(numel(Rin),1);

for i = 1:numel(Rin)
    x = Rin(i).(xField);
    y = Rin(i).(yField);
    if isempty(x) || isempty(y), continue; end
    good = isfinite(x) & isfinite(y);
    if ~any(good), continue; end
    xMin = min(xMin, min(x(good)));
    xMax = max(xMax, max(x(good)));
    keep(i) = true;
end

Rin = Rin(keep);
if isempty(Rin) || ~isfinite(xMin) || ~isfinite(xMax) || xMin == xMax
    return
end

xGrid = linspace(xMin, xMax, nGrid);
Yall = nan(numel(Rin), nGrid);

for i = 1:numel(Rin)
    x = Rin(i).(xField);
    y = Rin(i).(yField);

    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);

    if numel(unique(x)) < 2, continue; end
    [x, ia] = unique(x);
    y = y(ia);

    Yall(i, :) = interp1(x, y, xGrid, 'linear', nan);
end

goodRows = any(isfinite(Yall), 2);
Yall = Yall(goodRows, :);
if isempty(Yall), return; end

C.x = xGrid;
C.y = collapseFcn(Yall, 1, 'omitnan');
C.lb = nan(1, nGrid);
C.ub = nan(1, nGrid);

for i = 1:nGrid
    vals = Yall(:, i);
    vals = vals(isfinite(vals));
    if isempty(vals), continue; end
    [~, C.lb(i), C.ub(i)] = getCI(vals, 1, 1);
end
C.n = size(Yall,1);
end

function [C, nGridCounts] = fxn_collapse_curve_with_counts(Rin, xField, yField, wField, nGrid, collapseFcn)
C = fxn_collapse_curve(Rin, xField, yField, nGrid, collapseFcn);
nGridCounts = [];

if isempty(C.x) || isempty(Rin) || ~isfield(Rin, wField)
    return
end

xGrid = C.x(:)';
nG = numel(xGrid);
if nG < 1
    return
end

if nG >= 2
    edges = [-inf, (xGrid(1:end-1) + xGrid(2:end)) / 2, inf];
else
    edges = [-inf, inf];
end

nGridCounts = zeros(1, nG);
for i = 1:numel(Rin)
    x = Rin(i).(xField);
    w = Rin(i).(wField);
    if isempty(x) || isempty(w), continue; end

    x = x(:);
    w = w(:);
    n = min(numel(x), numel(w));
    x = x(1:n);
    w = w(1:n);

    good = isfinite(x) & isfinite(w) & w > 0;
    if ~any(good), continue; end
    x = x(good);
    w = w(good);

    ib = discretize(x, edges);
    for iBin = 1:nG
        nGridCounts(iBin) = nGridCounts(iBin) + sum(w(ib == iBin), 'omitnan');
    end
end

if ~any(isfinite(nGridCounts) & nGridCounts > 0)
    nGridCounts = [];
end
end

function Rout = fxn_attach_medians_from_iter(Rin, namesMetrics_behav)
% Reconstruct summary fields used by plotting from per-iteration values.
Rout = Rin;
if isempty(Rout)
    return
end

paramNames_IN = {'Nmul', 'Nadd', 'Nshared'};

for iR = 1:numel(Rout)
    % Template-level summaries
    if isfield(Rout, 'template_rmse_iter') && ~isempty(Rout(iR).template_rmse_iter)
        [Rout(iR).template_rmse, ~, ~] = getCI(Rout(iR).template_rmse_iter(:), 1, 1);
    else
        Rout(iR).template_rmse = nan;
    end
    if isfield(Rout, 'template_R2_iter') && ~isempty(Rout(iR).template_R2_iter)
        [Rout(iR).template_R2, ~, ~] = getCI(Rout(iR).template_R2_iter(:), 1, 1);
    else
        Rout(iR).template_R2 = nan;
    end

    if isfield(Rout, 'margORI_est_iter_norm') && ~isempty(Rout(iR).margORI_est_iter_norm)
        [Rout(iR).margORI_est_med, Rout(iR).margORI_est_lb, Rout(iR).margORI_est_ub] = getCI(Rout(iR).margORI_est_iter_norm, 1, 1);
    else
        Rout(iR).margORI_est_med = nan(1, numel(Rout(iR).margORI_true));
        Rout(iR).margORI_est_lb = nan(1, numel(Rout(iR).margORI_true));
        Rout(iR).margORI_est_ub = nan(1, numel(Rout(iR).margORI_true));
    end
    if isfield(Rout, 'margSF_est_iter_norm') && ~isempty(Rout(iR).margSF_est_iter_norm)
        [Rout(iR).margSF_est_med, Rout(iR).margSF_est_lb, Rout(iR).margSF_est_ub] = getCI(Rout(iR).margSF_est_iter_norm, 1, 1);
    else
        Rout(iR).margSF_est_med = nan(1, numel(Rout(iR).margSF_true));
        Rout(iR).margSF_est_lb = nan(1, numel(Rout(iR).margSF_true));
        Rout(iR).margSF_est_ub = nan(1, numel(Rout(iR).margSF_true));
    end

    % Data metric summaries
    [Rout(iR).pYES_data_med, Rout(iR).pYES_data_lb, Rout(iR).pYES_data_ub] = getCI(Rout(iR).pYES_data_iter(:), 1, 1);
    [Rout(iR).pC_data_med, Rout(iR).pC_data_lb, Rout(iR).pC_data_ub] = getCI(Rout(iR).pC_data_iter(:), 1, 1);
    [Rout(iR).pA_data_med, Rout(iR).pA_data_lb, Rout(iR).pA_data_ub] = getCI(Rout(iR).pA_data_iter(:), 1, 1);

    [Rout(iR).DV_allBins_med, Rout(iR).DV_allBins_lb, Rout(iR).DV_allBins_ub] = fxn_ci_cols_triplet(Rout(iR).DV_allBins_iter);
    [Rout(iR).nTrials_allBins_med, Rout(iR).nTrials_allBins_lb, Rout(iR).nTrials_allBins_ub] = fxn_ci_cols_triplet(Rout(iR).nTrials_allBins_iter);
    [Rout(iR).pYES_data_curve_med, Rout(iR).pYES_data_curve_lb, Rout(iR).pYES_data_curve_ub] = fxn_ci_cols_triplet(Rout(iR).pYES_data_curve_iter);
    [Rout(iR).pC_data_curve_med, Rout(iR).pC_data_curve_lb, Rout(iR).pC_data_curve_ub] = fxn_ci_cols_triplet(Rout(iR).pC_data_curve_iter);
    [Rout(iR).pA_data_curve_med, Rout(iR).pA_data_curve_lb, Rout(iR).pA_data_curve_ub] = fxn_ci_cols_triplet(Rout(iR).pA_data_curve_iter);

    % Fit metric summaries
    for iMetric = 1:numel(namesMetrics_behav)
        m = namesMetrics_behav{iMetric};
        rmseIterField = sprintf('%s_rmse_allIter', m);
        R2IterField = sprintf('%s_R2_allIter', m);

        if isfield(Rout, rmseIterField) && ~isempty(Rout(iR).(rmseIterField))
            [Rout(iR).(sprintf('%s_rmse', m)), ~, ~] = getCI(Rout(iR).(rmseIterField)(:), 1, 1);
        else
            Rout(iR).(sprintf('%s_rmse', m)) = nan;
        end

        if isfield(Rout, R2IterField) && ~isempty(Rout(iR).(R2IterField))
            [Rout(iR).(sprintf('%s_R2', m)), ~, ~] = getCI(Rout(iR).(R2IterField)(:), 1, 1);
        else
            Rout(iR).(sprintf('%s_R2', m)) = nan;
        end
    end

    [Rout(iR).pYES_pred_med, Rout(iR).pYES_pred_lb, Rout(iR).pYES_pred_ub] = getCI(Rout(iR).pYES_pred_iter(:), 1, 1);
    [Rout(iR).pC_pred_med, Rout(iR).pC_pred_lb, Rout(iR).pC_pred_ub] = getCI(Rout(iR).pC_pred_iter(:), 1, 1);
    [Rout(iR).pA_pred_med, Rout(iR).pA_pred_lb, Rout(iR).pA_pred_ub] = getCI(Rout(iR).pA_pred_iter(:), 1, 1);

    [Rout(iR).pYES_pred_curve_med, Rout(iR).pYES_pred_curve_lb, Rout(iR).pYES_pred_curve_ub] = fxn_ci_cols_triplet(Rout(iR).pYES_pred_curve_iter);
    [Rout(iR).pC_pred_curve_med, Rout(iR).pC_pred_curve_lb, Rout(iR).pC_pred_curve_ub] = fxn_ci_cols_triplet(Rout(iR).pC_pred_curve_iter);
    [Rout(iR).pA_pred_curve_med, Rout(iR).pA_pred_curve_lb, Rout(iR).pA_pred_curve_ub] = fxn_ci_cols_triplet(Rout(iR).pA_pred_curve_iter);

    if isfield(Rout, 'nLL_test_allIter') && ~isempty(Rout(iR).nLL_test_allIter)
        [Rout(iR).nLL_med, Rout(iR).nLL_lb, Rout(iR).nLL_ub] = getCI(Rout(iR).nLL_test_allIter(:), 1, 1);
    else
        Rout(iR).nLL_med = nan;
        Rout(iR).nLL_lb = nan;
        Rout(iR).nLL_ub = nan;
    end

    % Criterion summaries
    if isfield(Rout, 'criterion_DV_rmse_allIter') && ~isempty(Rout(iR).criterion_DV_rmse_allIter)
        [Rout(iR).criterion_DV_rmse, ~, ~] = getCI(Rout(iR).criterion_DV_rmse_allIter(:), 1, 1);
    else
        Rout(iR).criterion_DV_rmse = nan;
    end
    if isfield(Rout, 'criterion_DV_est_allIter') && ~isempty(Rout(iR).criterion_DV_est_allIter)
        [Rout(iR).criterion_DV_est_med, Rout(iR).criterion_DV_est_lb, Rout(iR).criterion_DV_est_ub] = getCI(Rout(iR).criterion_DV_est_allIter(:), 1, 1);
    else
        Rout(iR).criterion_DV_est_med = nan;
        Rout(iR).criterion_DV_est_lb = nan;
        Rout(iR).criterion_DV_est_ub = nan;
    end

    % Internal-noise parameter summaries
    for iP = 1:numel(paramNames_IN)
        pName = paramNames_IN{iP};
        rmseIterField = sprintf('%s_rmse_allIter', pName);
        estIterField = sprintf('%s_est_allIter', pName);

        if isfield(Rout, rmseIterField) && ~isempty(Rout(iR).(rmseIterField))
            [Rout(iR).(sprintf('%s_rmse', pName)), ~, ~] = getCI(Rout(iR).(rmseIterField)(:), 1, 1);
        else
            Rout(iR).(sprintf('%s_rmse', pName)) = nan;
        end

        if isfield(Rout, estIterField) && ~isempty(Rout(iR).(estIterField))
            [Rout(iR).(sprintf('%s_est_med', pName)), Rout(iR).(sprintf('%s_est_lb', pName)), Rout(iR).(sprintf('%s_est_ub', pName))] = getCI(Rout(iR).(estIterField)(:), 1, 1);
        else
            Rout(iR).(sprintf('%s_est_med', pName)) = nan;
            Rout(iR).(sprintf('%s_est_lb', pName)) = nan;
            Rout(iR).(sprintf('%s_est_ub', pName)) = nan;
        end
    end
end
end

function medRow = fxn_ci_cols(vals)
if isempty(vals)
    medRow = [];
    return
end

if isvector(vals)
    medRow = getCI(vals(:), 1, 1);
    medRow = medRow(:)';
    return
end

nCol = size(vals, 2);
medRow = nan(1, nCol);
for iCol = 1:nCol
    [medRow(iCol), ~, ~] = getCI(vals(:, iCol), 1, 1);
end
end

function [medRow, lbRow, ubRow] = fxn_ci_cols_triplet(vals)
if isempty(vals)
    medRow = [];
    lbRow = [];
    ubRow = [];
    return
end

if isvector(vals)
    [medTmp, lbTmp, ubTmp] = getCI(vals(:), 1, 1);
    medRow = medTmp(:)';
    lbRow = lbTmp(:)';
    ubRow = ubTmp(:)';
    return
end

nCol = size(vals, 2);
medRow = nan(1, nCol);
lbRow = nan(1, nCol);
ubRow = nan(1, nCol);
for iCol = 1:nCol
    [medRow(iCol), lbRow(iCol), ubRow(iCol)] = getCI(vals(:, iCol), 1, 1);
end
end

function fxn_plot_connected_dots_with_err(x, y, lb, ub, colorRGB, P)
good = isfinite(x) & isfinite(y) & isfinite(lb) & isfinite(ub);
if ~any(good), return; end

x = x(good);
y = y(good);
lb = lb(good);
ub = ub(good);

plot(x, y, '-', 'Color', colorRGB, 'LineWidth', P.lineWidth);
for iPt = 1:numel(x)
    line([x(iPt) x(iPt)], [lb(iPt) ub(iPt)], 'Color', colorRGB, 'LineWidth', P.errLineWidth);
end
plot(x, y, 'o', ...
    'MarkerSize', P.markerSize, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', colorRGB, ...
    'LineStyle', 'none', ...
    'LineWidth', P.lineWidth);
end

function [r2, rho] = fxn_curve_fit_metrics(xRef, yRef, xPred, yPred, w)
% Interpolate prediction onto reference x grid, then compute R² and Pearson's r.
% Optional w: trial counts per bin for weighted R² / rho.
yInterp = interp1(xPred, yPred, xRef, 'linear', nan);
good = isfinite(yRef) & isfinite(yInterp);
if sum(good) < 2
    r2 = nan; rho = nan; return
end
yR = yRef(good);
yP = yInterp(good);
if nargin >= 5 && ~isempty(w) && numel(w) == numel(xRef)
    wg = w(good); wg = wg(:); wg(~isfinite(wg) | wg < 0) = 0;
    wSum = sum(wg);
    if wSum > 0
        wg = wg / wSum;
        yR_col = yR(:);
        yP_col = yP(:);
        yR_wmean = sum(wg .* yR_col);
        ss_res = sum(wg .* (yR_col - yP_col).^2);
        ss_tot = sum(wg .* (yR_col - yR_wmean).^2);
        r2  = 1 - ss_res / max(ss_tot, eps);
        yP_wmean = sum(wg .* yP_col);
        cov_w = sum(wg .* (yR_col - yR_wmean) .* (yP_col - yP_wmean));
        varR_w = sum(wg .* (yR_col - yR_wmean).^2);
        varP_w = sum(wg .* (yP_col - yP_wmean).^2);
        rho = cov_w / sqrt(max(varR_w * varP_w, eps));
        return
    end
end
ss_res = sum((yR - yP).^2);
ss_tot = sum((yR - mean(yR)).^2);
r2  = 1 - ss_res / ss_tot;
rho = corr(yR(:), yP(:));
end

function rmse = fxn_curve_rmse(xRef, yRef, xPred, yPred, w)
yInterp = interp1(xPred, yPred, xRef, 'linear', nan);
good = isfinite(yRef) & isfinite(yInterp);
if sum(good) < 1
    rmse = nan;
    return
end

if nargin >= 5 && ~isempty(w) && numel(w) == numel(xRef)
    wg = w(good); wg = wg(:) / sum(wg(isfinite(wg)));
    rmse = sqrt(sum(wg' .* (yRef(good) - yInterp(good)).^2, 'omitnan'));
    assert(isscalar(rmse), 'ALERT: rmse is not a scalar but a vector!!')
else
    rmse = sqrt(mean((yRef(good) - yInterp(good)).^2));
end
end

function plot_ranked_categorical(vals, xLabelText, groups, groupCmap, showLegend)
% Parse vals → string labels, tracking valid entries with keep mask
if isnumeric(vals) || islogical(vals)
    vals = vals(:);
    keep = isfinite(vals);
    labels = string(vals(keep));
elseif iscell(vals)
    raw = string(vals(:));
    keep = strlength(raw) > 0;
    labels = raw(keep);
elseif isstring(vals) || iscategorical(vals)
    raw = string(vals(:));
    keep = strlength(raw) > 0;
    labels = raw(keep);
else
    error('Unsupported type for ranked categorical plot.');
end

if isempty(labels)
    text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
    xlim([0 1]); ylim([0 1]); return
end

hasGroups = nargin >= 3 && ~isempty(groups);
if nargin < 5 || isempty(showLegend)
    showLegend = false;
end
if hasGroups
    grp = string(groups(:));
    grp = grp(keep);
    uGrp = unique(grp);
    nGrp = numel(uGrp);
    if nargin < 4 || isempty(groupCmap)
        groupCmap = lines(nGrp);
    end
end

[uCat, ~, ic] = unique(labels);
counts = accumarray(ic, 1);
pct = counts / sum(counts) * 100;
[~, idxSort] = sort(pct, 'descend');
u_sorted = uCat(idxSort);
nCat = numel(u_sorted);

if ~hasGroups
    pct_sorted = pct(idxSort);
    bar(1:nCat, pct_sorted, 0.7);
    pctTop = pct_sorted;
else
    % Build count matrix: rows = categories (sorted), cols = groups
    countMat = zeros(nCat, nGrp);
    for iCat = 1:nCat
        for iGrp = 1:nGrp
            countMat(iCat, iGrp) = sum(labels == u_sorted(iCat) & grp == uGrp(iGrp));
        end
    end
    pctMat = countMat / sum(countMat(:)) * 100;
    pctTop = sum(pctMat, 2);

    b = bar(1:nCat, pctMat, 0.7, 'stacked');
    for iGrp = 1:nGrp
        b(iGrp).FaceColor = groupCmap(iGrp, :);
        b(iGrp).EdgeColor = 'none';
    end
    if showLegend
        legend(b, cellstr(uGrp), 'Location', 'northeast', 'FontSize', 7, 'Box', 'off');
    end
end

xticks(1:nCat);
xticklabels(u_sorted);
% xtickangle(45);
xlabel(xLabelText);
ylim([0, max(pctTop) * 1.15]);

for i = 1:nCat
    text(i, pctTop(i), sprintf(' %.0f%%', pctTop(i)), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 9);
end
end

function vals = fxn_field_values_to_string(Rin, fieldName)
sampleVal = Rin(1).(fieldName);
if isnumeric(sampleVal) || islogical(sampleVal)
    vals = string([Rin.(fieldName)]);
elseif isstring(sampleVal) || ischar(sampleVal) || iscategorical(sampleVal)
    vals = string({Rin.(fieldName)});
else
    vals = string({Rin.(fieldName)});
end
end

function P = fxn_collect_fig4_panel_points(Rsub, pName, gaborCST_unik, cSDT_unik, collapseFcn, nBins_Part4)
P = struct('x', [], 'y', [], 'lb', [], 'ub', [], 'n', []);
for iSig = 1:numel(gaborCST_unik)
    for iCz = 1:numel(cSDT_unik)
        S = fxn_collapse_param_byLevel( ...
            Rsub, pName, ...
            collapseFcn, nBins_Part4, gaborCST_unik(iSig), cSDT_unik(iCz));
        if isempty(S.x), continue; end
        P.x = [P.x; S.x(:)]; %#ok<AGROW>
        P.y = [P.y; S.y(:)]; %#ok<AGROW>
        P.lb = [P.lb; S.lb(:)]; %#ok<AGROW>
        P.ub = [P.ub; S.ub(:)]; %#ok<AGROW>
        P.n = [P.n; S.n(:)]; %#ok<AGROW>
    end
end
end

function [rho, slope, bias, rmse] = fxn_param_recovery_stats(x, y)
good = isfinite(x) & isfinite(y);
x = x(good);
y = y(good);
if numel(x) < 2
    rho = nan; slope = nan; bias = nan; rmse = nan; return
end

rho = corr(x(:), y(:));
p = polyfit(x(:), y(:), 1);
slope = p(1);
bias = mean(y(:) - x(:), 'omitnan');
rmse = sqrt(mean((y(:) - x(:)).^2, 'omitnan'));
end

function S = fxn_collapse_param_byLevel(Rin, pName, collapseFcn, nBins_Part4, signalVal, czVal)
% Keep the examined parameter levels separate.
% Collapse only across the other parameters.
% Filter by signalCST x cSDT_true when signalVal/czVal are supplied.

S = struct('x', [], 'y', [], 'lb', [], 'ub', [], 'n', []);

if nargin >= 6
    idx = [Rin.gaborCST] == signalVal & [Rin.cSDT_true] == czVal;
    Rsub = Rin(idx);
else
    Rsub = Rin;
end
if isempty(Rsub), return; end

switch pName
    case 'criterion_DV'
        S = fxn_bin_criterion_values(Rsub, nBins_Part4, collapseFcn);
        return
    otherwise
        xField  = sprintf('%s_true', pName);
        yField  = sprintf('%s_est_med', pName);
        lbField = sprintf('%s_est_lb', pName);
        ubField = sprintf('%s_est_ub', pName);
end

xvals_all  = [Rsub.(xField)]';
yvals_all  = [Rsub.(yField)]';
lbvals_all = [Rsub.(lbField)]';
ubvals_all = [Rsub.(ubField)]';

good = isfinite(xvals_all) & isfinite(yvals_all);
xvals_all  = xvals_all(good);
yvals_all  = yvals_all(good);
lbvals_all = lbvals_all(good);
ubvals_all = ubvals_all(good);

if isempty(xvals_all), return; end

% group by the examined parameter's true levels
xLevels = unique(xvals_all);

x_out  = nan(numel(xLevels), 1);
y_out  = nan(numel(xLevels), 1);
lb_out = nan(numel(xLevels), 1);
ub_out = nan(numel(xLevels), 1);
n_out  = nan(numel(xLevels), 1);

for iLev = 1:numel(xLevels)
    idxLev = xvals_all == xLevels(iLev);

    x_out(iLev)  = xLevels(iLev);
    y_out(iLev)  = collapseFcn(yvals_all(idxLev), 'omitnan');
    lb_out(iLev) = collapseFcn(lbvals_all(idxLev), 'omitnan');
    ub_out(iLev) = collapseFcn(ubvals_all(idxLev), 'omitnan');
    n_out(iLev)  = sum(idxLev);
end

S.x  = x_out;
S.y  = y_out;
S.lb = lb_out;
S.ub = ub_out;
S.n  = n_out;
end

function S = fxn_bin_criterion_values(Rsub, nBins_Part4, collapseFcn)
S = struct('x', [], 'y', [], 'lb', [], 'ub', [], 'n', []);

xvals_all  = [Rsub.criterion_DV_true]';
yvals_all  = [Rsub.criterion_DV_est_med]';
lbvals_all = [Rsub.criterion_DV_est_lb]';
ubvals_all = [Rsub.criterion_DV_est_ub]';

good = isfinite(xvals_all) & isfinite(yvals_all);
xvals_all  = xvals_all(good);
yvals_all  = yvals_all(good);
lbvals_all = lbvals_all(good);
ubvals_all = ubvals_all(good);

if isempty(xvals_all), return; end

xMin = min(xvals_all);
xMax = max(xvals_all);
if ~isfinite(xMin) || ~isfinite(xMax), return; end

if xMin == xMax
    S.x = round(xMin);
    S.y = collapseFcn(yvals_all, 'omitnan');
    S.lb = collapseFcn(lbvals_all, 'omitnan');
    S.ub = collapseFcn(ubvals_all, 'omitnan');
    S.n = numel(yvals_all);
    return
end

binEdges = linspace(xMin, xMax, nBins_Part4 + 1);
binCtrs = round((binEdges(1:end-1) + binEdges(2:end)) / 2);
binIdx = discretize(xvals_all, binEdges);
binIdx(xvals_all == xMax) = nBins_Part4;

uCtrs = unique(binCtrs(:));
x_out  = nan(numel(uCtrs), 1);
y_out  = nan(numel(uCtrs), 1);
lb_out = nan(numel(uCtrs), 1);
ub_out = nan(numel(uCtrs), 1);
n_out  = nan(numel(uCtrs), 1);

for iCtr = 1:numel(uCtrs)
    idxCtr = ismember(binIdx, find(binCtrs == uCtrs(iCtr)));
    if ~any(idxCtr), continue; end

    x_out(iCtr)  = uCtrs(iCtr);
    y_out(iCtr)  = collapseFcn(yvals_all(idxCtr), 'omitnan');
    lb_out(iCtr) = collapseFcn(lbvals_all(idxCtr), 'omitnan');
    ub_out(iCtr) = collapseFcn(ubvals_all(idxCtr), 'omitnan');
    n_out(iCtr)  = sum(idxCtr);
end

keep = isfinite(x_out) & isfinite(y_out);
S.x  = x_out(keep);
S.y  = y_out(keep);
S.lb = lb_out(keep);
S.ub = ub_out(keep);
S.n  = n_out(keep);
end
