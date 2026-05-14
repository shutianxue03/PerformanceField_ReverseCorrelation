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

%% settings
flag_whitenDV=0; % NEEDED for SX_RC1_setting!!
str_part = 'whitenDV0'; % <-- change if needed
iModelA_fit = [1]; %1=use data-derived template; 2=use true template
nBfit = 7;          % number of B-model variants to evaluate (each uses a different internal-noise structure)
nBins_Part4 = 3; % define bins for collapsing parameter recovery points; use 3 for main text, 5 for Supp

%--------------%
SX_RC1_setting;
%--------------%
nameFolder_Data = sprintf('%s/Data_%s', nameFolder_server, str_part);
nameFile_R = sprintf('%s/Outputs/R_A%d_%s.mat', nameFolder_server, iModelA_fit, str_part);
nameFolder_Data_OOD = sprintf('%s/Data_OOD_%d%d', nameFolder_Data, nORI, nSF);
nameFolder_Data_NOM_Trialwise = sprintf('%s/Data_NOM_Trialwise_%d%d', nameFolder_Data, nORI, nSF);

namesMetrics_behav = {'pYES','pC','pA'};

% find IO folders
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

        try
            data_compIV = load(fullfile(nameDir_compIV(1).folder, nameDir_compIV(1).name));
        catch
            nameDir_fitNOM_chk = dir(fullfile(info.folder_NOM, sprintf('n*_A%dB*.mat', iModelA_fit)));
            nameDir_fitNOM_chk = nameDir_fitNOM_chk(~contains({nameDir_fitNOM_chk.name}, {'min','compIV'}));
            if ~isempty(nameDir_fitNOM_chk)
                fprintf('corrupted compIV file, with fitNOM files present: %s\n', info.folder_NOM);
            else
                fprintf('corrupted compIV file, but NO fitNOM files: %s\n', info.folder_NOM);
            end
            continue
        end

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

        % =========================================================
        % 1. template recovery
        % =========================================================
        [template_rmse_iter, template_R2_iter] = fxn_templateRecovery_simPlot2Style_iter(  truth, data_compIV, nIter);

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
        vals = data_compIV.basisWidthScaleORI_tmpl_allIter(:);
        strVec = string(vals);
        [uniqueVals, ~, groupIdx] = unique(strVec);
        counts = accumarray(groupIdx, 1);
        [~, idxMax] = max(counts);
        fileBase.basisWidthScaleORI_mode     = uniqueVals(idxMax);
        fileBase.basisWidthScaleORI_mode_pct = counts(idxMax) / numel(strVec);

        % Modal SF width-scale setting across iterations.
        vals = data_compIV.basisWidthScaleSF_tmpl_allIter(:);
        strVec = string(vals);
        [uniqueVals, ~, groupIdx] = unique(strVec);
        counts = accumarray(groupIdx, 1);
        [~, idxMax] = max(counts);
        fileBase.basisWidthScaleSF_mode     = uniqueVals(idxMax);
        fileBase.basisWidthScaleSF_mode_pct = counts(idxMax) / numel(strVec);

        % Modal SF asymmetry setting across iterations.
        vals = data_compIV.asymSF_rightLeftRatio_tmpl_allIter(:);
        strVec = string(vals);
        [uniqueVals, ~, groupIdx] = unique(strVec);
        counts = accumarray(groupIdx, 1);
        [~, idxMax] = max(counts);
        fileBase.asymSF_rightLeftRatio_mode     = uniqueVals(idxMax);
        fileBase.asymSF_rightLeftRatio_mode_pct = counts(idxMax) / numel(strVec);

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
        fileBase.criterion_DV_true = truth.criterion_DV_true;

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
            Ri.nLL_test_allIter = data_fitNOM.nLL_test_allIter(:)'; % store full iter vector

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

            parfor iIter = 1:nIter_fit
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
    save(nameFile_R, 'R')

else
    % fprintf('Output file already exists, skipping compilation: %s\n', nameFile_R);
    % Load output file
    load(nameFile_R, 'R')
    fprintf('\n\n%s: Loading compiled record R with %d entries.\n\n', datetime('now'), numel(R));

end % if ~exist(nameFile_R)

% ---------- unique levels ----------
gaborCST_unik = unique([R.gaborCST]);
cSDT_unik = unique([R.cSDT_true]);
Nmul_unik = unique([R.Nmul_true]);
Nadd_unik = unique([R.Nadd_true]);
Nshared_unik = unique([R.Nshared_true]);
Bsim_unik = unique([R.iModelB_sim]);
Bfit_unik = unique([R.iModelB_fit]);

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
setting.fig4_axisBufferProp = 0.45;
setting.fig4_markerSizeMin = 7;
setting.fig4_markerSizeMax = 20;
setting.fig4_markerSizeExp = 0.5;
setting.fig2_top_axisBufferProp = 0.12;
setting.fig3_gapColsLeft = 0;           % number of blank columns between heatmap and curves
setting.fig3_curveTileSpacing = 'compact'; % spacing among curve panels: 'none'|'compact'|'loose'

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
setting.fig1_mode_fields = {'nBasisORI_mode', 'nBasisSF_mode', 'basisFxnORI_mode', 'basisFxnSF_mode',  'basisWidthScaleORI_mode', 'basisWidthScaleSF_mode', 'ridge_mode', 'asymSF_rightLeftRatio_mode'};
setting.fig1_mode_labels = {'ORI basis mode', 'SF basis mode', 'ORI basis family mode', 'SF basis family mode',  'ORI width-scale mode', 'SF width-scale mode', 'Ridge mode', 'SF asymmetry mode'};

setting.fig1_pct_fields = {'nBasisORI_mode_pct', 'nBasisSF_mode_pct', 'basisFxnORI_mode_pct', 'basisFxnSF_mode_pct',  'basisWidthScaleORI_mode_pct', 'basisWidthScaleSF_mode_pct', 'ridge_mode_pct', 'asymSF_rightLeftRatio_mode_pct'};
setting.fig1_pct_labels = {'ORI mode selection rate', 'SF mode selection rate', 'ORI basis family mode selection rate', 'SF basis family mode selection rate',  'ORI width-scale mode selection rate', 'SF width-scale mode selection rate', 'Ridge mode selection rate', 'SF asymmetry mode selection rate'};

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

%% Figure 0: condition composition overview from compiled R
% Keep one row per simulated condition to avoid duplication across fitted B models.
fprintf('\n%s: Fig 0: Condition composition...', string(datetime('now')))

T0 = struct2table(R);
keyVars_fig0 = {'nameIO', 'noiseCST', 'gaborCST', 'nTrials',  'Nmul_true', 'Nadd_true', 'Nshared_true', 'cSDT_true', 'iModelB_sim'};
keyVars_fig0 = keyVars_fig0(ismember(keyVars_fig0, T0.Properties.VariableNames));
if ~isempty(keyVars_fig0)
    [~, ia_fig0] = unique(T0(:, keyVars_fig0), 'rows', 'stable');
    T0 = T0(ia_fig0, :);
end

h = figure('Position', [100 100 1600 720]);
tiledlayout(2, 3, 'Padding', 'loose', 'TileSpacing', 'loose');

fig0Fields = {'iModelB_sim', 'gaborCST', 'cSDT_true', 'Nmul_true', 'Nadd_true', 'Nshared_true'};
fig0Labels = {'Generating model', 'signal contrast', 'SDT criterion', 'multiplicative variability', 'additive variability', 'shared variability'};

% B-model shading: B1 = black (0.0), B7 = light grey (0.78)
Bsim_all = unique(T0.iModelB_sim);
Bsim_all = Bsim_all(isfinite(Bsim_all));
nBsim_all = numel(Bsim_all);
Bsim_shades = [0, linspace(0.3, 0.78, max(nBsim_all-1, 2))];  % one shade per B level
Bsim_shades = Bsim_shades(1:nBsim_all);
bsim_color = @(b) repmat(Bsim_shades(Bsim_all == b), 1, 3);  % grey RGB for a given B index

for iPanel = 1:numel(fig0Fields)
    nexttile; hold on;

    vals = T0.(fig0Fields{iPanel});
    vals = vals(:);
    vals = vals(isfinite(vals));

    [uVals, ~, iU] = unique(vals);
    counts = accumarray(iU, 1);
    nTotal_panel = sum(counts);
    counts_prop = counts ./ max(nTotal_panel, 1);

    xPos = 1:numel(uVals);

    if strcmp(fig0Fields{iPanel}, 'iModelB_sim')
        % Bsim panel: one bar per B level, each with its own shade
        b = bar(xPos, counts_prop, 0.75);
        b.FaceColor = 'flat';
        b.EdgeColor = 'k';
        b.LineWidth = 0.5;
        for ib = 1:numel(uVals)
            b.CData(ib, :) = bsim_color(uVals(ib));
        end
    else
        % All other panels: stacked bars segmented by Bsim contribution
        % countMat: rows = x levels, cols = Bsim levels
        countMat = zeros(numel(uVals), nBsim_all);
        for ib = 1:nBsim_all
            mask_b = T0.iModelB_sim == Bsim_all(ib);
            vals_b = T0.(fig0Fields{iPanel})(mask_b);
            vals_b = vals_b(isfinite(vals_b));
            for ix = 1:numel(uVals)
                countMat(ix, ib) = sum(vals_b == uVals(ix));
            end
        end
        countMat_prop = countMat ./ max(nTotal_panel, 1);
        % Draw stacked bars and color each series by Bsim shade
        b = bar(xPos, countMat_prop, 0.75, 'stacked');
        for ib = 1:nBsim_all
            b(ib).FaceColor = bsim_color(Bsim_all(ib));
            b(ib).EdgeColor = 'w';
            b(ib).LineWidth = 0.75;
        end
    end

    % Add count labels above each bar segment
    yTop = counts_prop;
    yPad = max(0.015, 0.05 * max(yTop));
    for ix = 1:numel(xPos)
        text(xPos(ix), yTop(ix) + yPad, sprintf('%d', counts(ix)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 10);
    end

    % Set y-axis limit slightly above the tallest bar
    ymax_panel = max(yTop) + 3 * yPad;
    if ~isfinite(ymax_panel) || ymax_panel <= 0
        ymax_panel = 1;
    end
    ylim([0, min(1.05, ymax_panel)]);

    ylabel('Proportion of conditions');
    xticks(xPos);
    if strcmp(fig0Fields{iPanel}, 'iModelB_sim')
        xticklabels(compose('M%d', round(uVals)));
        xlabel(fig0Labels{iPanel});
    else
        xticklabels(compose('%.3g', uVals));
        xlabel(sprintf('True %s', fig0Labels{iPanel}));
    end

    box off;

    % Add legend only to the Bsim panel
    if iPanel == 1
        hLeg = gobjects(nBsim_all, 1);
        for ib = 1:nBsim_all
            clr = bsim_color(Bsim_all(ib));
            hLeg(ib) = patch(nan, nan, clr,  'EdgeColor', 'k',  'LineWidth', 0.75,  'DisplayName', sprintf('M%d %s', Bsim_all(ib), namesModelB{ib}));
        end
        legend(hLeg, 'Location', 'northeast',  'Box', 'off',  'FontSize', 10);
    end

end % for iPanel

sgtitle('Figure 0: Number of conditions by simulation variable');
set(findall(gcf,'-property','FontSize'), 'FontSize', 15);

saveas(h, fullfile(nameFolder_Figures_part4, 'FigS0_conditionComposition.png'));
close(h);
fprintf('DONE\n\n')

%% Figure 1: Basis selection
fprintf('\n%s: Fig 1: Basis selection rates...', string(datetime('now')))

R_fig1 = R([R.iModelB_sim] == setting.fig12_Bsim & [R.iModelB_fit] == setting.fig12_Bfit);

grpLabels_fig1 = arrayfun(@(r) sprintf('sig=%.3g', r.gaborCST),  R_fig1, 'UniformOutput', false);

% grpLabels_fig1 = arrayfun(@(r) sprintf('Nmul=%g x Nadd=%g x Nshared=%g', r.Nmul_true, r.Nadd_true, r.Nshared_true),  %     R_fig1, 'UniformOutput', false);

uGrp_fig1    = unique(grpLabels_fig1);
grpCmap_fig1 = setting.palette(1:numel(uGrp_fig1), :);

h = figure('Position', [100 100 2e3 1e3]);

panel_order = [1, 5, 2, 6, 3, 7, 4, 8];
panel_xlabel = {'# ORI basis functions', '# SF basis functions',  'ORI basis family', 'SF basis family',  'ORI width scale', 'SF width scale',  'Ridge penalty', 'asymSF_rightLeftRatio'};
panel_title = {'nBasisORI', 'nBasisSF',  'basisFxnORI', 'basisFxnSF',  'basisWidthScaleORI', 'basisWidthScaleSF',  'ridge', 'asymSF_rightLeftRatio'};

for iPanel = 1:numel(setting.fig1_mode_fields)
    subplot(2,4,panel_order(iPanel)); hold on;
    plot_ranked_categorical({R_fig1.(setting.fig1_mode_fields{iPanel})}, panel_title{iPanel}, grpLabels_fig1, grpCmap_fig1, iPanel == 1);
    ylabel('% iterations selected');
    xlabel(panel_xlabel{iPanel});
    title(panel_title{iPanel});
    fxn_style_ax(gca, setting);
end

sgtitle('Fig 1: rank of selected basis settings across iterations', 'FontWeight', 'bold');
% Set all font size in this figure
set(findall(gcf,'-property','FontSize'), 'FontSize', 12);

saveas(h, fullfile(nameFolder_Figures_part4, 'FigS1_selectedBasisRank.png'));
close(h);
fprintf('DONE\n\n')

%% Figure 2: Template recovery
%    Use one representative Bsim/Bfit because these are file-level

fprintf('\n%s: Fig 2: Template recovery...', string(datetime('now')))

R_fig2 = R([R.iModelB_sim] == setting.fig12_Bsim & [R.iModelB_fit] == setting.fig12_Bfit);

nVars = numel(setting.fig2_varFields);
setting_fig2 = setting;

fig2Cache = struct([]);
for iCol = 1:nVars
    vField  = setting.fig2_varFields{iCol};
    vName   = setting.fig2_varNames{iCol};
    vLevels = unique([R_fig2.(vField)]);

    fig2Cache(iCol).vField = vField; %#ok<AGROW>
    fig2Cache(iCol).vName = vName;
    fig2Cache(iCol).vLevels = vLevels;
    fig2Cache(iCol).rmse_all = fxn_collapse_scalar_by_x_iterCI(  R_fig2, vField, 'template_rmse_iter', vLevels, setting.scalarCollapseFcn);

    [xAxisORI, yTrueORI, curvesORI] = fxn_collapse_tuning_panel(R_fig2, vField, vLevels, 'ORI');
    [xAxisSF,  yTrueSF,  curvesSF]  = fxn_collapse_tuning_panel(R_fig2, vField, vLevels, 'SF');
    fig2Cache(iCol).tuning(1).xAxis = xAxisORI;
    fig2Cache(iCol).tuning(1).yTrue = yTrueORI;
    fig2Cache(iCol).tuning(1).curves = curvesORI;
    fig2Cache(iCol).tuning(2).xAxis = xAxisSF;
    fig2Cache(iCol).tuning(2).yTrue = yTrueSF;
    fig2Cache(iCol).tuning(2).curves = curvesSF;
end % for iCol

h = figure('Position', [100 100 2e3 1.2e3]);
tiledlayout(3, nVars, 'TileSpacing', 'loose', 'Padding', 'loose');
ax_rmse_last = [];

for iCol = 1:nVars
    vName = fig2Cache(iCol).vName;

    % ------ Row 1: RMSE ------
    ax_rmse = nexttile(iCol); hold on;
    ax_rmse_last = ax_rmse;

    % Always overlay grand mean.
    S = fig2Cache(iCol).rmse_all;
    if ~isempty(S.x)
        xTop = S.x;
        yTop = S.y;
        semTop = S.sem;

        if ~isempty(xTop)
            % Plot the connecting line
            plot(xTop, yTop, '-', 'Color', [0 0 0], 'LineWidth', setting.trueWidth);

            nLevTop = numel(xTop);
            if nLevTop >= 2
                grayValsTop = log10(logspace(0.55, 0.00, nLevTop)');
            else
                grayValsTop = 0.55;
            end

            % Plot dots
            for iData = 1:nLevTop
                cPt = repmat(grayValsTop(iData), 1, 3);
                plot(xTop(iData), yTop(iData), 'o', ...
                    'MarkerSize', 20, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', cPt, 'LineStyle', 'none', 'LineWidth', 3);
            end % iData

            % Plot errorbars
            errorbar(xTop, yTop, semTop, 'LineStyle', 'none', 'Color', [0 0 0], 'LineWidth', 3, 'CapSize', 0);

        end
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
    xlabel(vName);
    ylim([0.1 0.14]);
    yticks(0.1:0.01:0.14);
    title(vName);
    if iCol == 1
        ylabel('2D Template RMSE');
    end

    % ------ Rows 2 & 3: ORI and SF ------
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
                legStr{nLeg} = sprintf('%g | %.3f', curves(iLev).level, rmse);
            end
            legend(legH(1:nLeg), legStr(1:nLeg), 'Location', 'south', 'FontSize', setting_fig2.fontSize, 'Box', 'off');
        end

        fxn_style_ax(gca, setting_fig2);
        xticks(axisTicks_tuning{iRow}); xticklabels(axisTL_tuning{iRow});
        ylim([-.2, 1]);
        % title(vName);
        xlabel(namesFeature{iRow});
        if iCol == 1
            ylabel('Marg. weight (divided by max)');
        end
    end
end

sgtitle('Fig 2: Template RMSE + tuning recovery', 'FontWeight', 'bold', 'FontSize', setting_fig2.fontSize);
set(findall(gcf,'-property','FontSize'), 'FontSize', 20);

saveas(h, fullfile(nameFolder_Figures_part4, 'FigS2_TemplateRecovery.png'));
close(h);
fprintf('DONE\n\n')

%% Figure 3: Metric recovery
%    Leftmost column: heatmap of RMSE (fit × sim); remaining columns: curves per generating model.

fprintf('\n%s: Fig 3: Metric recovery...', string(datetime('now')))

nMetric_fig3 = numel(setting.metricFields);
nSim_fig3 = numel(Bsim_unik);
nFit_fig3 = numel(Bfit_unik);

% Precompute mean/SEM of RMSE across conditions for each (metric, sim, fit) triple.
rmseByFit_mean_fig3 = nan(nMetric_fig3, nSim_fig3, nFit_fig3);
rmseByFit_sem_fig3  = nan(nMetric_fig3, nSim_fig3, nFit_fig3);

for iMetric_fig3 = 1:nMetric_fig3
    mName_fig3 = setting.metricFields{iMetric_fig3};
    rmseIterField_fig3 = sprintf('%s_rmse_allIter', mName_fig3);

    for iSim_fig3 = 1:nSim_fig3
        simModel_fig3 = Bsim_unik(iSim_fig3);
        for iFit_fig3 = 1:nFit_fig3
            fitModel_fig3 = Bfit_unik(iFit_fig3);
            Rsub_fig3 = R([R.iModelB_sim] == simModel_fig3 & [R.iModelB_fit] == fitModel_fig3);

            medians_thisSim = [];
            for iRow_fig3 = 1:numel(Rsub_fig3)
                vals_fig3 = Rsub_fig3(iRow_fig3).(rmseIterField_fig3);
                vals_fig3 = vals_fig3(:);
                [med_iter, ~, ~] = getCI(vals_fig3, 1, 1);
                medians_thisSim = [medians_thisSim; med_iter]; %#ok<AGROW>
            end

            [rmseByFit_mean_fig3(iMetric_fig3, iSim_fig3, iFit_fig3), ~, ~, ...
                rmseByFit_sem_fig3(iMetric_fig3, iSim_fig3, iFit_fig3)] = getCI(medians_thisSim, 2, 1);
        end
    end
end

% Layout: nMetric rows × (1 heatmap + nGapCols + nSim curve columns)
h = figure('Position', [0 0 2e3 7e2]);
nGapCols_fig3 = max(0, round(setting.fig3_gapColsLeft));
tlo_fig3 = tiledlayout(nMetric_fig3, 1 + nGapCols_fig3 + nSim_fig3, ...
    'TileSpacing', setting.fig3_curveTileSpacing, 'Padding', 'compact');

% Axis tick labels for heatmap (same logic as Fig 5)
xTickLabels_fig3 = arrayfun(@(k) sprintf('M%d', k), Bsim_unik, 'UniformOutput', false);
yTickLabels_fig3 = arrayfun(@(k) sprintf('M%d', k), Bfit_unik, 'UniformOutput', false);
% if exist('namesModelB', 'var') && ~isempty(namesModelB)
%     for iLab = 1:numel(Bsim_unik)
%         k = Bsim_unik(iLab);
%         if k >= 1 && k <= numel(namesModelB)
%             xTickLabels_fig3{iLab} = char(string(namesModelB{k}));
%         end
%     end
%     for iLab = 1:numel(Bfit_unik)
%         k = Bfit_unik(iLab);
%         if k >= 1 && k <= numel(namesModelB)
%             yTickLabels_fig3{iLab} = char(string(namesModelB{k}));
%         end
%     end
% end

for iRow_metric = 1:nMetric_fig3

    mName    = setting.metricFields{iRow_metric};
    dataField = setting.metricCurveDataFields{iRow_metric};
    predField = setting.metricCurvePredFields{iRow_metric};

    switch mName
        case 'pYES';  metricTitle = 'Detection rate';
        case 'pC';    metricTitle = 'Accuracy';
        case 'pA';    metricTitle = 'Consistency rate';
        otherwise;    metricTitle = mName;
    end

    % ------------------------------------------------------------------ %
    % Leftmost column: heatmap of mean RMSE  (rows = fit model, cols = sim) %
    % ------------------------------------------------------------------ %
    M_heatmap_fig3 = squeeze(rmseByFit_mean_fig3(iRow_metric, :, :))'; % nFit × nSim

    axHeat = nexttile(tlo_fig3);
    hold(axHeat, 'on');

    finVals_h = M_heatmap_fig3(isfinite(M_heatmap_fig3));
    if ~isempty(finVals_h)
        q_h = prctile(finVals_h, [2.5 97.5]);
        cLo_h = q_h(1); cHi_h = q_h(2);
        if ~isfinite(cLo_h) || ~isfinite(cHi_h) || cLo_h >= cHi_h
            cLo_h = min(finVals_h); cHi_h = max(finVals_h);
        end
        if cLo_h == cHi_h; cLo_h = cLo_h - eps; cHi_h = cHi_h + eps; end
    else
        cLo_h = 0; cHi_h = 1;
    end
    M_heatmap_disp = min(max(M_heatmap_fig3, cLo_h), cHi_h);

    imagesc(axHeat, 1:nSim_fig3, 1:nFit_fig3, M_heatmap_disp);
    set(axHeat, 'YDir', 'normal', ...
        'XTick', 1:nSim_fig3, 'XTickLabel', xTickLabels_fig3, ...
        'YTick', 1:nFit_fig3, 'YTickLabel', yTickLabels_fig3, ...
        'XLim', [0.5, nSim_fig3+0.5], 'YLim', [0.5, nFit_fig3+0.5]);
    xlabel(axHeat, 'Generating model');
    ylabel(axHeat, 'Fitting model');
    title(axHeat, [metricTitle ' RMSE']);
    axis(axHeat, 'square');
    clim(axHeat, [cLo_h cHi_h]);
    cbar_h = colorbar(axHeat, 'Location', 'eastoutside');
    fxn_set_colorbar_ticks(cbar_h, cLo_h, cHi_h);

    % % Annotate heatmap cell values
    % for iFit_h = 1:nFit_fig3
    %     for iSim_h = 1:nSim_fig3
    %         v_h = M_heatmap_fig3(iFit_h, iSim_h);
    %         if ~isfinite(v_h), continue; end
    %         frac_h = max(0, min(1, (v_h - cLo_h) / (cHi_h - cLo_h)));
    %         if frac_h < 0.5; txtClr_h = [1 1 1]; else; txtClr_h = [0 0 0]; end
    %         text(axHeat, iSim_h, iFit_h, sprintf('%.3f', v_h), ...
    %             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    %             'FontSize', 8, 'Color', txtClr_h);
    %     end
    % end

    % Matching models: red solid boxes on the diagonal
    nDiag_h = min(nFit_fig3, nSim_fig3);
    for iDiag_h = 1:nDiag_h
        rectangle('Parent', axHeat, 'Position', [iDiag_h-0.5, iDiag_h-0.5, 1, 1], 'EdgeColor', 'r', 'LineWidth', 2.5, 'LineStyle', '-');
    end

    % Lowest RMSE per column: Black dashed boxes
    for iSim_h = 1:nSim_fig3
        colVals_h = M_heatmap_fig3(:, iSim_h);
        finMask_h = isfinite(colVals_h);
        if ~any(finMask_h), continue; end
        finIdx_h = find(finMask_h);
        [~, iLoc_h] = min(colVals_h(finMask_h));
        bestRow_h = finIdx_h(iLoc_h);
        rectangle('Parent', axHeat, 'Position', [iSim_h-0.5, bestRow_h-0.5, 1, 1], 'EdgeColor', 'k', 'LineWidth', 2.5, 'LineStyle', '--');
    end

    % Spacer tiles to separate the heatmap column from curve columns.
    for iGap_fig3 = 1:nGapCols_fig3
        axGap_fig3 = nexttile(tlo_fig3);
        axis(axGap_fig3, 'off');
    end

    % ------------------------------------------------------------------ %
    % RIGHT: Curve panels: one per generating model                              %
    % ------------------------------------------------------------------ %
    for iCol_Bsim = 1:nSim_fig3
        simModel = Bsim_unik(iCol_Bsim);

        R_data = R([R.iModelB_sim] == simModel & [R.iModelB_fit] == 1);

        axMain = nexttile(tlo_fig3);
        hold(axMain, 'on');
        isFirstPanel = (iCol_Bsim == 1) && (iRow_metric == 1);

        [Cdata, nDataBins] = fxn_collapse_curve_with_counts(R_data, 'DV_allBins_med', dataField, 'nTrials_allBins_med', setting.curveNGrid, @median);
        dot_nMin = nan; dot_nMax = nan;
        dot_sMin = max(4, setting.markerSize * 0.5);
        dot_sMax = max(dot_sMin + 2, setting.markerSize * 1.5);

        % First pass: collect predicted curves and compute RMSE for all fit models
        nFit = numel(Bfit_unik);
        CpredAll_fig3 = cell(1, nFit);
        rmseAll_fig3  = nan(1, nFit);
        for iFit = 1:nFit
            fitModel = Bfit_unik(iFit);
            R_pred = R([R.iModelB_sim] == simModel & [R.iModelB_fit] == fitModel);
            CpredAll_fig3{iFit} = fxn_collapse_curve(R_pred, 'DV_allBins_med', predField, setting.curveNGrid, @median);
            if ~isempty(CpredAll_fig3{iFit}.x) && ~isempty(Cdata.x)
                rmseAll_fig3(iFit) = fxn_curve_rmse(Cdata.x, Cdata.y, CpredAll_fig3{iFit}.x, CpredAll_fig3{iFit}.y, nDataBins);
            end
        end

        % Identify matching and best-RMSE fit models
        matchFitIdx_fig3 = find(Bfit_unik == simModel, 1);
        finRmseIdx_fig3  = find(isfinite(rmseAll_fig3));
        if ~isempty(finRmseIdx_fig3)
            [~, iLoc_best] = min(rmseAll_fig3(finRmseIdx_fig3));
            bestFitIdx_fig3 = finRmseIdx_fig3(iLoc_best);
        else
            bestFitIdx_fig3 = [];
        end

        % Second pass: grey curves (skip matching and best-non-matching)
        for iFit = 1:nFit
            if iFit == matchFitIdx_fig3, continue; end
            if ~isempty(bestFitIdx_fig3) && iFit == bestFitIdx_fig3 && bestFitIdx_fig3 ~= matchFitIdx_fig3
                continue;
            end
            Cpred = CpredAll_fig3{iFit};
            if isempty(Cpred.x), continue; end
            plot(Cpred.x, Cpred.y, 'Color', [.5, .5, .5], 'LineStyle', '-', 'LineWidth', setting.lineWidth);
        end

        % Plot matching model: red solid
        if ~isempty(matchFitIdx_fig3) && ~isempty(CpredAll_fig3{matchFitIdx_fig3}.x)
            plot(CpredAll_fig3{matchFitIdx_fig3}.x, CpredAll_fig3{matchFitIdx_fig3}.y, 'Color', [1, 0, 0], 'LineStyle', '-', 'LineWidth', setting.lineWidth*1.5);
        end

        % Plot best-RMSE model last: black dashed; label if not the matching model
        if ~isempty(bestFitIdx_fig3)
            Cbest = CpredAll_fig3{bestFitIdx_fig3};
            if ~isempty(Cbest.x)
                % plot the curve
                plot(Cbest.x, Cbest.y, 'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', setting.lineWidth*1.5);
                % Print the M
                % if bestFitIdx_fig3 ~= matchFitIdx_fig3
                %     bestModel = Bfit_unik(bestFitIdx_fig3);
                %     if exist('namesModelB', 'var') && ~isempty(namesModelB) && bestModel >= 1 && bestModel <= numel(namesModelB)
                %         bestLabel_fig3 = char(string(namesModelB{bestModel}));
                %     else
                %         bestLabel_fig3 = sprintf('M%d', bestModel);
                %     end
                %     [~, iRight] = max(Cbest.x);
                %     text(Cbest.x(iRight), Cbest.y(iRight), bestLabel_fig3, ...
                %         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                %         'FontSize', max(setting.fontSize - 2, 7), 'Color', [0 0 0]);
                % end
            end
        end

        % Plot data dots
        if ~isempty(Cdata.x)
            mSizes = setting.markerSize * ones(size(Cdata.x));
            if ~isempty(nDataBins) && numel(nDataBins) == numel(Cdata.x)
                nDot = nDataBins(:)';
                if any(isfinite(nDot))
                    nMin = min(nDot); nMax = max(nDot);
                    dot_nMin = nMin; dot_nMax = nMax;
                    if nMax > nMin
                        frac = (nDot - nMin) / (nMax - nMin);
                    else
                        frac = zeros(size(nDot));
                    end
                    mSizes = dot_sMin + (dot_sMax - dot_sMin) * (frac .^ setting.fig4_markerSizeExp);
                end
            end
            for iDot = 1:numel(Cdata.x)
                plot(Cdata.x(iDot), Cdata.y(iDot), 'ko', 'MarkerSize', mSizes(iDot), 'MarkerFaceColor', 'w', 'LineWidth', setting.lineWidth);
            end
        end

        fxn_style_ax(gca, setting);
        set(gca, 'YGrid', 'off');

        if iRow_metric == nMetric_fig3
            xlabel('Collapsed DV bins');
        end
        if iCol_Bsim == 1
            ylabel(metricTitle);
        end
        if iRow_metric == 1
            title(sprintf('M%d %s', iCol_Bsim, namesModelB{iCol_Bsim}));
        else
            title('');
        end
        axis square

        % Define ylimits
        switch mName
            case 'pYES';  ylim([0 1]); yticks(0:.2:1)
            case {'pC','pA'};  ylim([0.5 1]); yticks(0.5:.1:1)
        end

        % Add legend for dot size
        if isFirstPanel && ~isempty(Cdata.x)
            axPos = get(axMain, 'Position');
            axDot = axes('Parent', h, 'Position', axPos, 'Color', 'none', 'Visible', 'off', ...
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
            hDot1 = scatter(axDot, nan, nan, dotSizesRef(1)^2, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', setting.lineWidth);
            hDot2 = scatter(axDot, nan, nan, dotSizesRef(2)^2, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', setting.lineWidth);
            hDot3 = scatter(axDot, nan, nan, dotSizesRef(3)^2, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', setting.lineWidth);
            legDot = legend(axDot, [hDot1 hDot2 hDot3], {'10', '100', '1000'}, 'Location', 'northwest', 'FontSize', setting.fontSize-2, 'Box', 'off');
            title(legDot, 'Dot size');

        end
    end % iCol_Bsim
end % iRow_metric

sgtitle('Fig 3: metric recovery', 'FontWeight', 'bold');
% set(findall(gcf,'-property','FontSize'), 'FontSize', 5);

saveas(h, fullfile(nameFolder_Figures_part4, 'FigS3_metricRecovery.png'));
close(h);
fprintf('DONE\n\n')

%% Figure 4: Parameter recovery
%    Only matched pairs Bsim = Bfit

fprintf('\n%s: Fig 4: Parameter recovery...', string(datetime('now')))

setting_fig4 = setting;
setting_fig4.fontSize = setting.fontSize;
Pplot.markerSize = 50;
Pplot.lineWidth=1.5;
Pplot.errLineWidth = 1.5;

% Condition coding in Fig 4:
%   marker color hue                = cSDT_true
%   marker shading/lightness        = signalCST low→high

nCz_fig4  = numel(cSDT_unik);

% Grayscale ramp for signalCST: index 1 = lightest (0.78), last = black (0)
nSig_fig4 = numel(gaborCST_unik);

% Match Figure 6 color construction: hue by cSDT, shading by signalCST.
czMin_fig4 = min(cSDT_unik);
czMax_fig4 = max(cSDT_unik);
sigMin_fig4 = min(gaborCST_unik);
sigMax_fig4 = max(gaborCST_unik);

% Precompute RMSE heatmaps per parameter across (fit model x sim model).
nParam_fig4 = numel(setting.fig4_paramNames);
nFit_fig4 = numel(Bfit_unik);
nSim_fig4 = numel(Bsim_unik);
M_paramRMSE_fig4 = nan(nParam_fig4, nFit_fig4, nSim_fig4);
for iParam_fig4 = 1:nParam_fig4
    pName_fig4 = setting.fig4_paramNames{iParam_fig4};
    for iSim_fig4 = 1:nSim_fig4
        simModel_fig4 = Bsim_unik(iSim_fig4);
        for iFit_fig4 = 1:nFit_fig4
            fitModel_fig4 = Bfit_unik(iFit_fig4);

            % Skip parameters absent in reduced models.
            if strcmp(pName_fig4, 'Nshared') && ~any(strcmp(namesModelBparams_short{fitModel_fig4}, 'Nshared'))
                continue
            end
            if strcmp(pName_fig4, 'Nmul') && ~any(strcmp(namesModelBparams_short{fitModel_fig4}, 'Nmul'))
                continue
            end
            if strcmp(pName_fig4, 'Nadd') && ~any(strcmp(namesModelBparams_short{fitModel_fig4}, 'Nadd'))
                continue
            end

            Rsub_fig4hm = R([R.iModelB_sim] == simModel_fig4 & [R.iModelB_fit] == fitModel_fig4);
            if isempty(Rsub_fig4hm)
                continue
            end
            rmseVals_fig4hm = [Rsub_fig4hm.(sprintf('%s_rmse', pName_fig4))]';
            [M_paramRMSE_fig4(iParam_fig4, iFit_fig4, iSim_fig4), ~, ~] = getCI(rmseVals_fig4hm, 2, 1);
        end
    end
end

h = figure('Position', [0 0 2e3 9e2]);
tiledlayout(numel(setting.fig4_paramNames), 1 + numel(Bfit_unik), 'TileSpacing', 'compact', 'Padding', 'compact');
ax_fig4_last = gobjects(1);

xTickLabels_fig4 = arrayfun(@(k) sprintf('M%d', k), Bsim_unik, 'UniformOutput', false);
yTickLabels_fig4 = arrayfun(@(k) sprintf('M%d', k), Bfit_unik, 'UniformOutput', false);

for iRow_param = 1:numel(setting.fig4_paramNames)
    pName = setting.fig4_paramNames{iRow_param};
    pTitle = setting.fig4_paramTitles{iRow_param};

    % ------------------------------------------------------------------ %
    % Leftmost column: heatmap of mean RMSE
    % ------------------------------------------------------------------ %
    M_heatmap_fig4 = squeeze(M_paramRMSE_fig4(iRow_param, :, :));
    axHeat_fig4 = nexttile;
    hold(axHeat_fig4, 'on');

    finVals_fig4 = M_heatmap_fig4(isfinite(M_heatmap_fig4));
    if ~isempty(finVals_fig4)
        q_fig4 = prctile(finVals_fig4, [2.5 97.5]);
        cLo_fig4 = q_fig4(1);
        cHi_fig4 = q_fig4(2);
        if ~isfinite(cLo_fig4) || ~isfinite(cHi_fig4) || cLo_fig4 >= cHi_fig4
            cLo_fig4 = min(finVals_fig4);
            cHi_fig4 = max(finVals_fig4);
        end
        if cLo_fig4 == cHi_fig4
            cLo_fig4 = cLo_fig4 - eps;
            cHi_fig4 = cHi_fig4 + eps;
        end
    else
        cLo_fig4 = 0;
        cHi_fig4 = 1;
    end
    M_heatmap_disp_fig4 = min(max(M_heatmap_fig4, cLo_fig4), cHi_fig4);

    imagesc(axHeat_fig4, 1:nSim_fig4, 1:nFit_fig4, M_heatmap_disp_fig4);
    set(axHeat_fig4, 'YDir', 'normal', ...
        'XTick', 1:nSim_fig4, 'XTickLabel', xTickLabels_fig4, ...
        'YTick', 1:nFit_fig4, 'YTickLabel', yTickLabels_fig4, ...
        'XLim', [0.5, nSim_fig4+0.5], 'YLim', [0.5, nFit_fig4+0.5]);
    xlabel(axHeat_fig4, 'Generating model');
    ylabel(axHeat_fig4, 'Fitting model');
    title(axHeat_fig4, sprintf('%s RMSE', pTitle));
    axis(axHeat_fig4, 'square');
    clim(axHeat_fig4, [cLo_fig4 cHi_fig4]);
    cbar_fig4 = colorbar(axHeat_fig4, 'Location', 'eastoutside');
    fxn_set_colorbar_ticks(cbar_fig4, cLo_fig4, cHi_fig4);

    % Red box on diagonal (matched models)
    nDiag_fig4 = min(nFit_fig4, nSim_fig4);
    for iDiag_fig4 = 1:nDiag_fig4
        rectangle('Parent', axHeat_fig4, 'Position', [iDiag_fig4-0.5, iDiag_fig4-0.5, 1, 1], ...
            'EdgeColor', 'r', 'LineWidth', 2.5, 'LineStyle', '-');
    end

    % Black dashed box on lowest RMSE per sim-model column
    for iSim_hm = 1:nSim_fig4
        colVals_hm = M_heatmap_fig4(:, iSim_hm);
        finMask_hm = isfinite(colVals_hm);
        if ~any(finMask_hm), continue; end
        finIdx_hm = find(finMask_hm);
        [~, iLoc_hm] = min(colVals_hm(finMask_hm));
        bestRow_hm = finIdx_hm(iLoc_hm);
        rectangle('Parent', axHeat_fig4, 'Position', [iSim_hm-0.5, bestRow_hm-0.5, 1, 1], ...
            'EdgeColor', 'k', 'LineWidth', 2.5, 'LineStyle', '--');
    end

    % ------------------------------------------------------------------ %
    % Right columns: scatter plots
    % ------------------------------------------------------------------ %
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

        % Define x and y limit to draw unity line
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

        % plot unit line
        plot([mn mx], [mn mx], 'k--', 'LineWidth', 1);

        % Define xlim and ylim
        xlim([mn mx]); ylim([mn mx]);
        if ~isempty(tickVals)
            xticks(tickVals);
            yticks(tickVals);
        end

        % one series per (signalCST x cSDT_true)
        for iSig = 1:numel(gaborCST_unik)
            for iCz = 1:numel(cSDT_unik)
                S = fxn_collapse_param_byLevel(  Rsub, pName,  setting.scalarCollapseFcn, nBins_Part4, gaborCST_unik(iSig), cSDT_unik(iCz));

                if isempty(S.x), continue; end

                % marker hue = cSDT, shading = signalCST (same mapping as Fig 6)
                czVal_fig4 = cSDT_unik(iCz);
                sigVal_fig4 = gaborCST_unik(iSig);

                if czMax_fig4 > czMin_fig4
                    czNorm_fig4 = (czVal_fig4 - czMin_fig4) / (czMax_fig4 - czMin_fig4);
                else
                    czNorm_fig4 = 0.5;
                end
                if sigMax_fig4 > sigMin_fig4
                    sigNorm_fig4 = (sigVal_fig4 - sigMin_fig4) / (sigMax_fig4 - sigMin_fig4);
                else
                    sigNorm_fig4 = 0.5;
                end

                hCol_fig4 = (1 - czNorm_fig4) * (2 / 3);
                sCol_fig4 = 0.25 + 0.70 * sigNorm_fig4;
                vCol_fig4 = 0.98 - 0.13 * sigNorm_fig4;
                faceClr = hsv2rgb([hCol_fig4, sCol_fig4, vCol_fig4]);

                % Plot connecting line
                plot(S.x, S.y, '-',  'Color', faceClr,  'LineWidth', 1);

                for iPt = 1:numel(S.x)

                    % Plot dots
                    scatter(S.x(iPt), S.y(iPt), Pplot.markerSize,  'o', 'MarkerFaceColor', faceClr, 'MarkerEdgeColor', 'w', 'LineWidth', Pplot.lineWidth);

                    % Plot error bars
                    if isfield(S, 'sem') && numel(S.sem) >= iPt && isfinite(S.sem(iPt))
                        errorbar(S.x(iPt), S.y(iPt), S.sem(iPt), 'LineStyle', 'none', 'Color', faceClr, 'LineWidth', Pplot.errLineWidth, 'CapSize', 0);
                    end
                end % for iPt
            end % for iCz
        end % for iSig

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
        Ppanel = fxn_collect_fig4_panel_points(  Rsub, pName, gaborCST_unik, cSDT_unik, setting.scalarCollapseFcn, nBins_Part4);
        [rho_panel, ~, ~, rmse_panel] = fxn_param_recovery_stats(Ppanel.x, Ppanel.y);
        % title(sprintf('M(sim/gen)=%d %s \nr=%.2f | RMSE=%.3f', Bsimfit, namesModelB{Bsimfit}, rho_panel, rmse_panel));
        if iRow_param==1
            title(sprintf('M%d %s', Bsimfit, namesModelB{Bsimfit}));
        end
        ax_fig4_last = gca;
    end % iCol_Bsimfit
end % iRow_param

sgtitle('Fig 4: Parameter recovery | Matched model pairs', 'FontWeight', 'bold');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 13)
saveas(h, fullfile(nameFolder_Figures_part4, 'FigS4_paramRecovery.png'));
close(h);
fprintf('DONE\n\n')

%% Figure 5: Parameter-pair correlation summary by signalCST x cSDT
fprintf('\n%s: Fig 5: Parameter-pair correlation summary...', string(datetime('now')))

% Fixed 6 parameter-pair columns from the canonical set used in Fig 5.
pairParams_fig5 = {'Nmul', 'Nadd', 'Nshared', 'criterion_DV'};
pairLabels_fig5 = {'N_{mul}', 'N_{add}', 'N_{shared}', 'Criterion DV'};
pairIdx_fig5 = nchoosek(1:numel(pairParams_fig5), 2);
pairIdx_fig5 = pairIdx_fig5([1, 2, 4, 3, 5, 6], :); % rearrange so all xx-cDV pairs are last
nCols_fig5 = size(pairIdx_fig5, 1); % expected = 6

% One row per matched model (Bfit = Bsim), expected 7 rows.
Bdiag_fig5 = intersect(Bfit_unik(:)', Bsim_unik(:)');
nRows_fig5 = numel(Bdiag_fig5);

nSig_fig5 = numel(gaborCST_unik);
nCz_fig5 = numel(cSDT_unik);
nBars_fig5 = nSig_fig5 * nCz_fig5;
xBar_fig5 = 1:nBars_fig5;
nPerm_fig5 = 5000;

% For labeling and indexing bars: one per signalCST x cSDT combination
barLabels_fig5 = cell(1, nBars_fig5);
barSigIdx_fig5 = nan(1, nBars_fig5);
barCzIdx_fig5 = nan(1, nBars_fig5);
iBar_fig5 = 0;
for iCz = 1:nCz_fig5
    for iSig = 1:nSig_fig5
        iBar_fig5 = iBar_fig5 + 1;
        barSigIdx_fig5(iBar_fig5) = iSig;
        barCzIdx_fig5(iBar_fig5) = iCz;
        barLabels_fig5{iBar_fig5} = sprintf('G%g|C%g', gaborCST_unik(iSig), cSDT_unik(iCz));
    end
end

% Create a struct to hold all permutation results
permFile_fig5 = fullfile(nameFolder_Figures_part4, 'FigS5_paramCorr_summaryBars_permNull.mat');
cacheLoaded_fig5 = false;

% Precompute bar colors once because they only depend on signalCST and cSDT level.
czMin_fig5 = min(cSDT_unik);
czMax_fig5 = max(cSDT_unik);
sigMin_fig5 = min(gaborCST_unik);
sigMax_fig5 = max(gaborCST_unik);
barColors_fig5 = nan(nBars_fig5, 3);
for iBar = 1:nBars_fig5
    czVal = cSDT_unik(barCzIdx_fig5(iBar));
    sigVal = gaborCST_unik(barSigIdx_fig5(iBar));

    if czMax_fig5 > czMin_fig5
        czNorm = (czVal - czMin_fig5) / (czMax_fig5 - czMin_fig5);
    else
        czNorm = 0.5;
    end

    if sigMax_fig5 > sigMin_fig5
        sigNorm = (sigVal - sigMin_fig5) / (sigMax_fig5 - sigMin_fig5);
    else
        sigNorm = 0.5;
    end

    hCol = (1 - czNorm) * (2 / 3);
    sCol = 0.25 + 0.70 * sigNorm;
    vCol = 0.98 - 0.13 * sigNorm;
    barColors_fig5(iBar, :) = hsv2rgb([hCol, sCol, vCol]);
end

% load the cached permutation results if they exist, otherwise compute and save them
cacheLoaded_fig5 = false;
if exist(permFile_fig5, 'file')
    % fprintf('\n Loading cached permStats_fig5 from %s\n', permFile_fig5);
    load(permFile_fig5, 'permStats_fig5');
    % Validate cache dimensions against current data
    if size(permStats_fig5.observedMean, 3) ~= nBars_fig5
        fprintf('\n Cache nBars mismatch (%d vs %d). Recomputing...\n', size(permStats_fig5.observedMean, 3), nBars_fig5);
        clear permStats_fig5;
    else
        cacheLoaded_fig5 = true;
    end
end
if ~cacheLoaded_fig5
    fprintf('\n Computing permutation nulls for Fig 5 (nPerm = %d)...\n', nPerm_fig5);
    % Create a struct to hold all permutation results
    permStats_fig5 = struct();
    permStats_fig5.nPerm = nPerm_fig5;
    permStats_fig5.fitModels = Bdiag_fig5;
    permStats_fig5.pairIdx = pairIdx_fig5;
    permStats_fig5.pairParams = {pairParams_fig5};
    permStats_fig5.barSigIdx = barSigIdx_fig5;
    permStats_fig5.barCzIdx = barCzIdx_fig5;
    permStats_fig5.barSigVal = gaborCST_unik(barSigIdx_fig5);
    permStats_fig5.barCzVal = cSDT_unik(barCzIdx_fig5);
    permStats_fig5.observedMean = nan(nCols_fig5, nRows_fig5, nBars_fig5);
    permStats_fig5.observedSem = nan(nCols_fig5, nRows_fig5, nBars_fig5);
    permStats_fig5.nRep = nan(nCols_fig5, nRows_fig5, nBars_fig5);
    permStats_fig5.null68 = nan(nCols_fig5, nRows_fig5, nBars_fig5, 2);
    permStats_fig5.null95 = nan(nCols_fig5, nRows_fig5, nBars_fig5, 2);
    permStats_fig5.permP = nan(nCols_fig5, nRows_fig5, nBars_fig5);
    permStats_fig5.sig95 = false(nCols_fig5, nRows_fig5, nBars_fig5);

    % Precompute matched-model subsets and bar membership once because they do not depend on parameter pair.
    Rdiag_fig5 = cell(1, nRows_fig5);
    rowGroups_fig5 = cell(nRows_fig5, nBars_fig5);
    for iModel = 1:nRows_fig5
        Rdiag_fig5{iModel} = R([R.iModelB_sim] == Bdiag_fig5(iModel) & [R.iModelB_fit] == Bdiag_fig5(iModel)); % only matched models are included in the correlation summary
        if isempty(Rdiag_fig5{iModel})
            continue;
        end
        gaborVals_fig5 = [Rdiag_fig5{iModel}.gaborCST];
        cSDTVals_fig5 = [Rdiag_fig5{iModel}.cSDT_true];
        for iB = 1:nBars_fig5
            rowGroups_fig5{iModel, iB} = find(gaborVals_fig5 == gaborCST_unik(barSigIdx_fig5(iB)) &  cSDTVals_fig5 == cSDT_unik(barCzIdx_fig5(iB)));
        end
    end % iModel

    for iPair = 1:nCols_fig5   % row = param pair
        iP1 = pairIdx_fig5(iPair, 1);
        iP2 = pairIdx_fig5(iPair, 2);
        p1 = pairParams_fig5{iP1};
        p2 = pairParams_fig5{iP2};

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

        for iModel = 1:nRows_fig5   % col = model variant
            fprintf('Processing pair %d/%d, model %d/%d...\n', iPair, nCols_fig5, iModel, nRows_fig5);
            fitModel = Bdiag_fig5(iModel);
            Rsub = Rdiag_fig5{iModel};

            % Skip non-identifiable pair/model combos: if either parameter is not
            % part of this reduced model, the pairwise correlation is undefined.
            modelParamNames_fig5 = namesModelBparams_short{fitModel};
            hasP1_fig5 = strcmp(p1, 'criterion_DV') || any(strcmp(modelParamNames_fig5, p1));
            hasP2_fig5 = strcmp(p2, 'criterion_DV') || any(strcmp(modelParamNames_fig5, p2));

            if ~(hasP1_fig5 && hasP2_fig5)
                continue;
            end

            yMean = nan(nBars_fig5, 1);
            ySem  = nan(nBars_fig5, 1);
            nRep  = zeros(nBars_fig5, 1);
            null68_fig5 = nan(nBars_fig5, 2);
            null95_fig5 = nan(nBars_fig5, 2);
            permP_fig5 = nan(nBars_fig5, 1);

            for iB = 1:nBars_fig5 % loop over signalCST x cSDT bars

                % Get row indices for this signalCST x cSDT condition in this model
                idxRows = rowGroups_fig5{iModel, iB};

                rVals = nan(numel(idxRows), 1);
                rNull = nan(nPerm_fig5, numel(idxRows));

                % Loop over rows (matched model instances) for this signalCST x cSDT condition, collect r values and null distributions.
                parfor iRR = 1:numel(idxRows)
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

                    if numel(xIter) < 2 || numel(yIter) < 2
                        continue;
                    end

                    xi = xIter;
                    yi = yIter;
                    if numel(unique(xi)) < 2 || numel(unique(yi)) < 2
                        continue;
                    end
                    %---------------------%
                    Sperm_fig5 = fxn_perm_corr(xi, yi, nPerm_fig5);
                    %---------------------%

                    rVals(iRR) = Sperm_fig5.rObs;
                    rNull(:, iRR) = Sperm_fig5.rPerm;
                end % iRR

                % Compute observed mean r, SEM, and null distribution percentiles for this bar.
                [yMean(iB), ~, ~, ySem(iB)] = getCI(rVals, 2, 1);
                nRep(iB) = numel(rVals);

                nullMean_fig5 = mean(rNull, 2, 'omitnan');
                [~, null68_fig5(iB, 1), null68_fig5(iB, 2)] = getCI(nullMean_fig5, 1, 1);
                [~, null95_fig5(iB, 1), null95_fig5(iB, 2)] = getCI(nullMean_fig5, 1, 1, 0.95);
                permP_fig5(iB) = (sum(abs(nullMean_fig5) >= abs(yMean(iB))) + 1) / (numel(nullMean_fig5) + 1);
            end % for iB

            % Skip non-existent panels (e.g., reduced models)
            if all(~isfinite(yMean))
                continue;
            end
            sig95_fig5 = isfinite(yMean) & isfinite(null95_fig5(:, 1)) & isfinite(null95_fig5(:, 2)) &  (yMean < null95_fig5(:, 1) | yMean > null95_fig5(:, 2));

            permStats_fig5.observedMean(iPair, iModel, :) = yMean;
            permStats_fig5.observedSem(iPair, iModel, :) = ySem;
            permStats_fig5.nRep(iPair, iModel, :) = nRep;
            permStats_fig5.null68(iPair, iModel, :, :) = reshape(null68_fig5, [1 1 nBars_fig5 2]);
            permStats_fig5.null95(iPair, iModel, :, :) = reshape(null95_fig5, [1 1 nBars_fig5 2]);
            permStats_fig5.permP(iPair, iModel, :) = permP_fig5;
            permStats_fig5.sig95(iPair, iModel, :) = sig95_fig5;

        end % iModel
    end % iPair

    save(permFile_fig5, 'permStats_fig5');
    cacheLoaded_fig5 = true;
end % if exist(permFile_fig5)

% PLOTTING
if cacheLoaded_fig5
    h = figure('Position', [80 80 330*nRows_fig5 210*nCols_fig5]);
    tiledlayout(nCols_fig5, nRows_fig5, 'TileSpacing', 'compact', 'Padding', 'compact');

    for iPair = 1:nCols_fig5
        for iModel = 1:nRows_fig5
            nexttile((iPair-1)*nRows_fig5 + iModel); hold on;

            yMean = squeeze(permStats_fig5.observedMean(iPair, iModel, :));
            ySem = squeeze(permStats_fig5.observedSem(iPair, iModel, :));
            nRep = squeeze(permStats_fig5.nRep(iPair, iModel, :));
            sig95_fig5 = squeeze(permStats_fig5.sig95(iPair, iModel, :));

            if all(~isfinite(yMean))
                % text(0.5, 0.5, 'N/A', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                axis off;
                continue;
            end

            nRepGood = nRep;
            dotSize_fig5 = nan(nBars_fig5, 1);
            if any(nRepGood ~= nRepGood(1))
                dotSize_fig5 = 80 + 240 * (nRepGood - min(nRepGood)) / max(max(nRepGood) - min(nRepGood), eps);
            else
                dotSize_fig5(:) = 180;
            end

            % Plot data dots
            scatter(xBar_fig5, yMean, dotSize_fig5, barColors_fig5, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.1);

            % Plot error bars
            for iPt = 1:numel(xBar_fig5)
                errorbar(xBar_fig5(iPt), yMean(iPt), ySem(iPt), 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1.2, 'CapSize', 0);
            end

            % Add aestercks
            for iStar = find(sig95_fig5(:))'
                if yMean(iStar) >= 0
                    yStar = min(0.95, yMean(iStar) + max(ySem(iStar), 0.03) + 0.05);
                    starAlign = 'bottom';
                else
                    yStar = max(-0.95, yMean(iStar) - max(ySem(iStar), 0.03) - 0.05);
                    starAlign = 'top';
                end
                text(xBar_fig5(iStar), yStar, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', starAlign, 'FontSize', 15, 'FontWeight', 'bold', 'Color', 'k');
            end

            yline(0, 'k-', 'LineWidth', 0.8);
            xlim([0.4, nBars_fig5 + 0.6]);
            ylim([-1, 1]);
            if iModel==1
                yticks([-1,0,1])
                ylabel('Pearson''s r')
            else
                yticks([])
            end
            for iSep = 1:max(nCz_fig5 - 1, 0)
                xline(iSep * nSig_fig5 + 0.5, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.9);
            end
            xticks([]);
            fxn_style_ax(gca, setting);
            set(gca, 'YGrid', 'off');
        end
    end

    sgtitle('Fig 5: r by signalCST x cSDT | rows: param pair, cols: model variant (Bfit=Bsim)', 'FontWeight', 'bold');
    set(findall(h,'-property','FontSize'), 'FontSize', 18);
    saveas(h, fullfile(nameFolder_Figures_part4, 'FigS5_paramCorr_summaryBars.png'));
    close(h);
end

% ── Color legend: rows = signalCST, cols = cSDT ──────────────────────
czMin_leg = min(cSDT_unik);    czMax_leg = max(cSDT_unik);
sigMin_leg = min(gaborCST_unik); sigMax_leg = max(gaborCST_unik);

sqSz  = 52;   % pixels per square
sqGap = 0.18; % adjustable gap between squares (0 to <1, in cell units)
marg  = 60;   % margin for axis labels
figW  = marg + nCz_fig5  * sqSz + marg;
figH  = marg + nSig_fig5 * sqSz + marg;

hLeg = figure('Position', [100 100 figW figH]);
axLeg = axes('Parent', hLeg);
hold(axLeg, 'on');

for iSig = 1:nSig_fig5
    sigVal = gaborCST_unik(iSig);
    if sigMax_leg > sigMin_leg
        sigNorm = (sigVal - sigMin_leg) / (sigMax_leg - sigMin_leg);
    else
        sigNorm = 0.5;
    end

    for iCz = 1:nCz_fig5
        czVal = cSDT_unik(iCz);
        if czMax_leg > czMin_leg
            czNorm = (czVal - czMin_leg) / (czMax_leg - czMin_leg);
        else
            czNorm = 0.5;
        end

        hCol  = (1 - czNorm) * (2/3);          % blue -> red
        sCol  = 0.25 + 0.70 * sigNorm;
        vCol  = 0.98 - 0.13 * sigNorm;
        faceC = hsv2rgb([hCol, sCol, vCol]);

        % iSig=1 → top row; flip y so row 1 is at top
        yBot = nSig_fig5 - iSig + sqGap / 2;
        yTop = nSig_fig5 - iSig + 1 - sqGap / 2;
        xL   = iCz - 1 + sqGap / 2;
        xR   = iCz - sqGap / 2;
        patch(axLeg, [xL xR xR xL], [yBot yBot yTop yTop], faceC,  'EdgeColor', 'none');
    end
end

% Show tick labels, but hide tick marks themselves.
set(axLeg, 'XTick', (0.5 : 1 : nCz_fig5 - 0.5),  'XTickLabel', arrayfun(@num2str, cSDT_unik(:)', 'UniformOutput', false));
set(axLeg, 'YTick', (0.5 : 1 : nSig_fig5 - 0.5),  'YTickLabel', arrayfun(@num2str, fliplr(gaborCST_unik(:)'), 'UniformOutput', false));
axLeg.XAxisLocation = 'top';
xlabel(axLeg, 'SDT criterion');
ylabel(axLeg, 'Signal contrast');
% title(axLeg, 'Color key: bar fill color');
xlim(axLeg, [0, nCz_fig5]);
ylim(axLeg, [0, nSig_fig5]);
axis(axLeg, 'equal');
set(axLeg, 'TickLength', [0 0]);
axLeg.XRuler.Axle.Visible = 'off';
axLeg.YRuler.Axle.Visible = 'off';
box(axLeg, 'off');
% fxn_style_ax(axLeg, setting);
set(findall(hLeg, '-property', 'FontSize'), 'FontSize', 20);

saveas(hLeg, fullfile(nameFolder_Figures_part4, 'FigS5_paramCorr_legend.png'));
close(hLeg);
fprintf('DONE\n\n')

%% Figure 6: Model recovery
fprintf('\n%s: Fig 6: Model recovery...', string(datetime('now')))

% Precompute index lookups once to avoid repeated struct filtering.
nSim_fig6 = numel(Bsim_unik);
nFit_fig6 = numel(Bfit_unik);
simAll_fig6 = [R.iModelB_sim]; simAll_fig6 = simAll_fig6(:);
fitAll_fig6 = [R.iModelB_fit]; fitAll_fig6 = fitAll_fig6(:);

if isfield(R, 'nameIO')
    datasetIDs_fig6 = string({R.nameIO})';
else
    datasetIDs_fig6 = string(arrayfun(@(r) sprintf('Nm%g_Na%g_Ns%g_G%g_cSDT%g',  r.Nmul_true, r.Nadd_true, r.Nshared_true, r.gaborCST, r.cSDT_true), R, 'UniformOutput', false))';
end

idsPerSim_fig6 = cell(nSim_fig6, 1);
rowsPerSimCondFit_fig6 = cell(nSim_fig6, 1);
rowsPerSimFit_fig6 = cell(nSim_fig6, nFit_fig6);

for iSim_fig6 = 1:nSim_fig6
    simVal_fig6 = Bsim_unik(iSim_fig6);
    idxSimRows_fig6 = find(simAll_fig6 == simVal_fig6);
    idsSim_fig6 = unique(datasetIDs_fig6(idxSimRows_fig6));
    idsPerSim_fig6{iSim_fig6} = idsSim_fig6;

    rowsCell_fig6 = cell(numel(idsSim_fig6), nFit_fig6);
    for iID_fig6 = 1:numel(idsSim_fig6)
        idxCondRows_fig6 = idxSimRows_fig6(datasetIDs_fig6(idxSimRows_fig6) == idsSim_fig6(iID_fig6));
        fitCond_fig6 = fitAll_fig6(idxCondRows_fig6);
        for iFit_fig6 = 1:nFit_fig6
            rowsCell_fig6{iID_fig6, iFit_fig6} = idxCondRows_fig6(fitCond_fig6 == Bfit_unik(iFit_fig6));
        end
    end
    rowsPerSimCondFit_fig6{iSim_fig6} = rowsCell_fig6;

    for iFit_fig6 = 1:nFit_fig6
        rowsPerSimFit_fig6{iSim_fig6, iFit_fig6} = idxSimRows_fig6(fitAll_fig6(idxSimRows_fig6) == Bfit_unik(iFit_fig6));
    end
end

% Win rate: for each condition, compute fraction of iterations won by each fit model,
% then average across conditions.
W_bestNLL = nan(nFit_fig6, nSim_fig6);
if ~isempty(R) && isfield(R, 'nLL_test_allIter')
    for iSim_wr = 1:nSim_fig6
        % fprintf('Computing win rates for Bsim=%d/%d...\n', iSim_wr, nSim_fig6);
        rowsCell_wr = rowsPerSimCondFit_fig6{iSim_wr};
        nCond_wr = size(rowsCell_wr, 1);
        winRatePerCond_wr = nan(nFit_fig6, nCond_wr);

        for iID_wr = 1:nCond_wr
            nIter_wr = 0;
            nllRows_wr = cell(1, nFit_fig6);
            for iFit_wr = 1:nFit_fig6
                rows_wr = rowsCell_wr{iID_wr, iFit_wr};
                if isempty(rows_wr), continue; end
                v_wr = R(rows_wr(1)).nLL_test_allIter(:)'; % preserve original first-match behavior
                nllRows_wr{iFit_wr} = v_wr;
                nIter_wr = max(nIter_wr, numel(v_wr));
            end
            if nIter_wr == 0, continue; end

            nllMat_wr = nan(nFit_fig6, nIter_wr);
            for iFit_wr = 1:nFit_fig6
                v_wr = nllRows_wr{iFit_wr};
                if isempty(v_wr), continue; end
                nllMat_wr(iFit_wr, 1:numel(v_wr)) = v_wr;
            end

            [minVals_wr, bestIdx_wr] = min(nllMat_wr, [], 1);
            validCol_wr = isfinite(minVals_wr);
            if any(validCol_wr)
                wins_wr = accumarray(bestIdx_wr(validCol_wr)', 1, [nFit_fig6, 1], @sum, 0);
                winRatePerCond_wr(:, iID_wr) = wins_wr / sum(validCol_wr);
            end
        end

        validCols_wr = any(isfinite(winRatePerCond_wr), 1);
        if any(validCols_wr)
            W_bestNLL(:, iSim_wr) = mean(winRatePerCond_wr(:, validCols_wr), 2, 'omitnan');
        end
    end
end

% Placeholders for RMSE for metrics and params
M_curveRMSE = nan(nFit_fig6, nSim_fig6);
M_paramRMSE = nan(nFit_fig6, nSim_fig6);

for iSim = 1:nSim_fig6
    for iFit = 1:nFit_fig6
        fitModel = Bfit_unik(iFit);
        rowsSimFit_fig6 = rowsPerSimFit_fig6{iSim, iFit};
        Rsub = R(rowsSimFit_fig6);
        if isempty(Rsub), continue; end

        % metric RMSE
        rmseMetricVals = [[Rsub.pYES_rmse]'; [Rsub.pC_rmse]'; [Rsub.pA_rmse]'];
        % rmseMetricVals = [Rsub.pA_rmse]';
        [M_curveRMSE(iFit, iSim), ~, ~] = getCI(rmseMetricVals, 2, 1);

        % parameter RMSE
        rmseParamFields = strcat(namesModelBparams_short{fitModel}, '_rmse');
        rmseParamVals = [];
        for iFld = 1:numel(rmseParamFields)
            vals = [Rsub.(rmseParamFields{iFld})]';
            rmseParamVals = [rmseParamVals; vals]; %#ok<AGROW>
        end
        [M_paramRMSE(iFit, iSim), ~, ~] = getCI(rmseParamVals, 2, 1);
    end
end

% ---------- delta NLL matrix (median ΔnLL from best fit, aggregated across all conditions) ----------
M_deltaNLL = nan(nFit_fig6, nSim_fig6);

if ~isempty(R) && isfield(R, 'nLL_med')
    for iSim_dn = 1:nSim_fig6
        rowsCell_dn = rowsPerSimCondFit_fig6{iSim_dn};
        nCond_dn = size(rowsCell_dn, 1);
        deltaVals_dn = cell(nFit_fig6, 1);

        for iICond_dn = 1:nCond_dn
            scores_dn = nan(nFit_fig6, 1);
            for iFit_dn = 1:nFit_fig6
                rows_dn = rowsCell_dn{iICond_dn, iFit_dn};
                vals_dn = [R(rows_dn).nLL_med]';
                if isempty(vals_dn)
                    scores_dn(iFit_dn)=nan;
                else
                [scores_dn(iFit_dn), ~, ~] = getCI(vals_dn, 2, 1);
                end
            end

            if ~any(isfinite(scores_dn)), continue; end
            bestScore_dn = min(scores_dn);

            for iFit_dn = 1:nFit_fig6
                if isfinite(scores_dn(iFit_dn))
                    deltaVals_dn{iFit_dn}(end+1, 1) = scores_dn(iFit_dn) - bestScore_dn; %#ok<AGROW>
                end
            end
        end

        for iFit_dn = 1:nFit_fig6
            [M_deltaNLL(iFit_dn, iSim_dn), ~, ~] = getCI(deltaVals_dn{iFit_dn}, 2, 1);
        end
    end
end

% Layout: 2x2
%   (1,1) Win rate      (1,2) Metric RMSE
%   (2,1) Median ΔnLL   (2,2) Param RMSE
panels_fig6 = {  W_bestNLL,   'Win rate (% iters, avg over conds)',  [0 1];  M_deltaNLL,  'ΔNLL from the best model', [ 0 160];  M_curveRMSE, 'Metrics RMSE',            [0.05, 0.07];  M_paramRMSE, 'Parameters RMSE',             [1 2.8] };

% tile order in a 2x2 tiledlayout (row-major): (1,1),(1,2),(2,1),(2,2)
% desired: WinRate, MetricRMSE, DeltaNLL, ParamRMSE
tileOrder = [1, 3, 2, 4]; % index into panels_fig6

h = figure('Position', [100 100 900 700]);
tlo = tiledlayout(2, 2, 'TileSpacing', 'loose', 'Padding', 'loose');

% Use model names for axis ticks when available.
xTickLabels_fig6 = arrayfun(@(k) sprintf('M%d', k), Bsim_unik, 'UniformOutput', false);
yTickLabels_fig6 = arrayfun(@(k) sprintf('M%d', k), Bfit_unik, 'UniformOutput', false);
if exist('namesModelB', 'var') && ~isempty(namesModelB)
    for iLab = 1:numel(Bsim_unik)
        k = Bsim_unik(iLab);
        if k >= 1 && k <= numel(namesModelB)
            xTickLabels_fig6{iLab} = char(string(namesModelB{k}));
        end
    end
    for iLab = 1:numel(Bfit_unik)
        k = Bfit_unik(iLab);
        if k >= 1 && k <= numel(namesModelB)
            yTickLabels_fig6{iLab} = char(string(namesModelB{k}));
        end
    end
end

% Loop over tiles in specified order to control color limits and annotations per panel.
nTiltes=4;
% fprintf('Plotting %d tiles in order: %s\n', nTiltes);
for iTile = 1:nTiltes
    % fprintf('Plotting tile %d/%d: %s...\n', iTile, nTiltes, panels_fig6{tileOrder(iTile), 2});
    iPan = tileOrder(iTile);
    M_plot   = panels_fig6{iPan, 1};
    ttl_plot = panels_fig6{iPan, 2};
    clim_fix = panels_fig6{iPan, 3};
    isWinRatePanel = (iPan == 1);

    if isWinRatePanel
        M_plot = 100 * M_plot;
        ttl_plot = 'Win rate (%)';
    end

    % Use robust panel-wise clipping so outliers do not dominate colormap.
    finVals = M_plot(isfinite(M_plot));
    if ~isempty(finVals)
        q = prctile(finVals, [2.5 97.5]);
        cLo = q(1);
        cHi = q(2);

        if ~isfinite(cLo) || ~isfinite(cHi) || cLo >= cHi
            cLo = min(finVals);
            cHi = max(finVals);
        end
        if cLo == cHi
            cLo = cLo - eps;
            cHi = cHi + eps;
        end
    else
        cLo = 0;
        cHi = 1;
    end
    M_plot_disp = min(max(M_plot, cLo), cHi);

    nexttile; hold on;
    imagesc(1:numel(Bsim_unik), 1:numel(Bfit_unik), M_plot_disp);
    set(gca, 'YDir', 'normal',  'XTick', 1:numel(Bsim_unik), 'XTickLabel', xTickLabels_fig6,  'YTick', 1:numel(Bfit_unik), 'YTickLabel', yTickLabels_fig6,  'XLim', [0.5, numel(Bsim_unik)+0.5],  'YLim', [0.5, numel(Bfit_unik)+0.5]);
    xlabel('Generating model'); ylabel('Fitting model');
    title(ttl_plot);
    % fxn_style_ax(gca, setting);
    axis square;

    % Color axis uses clipped 95% interval range for this panel.
    clim([cLo cHi]);
    % show colorbar
    cbar_fig6 = colorbar('Location', 'eastoutside');
    fxn_set_colorbar_ticks(cbar_fig6, cLo, cHi);

    % Annotate values with adaptive font color (low → white, high → black)
    % fprintf('    Annotating values for tile %d/%d...\n', iTile, nTiltes);
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
            if isWinRatePanel
                txtVal = sprintf('%.1f%%', v);
            else
                txtVal = sprintf('%.3f', v);
            end
            text(iSim_ann, iFit_ann, txtVal,  'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',  'FontSize', 8, 'Color', txtClr);
        end
    end


    % Highlight diagonal cells (matched sim/fit model) with thick red borders
    nDiag = min(numel(Bfit_unik), numel(Bsim_unik));
    for iDiag = 1:nDiag
        rectangle('Position', [iDiag - 0.5, iDiag - 0.5, 1, 1], 'EdgeColor', 'r', 'LineWidth', 2.5, 'LineStyle', '-');
    end

    % For each column, highlight the best cell with a black box:
    %   Win rate panel (iPan==1): highest value is best.
    %   All other panels: lowest value is best.
    for iSim_box = 1:numel(Bsim_unik)
        colVals = M_plot(:, iSim_box);
        finMask = isfinite(colVals);

        if isWinRatePanel
            [~, bestRow] = max(colVals(finMask));
        else
            [~, bestRow] = min(colVals(finMask));
        end
        finIdx = find(finMask);
        bestRow = finIdx(bestRow);
        % Only draw black box if this cell is NOT on the diagonal
        % if bestRow ~= iSim_box
        rectangle('Position', [iSim_box - 0.5, bestRow - 0.5, 1, 1], 'EdgeColor', 'k', 'LineWidth', 2.5, 'LineStyle', '--');
        % end
    end

end % for iTile

sgtitle('Fig 6: Model recovery', 'FontWeight', 'bold');
% fprintf('Saving figure...\n');
nameFile_fig6 = fullfile(nameFolder_Figures_part4, 'FigS6_modelRecovery.png');
saveas(h, nameFile_fig6);
% fprintf('Closing figure...\n');
close(h);

fprintf('DONE\n\n')





% %% Figure 6B: parameter-pair correlations
% fprintf('\n%s: Fig 6B: Parameter-pair correlations\n', string(datetime('now')))
%
% for iFit = 1:numel(Bfit_unik)
%     fprintf('Processing fit model %d/%d...\n', iFit, numel(Bfit_unik));
%     fitModel = Bfit_unik(iFit);
%
%     % Use matched pairs to keep one coherent parameterization per model.
%     Rsub = R([R.iModelB_sim] == fitModel & [R.iModelB_fit] == fitModel);
%     if isempty(Rsub), continue; end
%
%     pNames_all = namesModelBparams_short{fitModel};
%     % Use per-iteration estimate arrays as source for dots
%     iterFields = cell(size(pNames_all));
%     for iP = 1:numel(pNames_all)
%         if strcmp(pNames_all{iP}, 'criterion_DV')
%             iterFields{iP} = 'criterion_DV_est_allIter';
%         else
%             iterFields{iP} = sprintf('%s_est_allIter', pNames_all{iP});
%         end
%     end
%
%     keepParam = false(size(pNames_all));
%     for iP = 1:numel(pNames_all)
%         if ~isfield(Rsub, iterFields{iP}), continue; end
%         vals = vertcat(Rsub.(iterFields{iP}));
%         keepParam(iP) = any(isfinite(vals(:)));
%     end
%
%     pNames = pNames_all(keepParam);
%     iterFields = iterFields(keepParam);
%     nParam = numel(pNames);
%     if nParam < 2, continue; end
%
%     sigVals_all = [Rsub.gaborCST]';
%     czVals_all  = [Rsub.cSDT_true]';
%     [isSig, sigIdx_all] = ismember(sigVals_all, gaborCST_unik);
%     [isCz,  czIdx_all]  = ismember(czVals_all,  cSDT_unik);
%     validGroupBase = isfinite(sigVals_all) & isfinite(czVals_all) & isSig & isCz;
%
%     groupRows = cell(numel(gaborCST_unik), numel(cSDT_unik));
%     for iSig = 1:numel(gaborCST_unik)
%         for iCz = 1:numel(cSDT_unik)
%             groupRows{iSig, iCz} = find(validGroupBase & sigIdx_all == iSig & czIdx_all == iCz);
%         end
%     end
%
%     % Grayscale color per condition based on sumsqr(Nmul_true, Nadd_true, Nshared_true).
%     % Black = highest sumsqr, light gray = lowest.
%     inFields_gray = {'Nmul_true', 'Nadd_true', 'Nshared_true'};
%     sumsqr_cond = zeros(numel(Rsub), 1);
%     for iGF = 1:numel(inFields_gray)
%         if isfield(Rsub, inFields_gray{iGF})
%             v = [Rsub.(inFields_gray{iGF})]';
%             v(~isfinite(v)) = 0;
%             sumsqr_cond = sumsqr_cond + v .^ 2;
%         end
%     end
%     sq_min = min(sumsqr_cond);
%     sq_max = max(sumsqr_cond);
%     if sq_max > sq_min
%         sq_norm = (sumsqr_cond - sq_min) / (sq_max - sq_min);
%     else
%         sq_norm = zeros(size(sumsqr_cond));
%     end
%     lightGrayVal = 0.78;
%     % normVal=1 → black (0), normVal=0 → light gray (lightGrayVal)
%     condGrayColor = (1 - sq_norm) * lightGrayVal .* ones(numel(Rsub), 3);
%
%     pairIdx = nchoosek(1:nParam, 2);
%     nPairs  = size(pairIdx, 1);
%     nCols   = min(3, nPairs);
%     nRows   = ceil(nPairs / nCols);
%
%     h = figure('Position', [100 100 420*nCols 320*nRows]);
%     tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
%
%     for iPair = 1:nPairs
%         iP1 = pairIdx(iPair, 1);
%         iP2 = pairIdx(iPair, 2);
%
%         nexttile; hold on;
%
%         % Collect all iter values pooled across conditions for overall r
%         xPool = [];
%         yPool = [];
%
%         for iSig = 1:numel(gaborCST_unik)
%             for iCz = 1:numel(cSDT_unik)
%                 idxRows = groupRows{iSig, iCz};
%                 if isempty(idxRows), continue; end
%
%                 mkr = setting.marker_signal{iSig};
%
%                 for iRow = 1:numel(idxRows)
%                     rIdx = idxRows(iRow);
%                     col  = condGrayColor(rIdx, :);
%                     xIter = Rsub(rIdx).(iterFields{iP1})(:);
%                     yIter = Rsub(rIdx).(iterFields{iP2})(:);
%                     xi = xIter;
%                     yi = yIter;
%                     xPool = [xPool; xi]; %#ok<AGROW>
%                     yPool = [yPool; yi]; %#ok<AGROW>
%
%                     % dots: one per iteration
%                     plot(xi, yi, mkr,  'LineStyle', 'none',  'MarkerSize', max(setting.markerSize - 2, 4), 'MarkerFaceColor', col,  'MarkerEdgeColor', 'w',  'LineWidth', 0.5);
%
%                     % regression line per condition
%                     if numel(xi) >= 2 && numel(unique(xi)) >= 2
%                         pfit_cond = polyfit(xi, yi, 1);
%                         xfit_cond = linspace(min(xi), max(xi), 60);
%                         plot(xfit_cond, polyval(pfit_cond, xfit_cond), '-',  'Color', col, 'LineWidth', 1.4);
%                     end
%                 end
%             end
%         end
%
%         if numel(xPool) >= 2
%             [r_val, p_val] = corr(xPool, yPool);
%             title(sprintf('%s vs %s | r=%.2f, p=%.3g', pNames{iP1}, pNames{iP2}, r_val, p_val));
%         else
%             title(sprintf('%s vs %s', pNames{iP1}, pNames{iP2}));
%         end
%
%         axis square;
%         xlabel(sprintf('Est. %s', pNames{iP1}));
%         ylabel(sprintf('Est. %s', pNames{iP2}));
%         fxn_style_ax(gca, setting);
%     end
%
%     % Figure-level legends: grayscale = sumsqr(IN params), marker = signalCST
%     axBase = gca;
%     pos = axBase.Position;
%
%     % Build a few representative grayscale legend entries
%     axColor = axes('Position', pos, 'Color', 'none', 'Visible', 'off',  'HitTest', 'off', 'HandleVisibility', 'off');
%     hold(axColor, 'on');
%     sq_unik_sorted = unique(sumsqr_cond);
%     nGrayLeg = min(numel(sq_unik_sorted), 5);
%     legIdxGray = round(linspace(1, numel(sq_unik_sorted), nGrayLeg));
%     hGray = gobjects(nGrayLeg, 1);
%     grayLegLabels = cell(nGrayLeg, 1);
%     for iGL = 1:nGrayLeg
%         sqVal = sq_unik_sorted(legIdxGray(iGL));
%         normGL = (sq_max > sq_min) * (sqVal - sq_min) / max(sq_max - sq_min, eps);
%         grayGL = (1 - normGL) * lightGrayVal;
%         hGray(iGL) = plot(axColor, nan, nan, 'o',  'LineStyle', 'none',  'MarkerFaceColor', grayGL * [1 1 1],  'MarkerEdgeColor', grayGL * [1 1 1],  'MarkerSize', 7);
%         grayLegLabels{iGL} = sprintf('%.3g', sqVal);
%     end
%     leg1 = legend(axColor, hGray, grayLegLabels,  'Location', 'northeast', 'Box', 'on', 'FontSize', max(setting.fontSize - 1, 8));
%     title(leg1, 'Total internal var.');
%
%     axSig = axes('Position', pos, 'Color', 'none', 'Visible', 'off',  'HitTest', 'off', 'HandleVisibility', 'off');
%     hold(axSig, 'on');
%     hSig = gobjects(numel(gaborCST_unik), 1);
%     for iSig = 1:numel(gaborCST_unik)
%         hSig(iSig) = plot(axSig, nan, nan,  'LineStyle', 'none',  'Marker', setting.marker_signal{iSig},  'MarkerFaceColor', 'w',  'MarkerEdgeColor', 'k',  'MarkerSize', 7);
%     end
%     leg2 = legend(axSig, hSig,  arrayfun(@(x) sprintf('%g', x), gaborCST_unik, 'UniformOutput', false),  'Location', 'southeast', 'Box', 'on', 'FontSize', max(setting.fontSize - 1, 8));
%     title(leg2, 'signalCST');
%
%     sgtitle(sprintf('Fig 6B: parameter-pair correlations | Bfit=Bsim=%d', fitModel), 'FontWeight', 'bold');
%     saveas(h, fullfile(nameFolder_Figures_part4, sprintf('FigS6B_paramCorr_B%d.png', fitModel)));
%     close(h);
% end % iModelB_fit
%
% %% Figure 7: Full vs. reduced models—Iteration-level ΔnLL distributions (Generating model is the full model B1)
% %    For each reduced model Br, shows the per-iteration distribution of
% %    ΔnLL(Br − B1) aggregated across conditions.
% %    A positive ΔnLL means B1 fits the held-out data better on that iteration.
% %    If distributions straddle 0, B1 barely wins and identifiability is weak.
%
% fprintf('\n%s: Fig 7: Iteration-level delta-nLL distributions\n', string(datetime('now')))
%
% Br_fig7 = setdiff(Bfit_unik, 1); % reduced models to compare against B1
% nReduced_f7 = numel(Br_fig7);
%
% if ~isfield(R, 'nLL_test_allIter') || ~any([R.iModelB_sim] == 1) || nReduced_f7 == 0
%     fprintf('Skipping Fig 7: nLL_test_allIter not available or no B1-simulated data.\n');
% else
%     if isfield(R, 'nameIO')
%         datasetIDs_f7 = string({R.nameIO});
%     else
%         datasetIDs_f7 = string(arrayfun(@(r) sprintf('Nm%g_Na%g_Ns%g_G%g_cSDT%g',  r.Nmul_true, r.Nadd_true, r.Nshared_true, r.gaborCST, r.cSDT_true), R, 'UniformOutput', false));
%     end
%     idxSim1_f7 = [R.iModelB_sim] == 1;
%     idsSim1_f7 = unique(datasetIDs_f7(idxSim1_f7));
%     nConds_f7  = numel(idsSim1_f7);
%
%     dnLL_iter_allConds = cell(nReduced_f7, nConds_f7);
%     for iCond = 1:nConds_f7
%         idxCond = idxSim1_f7 & datasetIDs_f7 == idsSim1_f7(iCond);
%         Rcond   = R(idxCond);
%         if isempty(Rcond), continue; end
%
%         idxB1 = find([Rcond.iModelB_fit] == 1, 1);
%         if isempty(idxB1) || isempty(Rcond(idxB1).nLL_test_allIter), continue; end
%         nll_B1 = Rcond(idxB1).nLL_test_allIter(:);
%
%         for iR = 1:nReduced_f7
%             Br      = Br_fig7(iR);
%             idxBr   = find([Rcond.iModelB_fit] == Br, 1);
%             if isempty(idxBr) || isempty(Rcond(idxBr).nLL_test_allIter), continue; end
%             nll_Br  = Rcond(idxBr).nLL_test_allIter(:);
%             nIt     = min(numel(nll_B1), numel(nll_Br));
%             dnLL_iter_allConds{iR, iCond} = nll_Br(1:nIt) - nll_B1(1:nIt);
%         end
%     end
%
%     h = figure('Position', [100 100 280*nReduced_f7 480]);
%     tiledlayout(1, nReduced_f7, 'TileSpacing', 'loose', 'Padding', 'loose');
%
%     for iR = 1:nReduced_f7
%         Br = Br_fig7(iR);
%         nexttile; hold on;
%
%         condMeds = nan(nConds_f7, 1);
%         for iCond = 1:nConds_f7
%             v = dnLL_iter_allConds{iR, iCond};
%             if ~isempty(v), condMeds(iCond) = median(v, 'omitnan'); end
%         end
%
%         idx = isfinite(condMeds);
%         if any(idx)
%             xJit = 1 + 0.12 * (rand(sum(idx), 1) - 0.5);
%             scatter(xJit, condMeds(idx), 28, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5);
%             [gMed, ~, ~, gSEM] = getCI(condMeds(idx), 1, 1);
%             errorbar(1, gMed, gSEM, 'ko',  'MarkerFaceColor', 'k', 'MarkerSize', 7, 'LineWidth', 1.8, 'CapSize', 0);
%             winRate_f7 = mean(condMeds(idx) > 0);
%             text(1.35, gMed, sprintf('win=%.0f%%', winRate_f7*100),  'FontSize', setting.fontSize, 'VerticalAlignment', 'middle');
%         end
%
%         yline(0, 'k--', 'LineWidth', 1.2);
%         xlim([0.5 1.8]); xticks([]);
%         ylabel('\Delta nLL (B_r - B_1) per condition median');
%         title(sprintf('B%d - B1', Br));
%         fxn_style_ax(gca, setting);
%     end
%
%     sgtitle('Fig 7: Full vs. reduced models\Delta nLL distributions | Bsim=Full model', 'FontWeight', 'bold');
%     saveas(h, fullfile(nameFolder_Figures_part4, 'FigS7_dnLL_distributions.png'));
%     close(h);
% end
%
% %% Figure 8: ΔnLL vs true omitted parameter value (Generating model is the full model B1)
% %    For each reduced model Br, plots per-condition median ΔnLL(Br − B1) on the
% %    y-axis against the true value of the parameter Br omits on the x-axis.
% %    Expected: ΔnLL should increase as the omitted parameter grows away from 0.
% %    Weak or flat relationship = identifiability problem even with large true parameter.
%
% fprintf('\n%s: Fig 8: delta-nLL vs true omitted parameter\n', string(datetime('now')))
%
% omitB_fig8      = [2,              3,            4,          5,          5,              6,          6,              7,          7        ];
% omitParam_fig8  = {'Nshared',     'Nmul',       'Nadd',     'Nadd',     'Nshared',      'Nmul',     'Nshared',      'Nmul',     'Nadd'   };
% trueField_fig8  = {'Nshared_true','Nmul_true',  'Nadd_true','Nadd_true','Nshared_true', 'Nmul_true','Nshared_true', 'Nmul_true','Nadd_true'};
%
% keepPair = ismember(omitB_fig8, Bfit_unik);
% omitB_fig8     = omitB_fig8(keepPair);
% omitParam_fig8 = omitParam_fig8(keepPair);
% trueField_fig8 = trueField_fig8(keepPair);
% nPairs_f8 = numel(omitB_fig8);
%
% if ~isfield(R, 'nLL_test_allIter') || ~any([R.iModelB_sim] == 1) || nPairs_f8 == 0
%     fprintf('Skipping Fig 8: nLL_test_allIter not available or no valid model pairs.\n');
% else
%     if isfield(R, 'nameIO')
%         datasetIDs_f8 = string({R.nameIO});
%     else
%         datasetIDs_f8 = string(arrayfun(@(r) sprintf('Nm%g_Na%g_Ns%g_G%g_cSDT%g',  r.Nmul_true, r.Nadd_true, r.Nshared_true, r.gaborCST, r.cSDT_true), R, 'UniformOutput', false));
%     end
%     idxSim1_f8 = [R.iModelB_sim] == 1;
%     idsSim1_f8 = unique(datasetIDs_f8(idxSim1_f8));
%
%     paletteSrc = lines(max(nPairs_f8, 3));
%
%     h = figure('Position', [100 100 380*nPairs_f8 420]);
%     tiledlayout(1, nPairs_f8, 'TileSpacing', 'loose', 'Padding', 'loose');
%
%     for iP = 1:nPairs_f8
%         Br          = omitB_fig8(iP);
%         tField      = trueField_fig8{iP};
%         pLabel      = omitParam_fig8{iP};
%
%         xvals = nan(numel(idsSim1_f8), 1);
%         yvals = nan(numel(idsSim1_f8), 1);
%
%         for iCond = 1:numel(idsSim1_f8)
%             idxCond = idxSim1_f8 & datasetIDs_f8 == idsSim1_f8(iCond);
%             Rcond   = R(idxCond);
%             if isempty(Rcond), continue; end
%
%             if isfield(Rcond, tField), xvals(iCond) = Rcond(1).(tField); end
%
%             idxB1 = find([Rcond.iModelB_fit] == 1, 1);
%             idxBr = find([Rcond.iModelB_fit] == Br, 1);
%             if isempty(idxB1) || isempty(idxBr), continue; end
%             nll_B1 = Rcond(idxB1).nLL_test_allIter(:);
%             nll_Br = Rcond(idxBr).nLL_test_allIter(:);
%             if isempty(nll_B1) || isempty(nll_Br), continue; end
%             nIt = min(numel(nll_B1), numel(nll_Br));
%             yvals(iCond) = median(nll_Br(1:nIt) - nll_B1(1:nIt), 'omitnan');
%         end
%
%         idx = isfinite(xvals) & isfinite(yvals);
%         nexttile; hold on;
%
%         scatter(xvals(idx), yvals(idx), 40, paletteSrc(iP, :),  'filled', 'MarkerFaceAlpha', 0.75);
%
%         if sum(idx) >= 3 && range(xvals(idx)) > eps(max(abs(xvals(idx)))) * numel(xvals(idx))
%             pCoef = polyfit(xvals(idx) - mean(xvals(idx)), yvals(idx), 1);
%             xFit  = linspace(min(xvals(idx)), max(xvals(idx)), 60);
%             plot(xFit, polyval(pCoef, xFit - mean(xvals(idx))), '-', 'Color', paletteSrc(iP,:), 'LineWidth', 1.8);
%             [rho_f8, ~] = corr(xvals(idx), yvals(idx));
%             text(0.05, 0.92, sprintf('r = %.2f', rho_f8), 'Units', 'normalized',  'FontSize', setting.fontSize, 'Color', paletteSrc(iP,:));
%         end
%
%         yline(0, 'k--', 'LineWidth', 1.2);
%         xlabel(sprintf('True %s', pLabel));
%         ylabel(sprintf('\\Delta nLL median (B%d - B1)', Br));
%         title(sprintf('B%d vs B1 | omitted: %s', Br, pLabel));
%         fxn_style_ax(gca, setting);
%     end
%
%     sgtitle('Fig 8: \Delta nLL vs true omitted parameter | Bsim=B1', 'FontWeight', 'bold');
%     saveas(h, fullfile(nameFolder_Figures_part4, 'FigS8_dnLL_vs_trueParam.png'));
%     close(h);
% end
%
% %% Figure 9: Correlate template RMSE with each parameter's RMSE across simulating conditions
%
% fprintf('\n%s: Fig 9: Template RMSE vs parameter RMSE correlations\n', string(datetime('now')))
%
% % Use matched pairs (Bsim == Bfit) only, so each row is a consistent condition.
% % For each panel: scatter template_rmse (x) vs that parameter's RMSE (y) and annotate r.
% % Panels 1-4: RMSE of fitted noise params + DV criterion (never excluded).
% % Panel 5:    True SDT criterion |cSDT_true| — shows if criterion bias predicts template error.
%
% % --- RMSE-based panels (one per fitted parameter) ---
% rmseParamNames_fig9  = {'Nmul',       'Nadd',       'Nshared',       'criterion_DV'};
% rmseParamLabels_fig9 = {'N_{mul}',    'N_{add}',    'N_{shared}',    'Criterion DV'};
%
% % Keep only those with at least one finite value; criterion_DV is always kept.
% Rmatch_fig9 = R([R.iModelB_sim] == [R.iModelB_fit]);
% keepRmse_fig9 = false(1, numel(rmseParamNames_fig9));
% for iP9 = 1:numel(rmseParamNames_fig9)
%     fld = sprintf('%s_rmse', rmseParamNames_fig9{iP9});
%     if strcmp(rmseParamNames_fig9{iP9}, 'criterion_DV')
%         keepRmse_fig9(iP9) = true;   % never exclude criterion DV
%     elseif isfield(Rmatch_fig9, fld) && any(isfinite([Rmatch_fig9.(fld)]))
%         keepRmse_fig9(iP9) = true;
%     end
% end
% rmseParamNames_fig9  = rmseParamNames_fig9(keepRmse_fig9);
% rmseParamLabels_fig9 = rmseParamLabels_fig9(keepRmse_fig9);
%
% nPanels_fig9 = numel(rmseParamNames_fig9) + 1;   % +1 for SDT criterion panel
% nCols_fig9   = min(nPanels_fig9, 3);
% nRows_fig9   = ceil(nPanels_fig9 / nCols_fig9);
% cmap_fig9    = lines(numel(Bfit_unik));
%
% h = figure('Position', [100 100 380*nCols_fig9 320*nRows_fig9]);
% tiledlayout(nRows_fig9, nCols_fig9, 'TileSpacing', 'compact', 'Padding', 'compact');
%
% % --- Panels: template RMSE vs each parameter RMSE ---
% for iP9 = 1:numel(rmseParamNames_fig9)
%     pName9    = rmseParamNames_fig9{iP9};
%     pLabel9   = rmseParamLabels_fig9{iP9};
%     rmseField = sprintf('%s_rmse', pName9);
%
%     nexttile; hold on;
%
%     hLines_fig9 = gobjects(numel(Bfit_unik), 1);
%     for iBfit9 = 1:numel(Bfit_unik)
%         fitModel9 = Bfit_unik(iBfit9);
%         Rsub9 = R([R.iModelB_sim] == fitModel9 & [R.iModelB_fit] == fitModel9);
%         if isempty(Rsub9) || ~isfield(Rsub9, rmseField), continue; end
%
%         x9 = [Rsub9.template_rmse]';
%         y9 = [Rsub9.(rmseField)]';
%         idx9 = isfinite(x9) & isfinite(y9);
%         if sum(idx9) < 2, continue; end
%
%         col9 = cmap_fig9(iBfit9, :);
%         hLines_fig9(iBfit9) = scatter(x9(idx9), y9(idx9), 28, col9,  'filled', 'MarkerFaceAlpha', 0.65);
%
%         if numel(unique(x9(idx9))) >= 2
%             pCoef9 = polyfit(x9(idx9), y9(idx9), 1);
%             xFit9  = linspace(min(x9(idx9)), max(x9(idx9)), 80);
%             plot(xFit9, polyval(pCoef9, xFit9), '-', 'Color', col9,  'LineWidth', 1.6, 'HandleVisibility', 'off');
%         end
%
%         [r9, ~] = corr(x9(idx9), y9(idx9));
%         hLines_fig9(iBfit9).DisplayName = sprintf('B%d: r = %.2f', fitModel9, r9);
%     end % iBfit9
%
%     xlabel('Template RMSE');
%     ylabel(sprintf('%s RMSE', pLabel9));
%     title(sprintf('Template RMSE vs %s RMSE', pLabel9));
%     fxn_style_ax(gca, setting);
%     set(gca, 'YGrid', 'off');
%     axis square;
%     validH9 = hLines_fig9(isgraphics(hLines_fig9));
%     if ~isempty(validH9)
%         legend(validH9, 'Location', 'best', 'Box', 'off', 'FontSize', max(setting.fontSize-1,7));
%     end
% end % iP9
%
% % --- Final panel: template RMSE vs true SDT criterion |cSDT_true| ---
% nexttile; hold on;
% hLines_sdt = gobjects(numel(Bfit_unik), 1);
% for iBfit9 = 1:numel(Bfit_unik)
%     fitModel9 = Bfit_unik(iBfit9);
%     Rsub9 = R([R.iModelB_sim] == fitModel9 & [R.iModelB_fit] == fitModel9);
%     if isempty(Rsub9) || ~isfield(Rsub9, 'cSDT_true'), continue; end
%
%     x9   = [Rsub9.template_rmse]';
%     y9   = abs([Rsub9.cSDT_true]');
%     idx9 = isfinite(x9) & isfinite(y9);
%     if sum(idx9) < 2, continue; end
%
%     col9 = cmap_fig9(iBfit9, :);
%     hLines_sdt(iBfit9) = scatter(x9(idx9), y9(idx9), 28, col9,  'filled', 'MarkerFaceAlpha', 0.65);
%
%     if numel(unique(x9(idx9))) >= 2
%         pCoef9 = polyfit(x9(idx9), y9(idx9), 1);
%         xFit9  = linspace(min(x9(idx9)), max(x9(idx9)), 80);
%         plot(xFit9, polyval(pCoef9, xFit9), '-', 'Color', col9,  'LineWidth', 1.6, 'HandleVisibility', 'off');
%     end
%
%     [r9, ~] = corr(x9(idx9), y9(idx9));
%     hLines_sdt(iBfit9).DisplayName = sprintf('B%d r=%.2f', fitModel9, r9);
% end
% xlabel('Template RMSE');
% ylabel('|SDT criterion| (|Cz_{true}|)');
% title('Template RMSE vs SDT criterion');
% fxn_style_ax(gca, setting);
% set(gca, 'YGrid', 'off');
% axis square;
% validH_sdt = hLines_sdt(isgraphics(hLines_sdt));
% if ~isempty(validH_sdt)
%     legend(validH_sdt, 'Location', 'best', 'Box', 'off', 'FontSize', max(setting.fontSize-1,7));
% end
%
% sgtitle('Fig 9: Template RMSE vs Parameter RMSE | Matched Bsim=Bfit', 'FontWeight', 'bold');
% saveas(h, fullfile(nameFolder_Figures_part4, 'FigS9_Corr_tempRMSE_paramRMSE.png'));
% close(h);

%% Local helpers  (only functions used 2+ times)

function info = fxn_parse_nameIO(nameIO)

info = struct();

tok = regexp(nameIO,  ['IO_cN([\d\.]+)_cG([\d\.]+)_nT([A-Za-z0-9\.]+)' ...
    '_Nm([\d\.]+)_Na([\d\.]+)_Ns([\d\.]+)_cSDT(-?[\d\.]+)' ...
    '_cont[\d\.]+_whiten[\d\.]+_R(\d+)_\d+\d+_B(\d+)'],  'tokens', 'once');

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
info.flag_regressType = str2double(tok{8});
info.iModelB_sim = str2double(tok{9});
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

    idx = isfinite(yData) & isfinite(yPred) & isfinite(nTrials) & (nTrials > 0);
    yData = yData(idx);
    yPred = yPred(idx);
    w = nTrials(idx);

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
set(ax, 'FontSize', P.fontSize, 'LineWidth', P.axisLineWidth,  'YGrid', 'on', 'GridAlpha', P.gridAlpha, 'TickDir', P.tickDir);
end

function fxn_set_colorbar_ticks(cbar, cLo, cHi)
ticks = fxn_nice_ticks(cLo, cHi, 5);
ticks(abs(ticks) < 1e-12) = 0;
tickLabels = arrayfun(@(x) fxn_format_tick_label(x), ticks, 'UniformOutput', false);
set(cbar, 'Ticks', ticks, 'TickLabels', tickLabels);
end

function ticks = fxn_nice_ticks(cLo, cHi, nTicks)
if nargin < 3 || isempty(nTicks)
    nTicks = 5;
end

if ~isfinite(cLo) || ~isfinite(cHi)
    ticks = linspace(0, 1, nTicks);
    return
end

if cHi < cLo
    tmp = cLo;
    cLo = cHi;
    cHi = tmp;
end

span = cHi - cLo;
if span <= eps(max(abs([cLo cHi 1])))
    step = fxn_nice_step(max(abs(cHi), 1) / (nTicks - 1));
    start = cLo - floor((nTicks - 1) / 2) * step;
    ticks = start + (0:nTicks-1) * step;
    return
end

desiredStep = span / (nTicks - 1);

% Try a few nearby nice steps and pick the best one that covers the range.
stepCandidates = unique([ ...
    fxn_nice_step(desiredStep * 0.5), ...
    fxn_nice_step(desiredStep), ...
    fxn_nice_step(desiredStep * 1.5), ...
    fxn_nice_step(desiredStep * 2.0), ...
    fxn_nice_step(desiredStep * 3.0)]);
stepCandidates = stepCandidates(isfinite(stepCandidates) & stepCandidates > 0);

bestScore = inf;
bestTicks = [];
midVal = 0.5 * (cLo + cHi);

for iStep = 1:numel(stepCandidates)
    step = stepCandidates(iStep);
    totalSpan = step * (nTicks - 1);

    kMin = ceil((cHi - totalSpan) / step);
    kMax = floor(cLo / step);
    if kMin > kMax
        continue;
    end

    startVals = (kMin:kMax) * step;
    for iStart = 1:numel(startVals)
        start = startVals(iStart);
        cand = start + (0:nTicks-1) * step;

        % Coverage (hard requirement)
        if cand(1) > cLo + 1e-12 || cand(end) < cHi - 1e-12
            continue;
        end

        % Score = small outside padding + centered around data mid + optional zero anchor bonus.
        pad = (cLo - cand(1)) + (cand(end) - cHi);
        midErr = abs(mean(cand) - midVal);
        zeroBonus = 0;

        if cLo >= 0 && cHi < 1
            % Prefer 0, 0.2, ... style for bounded rates.
            if abs(cand(1)) > 1e-12
                zeroBonus = 0.25 * step;
            end
        elseif cLo <= 0 && cHi >= 0
            if ~any(abs(cand) < 1e-12)
                zeroBonus = 0.20 * step;
            end
        end

        score = pad + 0.15 * midErr + zeroBonus;
        if score < bestScore
            bestScore = score;
            bestTicks = cand;
        end
    end
end

if ~isempty(bestTicks)
    ticks = bestTicks;
    return
end

% Robust fallback: exactly nTicks and guaranteed coverage.
step = fxn_nice_step(desiredStep);
while step * (nTicks - 1) < span
    step = fxn_nice_step(step * 1.01);
    if step <= 0
        step = desiredStep;
        break
    end
end
start = floor(cLo / step) * step;
ticks = start + (0:nTicks-1) * step;
end

function step = fxn_nice_step(rawStep)
if ~isfinite(rawStep) || rawStep <= 0
    step = 1;
    return
end

pow10 = 10 ^ floor(log10(rawStep));
base = rawStep / pow10;

if base <= 1
    niceBase = 1;
elseif base <= 2
    niceBase = 2;
elseif base <= 2.5
    niceBase = 2.5;
elseif base <= 5
    niceBase = 5;
else
    niceBase = 10;
end

step = niceBase * pow10;
end

function label = fxn_format_tick_label(x)
x = round(x, 2);
if abs(x) < 1e-12
    x = 0;
end
label = sprintf('%.2f', x);
label = regexprep(label, '([\.]\d*?)0+$', '$1');
label = regexprep(label, '\.$', '');
end

function S = fxn_perm_corr(A, B, nPerm)
% Computes the observed correlation between A and B, generates a null distribution of correlations by permuting B, and calculates confidence intervals and a two-sided p-value.
if nargin < 3 || isempty(nPerm)
    nPerm = 5000;
end

S = struct('rObs', nan, 'rPerm', nan(nPerm, 1),  'ci68', [nan nan], 'ci95', [nan nan], 'pTwoSided', nan);

% Ensure A and B are vectors of the same length, and remove any pairs with non-finite values.
A = A(:);
B = B(:);
idx = isfinite(A) & isfinite(B);
A = A(idx);
B = B(idx);

if numel(A) < 2 || numel(unique(A)) < 2 || numel(unique(B)) < 2
    return;
end

% Use centered/unit-norm vectors so correlation reduces to a dot product.
A = A - mean(A);
B = B - mean(B);
normA = norm(A);
normB = norm(B);
if normA <= eps || normB <= eps
    return;
end
A = A / normA;
B = B / normB;

S.rObs = A' * B;

% Compute permutation correlations in vectorized batches to avoid thousands of scalar corr calls.
nObs = numel(A);
batchSize = min(1000, nPerm);
iStart = 1;
while iStart <= nPerm
    iEnd = min(iStart + batchSize - 1, nPerm);
    nBatch = iEnd - iStart + 1;
    [~, permIdx] = sort(rand(nObs, nBatch), 1);
    Bperm = B(permIdx);
    S.rPerm(iStart:iEnd) = A' * Bperm;
    iStart = iEnd + 1;
end

idxPerm = isfinite(S.rPerm);
if ~any(idxPerm)
    return;
end
% Calculate confidence intervals and two-sided p-value from the null distribution.
rPermValid = S.rPerm(idxPerm);
S.ci68 = prctile(rPermValid, [16 84]);
S.ci95 = prctile(rPermValid, [2.5 97.5]);
S.pTwoSided = (sum(abs(rPermValid) >= abs(S.rObs)) + 1) / (numel(rPermValid) + 1);
end

function S = fxn_collapse_scalar_by_x_iterCI(Rin, xField, yIterField, xLevels, collapseFcn)
% Collapses scalar data with iteration-level CI across x-levels.
% Uses dual getCI calls: first on iteration arrays per record, then on
% the set of medians across records within each x-level to get mean and SEM.

S = struct('x', [], 'y', [], 'sem', [], 'lb', [], 'ub', [], 'n', []);

x = [];
y = [];
sem = [];
n = [];

for i = 1:numel(xLevels)
    idx = [Rin.(xField)] == xLevels(i);
    Rsub = Rin(idx);
    if isempty(Rsub), continue; end

    yRow = nan(numel(Rsub), 1);

    for iRow = 1:numel(Rsub)
        vals = Rsub(iRow).(yIterField);
        vals = vals(:);
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end
        % First getCI: median and CI of iteration-level values
        [yRow(iRow), ~, ~] = getCI(vals, 1, 1);
    end

    idx = isfinite(yRow);
    if ~any(idx), continue; end

    % Second getCI: apply to the set of per-record medians to get mean and SEM
    medians_goodRows = yRow(idx);
    [y_mean, ~, ~, y_sem] = getCI(medians_goodRows, 2, 1);

    x(end+1,1) = xLevels(i); %#ok<AGROW>
    y(end+1,1) = y_mean; %#ok<AGROW>
    sem(end+1,1) = y_sem; %#ok<AGROW>
    n(end+1,1) = sum(idx); %#ok<AGROW>
end

S.x = x;
S.y = y;
S.sem = sem;
S.n = n;
end


function [xAxis, yTrue, curves] = fxn_collapse_tuning_panel(Rin, vField, vLevels, domainName)
xAxis = [];
yTrue = [];
curves = struct('level', {}, 'y', {}, 'n', {});

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
    for iRowSub = 1:numel(Rsub)
        Yiter = Rsub(iRowSub).(fieldEstIter);
        if isempty(Yiter), continue; end

        for iCh = 1:nCh
            vals = Yiter(:, iCh);
            vals = vals(isfinite(vals));
            if isempty(vals), continue; end
            [yRow(iRowSub, iCh), ~, ~] = getCI(vals, 1, 1);
        end
    end

    y = nan(1, size(yRow, 2));
    for iCh = 1:size(yRow, 2)
        idxRows = isfinite(yRow(:, iCh));
        if any(idxRows)
            [y(iCh), ~, ~] = getCI(yRow(idxRows, iCh), 2, 1);
        end
    end

    curves(end+1).level = vLevels(i); %#ok<AGROW>
    curves(end).y = y;
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
    idx = isfinite(x) & isfinite(y);
    if ~any(idx), continue; end
    xMin = min(xMin, min(x(idx)));
    xMax = max(xMax, max(x(idx)));
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

    idx = isfinite(x) & isfinite(y);
    x = x(idx);
    y = y(idx);

    if numel(unique(x)) < 2, continue; end
    [x, ia] = unique(x);
    y = y(ia);

    Yall(i, :) = interp1(x, y, xGrid, 'linear', nan);
end

idxRows = any(isfinite(Yall), 2);
Yall = Yall(idxRows, :);
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

    idx = isfinite(x) & isfinite(w) & w > 0;
    if ~any(idx), continue; end
    x = x(idx);
    w = w(idx);

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
idx = isfinite(x) & isfinite(y) & isfinite(lb) & isfinite(ub);
if ~any(idx), return; end

x = x(idx);
y = y(idx);
lb = lb(idx);
ub = ub(idx);

if isfield(P, 'markerSizeVec') && ~isempty(P.markerSizeVec) && numel(P.markerSizeVec) == numel(idx)
    markerSizes = P.markerSizeVec(idx);
else
    markerSizes = P.markerSize * ones(size(x));
end

plot(x, y, '-', 'Color', colorRGB, 'LineWidth', P.lineWidth);
errLo = y - lb;
errHi = ub - y;
errorbar(x, y, errLo, errHi, 'LineStyle', 'none', 'Color', colorRGB, 'LineWidth', P.errLineWidth, 'CapSize', 0);
for iPt = 1:numel(x)
    plot(x(iPt), y(iPt), 'o',  'MarkerSize', markerSizes(iPt),  'MarkerFaceColor', 'w',  'MarkerEdgeColor', colorRGB,  'LineStyle', 'none',  'LineWidth', P.lineWidth);
end
end

function [r2, rho] = fxn_curve_fit_metrics(xRef, yRef, xPred, yPred, w)
% Interpolate prediction onto reference x grid, then compute R² and Pearson's r.
% Optional w: trial counts per bin for weighted R² / rho.
yInterp = interp1(xPred, yPred, xRef, 'linear', nan);
idx = isfinite(yRef) & isfinite(yInterp);
if sum(idx) < 2
    r2 = nan; rho = nan; return
end
yR = yRef(idx);
yP = yInterp(idx);
if nargin >= 5 && ~isempty(w) && numel(w) == numel(xRef)
    wg = w(idx); wg = wg(:); wg(~isfinite(wg) | wg < 0) = 0;
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
idx = isfinite(yRef) & isfinite(yInterp);
if sum(idx) < 1
    rmse = nan;
    return
end

if nargin >= 5 && ~isempty(w) && numel(w) == numel(xRef)
    wg = w(idx); wg = wg(:) / sum(wg(isfinite(wg)));
    rmse = sqrt(sum(wg' .* (yRef(idx) - yInterp(idx)).^2, 'omitnan'));
    assert(isscalar(rmse), 'ALERT: rmse is not a scalar but a vector!!')
else
    rmse = sqrt(mean((yRef(idx) - yInterp(idx)).^2));
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
    text(i, pctTop(i), sprintf(' %.0f%%', pctTop(i)),  'VerticalAlignment', 'bottom',  'HorizontalAlignment', 'left',  'FontSize', 9);
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
P = struct('x', [], 'y', [], 'sem', [], 'n', []);
for iSig = 1:numel(gaborCST_unik)
    for iCz = 1:numel(cSDT_unik)
        S = fxn_collapse_param_byLevel(  Rsub, pName,  collapseFcn, nBins_Part4, gaborCST_unik(iSig), cSDT_unik(iCz));
        if isempty(S.x), continue; end
        P.x = [P.x; S.x(:)]; %#ok<AGROW>
        P.y = [P.y; S.y(:)]; %#ok<AGROW>
        P.sem = [P.sem; S.sem(:)]; %#ok<AGROW>
        P.n = [P.n; S.n(:)]; %#ok<AGROW>
    end
end
end

function [rho, slope, bias, rmse] = fxn_param_recovery_stats(x, y)
idx = isfinite(x) & isfinite(y);
x = x(idx);
y = y(idx);
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
% Uses dual getCI calls: first on per-record estimates to get medians,
% then on the set of medians to get mean and SEM.

S = struct('x', [], 'y', [], 'sem', [], 'n', []);

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
end

xvals_all  = [Rsub.(xField)]';
yvals_all  = [Rsub.(yField)]';

idx = isfinite(xvals_all) & isfinite(yvals_all);
xvals_all  = xvals_all(idx);
yvals_all  = yvals_all(idx);

% group by the examined parameter's true levels
xLevels = unique(xvals_all);

x_out  = nan(numel(xLevels), 1);
y_out  = nan(numel(xLevels), 1);
sem_out = nan(numel(xLevels), 1);
n_out  = nan(numel(xLevels), 1);

for iLev = 1:numel(xLevels)
    idxLev = xvals_all == xLevels(iLev);

    % Apply second getCI to the set of medians within this level
    % to get mean and SEM (via 68% CI)
    medians_thisLev = yvals_all(idxLev);
    [y_mean, ~, ~, y_sem] = getCI(medians_thisLev, 2, 1);

    x_out(iLev)  = xLevels(iLev);
    y_out(iLev)  = y_mean;
    sem_out(iLev) = y_sem;
    n_out(iLev)  = sum(idxLev);
end

S.x  = x_out;
S.y  = y_out;
S.sem = sem_out;
S.n  = n_out;
end

function S = fxn_bin_criterion_values(Rsub, nBins_Part4, collapseFcn)
% Bin criterion values and aggregate using dual getCI calls.
% First getCI: per-record estimates → medians
% Second getCI: per-bin medians → mean and SEM

S = struct('x', [], 'y', [], 'sem', [], 'n', []);

xvals_all  = [Rsub.criterion_DV_true]';
yvals_all  = [Rsub.criterion_DV_est_med]';

idx = isfinite(xvals_all) & isfinite(yvals_all);
xvals_all  = xvals_all(idx);
yvals_all  = yvals_all(idx);

xMin = min(xvals_all);
xMax = max(xvals_all);

if xMin == xMax
    % Only one unique value: apply getCI to all medians
    [y_mean, ~, ~, y_sem] = getCI(yvals_all, 2, 1);
    S.x = round(xMin);
    S.y = y_mean;
    S.sem = y_sem;
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
sem_out = nan(numel(uCtrs), 1);
n_out  = nan(numel(uCtrs), 1);

for iCtr = 1:numel(uCtrs)
    idxCtr = ismember(binIdx, find(binCtrs == uCtrs(iCtr)));

    % Apply second getCI to the set of medians within this bin
    medians_thisBin = yvals_all(idxCtr);
    [y_mean, ~, ~, y_sem] = getCI(medians_thisBin, 2, 1);

    x_out(iCtr)  = uCtrs(iCtr);
    y_out(iCtr)  = y_mean;
    sem_out(iCtr) = y_sem;
    n_out(iCtr)  = sum(idxCtr);
end

keep = isfinite(x_out) & isfinite(y_out);
S.x  = x_out(keep);
S.y  = y_out(keep);
S.sem = sem_out(keep);
S.n  = n_out(keep);
end
