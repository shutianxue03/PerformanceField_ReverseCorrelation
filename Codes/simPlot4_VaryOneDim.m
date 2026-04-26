% Visualization of model-simulation outputs
% Assumes compiled struct R already exists in workspace.
%
% Figures:
%   Fig 1: template GoF
%   Fig 2: recovered tuning functions
%   Fig 3: metric recovery curves
%   Fig 4: parameter recovery scatter

clc; close all;
%% Shared plotting formats
iModelA_fit = 1 % 1=fitting using data-derived template; 2=true template
nBins_Part4 = 3; % define bins for collapsing parameter recovery points; use 3 for main text, 5 for Supp

% Load data
load(sprintf('%s/Outputs/R_A%d.mat', nameFolder_server, iModelA_fit), 'R')
setting = struct();

% ---------- filtering ----------
setting.flag_use_filtered = true;

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
setting.scalarCollapseFcn = @mean;   % mean across rows after filtering
setting.curveCollapseFcn  = @mean;   % mean across rows after collapsing
setting.curveNGrid = 9;              % redefine bins after collapsing

% ---------- figure-specific fixed choices ----------
setting.fig12_Bsim = 1;
setting.fig12_Bfit = 1;

% ---------- field names ----------
setting.fig1_mode_fields = {'nBasisORI_mode', 'nBasisSF_mode', 'basisFxnORI_mode', 'basisFxnSF_mode', 'ridge_mode'};
setting.fig1_mode_labels = {'ORI basis mode', 'SF basis mode', 'ORI basis family mode', 'SF basis family mode', 'Ridge mode'};

setting.fig1_pct_fields = {'nBasisORI_mode_pct', 'nBasisSF_mode_pct', 'basisFxnORI_mode_pct', 'basisFxnSF_mode_pct', 'ridge_mode_pct'};
setting.fig1_pct_labels = {'ORI mode selection rate', 'SF mode selection rate', 'ORI basis family mode selection rate', 'SF basis family mode selection rate', 'Ridge mode selection rate'};

setting.varFields = {'Nmul_true', 'Nadd_true', 'Nshared_true', 'gaborCST', 'Cz_true', 'lambda_whiten', 'C_contribution'};
setting.varNames  = {'Nmul', 'Nadd', 'Nshared', 'signalCST', 'Cz', 'lambda_white', 'C_contribution'};

setting.fig1_varFields = setting.varFields;
setting.fig1_varNames  = setting.varNames;

setting.fig2_varFields = setting.varFields;
setting.fig2_varNames  = setting.varNames;

setting.metricFields = {'pYES', 'pC', 'pA'};
setting.metricCurveDataFields = {'pYES_data_curve_med', 'pC_data_curve_med', 'pA_data_curve_med'};
setting.metricCurvePredFields = {'pYES_pred_curve_med', 'pC_pred_curve_med', 'pA_pred_curve_med'};

setting.fig4_paramNames = {'Nmul', 'Nadd', 'Nshared', 'criterion_DV'};
setting.fig4_paramTitles = {'Nmul', 'Nadd', 'Nshared', 'criterion DV'};

% ---------- unique levels ----------
gaborCST_unik = unique([R.gaborCST]);
Cz_unik = unique([R.Cz_true]);
Nmul_unik = unique([R.Nmul_true]);
Nadd_unik = unique([R.Nadd_true]);
Nshared_unik = unique([R.Nshared_true]);
Bsim_unik = unique([R.iModelB_sim]);
Bfit_unik = unique([R.iModelB_fit]);
lambda_whiten_unik = unique([R.lambda_whiten]);
C_contribution_unik = unique([R.C_contribution]);

% Expand plotting styles to match however many levels are present.
base_cmap_Cz = setting.cmap_Cz;
if numel(Cz_unik) <= size(base_cmap_Cz, 1)
    setting.cmap_Cz = base_cmap_Cz(1:numel(Cz_unik), :);
else
    nExtraCz = numel(Cz_unik) - size(base_cmap_Cz, 1);
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

%% Figure 1: Basis selection
fprintf('\n%s: Fig 1: Basis selection rates\n', string(datetime('now')))

R_fig1 = R([R.iModelB_sim] == setting.fig12_Bsim & [R.iModelB_fit] == setting.fig12_Bfit);

% group labels: signalCST × Cz, one per entry in R_fig1
grpLabels_fig1 = arrayfun(@(r) sprintf('sig=%.3g x Cz=%.3g', r.gaborCST, r.Cz_true), ...
    R_fig1, 'UniformOutput', false);

% grpLabels_fig1 = arrayfun(@(r) sprintf('Nmul=%g x Nadd=%g x Nshared=%g', r.Nmul_true, r.Nadd_true, r.Nshared_true), ...
%     R_fig1, 'UniformOutput', false);

uGrp_fig1    = unique(grpLabels_fig1);
grpCmap_fig1 = setting.palette(1:numel(uGrp_fig1), :);

h = figure('Position', [100 100 1600 800]);

subplot(2,3,1); hold on;
%=====================%
plot_ranked_categorical({R_fig1.(setting.fig1_mode_fields{1})}, 'nBasisORI', grpLabels_fig1, grpCmap_fig1, true);
%=====================%
ylabel('% iterations selected');
xlabel('# ORI basis functions');
%=====================%
fxn_style_ax(gca, setting);
%=====================%

subplot(2,3,2); hold on;
%=====================%
plot_ranked_categorical({R_fig1.(setting.fig1_mode_fields{2})}, 'nBasisSF', grpLabels_fig1, grpCmap_fig1, false);
%=====================%
ylabel('% iterations selected');
xlabel('# SF basis functions');
%=====================%
fxn_style_ax(gca, setting);
%=====================%

subplot(2,3,6); hold on;
%=====================%
plot_ranked_categorical({R_fig1.(setting.fig1_mode_fields{5})}, 'ridge', grpLabels_fig1, grpCmap_fig1, false);
%=====================%
ylabel('% iterations selected');
xlabel('Ridge penalty');
%=====================%
fxn_style_ax(gca, setting);
%=====================%

subplot(2,3,4); hold on;
%=====================%
plot_ranked_categorical({R_fig1.(setting.fig1_mode_fields{3})}, 'basisFxnORI', grpLabels_fig1, grpCmap_fig1, false);
%=====================%
ylabel('% iterations selected');
xlabel('ORI basis family');
%=====================%
fxn_style_ax(gca, setting);
%=====================%

subplot(2,3,5); hold on;
%=====================%
plot_ranked_categorical({R_fig1.(setting.fig1_mode_fields{4})}, 'basisFxnSF', grpLabels_fig1, grpCmap_fig1, false);
%=====================%
ylabel('% iterations selected');
xlabel('SF basis family');
%=====================%
fxn_style_ax(gca, setting);
%=====================%

sgtitle('Fig 1: rank of selected basis settings across iterations', 'FontWeight', 'bold');

saveas(h, fullfile(nameFolder_Figures_part4, 'FigS1_selectedBasisRank.png'));
close(h);

%% Figure 1B: Basis selection co-occurrence
% fprintf('\n%s: Fig 1B\n', string(datetime('now')))

% h = figure('Position', [100 100 1500 900]);
% tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% cooccPairs = {
%     'nBasisORI_mode', 'basisFxnORI_mode', 'nBasisORI', 'basisFxnORI';
%     'nBasisSF_mode',  'basisFxnSF_mode',  'nBasisSF',  'basisFxnSF';
%     'nBasisORI_mode', 'ridge_mode',       'nBasisORI', 'ridge';
%     'nBasisSF_mode',  'ridge_mode',       'nBasisSF',  'ridge';
%     'basisFxnORI_mode', 'ridge_mode',     'basisFxnORI', 'ridge';
%     'basisFxnSF_mode',  'ridge_mode',     'basisFxnSF',  'ridge'};

% for iPair = 1:size(cooccPairs, 1)
%     nexttile; hold on;
%     %=====================%
%     valsA_cc = fxn_field_values_to_string(R_fig1, cooccPairs{iPair, 1});
%     valsB_cc = fxn_field_values_to_string(R_fig1, cooccPairs{iPair, 2});
%     uA_cc = unique(valsA_cc);
%     uB_cc = unique(valsB_cc);
%     M_cc = zeros(numel(uB_cc), numel(uA_cc));
%     for iA_cc = 1:numel(uA_cc)
%         for iB_cc = 1:numel(uB_cc)
%             M_cc(iB_cc, iA_cc) = sum(valsA_cc == uA_cc(iA_cc) & valsB_cc == uB_cc(iB_cc));
%         end
%     end
%     if sum(M_cc(:)) > 0
%         M_cc = M_cc / sum(M_cc(:)) * 100;
%     end
%     imagesc(1:numel(uA_cc), 1:numel(uB_cc), M_cc);
%     set(gca, 'YDir', 'normal', ...
%         'XTick', 1:numel(uA_cc), 'XTickLabel', cellstr(uA_cc), ...
%         'YTick', 1:numel(uB_cc), 'YTickLabel', cellstr(uB_cc));
%     xtickangle(45);
%     xlabel(cooccPairs{iPair, 3});
%     ylabel(cooccPairs{iPair, 4});
%     title(sprintf('%s x %s', cooccPairs{iPair, 3}, cooccPairs{iPair, 4}));
%     colorbar;
%     %=====================%
%     fxn_style_ax(gca, setting);
%     %=====================%
% end

% sgtitle(sprintf('Fig 1B: basis-setting co-occurrence | filtered | Bsim = %d | Bfit = %d', ...
%     setting.fig12_Bsim, setting.fig12_Bfit), 'FontWeight', 'bold');
% saveas(h, fullfile(nameFolder_Figures_part4, 'FigS1B_basisCooccurrence.png'));
% close(h);

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
    xlabel(vName); ylabel('Template RMSE'); ylim([0 0.3]); title(vName);

    % ------ Rows 2 & 3: ORI then SF tuning ------
    for iRow = 1:2
        nexttile(iRow*nVars + iCol); hold on;

        yTrue = fig2Cache(iCol).tuning(iRow).yTrue;
        curves = fig2Cache(iCol).tuning(iRow).curves;
        if ~isempty(yTrue)
            hTrue = plot(axis_tuning{iRow}, yTrue, 'r-', 'LineWidth', setting.trueWidth);
            yTrueNorm = yTrue / max(yTrue);
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
                legH(nLeg) = plot(axis_tuning{iRow}, curves(iLev).y, 'LineWidth', setting.lineWidth, 'Color', cLev);
                rmse = fxn_curve_rmse(xIdx, yTrueNorm, xIdx, curves(iLev).y);
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

h = figure('Position', [100 100 1300 900]);
tiledlayout(numel(Bsim_unik), numel(setting.metricFields), 'TileSpacing', 'loose', 'Padding', 'loose');

for iCol_Bsim = 1:numel(Bsim_unik)
    simModel = Bsim_unik(iCol_Bsim);

    % simulated data are duplicated across Bfit; use Bfit = 1
    R_data = R([R.iModelB_sim] == simModel & [R.iModelB_fit] == 1);

    for iCol_metric = 1:numel(setting.metricFields)
        mName = setting.metricFields{iCol_metric};
        dataField = setting.metricCurveDataFields{iCol_metric};
        predField = setting.metricCurvePredFields{iCol_metric};

        nexttile; hold on;
        isFirstPanel = (iCol_Bsim == 1) && (iCol_metric == 1);

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
        if iCol_Bsim == numel(Bsim_unik)
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
    end % iCol_metric
end % iCol_Bsim

sgtitle('Fig 3: metric recovery | filtered', 'FontWeight', 'bold');
saveas(h, fullfile(nameFolder_Figures_part4, 'FigS3_metricRecovery.png'));
close(h);

%% Figure 4: parameter recovery
%    Only matched pairs Bsim = Bfit

fprintf('\n%s: Fig 4: Parameter recovery\n', string(datetime('now')))

setting_fig4 = setting;
setting_fig4.fontSize = setting.fontSize;

% Condition coding in Fig 4: color = C_contribution, marker = lambda_white.
cmap_cont_fig4 = setting.palette(1:numel(C_contribution_unik), :);
marker_lambda_fig4 = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h', 'x', '+'};
if numel(lambda_whiten_unik) > numel(marker_lambda_fig4)
    nRep_mk = ceil(numel(lambda_whiten_unik) / numel(marker_lambda_fig4));
    marker_lambda_fig4 = repmat(marker_lambda_fig4, 1, nRep_mk);
end
marker_lambda_fig4 = marker_lambda_fig4(1:numel(lambda_whiten_unik));

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

        % one series per lambda_white x C_contribution, preserving parameter levels
        for iLW = 1:numel(lambda_whiten_unik)
            for iCont = 1:numel(C_contribution_unik)

                S = fxn_collapse_param_byLevel(Rsub, pName, lambda_whiten_unik(iLW), C_contribution_unik(iCont), setting.scalarCollapseFcn, nBins_Part4);

                if isempty(S.x), continue; end

                % optional: connect points from the same signalCST × Cz group
                good = isfinite(S.x) & isfinite(S.y);
                if any(good)
                    plot(S.x(good), S.y(good), '-', ...
                        'Color', cmap_cont_fig4(iCont,:), ...
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
                        Pplot.markerSize = setting.fig4_markerSizeMin + ...
                            (setting.fig4_markerSizeMax - setting.fig4_markerSizeMin) * (sc - 1) / max(sc, 1);
                        Pplot.markerSize = min(max(Pplot.markerSize, setting.fig4_markerSizeMin), setting.fig4_markerSizeMax);
                    end

                    %=====================%
                    line([S.x(iPt) S.x(iPt)], [S.lb(iPt) S.ub(iPt)], ...
                        'Color', cmap_cont_fig4(iCont,:), 'LineWidth', Pplot.errLineWidth);
                    plot(S.x(iPt), S.y(iPt), marker_lambda_fig4{iLW}, ...
                        'MarkerSize', Pplot.markerSize, ...
                        'MarkerFaceColor', 'w', ...
                        'MarkerEdgeColor', cmap_cont_fig4(iCont,:), ...
                        'LineStyle', 'none', ...
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
        Ppanel = fxn_collect_fig4_panel_points(Rsub, pName, lambda_whiten_unik, C_contribution_unik, setting.scalarCollapseFcn, nBins_Part4);
        [rho_panel, ~, ~, rmse_panel] = fxn_param_recovery_stats(Ppanel.x, Ppanel.y);
        title(sprintf('Bfit=Bsim=%d \nr=%.2f | RMSE=%.3f', Bsimfit, rho_panel, rmse_panel));
        ax_fig4_last = gca;
    end % iCol_Bsimfit
end % iRow_param

sgtitle('Fig 4: Parameter recovery | Matched model pairs', 'FontWeight', 'bold');

if isgraphics(ax_fig4_last)
    pos_fig4leg = ax_fig4_last.Position;
    axColor_fig4 = axes('Position', pos_fig4leg, 'Color', 'none', 'Visible', 'off', ...
        'HitTest', 'off', 'HandleVisibility', 'off');
    hold(axColor_fig4, 'on');
    hCont_fig4 = gobjects(numel(C_contribution_unik), 1);
    for i_leg = 1:numel(C_contribution_unik)
        hCont_fig4(i_leg) = plot(axColor_fig4, nan, nan, '-', 'Color', cmap_cont_fig4(i_leg,:), 'LineWidth', 2);
    end
    leg1_fig4 = legend(axColor_fig4, hCont_fig4, ...
        arrayfun(@(x) sprintf('%g', x), C_contribution_unik, 'UniformOutput', false), ...
        'Location', 'northeast', 'Orientation', 'vertical', 'Box', 'on', ...
        'FontSize', max(ax_fig4_last.FontSize - 1, 7));
    title(leg1_fig4, 'C_contribution');
    axSig_fig4 = axes('Position', pos_fig4leg, 'Color', 'none', 'Visible', 'off', ...
        'HitTest', 'off', 'HandleVisibility', 'off');
    hold(axSig_fig4, 'on');
    hSig_fig4 = gobjects(numel(lambda_whiten_unik), 1);
    for i_leg = 1:numel(lambda_whiten_unik)
        hSig_fig4(i_leg) = plot(axSig_fig4, nan, nan, ...
            'LineStyle', 'none', ...
            'Marker', marker_lambda_fig4{i_leg}, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', 'k');
    end
    leg2_fig4 = legend(axSig_fig4, hSig_fig4, ...
        arrayfun(@(x) sprintf('%g', x), lambda_whiten_unik, 'UniformOutput', false), ...
        'Location', 'southeast', 'Orientation', 'vertical', 'Box', 'on', ...
        'FontSize', max(ax_fig4_last.FontSize - 1, 7));
    title(leg2_fig4, 'lambda_white');
end
saveas(h, fullfile(nameFolder_Figures_part4, 'FigS4_paramRecovery.png'));
close(h);

%% Figure 5: model recovery
fprintf('\n%s: Fig 5: Model recovery\n', string(datetime('now')))

W_bestNLL = nan(numel(Bfit_unik), numel(Bsim_unik));
if ~isempty(R)
    if isfield(R, 'nameIO')
        datasetIDs_wr = string({R.nameIO});
    else
        datasetIDs_wr = string(arrayfun(@(r) sprintf('Nm%g_Na%g_Ns%g_G%g_Cz%g', ...
            r.Nmul_true, r.Nadd_true, r.Nshared_true, r.gaborCST, r.Cz_true), R, 'UniformOutput', false));
    end
    for iSim_wr = 1:numel(Bsim_unik)
        simModel_wr = Bsim_unik(iSim_wr);
        idxSim_wr   = [R.iModelB_sim] == simModel_wr;
        idsSim_wr   = unique(datasetIDs_wr(idxSim_wr));
        wins_wr     = zeros(numel(Bfit_unik), 1);
        nValid_wr   = 0;
        for iID_wr = 1:numel(idsSim_wr)
            idxID_wr = idxSim_wr & datasetIDs_wr == idsSim_wr(iID_wr);
            Rsub_wr  = R(idxID_wr);
            if isempty(Rsub_wr), continue; end
            nllVals_wr = nan(numel(Bfit_unik), 1);
            for iFit_wr = 1:numel(Bfit_unik)
                idxFit_wr = [Rsub_wr.iModelB_fit] == Bfit_unik(iFit_wr);
                vals_wr   = [Rsub_wr(idxFit_wr).nLL_med]';
                vals_wr   = vals_wr(isfinite(vals_wr));
                if ~isempty(vals_wr)
                    nllVals_wr(iFit_wr) = mean(vals_wr, 'omitnan');
                end
            end
            if ~any(isfinite(nllVals_wr)), continue; end
            [~, idxBest_wr] = min(nllVals_wr);
            wins_wr(idxBest_wr) = wins_wr(idxBest_wr) + 1;
            nValid_wr = nValid_wr + 1;
        end
        if nValid_wr > 0
            W_bestNLL(:, iSim_wr) = wins_wr / nValid_wr;
        end
    end
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
        datasetIDs_dn = string(arrayfun(@(r) sprintf('Nm%g_Na%g_Ns%g_G%g_Cz%g', ...
            r.Nmul_true, r.Nadd_true, r.Nshared_true, r.gaborCST, r.Cz_true), R, 'UniformOutput', false));
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
    W_bestNLL,   'Best model win rate',      [0 1]; ...
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

    % Increase spacing by shrinking each axes box within its tile.
    axPos = get(gca, 'Position');
    shrink = 0.86;
    newW = axPos(3) * shrink;
    newH = axPos(4) * shrink;
    newX = axPos(1) + (axPos(3) - newW) / 2;
    newY = axPos(2) + (axPos(4) - newH) / 2;
    set(gca, 'Position', [newX, newY, newW, newH]);

    % color axis
    finVals = M_plot(isfinite(M_plot));
    if ~isempty(clim_fix) && numel(clim_fix) == 2
        cLo = clim_fix(1); cHi = clim_fix(2);
    elseif ~isempty(finVals)
        cLo = min(finVals); cHi = max(finVals);
        if cLo == cHi, cLo = cLo - eps; cHi = cHi + eps; end
    else
        cLo = 0; cHi = 1;
    end
    clim([cLo cHi]);
    % colorbar hidden for cleaner Figure 5 panels

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
    estFields = cell(size(pNames_all));
    for iP = 1:numel(pNames_all)
        if strcmp(pNames_all{iP}, 'criterion_DV')
            estFields{iP} = 'criterion_DV_est_med';
        else
            estFields{iP} = sprintf('%s_est_med', pNames_all{iP});
        end
    end

    keepParam = false(size(pNames_all));
    for iP = 1:numel(pNames_all)
        if ~isfield(Rsub, estFields{iP}), continue; end
        vals = [Rsub.(estFields{iP})]';
        keepParam(iP) = any(isfinite(vals));
    end

    pNames = pNames_all(keepParam);
    estFields = estFields(keepParam);
    nParam = numel(pNames);
    if nParam < 2, continue; end

    nObs = numel(Rsub);
    paramVals = nan(nObs, nParam);
    for iP = 1:nParam
        paramVals(:, iP) = [Rsub.(estFields{iP})]';
    end

    sigVals_all = [Rsub.gaborCST]';
    czVals_all = [Rsub.Cz_true]';
    [isSig, sigIdx_all] = ismember(sigVals_all, gaborCST_unik);
    [isCz, czIdx_all] = ismember(czVals_all, Cz_unik);
    validGroupBase = isfinite(sigVals_all) & isfinite(czVals_all) & isSig & isCz;

    groupRows = cell(numel(gaborCST_unik), numel(Cz_unik));
    for iSig = 1:numel(gaborCST_unik)
        for iCz = 1:numel(Cz_unik)
            groupRows{iSig, iCz} = find(validGroupBase & sigIdx_all == iSig & czIdx_all == iCz);
        end
    end

    pairIdx = nchoosek(1:nParam, 2);
    nPairs = size(pairIdx, 1);
    nCols = min(3, nPairs);
    nRows = ceil(nPairs / nCols);
    czLineStyles = {'-', '--', ':', '-.'};
    markerIsCircle = contains(setting.marker_signal, 'o');

    h = figure('Position', [100 100 420*nCols 320*nRows]);
    tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for iPair = 1:nPairs
        iP1 = pairIdx(iPair, 1);
        iP2 = pairIdx(iPair, 2);

        xAll = paramVals(:, iP1);
        yAll = paramVals(:, iP2);
        good = validGroupBase & isfinite(xAll) & isfinite(yAll);
        x = xAll(good);
        y = yAll(good);

        nexttile; hold on;
        if isempty(x)
            text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
            xlim([0 1]); ylim([0 1]);
            title(sprintf('%s vs %s', pNames{iP1}, pNames{iP2}));
            fxn_style_ax(gca, setting);
            continue;
        end

        for iSig = 1:numel(gaborCST_unik)
            for iCz = 1:numel(Cz_unik)
                idxRows = groupRows{iSig, iCz};
                if isempty(idxRows), continue; end

                idxRows = idxRows(isfinite(xAll(idxRows)) & isfinite(yAll(idxRows)));
                if isempty(idxRows), continue; end

                xGrp = xAll(idxRows);
                yGrp = yAll(idxRows);

                % Plot individual data points
                plot(xGrp, yGrp, setting.marker_signal{iSig}, ...
                    'LineStyle', 'none', ...
                    'MarkerSize', max(setting.markerSize, 6), ...
                    'MarkerFaceColor', setting.cmap_Cz(iCz,:), ...
                    'MarkerEdgeColor', 'w', ...
                    'LineWidth', 1.0);

                % Fit and plot a linear regression line for the group
                if numel(xGrp) >= 2
                    pfit_grp = polyfit(xGrp, yGrp, 1);
                    xfit_grp = linspace(min(xGrp), max(xGrp), 60);
                    yfit_grp = polyval(pfit_grp, xfit_grp);
                    lineStyle_grp = czLineStyles{mod(iCz-1, numel(czLineStyles)) + 1};
                    if markerIsCircle(iSig)
                        lineStyle_grp = '-';
                    else
                        lineStyle_grp = '--';
                    end
                    plot(xfit_grp, yfit_grp, ...
                        'Color', setting.cmap_Cz(iCz,:), ...
                        'LineStyle', lineStyle_grp, ...
                        'LineWidth', 1.4);
                end
            end
        end

        if numel(x) >= 2
            [r, p] = corr(x, y);
            title(sprintf('%s vs %s | r=%.2f, p=%.3g', pNames{iP1}, pNames{iP2}, r, p));
        else
            title(sprintf('%s vs %s', pNames{iP1}, pNames{iP2}));
        end

        axis square;
        xlabel(sprintf('Est. %s', pNames{iP1}));
        ylabel(sprintf('Est. %s', pNames{iP2}));
        fxn_style_ax(gca, setting);
    end

    % Figure-level legends (color=Cz, marker=signalCST)
    axBase = gca;
    pos = axBase.Position;

    axColor = axes('Position', pos, 'Color', 'none', 'Visible', 'off', ...
        'HitTest', 'off', 'HandleVisibility', 'off');
    hold(axColor, 'on');
    hCz = gobjects(numel(Cz_unik), 1);
    for iCz = 1:numel(Cz_unik)
        hCz(iCz) = plot(axColor, nan, nan, 'o', ...
            'LineStyle', 'none', ...
            'MarkerFaceColor', setting.cmap_Cz(iCz,:), ...
            'MarkerEdgeColor', setting.cmap_Cz(iCz,:), ...
            'MarkerSize', 7);
    end
    leg1 = legend(axColor, hCz, ...
        arrayfun(@(x) sprintf('%g', x), Cz_unik, 'UniformOutput', false), ...
        'Location', 'northeast', 'Box', 'on', 'FontSize', max(setting.fontSize - 1, 8));
    title(leg1, 'Cz');

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

%% Local helpers  (only functions used 2+ times)
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

function P = fxn_collect_fig4_panel_points(Rsub, pName, lambda_whiten_unik, C_contribution_unik, collapseFcn, nBins_Part4)
P = struct('x', [], 'y', [], 'lb', [], 'ub', [], 'n', []);
for iLW = 1:numel(lambda_whiten_unik)
    for iCont = 1:numel(C_contribution_unik)
        S = fxn_collapse_param_byLevel(Rsub, pName, lambda_whiten_unik(iLW), C_contribution_unik(iCont), collapseFcn, nBins_Part4);
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

function S = fxn_collapse_param_byLevel(Rin, pName, lambdaVal, contributionVal, collapseFcn, nBins_Part4)
% Keep the examined parameter levels separate.
% Collapse only across the other parameters.

S = struct('x', [], 'y', [], 'lb', [], 'ub', [], 'n', []);

% restrict to one lambda_whiten x C_contribution combination
idx = [Rin.lambda_whiten] == lambdaVal & [Rin.C_contribution] == contributionVal;
Rsub = Rin(idx);
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
