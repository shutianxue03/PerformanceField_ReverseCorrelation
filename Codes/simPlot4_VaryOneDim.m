% Visualization of model-simulation outputs
% Assumes compiled struct R already exists in workspace.
%
% Figures:
%   Fig 1: template GoF
%   Fig 2: recovered tuning functions
%   Fig 3: metric recovery curves
%   Fig 4: parameter recovery scatter

clc; close all;
iModelA_fit=1 % 1=fitting using data-derived template; 2=true template
load(sprintf('%s/Outputs/R_A%d.mat', nameFolder_server, iModelA_fit), 'R')

%% Shared plotting formats

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
setting.cmap_Cz = [1 0 0; 0 0 0; 0 0 1];   % low/mid/high Cz
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

% ---------- collapse rules ----------
setting.scalarCollapseFcn = @mean;   % mean across rows after filtering
setting.curveCollapseFcn  = @mean;   % mean across rows after collapsing
setting.curveNGrid = 9;              % redefine bins after collapsing

% ---------- figure-specific fixed choices ----------
setting.fig12_Bsim = 1;
setting.fig12_Bfit = 1;

% ---------- field names ----------
setting.fig0_mode_fields = {'nBasisORI_mode', 'nBasisSF_mode'};
setting.fig0_mode_labels = {'ORI basis mode', 'SF basis mode'};

setting.fig0_pct_fields = {'nBasisORI_mode_pct', 'nBasisSF_mode_pct'};
setting.fig0_pct_labels = {'ORI mode selection rate', 'SF mode selection rate'};

setting.fig1_yfields = {'template_rmse', 'template_R2'};
setting.fig1_ylabels = {'Template RMSE', 'Template R^2'};

setting.varFields = {'Nmul_true', 'Nadd_true', 'Nshared_true', 'gaborCST', 'Cz_true'};
setting.varNames  = {'Nmul', 'Nadd', 'Nshared', 'signalCST', 'Cz'};

setting.fig0_varFields = setting.varFields;
setting.fig0_varNames = setting.varNames;

setting.fig1_varFields = setting.varFields;
setting.fig1_varNames = setting.varNames;

setting.fig2_varFields = setting.varFields;
setting.fig2_varNames  = setting.varNames;

setting.metricFields = {'pYES', 'pC', 'pA'};
setting.metricCurveDataFields = {'pYES_data_curve_med', 'pC_data_curve_med', 'pA_data_curve_med'};
setting.metricCurvePredFields = {'pYES_pred_curve_med', 'pC_pred_curve_med', 'pA_pred_curve_med'};

setting.fig4_paramNames = {'Nmul', 'Nadd', 'Nshared', 'criterion_DV'};
setting.fig4_paramTitles = {'Nmul', 'Nadd', 'Nshared', 'criterion DV'};

%% ---------------- Figure 0A: Basis counts ----------------
fprintf('\n%s: Fig 0A\n', string(datetime('now')))

R_fig0 = R([R.iModelB_sim] == setting.fig12_Bsim & [R.iModelB_fit] == setting.fig12_Bfit);

h = figure('Position', [100 100 1800 650]);
tiledlayout(2, numel(setting.fig0_varFields), 'TileSpacing', 'compact', 'Padding', 'compact');

for iCol = 1:numel(setting.fig0_varFields)
    xField = setting.fig0_varFields{iCol};
    xName  = setting.fig0_varNames{iCol};
    xLevels = unique([R_fig0.(xField)]);

    for iRow=1:2 % 1=ORI; 2=SF
        nexttile((iRow-1)*numel(setting.fig0_varFields) + iCol); hold on;

        S = fxn_collapse_scalar_by_x_simple(R_fig0, ...
            xField, setting.fig0_mode_fields{iRow}, ...
            xLevels, setting.scalarCollapseFcn);

        if ~isempty(S.x)
            fxn_plot_line_with_band(S.x, S.y, S.lb, S.ub, [0 0 0], '-', setting);
        end

        box on;
        set(gca, 'FontSize', setting.fontSize, 'LineWidth', setting.axisLineWidth);
        xlabel(xName);
        ylabel(setting.fig0_mode_labels{iRow});
        yticks([2:9]), ylim([2,9])
        title(xName);
    end % iRow
end % iItem

sgtitle(sprintf('Fig 0A: modal basis counts | filtered | Bsim = %d | Bfit = %d', ...
    setting.fig12_Bsim, setting.fig12_Bfit), 'FontWeight', 'bold');

saveas(h, fullfile(nameFolder_Figures_part4, 'Fig0A_basisMode.png'));
close(h);

%% ---------------- Figure 0B: Basis counts percentage ----------------
fprintf('\n%s: Fig 0B\n', string(datetime('now')))

h = figure('Position', [100 100 1800 650]);
tiledlayout(2, numel(setting.fig0_varFields), 'TileSpacing', 'compact', 'Padding', 'compact');

for iCol = 1:numel(setting.fig0_varFields)
    xField = setting.fig0_varFields{iCol};
    xName  = setting.fig0_varNames{iCol};
    xLevels = unique([R_fig0.(xField)]);

    for iRow=1:2
        % ----- row 1: nBasisORI_mode_pct -----
        nexttile((iRow-1)*numel(setting.fig0_varFields) + iCol); hold on;

        S = fxn_collapse_scalar_by_x_simple(R_fig0, ...
            xField, setting.fig0_pct_fields{iRow}, ...
            xLevels, setting.scalarCollapseFcn);

        if ~isempty(S.x)
            fxn_plot_line_with_band(S.x, S.y, S.lb, S.ub, [0 0 0], '-', setting);
        end

        box on;
        set(gca, 'FontSize', setting.fontSize, 'LineWidth', setting.axisLineWidth);
        xlabel(xName);
        ylabel(setting.fig0_pct_labels{iRow});
        ylim([0 1]);
        title(xName);
        ytickformat('percentage')
    end % iRow
end % iCol

sgtitle(sprintf('Fig 0B: modal basis selection rate | filtered | Bsim = %d | Bfit = %d', ...
    setting.fig12_Bsim, setting.fig12_Bfit), 'FontWeight', 'bold');

saveas(h, fullfile(nameFolder_Figures_part4, 'Fig0B_basisModePct.png'));
close(h);

%% Figure 1: template GoF
%    Use one representative Bsim/Bfit because these are file-level

fprintf('\n%s: Fig 1\n', string(datetime('now')))

R_fig1 = R([R.iModelB_sim] == setting.fig12_Bsim & [R.iModelB_fit] == setting.fig12_Bfit);

h = figure('Position', [100 100 1800 650]);
tiledlayout(2, numel(setting.fig1_varFields), 'TileSpacing', 'compact', 'Padding', 'compact');

for iCol = 1:numel(setting.fig1_varFields)
    xField = setting.fig1_varFields{iCol};
    xName  = setting.fig1_varNames{iCol};
    xLevels = unique([R_fig1.(xField)]);

    for iRow = 1:2 % 1=RMSE; 2=R2
        % ----- row 1: RMSE -----
        nexttile((iRow-1)*numel(setting.fig1_varFields) + iCol); hold on;

        S = fxn_collapse_scalar_by_x_simple(R_fig1, ...
            xField, setting.fig1_yfields{iRow}, ...
            xLevels, setting.scalarCollapseFcn);

        if ~isempty(S.x)
            fxn_plot_line_with_band(S.x, S.y, S.lb, S.ub, [0 0 0], '-', setting);
        end

        box on;
        set(gca, 'FontSize', setting.fontSize, 'LineWidth', setting.axisLineWidth);
        xlabel(xName);
        ylabel(setting.fig1_ylabels{iRow});
        ylim([0, .5]*iRow)
        title(xName);

    end % iRow
end % iCol

sgtitle(sprintf('Fig 1: template GoF | filtered | Bsim = %d | Bfit = %d', ...
    setting.fig12_Bsim, setting.fig12_Bfit), 'FontWeight', 'bold');

saveas(h, fullfile(nameFolder_Figures_part4, 'Fig1_templateGoF.png'));
close(h);

%% Figure 2: recovered tuning function
%    Use one representative Bsim/Bfit because these are file-level

fprintf('\n%s: Fig 2\n', string(datetime('now')))

R_fig2 = R([R.iModelB_sim] == setting.fig12_Bsim & [R.iModelB_fit] == setting.fig12_Bfit);

h = figure('Position', [100 100 1800 700]);
tiledlayout(2, numel(setting.fig2_varFields), 'TileSpacing', 'compact', 'Padding', 'compact');

for iCol = 1:numel(setting.fig2_varFields)
    vField = setting.fig2_varFields{iCol};
    vName  = setting.fig2_varNames{iCol};
    vLevels = unique([R_fig2.(vField)]);

    for iRow=1:2 %1=ORI; 2=SF

        nexttile((iRow-1)*numel(setting.fig2_varFields) + iCol); hold on;
        [xORI, yTrueORI, curvesORI] = fxn_collapse_tuning_panel(R_fig2, vField, vLevels, namesFeature{iRow});
        if ~isempty(xORI)
            plot(axis_tuning{iRow}, yTrueORI, 'k-', 'LineWidth', setting.trueWidth);
            for iLev = 1:numel(curvesORI)
                if isempty(curvesORI(iLev).y), continue; end
                plot(axis_tuning{iRow}, curvesORI(iLev).y, 'LineWidth', setting.lineWidth, 'DisplayName', sprintf('%g', curvesORI(iLev).level));
            end
        end
        box on; set(gca, 'FontSize', setting.fontSize, 'LineWidth', setting.axisLineWidth);
        xticks(axisTicks_tuning{iRow})
        ylim([-.2, 1])
        title(vName); ylabel('Marg. weight (divided by max)');
        xlabel(namesFeature{iRow})
        l = legend('show', 'Location', 'best', 'FontSize', setting.fontSize-1); title(l, vField)
    end % iRow

end % iCol

sgtitle(sprintf('Fig 2: recovered tuning | filtered | Bsim = %d | Bfit = %d', ...
    setting.fig12_Bsim, setting.fig12_Bfit), 'FontWeight', 'bold');
saveas(h, fullfile(nameFolder_Figures_part4, 'Fig2_tuning.png'));
close(h);

%% Figure 3: metric recovery
%    Simulated curve from Bfit = 1 only; predictions from all Bfit

fprintf('\n%s: Fig 3\n', string(datetime('now')))

h = figure('Position', [100 100 1300 900]);
tiledlayout(numel(Bsim_unik), numel(setting.metricFields), 'TileSpacing', 'compact', 'Padding', 'compact');

for iCol_Bsim = 1:numel(Bsim_unik)
    simModel = Bsim_unik(iCol_Bsim);

    % simulated data are duplicated across Bfit; use Bfit = 1
    R_data = R([R.iModelB_sim] == simModel & [R.iModelB_fit] == 1);

    for iCol_metric = 1:numel(setting.metricFields)
        mName = setting.metricFields{iCol_metric};
        dataField = setting.metricCurveDataFields{iCol_metric};
        predField = setting.metricCurvePredFields{iCol_metric};

        nexttile; hold on;

        % collapsed simulated curve
        Cdata = fxn_collapse_curve(R_data, 'DV_allBins_med', dataField, setting.curveNGrid, setting.curveCollapseFcn);
        if ~isempty(Cdata.x)
            plot(Cdata.x, Cdata.y, 'ko', ...
                'MarkerSize', setting.markerSize, ...
                'MarkerFaceColor', 'k');
            fxn_plot_band_only(Cdata.x, Cdata.lb, Cdata.ub, [0 0 0], setting.bandAlpha);
        end

        % predicted curves across Bfit
        for iFit = 1:numel(Bfit_unik)
            fitModel = Bfit_unik(iFit);
            R_pred = R([R.iModelB_sim] == simModel & [R.iModelB_fit] == fitModel);

            Cpred = fxn_collapse_curve(R_pred, 'DV_allBins_med', predField, setting.curveNGrid, setting.curveCollapseFcn);
            if isempty(Cpred.x), continue; end

            fxn_plot_line_with_band(Cpred.x, Cpred.y, Cpred.lb, Cpred.ub, ...
                [0 0 0], setting.lineStyles_fit{iFit}, setting);
        end

        box on; set(gca, 'FontSize', setting.fontSize, 'LineWidth', setting.axisLineWidth);
        if iCol_Bsim == numel(Bsim_unik)
            xlabel('Collapsed DV bins');
        end
        if iCol_metric == 1
            ylabel(sprintf('Bsim = %d', simModel));
        end
        title(mName);

        switch mName
            case 'pYES'
                ylim([0 1]);
            case {'pC','pA'}
                ylim([0.5 1]);
        end
    end % iCol_metric
end % iCol_Bsim

sgtitle('Fig 3: metric recovery | filtered', 'FontWeight', 'bold');
fxn_addLegend_fitOnly(Bfit_unik, setting.lineStyles_fit);
saveas(h, fullfile(nameFolder_Figures_part4, 'Fig3_metricRecovery.png'));
close(h);

%% Figure 4: parameter recovery
%    Only matched pairs Bsim = Bfit

fprintf('\n%s: Fig 4\n', string(datetime('now')))

h = figure('Position', [100 100 1300 900]);
tiledlayout(numel(setting.fig4_paramNames), numel(Bfit_unik), 'TileSpacing', 'compact', 'Padding', 'compact');

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

        % one dot per signalCST × Cz after collapsing across IN parameters
        for iSig = 1:numel(gaborCST_unik)
            for iCz = 1:numel(Cz_unik)
                S = fxn_collapse_param_dot(Rsub, pName, gaborCST_unik(iSig), Cz_unik(iCz), setting.scalarCollapseFcn);

                if isempty(S.x) || isnan(S.y), continue; end

                fxn_plot_dot_with_err(S.x, S.y, S.lb, S.ub, ...
                    setting.cmap_Cz(iCz,:), setting.marker_signal{iSig}, setting);
            end
        end

        % unity line
        xl = xlim; yl = ylim;
        mn = min([xl(1), yl(1)]);
        mx = max([xl(2), yl(2)]);
        plot([mn mx], [mn mx], 'k--', 'LineWidth', 1);
        xlim([mn mx]); ylim([mn mx]);

        box on; axis square;
        set(gca, 'FontSize', setting.fontSize, 'LineWidth', setting.axisLineWidth);
        if iRow_param == numel(setting.fig4_paramNames)
            xlabel('True');
        end
        if iCol_Bsimfit == 1
            ylabel(pTitle);
        end
        title(sprintf('Bfit=Bsim=%d', Bsimfit));
    end % iCol_Bsimfit
end % iRow_param

sgtitle('Fig 4: parameter recovery | filtered | matched pairs only', 'FontWeight', 'bold');
fxn_addLegend_CzSignal_markers(setting.cmap_Cz, Cz_unik, setting.marker_signal, gaborCST_unik);
saveas(h, fullfile(nameFolder_Figures_part4, 'Fig4_paramRecovery.png'));
close(h);

fprintf('\n%s: simPlot4 done.\n', string(datetime('now')))

%% Local helpers
function S = fxn_collapse_scalar_by_x_simple(Rin, xField, yField, xLevels, collapseFcn)
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

    vals = [Rsub.(yField)]';
    vals = vals(isfinite(vals));
    if isempty(vals), continue; end

    x(end+1,1) = xLevels(i); %#ok<AGROW>
    y(end+1,1) = collapseFcn(vals, 'omitnan'); %#ok<AGROW>
    [~, lb_i, ub_i] = getCI(vals, 1, 1);
    lb(end+1,1) = lb_i; %#ok<AGROW>
    ub(end+1,1) = ub_i; %#ok<AGROW>
    n(end+1,1) = numel(vals); %#ok<AGROW>
end

S.x = x;
S.y = y;
S.lb = lb;
S.ub = ub;
S.n = n;
end

function S = fxn_collapse_scalar_by_x(Rin, xField, yField, fixedFields, fixedVals, xLevels, collapseFcn)
S = struct('x', [], 'y', [], 'lb', [], 'ub', [], 'n', []);
x = []; y = []; lb = []; ub = []; n = [];
for i = 1:numel(xLevels)
    idx = [Rin.(xField)] == xLevels(i);
    for j = 1:numel(fixedFields)
        idx = idx & ([Rin.(fixedFields{j})] == fixedVals(j));
    end
    Rsub = Rin(idx);
    if isempty(Rsub), continue; end
    vals = [Rsub.(yField)]';
    vals = vals(isfinite(vals));
    if isempty(vals), continue; end

    x(end+1,1) = xLevels(i); %#ok<AGROW>
    y(end+1,1) = collapseFcn(vals, 'omitnan'); %#ok<AGROW>
    [~, lb_i, ub_i] = getCI(vals, 1, 1);
    lb(end+1,1) = lb_i; %#ok<AGROW>
    ub(end+1,1) = ub_i; %#ok<AGROW>
    n(end+1,1) = numel(vals); %#ok<AGROW>
end
S.x = x; S.y = y; S.lb = lb; S.ub = ub; S.n = n;
end

function [xAxis, yTrue, curves] = fxn_collapse_tuning_panel(Rin, vField, vLevels, domainName)
xAxis = [];
yTrue = [];
curves = struct('level', {}, 'y', {}, 'lb', {}, 'ub', {}, 'n', {});

if isempty(Rin), return; end

switch domainName
    case 'ORI'
        fieldTrue = 'margORI_true';
        fieldEst  = 'margORI_est_med';
        nCh = numel(Rin(1).margORI_true);
    case 'SF'
        fieldTrue = 'margSF_true';
        fieldEst  = 'margSF_est_med';
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

    Y = vertcat(Rsub.(fieldEst));
    y = mean(Y, 1, 'omitnan');
    
    % max out
    y = y/max(y);

    lb = nan(1, nCh);
    ub = nan(1, nCh);
    for iCh = 1:nCh
        [~, lb(iCh), ub(iCh)] = getCI(Y(:, iCh), 1, 1);
    end

    curves(end+1).level = vLevels(i); %#ok<AGROW>
    curves(end).y = y;
    curves(end).lb = lb;
    curves(end).ub = ub;
    curves(end).n = size(Y,1);
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

function S = fxn_collapse_param_dot(Rin, pName, sigVal, czVal, collapseFcn)
S = struct('x', [], 'y', nan, 'lb', nan, 'ub', nan, 'n', 0);

idx = [Rin.gaborCST] == sigVal & [Rin.Cz_true] == czVal;
Rsub = Rin(idx);
if isempty(Rsub), return; end

switch pName
    case 'criterion_DV'
        xvals = [Rsub.criterion_DV_true]';
        yvals = [Rsub.criterion_DV_est_med]';
        lbvals = [Rsub.criterion_DV_est_lb]';
        ubvals = [Rsub.criterion_DV_est_ub]';
    otherwise
        xvals = [Rsub.(sprintf('%s_true', pName))]';
        yvals = [Rsub.(sprintf('%s_est_med', pName))]';
        lbvals = [Rsub.(sprintf('%s_est_lb', pName))]';
        ubvals = [Rsub.(sprintf('%s_est_ub', pName))]';
end

good = isfinite(xvals) & isfinite(yvals);
xvals = xvals(good);
yvals = yvals(good);
lbvals = lbvals(good);
ubvals = ubvals(good);

if isempty(xvals), return; end

S.x = collapseFcn(xvals, 'omitnan');
S.y = collapseFcn(yvals, 'omitnan');
S.lb = collapseFcn(lbvals, 'omitnan');
S.ub = collapseFcn(ubvals, 'omitnan');
S.n = numel(yvals);
end

function fxn_plot_line_with_band(x, y, lb, ub, colorRGB, lineStyle, P)
if isempty(x), return; end
fxn_plot_band_only(x, lb, ub, colorRGB, P.bandAlpha);
plot(x, y, 'Color', colorRGB, 'LineStyle', lineStyle, 'LineWidth', P.lineWidth);
end

function fxn_plot_band_only(x, lb, ub, colorRGB, alphaVal)
good = isfinite(x) & isfinite(lb) & isfinite(ub);
if ~any(good), return; end
x = x(good); lb = lb(good); ub = ub(good);
patch([x(:); flipud(x(:))], [lb(:); flipud(ub(:))], colorRGB, ...
    'FaceAlpha', alphaVal, 'EdgeColor', 'none');
end

function fxn_plot_dot_with_err(x, y, lb, ub, colorRGB, markerSym, P)
if isempty(x) || ~isfinite(y), return; end
line([x x], [lb ub], 'Color', colorRGB, 'LineWidth', P.errLineWidth);
plot(x, y, markerSym, ...
    'MarkerSize', P.markerSize, ...
    'MarkerFaceColor', colorRGB, ...
    'MarkerEdgeColor', colorRGB, ...
    'LineStyle', 'none');
end

function fxn_addLegend_CzSignal(cmap_Cz, Cz_unik, lineStyles_signal, gaborCST_unik)
ax1 = gca;

hCz = gobjects(numel(Cz_unik),1);
for i = 1:numel(Cz_unik)
    hCz(i) = plot(nan, nan, 'Color', cmap_Cz(i,:), 'LineWidth', 2);
end
leg1 = legend(ax1, hCz, ...
    arrayfun(@(x) sprintf('Cz = %g', x), Cz_unik, 'UniformOutput', false), ...
    'Location', 'northoutside', 'Orientation', 'horizontal');
title(leg1, 'Color');

ax2 = axes('Position', ax1.Position, 'Color', 'none', 'Visible', 'off');
hold(ax2, 'on');
hSig = gobjects(numel(gaborCST_unik),1);
for i = 1:numel(gaborCST_unik)
    hSig(i) = plot(ax2, nan, nan, 'k', 'LineStyle', lineStyles_signal{i}, 'LineWidth', 1.5);
end
leg2 = legend(ax2, hSig, ...
    arrayfun(@(x) sprintf('signal = %g', x), gaborCST_unik, 'UniformOutput', false), ...
    'Location', 'eastoutside');
title(leg2, 'Line style');
end

function fxn_addLegend_CzSignal_markers(cmap_Cz, Cz_unik, marker_signal, gaborCST_unik)
ax1 = gca;

hCz = gobjects(numel(Cz_unik),1);
for i = 1:numel(Cz_unik)
    hCz(i) = plot(nan, nan, 'Color', cmap_Cz(i,:), 'LineWidth', 2);
end
leg1 = legend(ax1, hCz, ...
    arrayfun(@(x) sprintf('Cz = %g', x), Cz_unik, 'UniformOutput', false), ...
    'Location', 'northoutside', 'Orientation', 'horizontal');
title(leg1, 'Color');

ax2 = axes('Position', ax1.Position, 'Color', 'none', 'Visible', 'off');
hold(ax2, 'on');
hSig = gobjects(numel(gaborCST_unik),1);
for i = 1:numel(gaborCST_unik)
    hSig(i) = plot(ax2, nan, nan, ...
        'LineStyle', 'none', ...
        'Marker', marker_signal{i}, ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k');
end
leg2 = legend(ax2, hSig, ...
    arrayfun(@(x) sprintf('signal = %g', x), gaborCST_unik, 'UniformOutput', false), ...
    'Location', 'eastoutside');
title(leg2, 'Marker');
end

function fxn_addLegend_fitOnly(Bfit_unik, lineStyles_fit)
ax = gca;
h = gobjects(numel(Bfit_unik),1);
for i = 1:numel(Bfit_unik)
    h(i) = plot(nan, nan, 'k', 'LineStyle', lineStyles_fit{i}, 'LineWidth', 1.5);
end
legend(ax, h, ...
    arrayfun(@(x) sprintf('Bfit = %d', x), Bfit_unik, 'UniformOutput', false), ...
    'Location', 'eastoutside');
end