%% [NOM] Metrics vs. DV; group averages x single locations
clc, fprintf('\n\n 14/24 Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'IV_allCond', 'nTrials_allCond', 'metric_*_allCond', 'R2_NOM_allCond')

% Define folder for saving figures for each model A and model B
nameFolder_Fig_NOM_metrics = sprintf('%s/NOMpred_A%d', nameFolder_Fig_NOM_Trialwise, iModelA_plot);
if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end

x_base   = 0.60;   % move a bit further left to make space for name column
x_width  = 0.4;   % total width of the table block
y_base   = 0.10;
y_height = 0.25;

sz_font = 15;
sz_axis = 18;
lineStyle_all = {'-', '-', '--', ':', ':', '-', '--', ':', '-.'};
namesPlotMode = {'Full', 'NoShared', 'NoMul', 'NoAdd', 'AllModels', 'NoPred'};

namesModelB_plot = {'Full model', 'No $\sigma_{shared}$', 'No $\sigma_{mul}$', 'No $\sigma_{add}$'};
szScaling = 10; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end
szBase = 5;

for iSet = 1:numel(iLocSingle_allSets)

    iLocSingle_perSet = iLocSingle_allSets{iSet};

    fprintf('\nL%s...', strjoin(string(iLocSingle_perSet), ''))

    for iPlotMode = 1:numel(namesPlotMode)

        % --------------------------
        % Decide what to plot
        % --------------------------
        switch iPlotMode
            case {1,2,3,4}
                % Plot ONE model per figure (B1..B4)
                iModelB_all_plot = iPlotMode;   % assumes your ModelB indices are 1..4
                flag_plotPred = true;
                flag_showR2table = true;

            case 5
                % Plot ALL models together
                iModelB_all_plot = iModelB_all; % e.g., [1 2 3 4]
                flag_plotPred = true;
                flag_showR2table = true;

            case 6
                % Plot DATA only (no predictions)
                iModelB_all_plot = 1;           % use a single source for data
                flag_plotPred = false;
                flag_showR2table = false;       % set true if you still want table (would be 1-row)
            otherwise
                error('Unknown iPlotMode=%d', iPlotMode);
        end

        nModelB_plot = numel(iModelB_all_plot);

        %%% --- PLOT ---
        for iMetric_prob = 1:nMetrics_prob
            figure('Position', [0 0 500 500]); hold on

            % R2 table: rows = plotted ModelB, cols = locations in this set
            R2_tab = nan(nModelB_plot, numel(iLocSingle_perSet));

            for rModel = 1:nModelB_plot
                iModelB = iModelB_all_plot(rModel);

                for iiLoc = 1:numel(iLocSingle_perSet)
                    locID = iLocSingle_perSet(iiLoc);

                    % Compute medians and CIs
                    [IV_ave, ~, ~, IV_SEM] = getCI(getCI(IV_allCond(iModelA_plot, iModelB, locID, :, :, :), 1, 5), 2, 1);
                    [nTrials_ave, ~, ~, nTrials_SEM] = getCI(getCI(nTrials_allCond(iModelA_plot, iModelB, locID, :, :, :), 1, 5), 2, 1);
                    [data_ave, ~, ~, data_SEM] = getCI(getCI(metric_data_allCond(iModelA_plot, iModelB, locID, iMetric_prob, :, :, :), 1, 6), 2, 1);

                    % pred exists but may not be plotted
                    [pred_ave, ~, ~, pred_SEM] = getCI(getCI(metric_pred_allCond(iModelA_plot, iModelB, locID, iMetric_prob, :, :, :), 1, 6), 2, 1);
                    [R2_NOM_ave, ~, ~, R2_NOM_SEM] = getCI(getCI(R2_NOM_allCond(iModelA_plot, iModelB, locID, iMetric_prob, :, :), 1, 6), 2, 1);

                    % store R2 mean
                    R2_tab(rModel, iiLoc) = R2_NOM_ave;

                    % ---------
                    % Styling
                    % ---------
                    cLoc = colors_comb(locID, :);

                    % For "all models together", differentiate models by linestyle/linewidth
                    % For "single model", keep clean.
                    if numel(iModelB_all_plot) > 1
                        ls = lineStyle_all{iModelB}; 
                        if iModelB == 1, lw = 3; else, lw = 1.5; end
                    else
                        ls = lineStyle_all{iModelB}; 
                        lw = 3;
                    end

                    % --------------------------
                    % Plot prediction (optional)
                    % --------------------------
                    if flag_plotPred
                        patch([IV_ave; flip(IV_ave)], ...
                            [pred_ave-pred_SEM; flip(pred_ave+pred_SEM)], ...
                            cLoc, 'FaceAlpha', .15, 'LineStyle', 'none', 'HandleVisibility', 'off');

                        plot(IV_ave, pred_ave, 'LineStyle', ls, 'Color', cLoc, 'LineWidth', lw, ...
                            'HandleVisibility', 'on');
                    end

                    % --------------------------
                    % Plot measurement (always)
                    % --------------------------
                    errorbar(IV_ave, data_ave, data_SEM, '.', 'vertical', 'CapSize', 0, ...
                        'Color', cLoc, 'HandleVisibility', 'off', 'LineWidth', max(lw/1.5, 1));

                    % Plot averaged data of each bin (dot size indicates number of trials)
                    for iBin = 1:nBins
                        plot(IV_ave(iBin), data_ave(iBin), 'o', ...
                            'MarkerEdgeColor', cLoc, 'MarkerFaceColor', 'w', ...
                            'MarkerSize', nTrials_ave(iBin) / szScaling + szBase, ...
                            'LineWidth', max(lw/1.5, 1), 'LineStyle', 'none', 'HandleVisibility', 'off');
                    end
                end % iiLoc
            end % rModel

            % reference line and y formatting
            yline(.5, '--', 'LineWidth', 2, 'Color', ones(1,3)/2);

            switch iMetric_prob
                case 1, ylim([0, 1]);    yticks(0:.2:1)
                case 2, ylim([.45, 1]);  yticks(.5:.1:1)
            end
            ylabel(namesMetrics_prob_full{iMetric_prob})

            if any(iLocSingle_perSet == 1), x_ticks = linspace(0, 180, 5);
            else,                          x_ticks = linspace(0, 80, 5);
            end
            xticks(x_ticks); xlim(x_ticks([1, end]))
            xlabel('Binned decision variable')

            % --------------------------
            % Legend: only when pred plotted
            % --------------------------
            if flag_plotPred
                if numel(iModelB_all_plot) == 1
                    legend(namesModelB{iModelB_all_plot(1)}, 'Location', 'best');
                else
                    legend(namesModelB(iModelB_all_plot), 'Location', 'best');
                end
            end

            % --------------------------
            % R2 table (kept)
            % --------------------------
            if flag_showR2table
                ax = gca;

                nRows = nModelB_plot;
                nCols = 1 + numel(iLocSingle_perSet); % [ModelName | one col per loc]

                x_pos   = linspace(x_base, x_base + x_width, nCols + 1);
                x_cells = x_pos(1:nCols) + diff(x_pos(1:2))/2;

                y_pos   = linspace(y_base, y_base + y_height, nRows + 1);
                y_cells = fliplr(y_pos(1:nRows)) + diff(y_pos(1:2))/2;

                % Column 1: model names
                for r = 1:nRows
                    iModelB = iModelB_all_plot(r);
                    name_str = sprintf('%s', namesModelB_plot{iModelB}); % your existing namesModelB_plot
                    text(x_cells(1), y_cells(r), name_str, ...
                        'Units','normalized', 'HorizontalAlignment','left', ...
                        'VerticalAlignment','middle', 'FontSize', sz_font, ...
                        'FontWeight','bold', 'Interpreter','latex');
                end

                % Columns 2..end: R2 values
                for r = 1:nRows
                    for c = 1:numel(iLocSingle_perSet)
                        R2_val = R2_tab(r, c);
                        if isnan(R2_val), continue; end
                        str_cell = sprintf('%.0f%%', R2_val * 100);
                        text(x_cells(1+c), y_cells(r), str_cell, ...
                            'Units','normalized', 'HorizontalAlignment','center', ...
                            'VerticalAlignment','middle', 'FontSize', sz_font);
                    end
                end
            end

            % axis cosmetics
            ax = gca;
            ax.FontSize = sz_axis;
            ax.LineWidth = 1.5;

            title(sprintf('n=%d [A%d] [L%s] [nIter=%d] %s | %s', ...
                nsubj, iModelA_plot, strjoin(string(iLocSingle_perSet), ''), nIterxJob, ...
                namesMetrics_prob{iMetric_prob}, namesPlotMode{iPlotMode}), ...
                'FontSize', 10);

            % Save
            saveas(gcf, sprintf('%s/n%d_L%s_%s_group_%s.png', ...
                nameFolder_Fig_NOM_metrics, nsubj, strjoin(string(iLocSingle_perSet), ''), ...
                namesMetrics_prob{iMetric_prob}, namesPlotMode{iPlotMode}));
            close(gcf)

        end % iMetric_prob
    end % iPlotMode
end % iSet