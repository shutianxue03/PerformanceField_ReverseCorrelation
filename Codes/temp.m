%% Behavioral metrics: single locations
clc, fprintf('\n\n Plotting STARTED......\n\n')
% Load data
load(nameFolder_Data_SaveCompile, 'metrics_allCond')

% Define folder for saving figures
nameFolder_Fig_behav = sprintf('%s/Behav', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_behav)), mkdir(nameFolder_Fig_behav), end

flag_plotIDVD = 1;
flag_plotDiff = 1;
namesMetrics_behav = {'CS', 'pA', 'dprime', 'criterion', 'pC', 'RT'}; nMetrics_behav = length(namesMetrics_behav);
[d70,~] = SX_sim06_SDT(.7, .3);

sz_wd_perBar = 70;

for iSet = 1:numel(iLocSingle_allSets)

    iLocSingle_all = iLocSingle_allSets{iSet};

    nBars = numel(iLocSingle_all);
    sz_fig = [nBars*sz_wd_perBar, 350];

    for iMetric_prob=1:nMetrics_behav
        switch iMetric_prob
            case 1, iMetric_vec = 10; x_ticks = linspace(1.6, 3.6, 5); flag_plotIDVD = 1; ref=nan;
            case 2, iMetric_vec = 6; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 1; ref=nan;
            case 3, iMetric_vec = 1; x_ticks = linspace(0, 1.6, 5); flag_plotIDVD = 0; ref=d70;
            case 4, iMetric_vec = 2; x_ticks = linspace(-1,1, 5); flag_plotIDVD = 0; ref=0;
            case 5, iMetric_vec = 3; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 0; ref=nan;
            case 6, iMetric_vec = 11; x_ticks = linspace(0, .2, 5); flag_plotIDVD = 0; ref=nan;
        end
        x_ticks = round(x_ticks, 2);

        data_allIter_allSubj = squeeze(metrics_allCond(iModelA_plot, iLocSingle_all, :, :, iDataset_plotNOM, iMetric_vec));

        str_title = sprintf('%s nIter=%d L%s', namesMetrics_behav{iMetric_prob}, nIterxJob, strjoin(string(iLocSingle_all), ''));
        %------------------------------%
        basicFxn_drawBars_permutation(data_allIter_allSubj, ref, colors_comb(iLocSingle_all, :), namesLocComb(iLocSingle_all), x_ticks, x_ticks, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj)
        % ------------------------------%
        % ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily, 2}{iTunC}))
        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%s_%s.png', nameFolder_Fig_behav, nsubj, strjoin(string(iLocSingle_all), ''), namesMetrics_behav{iMetric_prob}))
        close(gcf)
    end % iMetric
end % iSet
clear metrics_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Behavioral metrics: paired locations
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'metrics_allCond')

sz_wd_perBar = 200;
nBars = 2;
sz_fig = [nBars*sz_wd_perBar, 250];

for iGroup = 1:nGroups
    iLocPair_all = iLocGroups_all{iGroup};
    for iMetric_prob = 1:nMetrics_behav
        switch iMetric_prob
            case 1, iMetric_vec = 10; x_ticks = linspace(1.6, 3.6, 5); flag_plotIDVD = 1; ref=nan;
            case 2, iMetric_vec = 6; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 1; ref=nan;
            case 3, iMetric_vec = 1; x_ticks = linspace(0, 1.6, 5); flag_plotIDVD = 0; ref=d70;
            case 4, iMetric_vec = 2; x_ticks = linspace(-1,1, 5); flag_plotIDVD = 0; ref=0;
            case 5, iMetric_vec = 3; x_ticks = linspace(.5, .9, 5); flag_plotIDVD = 0; ref=nan;
            case 6, iMetric_vec = 11; x_ticks = linspace(0, .2, 5); flag_plotIDVD = 0; ref=nan;
        end
        x_ticks = round(x_ticks, 2);

        data_allIter_allSubj = squeeze(metrics_allCond(iModelA_plot, iLocPair_all, :, :, iDataset_plotNOM, iMetric_vec));

        str_title = sprintf('%s nIter=%d L%d%d', namesMetrics_behav{iMetric_prob}, nIterxJob, iLocPair_all);
        %------------------------------%
        basicFxn_drawBars_permutation(data_allIter_allSubj, ref, colors_comb(iLocPair_all, :), namesLocComb(iLocPair_all), x_ticks, x_ticks, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj)
        %------------------------------%
        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_%s.png', nameFolder_Fig_behav, nsubj, iLocPair_all, namesMetrics_behav{iMetric_prob}))
        close(gcf)
    end % iMetric
end % iGroup
clear metrics_allCond
fprintf('\n\n Plotting DONE\n\n')


%% Plot templates for IDVD and group averages
clc, fprintf('\n\n Plotting STARTED......\n\n')
% Load data
load(nameFolder_Data_SaveCompile, 'template_*_allCond')

pixMin = 0;
pixMax = 0;
% Define folder for saving figures
nameFolder_Fig_NOM_Template = sprintf('%s/Template', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_NOM_Template)), mkdir(nameFolder_Fig_NOM_Template), end

for iDataset = 1:nDatasets
    switch iDataset % needs to matchOOD_xx_compIV ("for iDataset = 1:2")
        case 1
            template_allSubj = squeeze(template_tmpl_allCond(iModelA_plot, iModelB_plot, :, :, :, :)); % ... nLoc x nSubj x nORI x nSF
        case 2
            template_allSubj = squeeze(template_full_allCond(iModelA_plot, iModelB_plot, :, :, :, :)); % ... nLoc x nSubj x nORI x nSF
    end
    str_dataset = namesDataset{iDataset};

    %%%% Group ave %%%%%
    for iLocComb = iLocComb_all
        figure('Position', [0 0 1e3 1e3]), hold on
        e2D_ave = getCI(template_allSubj(iLocComb, :, :, :), 2, 2);
        e2D_ave = e2D_ave';
        fprintf('\n%s L%d: Min = %.3f, Max=%.3f\n', str_dataset, iLocComb, min(e2D_ave(:)), max(e2D_ave(:)))

        if iLocComb==1
            cLim = [-.03, .26];
        else
            cLim = [-.02, .14];
            pixMin = min(pixMin, min(e2D_ave(:)));
            pixMax = max(pixMax, max(e2D_ave(:)));
        end
        RCplot_2Dkernel(e2D_ave, cLim)
        title(sprintf('n=%d nIter=%d L%d %s (%s)', nsubj, nIterxJob, iLocComb, namesLocComb{iLocComb}, str_dataset))
        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_group_L%d_A%dB%d_%s.png', nameFolder_Fig_NOM_Template, nsubj, iLocComb, iModelA_plot, iModelB_plot, str_dataset))
        close(gcf)
    end % iiLoc

    %%%%% Idvd data in one figure, per loc %%%%%
    % for iiLoc = 1:nLocComb8
    %
    %     figure('Position', [0, 0, 2e3, 1.8e3])
    %     for isubj = 1:nsubj
    %
    %         e2D = squeeze(template_allSubj(iLocComb_all(iiLoc), isubj, :,:))';
    %
    %         subplot(nRows, nCols, isubj), hold on
    %         RCplot_2Dkernel(e2D)
    %         title(subjList{isubj})
    %
    %     end % isubj
    %     set(findall(gcf, '-property', 'FontSize'), 'FontSize',12)
    %     sgtitle(sprintf('L%d %s [A%dB%d] (%s)', iLocComb_all(iiLoc), namesLocComb{iLocComb_all(iiLoc)}, iModelA_plot, iModelB_plot, str_title), 'FontSize',20)
    %
    %     % save
    %     drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d_A%dB%d_%s.png', nameFolder_Fig_NOM_Template, nsubj, iLocComb_all(iiLoc), iModelA_plot, iModelB_plot, str_title))
    %     close(gcf)
    %
    % end % iiLoc
end % iDataset

% Print min and max pixel value
fprintf('\n Summary: Min = %.3f, Max=%.3f\n', pixMin, pixMax)

clear template_*_allCond
fprintf('\n\n Plotting DONE\n\n')


%% Plot separability (use basicFxn_drawBars_permutation only)
clc; fprintf('\n\n Plotting STARTED......\n\n')
str_dataset = namesDataset{iDataset_plotRC};

% Load compiled separability: sep_allCond
load(nameFolder_Data_SaveCompile, 'sep_allCond')

% Folder for saving figures
nameFolder_Fig_Sep = sprintf('%s/Separability', nameFolder_Fig_NOM_Trialwise);
if ~exist(nameFolder_Fig_Sep, 'dir'), mkdir(nameFolder_Fig_Sep); end

% ---- Loop across sets of locations ----
for iSet = 1:numel(iLocSingle_allSets)

    % Locations included in this figure
    iLocSingle_all = iLocSingle_allSets{iSet};
    nLoc = numel(iLocSingle_all);

    % Extract separability for this set
    % sep_allIter_allSubj: [nLoc x nSubj x nIter]
    sep_allIter_allSubj = squeeze(sep_allCond(iModelA_plot, iModelB_plot, iLocSingle_all, :, :, iDataset_plotRC));

    % Reformat to match basicFxn_drawBars_permutation:
    % data_allIter_allSubj: [nIter x nSubj x nCond], where nCond = nLoc
    data_allIter_allSubj = permute(sep_allIter_allSubj, [3 2 1]);  % [nIter x nSubj x nLoc]

    % Colors for each bar/condition (each location)
    colors = colors_comb(iLocSingle_all, :); % [nLoc x 3]

    % Condition labels (not displayed on x-axis inside the function, but used in title strings)
    x_ticks = cell(1, nLoc);
    for iiLoc = 1:nLoc
        x_ticks{iiLoc} = sprintf('L%d', iLocSingle_all(iiLoc));
    end

    % Plot settings for this figure
    ref = 0.5;                 % reference line used in your old plot
    y_ticks = [0.5 0.6 0.7 0.8 0.9 1.0];
    y_ticklabels = y_ticks;

    flag_plotIDVD = 1;         % show subject-level lines (recommended)
    flag_plotDiff = 0;         % only meaningful when nCond==2 (turn on if you want)
    sz_fig = [800 600];

    str_loc = strjoin(string(iLocSingle_all), '');
    str_title = sprintf('n%d nIter=%d [A%dB%d] L%s %s', ...
        nsubj, nIterxJob, iModelA_plot, iModelB_plot, str_loc, str_dataset);

    % ---- draw ----
    basicFxn_drawBars_permutation( ...
        data_allIter_allSubj, ...
        ref, colors, x_ticks, ...
        y_ticks, y_ticklabels, ...
        flag_plotIDVD, flag_plotDiff, ...
        str_title, sz_fig, ...
        nIterxJob, markers_allSubj);

    % ---- save ----
    outName = sprintf('%s/n%d_L%s_A%dB%d_%s.png', ...
        nameFolder_Fig_Sep, nsubj, str_loc, iModelA_plot, iModelB_plot, str_dataset);
    saveas(gcf, outName);
    close(gcf);

end % iSet

clear sep_allCond
fprintf('\n\n Plotting DONE\n\n')


%% Plot separability
% clc; fprintf('\n\n Plotting STARTED......\n\n')
% str_dataset = namesDataset{iDataset_plotRC};
% 
% % Load compiled separability: sep_allCond
% % Expected dims (based on your squeeze usage below):
% %   sep_allCond: [nModelsA x nModelsB x nLoc x nSubj x nIter x nSet]
% load(nameFolder_Data_SaveCompile, 'sep_allCond')
% 
% % Folder for saving figures
% nameFolder_Fig_Sep = sprintf('%s/Separability', nameFolder_Fig_NOM_Trialwise);
% if ~exist(nameFolder_Fig_Sep, 'dir'), mkdir(nameFolder_Fig_Sep); end
% 
% wd_border = 2;
% sz_title  = 10;
% sz_label  = 25;     % 40
% sz_ticks  = 35;
% wd_all    = 3;
% sz_marker = 20;
% step      = 0.1;    % spacing between locations (x-jitter across locs)
% 
% % ---- Loop across sets of locations ----
% for iSet = 1:numel(iLocSingle_allSets)
% 
%     % location indices included in this figure (e.g., [1 3 5] etc.)
%     iLocSingle_all = iLocSingle_allSets{iSet};
%     nLoc = numel(iLocSingle_all);
% 
%     % Extract separability for this set
%     % After squeeze: sep_allIter_allSubj is [nLoc x nSubj x nIter]
%     sep_allIter_allSubj = squeeze(sep_allCond(iModelA_plot, iModelB_plot, iLocSingle_all, :, :, iDataset_plotRC));
% 
%     % ==========================
%     %  STATS (ANOVA + effect size CI + permutation p)
%     % ==========================
% 
%     % indLoc labels for one-way ANOVA on vectorized data:
%     % indLoc is [nLoc x nSubj], with values 1..nLoc repeated across subjects
%     indLoc = repmat(1:nLoc, nsubj, 1)';   % [nLoc x nSubj]
% 
%     % --- Per-iteration F and eta^2 (partial) for CI across iters ---
%     eta2_allIter = nan(nIterxJob, 1);
%     F_allIter    = nan(nIterxJob, 1);
% 
%     for iIter = 1:nIterxJob
%         % sep_perIter: [nLoc x nSubj] at this iteration
%         sep_perIter = squeeze(sep_allIter_allSubj(:, :, iIter));
% 
%         % One-way ANOVA with factor Loc (your helper prints + returns table)
%         % tbl indexing: you already use tbl{2,7} elsewhere for p;
%         % Here we assume (row2=Loc) and columns:
%         %   col3=df, col5=F, row3 col3=df_error
%         [~, tbl] = print_nANOVA({'Loc'}, sep_perIter(:), {indLoc(:)}, nsubj);
% 
%         F   = tbl{2,5};
%         df1 = tbl{2,3};
%         df2 = tbl{3,3};
% 
%         F_allIter(iIter)    = F;
%         eta2_allIter(iIter) = (F*df1) / (F*df1 + df2);  % partial eta^2
%     end % iIter
% 
%     % --- CI across iters (median + CI) ---
%     [F_med,   F_lb,   F_ub]   = getCI(F_allIter,    1, 1);
%     [eta_med, eta_lb, eta_ub] = getCI(eta2_allIter, 1, 1);
% 
%     % --- Median separability map across iters + CI for plotting ---
%     % sep_med: [nLoc x nSubj]
%     % sep_neg/sep_pos: lower/upper error bars per point (same dims)
%     [sep_med, ~, ~, sep_neg, sep_pos] = getCI(sep_allIter_allSubj, 1, 3);
% 
%     % --- Permutation test: RM-null by shuffling condition labels WITHIN subject ---
%     % This corresponds to: under H0, location labels are exchangeable within each subject.
%     nPerm = 1e4;
%     Fvalue_allPerm = nan(nPerm, 1);
% 
%     % Observed F-value computed on the *median* separability map
%     Fvalue_obs = rm_oneway(sep_med);   % expects [nLoc x nSubj]
% 
%     for iPerm = 1:nPerm
%         data_shuffled = sep_med;
% 
%         % Shuffle location labels within each subject (column-wise permutation)
%         % This preserves each subject’s marginal distribution.
%         for iSubj = 1:nsubj
%             data_shuffled(:, iSubj) = data_shuffled(randperm(nLoc), iSubj);
%         end
% 
%         % Compute permuted F
%         Fvalue_allPerm(iPerm) = rm_oneway(data_shuffled);
%     end % iPerm
% 
%     % Right-tailed p-value (F >= 0)
%     p_ANOVA_perm = (1 + sum(Fvalue_allPerm >= Fvalue_obs)) / (nPerm + 1);
% 
%     % Summary string (kept your style)
%     str_ANOVA = sprintf(['ANOVA: F=%.3f [%.3f, %.3f], perm p=%.3f, ' ...
%                          'eta2=%.3f [%.3f, %.3f]\n'], ...
%                          F_med, F_lb, F_ub, p_ANOVA_perm, eta_med, eta_lb, eta_ub);
% 
%     % ==========================
%     %  PLOT
%     % ==========================
%     figure('Position', [0 0 1e3 600]); hold on; box on
%     yline(.5, '--', 'LineWidth', wd_all);
% 
%     str_sepValues = "";
% 
%     % Loop across locations (jittered by location)
%     for iiLoc = 1:nLoc
%         locID = iLocSingle_all(iiLoc);
% 
%         % Centered jitter across locations around each subject index
%         offset = (iiLoc - (nLoc + 1)/2) * step;
%         x = (1:nsubj) + offset;
% 
%         % Plot per-subject point + CI at this location
%         for isubj = 1:nsubj
%             c = colors_comb(locID, :);
% 
%             % CI bars per subject/location (asymmetric error: sep_neg/sep_pos)
%             errorbar(x(isubj), sep_med(iiLoc, isubj), ...
%                 sep_neg(iiLoc, isubj), sep_pos(iiLoc, isubj), ...
%                 'Color', c, 'CapSize', 0);
% 
%             % point
%             plot(x(isubj), sep_med(iiLoc, isubj), markers_allSubj{isubj}, ...
%                 'MarkerFaceColor', c, 'MarkerEdgeColor', c, ...
%                 'MarkerSize', sz_marker, 'LineWidth', 3);
%         end
% 
%         % Append brief summary per location: mean (SEM)
%         mu  = mean(sep_med(iiLoc, :), 'omitnan');
%         sem = std(sep_med(iiLoc, :), 'omitnan') / sqrt(nsubj);
%         str_sepValues = str_sepValues + sprintf('L%d: %.2f (%.2f); ', locID, mu, sem);
%     end
% 
%     % Axes labels/ticks
%     ylabel('Correlation', 'FontSize', sz_label);
%     yticks([0, .5:.1:1]);
%     ylim([.5, 1]);
% 
%     xticks(1:nsubj);
%     xticklabels(1:nsubj);       % keep numeric IDs
%     % xticklabels(subjList);    % (optional) toggle later
%     xlim([0, nsubj + 1]);
%     xlabel('Observer #');
% 
%     % Style
%     ax = gca;
%     ax.XAxis.FontSize = sz_ticks;
%     ax.YAxis.FontSize = sz_ticks;
%     ax.LineWidth      = wd_border;
% 
%     set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 2);
%     set(findall(gcf, '-property', 'FontSize'),  'FontSize', 30);
% 
%     % Title + save
%     str_loc = strjoin(string(iLocSingle_all), '');
%     title(sprintf('n%d nIter=%d [A%dB%d] L%s %s\n%s\n%s', nsubj, nIterxJob, iModelA_plot, iModelB_plot, str_loc, str_dataset, str_sepValues, str_ANOVA), 'FontSize', sz_title);
% 
%     drawnow; pause(1);
% 
%     outName = sprintf('%s/n%d_L%s_A%dB%d_%s.png', ...
%         nameFolder_Fig_Sep, nsubj, str_loc, iModelA_plot, iModelB_plot, str_dataset);
%     saveas(gcf, outName);
%     close(gcf);
% 
% end % iSet
% 
% clear sep_allCond
% fprintf('\n\n Plotting DONE\n\n')


%% Tuning functions for group averages
clc, fprintf('\n\n Plotting STARTED......\n\n')
% Load data
load(nameFolder_Data_SaveCompile, 'marg*_allCond', 'margPred*_allCond', 'margR2*_allCond')

% Define folder for saving figures
nameFolder_Fig_NOM_Tuning = sprintf('%s/Tuning', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_NOM_Tuning)), mkdir(nameFolder_Fig_NOM_Tuning), end

wd_border = 4; % default 5
sz_ticks = 35;% default 35

iplots = reshape((1:nGroups*nFeatures)', nGroups, nFeatures)';
for iDataset=1:nDatasets
    switch iDataset
        case 1
            markerStyle = 'o';
            lineStyle = '-';
        case 2
            markerStyle = 's';
            lineStyle = '--';
    end

    for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
        iLocPair_all = iLocGroups_all{iGroup};

        for iFeature = 1:nFeatures

            xaxis = axis_tuning{iFeature};

            figure('Position', [0 0 1.1e3 8e2]) % default 8e2
            hold on
            for iLoc = iLocPair_all

                % Set color
                color_comb = colors_comb(iLoc, :);

                % Set y-ticks
                if iFeature==1 % ORI tuning fxn
                    marg_allCond = margORI_allCond;
                    margPred_allCond = margPred_ORI_allCond;
                    margR2_allCond = margR2_ORI_allCond;
                    if find(iLocPair_all==1), yticks_ = [-.03, linspace(0, .24, 4)]; % fov vs. peri, higher ub
                    else, yticks_ = [-.03, linspace(0, .12, 4)];
                    end
                else % SF tuning fxn
                    marg_allCond = margSF_allCond;
                    margPred_allCond = margPred_SF_allCond;
                    margR2_allCond = margR2_SF_allCond;
                    if find(iLocPair_all==1), yticks_ = [-.01, linspace(0, .12, 4)];
                    else, yticks_ = [-.01, linspace(0, .06, 4)]; %[-.02, 0, .02, .04, .06];
                    end
                end

                ymax = max(yticks_);
                ymin = min(yticks_);

                % Obtain data and prediction
                [marg_ave, ~, ~, marg_sem] = getCI(getCI(marg_allCond(iModelA_plot, iModelB_plot, iLoc, :, :, iDataset, :), 1, 5), 2, 1);
                [margPred_ave, ~, ~, margPred_sem] = getCI(getCI(margPred_allCond(iModelA_plot, iModelB_plot, iLoc, :, :, iDataset, :), 1, 5), 2, 1);
                [R2_ave, ~, ~, R2_sem] = getCI(getCI(margR2_allCond(iModelA_plot, iModelB_plot, iLoc, :, :, iDataset), 1, 5), 2, 1);

                % Data (dots + errorbars)
                errorbar(xaxis, marg_ave, marg_sem, markerStyle, 'Color', color_comb, 'CapSize',0)
                plot(xaxis, marg_ave, 'o', 'color', color_comb, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'HandleVisibility','off')

                % Prediction (lines + bands)
                patch([xaxis, flip(xaxis)], [margPred_ave-margPred_sem, flip(margPred_ave+margPred_sem)], color_comb, 'FaceAlpha', .3, 'linestyle', 'none')
                plot(xaxis, margPred_ave, '-', 'color', color_comb)

                % Draw reference lines
                yline(0, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);
                xline(iFeature-1, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);

                % Set ticks, labels and limits (for EACH loc, to print R2 at the right loc)
                xlabel(namesFeature_axis{iFeature})
                ylabel(namesFeature_axis_Tuning{iFeature})

                xlim(axisLim{iFeature})
                ylim([ymin, ymax])
                yticks(yticks_)
                %     xlabel(xlabels_tuning{ifeature}, 'FontSize', sz_label)
                xticks(axisTicks_tuning{iFeature})
                xticklabels(axisTL_tuning{iFeature})
                if iFeature==2, xticklabels(round(axisTL_tuning{iFeature}, 2)), end

                % Print R2 in the figure
                xlims = axisLim{iFeature};
                switch iFeature
                    case 1, x_R2 = xlims(1) + 0.12 * range(xlims);   % slightly right of left boundary
                    case 2, x_R2 = xlims(2) - 0.3 * range(xlims);   % slightly left of rightboundary
                end
                y_R2 = ymax - 0.01-0.08*(find(iLoc == iLocPair_all)-1) * (ymax-ymin);   % slightly below top boundary

                text(x_R2, y_R2, ...
                    sprintf('$R^2 = %.2f \\pm %.2f$', R2_ave, R2_sem), ...
                    'Interpreter', 'latex', ...
                    'FontSize', 30, ...
                    'Color', color_comb, ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'top');
            end % iLoc

            ax = gca;
            ax.XAxis.FontSize = sz_ticks;
            ax.YAxis.FontSize = sz_ticks;
            ax.LineWidth = wd_border/1.5;

            [R2_ave, ~, ~, R2_sem] = getCI(getCI(margR2_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset), 1, 5), 2, 2);

            title(sprintf('n=%d L%d%d [A%dB%d] %s tuning (%s)', nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, namesDataset{iDataset}))

            set(findall(gcf, '-property', 'linewidth'), 'linewidth', 2)

            % Save the figure
            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_group_L%d%d_A%dB%d_%s_%s.png', nameFolder_Fig_NOM_Tuning, nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, namesDataset{iDataset}))
            close(gcf)
        end % iFeature
    end % iGroup
end % iDataset = 1:2
clear marg*_allCond margPred*_allCond margR2*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Tuning functions for idvd (LARGE FIGURES!!)
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'marg*_allCond', 'margPred*_allCond', 'margParams*_allCond')

% Define folder for saving figures
nameFolder_Fig_NOM_Tuning = sprintf('%s/TuningFxns', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_NOM_Tuning)), mkdir(nameFolder_Fig_NOM_Tuning), end

wd_border = 4; % default 5
sz_ticks = 30;% default 35

for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
    iLocPair_all = iLocGroups_all{iGroup};

    for iFeature = 2%1:nFeatures

        xaxis = axis_tuning{iFeature};

        for iDataset = 1:nDatasets
            switch iDataset
                case 1 % TmplSet
                    markerStyle = 'o';
                    lineStyle = '-';
                case 2 % FullSet
                    markerStyle = 's';
                    lineStyle = '--';
            end

            figure('Position', [0 0 2e3 1.5e3])

            for isubj = 1:nsubj
                subplot(nRows_subj, nCols_subj, isubj), hold on

                subjName = subjList{isubj};

                str_TunParams ='';

                for iLoc = iLocPair_all
                    str_TunParams = sprintf('%s\nL%d', str_TunParams, iLoc);

                    % fprintf('\nL%d...', iLoc)
                    % Set color
                    color_comb = colors_comb(iLoc, :);

                    % Set y-ticks
                    if iFeature==1
                        marg_allCond = margORI_allCond;
                        margPred_allCond = margPred_ORI_allCond;
                        margParams_allCond = margParams_ORI_allCond;
                        if iLoc==1, yticks_ = [-.03, 0, .05, .10, .15]; % fov vs. peri, higher ub
                        else, yticks_ = [-.02, linspace(0, .12, 4)];
                        end
                    else
                        marg_allCond = margSF_allCond;
                        margPred_allCond = margPred_SF_allCond;
                        margParams_allCond = margParams_SF_allCond;
                        if iLoc==1, yticks_ = [-.01, linspace(0, .08, 4)];
                        else, yticks_ = [-.01, linspace(0, .06, 4)]; %[-.02, 0, .02, .04, .06];
                        end
                    end

                    ymax = max(yticks_);
                    ymin = min(yticks_);

                    % Obtain data, prediction, and params
                    [marg_med, ~, ~, marg_lb, marg_ub] = getCI(marg_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :), 1, 5);
                    [margPred_med, margPred_lb, margPred_ub] = getCI(margPred_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :), 1, 5);
                    [margParam_med, margParam_lb, margParam_ub] = getCI(margParams_allCond(iModelA_plot, iModelB_plot, iLoc, isubj, :, iDataset, :), 1, 5);

                    % Data (dots + errorbars)
                    errorbar(xaxis, marg_med, marg_lb, marg_ub, markerStyle, 'Color', color_comb, 'CapSize',0)
                    plot(xaxis, marg_med, markerStyle, 'color', color_comb, 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'HandleVisibility','off')

                    % Prediction (lines + bands)
                    patch([xaxis, flip(xaxis)], [margPred_lb', flip(margPred_ub')], color_comb, 'FaceAlpha', .3, 'linestyle', 'none', 'linewidth', 2)
                    plot(xaxis, margPred_med, '-', 'color', color_comb)

                    % Print estimated parameters
                    % each loc has 3-4 params, think of how to arrange (maybe drop CI)
                    for iTunParam = 1:length(margParam_med)
                        str_TunParams = sprintf('%s | %.2f', str_TunParams, margParam_med(iTunParam));
                    end % iTunParam
                    % str_TunParams = sprintf('%s\n', str_TunParams);
                end % iLoc

                % Draw reference lines
                yline(0, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);
                xline(iFeature-1, 'handlevisibility', 'off', 'linewidth', wd_border, 'color', [.7, .7, .7]);

                xlabel('Orientation (deg)')
                ylabel('Amplitude')
                title(sprintf('%s\n%s\n', subjName, str_TunParams))
                % legend('Location', 'best')

            end % isubj
            sgtitle(sprintf('n=%d L%d%d [A%dB%d] %s tuning (%s)', nsubj, iLocPair_all,iModelA_plot, iModelB_plot, namesFeature{iFeature}, namesDataset{iDataset}))
            set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10)

            % Save the figure
            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_A%dB%d_%s_%s.png', nameFolder_Fig_NOM_Tuning, nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, namesDataset{iDataset}))
            close(gcf)
        end % i=1:2
    end % iFeature
end % iGroup
clear marg*_allCond margPred*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Tuning characteristics
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'margTunC_*_allCond')

% Define folder for saving figures
nameFolder_Fig_tunC = sprintf('%s/TuningCs', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_tunC)), mkdir(nameFolder_Fig_tunC), end

flag_plotDist = 0;
flag_plotIDVD = 1;
flag_plotDiff = 1;
paramMode = 2;
sz_fig = [350 350]; % size of the figure canvas

for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
    iLocPair_all = iLocGroups_all{iGroup};

    % Plotting settings
    colors = colors_comb(iLocPair_all, :);
    x_ticks = namesLocComb(iLocPair_all);

    % Loop through each feature
    for iFeature = 1:nFeatures

        % Define the x-axis and parameters for the current feature
        xaxis = axis_tuning{iFeature};
        nfilters = length(xaxis);
        iFamily = iFamily_perF(iFeature);
        namesTunC = namesTunC_unit_perF{iFamily, paramMode};
        nTunC_full = length(namesTunC);

        y_ticks_all = [];

        % Obtain y_ticks of ORI/SF domain
        switch flag_plotDist
            case 0
                switch iFamily_perF(iFeature)
                    case 1 % ORI tuniningC | Scaled Gaussian
                        if flag_plotIDVD
                            if find(iLocPair_all==1)
                                y_ticks_all{1} = linspace(0, .4, 5); % ORI peak amp
                            else
                                y_ticks_all{1} = linspace(0, .2, 5); % ORI peak amp
                            end
                            y_ticks_all{2} = linspace(0, 100, 5); % ORI band
                            y_ticks_all{3} = linspace(-.05, .03, 5); % ORI baseline
                        else
                            if find(iLocPair_all==1)
                                y_ticks_all{1} = linspace(0, .24, 5); % ORI peak amp
                            else
                                y_ticks_all{1} = linspace(0, .12, 5); % ORI peak amp
                            end

                            y_ticks_all{2} = linspace(0, 100, 5); % ORI band
                            y_ticks_all{3} = linspace(-.08, .08, 5); % ORI baseline
                        end

                    case 2 % SF tuniningC | log parabola
                        if flag_plotIDVD

                            y_ticks_all{1} = linspace(0, 2, 5); % peak SF
                            if find(iLocPair_all==1)
                                y_ticks_all{2} = linspace(.01, .2, 5); % SF peak amp
                            else
                                y_ticks_all{2} = linspace(.01, .09, 5); % SF peak amp
                            end
                            y_ticks_all{3} = linspace(0, 2, 5); % SF bandwidth
                            y_ticks_all{4} = linspace(-.1, .1, 5); % SF baseline
                        else
                            y_ticks_all{1} = linspace(0, 2, 5); % peak SF

                            if find(iLocPair_all==1)
                                y_ticks_all{2} = linspace(.02, .12, 5); % SF peak amp
                            else
                                y_ticks_all{2} = linspace(.02, .06, 5); % SF peak amp
                            end
                            y_ticks_all{3} = linspace(.3, 2, 5); % SF bandwidth
                            y_ticks_all{4} = linspace(-.08, 0, 5); % SF baseline
                        end
                end
            case 1 % plot distribution of group averages, not bars
                switch iFamily_perF(iFeature)
                    case 1 % ORI tuniningC | Scaled Gaussian
                        if find(iLocPair_all==1)
                            y_ticks_all{1} = linspace(0, .3, 5); % ORI peak amp (fovea)
                        else
                            y_ticks_all{1} = linspace(.06, .14, 5); % ORI peak amp
                        end
                        y_ticks_all{2} = linspace(10, 50, 5); % ORI band
                        y_ticks_all{3} = linspace(-.03, .01, 5); % ORI baseline

                    case 2 % SF tuniningC | log parabola
                        y_ticks_all{1} = linspace(log2(1.4), log2(2.8), 5); % peak SF
                        if find(iLocPair_all==1)
                            y_ticks_all{2} = linspace(.01, .13, 5); % SF peak amp (fovea)
                        else
                            y_ticks_all{2} = linspace(.03, .07, 5); % SF peak amp
                        end
                        y_ticks_all{3} = linspace(.5, .9, 5); % SF bandwidth
                        y_ticks_all{4} = linspace(-.04, .04, 5); % SF baseline

                end
        end

        % Loop through each tuning characteristic
        for iTunC = 1:nTunC_full

            % Compute median for each observer
            % margTunC_ORI_allCond: nModelA x nModelB x nLoc_pair x nsubj x nIter x nDataset x nTunC
            % tunC_allSubj_allIter: nLoc_pair x nsubj x nIter
            switch iFeature
                case 1
                    data_allIter_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                    data_obs_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, 1, 2, iTunC)); % full data, without resampling
                case 2
                    data_allIter_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                    data_obs_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, 1, 2, iTunC));
            end

            % Obtain ytick labels and ref
            y_ticklabels = round(y_ticks_all{iTunC}, 2);
            ref = nan;
            switch iFamily_perF(iFeature)
                case 1, if find(iTunC==3), ref = 0; end

                case 2 % log parabola
                    switch iTunC
                        case 1, y_ticklabels = round(2.^y_ticks_all{iTunC}, 1); ref = log2(2);
                            data_allIter_allSubj = log2(data_allIter_allSubj); % pref SF
                        case 4, ref = 0; % baseline
                    end
            end

            str_title = sprintf('n=%d nIter=%d L%d%d [A%dB%d] [%s] | %s %s', nsubj, nIterxJob, iLocPair_all, iModelA_plot, iModelB_plot, namesDataset{iDataset_plotRC}, namesFeature{iFeature}, namesTunC{iTunC});

            switch flag_plotDist
                case 0
                    %------------------------------%
                    basicFxn_drawBars_permutation(data_allIter_allSubj, ref, colors, x_ticks, y_ticks_all{iTunC}, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj)
                    % ------------------------------%
                    ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily, 2}{iTunC}))
                case 1
                    %------------------------------%
                    basicFxn_drawDist_permutation(data_allIter_allSubj, ref, colors, x_ticks, y_ticks_all{iTunC}, y_ticklabels, str_title, sz_fig, nIterxJob, nsubj)
                    %------------------------------%
                    xlabel(sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily, 2}{iTunC}))
                    ylabel('Probabillity')
            end
            % Save the figure
            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_A%dB%d_%s%d.png', nameFolder_Fig_tunC, nsubj, iLocPair_all, iModelA_plot, iModelB_plot, namesFeature{iFeature}, iTunC))
            close(gcf)

        end % end of iTunC
    end % end of iFeature
end % end of iGroup
clear margTunC_*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Corr0: Corr between CS and tunParams (just to check bound-hitting)
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'margParams*_allCond')

nameVarX = 'CS';
nameVarY = 'tunParam';
sz_label = 55;
sz_labelOffset = .05;
sz_axOffset = 0.05; % extra breathing room

nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

flag_zeroMean = 0;
flag_plotIdvdCI = 0;
flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles
type_corr = 'pearson';
type_tail = 'both';

for iSet = 1:length(iLocSingle_allSets)

    iLocCorr_all = iLocSingle_allSets{iSet};

    % Obtain CS (x-axis)
    X_allSubj = CS_allSubj(:, iLocCorr_all);
    x_ticks = linspace(1.5, 3.5, 5); % EE
    X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);

    for iFeature = 1:nFeatures
        iFamily = iFamily_perF(iFeature);
        namesTunCs = namesTunC_unit_perF{iFamily, 1};
        nTunCs_full = length(namesTunCs);

        for iTunC = 1:nTunCs_full
            switch iFeature
                case 1, NOMp_allIter_allSubj = squeeze(margParams_ORI_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
                case 2, NOMp_allIter_allSubj = squeeze(margParams_SF_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
            end

            switch iFeature
                case 1, y_ticks_allTunC_lb = [0, 0, -.1]; y_ticks_allTunC_ub = [.36, 90, .1];% 3 values are ORI peak amplitude, width, baseline
                case 2, y_ticks_allTunC_lb = [0, 0, 0, -.1]; y_ticks_allTunC_ub = [4, .24, 1, .1]; % 4 values are SF peak, peak amplitude, width, baseline
            end

            % y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);
            y_ticks = [];

            x_ticklabels = x_ticks;
            y_ticklabels = y_ticks;

            nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunCs{iTunC});
            nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

            str_title = sprintf('n=%d, nIter=%d %s vs. %s [L%s]', nsubj, nIterxJob, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

            basicFxn_drawCorr_permutation(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, type_corr, type_tail, str_title, markers_allSubj, nIterxJob);
            % plot ub and lb when fitting
            yline(ub_full_all{iFamily}(iTunC), 'k--');
            yline(lb_full_all{iFamily}(iTunC), 'k--');

            xlabel('Contrast sensitivity (1/contrast)', 'fontsize', sz_label)
            ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunCs{iTunC}), 'fontsize', sz_label)

            % Adjust distance between components
            ax = gca;
            ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
            ax.XLabel.Units = 'normalized';
            ax.YLabel.Units = 'normalized';

            ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
            ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left

            ax.Position = [ ...
                ti(1) + sz_axOffset, ...
                ti(2) + sz_axOffset, ...
                1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                1 - ti(2) - ti(4) - 2*sz_axOffset];

            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_%s_L%s.png', nameFolder_Fig_NOM_corr, nsubj, nameVarY_fileTitle, strjoin(string(iLocCorr_all), '')))
            close(gcf)

        end % iTunC
    end % iFeature
end % iSet
clear margParams*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Corr1: Corr between CS and tunC
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'margTunC*_allCond')

nameVarX = 'CS';
nameVarY = 'tunC';
sz_label = 55;
sz_labelOffset = .05;
sz_axOffset = 0.05; % extra breathing room

nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

flag_zeroMean = 0;
flag_plotIdvdCI = 0;
flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles
type_corr = 'pearson';
type_tail = 'both';

for iSet = 1:numel(iLocSingle_allSets)

    iLocCorr_all = iLocSingle_allSets{iSet};

    % Obtain CS (x-axis)
    X_allSubj = CS_allSubj(:, iLocCorr_all);
    x_ticks = linspace(1.5, 3.5, 5); % EE
    X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);

    for iFeature = 1:nFeatures
        iFamily = iFamily_perF(iFeature);
        namesTunCs = namesTunC_unit_perF{iFamily, 2};
        nTunCs_full = length(namesTunCs);

        for iTunC = 1:nTunCs_full
            switch iFeature
                case 1, NOMp_allIter_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
                case 2, NOMp_allIter_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iDataset_plotRC, iTunC));
            end

            switch iFeature
                case 1, y_ticks_allTunC_lb = [0, 0, -.1]; y_ticks_allTunC_ub = [.36, 90, .1];% 3 values are ORI peak amplitude, width, baseline
                case 2, y_ticks_allTunC_lb = [0, 0, 0, -.1]; y_ticks_allTunC_ub = [4, .24, 2.8, .1]; % 4 values are SF peak, peak amplitude, width, baseline
            end

            y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);

            x_ticklabels = x_ticks;
            y_ticklabels = y_ticks;

            nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily_perF(iFeature), 2}{iTunC});
            nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

            str_title = sprintf('n=%d, nIter=%d %s vs. %s [L%s]', nsubj, nIterxJob, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

            %======================%
            basicFxn_drawCorr_permutation(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, type_corr, type_tail, str_title, markers_allSubj, nIterxJob);
            %======================%

            xlabel('Contrast sensitivity (1/contrast)', 'fontsize', sz_label)
            ylabel(sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily_perF(iFeature), 2}{iTunC}), 'fontsize', sz_label)

            % Adjust distance between components
            ax = gca;
            ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
            ax.XLabel.Units = 'normalized';
            ax.YLabel.Units = 'normalized';

            ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
            ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left

            ax.Position = [ ...
                ti(1) + sz_axOffset, ...
                ti(2) + sz_axOffset, ...
                1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                1 - ti(2) - ti(4) - 2*sz_axOffset];

            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_%s_L%s.png', nameFolder_Fig_NOM_corr, nsubj, nameVarY_fileTitle, strjoin(string(iLocCorr_all), '')))
            close(gcf)

        end % iTunC
    end % iFeature
end % iSet
clear margTunC*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Corr2: Corr between CS and NOM params
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond')

nameVarX = 'CS';
nameVarY = 'NOMparams';

sz_label = 55;
sz_labelOffset = .05;
sz_axOffset = 0.05; % extra breathing room

nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

flag_zeroMean = 0;
flag_plotIdvdCI = 0;
flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles
type_corr = 'pearson';
type_tail = 'both';

for iSet = 1:numel(iLocSingle_allSets)

    iLocCorr_all = iLocSingle_allSets{iSet};

    % Obtain CS (x-axis)
    X_allSubj = CS_allSubj(:, iLocCorr_all);
    x_ticks = linspace(1.5, 3.5, 5); % EE
    X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);
    % X_med_allSubj = X_allSubj;

    for iNOMparam = 1:nNOMparams
        NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iNOMparam));

        Y_med_allSubj = getCI(NOMp_allIter_allSubj, 1, 3)'; % rotate to match the format needed by basicFxn_drawCorr

        y_ticks_lb = [0, 0, 0]; y_ticks_ub = [.8, 40, 40];

        y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

        x_ticklabels = x_ticks;
        y_ticklabels = y_ticks;

        nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
        nameVarY_fileTitle = nameVarY_figTitle;

        str_title = sprintf('n=%d, nIter=%d %s vs. %s [L%s]', nsubj, nIterxJob, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

        basicFxn_drawCorr_permutation(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, type_corr, type_tail, str_title, markers_allSubj, nIterxJob);

        xlabel('Contrast sensitivity (1/contrast)', 'fontsize', sz_label)
        ylabel(nameVarY_figTitle, 'fontsize', sz_label)

        % Adjust distance between components
        ax = gca;
        ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
        ax.XLabel.Units = 'normalized';
        ax.YLabel.Units = 'normalized';

        ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
        ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left

        ax.Position = [ ...
            ti(1) + sz_axOffset, ...
            ti(2) + sz_axOffset, ...
            1 - ti(1) - ti(3) - 2*sz_axOffset, ...
            1 - ti(2) - ti(4) - 2*sz_axOffset];

        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_NOM%d_L%s.png', nameFolder_Fig_NOM_corr, nsubj, iNOMparam, strjoin(string(iLocCorr_all), '')))
        close(gcf)
    end % iNOMparam
end % iSet
clear params_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Corr3: Corr between CS and pA
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'metric_data_allCond')

nameVarX = 'CS';
nameVarY = 'pA'; iMetric_pA = 2;

sz_label = 55;
sz_labelOffset = .05;
sz_axOffset = 0.05; % extra breathing room

nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

flag_zeroMean = 0;
flag_plotIdvdCI = 0;
flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles
type_corr = 'pearson';
type_tail = 'both';

for iSet = 1:numel(iLocSingle_allSets)

    iLocCorr_all = iLocSingle_allSets{iSet};

    % Obtain CS (x-axis)
    X_allSubj = CS_allSubj(:, iLocCorr_all);
    x_ticks = linspace(1.5, 3.5, 5); % EE
    X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);

    NOMp_allIter_allSubj = getCI(metric_data_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, iMetric_pA, :, :, :), 2, 7);

    y_ticks = round(linspace(.6, .9, 5), 2);

    x_ticklabels = x_ticks;
    y_ticklabels = y_ticks;

    nameVarY_figTitle = nameVarY;
    nameVarY_fileTitle = nameVarY_figTitle;

    str_title = sprintf('n=%d, nIter=%d%s vs. %s [L%s]', nsubj, nIterxJob, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

    basicFxn_drawCorr_permutation(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, type_corr, type_tail, str_title, markers_allSubj, nIterxJob);

    xlabel('Contrast sensitivity (1/contrast)', 'fontsize', sz_label)
    ylabel(nameVarY_figTitle, 'fontsize', sz_label)

    % Adjust distance between components
    ax = gca;
    ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
    ax.XLabel.Units = 'normalized';
    ax.YLabel.Units = 'normalized';

    ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
    ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left

    ax.Position = [ ...
        ti(1) + sz_axOffset, ...
        ti(2) + sz_axOffset, ...
        1 - ti(1) - ti(3) - 2*sz_axOffset, ...
        1 - ti(2) - ti(4) - 2*sz_axOffset];

    drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%s.png', nameFolder_Fig_NOM_corr, nsubj, strjoin(string(iLocCorr_all), '')))
    close(gcf)

end % iSet
clear metric_data_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Corr4: Corr between pA and NOM params
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond')
nameVarX = 'pA';
nameVarY = 'NOMparams';

sz_label = 55;
sz_labelOffset = .05;
sz_axOffset = 0.05; % extra breathing room

nameFolder_Fig_NOM_corr = sprintf('%s/Corr/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_corr)); mkdir(nameFolder_Fig_NOM_corr), end

flag_zeroMean = 0;
flag_plotIdvdCI = 0;
flag_plotUnikSymbol = 0; % 1=each subj has a unique marker; 0=all are circles
type_corr = 'pearson';
type_tail = 'both';

for iSet = 1:numel(iLocSingle_allSets)

    iLocCorr_all = iLocSingle_allSets{iSet};

    % Obtain x-axis values
    X_allSubj = pA_allSubj(:, iLocCorr_all);
    x_ticks = linspace(.6, .8, 5); % EE
    X_allIter_allSubj = repmat(X_allSubj, 1, 1, nIterxJob);

    for iNOMparam = 1:nNOMparams
        NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocCorr_all, :, :, iNOMparam));

        y_ticks_lb = [0, 5, 0]; y_ticks_ub = [.8, 45, 40];

        y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

        x_ticklabels = x_ticks;
        y_ticklabels = y_ticks;

        nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
        nameVarY_fileTitle = nameVarY_figTitle;

        str_title = sprintf('n=%d, nIter=%d%s vs. %s [L%s]', nsubj, nIterxJob, nameVarX, nameVarY_figTitle, strjoin(string(iLocCorr_all), ''));

        basicFxn_drawCorr_permutation(X_allIter_allSubj, NOMp_allIter_allSubj, colors_comb(iLocCorr_all, :), x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, flag_plotUnikSymbol, type_corr, type_tail, str_title, markers_allSubj, nIterxJob);

        xlabel(nameVarX, 'fontsize', sz_label)
        ylabel(nameVarY_figTitle, 'fontsize', sz_label)

        % Adjust distance between components
        ax = gca;
        ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
        ax.XLabel.Units = 'normalized';
        ax.YLabel.Units = 'normalized';

        ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
        ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left

        ax.Position = [ ...
            ti(1) + sz_axOffset, ...
            ti(2) + sz_axOffset, ...
            1 - ti(1) - ti(3) - 2*sz_axOffset, ...
            1 - ti(2) - ti(4) - 2*sz_axOffset];

        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_NOM%d_L%s.png', nameFolder_Fig_NOM_corr, nsubj, iNOMparam, strjoin(string(iLocCorr_all), '')))
        close(gcf)
    end % iNOMparam
end % iSet
clear params_allCond
fprintf('\n\n Plotting DONE\n\n')

%% CompAsym1: (CS and tunC): bin and compare between extents of EE/HVA/VMA
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'margTunC_*_allCond')

nameVarX = 'CS';
nameVarY = 'tunC';

nameFolder_Fig_NOM_CompAsym = sprintf('%s/CompAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_CompAsym)); mkdir(nameFolder_Fig_NOM_CompAsym), end

sz_label = 15;
sz_labelOffset = 0.01;
sz_axOffset = 0.05; % extra breathing room
nBinsCompAsym = 2;
sz_wd_perBar = 200;
sz_fig = [1e3, 500];

for iGroup = 1:nGroups
    iLocPair_all = iLocGroups_all{iGroup};

    switch iLocPair_all(1)
        case 1, nameAsymX = 'Ecc. effect'; nameAsymY = 'Ecc. effect';
        case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
        case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
    end

    % X-axis
    asymX_allSubj = (CS_allSubj(:, iLocPair_all(1))-CS_allSubj(:, iLocPair_all(2)))./(CS_allSubj(:, iLocPair_all(1))+CS_allSubj(:, iLocPair_all(2)));
    switch iLocPair_all(1)
        case 1, x_ticks = linspace(0, 20, 5); % EE
        case 6, x_ticks = linspace(-5, 15, 5); % HVA
        case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
    end
    asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob); % [nsubj x nIterxJob]

    for iFeature = 1:nFeatures
        iFamily = iFamily_perF(iFeature);
        namesTunCs = namesTunC_unit_perF{iFamily, 2};
        nTunCs_full = length(namesTunCs);

        for iTunC = 1:nTunCs_full
            % Y-axis
            switch iFeature
                case 1, NOMp_allIter_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                case 2, NOMp_allIter_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
            end

            % Y-value: [nsubj x nIterxJob]
            asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

            switch iLocPair_all(1)
                case 1 % EE
                    switch iFamily
                        case 1, y_ticks_allTunC_lb = -[0, 32, 100]; y_ticks_allTunC_ub = [72, 20, 100];
                            % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                        case 2, y_ticks_allTunC_lb = -[50, 0, 50, 145]; y_ticks_allTunC_ub = [50, 80, 50, 265];
                    end
                case 6 % HVA
                    switch iFamily
                        case 1, y_ticks_allTunC_lb = -[40, 30, 100]; y_ticks_allTunC_ub = [40, 50, 100];
                            % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                        case 2, y_ticks_allTunC_lb = -[50, 20, 60, 145]; y_ticks_allTunC_ub = [50, 60, 60, 265];
                    end
                case 5 % VMA
                    switch iFamily
                        case 1, y_ticks_allTunC_lb = -[30, 70, 50]; y_ticks_allTunC_ub = [90, 30, 50];
                            % case 8, y_ticks_allTunC_lb = -[100, 40, 300, 320, 50, 320]; y_ticks_allTunC_ub = [100, 44, 300, 320, 50, 320];
                        case 2, y_ticks_allTunC_lb = -[50, 40, 60, 250]; y_ticks_allTunC_ub = [50, 100, 40, 250];
                    end
            end
            y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);

            x_ticklabels = nan;
            y_ticklabels = nan;

            nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily_perF(iFeature), 2}{iTunC});
            nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

            str_title = sprintf('n=%d, nIter=%d %s (%s) vs. %s (%s)', nsubj, nIterxJob, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

            % remove later!!s
            asymX_allIter_allSubj_ = repmat(asymX_allIter_allSubj, 5, 1);
            asymX_allIter_allSubj_ = asymX_allIter_allSubj_+randn(size(asymX_allIter_allSubj_))*mean(asymX_allIter_allSubj_(:))/20;
            asymY_allIter_allSubj_ = repmat(asymY_allIter_allSubj, 5, 1);
            asymY_allIter_allSubj_ = asymY_allIter_allSubj_+randn(size(asymY_allIter_allSubj_))*mean(asymY_allIter_allSubj_(:))/20;

            str_ylabel = sprintf('\\Delta %s %s (%%)', namesFeature{iFeature}, namesTunC_noUnit{iFamily_perF(iFeature), 2}{iTunC});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            basicFxn_compAsym_permutation(asymX_allIter_allSubj_*100, asymY_allIter_allSubj_*100, nBinsCompAsym, y_ticks, sz_fig, str_title, str_ylabel)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % xlabel(sprintf('\\Delta contrast sensitivity (%%)'), 'FontSize', sz_label);
            % ylabel(str_ylabel, 'fontsize', sz_label)

            % Adjust distance between components
            % ax = gca;
            % ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
            %
            % ax.XLabel.Units = 'normalized';
            % ax.YLabel.Units = 'normalized';
            %
            % ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
            % ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left
            % ax.Position = [ ...
            %     ti(1) + sz_axOffset, ...
            %     ti(2) + sz_axOffset, ...
            %     1 - ti(1) - ti(3) - 2*sz_axOffset, ...
            %     1 - ti(2) - ti(4) - 2*sz_axOffset];

            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_%s.png', nameFolder_Fig_NOM_CompAsym, nsubj, iLocPair_all, nameVarY_fileTitle))
            close(gcf)

        end % iTunC
    end % iFeature
end % iGroup
clear margTunC_*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% CompAsym2: CorrAsym (CS and NOMparams): corr between extents of EE/HVA/VMA
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond')

nameVarX = 'CS';
nameVarY = 'NOMparams';

nameFolder_Fig_NOM_CompAsym = sprintf('%s/CompAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_CompAsym)); mkdir(nameFolder_Fig_NOM_CompAsym), end

sz_label = 15;
sz_labelOffset = 0.05;
sz_axOffset = 0.05; % extra breathing room

for iGroup = 1:nGroups
    iLocPair_all = iLocGroups_all{iGroup};

    switch iLocPair_all(1)
        case 1, nameAsymX = 'Ecc. effect'; nameAsymY = 'Ecc. effect';
        case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
        case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
    end

    asymX_allSubj = (CS_allSubj(:, iLocPair_all(1))-CS_allSubj(:, iLocPair_all(2)))./(CS_allSubj(:, iLocPair_all(1))+CS_allSubj(:, iLocPair_all(2)));
    switch iLocPair_all(1)
        case 1, x_ticks = linspace(0, 20, 5); % EE
        case 6, x_ticks = linspace(-5, 15, 5); % HVA
        case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
    end
    asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';

    for iNOMparam = 1:nNOMparams
        NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iNOMparam));

        asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

        y_ticks_lb = -[30, 20, 20]; y_ticks_ub = [50 40 60];
        y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

        x_ticklabels = nan;
        y_ticklabels = nan;

        nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
        nameVarY_fileTitle = nameVarY_figTitle;

        str_title = sprintf('n=%d, nIter=%d %s (%s) vs. %s (%s)', nsubj, nIterxJob, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

        % remove later!!s
        asymX_allIter_allSubj_ = repmat(asymX_allIter_allSubj', 5, 1);
        asymX_allIter_allSubj_ = asymX_allIter_allSubj_+randn(size(asymX_allIter_allSubj_))*mean(asymX_allIter_allSubj_(:))/20;
        asymY_allIter_allSubj_ = repmat(asymY_allIter_allSubj, 5, 1);
        asymY_allIter_allSubj_ = asymY_allIter_allSubj_+randn(size(asymY_allIter_allSubj_))*mean(asymY_allIter_allSubj_(:))/20;

        str_ylabel = sprintf('\\Delta %s (%%)', nameVarY_figTitle);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        basicFxn_compAsym_permutation(asymX_allIter_allSubj_*100, asymY_allIter_allSubj_*100, nBinsCompAsym, y_ticks, sz_fig, str_title, str_ylabel)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % xlabel(sprintf('\\Delta contrast sensitivity (%%)'), 'FontSize', sz_label);
        % ylabel(, 'fontsize', sz_label)

        % Adjust distance between components
        % ax = gca;
        % ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks
        %
        % ax.XLabel.Units = 'normalized';
        % ax.YLabel.Units = 'normalized';
        %
        % ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
        % ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left
        % ax.Position = [ ...
        %     ti(1) + sz_axOffset, ...
        %     ti(2) + sz_axOffset, ...
        %     1 - ti(1) - ti(3) - 2*sz_axOffset, ...
        %     1 - ti(2) - ti(4) - 2*sz_axOffset];

        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_NOM%d.png', nameFolder_Fig_NOM_CompAsym, nsubj, iLocPair_all, iNOMparam))
        close(gcf)

    end % iNOMparam
end % iGroup
clear params_allCond
fprintf('\n\n Plotting DONE\n\n')


%% CorrAsym1: (CS and tunC): corr between extents of EE/HVA/VMA
flag_plotIdvdCI = 1;
flag_plotUnikSymbol = 0;

clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'margTunC_*_allCond')

nameVarX = 'CS';
nameVarY = 'tunC';

nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

sz_label = 55;
sz_labelOffset = 0.01;
sz_axOffset = 0.05; % extra breathing room

for iGroup = 1:nGroups
    iLocPair_all = iLocGroups_all{iGroup};

    switch iLocPair_all(1)
        case 1, nameAsymX = 'Ecc. effect'; nameAsymY = 'Ecc. effect';
        case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
        case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
    end

    % X-axis
    asymX_allSubj = (CS_allSubj(:, iLocPair_all(1))-CS_allSubj(:, iLocPair_all(2)))./(CS_allSubj(:, iLocPair_all(1))+CS_allSubj(:, iLocPair_all(2)));
    switch iLocPair_all(1)
        case 1, x_ticks = linspace(0, 20, 5); % EE
        case 6, x_ticks = linspace(-5, 15, 5); % HVA
        case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
    end
    asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';

    for iFeature = 1:nFeatures
        iFamily = iFamily_perF(iFeature);
        namesTunCs = namesTunC_unit_perF{iFamily, 2};
        nTunCs_full = length(namesTunCs);

        for iTunC = 1:nTunCs_full
            % Y-axis
            switch iFeature
                case 1, NOMp_allIter_allSubj = squeeze(margTunC_ORI_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
                case 2, NOMp_allIter_allSubj = squeeze(margTunC_SF_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iDataset_plotRC, iTunC));
            end

            asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

            switch iLocPair_all(1)
                case 1 % EE
                    switch iFamily
                        case 1, y_ticks_allTunC_lb = -[0, 32, 100]; y_ticks_allTunC_ub = [72, 20, 100];
                            % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                        case 2, y_ticks_allTunC_lb = -[50, 0, 50, 145]; y_ticks_allTunC_ub = [50, 80, 50, 265];
                    end
                case 6 % HVA
                    switch iFamily
                        case 1, y_ticks_allTunC_lb = -[40, 30, 100]; y_ticks_allTunC_ub = [40, 50, 100];
                            % case 8, y_ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; y_ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                        case 2, y_ticks_allTunC_lb = -[50, 20, 60, 145]; y_ticks_allTunC_ub = [50, 60, 60, 265];
                    end
                case 5 % VMA
                    switch iFamily
                        case 1, y_ticks_allTunC_lb = -[30, 70, 50]; y_ticks_allTunC_ub = [90, 30, 50];
                            % case 8, y_ticks_allTunC_lb = -[100, 40, 300, 320, 50, 320]; y_ticks_allTunC_ub = [100, 44, 300, 320, 50, 320];
                        case 2, y_ticks_allTunC_lb = -[50, 40, 60, 250]; y_ticks_allTunC_ub = [50, 100, 40, 250];
                    end
            end
            y_ticks = linspace(y_ticks_allTunC_lb(iTunC), y_ticks_allTunC_ub(iTunC), 5);

            x_ticklabels = nan;
            y_ticklabels = nan;

            nameVarY_figTitle = sprintf('%s %s', namesFeature{iFeature}, namesTunC_unit_perF{iFamily_perF(iFeature), 2}{iTunC});
            nameVarY_fileTitle = sprintf('%s%d', namesFeature{iFeature}, iTunC);

            str_title = sprintf('n=%d, nIter=%d %s (%s) vs. %s (%s)', nsubj, nIterxJob, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            basicFxn_drawCorrAsym_permutation(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            xlabel(sprintf('\\Delta contrast sensitivity (%%)'), 'FontSize', sz_label);
            ylabel(sprintf('\\Delta %s %s (%%)', namesFeature{iFeature}, namesTunC_noUnit{iFamily_perF(iFeature), 2}{iTunC}), 'fontsize', sz_label)

            % Adjust distance between components
            ax = gca;
            ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks

            ax.XLabel.Units = 'normalized';
            ax.YLabel.Units = 'normalized';

            ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
            ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left
            ax.Position = [ ...
                ti(1) + sz_axOffset, ...
                ti(2) + sz_axOffset, ...
                1 - ti(1) - ti(3) - 2*sz_axOffset, ...
                1 - ti(2) - ti(4) - 2*sz_axOffset];

            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_%s.png', nameFolder_Fig_NOM_CorrAsym, nsubj, iLocPair_all, nameVarY_fileTitle))
            close(gcf)

        end % iTunC
    end % iFeature
end % iGroup
clear margTunC_*_allCond
fprintf('\n\n Plotting DONE\n\n')

%% CorrAsym2: CorrAsym (CS and NOMparams): corr between extents of EE/HVA/VMA
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond')

nameVarX = 'CS';
nameVarY = 'NOMparams';

nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

sz_label = 60;
sz_labelOffset = 0.05;
sz_axOffset = 0.05; % extra breathing room

for iGroup = 1:nGroups
    iLocPair_all = iLocGroups_all{iGroup};

    switch iLocPair_all(1)
        case 1, nameAsymX = 'Ecc. effect'; nameAsymY = 'Ecc. effect';
        case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
        case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
    end

    asymX_allSubj = (CS_allSubj(:, iLocPair_all(1))-CS_allSubj(:, iLocPair_all(2)))./(CS_allSubj(:, iLocPair_all(1))+CS_allSubj(:, iLocPair_all(2)));
    switch iLocPair_all(1)
        case 1, x_ticks = linspace(0, 20, 5); % EE
        case 6, x_ticks = linspace(-5, 15, 5); % HVA
        case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
    end
    asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';

    for iNOMparam = 1:nNOMparams
        NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iNOMparam));

        asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

        y_ticks_lb = -[30, 20, 20]; y_ticks_ub = [50 40 40];
        y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

        x_ticklabels = nan;
        y_ticklabels = nan;

        nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
        nameVarY_fileTitle = nameVarY_figTitle;

        str_title = sprintf('n=%d, nIter=%d %s (%s) vs. %s (%s)', nsubj, nIterxJob, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        basicFxn_drawCorrAsym_permutation(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        xlabel(sprintf('\\Delta contrast sensitivity (%%)'), 'FontSize', sz_label);
        ylabel(sprintf('\\Delta %s (%%)', nameVarY_figTitle), 'fontsize', sz_label)

        % Adjust distance between components
        ax = gca;
        ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks

        ax.XLabel.Units = 'normalized';
        ax.YLabel.Units = 'normalized';

        ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
        ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left
        ax.Position = [ ...
            ti(1) + sz_axOffset, ...
            ti(2) + sz_axOffset, ...
            1 - ti(1) - ti(3) - 2*sz_axOffset, ...
            1 - ti(2) - ti(4) - 2*sz_axOffset];

        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_NOM%d.png', nameFolder_Fig_NOM_CorrAsym, nsubj, iLocPair_all, iNOMparam))
        close(gcf)

    end % iNOMparam
end % iGroup
clear params_allCond
fprintf('\n\n Plotting DONE\n\n')

%% CorrAsym3: (CS and pA): corr between extents of EE/HVA/VMA
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'metric_data_allCond')

nameVarX = 'CS';
nameVarY = 'pA';

iMetric_pA = 2;

nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

sz_label = 60;
sz_labelOffset = 0.05;
sz_axOffset = 0.05; % extra breathing room

for iGroup = 1:nGroups
    iLocPair_all = iLocGroups_all{iGroup};

    switch iLocPair_all(1)
        case 1, nameAsymX = 'EE'; nameAsymY = 'EE';
        case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
        case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
    end

    asymX_allSubj = (CS_allSubj(:, iLocPair_all(1))-CS_allSubj(:, iLocPair_all(2)))./(CS_allSubj(:, iLocPair_all(1))+CS_allSubj(:, iLocPair_all(2)));
    switch iLocPair_all(1)
        case 1, x_ticks = linspace(0, 20, 5); % EE
        case 6, x_ticks = linspace(-5, 15, 5); % HVA
        case 5, x_ticks = linspace(0, 12, 5); % VMA (extent is smaller)
    end
    asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';

    NOMp_allIter_allSubj = getCI(metric_data_allCond(iModelA_plot, iModelB_plot, iLocPair_all, iMetric_pA, :, :, :), 2, 7);

    asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

    % switch iLocPair_all(1)
    %     case 1 % EE
    %         y_ticks_lb = -[10, 70, 12]; y_ticks_ub = [6, 40, 30];
    %     case 6 % HVA
    %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
    %     case 5 % VMA
    %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
    % end

    y_ticks = linspace(-10, 10, 5);

    x_ticklabels = nan;
    y_ticklabels = nan;

    nameVarY_figTitle = nameVarY;
    nameVarY_fileTitle = nameVarY_figTitle;

    str_title = sprintf('n=%d, nIter=%d %s (%s) vs. %s (%s)', nsubj, nIterxJob, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    basicFxn_drawCorrAsym_permutation(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xlabel(sprintf('%s of contrast sensitivity (%%)', nameAsymX), 'fontsize', sz_label)
    ylabel(sprintf('%s of %s (%%)', nameAsymX, nameVarY_figTitle), 'fontsize', sz_label)

    % Adjust distance between components
    ax = gca;
    ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks

    ax.XLabel.Units = 'normalized';
    ax.YLabel.Units = 'normalized';

    ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
    ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left
    ax.Position = [ ...
        ti(1) + sz_axOffset, ...
        ti(2) + sz_axOffset, ...
        1 - ti(1) - ti(3) - 2*sz_axOffset, ...
        1 - ti(2) - ti(4) - 2*sz_axOffset];

    drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d.png', nameFolder_Fig_NOM_CorrAsym, nsubj, iLocPair_all))
    close(gcf)

end % iGroup
clear metric_data_allCond
fprintf('\n\n Plotting DONE\n\n')

%% CorrAsym4: (pA and NOMparams): corr between extents of EE/HVA/VMA
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond')

nameVarX = 'pA';
nameVarY = 'NOMparams';

nameFolder_Fig_NOM_CorrAsym = sprintf('%s/CorrAsym/%s_%s', nameFolder_Fig_NOM_Trialwise, nameVarX, nameVarY);
if isempty(dir(nameFolder_Fig_NOM_CorrAsym)); mkdir(nameFolder_Fig_NOM_CorrAsym), end

sz_label = 60;
sz_labelOffset = 0.05;
sz_axOffset = 0.05; % extra breathing room

for iGroup = 1:nGroups
    iLocPair_all = iLocGroups_all{iGroup};

    switch iLocPair_all(1)
        case 1, nameAsymX = 'EE'; nameAsymY = 'EE';
        case 6, nameAsymX = 'HVA'; nameAsymY = 'HVA';
        case 5, nameAsymX = 'VMA'; nameAsymY = 'VMA';
    end

    % X-axis (pA)
    asymX_allSubj = (pA_allSubj(:, iLocPair_all(1))-pA_allSubj(:, iLocPair_all(2)))./(pA_allSubj(:, iLocPair_all(1))+pA_allSubj(:, iLocPair_all(2)));
    asymX_allIter_allSubj = repmat(asymX_allSubj, 1, nIterxJob)';
    % switch iLocPair_all(1)
    %     case 1, x_ticks = linspace(0, 20, 5); % EE
    %     case 6, x_ticks = linspace(-5, 15, 5); % HVA
    %     case 5, x_ticks = linspace(0, 16, 5); % VMA (extent is smaller)
    % end
    x_ticks = linspace(-6, 10, 5);

    for iNOMparam = 1:nNOMparams
        NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB_plot, iLocPair_all, :, :, iNOMparam));

        asymY_allIter_allSubj = squeeze((NOMp_allIter_allSubj(1, :, :)-NOMp_allIter_allSubj(2, :, :))./(NOMp_allIter_allSubj(1, :, :)+NOMp_allIter_allSubj(2, :, :)));

        % switch iLocPair_all(1)
        %     case 1 % EE
        %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        %     case 6 % HVA
        %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        %     case 5 % VMA
        %         y_ticks_lb = -[70, 70, 12]; y_ticks_ub = [30, 40, 30];
        % end
        y_ticks_lb = -[60, 40, 20]; y_ticks_ub = [40, 40, 60];
        y_ticks = linspace(y_ticks_lb(iNOMparam), y_ticks_ub(iNOMparam), 5);

        x_ticklabels = nan;
        y_ticklabels = nan;

        nameVarY_figTitle = namesModelBparams{iModelB_plot}{iNOMparam};
        nameVarY_fileTitle = nameVarY_figTitle;

        str_title = sprintf('n=%d, nIter=%d %s (%s) vs. %s (%s)', nsubj, nIterxJob, nameVarX, nameAsymX, nameVarY_figTitle, nameAsymY);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        basicFxn_drawCorrAsym_permutation(asymX_allIter_allSubj*100, asymY_allIter_allSubj*100, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_plotIdvdCI, flag_plotUnikSymbol, str_title, markers_allSubj)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        xlabel(sprintf('%s of %s (%%)', nameAsymX, nameVarX), 'fontsize', sz_label)
        ylabel(sprintf('%s of %s (%%)', nameAsymX, nameVarY_figTitle), 'fontsize', sz_label)

        % Adjust distance between components
        ax = gca;
        ti = ax.TightInset;   % [left bottom right top] padding needed for labels/ticks

        ax.XLabel.Units = 'normalized';
        ax.YLabel.Units = 'normalized';

        ax.XLabel.Position(2) = ax.XLabel.Position(2) - sz_labelOffset;   % move label down
        ax.YLabel.Position(1) = ax.YLabel.Position(1) - sz_labelOffset;   % move label left
        ax.Position = [ ...
            ti(1) + sz_axOffset, ...
            ti(2) + sz_axOffset, ...
            1 - ti(1) - ti(3) - 2*sz_axOffset, ...
            1 - ti(2) - ti(4) - 2*sz_axOffset];

        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_NOM%d.png', nameFolder_Fig_NOM_CorrAsym, nsubj, iLocPair_all, iNOMparam))
        close(gcf)

    end % iNOMparam
end % iGroup
clear params_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Plot metrics vs. IV for group averages (single locations)
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'IV_allCond', 'nTrials_allCond', 'metric_*_allCond', 'R2_NOM_allCond')

% Define folder for saving figures
nameFolder_Fig_NOM_Trialwise = sprintf('%s/NOM_Trialwise_%d%d', nameFolder_Figures, nORI, nSF);
if isempty(dir(nameFolder_Fig_NOM_Trialwise)), mkdir(nameFolder_Fig_NOM_Trialwise), end
namesModelB_plot = {'Full model', 'No $\sigma_{shared}$', 'No $\sigma_{mul}$', 'No $\sigma_{private}$'};

x_base   = 0.60;   % move a bit further left to make space for name column
x_width  = 0.4;   % total width of the table block
y_base   = 0.10;
y_height = 0.25;

sz_font = 15;
sz_axis = 18;
lineStyle_all = {'-', '-', '--', ':', ':', '-', '--', ':', '-.'};
namesPlotMode = {'NoPred', 'FullOnly', 'AllModels'};
szScaling = 10; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end
szBase = 5;

for iModelA = 1%iModelA_all % just plot pred of the model using RC-derived template (happens to be yhe best model)

    % Define folder for saving figures for each model A and model B
    nameFolder_Fig_NOM_metrics = sprintf('%s/A%d', nameFolder_Fig_NOM_Trialwise, iModelA);
    if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end


    for iSet = 1:numel(iLocSingle_allSets)

        iLocSingle_all = iLocSingle_allSets{iSet};

        fprintf('\nL%s...', strjoin(string(iLocSingle_all), ''))

        for iPlotMode = 1:length(namesPlotMode)

            % Decide the number of models to plot
            if iPlotMode == 2
                iModelB_all_plot = 1;
            else
                % iModelB_all_plot = 1:4;
                iModelB_all_plot = iModelB_all;
            end
            nModelB_plot = length(iModelB_all_plot);


            %%% --- PLOT ---
            figure('Position', [0 0 500 1e3])

            for iMetric_prob = 1:nMetrics_prob

                subplot(nMetrics_prob, 1, iMetric_prob), hold on

                R2_tab = nan(nModelB_plot, 2);   % [row = model, col = loc within pair]

                for iModelB = iModelB_all_plot

                    if iMetric_prob==1
                        fprintf('%s: %s | ', namesModelB{iModelB_all_plot(iModelB)}, lineStyle_all{iModelB})
                    end

                    for iiLoc = 1:numel(iLocSingle_all)
                        % Compute medians and CIs
                        [IV_ave, ~, ~, IV_SEM] = getCI(getCI(IV_allCond(iModelA, iModelB, iLocSingle_all(iiLoc), :, :, :), 1, 5), 2, 1);
                        [nTrials_ave, ~, ~, nTrials_SEM] = getCI(getCI(nTrials_allCond(iModelA, iModelB, iLocSingle_all(iiLoc), :, :, :), 1, 5), 2, 1);
                        [data_ave, ~, ~, data_SEM] = getCI(getCI(metric_data_allCond(iModelA, iModelB, iLocSingle_all(iiLoc), iMetric_prob, :, :, :), 1, 6), 2, 1);
                        [pred_ave, ~, ~, pred_SEM] = getCI(getCI(metric_pred_allCond(iModelA, iModelB, iLocSingle_all(iiLoc), iMetric_prob, :, :, :), 1, 6), 2, 1);
                        [R2_NOM_ave, ~, ~, R2_NOM_SEM] = getCI(getCI(R2_NOM_allCond(iModelA, iModelB, iLocSingle_all(iiLoc), iMetric_prob, :, :), 1, 6), 2, 1);

                        % ---- Store R2 mean into table ----
                        rowIdx = find(iModelB_all_plot == iModelB);  % which row for this model
                        colIdx = iiLoc;                              % col 1 or 2 for the pair
                        R2_tab(rowIdx, colIdx) = R2_NOM_ave;        % store mean R2 (not SEM)

                        if iModelB==1, color_pred = colors_comb(iLocSingle_all(iiLoc), :); lw = 3;
                        else, color_pred='k'; lw = 1.5;
                        end

                        if iPlotMode ~= 1
                            % Plot prediction
                            patch([IV_ave; flip(IV_ave)], [pred_ave-pred_SEM; flip(pred_ave+pred_SEM)], ones(1, 3) / 2, 'FaceAlpha', .3, 'linestyle', 'none', 'handlevisibility', 'off')
                            plot(IV_ave, pred_ave, 'lineStyle', lineStyle_all{iModelB}, 'color', colors_comb(iLocSingle_all(iiLoc), :), 'linewidth', lw)
                        end

                        % Plot measurement
                        % errorbar(IV_ave, data_ave, IV_SEM, '.', 'horizontal', 'CapSize', 0, 'color', colors_comb(iLoc_Sep_all(iiLoc), :), 'handlevisibility', 'off')
                        errorbar(IV_ave, data_ave, data_SEM, '.', 'vertical', 'CapSize', 0, 'color', colors_comb(iLocSingle_all(iiLoc), :), 'handlevisibility', 'off', 'linewidth', lw/1.5)

                        % Plot averaged data of each bin (dot size indicates number of trials)
                        for iBin = 1:nBins
                            facecolor = 'w';
                            plot(IV_ave(iBin), data_ave(iBin), 'o', 'markeredgecolor', colors_comb(iLocSingle_all(iiLoc), :), ...
                                'markerfacecolor', facecolor, 'MarkerSize', nTrials_ave(iBin) / szScaling + szBase, 'LineWidth', lw/1.5, 'LineStyle', 'none', 'handlevisibility', 'off')
                        end
                    end % iiLoc

                    yline(.5, '--', 'linewidth', 2, 'color', ones(1,3)/2);
                    ylim([0, 1])
                    yticks(0:.2:1)
                    ylabel(namesMetrics_prob_full{iMetric_prob})

                    if find(iLocSingle_all==1), x_ticks = linspace(0, 180, 5);
                    else, x_ticks = linspace(0, 100, 5);
                    end
                    xticks(x_ticks)
                    xlim(x_ticks([1, end]))
                    xlabel('Binned decision variable')
                    if isubj == 1, legend(namesModelB(iModelB_all_plot), 'Location', 'best'), end

                end % iModelB

                % ==== Print R2 as a 4x3 table (left column = model name) ====
                ax = gca;

                nRows_subj = nModelB_plot;
                nCols= 3;   % now: [ModelName | Loc1 | Loc2]

                % Create normalized coordinates for each cell
                x_pos = linspace(x_base, x_base + x_width, nCols + 1);
                x_cells = x_pos(1:nCols) + diff(x_pos(1:2))/2;   % center of each column

                y_pos = linspace(y_base, y_base + y_height, nRows_subj + 1);
                y_cells = fliplr(y_pos(1:nRows_subj)) + diff(y_pos(1:2))/2;

                if iPlotMode ~= 1
                    % ---- Column 1: model names (B1, B2, ...) ----
                    for r = 1:nRows_subj
                        name_str = sprintf('%s', namesModelB_plot{iModelB_all_plot(r)});

                        text(x_cells(1), y_cells(r), name_str, 'Units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize', sz_font, 'FontWeight', 'bold', 'interpreter', 'latex');
                        % text(x_cells(1), y_cells(r), '$R^2$=', 'Units','normalized', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', sz_font, 'FontWeight', 'bold', 'interpreter', 'latex');
                    end

                    % ---- Columns 2–3: R2 values ----
                    for r = 1:nRows_subj
                        for c = 1:2   % two locations
                            R2_val = R2_tab(r, c);
                            if isnan(R2_val), continue; end

                            str_cell = sprintf('%.0f%%', R2_val * 100);

                            % text(x_cells(c+1), y_cells(r), str_cell, 'Units','normalized', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', sz_font);
                        end
                    end
                end % if iPlotMode==2
                ax = gca;
                ax.FontSize=sz_axis;
                ax.LineWidth = 1.5; % Adjust the value as desired
            end % iMetric

            sgtitle(sprintf('n=%d [A%d] [L%s] [nIter=%d] %s', nsubj, iModelA, strjoin(string(iLocSingle_all), ''), nIterxJob, namesMetrics_prob{iMetric_prob}))

            % Save the figure
            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%s_A%d_group_%s.png', nameFolder_Fig_NOM_metrics, nsubj, strjoin(string(iLocSingle_all), ''), iModelA_plot, namesPlotMode{iPlotMode}))
            close(gcf)
            fprintf('DONE\n')
        end % iPlotMode
    end % iModelA
end % iSet
clear *_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Plot metrics vs. IV for group averages (paired locations)
clc, fprintf('\n\n Plotting STARTED......\n\n')

% Load data
load(nameFolder_Data_SaveCompile, 'IV_allCond', 'nTrials_allCond', 'metric_*_allCond', 'R2_NOM_allCond')

% Define folder for saving figures
nameFolder_Fig_NOM_Trialwise = sprintf('%s/NOM_Trialwise_%d%d', nameFolder_Figures, nORI, nSF);
if isempty(dir(nameFolder_Fig_NOM_Trialwise)), mkdir(nameFolder_Fig_NOM_Trialwise), end
namesModelB_plot = {'Full model', 'No $\sigma_{shared}$', 'No $\sigma_{mul}$', 'No $\sigma_{private}$'};

x_base   = 0.60;   % move a bit further left to make space for name column
x_width  = 0.4;   % total width of the table block
y_base   = 0.10;
y_height = 0.25;

sz_font = 15;
sz_axis = 18;
lineStyle_all = {'-', '-', '--', ':', ':', '-', '--', ':', '-.'};
namesPlotMode = {'NoPred', 'FullOnly', 'AllModels'};

for iPlotMode = 1:length(namesPlotMode)

    % Decide the number of models to plot
    if iPlotMode == 2
        iModelB_all_plot = 1;
    else
        % iModelB_all_plot = 1:4;
        iModelB_all_plot = iModelB_all;
    end
    nModelB_plot = length(iModelB_all_plot);

    for iModelA = 1%iModelA_all % just plot pred of the model using RC-derived template (happens to be yhe best model)

        % Define folder for saving figures for each model A and model B
        nameFolder_Fig_NOM_metrics = sprintf('%s/A%d', nameFolder_Fig_NOM_Trialwise, iModelA);
        if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end

        % for iLocComb = iLocComb_all
        for iGroup = 1:nGroups
            iLocPair_all = iLocGroups_all{iGroup};

            fprintf('\nL%d%d...', iLocPair_all)
            szScaling = 80; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end

            figure('Position', [0 0 500 1e3])

            for iMetric_prob = 1:nMetrics_prob

                subplot(nMetrics_prob, 1, iMetric_prob), hold on

                R2_tab = nan(nModelB_plot, 2);   % [row = model, col = loc within pair]

                for iModelB = iModelB_all_plot

                    if iMetric_prob==1
                        fprintf('%s: %s | ', namesModelB{iModelB_all_plot(iModelB)}, lineStyle_all{iModelB})
                    end

                    for iiLoc = 1:2
                        % Compute medians and CIs
                        [IV_ave, ~, ~, IV_SEM] = getCI(getCI(IV_allCond(iModelA, iModelB, iLocPair_all(iiLoc), :, :, :), 1, 5), 2, 1);
                        [nTrials_ave, ~, ~, nTrials_SEM] = getCI(getCI(nTrials_allCond(iModelA, iModelB, iLocPair_all(iiLoc), :, :, :), 1, 5), 2, 1);
                        [data_ave, ~, ~, data_SEM] = getCI(getCI(metric_data_allCond(iModelA, iModelB, iLocPair_all(iiLoc), iMetric_prob, :, :, :), 1, 6), 2, 1);
                        [pred_ave, ~, ~, pred_SEM] = getCI(getCI(metric_pred_allCond(iModelA, iModelB, iLocPair_all(iiLoc), iMetric_prob, :, :, :), 1, 6), 2, 1);
                        [R2_NOM_ave, ~, ~, R2_NOM_SEM] = getCI(getCI(R2_NOM_allCond(iModelA, iModelB, iLocPair_all(iiLoc), iMetric_prob, :, :), 1, 6), 2, 1);

                        % ---- Store R2 mean into table ----
                        rowIdx = find(iModelB_all_plot == iModelB);  % which row for this model
                        colIdx = iiLoc;                              % col 1 or 2 for the pair
                        R2_tab(rowIdx, colIdx) = R2_NOM_ave;        % store mean R2 (not SEM)

                        if iModelB==1, color_pred = colors_comb(iLocPair_all(iiLoc), :); lw = 3;
                        else, color_pred='k'; lw = 1.5;
                        end

                        if iPlotMode ~= 1
                            % Plot prediction
                            patch([IV_ave; flip(IV_ave)], [pred_ave-pred_SEM; flip(pred_ave+pred_SEM)], ones(1, 3) / 2, 'FaceAlpha', .3, 'linestyle', 'none', 'handlevisibility', 'off')
                            plot(IV_ave, pred_ave, 'lineStyle', lineStyle_all{iModelB}, 'color', colors_comb(iLocPair_all(iiLoc), :), 'linewidth', lw)
                        end

                        % Plot measurement
                        % errorbar(IV_ave, data_ave, IV_SEM, '.', 'horizontal', 'CapSize', 0, 'color', colors_comb(iLocPair_all(iiLoc), :), 'handlevisibility', 'off')
                        errorbar(IV_ave, data_ave, data_SEM, '.', 'vertical', 'CapSize', 0, 'color', colors_comb(iLocPair_all(iiLoc), :), 'handlevisibility', 'off', 'linewidth', lw/1.5)

                        % Plot averaged data of each bin
                        for iBin = 1:nBins
                            facecolor = 'w';
                            plot(IV_ave(iBin), data_ave(iBin), 'o', 'markeredgecolor', colors_comb(iLocPair_all(iiLoc), :), ...
                                'markerfacecolor', facecolor, 'MarkerSize', nTrials_ave(iBin) / szScaling + 5, 'LineWidth', lw/1.5, 'LineStyle', 'none', 'handlevisibility', 'off')
                        end
                    end % iiLoc

                    yline(.5, '--', 'linewidth', 2);

                    % xlim([,150])
                    ylim([0, 1])
                    yticks(0:.2:1)
                    ylabel(namesMetrics_prob_full{iMetric_prob})

                    switch iLocPair_all(1)
                        case 1, xlim([0, 200])
                        case 6, xlim([0, 100])
                        case 5, xlim([0, 80])
                    end
                    xlabel('Binned decision variable')
                    if isubj == 1, legend(namesModelB(iModelB_all_plot), 'Location', 'best'), end

                end % iModelB

                % ==== Print R2 as a 4x3 table (left column = model name) ====
                ax = gca;

                nRows_subj = nModelB_plot;
                nCols= 3;   % now: [ModelName | Loc1 | Loc2]

                % Create normalized coordinates for each cell
                x_pos = linspace(x_base, x_base + x_width, nCols + 1);
                x_cells = x_pos(1:nCols) + diff(x_pos(1:2))/2;   % center of each column

                y_pos = linspace(y_base, y_base + y_height, nRows_subj + 1);
                y_cells = fliplr(y_pos(1:nRows_subj)) + diff(y_pos(1:2))/2;

                if iPlotMode ~= 1
                    % ---- Column 1: model names (B1, B2, ...) ----
                    for r = 1:nRows_subj
                        name_str = sprintf('%s', namesModelB_plot{iModelB_all_plot(r)});

                        text(x_cells(1), y_cells(r), name_str, 'Units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize', sz_font, 'FontWeight', 'bold', 'interpreter', 'latex');
                        % text(x_cells(1), y_cells(r), '$R^2$=', 'Units','normalized', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', sz_font, 'FontWeight', 'bold', 'interpreter', 'latex');
                    end

                    % ---- Columns 2–3: R2 values ----
                    for r = 1:nRows_subj
                        for c = 1:2   % two locations
                            R2_val = R2_tab(r, c);
                            if isnan(R2_val), continue; end

                            str_cell = sprintf('%.0f%%', R2_val * 100);

                            % text(x_cells(c+1), y_cells(r), str_cell, 'Units','normalized', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', sz_font);
                        end
                    end
                end % if iPlotMode==2
                ax = gca;
                ax.FontSize=sz_axis;
                ax.LineWidth = 1.5; % Adjust the value as desired
            end % iMetric

            sgtitle(sprintf('n=%d [A%d] [L%d%d] [nIter=%d] %s', nsubj, iModelA, iLocPair_all, nIterxJob, namesMetrics_prob{iMetric_prob}))
            % set(findall(gcf, '-property', 'fontsize'), 'fontsize', 15)
            % set(findall(gcf, '-property', 'linewidth'), 'linewidth', 2)

            % Save the figure
            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_A%d_group_%s.png', nameFolder_Fig_NOM_metrics, nsubj, iLocPair_all, iModelA_plot, namesPlotMode{iPlotMode}))
            close(gcf)
            fprintf('DONE\n')
        end % iLocComb
    end % iModelA
end % iPlotMode

clear *_allCond
fprintf('\n\n Plotting DONE\n\n')

%% Plot metrics vs. IV for each idvd
clc, fprintf('\n\n Plotting STARTED......\n\n')
lineStyle_all = {'-', '-', '--', ':', '-.', '-', '--', ':', '-.'};

% iModelA_plot = 1;%iModelA_all % just plot pred of the model using RC-derived template (happens to be yhe best model)
% iModelA = iModelA_plot;
iModelB_all_plot = 1;

% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond', 'IV_allCond', 'nTrials_allCond', 'metric_data_allCond', 'metric_pred_allCond', 'R2_NOM_allCond')
% Define folder for saving figures for each model A and model B
nameFolder_Fig_NOM_metrics = sprintf('%s/A%d', nameFolder_Fig_NOM_Trialwise, iModelA_plot);
if isempty(dir(nameFolder_Fig_NOM_metrics)), mkdir(nameFolder_Fig_NOM_metrics), end

for iLocSingle = iLocComb_all
    fprintf('\nL%d...', iLocSingle)
    szScaling = 80; %if any(iLocComb==[6,7]), scalingF=100; elseif iLocComb==8, scalingF=200; else, scalingF=50; end

    %%%%%%%%%%%%
    for iMetric_prob = 1:nMetrics_prob
        figure('Position', [0 0 2e3 1.5e3])

        for isubj = 1:nsubj
            subplot(nRows_subj, nCols_subj, isubj), hold on

            subjName = subjList{isubj};

            % Create a string to store R2
            str_R2 = '|';

            % Create a string to store parameter estimates for display
            str_est = [];
            params_allCond(isnan(params_allCond))=0;
            for iModelB = iModelB_all_plot
                vals = strjoin(string(round(getCI(params_allCond(iModelA_plot, iModelB, iLocSingle, isubj, :, :), 1, 5)',2)), ", ");
                str_est = [str_est, sprintf('B%d: %s\n', iModelB, vals)];
            end

            for iModelB = iModelB_all_plot

                % Compute medians and CIs
                [IV_allBins, ~, ~, IV_allBins_neg, IV_allBins_pos] = getCI(IV_allCond(iModelA_plot, iModelB, iLocSingle, isubj, :, :), 1, 5);
                [nTrials_allBins, nTrials_allBins_lb, nTrials_allBins_ub] = getCI(nTrials_allCond(iModelA_plot, iModelB, iLocSingle, isubj, :, :), 1, 5);
                [data_allBins, ~, ~, data_allBins_neg, data_allBins_pos] = getCI(metric_data_allCond(iModelA_plot, iModelB, iLocSingle, iMetric_prob, isubj, :, :), 1, 6);
                [pred_allBins, pred_allBins_lb, pred_allBins_ub] = getCI(metric_pred_allCond(iModelA_plot, iModelB, iLocSingle, iMetric_prob, isubj, :, :), 1, 6);
                [params_allBins, params_allBins_lb, params_allBins_ub] = getCI(params_allCond(iModelA_plot, iModelB, iLocSingle, isubj, :, :), 1, 5);
                [R2_NOM_med, R2_NOM_lb, R2_NOM_ub] = getCI(R2_NOM_allCond(iModelA_plot, iModelB, iLocSingle, iMetric_prob, isubj, :), 1, 6, .95);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Plot prediction
                if iModelB==1, color_pred = colors_comb(iLocSingle, :); lw = 3;
                else, color_pred='k'; lw = 1;
                end
                patch([IV_allBins; flip(IV_allBins)], [pred_allBins_lb; flip(pred_allBins_ub)], ones(1, 3) / 2, 'FaceAlpha', .3, 'linestyle', 'none', 'handlevisibility', 'off')
                plot(IV_allBins, pred_allBins, 'lineStyle', lineStyle_all{iModelB}, 'color', color_pred, 'linewidth', lw)

                % Plot measurement
                errorbar(IV_allBins, data_allBins, IV_allBins_neg, IV_allBins_pos, '.', 'horizontal', 'CapSize', 0, 'color', colors_comb(iLocSingle, :), 'handlevisibility', 'off', 'linewidth', lw)
                errorbar(IV_allBins, data_allBins, data_allBins_neg, data_allBins_pos, '.', 'vertical', 'CapSize', 0, 'color', colors_comb(iLocSingle, :), 'handlevisibility', 'off', 'linewidth', lw)

                yline(.5, 'k--');

                % Plot averaged data of each bin
                for iBin = 1:nBins
                    % if isubj>nMarkersMax, facecolor=colors_comb(iLocComb, :); else, facecolor='w'; end
                    facecolor = 'w';
                    plot(IV_allBins(iBin), data_allBins(iBin), 'o', 'markeredgecolor', colors_comb(iLocSingle, :), ...
                        'markerfacecolor', facecolor, 'MarkerSize', nTrials_allBins(iBin) / szScaling + 5, 'LineWidth', 1, 'LineStyle', 'none', 'handlevisibility', 'off')
                end

                % xlim([,150])
                ylim([0, 1])
                ylabel(namesMetrics_prob{iMetric_prob})
                xlabel('Binned decision variable')
                if isubj == 1, legend(namesModelB(iModelB_all), 'Location', 'best'), end

                % R2 median and CI
                str_R2 = sprintf('%s B%d: %.2f [%.2f, %.2f]', str_R2, iModelB, R2_NOM_med, R2_NOM_lb, R2_NOM_ub);
            end % iModelB

            % print estimates
            text(100, .2, str_est, 'FontSize', 6)
            title(sprintf('[%s] R^2: %s', subjName, str_R2))
        end % isubj

        sgtitle(sprintf('n=%d [A%d] [L%d] [nIter=%d] %s', nsubj, iModelA_plot, iLocSingle, nIterxJob, namesMetrics_prob{iMetric_prob}))
        set(findall(gcf, '-property', 'fontsize'), 'fontsize', 12)

        % Save the figure
        drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d_A%d_%s.png', nameFolder_Fig_NOM_metrics, nsubj, iLocSingle, iModelA_plot, namesMetrics_prob{iMetric_prob}))
        close(gcf)
    end % iMetric

    fprintf('DONE\n')
end % iLocComb

% end % iModelA

fprintf('\n\n Plotting DONE\n\n')

%% Plot GoF for ModelA x ModelB x Loc (use basicFxn_drawBars_permutation)
clc
clc, fprintf('\n\n Plotting STARTED......\n\n')
% Load data
load(nameFolder_Data_SaveCompile, 'nLL_allCond')

% Define folder for saving GoF figures
nameFolder_Fig_NOM_GoF = fullfile(nameFolder_Fig_NOM_Trialwise, 'GoF');
if isempty(dir(nameFolder_Fig_NOM_GoF)), mkdir(nameFolder_Fig_NOM_GoF); end

% iModelB_selected = [1,3,2,4]; flag_plotIDVD = 0; flag_plotDiff = 0;
iModelB_selected = [1,2,3,4]; flag_plotIDVD = 0; flag_plotDiff = 0;
iModelB_selected = [1,3,4]; flag_plotIDVD = 1; flag_plotDiff = 1; %  {'FullModel'}  {'NoMultiN'}    {'NoPrivN'}
iModelB_selected = [1,2]; flag_plotIDVD = 1; flag_plotDiff = 1; %  {'FullModel'}  {'NoMultiN'}

% Goodness-of-fit (GoF) measures to plot
namesGoF = {'nLL', 'AIC-nLL', 'AICc-nLL', 'BIC-nLL'}; % no raw nLL
nGoFs = numel(namesGoF);

% Common y-limit for all GoF measures (already positive deltas)
y_lim = [0, 16];
sz_label = 20;
% Loop over GoF measures
for iGoF = 1% 1:nGoFs % just plot nLL

    for iModelA = 1%iModelA_all

        for iLocSingle = 1:nLocComb8

            sz_wd_perBar = 70;
            nBars = numel(iModelB_all(iModelB_selected));
            sz_fig = [nBars*sz_wd_perBar, 250];

            % data: [ModelA x ModelB x LocComb x Metric x Subj x Iteration]
            nLL_allIter_allSubj = nLL_allCond(iModelA, iModelB_selected, iLocSingle, :, :);
            dnLL_allIter_allSubj = nan(size(nLL_allIter_allSubj));

            % For each subject and each iteration, subtract the minimum across models
            for isubj = 1:nsubj
                parfor iIter = 1:nIterxJob
                    d = nLL_allIter_allSubj(:, :, :, isubj, iIter);
                    d_min = min(d(:));
                    dnLL_allIter_allSubj(:, :, :, isubj, iIter) = d - d_min;
                end
            end

            dnLL_allIter_allSubj = squeeze(dnLL_allIter_allSubj);

            str_title = sprintf('n=%d nIter=%d %s [L%d] ModelB [%s]', nsubj, nIterxJob, namesGoF{iGoF}, iLocSingle, strjoin(string(iModelB_selected), ' '));
            x_ticks = nan;
            ref = nan;

            %------------------------------%
            basicFxn_drawBars_permutation(dnLL_allIter_allSubj, ref, repmat(colors_comb(iLocSingle, :), nBars, 1), namesModelB(iModelB_all(iModelB_selected)), x_ticks, x_ticks, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj)
            %------------------------------%
            % ylim(y_lim)
            ylabel(sprintf('\\Delta nLL'), 'FontSize', sz_label);

            drawnow, pause(1), saveas(gcf, fullfile(nameFolder_Fig_NOM_GoF, sprintf('n%d_%s_B%s_L%d.png', nsubj, namesGoF{iGoF}, strjoin(string(iModelB_selected), ''), iLocSingle)));
            close(gcf)

        end % iLocSingle
    end % iModelA
end % iGoF

clear nLL_allCond
fprintf('\n\n Plotting DONE\n\n')



%% Compare NOM parameters across locations (use basicFxn)
clc, fprintf('\n\n Plotting STARTED......\n\n')
% Load data
load(nameFolder_Data_SaveCompile, 'params_allCond')

% Define folder for saving comparison figures
nameFolder_Fig_NOM_compParams = fullfile(nameFolder_Fig_NOM_Trialwise, 'CompNOMparams');
if isempty(dir(nameFolder_Fig_NOM_compParams))
    mkdir(nameFolder_Fig_NOM_compParams);
end

% nParamsMax = 3; % max three params
sz_wd_perBar = 120;
nBars = 2;
sz_fig = [nBars*sz_wd_perBar, 250];
flag_plotIDVD = 1;
flag_plotDiff = 1;

for iModelB = iModelB_plot %iModelB_all % just plot the best model

    % Number of parameters for this Model B
    nParamsB = numel(namesModelBparams{iModelB});

    for iGroup = 1:nGroups % e.g., {[1, 8]} or more pairs if desired
        iLocPair_all = iLocGroups_all{iGroup};

        for iParam = 1:nParamsB

            NOMp_allIter_allSubj = squeeze(params_allCond(iModelA_plot, iModelB, iLocPair_all, :, :, iParam));

            x_ticks = namesLocComb(iLocPair_all);
            switch iModelB
                case 1
                    switch iParam
                        case 1, y_ticks = linspace(NOMp1_lb, .8, 5);
                        case 2, y_ticks = linspace(NOMp2_lb, 40, 5);
                        case 3, y_ticks = linspace(NOMp3_lb, 40, 5);
                    end
            end
            y_ticklabels = round(y_ticks, 2);

            str_title = sprintf('n=%d nIter=%d L%d%d [A%dB%d] %s', nsubj, nIterxJob, iLocPair_all, iModelA_plot, iModelB_plot, namesModelBparams{iModelB}{iParam});
            %------------------------------%
            basicFxn_drawBars_permutation(NOMp_allIter_allSubj, ref, colors_comb(iLocPair_all, :), x_ticks, y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIterxJob, markers_allSubj)
            %------------------------------%
            ylabel(sprintf('%s', namesModelBparams{iModelB}{iParam}))
            % Save the figure
            drawnow, pause(1), saveas(gcf, sprintf('%s/n%d_L%d%d_A%dB%d_NOMp%d.png', nameFolder_Fig_NOM_compParams, nsubj, iLocPair_all, iModelA_plot, iModelB_plot, iParam))
            close(gcf);
        end % iParam
    end % iGroup
end % iModelB
clear params_allCond
fprintf('\n\n Plotting DONE\n\n')
