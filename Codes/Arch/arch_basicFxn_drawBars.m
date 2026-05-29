
function flag_sig = basicFxn_drawBars(med_allSubj, ref, colors, x_ticks, y_ticks, y_ticklabels, flag_plotIDVD, flag_plotDiff, text_title, flag_pairwiseComp, sz_fig, nPairs, str_tail)

% Inputs:
%    med_allSubj: nsubj x nBars, median of boostrapping
%    ref: one reference line, ylines
%    colors: nBars x 3, each rows indicates a bar
%    x_ticks: cell, containing strings indicating location names
%    y_ticks: vector, containing 3 values
%    y_ticklabels: vector, containing 3 value
%    flag_plotIdvd: 1=plot idvd data; 0=NOT
%    flag_plotDiff: ~isnan=plot SEM of difference
%    text_title: string
%    flag_pairwiseComp: 1=print pairwise comparison in the title
%    sz_fig

%% define sizes
sz_marker_idvd = 10;
sz_marker_ave = 30;
fsz_ticks = 15; % tunC: 20; BEHAV: xx; NOM: 15
fsz_title = 10; % font size of titlte
wd = 2; % line width of axis
wd_bar = .5; % the width of the bar (not the bar edge!!)
wd_ref = wd; % line width of the reference line
% markers_allSubj = {'o', 's', 'd', '^','v',  '<', '+','p', 'h', 'x', '>',      'o', 's', 'd', '^'}; % for each subj
interval_diffBar = 8; % higher, closer the comparison bar is to the top of the figure
% markers_allSubj = {'o', 'o', 'o', 'o','o',  'o', 'o','o', 'o', 'o', 'o',      'o', 'o', 'o'};
nsubj_max = 11;
fprintf('\n\n *** %s ***\n', text_title)

%% extract nsubj and nBars and make assertion
[nsubj, nBars] = size(med_allSubj);
% assert(length(markers_allSubj) >= nsubj)
assert(length(x_ticks) == nBars)

%% get ave/SEM
if nsubj>100
    [ave, ~, ~, sem] = getCI(med_allSubj, 1, 1);
else
    [ave, ~, ~, sem] = getCI(med_allSubj, 2, 1);
end

figure('Position', [0, 200, sz_fig])
hold on

%% plot bars/disks and errorbars
for iBar = 1:nBars
    % normal
    %     if flag_plotIdvd
    %         plot(iBar, ave(iBar),'o', 'MarkerFaceColor', colors(iBar, :),  'MarkerEdgeColor', 'w', 'MarkerSize', sz_marker_ave, 'linewidth', wd, 'HandleVisibility', 'off')
    %     else
    bar(iBar, ave(iBar), 'FaceColor', colors(iBar, :),  'EdgeColor', colors(iBar, :), 'barwidth', wd_bar, 'HandleVisibility', 'off')
    %     end
    errorbar(iBar, ave(iBar), sem(iBar), '.', 'color', colors(iBar, :), 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off')

    %     if ~flag_plotIdvd
    if ave(iBar)>0
        errorbar(iBar, ave(iBar), sem(iBar), 0, '.w', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off')
    else
        errorbar(iBar, ave(iBar), 0, sem(iBar), '.w', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off')
    end
    %     end

    % if plot HVA vs. VMA
    %     bar(iBar, ave(iBar), 'FaceColor', 'w',  'EdgeColor', 'k', 'barwidth', wd_bar, 'linewidth', wd)
    %     errorbar(iBar, ave(iBar), sem(iBar), '.', 'color', 'k', 'CapSize', 0, 'linewidth', wd)
end % end of iiLoc

%% idvd data
buffer = .2;
if flag_plotIDVD

    if nBars == 2 % for nBars=2, plot idvd data between two bars
        x = [1+buffer, 2-buffer];
    else % for nBars>2, plot idvd data at the center of each bar
        x = 1:nBars;
    end

    for isubj = 1:nsubj
        if isubj<=nsubj_max, c = 'w'; else, c= ones(1,3)/2; end
        if nBars==2

            %             plot(x, med_allSubj(isubj, :), [markers_allSubj{isubj}, '-'],  ...
            plot(x, med_allSubj(isubj, :), '-',  ...
                'color',ones(1,3)*.7, 'markerfacecolor', 'w', 'markeredgecolor', ones(1,3)*.7, ...
                'markersize', sz_marker_idvd, 'linewidth', wd)
        else
            %                     plot([1,2], med_allSubj(isubj, [1,2]), [markers_allSubj{isubj}, '-'],  ...
            %             'color',ones(1,3)*.7, 'markerfacecolor', 'w', 'markeredgecolor', ones(1,3)*.7, ...
            %             'markersize', sz_marker_idvd, 'linewidth', wd_bar, 'linewidth', 2)
            %
            %             plot([1,3], med_allSubj(isubj, [1,3]), [markers_allSubj{isubj}, '-'],  ...
            plot(1:nBars+buffer, med_allSubj(isubj, :), '-',  ...
                'color',ones(1,3)*.7, 'markerfacecolor', 'w', 'markeredgecolor', ones(1,3)*.7, ...
                'markersize', sz_marker_idvd, 'linewidth', wd)
        end
    end
end

%% draw ref and compare bars with ref (print in command window)
if ~isnan(ref)
    if length(ref)>1; error('ALERT: there are more than one REF!!!'), end
    yline(ref, 'color', ones(1,3)/2, 'linewidth', wd_ref, 'HandleVisibility', 'off');
    for iBar = 1:nBars
        [~, p, ~, stats] = ttest(med_allSubj(:, iBar)-ref);
        cohenD = fxn_getES(med_allSubj(:, iBar), ref);
        fprintf('   %s (%.2f) vs. ref (%.1f): t=%.2f, p=%.3f, d=%.2f (%d/%d)\n', ...
            x_ticks{iBar}, mean(med_allSubj(:, iBar)), ref, stats.tstat, p, cohenD, sum(med_allSubj(:, iBar)> ref), nsubj);
    end
end

%% conduct ANOVA (regardless of nBars)
% indBar = repmat(1:nBars, nsubj, 1);
% text_ANOVA = print_nANOVA({'Loc'}, med_allSubj(:), {indBar(:)}, nsubj, 1);

%% if nBars>2, compare every pair
text_testPairs = '';
if nBars>2 && flag_pairwiseComp
    indPairs = nchoosek(1:nBars, 2); % all possible pairs across columns
    %     npairs = size(indPairs , 1);
    for iPair = 1:nPairs
        x1 = med_allSubj(:, indPairs(iPair, 1));
        x2 = med_allSubj(:, indPairs(iPair, 2));
        [~, p,~, stats] = ttest(x1, x2, 'tail', str_tail);
        p=p*nPairs; %
        cohenD = fxn_getES(x1, x2);
        text_testPairs = [text_testPairs, sprintf('%s vs. %s: t=%.2f, p=%.3f, d=%.2f (%d/%d)\n', ...
            x_ticks{indPairs(iPair, 1)}, x_ticks{indPairs(iPair, 2)}, stats.tstat, p, cohenD, sum(med_allSubj(:, indPairs(iPair, 1))> med_allSubj(:, indPairs(iPair, 2))), nsubj)];
    end
end

%% compare two loc (ttest, draw sem of diff)
if nBars==2
    x1 = med_allSubj(:, 1);
    x2 = med_allSubj(:, 2);
    [~, p,~, stats] = ttest(x1, x2);
    cohenD = fxn_getES(x1, x2);
    flag_sig=''; if p<.05, flag_sig='_sig'; elseif p<.1, flag_sig='_mg'; end
    text_testPairs = sprintf('%s vs. %s: t(%d)=%.2f, p=%.3f, d=%.2f (%d/%d)\n', ...
        x_ticks{1}, x_ticks{2},stats.df, stats.tstat, p, cohenD, sum(med_allSubj(:, 1)> med_allSubj(:, 2)), nsubj);

    if flag_plotDiff
        yDiffSEM = y_ticks(end) - (y_ticks(end) - y_ticks(1))/interval_diffBar;
        diff_sem = std(x1 - x2)/sqrt(nsubj);
        errorbar(1.5, yDiffSEM, diff_sem, 'k', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off')
        plot([1,2], [yDiffSEM, yDiffSEM], 'k-', 'linewidth', wd, 'HandleVisibility', 'off')
        %     string_s = getString_starts(p);
    end
end

%% ticks, limits and labels
if ~isnan(y_ticks), yticks(y_ticks), ylim(y_ticks([1, end])), end
if ~isnan(y_ticklabels), yticklabels(y_ticklabels), end
xticks(1:nBars), xticklabels([]) % manually add tick labels /symbols on the slide
% if ~isnan(x_label), xlabel(x_label), end
% if ~isnan(y_label), ylabel(y_label), end

% if plot idvd data, leave more space for the middle
if flag_plotIDVD, buffer = .6; else, buffer = .5; end
xlim([1-buffer, nBars+buffer])

%% size
ax = gca;
ax.XAxis.FontSize = fsz_ticks;
ax.YAxis.FontSize = fsz_ticks;
ax.LineWidth = wd;

%% title
% title(sprintf('%s\n%s%s', text_title, text_ANOVA, text_testPairs), 'fontsize', fsz_title)

title(sprintf('%s\n%s', text_title, text_testPairs), 'fontsize', fsz_title)


