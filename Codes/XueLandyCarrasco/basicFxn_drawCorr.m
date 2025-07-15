function  flag_sig = basicFxn_drawCorr(x_med_allSubj, y_med_allSubj, colors, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, type_corr, type_tail, text_title, markers_allSubj)

% flag_sig = basicFxn_drawCorr(x_med_allSubj, y_med_allSubj, colors, x_ticks, y_ticks, x_ticklabels, y_ticklabels, flag_zeroMean, type_corr, type_tail, text_title, markers_allSubj)
% to plot and compare two values across locations for (1) performance and (2) tuning characteristics
% and assess partial correlation while controlling for location
% Inputs:
%    x_med_allSubj: nsubj x nLoc, median of boostrapping, to plot on the x-axis
%    y_med_allSubj: nsubj x nLoc, median of boostrapping, to plot on the y-axis
%    colors: nLoc x 3, each rows indicates a location
%    x_ticks: vector, containing 5 values
%    y_ticks: vector, containing 5 values
%    x_ticklabels: vector, containing 5 values
%    y_ticklabels: vector, containing 5 values
%    flag_zeroMean: 1=subtract mean; 0=NOT
%    typeCorr: string, 'pearson'=Pearson's r; 'spearman'=spearman's rho, 'kendall'=Kendall's tau
%    type_tail: 'both' (two tail), 'left', 'right'
%    text_title
%    markers_allSubj

%% figure setting
wd_border = 4; % 4
sz_ticks = 40; % 30
sz_marker = 30;
nMarkerMax = 11;

%% extract nsubj and nLoc and make assertion
[nsubj, nLoc] = size(x_med_allSubj);
assert(length(markers_allSubj) >= nsubj)

%% zero-mean
switch flag_zeroMean
    case 1 % control for VAR 1
        x_med_0mean_allSubj = x_med_allSubj - mean(x_med_allSubj, 1);
        y_med_0mean_allSubj = y_med_allSubj - mean(y_med_allSubj, 1);
    case 2 % control for VAR 2
        x_med_0mean_allSubj = x_med_allSubj - mean(x_med_allSubj, 2);
        y_med_0mean_allSubj = y_med_allSubj - mean(y_med_allSubj, 2);
    case 0 % no zero-mean
        x_med_0mean_allSubj = x_med_allSubj;
        y_med_0mean_allSubj = y_med_allSubj;
end

% %% get ave and SEM of x and y
% [x_ave, ~, ~, x_sem] = getCI(x_med_0mean_allSubj, 2, 1);
% [y_ave, ~, ~, y_sem] = getCI(y_med_0mean_allSubj, 2, 1);

%%
figure('Position', [0 200 1e3 1e3]); hold on, box on

%% idvd data
for iLoc = 1:nLoc
    for isubj = 1:nsubj
        %         if isubj<=nMarkerMax,
        faceColor = 'w';
        %else, faceColor = colors(iLoc, :); end
        plot(x_med_0mean_allSubj(isubj, iLoc), y_med_0mean_allSubj(isubj, iLoc),  markers_allSubj{isubj}, 'markerfacecolor',faceColor, ...
            'markeredgecolor', colors(iLoc, :), 'markersize', sz_marker, 'linewidth', wd_border)
    end % isubj
end % iLoc


%% get partial corr
ANOVA_indLoc = repmat(1:nLoc, nsubj, 1);
ANOVA_indSubj = repmat((1:nsubj)', 1, nLoc);
switch flag_zeroMean
    case 0
        [r_partial, p_partial] = corr(x_med_allSubj(:), y_med_allSubj(:), 'type', type_corr, 'tail', type_tail);
    case 1
        [r_partial, p_partial] = partialcorr(x_med_allSubj(:), y_med_allSubj(:), ANOVA_indLoc(:), ...
            'type', type_corr, 'tail', type_tail);
    case 2
        [r_partial, p_partial] = partialcorr(x_med_allSubj(:), y_med_allSubj(:), ANOVA_indSubj(:), ...
            'type', type_corr, 'tail', type_tail);
end
flag_sig=''; if p_partial<.05, flag_sig='_sig'; elseif p_partial<.1, flag_sig = '_mg'; end

%% corr for each loc
text_corr_perL = [];
for iLoc = 1:nLoc
    [r, p] = corr(x_med_allSubj(:, iLoc), y_med_allSubj(:, iLoc));
    text_corr_perL = [text_corr_perL, sprintf('L%d: r=%.2f, p=%.3f\n', iLoc, r, p)];
end

%% linear regression
lm = polyfit(x_med_0mean_allSubj(:), y_med_0mean_allSubj(:), 1);
x_lm2 = linspace(min(x_med_0mean_allSubj(:)), max(x_med_0mean_allSubj(:)), 2);
yfit = polyval(lm, x_lm2);
eta2 = var(polyval(lm, x_med_0mean_allSubj(:)))/var(y_med_0mean_allSubj(:));
if p_partial<.1
    plot(x_lm2, yfit,'-', 'color', ones(1,3)*.4, 'handlevisibility', 'off', 'linewidth', wd_border * 1.5);
end

%% ticks and limits
if ~isnan(x_ticks), xticks(x_ticks), xlim(x_ticks([1, end])), end
if ~isnan(y_ticks), yticks(y_ticks), ylim(y_ticks([1, end])), end
if ~isnan(x_ticklabels), xticks(x_ticklabels),  end
if ~isnan(y_ticklabels), yticks(y_ticklabels),  end

%% figure format
axis square
ax = gca;
ax.XAxis.FontSize = sz_ticks;
ax.YAxis.FontSize = sz_ticks;
ax.LineWidth = wd_border;

%% title
title(sprintf('%s\nPartial r = %.2f (p = %.3f) eta^2=%.2f\n%s', ...
    text_title, r_partial, p_partial, eta2, text_corr_perL))



