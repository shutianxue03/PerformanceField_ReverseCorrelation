function flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj, asymY_med_allSubj, x_ticks, y_ticks, x_ticklabels, y_ticklabels, text_title, markers_allSubj)

% to plot and compare asymmetries between two measurements
% Inputs:
%    asym1_med_allSubj: nsubj x 1, median of boostrapping, to plot on the x-axis
%    asymY_med_allSubj: nsubj x 1, median of boostrapping, to plot on the y-axis
%    x_ticks: vector, containing 5 values
%    y_ticks: vector, containing 5 values
%    x_ticklabels: vector, containing 5 values
%    y_ticklabels: vector, containing 5 values
%    text_title

%% figure setting
wd_ref = 3; % line width of the reference line
wd_border = 4; % default =4
fsz_ticks = 45; % TunC: 45; pA; NOM
sz_marker = 40;
nMarkerMax = 11;
% markers_allSubj = {'o', 's', 'd', '^','v',  '<', '+','p', 'h', 'x', '>',     'o', 's', 'd', '^'}; % for each subjclc

%% extract nsubj and nLoc and make assertion
nsubj = length(asymX_med_allSubj);
assert(length(markers_allSubj) >= nsubj)

%%
figure('Position', [0 200 1e3 1e3]); hold on, box on

%% idvd data
for isubj = 1:nsubj
%     if isubj<=nMarkerMax, 
        faceColor = 'w'; edgeColor = 'k';
%     else, 
        faceColor = 'k'; edgeColor = 'w';
%     end
    plot(asymX_med_allSubj(isubj), asymY_med_allSubj(isubj),  markers_allSubj{isubj}, ...
        'markerfacecolor', faceColor, 'markeredgecolor',edgeColor, 'markersize', sz_marker, 'linewidth', wd_border)
end

%% get corr of different types
[r2, p_r2] = corr(asymX_med_allSubj, asymY_med_allSubj); % left: assume slope <0
[rL, p_rL] = corr(asymX_med_allSubj, asymY_med_allSubj, 'Tail','left'); % left: assume slope <0
[rR, p_rR] = corr(asymX_med_allSubj, asymY_med_allSubj, 'Tail','right'); % left: assume slope <0

[rho2, p_rho2] = corr(asymX_med_allSubj, asymY_med_allSubj, 'Type', 'spearman'); % left: assume slope <0
[rhoL, p_rhoL] = corr(asymX_med_allSubj, asymY_med_allSubj, 'Type', 'spearman', 'Tail','left'); % left: assume slope <0
[rhoR, p_rhoR] = corr(asymX_med_allSubj, asymY_med_allSubj, 'Type', 'spearman', 'Tail','right'); % right: assume slope>0

[tau2, p_tau2] = corr(asymX_med_allSubj, asymY_med_allSubj, 'Type', 'Kendall'); % left: assume slope <0
[tauL, p_tauL] = corr(asymX_med_allSubj, asymY_med_allSubj, 'Type', 'Kendall', 'Tail','left'); % left: assume slope <0
[tauR, p_tauR] = corr(asymX_med_allSubj, asymY_med_allSubj,'Type', 'Kendall', 'Tail','right'); % right: assume slope>0

text_corr_all = sprintf('r=%.2f (%.3f) | rho=%.2f (%.3f)  | tau=%.2f (%.3f)\n', ...
    r2, p_r2,    rho2, p_rho2,    tau2, p_tau2);

flag_sig=''; 
% if p_tR<.05, flag_sig='_sig'; elseif p_tR<.1, flag_sig = '_mg'; end

%% linear regression
% if any([p_r2/2, p_rho2/2, p_tau2/2]<=.1)
% % if p_tR<.1
%     lm = polyfit(asymX_med_allSubj, asymY_med_allSubj, 1);
%     x_lm2 = linspace(min(asymX_med_allSubj), max(asymX_med_allSubj), 2);
%     yfit = polyval(lm, x_lm2);
%     plot(x_lm2, yfit,'-', 'color', ones(1,3)*.4, 'handlevisibility', 'off', 'linewidth', wd_border * 1.5);
%     % eta2 = var(polyval(lm, HVA_med_allSubj))/var(VMA_med_allSubj);
% end

%% ref (i.e., at 0)
xline(0, 'color', ones(1,3)/2, 'linewidth', wd_ref);
yline(0, 'color', ones(1,3)/2, 'linewidth', wd_ref);

%% ticks and limits
if ~isnan(x_ticks), xticks(x_ticks), xlim(x_ticks([1, end])), end
if ~isnan(y_ticks), yticks(y_ticks), ylim(y_ticks([1, end])), end
if ~isnan(x_ticklabels), xticklabels(x_ticklabels),  end
if ~isnan(y_ticklabels), yticklabels(y_ticklabels),  end

%% Print correlation coefficient

%% figure format
axis square
ax = gca;
ax.XAxis.FontSize = fsz_ticks;
ax.YAxis.FontSize = fsz_ticks;
ax.LineWidth = wd_border;

%% title
title(sprintf('%s\n%s', text_title, text_corr_all))



