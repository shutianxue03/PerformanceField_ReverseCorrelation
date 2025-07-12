
function p = fxn_drawCorr_Asym(x_med_allSubj, y_med_allSubj, xticks_, yticks_, namesType, namesLocComb, iLocComb_all, flag_plotCI, markers_allSubj, itype)

getAsym = @(a,b) (a-b)./(a+b);
nsubj = size(x_med_allSubj, 1);

%% figure setting
wd_border = 6;
sz_title = 20;
sz_ticks = 50; % 30
sz_marker = 40;


%     x_med = getCI(x_med_allSubj, 1, 2);
%     y_med = getCI(y_med_allSubj(:, :, :, ip), 1, 2);
    asym_x_med_allSubj = getAsym(x_med_allSubj(:, 1), x_med_allSubj(:, 2));
    asym_y_med_allSubj = getAsym(y_med_allSubj(:, 1), y_med_allSubj(:, 2));


%% plot

for isubj = 1:nsubj
    plot(asym_x_med_allSubj(isubj), asym_y_med_allSubj(isubj),  markers_allSubj{isubj}, ...
        'markerfacecolor', 'w', 'markeredgecolor', 'k', 'markersize', sz_marker, 'linewidth', wd_border)
end

%%  linear regression
lm = polyfit(asym_x_med_allSubj, asym_y_med_allSubj, 1);
x_lm2 = linspace(min(asym_x_med_allSubj), max(asym_x_med_allSubj), 2);
yfit = polyval(lm, x_lm2);
plot(x_lm2, yfit, '', 'color', ones(1,3)*.4, 'handlevisibility', 'off', 'linewidth', wd_border*1.5);

%% extra lines
% if (ifeature == 1) && (iparam == 3), yline(0, 'color', ones(1,3)*.5); end

%% get corr
[r2, p_p2] = corr(asym_x_med_allSubj, asym_y_med_allSubj); % left: assume slope <0
[rL, p_pL] = corr(asym_x_med_allSubj, asym_y_med_allSubj, 'Tail','left'); % left: assume slope <0
[rR, p_pR] = corr(asym_x_med_allSubj, asym_y_med_allSubj, 'Tail','right'); % left: assume slope <0

[rho2, p_r2] = corr(asym_x_med_allSubj, asym_y_med_allSubj); % left: assume slope <0
[rhoL, p_rL] = corr(asym_x_med_allSubj, asym_y_med_allSubj, 'Type', 'spearman', 'Tail','left'); % left: assume slope <0
[rhoR, p_rR] = corr(asym_x_med_allSubj, asym_y_med_allSubj, 'Type', 'spearman', 'Tail','right'); % left: assume slope <0

[tau2, p_t2] = corr(asym_x_med_allSubj, asym_y_med_allSubj); % left: assume slope <0
[tauL, p_tL] = corr(asym_x_med_allSubj, asym_y_med_allSubj, 'Type', 'Kendall', 'Tail','left'); % left: assume slope <0
[tauR, p_tR] = corr(asym_x_med_allSubj, asym_y_med_allSubj,'Type', 'Kendall', 'Tail','right'); % left: assume slope <0

tail=["2-tail"; "R-tail"; "L-tail"];

Pearson=round([r2;rL;rR],2);
p1=round([p_p2;p_pL;p_pR],3);
Spearman=round([rho2;rhoL;rhoR],2);
p2=round([p_r2;p_rL;p_rR],3);
Kendall=round([tau2;tauL;tauR],2);
p3=round([p_t2;p_tL;p_tR],3);

table(tail, Pearson, p1, Spearman, p2, Kendall, p3)

% report Tau
r = tauR; 
p = p_tR;
text_corr = sprintf('r = %.2f (p = %.3f) eta^2=%.2f', ...
    round(r,2), round(p,3), round(var(polyval(lm, asym_x_med_allSubj))/var(asym_y_med_allSubj), 2));

%% xticks
% xlim_ = xticks_([1, end]);
% ylim_ = yticks_([1, end]);
% xlim(xlim_), xticks(xticks_) % temporaily hard-coded, since x is usually CS
% ylim(ylim_), yticks(yticks_)
xlabel([])
axis square
% title(sprintf('%s vs. CS\nPartial r = %.2f (p = %.3f)', namesYaxis{ip}, round(r,2), round(p,3)))

ax = gca;
ax.XAxis.FontSize = sz_ticks;
ax.YAxis.FontSize = sz_ticks;
ax.LineWidth = wd_border;

yline(0, 'color', ones(1,3)*.5, 'linewidth', wd_border);
xline(0, 'color', ones(1,3)*.5, 'linewidth', wd_border);

axis square

title(sprintf('[%s] %s vs. %s %s', ...
    namesType{itype}, namesLocComb{iLocComb_all(1)}, namesLocComb{iLocComb_all(2)}, text_corr), 'fontsize', sz_title)

