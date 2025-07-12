
function p = fxn_drawCorr_EEI(ip, x, y, xticks_, yticks_, namesType, namesLocComb, iLocComb_all, flag_plotCI, markers_allSubj, itype)

getEEI = @(a,b) (a-b)./(a+b);
% getEEI = @(a,b) (a-b);
% getEEI = @(a,b) (a./b)-1;
nsubj = size(x, 1);

%% figure setting
wd_border = 6;
sz_title = 20;
sz_ticks = 50; % 30
sz_marker = 40;

%% get EEI
if flag_plotCI
    % way 1: get EEI of all boostraps, then take median
    EEI_x_allSubj = nan(nsubj, nB);
    EEI_y_allSubj = EEI_x_allSubj;
    for isubj = 1:nsubj
        x_ = squeeze(x(isubj, :, :)); if size(x_, 1) ~= nB, x_ = x_'; end
        y_ = squeeze(y(isubj, :, :)); if size(y_, 1) ~= nB, y_ = y_'; end
        % to ensure all values are non-negative (no need for x, since CST and pA are always positive)
        y_shifted = y_;% - min(y_(:)); % some tunC could be negative (baseline)
        x_shifted = x_;
        EEI_x = getEEI(x_shifted(:, 1), x_shifted(:, 2));
        EEI_y = getEEI(y_shifted(:, 1), y_shifted(:, 2));
        EEI_x_allSubj(isubj,:) = EEI_x;
        EEI_y_allSubj(isubj,:) = EEI_y;
    end
    
    % get median
    [EEI_x_med_allSubj, ~, ~, EEI_x_neg_allSubj, EEI_x_pos_allSubj] = getCI(EEI_x_allSubj, 1, 2);
    [EEI_y_med_allSubj, ~, ~, EEI_y_neg_allSubj, EEI_y_pos_allSubj] = getCI(EEI_y_allSubj, 1, 2);
else
    % way 2: get EEI of the median of all bootstraps
    % EEI_x_med_allSubj = getEEI(squeeze(median(cs_allSubj_(:, :, 1), 2)), squeeze(median(cs_allSubj_(:, :, 2), 2)));
    % EEI_y_med_allSubj = getEEI(squeeze(NOM_params_med_allSubj(:, 1, ip)), squeeze(NOM_params_med_allSubj(:, 2, ip)));
    
    x_med = getCI(x, 1, 2);
    y_med = getCI(y(:, :, :, ip), 1, 2);
    EEI_x_med_allSubj = getEEI(x_med(:, 1), x_med(:, 2));
    EEI_y_med_allSubj = getEEI(y_med(:, 1), y_med(:, 2));
end

%% plot
if flag_plotCI
    errorbar(EEI_x_med_allSubj, EEI_y_med_allSubj, EEI_y_neg_allSubj, EEI_y_pos_allSubj, EEI_x_neg_allSubj, EEI_x_pos_allSubj, 'k.', 'CapSize', 0, 'linewidth', wd_border/2)
end

for isubj = 1:nsubj
    plot(EEI_x_med_allSubj(isubj), EEI_y_med_allSubj(isubj),  markers_allSubj{isubj}, ...
        'markerfacecolor', 'w', 'markeredgecolor', 'k', 'markersize', sz_marker, 'linewidth', wd_border)
end

%%  linear regression
lm = polyfit(EEI_x_med_allSubj, EEI_y_med_allSubj, 1);
x_lm2 = linspace(min(EEI_x_med_allSubj), max(EEI_x_med_allSubj), 2);
yfit = polyval(lm, x_lm2);
plot(x_lm2, yfit, '', 'color', ones(1,3)*.4, 'handlevisibility', 'off', 'linewidth', wd_border*1.5);

%% extra lines
% if (ifeature == 1) && (iparam == 3), yline(0, 'color', ones(1,3)*.5); end

%% get corr
[r2, p_p2] = corr(EEI_x_med_allSubj, EEI_y_med_allSubj); % left: assume slope <0
[rL, p_pL] = corr(EEI_x_med_allSubj, EEI_y_med_allSubj, 'Tail','left'); % left: assume slope <0
[rR, p_pR] = corr(EEI_x_med_allSubj, EEI_y_med_allSubj, 'Tail','right'); % left: assume slope <0

[rho2, p_r2] = corr(EEI_x_med_allSubj, EEI_y_med_allSubj); % left: assume slope <0
[rhoL, p_rL] = corr(EEI_x_med_allSubj, EEI_y_med_allSubj, 'Type', 'spearman', 'Tail','left'); % left: assume slope <0
[rhoR, p_rR] = corr(EEI_x_med_allSubj, EEI_y_med_allSubj, 'Type', 'spearman', 'Tail','right'); % left: assume slope <0

[tau2, p_t2] = corr(EEI_x_med_allSubj, EEI_y_med_allSubj); % left: assume slope <0
[tauL, p_tL] = corr(EEI_x_med_allSubj, EEI_y_med_allSubj, 'Type', 'Kendall', 'Tail','left'); % left: assume slope <0
[tauR, p_tR] = corr(EEI_x_med_allSubj, EEI_y_med_allSubj,'Type', 'Kendall', 'Tail','right'); % left: assume slope <0

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
    round(r,2), round(p,3), round(var(polyval(lm, EEI_x_med_allSubj))/var(EEI_y_med_allSubj), 2));

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

