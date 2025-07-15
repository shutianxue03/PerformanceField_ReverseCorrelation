
wd_border = 3;
sz_title = 20;
sz_label = 40; %40
sz_ticks = 40; % 30
sz_marker = 22;

getEEI = @(a,b) (a-b)./(a+b);

switch itype, case 1, itype_pA = 7; case 2, itype_pA = 8; case 3, itype_pA = 6; end
itype_pA = 6;

x = cs_allSubj_; % cs (1/cst)

y = squeeze(metrics_allSubj(:, :,:, itype_pA));

figure('Position', [0 0 1200 1000]);
hold on, box on

% get EEI at each iteEEIn of boot
EEI_x_allSubj = nan(nsubj, nB); EEI_y_allSubj = EEI_x_allSubj;
for isubj = 1:nsubj
    x_ = squeeze(x(isubj, :, :));
    y_ = squeeze(y(isubj, :, :));
    % to ensure all values are non-negative (no need for x, since CST and pA are always positive)
    y_shifted = y_ - min(y_(:)); % some tunC could be negative (baseline)
    x_shifted = x_;
    EEI_x = getEEI(x_shifted(:, 1), x_shifted(:, 2));
    EEI_y = getEEI(y_shifted(:, 1), y_shifted(:, 2));
    EEI_x_allSubj(isubj,:) = EEI_x;
    EEI_y_allSubj(isubj,:) = EEI_y;
end

[EEI_x_med_allSubj, ~, ~, EEI_x_neg_allSubj, EEI_x_pos_allSubj] = getCI(EEI_x_allSubj, 1, 2);
[EEI_y_med_allSubj, ~, ~, EEI_y_neg_allSubj, EEI_y_pos_allSubj] = getCI(EEI_y_allSubj, 1, 2);

errorbar(EEI_x_med_allSubj, EEI_y_med_allSubj, EEI_y_neg_allSubj, EEI_y_pos_allSubj, EEI_x_neg_allSubj, EEI_x_pos_allSubj, 'k.', 'CapSize', 0, 'linewidth', wd_border/2)
for isubj = 1:nsubj
    plot(EEI_x_med_allSubj(isubj), EEI_y_med_allSubj(isubj), markers_allSubj{isubj}, ...
        'MarkerSize', sz_marker, 'MarkerFaceColor', 'w', 'MarkerEdgeColor',  'k', 'linewidth', wd_border/1.5)
end

%  linear regression
lm = polyfit(EEI_x_med_allSubj, EEI_y_med_allSubj, 1);
x_lm2 = linspace(min(EEI_x_med_allSubj) - std(EEI_x_med_allSubj), max(EEI_x_med_allSubj) + std(EEI_x_med_allSubj), 2);
yfit = polyval(lm, x_lm2);
plot(x_lm2, yfit, 'k', 'handlevisibility', 'off', 'linewidth', wd_border*1.5);

% plot the CI of ln regression
%                 mdl = fitlm(table(x_med(:), y_med(:)));
%                 lr = plot(mdl,'handlevisibility', 'off');
%                 lr(4).Color='k'; lr(4).LineWidth=1.5;
%                 lr(3).Color='k'; lr(3).LineWidth=1.5;

%             % plot linear regression
%             plot(x_lm2, yfit, 'k-', 'handlevisibility', 'off', 'linewidth', 2,'handlevisibility', 'off')

% extra lines
% if (ifeature == 1) && (iparam == 3), yline(0, 'color', ones(1,3)*.5); end

ax = gca;
ax.XAxis.FontSize = sz_ticks;
ax.YAxis.FontSize = sz_ticks;
ax.LineWidth = wd_border;

xlabel([])
ylabel(sprintf('EEI for pA (%s)', namesType{itype}), 'fontsize', sz_label),

% get corr
[r, p] = corr(EEI_x_med_allSubj, EEI_y_med_allSubj);
if p<0.05, title_sig='_sig'; else, title_sig=''; end
text_corr = sprintf('r = %.2f (p = %.3f)', round(r,2), round(p,3));

%             xline(0, 'color', ones(1,3)*.5, 'linewidth', 1.5);
yline(0, 'color', ones(1,3)*.5, 'linewidth', wd_border);

axis square

title(sprintf('[%s] %s vs. %s %s', namesType{itype}, namesLocComb{iLocComb_all(1)}, namesLocComb{iLocComb_all(2)}, text_corr), 'fontsize', sz_title)
%
% folderName = sprintf('%s/%s/corrEEI_pA/%s/', nameFigFolder, nameEnergySource, nameFileLoc);
% folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
% 
% saveas(gcf, sprintf('%sn%d_%s.jpg', folderName, nsubj, namesType{itype}))


