function RCplot_2Dkernel(data2D, caxisLim, outline_pos, outline_neg)
% RCplot_2Dkernel(data2D, caxisLim, outline)

%----------------%
SX_RC1_setting
%----------------%

wd_ref = 6;
wd_contour = 4;
wd_border = 2;
sz_ticks = 50;
sz_colorbar = 40;

hold on
imagesc(axis_tuning{1}, axis_tuning{2}, data2D)
xline(axisTicks_tuning{1}(3), 'r-', 'linewidth', wd_ref);
yline(axisTicks_tuning{2}(3), 'r-', 'linewidth', wd_ref);

xticks(axisTicks_tuning{1})
xticklabels(axisTL_tuning{1})
yticks(axisTicks_tuning{2})
yticklabels(axisTL_tuning{2})

xlim(axisTicks_tuning{1}([1,end]))
ylim(axisTicks_tuning{2}([1,end]))

ch = colorbar;
axis square

% Delineate the pixels for multiple comparison
if nargin >= 3
    contour(axis_tuning{1}, axis_tuning{2}, outline_pos, 1, 'r', 'linewidth', wd_contour)
end

if nargin == 4
    contour(axis_tuning{1}, axis_tuning{2}, outline_neg, 1, 'k', 'linewidth', wd_contour)
end

if ~isnan(caxisLim)
    caxis(caxisLim) % make sure this is at the end!!
    ch.Ticks = round(linspace(caxisLim(1), caxisLim(2), 4), 2);
end
ch.FontSize = sz_colorbar;

ax = gca;
ax.XAxis.FontSize = sz_ticks;
ax.YAxis.FontSize = sz_ticks;
% ax.LineWidth = wd_border;
