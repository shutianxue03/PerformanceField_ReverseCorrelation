
if diffFlag, names = PFnames; iname = ititle; else, names = locNames; iname = iLoc; end

%% additional lines
yline(0, 'handlevisibility', 'off');
xline(line_extra(ikernel), 'handlevisibility', 'off');

%% xaxis
xlim(xticks_{ikernel}([1, end]) + ybuffer{ikernel})
xlabel(xlabels{ikernel}, 'FontSize', sz_label)
xticks(xticks_{ikernel})
xticklabels(xticklabels_{ikernel})

%% yaxis
ylim([yrange{ikernel}(1), yrange{ikernel}(end)])
ylabel('kernel (a.u.)', 'FontSize', sz_label)
yticks(yrange{ikernel})

ax = gca; ax.FontSize = sz_ticks;

% %% legend
% if nsubj == 1, legend('kernel', 'model')
% else, legend('average', '\pm 1 SEM', 'model')
% end
% 
% %% title
% title(names{iname}, 'FontSize', sz_title)
% 
% %%
% set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 2)



