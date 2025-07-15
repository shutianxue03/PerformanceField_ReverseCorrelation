

wd = 5;
wd_bar = 2;
sz_marker_ave = 30;
sz_marker_idvd = 12;
sz_ticks = 20;

[s_med_allSubj, ~, ~, s_neg_allSubj, s_pos_allSubj] = getCI(similarityAfterMirroring_allSubj(:,:, itype, :), 1, 2);
[s_ave, ~, ~, s_sem] = getCI(s_med_allSubj);

%% group ave
figure('Position', [0 200 400 300]),hold on

for iLoc = 1:nLoc
    errorbar(iLoc, s_ave(iLoc), s_sem(iLoc), '.', ...
        'color', colors_comb(iLocComb_all(iLoc), :), 'CapSize', 0, 'linewidth', wd)
    plot(iLoc, s_ave(iLoc), 'o', ...
        'MarkerEdgeColor', colors_comb(iLocComb_all(iLoc), :), 'MarkerSize', sz_marker_ave, 'linewidth', wd)
end % end of iiLoc

%% idvd
buffer = .2;
for isubj = 1:nsubj
    plot([1+buffer, 2-buffer], s_med_allSubj(isubj, :),  [markers_allSubj{isubj}, '-'], ...
        'color',ones(1,3)*.7, 'markerfacecolor', 'w', 'markeredgecolor', ones(1,3)*.7, ...
        'markersize', sz_marker_idvd, 'linewidth', wd_bar)
end

xticks([])
buffer = .3;
xlim([1-buffer, nLoc+buffer])
ylim([.5, 1])

ax = gca;
ax.XAxis.FontSize = sz_ticks;
ax.YAxis.FontSize = sz_ticks;
ax.LineWidth = wd;
ylabel('Similarity')
title(sprintf('%s: %.2f (%.2f)\n%s: %.2f (%.2f)', namesLocComb{iLocComb_all(1)}, s_ave(1), s_sem(1), namesLocComb{iLocComb_all(2)}, s_ave(2), s_sem(2)))

folderName = sprintf('%s/%s/kernels2D/', nameFigFolder, name_numFilters_Fitting);
folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
saveas(gcf, sprintf('%sn%d_L%d%d_similarity.jpg', folderName, nsubj, iLocComb_all))
