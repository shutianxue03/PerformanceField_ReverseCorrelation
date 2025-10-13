function ModelPlot_hist(iLoc, nLoc, IV_PRS, IV_ABS, params_est)
% ModelPlot_hist(iLoc, nLoc, IV_PRS, IV_ABS, params_est)
namesLocComb = {'Fovea', 'LHM', 'UVM', 'RHM', 'LVM', 'HM', 'VM', 'Peri'};

if ndims(params_est)==3
    params_est_ave = squeeze(median(params_est(:, iLoc, :), 1));
    params_est_std = squeeze(std(params_est(:, iLoc, :), [], 1));
else
    params_est_ave = params_est(iLoc, :);
    params_est_std = params_est(iLoc, :);
end

nparams_model = length(params_est_ave);
if nLoc == 1, figure, else, subplot(1,nLoc,iLoc), end

hold on
p_prs = histcounts(IV_PRS, 'Normalization', 'probability');
p_abs = histcounts(IV_ABS, 'Normalization', 'probability');
% histogram
histogram(IV_PRS, 'FaceColor', 'r', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(IV_ABS, 'FaceColor', 'b', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability');
% distribution mean
xline(mean(IV_PRS), 'r-', 'linewidth', 1.5);
xline(mean(IV_ABS), 'b-', 'linewidth', 1.5);
% distribution median
xline(median(IV_PRS), 'r--', 'linewidth', 1.5);
xline(median(IV_ABS), 'b--', 'linewidth', 1.5);
% width - external noise SD
errorbar(mean(IV_PRS), max(p_prs)/2, std(IV_PRS), 'r-', 'horizontal', 'linewidth', 1.5)
errorbar(mean(IV_ABS), max(p_abs)/2, std(IV_ABS), 'b-', 'horizontal', 'linewidth', 1.5)


switch nparams_model
    case 1
        % width - internal + external noise SD
        errorbar(mean(IV_PRS), max(p_prs)/2, std(IV_PRS)*(1+ params_est_ave), 'r', 'horizontal', 'CapSize',0, 'linewidth', 1.5)
        errorbar(mean(IV_ABS), max(p_abs)/2, std(IV_ABS)*(1+ params_est_ave), 'b', 'horizontal', 'CapSize',0, 'linewidth', 1.5)
        % estimated threshold
        xline(0, 'k', 'linewidth', 2);
    case 2
        % width - internal + external noise SD
        errorbar(mean(IV_PRS), max(p_prs)/2, std(IV_PRS)*(1+ params_est_ave(1)), 'r', 'horizontal', 'CapSize',0, 'linewidth', 1.5)
        errorbar(mean(IV_ABS), max(p_abs)/2, std(IV_ABS)*(1+ params_est_ave(1)), 'b', 'horizontal', 'CapSize',0, 'linewidth', 1.5)
        % estimated threshold
        xline(params_est_ave(2), 'k', 'linewidth', 2);
        errorbar(params_est_ave(2), .04, params_est_std(2), 'k', 'horizontal', 'linewidth', 2)
    case 3
        errorbar(mean(IV_PRS), max(p_prs)/2, std(IV_PRS)*(1+ params_est_ave(1)), 'r', 'horizontal', 'CapSize',0, 'linewidth', 1.5)
        errorbar(mean(IV_ABS), max(p_abs)/2, std(IV_ABS)*(1+ params_est_ave(2)), 'b', 'horizontal', 'CapSize',0, 'linewidth', 1.5)
        % estimated threshold
        xline(params_est_ave(3), 'k', 'linewidth', 2);
        errorbar(params_est_ave(3), .04, params_est_std(2), 'k', 'horizontal', 'linewidth', 2)
end
title(namesLocComb{iLoc})
end