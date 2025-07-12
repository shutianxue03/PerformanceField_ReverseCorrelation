function quickPlot_hist(IV_PRS, IV_ABS, IV_PRS_noisy, IV_ABS_noisy, IV_PRS_noisy_sim, IV_ABS_noisy_sim, mu_PRS, mu_ABS, sigma_IE_PRS, sigma_PRS, sigma_IE_ABS,  sigma_ABS)

%% histogram
subplot(2,2,[1,3]), hold on
histogram(IV_PRS, 50, 'DisplayStyle', 'stair', 'EdgeColor', 'r', 'linewidth', 2, 'Normalization', 'Probability')
histogram(IV_ABS, 50, 'DisplayStyle', 'stair', 'EdgeColor', 'b', 'linewidth', 2, 'Normalization', 'Probability')
histogram(IV_PRS_noisy, 50, 'DisplayStyle', 'bar', 'linestyle', 'none', 'FaceColor', 'r', 'Normalization', 'Probability')
histogram(IV_ABS_noisy, 50, 'DisplayStyle', 'bar', 'linestyle', 'none', 'FaceColor', 'b', 'Normalization', 'Probability')
histogram(IV_PRS_noisy_sim, 50, 'DisplayStyle', 'stair', 'EdgeColor', 'k', 'linewidth', 1.5, 'Normalization', 'Probability')
histogram(IV_ABS_noisy_sim, 50, 'DisplayStyle', 'stair', 'EdgeColor', 'k', 'linewidth', 1.5, 'Normalization', 'Probability')
legend({'PRS raw', 'ABS raw', 'Noisy PRS [math]', 'Noisy ABS [math]', 'Noisy PRS [sim]', 'Noisy ABS [sim]'}, ...
    'Location', 'southoutside')

%% mean and std
subplot(2,2,2), hold on
bar(1:2, [median(IV_PRS), median(IV_PRS_noisy), median(IV_PRS_noisy_sim);...
    median(IV_ABS), median(IV_ABS_noisy), median(IV_ABS_noisy_sim)])
yline(mu_PRS, 'k-');
yline(mu_ABS, 'k-');
xticks(1:2), xticklabels({'PRS median', 'ABS median'})
legend({'raw', 'noisy', 'noisy sim'}, 'Location', 'south')
title('Medians (should be the same)')

% sigma_IE_PRS, sigma_PRS, sigma_IE_ABS,  sigma_ABS
subplot(2,2,4), hold on
b = bar(1:2, [std(IV_PRS), std(IV_PRS_noisy), std(IV_PRS_noisy_sim);...
             std(IV_ABS), std(IV_ABS_noisy), std(IV_ABS_noisy_sim)]);
plot(1,sigma_IE_PRS, 'ko')
plot(.77,sigma_PRS, 'ko')
plot(2,sigma_IE_ABS, 'ko')
plot(1.77,sigma_ABS, 'ko')

xticks(1:2), xticklabels({'PRS SD', 'ABS SD'})
title('STDs (noisy should = noisy sim and diff from raw)')



