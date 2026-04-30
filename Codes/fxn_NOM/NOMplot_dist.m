

figure('Position', [0 200 2e3 1e3]);
% 1. Distribution of DVs
subplot(2,2,1); hold on;
histogram(DVnoisy_sim_allT(iPRS_allT == 1), 'FaceColor', 'r', 'DisplayName', 'Signal Present', 'normalization', 'probability');
histogram(DVnoisy_sim_allT(iPRS_allT == 0), 'FaceColor', 'b', 'DisplayName', 'Signal Absent', 'normalization', 'probability');
xline(criterion_DV_true, 'LineWidth', 2, 'DisplayName', 'True criterion');
xlim([min(DVnoisy_sim_allT), max(DVnoisy_sim_allT)]);
ylabel('Proportion');
legend('show', 'location', 'best');

% 2. Response and  pYES as a fxn of clean DV
subplot(2,2,3); hold on;
plot(DVclean_sim_allT, resp_allT, 'ro', 'DisplayName', 'Binary response');
plot(DVclean_sim_allT, pYES_pred_allT, 'k+', 'DisplayName', 'Pred pYES');
xline(criterion_DV_true, 'LineWidth', 2, 'DisplayName', 'True criterion');
xlabel('IV'); ylabel('pYES');
yline(0.5, 'k--');
xlim([min(DVnoisy_sim_allT), max(DVnoisy_sim_allT)]);
metrics_sim_ = metrics_sim; metrics_sim_(3:end) = metrics_sim_(3:end)*100;

% 3. Correlation between DV
subplot(2,2,2), hold on 
xPRS = DVnoisy_sim_allT(iPass_allT==1 & iPRS_allT==1);
yPRS = DVnoisy_sim_allT(iPass_allT==2 & iPRS_allT==1);
scatter(xPRS, yPRS, 'ro')
xABS = DVnoisy_sim_allT(iPass_allT==1 & iPRS_allT==0);
yABS = DVnoisy_sim_allT(iPass_allT==2 & iPRS_allT==0);
scatter(xABS, yABS, 'bo')
% unity line
plot([min(DVnoisy_sim_allT), max(DVnoisy_sim_allT)], [min(DVnoisy_sim_allT), max(DVnoisy_sim_allT)], 'k-')
axis square
xlim([min(DVnoisy_sim_allT), max(DVnoisy_sim_allT)])
ylim([min(DVnoisy_sim_allT), max(DVnoisy_sim_allT)])
xlabel('DV of Pass 1');
ylabel('DV of Pass 2');

corrPass_PRS = corr(xPRS, yPRS);
corrPass_ABS = corr(xABS, yABS);
x = DVnoisy_sim_allT(iPass_allT==1);
y = DVnoisy_sim_allT(iPass_allT==2);
z = iPRS_allT(iPass_allT==1);
corrPass_partial = partialcorr(x, y, z);
title(sprintf('Correlation between noisy DV of two passes (r=%.2f)\nPRS: r=%.2f; ABS: r=%.2f', corrPass_partial, corrPass_PRS, corrPass_ABS))

sgtitle(sprintf([ ...
    'Simulated data & metrics' ...
    '        \nMedian=%.2f, SD=%.2f, CI_{95}=[%.2f, %.2f], \n[TRUE] gN=%.2f, gC=%.2f, Nmul=%.1f, Nadd=%.1f, Nshared=%.1f' ...
    '       \n *** Criterion: True: %.2f (%.2f) | Simulated: %.2f ***'...
    '       \n[MEASURED] pYES=%.2f, pC=%.2f, pHit=%.2f, pFA=%.2f, pA=%.2f'], ...
    std(DVnoisy_sim_allT), median(DVnoisy_sim_allT), quantile(DVnoisy_sim_allT, [.05, .95]), noiseCST, gaborCST, ...
    Nmul_true, Nadd_true, Nshared_true, ...
    cSDT_true, criterion_DV_true, cSDT_sim, ...
    pYES_sim, pC_sim, pHit_sim, pFA_sim, pA_sim), ...
    'fontsize', 12);

saveas(gcf, sprintf('%s/0PerfHist_pC%.0f.jpg', nameFolder_Figures_perSubj, pC_sim*100));
close all;
fprintf('%s: Distribution of DV drawn.\n\n', datetime('now'))