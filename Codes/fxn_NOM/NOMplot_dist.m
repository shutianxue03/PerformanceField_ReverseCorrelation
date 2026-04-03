

figure('Position', [0 200 600 1e3]);
subplot(2,1,1); hold on;
histogram(DVnoisy_sim_allT(iPRS_allT == 1), 'FaceColor', 'r', 'DisplayName', 'Signal Present', 'normalization', 'probability');
histogram(DVnoisy_sim_allT(iPRS_allT == 0), 'FaceColor', 'b', 'DisplayName', 'Signal Absent', 'normalization', 'probability');
xline(criterion_DV_true, 'LineWidth', 2, 'DisplayName', 'True criterion');
xlim([min(DVnoisy_sim_allT), max(DVnoisy_sim_allT)]);
ylabel('Proportion');
legend('show', 'location', 'best');

subplot(2,1,2); hold on;
plot(DVclean_sim_allT, resp_allT, 'ro', 'DisplayName', 'Binary response');
plot(DVclean_sim_allT, pYES_pred_allT, 'k+', 'DisplayName', 'Pred pYES');
xline(criterion_DV_true, 'LineWidth', 2, 'DisplayName', 'True criterion');
xlabel('IV'); ylabel('pYES');
yline(0.5, 'k--');
xlim([min(DVnoisy_sim_allT), max(DVnoisy_sim_allT)]);
metrics_sim_ = metrics_sim; metrics_sim_(3:end) = metrics_sim_(3:end)*100;

sgtitle(sprintf([ ...
    'Simulated data & metrics' ...
    '        \nCI_{95}=[%.2f, %.2f], Median=%.2f\n[TRUE] gN=%.2f, gC=%.2f, Nmul=%s, Nadd=%s, Nshared=%s' ...
    '       \n *** Criterion: True: %.2f (%.2f) | Simulated: %.2f ***'...
    '       \n[MEASURED] pYES=%.2f, pC=%.2f, pHit=%.2f, pFA=%.2f, pA=%.2f, corrPass=%.2f'], ...
    round(quantile(DVnoisy_sim_allT, [.05, .95, .5]), 1), noiseCST, gaborCST, ...
    format_num2exp(Nmul_true), format_num2exp(Nadd_true), format_num2exp(Nshared_true), ...
    Cz_true, criterion_DV_true, Cz_sim, ...
    pYES_sim, pC_sim, pHit_sim, pFA_sim, pA_sim, corrPass), ...
    'fontsize', 12);

saveas(gcf, sprintf('%s/0PerfHist.jpg', nameFolder_Figures_perSubj));
close all;
fprintf('%s: Distribution of DV drawn.\n\n', datetime('now'))