
figure('Position', [0 300 800 800])

subplot(3,2,1), hold on
histogram(DV_PRS, 'normalization', 'probability', 'EdgeColor', 'r', 'DisplayStyle', 'stairs')
histogram(DV_ABS, 'normalization', 'probability', 'EdgeColor', 'b', 'DisplayStyle', 'stairs')
xlabel(sprintf('DV (gabor cst=%.1f)', gaborCST_all(end)))

subplot(3,2,3), hold on
plot(gaborCST_all, DV_PRS_mean, 'ro-')
plot(gaborCST_all, DV_ABS_mean, 'b-')
xlabel('Tested Gabor CST'), ylabel('Mean of DV')

subplot(3,2,4), hold on
plot(gaborCST_all, DV_PRS_var, 'r--')
plot(gaborCST_all, DV_ABS_var, 'b--')
xlabel('Tested Gabor CST'), ylabel('Var of DV')

subplot(3,2,5), hold on
histogram(DV_PRS_mean, 20, 'normalization', 'probability', 'EdgeColor', 'r', 'DisplayStyle', 'stairs')
histogram(DV_ABS_mean, 20, 'normalization', 'probability', 'EdgeColor', 'b', 'DisplayStyle', 'stairs')
xlabel('Mean of DV')
legend(sprintf('PRS (mean=%.2f, var=%.2f)', mean(DV_PRS_mean), var(DV_PRS_mean)), ...
    sprintf('ABS (mean=%.2f, var=%.2f)', mean(DV_ABS_mean), var(DV_ABS_mean)))

subplot(3,2,6), hold on
histogram(DV_PRS_var, 10, 'normalization', 'probability', 'EdgeColor', 'r', 'DisplayStyle', 'stairs')
histogram(DV_ABS_var, 10, 'normalization', 'probability', 'EdgeColor', 'b', 'DisplayStyle', 'stairs')
xlabel('Var of DV')
legend(sprintf('PRS (mean=%.2f, var=%.2f)', mean(DV_PRS_var), var(DV_PRS_var)), ...
    sprintf('ABS (mean=%.2f, var=%.2f)', mean(DV_ABS_var), var(DV_ABS_var)))

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
sgtitle([title_, ' - Mean and Var of DV'], 'FontSize',25)
