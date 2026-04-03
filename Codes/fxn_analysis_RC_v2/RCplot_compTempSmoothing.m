
template_ideal_norm = template_ideal/sum(template_ideal(:));
template_noSmoothing_norm = template_noSmoothing/sum(template_noSmoothing(:));
template_smooth_norm = template_smooth/sum(template_smooth(:));

%%
figure('Position', [0 0 2e3 1e3]);
subplot(1,3,1)
imagesc(axis_tuning{2}, axis_tuning{1}, template_ideal_norm);
axis xy;
xlabel('SF');
ylabel('ORI');
colorbar;
axis square
yline(axis_tuning{1}((nORI+1)/2), 'r-', 'LineWidth', 2);
xline(axis_tuning{2}((nSF+1)/2), 'r-', 'LineWidth', 2);
title('Ideal template');

subplot(1,3,2)
imagesc(axis_tuning{2}, axis_tuning{1}, template_noSmoothing_norm);
axis xy;
xlabel('SF');
ylabel('ORI');
colorbar;
axis square
yline(axis_tuning{1}((nORI+1)/2), 'r-', 'LineWidth', 2);
xline(axis_tuning{2}((nSF+1)/2), 'r-', 'LineWidth', 2);
title('No smoothening RC template');

subplot(1,3,3)
imagesc(axis_tuning{2}, axis_tuning{1}, template_smooth_norm);
axis xy;
xlabel('SF');
ylabel('ORI');
colorbar;
axis square
yline(axis_tuning{1}((nORI+1)/2), 'r-', 'LineWidth', 2);
xline(axis_tuning{2}((nSF+1)/2), 'r-', 'LineWidth', 2);
title('Smooth-basis RC template');

%%
figure
subplot(1,2,1), hold on, plot(template_noSmoothing_norm, template_ideal_norm, 'o'),
axis_max = max([template_noSmoothing_norm(:); template_ideal_norm(:)]);
plot([0, axis_max], [0, axis_max], 'k-')
xlim([0, axis_max])
ylim([0, axis_max])

subplot(1,2,2), hold on, plot(template_smooth_norm, template_ideal_norm, 'o'),
axis_max = max([template_smooth_norm(:); template_ideal_norm(:)]);
plot([0, axis_max], [0, axis_max], 'k-')
xlim([0, axis_max])
ylim([0, axis_max])