

xlim_allM = [.7, 1.6; -1, 1; .5, 1; .5, 1; .1, .5; .5, .9; .5, .9;.5, .9];
figure('Position', [0 0 1000 600])

for im = 1:nmetrics
    subplot(2,4,im), hold on
    plot(metrics2_test_allB(:, im), metrics2_pred_allB(:, im), 'o')
    plot([xlim_allM(im, 1), xlim_allM(im, 2)], [xlim_allM(im, 1), xlim_allM(im, 2)], 'k-')
    xlim(xlim_allM(im, :))
    ylim(xlim_allM(im, :))
    axis square
    title(namesMetrics{im})
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
sgtitle(sprintf('%s-%s', subjName, namesLocComb{iLocComb}))

%%
if flagIncludePA, namesParamsMuSigma = {'Thresh', 'Mu PRS', 'sigma PRS', 'Mu ABS', 'sigma ABS', 'pA buffer'};
else, namesParamsMuSigma = {'Thresh', 'Mu PRS', 'sigma PRS', 'Mu ABS', 'sigma ABS'};
end
np = length(namesParamsMuSigma);

figure('Position', [0 0 1200 200])
for ip = 1:np
    subplot(1,np,ip), hold on
    histogram(p_est_allB(:, ip))
    title(sprintf('%s\n%.2f', namesParamsMuSigma{ip}, median(p_est_allB(:, ip))))
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)