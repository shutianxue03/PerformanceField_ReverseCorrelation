
itest = 4; % IVs are boxcox-transformed
figure('Position', [3e3 0 1e3 500]), hold on

%% histogram of IV
h1 = histfit(IV_PRS, 20); h1(1).FaceColor = 'r'; h1(1).LineWidth = .5; h1(1).FaceAlpha = .5; h1(2).Color = 'r'; h1(2).LineWidth = 2;
h2 = histfit(IV_ABS, 20); h2(1).FaceColor = 'b'; h2(1).LineWidth = .5; h2(1).FaceAlpha = .5; h2(2).Color = 'b'; h2(2).LineWidth = 2;

%% threshold
xline(params_est(3), 'k', 'linewidth', 5);

%% theoretical function of noisyIV
N_mul = params_est(1);
sd_add = params_est(2);
sigma_IE_PRS = sqrt((1+N_mul^2) * sigma_PRS^2 + sd_add^2);
sigma_IE_ABS = sqrt((1+N_mul^2) * sigma_ABS^2 + sd_add^2);
x = linspace(min([IV_PRS;IV_ABS]), max([IV_PRS;IV_ABS]), 1e3);
noisyIV_PRS = normpdf(x, mu_PRS, sigma_IE_PRS); noisyIV_PRS = noisyIV_PRS/max(noisyIV_PRS) * max(h1(2).YData);
noisyIV_ABS = normpdf(x, mu_ABS, sigma_IE_ABS);noisyIV_ABS = noisyIV_ABS/max(noisyIV_ABS) * max(h2(2).YData);
plot(x, noisyIV_PRS, 'r*-')
plot(x, noisyIV_ABS, 'b*-')

%%
title(sprintf('%s\nIV: [PRS]~N(%.2f, %.3f) [ABS]~N(%.2f, %.3f)\nFitted: [PRS]~N(%.2f, %.3f) [ABS]~N(%.2f, %.3f)\nEst: %.3f, %.3f, %.3f (lb: %.3f, ub: %.3f)', ...
    namesNormalityTest{itest}, mu_PRS, sigma_PRS, mu_ABS, sigma_ABS, ...
    f1.mu,  f1.sigma, f2.mu,  f2.sigma, ...
    params_est, params_lb(end), params_ub(end)))
set(findall(gcf, '-property', 'fontsize'), 'fontsize',20)

