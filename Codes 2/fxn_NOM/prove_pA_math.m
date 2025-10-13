

%% params to vary
ntrials_pA = 1e4;

mu_PRS = 2;
sigma_PRS = 3;
mu_ABS = -2;
sigma_ABS = 2;
thresh = 0;

%%
PDF_PRS = @(x) normpdf(x, mu_PRS, sigma_PRS);
PDF_ABS = @(x) normpdf(x, mu_ABS, sigma_ABS);
CDF_PRS = normcdf(thresh, mu_PRS, sigma_PRS);
CDF_ABS = normcdf(thresh, mu_ABS, sigma_ABS);

fxn_PRS = @(x) PDF_PRS(x) .* ((CDF_PRS).^2 + (1-CDF_PRS).^2);
fxn_ABS = @(x) PDF_ABS(x) .* ((CDF_ABS).^2 + (1-CDF_ABS).^2);
pA_PRS_math = integral(fxn_PRS, -inf, inf);
pA_ABS_math = integral(fxn_ABS, -inf, inf);

%% simulate trials to calculate pA
nrep = 1e3;
pA_PRS_sim = nan(1,nrep);
pA_ABS_sim = nan(1,nrep);
parfor irep = 1:nrep
    IV_PRS_pA = randn(1, ntrials_pA) * sigma_PRS + mu_PRS;
    IV_ABS_pA = randn(1, ntrials_pA) * sigma_ABS + mu_ABS;
    
    %get internal responses (1=YES, 0=NO)
    IR_PRS = IV_PRS_pA >= thresh; % so that 1 is the correct answer
    IR_ABS = IV_ABS_pA >= thresh; % so that 0 is the correct answer
    
    % split trials into two passes
    IR_PRS_pass = reshape(IR_PRS, [2,ntrials_pA/2]);
    IR_ABS_pass = reshape(IR_ABS, [2,ntrials_pA/2]);
    
    % predict pA3
    pA_PRS_sim(irep) = mean(IR_PRS_pass(1,:) == IR_PRS_pass(2,:));
    pA_ABS_sim(irep) = mean(IR_ABS_pass(1,:) == IR_ABS_pass(2,:));
    % pA_sim = mean([pA_PRS_sim, pA_ABS_sim]);
end
fprintf('DONE\n')

%% PLOT
figure('Position', [200 200 900 300])
subplot(1,3,1), hold on
x = linspace(-10,10, 1e4);
plot(x, normpdf(x, mu_PRS, sigma_PRS), 'r')
plot(x, normpdf(x, mu_ABS, sigma_ABS), 'b')
xline(thresh, 'k-');
title('Variable Distribution')

subplot(1,3,2), hold on
histogram(pA_PRS_sim, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'Normalization', 'Probability')
xline(median(pA_PRS_sim), 'k--', 'linewidth', 1);
xline(pA_PRS_math, 'r', 'linewidth', 2);
legend({'Simulated pA', 'Median', 'Theoretical pA'}, 'Location', 'south')
% ylim([0, .1])
xlabel('pA [PRS]')
title(sprintf('[PRS] Distribution of simulated pA\nstd=%.4f', std(pA_PRS_sim)))

subplot(1,3,3), hold on
histogram(pA_ABS_sim, 'DisplayStyle', 'stair','EdgeColor', 'b', 'Normalization', 'Probability')
xline(median(pA_ABS_sim), 'k--', 'linewidth', 1);
xline(pA_ABS_math, 'b', 'linewidth', 2);
% ylim([0, .1])
xlabel('pA [ABS]')
title(sprintf('[ABS] Distribution of simulated pA\nstd=%.4f', std(pA_ABS_sim)))

sgtitle(sprintf('[1e%d trials] thresh = %.2f\nPRS ~ N(%.2f, %.2f) ABS ~ N(%.2f, %.2f)', round(log10(ntrials_pA)), thresh, mu_PRS, sigma_PRS, mu_ABS, sigma_ABS))

set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)