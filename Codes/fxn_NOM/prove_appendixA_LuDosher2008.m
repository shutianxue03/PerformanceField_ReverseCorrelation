% SIMULATE to test Eq A8&9in Appendix A (p32) of Lu % Dosher 2008
% calculate pC given a 2AFC discrm task
% Two ways to calculate the theoretical pC in a 2 AFC task produced the SAME result

mu1 = -1;
mu2 = 0;
sigma1 = 1;
sigma2 = 1.2;

% Way 1 (Eq A8)
pC1 = 1-normcdf(0, (mu2-mu1)/sqrt(sigma1^2+sigma2^2), 1);

% Way 2 (Eq A9)
PDF = @(x) normpdf(x, mu2, sigma2);
CDF = @(x) normcdf(x, mu1, sigma1);
fxn_pC2 = @(x) PDF(x) .* CDF(x);
pC2 = integral(fxn_pC2, -inf, inf);

if abs(pC1-pC2) < 1e-10, fprintf('pC calculated in two ways are the SAME.\n'), end

%% understand why Way 2 is right by simulation
clc
ntrials = 1e4;
ni = 1e4;
pC_sim = nan(1,ni);
parfor ii = 1:ni
    trials1 = randn(1, ntrials) * sigma1+mu1;
    trials2 = randn(1, ntrials) * sigma2+mu2;
    
    % the estimated pC, should be similar to pC
    pC_sim(ii) = mean(trials2>trials1);
end

figure, hold on
histogram(pC_sim)
xline(pC2, 'r-', 'linewidth', 2);
legend({'Simulated pC', 'Math pC'})





