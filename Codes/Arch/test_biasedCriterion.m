%% Simulation aligned with RC method (signal-absent only, probit GLM, z-scored energy)
clear; clc; rng(0);

N = 5000;                      % # signal-absent trials (per location/contrast)
E = randn(N,1);                % z-scored energy: mean 0, sd 1 (matches \tilde{E})

beta1_true = 1;              % "representation strength" (slope)
beta0_unbiased  = 0.0;              % less liberal / closer to unbiased
beta0_liberal  = 1.0;              % more liberal criterion (higher P(yes) overall)

% Generate responses under probit model: P(yes)=Phi(beta0 + beta1*E)
p_unb = normcdf(beta0_unbiased + beta1_true*E);
p_lib = normcdf(beta0_liberal + beta1_true*E);
R_unb = rand(N,1) < p_unb;
R_lib = rand(N,1) < p_lib;

fprintf('Yes-rate unbiased: %.3f\n', mean(R_unb));
fprintf('Yes-rate liberal : %.3f\n', mean(R_lib));

%% Fit GLM WITH intercept (matches your Eq. 4)
B_unb = glmfit(E, R_unb, 'binomial', 'link', 'probit');  % B(1)=beta0, B(2)=beta1
B_lib = glmfit(E, R_lib, 'binomial', 'link', 'probit');

fprintf('\nGLM with intercept (your method):\n');
fprintf('  True beta0: unb %.2f, lib %.2f\n', beta0_unbiased, beta0_liberal);
fprintf('  Hat  beta0: unb %.3f, lib %.3f\n', B_unb(1), B_lib(1));
fprintf('  True beta1: %.2f\n', beta1_true);
fprintf('  Hat  beta1: unb %.3f, lib %.3f\n', B_unb(2), B_lib(2));

%% Plot simulated binary data and GLM fits
figure; hold on;
scatter(E, R_unb, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.1);
scatter(E, R_lib, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.1);
Efit = linspace(-3,3,100)';
pfit_unb = normcdf(B_unb(1) + B_unb(2)*Efit);
pfit_lib = normcdf(B_lib(1) + B_lib(2)*Efit);
plot(Efit, pfit_unb, 'b-', 'LineWidth', 2);
plot(Efit, pfit_lib, 'r-', 'LineWidth', 2);
xlabel('z-scored energy E'); ylabel('P(yes)');
title('Simulated data and GLM fits with intercept');
legend({sprintf('Data unbiased (slope=%.2f, criterion=%.2f)', B_unb(2), B_unb(1)), sprintf('Data liberal (slope=%.2f, criterion=%.2f)', B_lib(2), B_lib(1)), ...
    'GLM unbiased', 'GLM liberal'}, 'Location','best');
grid on;
%% Optional: repeat many times to show distribution of beta1 estimates
nRep = 200;
b1hat_unb = nan(nRep,1);
b1hat_lib = nan(nRep,1);
for r = 1:nRep
    E = randn(N,1);
    p_unb = normcdf(beta0_unbiased + beta1_true*E);
    p_lib = normcdf(beta0_liberal + beta1_true*E);
    R_unb = rand(N,1) < p_unb;
    R_lib = rand(N,1) < p_lib;

    B_unb = glmfit(E, R_unb, 'binomial', 'link', 'probit');
    B_lib = glmfit(E, R_lib, 'binomial', 'link', 'probit');
    b1hat_unb(r) = B_unb(2);
    b1hat_lib(r) = B_lib(2);
end

fprintf('\nAcross %d repetitions:\n', nRep);
fprintf('  mean(beta1hat): unb %.3f, lib %.3f\n', mean(b1hat_unb), mean(b1hat_lib));
fprintf('  sd(beta1hat)  : unb %.3f, lib %.3f\n', std(b1hat_unb), std(b1hat_lib));

figure; hold on;
histogram(b1hat_unb, 25); histogram(b1hat_lib, 25);
legend({'Unbiased','Liberal'}, 'Location','best'); grid on;
xlabel('\beta_1 estimate'); ylabel('count');
title('Probit GLM slope (\beta_1): should not systematically drop with liberal criterion when intercept is included');
