
n = 1e5;



% prove that adding internal noise is to change the sigma from sigma_E to
% sqrt(1+alpha^2) * sigma_E

mu_E = 0;
sigma_E = 2;
alpha = .5; % internal-external ratio
sigma_I = sigma_E*alpha;
sigma_IE_correct = sqrt(1+alpha^2) * sigma_E;
sigma_IE_wrong = (alpha+1) * sigma_E;

IV = randn(1,n) * sigma_E + mu_E;
IV_int_sim = IV + randn(1,n) * sigma_I;
IV_int_correct = randn(1,n) * sigma_IE_correct + mu_E;
IV_int_wrong = randn(1,n) * sigma_IE_wrong + mu_E;

figure, hold on
histogram(IV, 'FaceColor', 'k', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability');
% histogram(IV_int_sim, 'FaceColor', 'r', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(IV_int_correct, 'FaceColor', 'b', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability');
% histogram(IV_int_wrong, 'FaceColor', 'k', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability');
% legend({'IV + sampled internal noise', 'sigma_IE= sqrt(1+alpha^2) * sigma_E', 'sigma_IE = (alpha+1) * sigma_E'})


