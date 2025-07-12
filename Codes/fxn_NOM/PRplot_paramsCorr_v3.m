function PRplot_paramsCorr_v3(params_est_allB, params_true, labels, ni, params_lb, params_ub)

ind = [1,2];
x = params_est_allB(:,ind(1));
y = params_est_allB(:,ind(2));

figure('Position', [2000 0 500 400])
hold on
plot(x, y,'o')
xline(params_true(ind(1)), 'r-');
yline(params_true(ind(2)), 'r-');
xlabel(labels{ind(1)})
ylabel(labels{ind(2)})
[r,p] = corr(x, y);
if nargin == 6
    xlim([params_lb(ind(1)), params_ub(ind(1))])
    ylim([params_lb(ind(2)), params_ub(ind(2))])
end
title(sprintf('Correlation between estimated params (nB = %d)\nr=%.3f, p=%.3f', ni, r, p))

set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
