% function PRplot_paramsCorr(params_est_allB, params_true, labels, ni, params_lb, params_ub)

nparams = length(params_true);
if nparams == 3, cc_all = [1,2;2,3;1,3];
else, cc_all = [1,2];
end
ncc = size(cc_all, 1);
figure('Position',[0 0 400*ncc 400])
for ic = 1:ncc
    ind = cc_all(ic, :);
    x = params_est_allB(:,ind(1));
    y = params_est_allB(:,ind(2));
    
    subplot(1,ncc,ic), hold on
    plot(x, y,'k.')
    xline(params_true(ind(1)), 'r-', 'linewidth', 2);
    yline(params_true(ind(2)), 'r-', 'linewidth', 2);
    
    [x_med, ~, ~, x_neg, x_pos] = getCI(x);
    [y_med, ~, ~, y_neg, y_pos] = getCI(y);
    errorbar(x_med, y_med, y_neg, y_pos,x_neg, x_pos,  'k', 'CapSize', 0, 'Linewidth', 2)
    
    xlabel(namesParamsModel{ind(1)})
    ylabel(namesParamsModel{ind(2)})
    [r,p] = corr(x, y);
%     if nargin == 6
%     xlim([params_lb(ind(1)), params_ub(ind(1))])
%     ylim([params_lb(ind(2)), params_ub(ind(2))])
%     end
    axis square 
    title(sprintf('r=%.3f, p=%.3f',r, p))
    
end
sgtitle(sprintf('Correlation between estimated params (nB = %d)', ni))

