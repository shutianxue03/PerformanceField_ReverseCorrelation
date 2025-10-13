
function [y_pred, params_est, Rsquared, nLL] = SX_sim08_fit(x, y, imodel, ub, lb, fitMode)

% at each each time of fitting, param0 should be sample randomly from [lb, ub]
sampParam0 = @(lb, ub) rand(1,length(ub)) .* (ub-lb) + lb;

% options = optimoptions('fmincon','Algorithm','interior-point','MaxIter',5000,'Display','off');
options = optimoptions('fmincon','MaxIterations', 5000,'Display','off');

params0 = sampParam0(lb, ub); % set param0

kernelFit_fxn = @(kernelParams) kernelFit(kernelParams, x, y, imodel, fitMode);

[params_est, nLL] = fmincon(kernelFit_fxn, params0, [], [], [], [], lb, ub, [], options);
y_pred = predSFkernel(x, imodel, params_est, 0);
Rsquared = 1-sumsqr(y-y_pred)/sumsqr(y-mean(y));

