% function pred = PR_pred_LuDosher2008(params_est, data)
% data.cst
% sigma_ext = .2; % or std(data.IV()
% beta = params_est(1);
% gamma = params_est(2);
% N_mul = params_est(3);
% sd_add = params_est(4);
% 
% 
% [pC, pA] = getpC_pA(cst, sigma_ext, sigma_int);

%%
 [pC, pA] = getpC_pA(.1, .2, .2);
 
%% helper fxn - getpC_pA
function [pC, pA] = getpC_pA(cst, sigma_ext, sigma_int)
% SX_normPDF = @(x,mu,sigma) exp(-(x-mu).^2/(2*sigma^2));

% get x
x = (-1e4:1e4)/100;
sigma_int_PRS = sigma_int;
sigma_int_ABS = sigma_int;

% combine internal & external noise
sigma_IE = sqrt(sigma_ext^2+(sigma_int_PRS^2+sigma_int_ABS^2)/2);
t = cst + x * sigma_IE;

% get pC
pC = sum((sigma_IE)*normpdf(t,cst,sigma_int_PRS).*normcdf(t,0,sigma_int_ABS)/100);

sigma_ext2 = sigma_ext * sqrt(2);
sigma_int2 = sqrt(sigma_int_PRS^2+sigma_int_ABS^2);
t = cst + x * sigma_ext2;

% get pA
pA = sum((sigma_ext2)*normpdf(t,cst,sigma_ext2).*(normcdf(t,0,sigma_int2).^2+(1-normcdf(t,0,sigma_int2)).^2)/100);

end
