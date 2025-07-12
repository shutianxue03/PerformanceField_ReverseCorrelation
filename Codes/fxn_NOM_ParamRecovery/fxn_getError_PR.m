
function nLL = fxn_getError_PR(iModelB, params_est, data)

% implementing the alternative trial-wise model suggested by Mike Landy
switch iModelB
    case 1 % multiplicative, additive noise, thresh
        lambda = params_est(1);
        N_mul = params_est(2);
        SD_add = params_est(3);
        crit = params_est(4);
    case 2 % additive noise, thresh
        N_mul = 0;
        lambda = params_est(1);
        SD_add = params_est(2);
        crit = params_est(3);
    case 3 % multiplicative noise, thresh
        SD_add = 0;
        lambda = params_est(1);
        N_mul = params_est(2);
        crit = params_est(3);
    case 4
        N_mul=0;
        SD_add=0;
        lambda = params_est(1);
        crit = params_est(2);
end

%% get IV and resp, equivalent to our 'DATA'
IV_PRS = data.IV_PRS;
IV_ABS = data.IV_ABS;
% IV_noisy_PRS = data.IV_noisy_PRS;
% IV_noisy_ABS = data.IV_noisy_ABS;
resp_PRS = data.resp_PRS;
resp_ABS = data.resp_ABS;

% predict the width of the new distribution
sigma_pred_PRS = sqrt((IV_PRS*N_mul).^2 + SD_add^2);
sigma_pred_ABS = sqrt((IV_ABS*N_mul).^2 + SD_add^2);

% predict probablity of responding yes
pYES_pred_PRS = lambda/2+(1-lambda)*(1-normcdf(crit, IV_PRS, sigma_pred_PRS));
pYES_pred_ABS = lambda/2+(1-lambda)*(1-normcdf(crit, IV_ABS, sigma_pred_ABS));

% calculate nLL by fitting the predicted p(YES) to binary responses
% fit p(YES) to reasured response
nLL_PRS = -sum( resp_PRS.* log(pYES_pred_PRS) + (1 - resp_PRS) .* log(1 - pYES_pred_PRS));
nLL_ABS = -sum( resp_ABS.* log(pYES_pred_ABS) + (1 - resp_ABS) .* log(1 - pYES_pred_ABS));

nLL = nLL_PRS+nLL_ABS;


%%
% figure, hold on
% plot(IV_PRS(resp_PRS==1), pYES_pred_PRS(resp_PRS==1), 'ro')
% plot(IV_PRS(resp_PRS==0), pYES_pred_PRS(resp_PRS==0), 'bo')
% xline(crit, 'k-')



