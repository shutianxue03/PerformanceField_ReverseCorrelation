
function nLL = fxn_getError_v6(iModelB, params_est, data)

% implementing the alternative trial-wise model suggested by Mike Landy

switch iModelB
    case 1 % multiplicative, additive noise, thresh
        lambda = params_est(1);
        sigma_template = params_est(2); %N_mul = params_est(2);
        SD_add = params_est(3);
        thresh = params_est(4);
        N_mul=0;
%     case 2 % additive noise, thresh
%         N_mul = 0;
%         lambda = params_est(1);
%         SD_add = params_est(2);
%         thresh = params_est(3);
%     case 3 % multiplicative noise, thresh
%         SD_add = 0;
%         lambda = params_est(1);
%         N_mul = params_est(2);
%         thresh = params_est(3);
%     case 4
%         N_mul=0;
%         SD_add=0;
%         lambda = params_est(1);
%         thresh = params_est(2);
end

% load data
IV_noNoise = data.IV;
e_noisy = data.e_noisy;
resp_data = data.resp;

% calculate noisy IV (=energy * noisy template)
IV_noisy = IV_noNoise + sigma_template*e_noisy;

% calculate the new sigma of each trial
sigma_pred = sqrt((IV_noisy*N_mul).^2 + SD_add^2);

% calculate the probablity of responding yes
pYES_pred = lambda/2+(1-lambda)*(1-normcdf(thresh, IV_noisy, sigma_pred));

% calculate nLL by fitting predicted p(YES) to binary responses
nLL_YES = -sum(resp_data.* log(pYES_pred) + (1 - resp_data) .* log(1 - pYES_pred));

nLL = nLL_YES;
