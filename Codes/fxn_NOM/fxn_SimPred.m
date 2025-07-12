%% helper fxn - fxn_SimPred
function metrics = fxn_SimPred(ntrials_sim, mu_IV_PRS, std_IV_PRS, mu_IV_ABS, std_IV_ABS, params, plotFlag)
% simulate double-pass trials and predict metrics

nparams_model = length(params);
switch nparams_model
    case 1
        alpha_PRS = params(1);
        alpha_ABS = alpha_PRS;
        thresh = 0;
    case 2,
        alpha_PRS = params(1);
        alpha_ABS = alpha_PRS;
        thresh = params(2);
    case 3
        alpha_PRS = params(1);
        alpha_ABS = params(2);
        thresh = params(3);
end

%% derive the width of the new IV distribution, with internal noise added
std_IE_PRS = sqrt(1+alpha_PRS^2) * std_IV_PRS;
std_IE_ABS = sqrt(1+alpha_ABS^2) * std_IV_ABS;

%% simulate trials
IV_PRS = randn(1, ntrials_sim) * std_IV_PRS + mu_IV_PRS;
IV_ABS = randn(1, ntrials_sim) * std_IV_ABS + mu_IV_ABS;

%% adding internal noise
IV_PRS_noisy = randn(1, ntrials_sim) * std_IE_PRS + mu_IV_PRS;
IV_ABS_noisy = randn(1, ntrials_sim) * std_IE_ABS + mu_IV_ABS;

IV_PRS_noisy_sim = IV_PRS + randn(1, ntrials_sim)*std_IV_PRS*alpha_PRS;
IV_ABS_noisy_sim = IV_ABS + randn(1, ntrials_sim)*std_IV_ABS*alpha_ABS;

%% get the PRS-ABS difference 
% IV_diff = IV_PRS-IV_ABS;
% IV_diff_mean = mean(IV_diff);
% IV_diff_std = std(IV_diff);
% IV_diff_ = (IV_diff - IV_diff_mean)/IV_diff_std;

% pCDF = cdf(x, )
% pCDF = mean(IV_diff_>thresh);
% pA_math = 

%% standardize noisy IV
[IV_PRS_norm, IV_ABS_norm] = fxn_standardize(IV_PRS, IV_ABS);
[IV_PRS_noisy_norm, IV_ABS_noisy_norm] = fxn_standardize(IV_PRS_noisy, IV_ABS_noisy);
[IV_PRS_noisy_sim_norm, IV_ABS_noisy_sim_norm] = fxn_standardize(IV_PRS_noisy_sim, IV_ABS_noisy_sim);

%% quickPlot
if plotFlag
figure('Position', [500 200 800 400]), quickPlot_hist(IV_PRS, IV_ABS, IV_PRS_noisy, IV_ABS_noisy, IV_PRS_noisy_sim, IV_ABS_noisy_sim), sgtitle('Unstandardized')
figure('Position', [0 200 800 400]), quickPlot_hist(IV_PRS_norm, IV_ABS_norm, IV_PRS_noisy_norm, IV_ABS_noisy_norm, IV_PRS_noisy_sim_norm, IV_ABS_noisy_sim_norm), sgtitle('Standardized')
end
%% get internal responses (1=YES, 0=NO)
IR_PRS = IV_PRS_noisy_ >= thresh; % so that 1 is the correct answer
IR_ABS = IV_ABS_noisy_ >= thresh; % so that 0 is the correct answer

% split trials into two passes
IR_PRS_pass = reshape(IR_PRS, [2,ntrials_sim/2]);
IR_ABS_pass = reshape(IR_ABS, [2,ntrials_sim/2]);

% make predicts
% predict pHit/pFA/pC3
pHit = mean(IR_PRS);
pFA = mean(IR_ABS);
pC = (pHit+1-pFA)/2;

% predict dprime/criterion
[dprime , criterion] = SX_sim06_SDT(pHit, pFA, 1e5);

% predict pA3
pA_PRS = mean(IR_PRS_pass(1,:) == IR_PRS_pass(2,:));
pA_ABS = mean(IR_ABS_pass(1,:) == IR_ABS_pass(2,:));
pA = mean([pA_PRS, pA_ABS]);

% compile
metrics.IV_PRS = IV_PRS;
metrics.IV_ABS = IV_ABS;
metrics.dprime = dprime;
metrics.criterion = criterion;
metrics.pC = pC;
metrics.pHit = pHit;
metrics.pFA = pFA;
metrics.pA = pA;
metrics.pA_PRS = pA_PRS;
metrics.pA_ABS = pA_ABS;

end