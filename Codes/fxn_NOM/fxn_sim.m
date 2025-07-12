

function metrics_sim = fxn_sim(ntrials_sim, mu_IV_PRS, std_IV_PRS, mu_IV_ABS, std_IV_ABS, params, plotFlag)
% simulate trials, then calculate ALL metrics 

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

% simulate raw trials (just for testing purpose)
IV_PRS = randn(1, ntrials_sim) * std_IV_PRS + mu_IV_PRS;
IV_ABS = randn(1, ntrials_sim) * std_IV_ABS + mu_IV_ABS;

% adding internal noise (just for testing purpose)
IV_PRS_noisy_sim = IV_PRS + randn(1, ntrials_sim)*std_IV_PRS*alpha_PRS;
IV_ABS_noisy_sim = IV_ABS + randn(1, ntrials_sim)*std_IV_ABS*alpha_ABS;

% adding derive new distribution given alpha
std_IE_PRS = sqrt(1+alpha_PRS^2) * std_IV_PRS;
std_IE_ABS = sqrt(1+alpha_ABS^2) * std_IV_ABS;
IV_PRS_noisy = randn(1, ntrials_sim) * std_IE_PRS + mu_IV_PRS;
IV_ABS_noisy = randn(1, ntrials_sim) * std_IE_ABS + mu_IV_ABS;

% standardize noisy IV (testing purpose)
[IV_PRS_norm, IV_ABS_norm] = fxn_standardize(IV_PRS, IV_ABS);
[IV_PRS_noisy_sim_norm, IV_ABS_noisy_sim_norm] = fxn_standardize(IV_PRS_noisy_sim, IV_ABS_noisy_sim);
[IV_PRS_noisy_norm, IV_ABS_noisy_norm] = fxn_standardize(IV_PRS_noisy, IV_ABS_noisy);

% quickPlot (testing purpose)
if plotFlag
figure('Position', [500 200 800 400]), quickPlot_hist(IV_PRS, IV_ABS, IV_PRS_noisy, IV_ABS_noisy, IV_PRS_noisy_sim, IV_ABS_noisy_sim), sgtitle('Unstandardized')
figure('Position', [0 200 800 400]), quickPlot_hist(IV_PRS_norm, IV_ABS_norm, IV_PRS_noisy_norm, IV_ABS_noisy_norm, IV_PRS_noisy_sim_norm, IV_ABS_noisy_sim_norm), sgtitle('Standardized')
end

% get internal responses (1=YES, 0=NO)
IR_PRS = IV_PRS_noisy_sim_norm >= thresh; % so that 1 is the correct answer
IR_ABS = IV_ABS_noisy_sim_norm >= thresh; % so that 0 is the correct answer

%% calculate metrics given simulated responses and threshold
pHit = mean(IR_PRS);
pFA = mean(IR_ABS);
pC = (pHit+1-pFA)/2;
 
% predict dprime/criterion
[dprime , criterion] = SX_sim06_SDT(pHit, pFA, 1e5);

% split trials into two passes
IR_PRS_pass = reshape(IR_PRS, [2,ntrials_sim/2]);
IR_ABS_pass = reshape(IR_ABS, [2,ntrials_sim/2]);

% predict pA3
pA_PRS = mean(IR_PRS_pass(1,:) == IR_PRS_pass(2,:));
pA_ABS = mean(IR_ABS_pass(1,:) == IR_ABS_pass(2,:));
pA = mean([pA_PRS, pA_ABS]);

% repeat the random sampling of trials to create the double pass 
% and generate the distribution of pA_PRS/ABS
% sim_randPass
% observation: the simulated pA is very consistent

% compile
metrics_sim.IV_PRS = IV_PRS; % the unstandardized IV distribution
metrics_sim.IV_ABS = IV_ABS;
metrics_sim.dprime = dprime;
metrics_sim.criterion = criterion;
metrics_sim.pC = pC;
metrics_sim.pHit = pHit;
metrics_sim.pFA = pFA;
metrics_sim.pA = pA;
metrics_sim.pA_PRS = pA_PRS;
metrics_sim.pA_ABS = pA_ABS;

end