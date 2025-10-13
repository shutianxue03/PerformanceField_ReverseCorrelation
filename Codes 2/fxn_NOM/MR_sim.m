

function sim = MR_sim(ntrials_sim, truth, params, plotFlag)

% extract params and truth
nparams_model = length(params);
switch nparams_model
    case 1
        alpha_PRS = params(1);
        alpha_ABS = alpha_PRS;
        thresh = 0;
    case 2
        alpha_PRS = params(1);
        alpha_ABS = alpha_PRS;
        thresh = params(2);
    case 3
        alpha_PRS = params(1);
        alpha_ABS = params(2);
        thresh = params(3);
end

mu_PRS = truth.mu_PRS;
sigma_PRS = truth.sigma_PRS;
mu_ABS = truth.mu_ABS;
sigma_ABS = truth.sigma_ABS;

% simulate raw trials 
IV_PRS = randn(1, ntrials_sim) * sigma_PRS + mu_PRS;
IV_ABS = randn(1, ntrials_sim) * sigma_ABS + mu_ABS;

% adding internal noise
IV_PRS_noisy_sim = IV_PRS + randn(1, ntrials_sim) * sigma_PRS * alpha_PRS;
IV_ABS_noisy_sim = IV_ABS + randn(1, ntrials_sim) * sigma_ABS * alpha_ABS;

% derive new distribution given alpha
sigma_IE_PRS = sqrt(1+alpha_PRS^2) * sigma_PRS;
sigma_IE_ABS = sqrt(1+alpha_ABS^2) * sigma_ABS;
IV_PRS_noisy = randn(1, ntrials_sim) * sigma_IE_PRS + mu_PRS;
IV_ABS_noisy = randn(1, ntrials_sim) * sigma_IE_ABS + mu_ABS;

% standardize noisy IV (for testing purpose)
[IV_PRS_norm, IV_ABS_norm] = fxn_standardize(IV_PRS, IV_ABS);
[IV_PRS_noisy_sim_norm, IV_ABS_noisy_sim_norm] = fxn_standardize(IV_PRS_noisy_sim, IV_ABS_noisy_sim);
[IV_PRS_noisy_norm, IV_ABS_noisy_norm] = fxn_standardize(IV_PRS_noisy, IV_ABS_noisy);

% get internal responses (1=YES, 0=NO)
% (use IV_PRS/ABS_noisy_sim since IV_PRS/ABS are saved for fitting)
IR_PRS = IV_PRS_noisy_sim >= thresh; % so that 1 is the correct answer
IR_ABS = IV_ABS_noisy_sim >= thresh; % so that 0 is the correct answer

% quickPlot (testing purpose)
if plotFlag
    % histogram of data 
    figure('Position', [500 200 800 400])
    quickPlot_hist(IV_PRS, IV_ABS, IV_PRS_noisy, IV_ABS_noisy, IV_PRS_noisy_sim, IV_ABS_noisy_sim, mu_PRS, mu_ABS, sigma_IE_PRS, sigma_PRS, sigma_IE_ABS,  sigma_ABS)
    % resp vs. IV
    figure('Position', [500 200 600 200])
    subplot(1,2,1), hold on
    plot(IV_PRS_noisy_sim, IR_PRS, 'ro')
    [slope, intercept, R2, R2_Tjur, pValues, pCat, yfit] = SX_sim07_RC(1, 1, IV_PRS_noisy_sim', IR_PRS', 'logit');
    plot(IV_PRS_noisy_sim, yfit, 'k.')
    title(sprintf('slope = %.2f, itcpt = %.2f, R^2 = %.2f', slope, intercept, R2))
    subplot(1,2,2), hold on
    plot(IV_ABS_noisy_sim, IR_ABS, 'b.'), 
    %     figure('Position', [0 200 800 400]), quickPlot_hist(IV_PRS_norm, IV_ABS_norm, IV_PRS_noisy_norm, IV_ABS_noisy_norm, IV_PRS_noisy_sim_norm, IV_ABS_noisy_sim_norm), sgtitle('Standardized')
end

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

pred_sim = [dprime, criterion, pC, pHit, pFA, pA, pA_PRS, pA_ABS];

%% [test] obtain all 'simulated' metrics by math
pred_math = MR_predMetricsMath(mu_PRS, mu_ABS, sigma_PRS, sigma_ABS, alpha_PRS, alpha_ABS, thresh);

%% compile
sim.IV= [IV_PRS; IV_ABS];
sim.resp = [IR_PRS; IR_ABS]; % the unstandardized IV distribution
sim.metrics_sim = pred_sim;
sim.metrics_math = pred_math;

end