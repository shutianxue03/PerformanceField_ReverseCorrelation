%% helper fxn - fxn_SimPred
function metrics_pred = fxn_pred(ntrials_pA, mu_IV_PRS, std_IV_PRS, mu_IV_ABS, std_IV_ABS, params)

% mathematically calculate pHit, pFA, pC, dprime, criterion by SDT theory
% calculate pA3 by simulation (at this point)

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

%% derive the width of the new IV distribution, with internal noise added
std_IE_PRS = sqrt(1+alpha_PRS^2) * std_IV_PRS;
std_IE_ABS = sqrt(1+alpha_ABS^2) * std_IV_ABS;

%% make predicts
% predict pHit/pFA/pC3
pHit_math = 1-normcdf(thresh, mu_IV_PRS, std_IE_PRS);
pFA_math = 1-normcdf(thresh, mu_IV_ABS, std_IE_ABS);
pC_math = (pHit_math+1-pFA_math)/2;
% predict dprime/criterion
[dprime_math , criterion_math] = SX_sim06_SDT(pHit_math, pFA_math, 1e5);

%% simulate trials to calculate pA
IV_PRS_pA = randn(1, ntrials_pA) * std_IE_PRS + mu_IV_PRS;
IV_ABS_pA = randn(1, ntrials_pA) * std_IE_PRS + mu_IV_ABS;

%get internal responses (1=YES, 0=NO)
IR_PRS = IV_PRS_pA >= thresh; % so that 1 is the correct answer
IR_ABS = IV_ABS_pA >= thresh; % so that 0 is the correct answer

% split trials into two passes
IR_PRS_pass = reshape(IR_PRS, [2,ntrials_pA/2]);
IR_ABS_pass = reshape(IR_ABS, [2,ntrials_pA/2]);

% predict pA3
pA_PRS_sim = mean(IR_PRS_pass(1,:) == IR_PRS_pass(2,:));
pA_ABS_sim = mean(IR_ABS_pass(1,:) == IR_ABS_pass(2,:));
pA_sim = mean([pA_PRS_sim, pA_ABS_sim]);

%% compile
metrics_pred.dprime = dprime_math;
metrics_pred.criterion = criterion_math;
metrics_pred.pC = pC_math;
metrics_pred.pHit = pHit_math;
metrics_pred.pFA = pFA_math;
metrics_pred.pA = pA_sim;
metrics_pred.pA_PRS = pA_PRS_sim;
metrics_pred.pA_ABS = pA_ABS_sim;

end