
function sim = PR_sim(iModelB, flagIncludePA, params_est, truth)
% function metrics_pred = MR_pred(mu_PRS, sigma_PRS, mu_ABS, sigma_ABS, params)

switch iModelB
    case 1 % multiplicative, additive noise, thresh
        N_mul = params_est(1);
        std_add = params_est(2);
        thresh = params_est(3);
        %         gamma = params_est(4);
    case 2 % additive noise, thresh
        N_mul = 0;
        std_add = params_est(1);
        thresh = params_est(2);
    case 3 % multiplicative noise, thresh
        N_mul = params_est(1);
        std_add = 0;
        thresh = params_est(2);
    case 4 % ratio, thresh
        ratio = params_est(1);
        thresh = params_est(2);
end

%% extract data
mu_PRS = truth.mu_PRS;
sigma_PRS = truth.sigma_PRS;
nPRS = truth.ndata_PRS;

mu_ABS = truth.mu_ABS;
sigma_ABS = truth.sigma_ABS;
nABS = truth.ndata_ABS;

%% predicts responses
IV_PRS = randn(nPRS, 1) * sigma_PRS + mu_PRS;
IV_ABS = randn(nABS, 1) * sigma_ABS + mu_ABS;
if iModelB==4
    IV_PRS_noiseAdd = IV_PRS + (ratio .* IV_PRS);
    IV_ABS_noiseAdd = IV_ABS + (ratio .* IV_ABS);
else
    %     IV_PRS_noiseAdd = IV_PRS + randn(nPRS, 1) * N_mul * sigma_PRS + randn(nPRS, 1) * std_add;
    %     IV_ABS_noiseAdd = IV_ABS + randn(nABS, 1) * N_mul * sigma_ABS + randn(nABS, 1) * std_add;
    
    IV_PRS_noiseAdd = IV_PRS + randn(nPRS, 1) * N_mul * sigma_PRS + randn(nPRS, 1) * std_add;
    IV_ABS_noiseAdd = IV_ABS + randn(nABS, 1) * N_mul * sigma_ABS + randn(nABS, 1) * std_add;
end

% quickPlot_IV_noiseAdded
%% get resp
resp = [ IV_PRS_noiseAdd > thresh, IV_ABS_noiseAdd > thresh];

% get measured metrics
pHit = mean(resp(:, 1));
pFA = mean(resp(:, 2));
pC = (pHit + 1-pFA)/2;
[d,c] = SX_sim06_SDT(pHit, pFA);

% split trials into two passes
resp_PRS_2pass = reshape(resp(:, 1), [2,nPRS/2]);
resp_ABS_2pass = reshape(resp(:, 2), [2,nABS/2]);

% predict pA3
respC_PRS = resp_PRS_2pass(1,:) == resp_PRS_2pass(2,:); pA_PRS = mean(respC_PRS);
respC_ABS = resp_ABS_2pass(1,:) == resp_ABS_2pass(2,:); pA_ABS = mean(respC_ABS);
pA = mean([pA_PRS, pA_ABS]);

sim_metrics = [d, c, pC, pHit, pFA, pA, pA_PRS, pA_ABS];

%% compile
sim.IV = [IV_PRS, IV_ABS];
sim.IV_noisy = [IV_PRS_noiseAdd, IV_ABS_noiseAdd];
sim.metrics_sim = sim_metrics;
sim.resp = resp;
sim.m_pA = (nPRS+nABS)/2; % number of all pairs
sim.n_pA = sum(respC_PRS) + sum(respC_ABS); % number of correct pairs
