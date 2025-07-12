
function pred = PR_pred_v5(iModelB, nBins, params_est, data, flag_plot)

% INPUT:
%   iModelB - Model type for internal noise structure and decision criterion
%   params_est - Estimated parameters for model
%   data - Data structure with fields for internal variables (IV) and responses
%   flag_plot - Flag to enable or disable plotting (unused in current code)

%% Model setup: Define parameters based on model type (iModelB)
switch iModelB
    case 1 % lapse rate, additive noise, criterion
        N_mul = 0;
        lambda = params_est(1);
        SD_add = params_est(2);
        criterion = params_est(3);
    case 2 % lapse rate, multiplicative noise, criterion
        SD_add = 0;
        lambda = params_est(1);
        N_mul = params_est(2);
        criterion = params_est(3);
    case 3 % lapse rate, criterion (no internal noise)
        N_mul=0;
        SD_add=0;
        lambda = params_est(1);
        criterion = params_est(2);
    case 4 % additive noise, criterion
        N_mul=0;
        SD_add=params_est(1);
        criterion = params_est(2);
    case 5 % multiplicative noise, criterion
        N_mul=params_est(1);
        SD_add=0;
        criterion = params_est(2);
end

%% Extract data variables
IV = data.IV;
resp_data = data.resp;
iPair = data.iPair;
iPRS = data.iPRS;
correctness = iPRS==resp_data; % correctness

%% Predict metrics
% Compute predicted standard deviation (sigma) for each trial
sigma_pred = sqrt((IV*N_mul).^2 + SD_add^2);

% Calculate probability of responding "Present" (pYES)
if iModelB<=3
    pYES_pred = lambda/2+(1-lambda)*(1-normcdf(criterion, IV, sigma_pred));
else
    pYES_pred = 1-normcdf(criterion, IV, sigma_pred);
end

% convert pYES to pC (accuracy): for ABS trials (iPRS==0), pC=1-pYES
pC_pred = pYES_pred;
pC_pred(iPRS==0) = 1-pYES_pred(iPRS==0);

% Predict probablity of responding consistently (pA)
pA_pred = pYES_pred.^2 + (1-pYES_pred).^2;

%% Bin data to calculate predicted and measured pC and pA for each bin
% load binning info
nTrials_allBins = data.nTrials_allBins;
iTrial4Bin = data.iTrial4Bin;

% Verify if paired trials are in the same bin
for iind=1:length(iTrial4Bin)
    iBin_A=iTrial4Bin(iind);
    iTrialB=find(iPair==iPair(iind));
    iTrialB=iTrialB(2);
    iBin_B=iTrial4Bin(iTrialB);
    if iBin_A~=iBin_B
        fprintf('Trial%d: Bin%d and %d\n', iind, iBin_A, iBin_B) % confirming that two trials of one pair are in the same bins
    end
end

% Initialize containers for binned values
IV_allBins = nan(nBins, 1);
pC_pred_allBins = nan(nBins, 1);
pC_data_allBins = pC_pred_allBins;
pA_pred_allBins = pC_pred_allBins;
pA_data_allBins = pC_pred_allBins;
pYES_pred_allBins = pC_pred_allBins;
pYES_data_allBins = pC_pred_allBins;

% Calculate binned values
for iBin = 1:nBins
    IV_allBins(iBin) = mean(IV(iTrial4Bin==iBin)); % Average IV for bin (should this be weighted by nTrials?)
    
    pC_pred_allBins(iBin) = mean(pC_pred(iTrial4Bin==iBin)); % Predicted accuracy
    pC_data_allBins(iBin) = mean(correctness(iTrial4Bin==iBin)); % Measured accuracy
    
    pYES_pred_allBins(iBin) = mean(pYES_pred(iTrial4Bin==iBin)); % Predicted pYES
    pYES_data_allBins(iBin) = mean(resp_data(iTrial4Bin==iBin)); % Measured pYES
    
    pA_pred_allBins(iBin) = mean(pA_pred(iTrial4Bin==iBin)); % Predicted pA
    
    %------------ Calculate measured pA for each bin ------------
    resp_perBin = resp_data(iTrial4Bin==iBin);
    iPair_perBin = iPair(iTrial4Bin==iBin);
    iPair_perBin_unik = unique(iPair_perBin);
    nUnik = length(iPair_perBin_unik);
    respC_perBin = nan(nUnik, 1);
    for iUnik=1:nUnik
        iTrialAB=find(iPair_perBin==iPair_perBin_unik(iUnik));
        respC_perBin(iUnik) = resp_perBin(iTrialAB(1))==resp_perBin(iTrialAB(2));
    end
    pA_data_allBins(iBin) = mean(respC_perBin);
    %------------------------------------------
end % iBin

%% save
pred.metrics.nTrials_allBins = nTrials_allBins;
pred.metrics.IV_allBins = IV_allBins;

pred.metrics.pC_pred_allBins = pC_pred_allBins;
pred.metrics.pC_data_allBins = pC_data_allBins;
pred.R2_pC = getR2(pC_data_allBins, pC_pred_allBins, nan);
pred.R2_weighted_pC = getR2(pC_data_allBins, pC_pred_allBins, nTrials_allBins);

pred.metrics.pA_pred_allBins = pA_pred_allBins;
pred.metrics.pA_data_allBins = pA_data_allBins;
pred.R2_pA = getR2(pA_data_allBins, pA_pred_allBins, nan);
pred.R2_weighted_pA = getR2(pA_data_allBins, pA_pred_allBins, nTrials_allBins);

pred.metrics.pYES_pred_allBins = pYES_pred_allBins;
pred.metrics.pYES_data_allBins = pYES_data_allBins;
pred.R2_pYES = getR2(pYES_data_allBins, pYES_pred_allBins, nan);
pred.R2_weighted_pYES = getR2(pYES_data_allBins, pYES_pred_allBins, nTrials_allBins);

%%
if flag_plot
    sz_scale = 80;
    figure('Position', [0 200 800 300 ])
    subplot(1,3,1), hold on
    plot(IV_allBins, pC_pred_allBins, 'k-')
    for iBin=1:nBins, plot(IV_allBins(iBin), pC_data_allBins(iBin), 'ko', 'MarkerSize', nTrials_allBins(iBin)/sz_scale+5), end
    title('Accuracy')
    subplot(1,3,2), hold on
    plot(IV_allBins, pYES_pred_allBins, 'k-')
    for iBin=1:nBins, plot(IV_allBins(iBin), pYES_data_allBins(iBin), 'ko', 'MarkerSize', nTrials_allBins(iBin)/sz_scale+5), end
    title('pYES')
    subplot(1,3,3), hold on
    plot(IV_allBins, pA_pred_allBins, 'k-')
    for iBin=1:nBins, plot(IV_allBins(iBin), pA_data_allBins(iBin), 'ko', 'MarkerSize', nTrials_allBins(iBin)/sz_scale+5), end
    title('Resp. consistency')
end
end