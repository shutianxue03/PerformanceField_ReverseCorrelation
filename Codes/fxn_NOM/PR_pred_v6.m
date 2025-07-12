function pred = PR_pred_v6(iModelB, nBins, params_est, data, flag_plot)
% copied from fxn_getError_v5

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

%% predict pYES, pC and pA of each trial
IV_noNoise = data.IV;
e_noisy = data.e_noisy;
IV_noisy = IV_noNoise + sigma_template*e_noisy;
resp_data = data.resp;
iPRS = data.iPRS;
iPair = data.iPair;
corr = iPRS==resp_data;

% new sigma of each trial
sigma_pred = sqrt((IV_noisy*N_mul).^2 + SD_add^2);

% probablity of responding yes
pYES_pred = lambda/2+(1-lambda)*(1-normcdf(thresh, IV_noisy, sigma_pred));

% convert pYES to pC (accuracy): for ABS trials (iPRS==0), pC=1-pYES
pC_pred = pYES_pred;
pC_pred(iPRS==0) = 1-pYES_pred(iPRS==0);

% prob of consistent resp
pA_pred = pYES_pred.^2 + (1-pYES_pred).^2;

%% bin data to calculate predicted and measured pC and pA of each bin
% load binning info
nTrials_allBins = data.nTrials_allBins;
iTrial4Bin = data.iTrial4Bin;

for iind=1:length(iTrial4Bin)
    iBin_A=iTrial4Bin(iind);
    iTrialB=find(iPair==iPair(iind));
    iTrialB=iTrialB(2);
    iBin_B=iTrial4Bin(iTrialB);
    if iBin_A~=iBin_B
        fprintf('Trial%d: Bin%d and %d\n', iind, iBin_A, iBin_B) % confirming that two trials of one pair are in the same bins
    end
end

IV_allBins = nan(nBins, 1);
pC_pred_allBins = nan(nBins, 1);
pC_data_allBins = pC_pred_allBins;
pA_pred_allBins = pC_pred_allBins;
pA_data_allBins = pC_pred_allBins;
pYES_pred_allBins = pC_pred_allBins;
pYES_data_allBins = pC_pred_allBins;

for iBin=1:nBins
    IV_allBins(iBin) = mean(IV_noNoise(iTrial4Bin==iBin));
    
    pC_pred_allBins(iBin) = mean(pC_pred(iTrial4Bin==iBin));
    pC_data_allBins(iBin) = mean(corr(iTrial4Bin==iBin));
    
    pYES_pred_allBins(iBin) = mean(pYES_pred(iTrial4Bin==iBin));
    pYES_data_allBins(iBin) = mean(resp_data(iTrial4Bin==iBin));
    
    pA_pred_allBins(iBin) = mean(pA_pred(iTrial4Bin==iBin));
    
    %------------ get measured pA ------------
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