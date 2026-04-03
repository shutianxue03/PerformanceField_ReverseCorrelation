function pred = PR_pred_v8(iModelB, params_est, data, nBins, flag_plotPerIter)

% Predict pYES, pC, and pA (and their binned versions) for the trial-wise noisy observer model.

%% -------- 1. Extract basic data --------
IV_allT = data.IV(:); % trial-wise IV
resp_allT = data.resp(:); % trial-wise response (1=YES, 0=NO)
iPair_allT = data.iPair(:); % pair ID per trial
iPRS_allT = data.iPRS(:); % 1=PRS, 0=ABS
correctness_allT = (iPRS_allT == resp_allT); % trial-wise accuracy

nTrials_allBins = data.nTrials_allBins(:);
iTrial4Bin = data.iTrial4Bin(:);

nTrials = numel(IV_allT);
if any([numel(resp_allT), numel(iPair_allT), numel(iTrial4Bin)] ~= nTrials)
    error('PR_pred_v7: IV, resp, iPair, iTrial4Bin must have same length.');
end

%% -------- 2. Get pYES (trial-wise) and pA (pair-wise) from core function --------
% fxn_predMetrics does all the noisy-observer math (sigma, criterion, bivariate DV, etc.)
%------------------------------%
[pYES_pred_allT, pA_pred_allPairs, ~] = fxn_predMetrics_v3(iModelB, params_est, data);
%------------------------------%

% Sanity check
if numel(pYES_pred_allT) ~= nTrials
    error('PR_pred_v7: pYES_pred must be trial-wise (nTrials x 1).');
end

%% -------- 3. Map pair-wise pA back to trials --------
% pA_pair is defined over unique(iPair). Here we convert it to trial-wise pA_trial
% so that we can average pA within IV bins using trial membership.

iPairUnik = unique(iPair_allT);
nPairs = numel(iPairUnik);

if numel(pA_pred_allPairs) ~= nPairs
    error('PR_pred_v7: pA_pair must have length equal to numel(unique(iPair)).');
end

pA_trial = nan(nTrials, 1);
for iPair = 1:nPairs 
    thisPairID = iPairUnik(iPair);
    idx = (iPair_allT == thisPairID);
    if sum(idx) ~= 2
        error('PR_pred_v7: each iPair must contain exactly 2 trials.');
    end
    pA_trial(idx) = pA_pred_allPairs(iPair);
end

%% -------- 4. Convert pYES to accuracy pC --------
pC_pred_allT = nan(size(pYES_pred_allT));
pC_pred_allT(iPRS_allT == 1) = pYES_pred_allT(iPRS_allT == 1); % PRS trials
pC_pred_allT(iPRS_allT == 0) = 1 - pYES_pred_allT(iPRS_allT == 0); % ABS trials

%% -------- 5. Bin metrics as a function of IV --------
IV_allBins = nan(nBins, 1);
pYES_pred_allBins = nan(nBins, 1);
pYES_data_allBins = nan(nBins, 1);
pC_pred_allBins = nan(nBins, 1);
pC_data_allBins = nan(nBins, 1);
pA_pred_allBins = nan(nBins, 1);
pA_data_allBins = nan(nBins, 1);

for iBin = 1:nBins
    indTrial = (iTrial4Bin == iBin);
    if nTrials_allBins(iBin) ~= sum(indTrial)
        error('PR_pred_v7: inconsistent number of trials in bin %d.', iBin);
    end

    IV_allBins(iBin) = mean(IV_allT(indTrial));

    pYES_pred_allBins(iBin) = mean(pYES_pred_allT(indTrial));
    pYES_data_allBins(iBin) = mean(resp_allT(indTrial));

    pC_pred_allBins(iBin) = mean(pC_pred_allT(indTrial));
    pC_data_allBins(iBin) = mean(correctness_allT(indTrial));

    % Predicted pA: defined per trial (constant within each pair)
    pA_pred_allBins(iBin) = mean(pA_trial(indTrial));

    % Measured pA: compute agreement per pair within this bin
    resp_perBin = resp_allT(indTrial);
    iPair_perBin = iPair_allT(indTrial);
    iPair_perBin_unik = unique(iPair_perBin);
    nUnik = numel(iPair_perBin_unik);
    respC_perBin = nan(nUnik, 1);

    for iUnik = 1:nUnik
        idxPair = find(iPair_perBin == iPair_perBin_unik(iUnik)); 
        assert(numel(idxPair) == 2, 'PR_pred_v7: each pair in a bin must have 2 trials.');
        respC_perBin(iUnik) = (resp_perBin(idxPair(1)) == resp_perBin(idxPair(2)));
    end
    pA_data_allBins(iBin) = mean(respC_perBin);
end

%% -------- 6. Optional plotting --------
if flag_plotPerIter
    % Histograms per bin
    figure('Position', [0,0, 2000, 2000])
    for iBin = 1:nBins
        indTrial = (iTrial4Bin == iBin);

        subplot(nBins,3,(iBin-1)*3+1); hold on
        histogram(pYES_pred_allT(indTrial), 'Normalization','probability');
        xline(median(pYES_pred_allT(indTrial)), 'r-', 'LineWidth', 2);
        xline(mean(pYES_pred_allT(indTrial)), 'g-', 'LineWidth', 2);
        if iBin == nBins, xlabel('pYES'); end
        ylabel('Probability');
        ylim([0, .5]); xlim([0, 1]);

        subplot(nBins,3,(iBin-1)*3+2); hold on
        histogram(pC_pred_allT(indTrial), 'Normalization','probability');
        xline(median(pC_pred_allT(indTrial)), 'r-', 'LineWidth', 2);
        xline(mean(pC_pred_allT(indTrial)), 'g-', 'LineWidth', 2);
        if iBin == nBins, xlabel('pC'); end
        ylabel('Probability');
        ylim([0, .5]); xlim([0, 1]);

        subplot(nBins,3,(iBin-1)*3+3); hold on
        histogram(pA_trial(indTrial), 'Normalization','probability');
        xline(median(pA_trial(indTrial)), 'r-', 'LineWidth', 2);
        xline(mean(pA_trial(indTrial)), 'g-', 'LineWidth', 2);
        if iBin == nBins, xlabel('pA'); end
        ylabel('Probability');
        ylim([0, .5]); xlim([0.5, 1]);
    end

    % Trialwise metrics vs IV
    figure('Position', [0,0, 2000, 400])

    subplot(1,3,1); hold on
    plot(IV_allT(iPRS_allT==1), pYES_pred_allT(iPRS_allT==1), 'r.');
    plot(IV_allT(iPRS_allT==0), pYES_pred_allT(iPRS_allT==0), 'b.');
    plot(IV_allT(iPRS_allT==1), resp_allT(iPRS_allT==1), 'ro');
    plot(IV_allT(iPRS_allT==0), resp_allT(iPRS_allT==0), 'bo');
    xlabel('Trialwise IV'); ylabel('Probability / response');
    ylim([0, 1]); title('pYES');

    subplot(1,3,2); hold on
    plot(IV_allT(iPRS_allT==1), pC_pred_allT(iPRS_allT==1), 'r.');
    plot(IV_allT(iPRS_allT==0), pC_pred_allT(iPRS_allT==0), 'b.');
    plot(IV_allT(iPRS_allT==1), correctness_allT(iPRS_allT==1), 'ro');
    plot(IV_allT(iPRS_allT==0), correctness_allT(iPRS_allT==0), 'bo');
    xlabel('Trialwise IV'); ylabel('Probability / response');
    ylim([0, 1]); title('pC');

    subplot(1,3,3); hold on
    plot(IV_allT(iPRS_allT==1), pA_trial(iPRS_allT==1), 'r.');
    plot(IV_allT(iPRS_allT==0), pA_trial(iPRS_allT==0), 'b.');
    xlabel('Trialwise IV'); ylabel('Predicted pA');
    ylim([0, 1]); title('pA');

    sgtitle('Trialwise metrics as a function of IV');
end

%% -------- 7. Pack outputs --------
pred.metrics.nTrials_allBins = nTrials_allBins;
pred.metrics.IV_allBins = IV_allBins;

pred.metrics.pYES_pred_allBins = pYES_pred_allBins;
pred.metrics.pYES_data_allBins = pYES_data_allBins;
pred.R2_pYES = getR2(pYES_data_allBins, pYES_pred_allBins, nan);
pred.R2_weighted_pYES = getR2(pYES_data_allBins, pYES_pred_allBins, nTrials_allBins);

pred.metrics.pC_pred_allBins = pC_pred_allBins;
pred.metrics.pC_data_allBins = pC_data_allBins;
pred.R2_pC = getR2(pC_data_allBins, pC_pred_allBins, nan);
pred.R2_weighted_pC = getR2(pC_data_allBins, pC_pred_allBins, nTrials_allBins);

pred.metrics.pA_pred_allBins = pA_pred_allBins;
pred.metrics.pA_data_allBins = pA_data_allBins;
pred.R2_pA = getR2(pA_data_allBins, pA_pred_allBins, nan);
pred.R2_weighted_pA = getR2(pA_data_allBins, pA_pred_allBins, nTrials_allBins);

end