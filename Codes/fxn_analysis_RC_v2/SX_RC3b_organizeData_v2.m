
%% [Behav] measurements of each session
% empty containers
dprime_perSess_perLoc = nan(nSess, nLoc5);
criterion_perSess_perLoc = dprime_perSess_perLoc;
RT_perSess_perLoc = dprime_perSess_perLoc;
cst_perSess_perLoc = dprime_perSess_perLoc;
pA3_perSess_perLoc = nan(nSess, nLoc5, nTypes);
pC3_perSess_perLoc = pA3_perSess_perLoc; 

iCuedLoc = iCuedLoc_allT(:, 1);
for iLoc = 1:nLoc5
    dprime_perSess_perLoc(:, iLoc) = dprime_perBlk(iCuedLoc==iLoc);
    criterion_perSess_perLoc(:, iLoc) = criterion_perBlk(iCuedLoc==iLoc);
    RT_perSess_perLoc(:, iLoc) = RT_aveT_perBlk(iCuedLoc==iLoc);
    cst_perSess_perLoc(:, iLoc) = cst_perBlk(iCuedLoc == iLoc);
    pA3_perSess_perLoc(:, iLoc, :) = pA3_perBlk(:, iCuedLoc==iLoc)';
    pC3_perSess_perLoc(:, iLoc, :) = pC3_perBlk(:, iCuedLoc==iLoc)';
end

metrics_perSess_perLoc = cat(3, dprime_perSess_perLoc, criterion_perSess_perLoc, RT_perSess_perLoc, cst_perSess_perLoc, pC3_perSess_perLoc, pA3_perSess_perLoc);

%% [Behav] derive metrics per combined loc
nmetrics = 8;
metrics_allT_perComb = nan(nLoc8, nmetrics);
for iLoc = 1:nLoc8
    if iLoc > 5, switch iLoc, case 6, iLoc_all = [2,4]; case 7, iLoc_all = [3,5]; case 8, iLoc_all = 2:5; end
    else, iLoc_all = iLoc;
    end
    nLoc_ = length(iLoc_all);
    iPRS = []; resp = []; iPair = [];
    for iiLoc = 1:nLoc_
        iPRS = [iPRS; dataMatrix(dataMatrix(:, 5) == iLoc_all(iiLoc), 6)];
        resp = [resp; dataMatrix(dataMatrix(:, 5) == iLoc_all(iiLoc), 9)];
        iPair = [iPair; dataMatrix(dataMatrix(:, 5) == iLoc_all(iiLoc), 8)];
    end
    npairs = length(iPair); iPair_unik = unique(iPair); npairs_unik = length(iPair_unik);
    % pC
    pHit = mean(resp & iPRS)*2; 
    pFA = mean(resp & 1-iPRS)*2;
    pC = (pHit+1-pFA)/2;
    
    % dprime, criterion
    [d,c] = SX_sim06_SDT(pHit, pFA);
    
    % pA
    respC = nan(npairs_unik, 3); % to distinguish from resp
    for ipair_unik = 1:npairs_unik
        resp_PRS = resp(iPair == iPair_unik(ipair_unik) & iPRS); if isempty(resp_PRS), respC_PRS = nan; else, respC_PRS = resp_PRS(1)==resp_PRS(2); end
        resp_ABS = resp(iPair == iPair_unik(ipair_unik) & 1-iPRS);if isempty(resp_ABS), respC_ABS = nan; else, respC_ABS = resp_ABS(1)==resp_ABS(2); end
        resp_both = resp(iPair == iPair_unik(ipair_unik) );if isempty(resp_both), respC_both = nan; else, respC_both = resp_both(1)==resp_both(2); end
        respC(ipair_unik,:) = [respC_both, respC_PRS, respC_ABS];
    end
   
    metrics_allT_perComb(iLoc, :) = [d, c, pC, pHit, pFA, nanmean(respC)];
end

%% [Behav] get ave/sem_perComb from _perBlk
metrics_perBlk = [dprime_perBlk,criterion_perBlk, RT_aveT_perBlk, cst_perBlk, pC3_perBlk', pA3_perBlk'];

%% [Behav] RT median
% three categorizations
% 1. PRS vs. ABS (#Col in the dataMatrix is 6)
% 2. correct vs. wrong (#Col in the dataMatrix is dataMatrix(:, 6) == dataMatrix(:, 9))
% 3. say yes vs. no (#Col in the dataMatrix is 9)

RT_log_allT_perLoc = cell(3, 2, nLoc5);
RT_log_mean_perLoc = nan(3, 2, nLoc5);

RT_med_perLoc = nan(3, nLoc5, 2);
RT_std_perLoc = nan(3, nLoc5, 2);
for icat = 1:3
    switch icat
        case 1, ind = dataMatrix(:, 6);
        case 2, ind = dataMatrix(:, 6) == dataMatrix(:, 9);
        case 3, ind = dataMatrix(:, 9);
    end
    for iLoc = 1:nLoc5
        RT_perLoc1 = dataMatrix((dataMatrix(:, 5) == iLoc) & ind == 1, 10);
        RT_perLoc2 = dataMatrix((dataMatrix(:, 5) == iLoc) & ind == 0, 10);
        
        RT_med_perLoc(icat, iLoc, :) = [median(RT_perLoc1),median(RT_perLoc2)];
        RT_std_perLoc(icat, iLoc, :) = [std(RT_perLoc1),std(RT_perLoc2)];
    end
end

%% [patch-output] get CI
% if flagGetPatch
%     % load noise_allT by load(nameNoisePatch)
%     CI_mean = cell(1,nLoc8);
%     CI_var = CI_mean;
%     
%     for iLoc = 1:nLoc8
%         if iLoc>5, switch iLoc, case 6, iLoc_all = [2,4]; case 7, iLoc_all = [3,5]; case 8, iLoc_all = 2:5; end
%         else, iLoc_all = iLoc;
%         end
%         nLoc_ = length(iLoc_all);
%         iPRS = []; resp = []; iPair = [];
%         noise_Hit = []; noise_FA=[];noise_CR=[];noise_Miss=[];
%         for iiLoc = 1:nLoc_
%             % 6th column: iPRS; 9th column: resp
%             % noise_allT: a column of cells
%             noise_perLoc = noise_allT(dataMatrix(:,5)==iLoc_all(iiLoc));
%             iPRS_perLoc = dataMatrix(dataMatrix(:,5)==iLoc_all(iiLoc), 6);
%             iResp_perLoc = dataMatrix(dataMatrix(:,5)==iLoc_all(iiLoc), 9);
%             
%             noise_Hit = [noise_Hit; noise_perLoc(iPRS_perLoc & iResp_perLoc)];
%             noise_FA = [noise_FA; noise_perLoc(1-iPRS_perLoc & iResp_perLoc)];
%             noise_CR = [noise_CR; noise_perLoc(1-iPRS_perLoc & 1-iResp_perLoc)];
%             noise_Miss = [noise_Miss; noise_perLoc(iPRS_perLoc & 1-iResp_perLoc)];
%         end
%         
%         noise_Hit_ = nan([length(noise_Hit), size(noise_Hit{1})]); for ii = 1:length(noise_Hit), noise_Hit_(ii, :, :) = noise_Hit{ii};end
%         noise_FA_ = nan([length(noise_FA), size(noise_FA{1})]); for ii = 1:length(noise_FA), noise_FA_(ii, :, :) = noise_FA{ii};end
%         noise_CR_ = nan([length(noise_CR), size(noise_CR{1})]); for ii = 1:length(noise_CR), noise_CR_(ii, :, :) = noise_CR{ii};end
%         noise_Miss_ = nan([length(noise_Miss), size(noise_Miss{1})]); for ii = 1:length(noise_Miss), noise_Miss_(ii, :, :) = noise_Miss{ii};end
%         
%         % get CI
%         CI_mean{iLoc} = squeeze(mean(noise_Hit_,1) + mean(noise_FA_,1) - mean(noise_CR_,1) - mean(noise_Miss_,1));
%         CI_var{iLoc} = squeeze(var(noise_Hit_,[],1) + var(noise_FA_,[],1) - var(noise_CR_,[],1) - var(noise_Miss_,[],1));
%         
%     end % iLoc
% end

%% [patch-output] check power spectrum of patches
% if flagGetPatch
%     RC_checkNoiseFFT % output is fft1D_mean_both, 2x96 vector, should take 2:13 for SF = 1-4 cpd
% end

%% save patch & related outputs
% if flagGetPatch
%     fprintf(' Saving patches...')
%     save(namePatchOutputs, 'subjName', 'CI*', 'fft1D_mean_both') % search for x_record
%     fprintf('  DONE\n')
% end

%%
clear noise* patch* pp* nn* nfilters*
clear filters* kernels* nparams*  pV* R2* % do NOT delete isubj, subjList, nBlocks (thou unwanted)
clear repInd_allT answer_allT correctness_allT tgtPrs_allT 

%% save data
nameFile_BehavMeas = sprintf('%s/Data_OOD_%d%d/%s%d/%s_behavMeas.mat', nameFolder_Data, nORI, nSF, subjName, nBlocks, subjName);
fprintf('  Updating behav data...')
tic
save(nameFile_BehavMeas, '*_perLoc', '*_perBlk', '*_perComb', '-append')
dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)



